"""Unloaded-Q extraction for a two-port transmission resonator (LGR-A meta-atom in its can).

Canonical copy (the measurement chain imports it from here).

Model of the transmission response near an isolated resonance (Kajfez, "Q Factor Measurements
Using MATLAB", 2011; Petersan & Anlage, J. Appl. Phys. 84, 3392, 1998):
    S21(f) = A * exp(j*phi) / (1 + 2j*Q_L*(f - f0)/f0) + (b0 + b1*(f - f0))     (complex)
The complex background (b0, b1) absorbs cable/leakage coupling between the loops. The fit is a
complex nonlinear least squares on the measured S21 (real and imaginary parts).

Coupling and unloaded Q:
  - Symmetric coupling (equal loops):         Qu = Q_L / (1 - |S21(f0)|)
  - General, using |S11|, |S22| at resonance:  beta1 = (1 - |S11|)/(|S11| + |S22|),
                                               beta2 = (1 - |S22|)/(|S11| + |S22|),
                                               Qu = Q_L * (1 + beta1 + beta2)
    (magnitude formulas; assume the detuned reflection magnitude is ~1, which holds for the
    weak coupling used here, |S21(f0)| <= -30 dB).
Uncertainty: 1-sigma from the fit Jacobian (residual-scaled covariance), plus a residual bootstrap.

CLI:  python q_extraction.py --csv meas.csv [--span-lw 8]   (columns f_Hz,S11_re,S11_im,S21_re,S21_im[,S22_re,S22_im])
      python q_extraction.py --selftest
"""
import argparse
import math
import sys

import numpy as np
from scipy.optimize import least_squares


def model(p, f):
    f0, ql, a, phi, b0r, b0i, b1r, b1i = p
    x = 2 * ql * (f - f0) / f0
    return a * np.exp(1j * phi) / (1 + 1j * x) + (b0r + 1j * b0i) + (b1r + 1j * b1i) * (f - f0)


def initial_guess(f, s21):
    i = int(np.argmax(np.abs(s21)))
    f0 = f[i]
    peak = np.abs(s21[i])
    half = peak / math.sqrt(2)
    above = np.flatnonzero(np.abs(s21) >= half)
    bw = max(f[above[-1]] - f[above[0]], f[1] - f[0]) if above.size > 1 else (f[-1] - f[0]) / 20
    return [f0, f0 / bw, peak, np.angle(s21[i]), 0, 0, 0, 0]


def fit_s21(f, s21, span_lw=None):
    f = np.asarray(f, float)
    s21 = np.asarray(s21, complex)
    p0 = initial_guess(f, s21)
    if span_lw:  # restrict the fit to +/- span_lw linewidths around the peak
        lw = p0[0] / p0[1]
        m = np.abs(f - p0[0]) <= span_lw * lw
        f, s21 = f[m], s21[m]
    scale = np.abs(s21).max()

    def resid(p):
        r = (model(p, f) - s21) / scale
        return np.concatenate([r.real, r.imag])

    lo = [f.min(), 1.0, 0, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf]
    hi = [f.max(), 1e7, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]
    sol = least_squares(resid, p0, bounds=(lo, hi), x_scale="jac", max_nfev=20000)
    J = sol.jac
    dof = max(1, 2 * f.size - len(p0))
    s2 = 2 * sol.cost / dof
    try:
        cov = np.linalg.inv(J.T @ J) * s2
        sig = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        sig = np.full(len(p0), np.nan)
    f0, ql, a = sol.x[0], sol.x[1], sol.x[2]
    return dict(f0=f0, QL=ql, S21_peak=a, sig_f0=sig[0], sig_QL=sig[1], sig_peak=sig[2], params=sol.x,
                f=f, s21=s21, rms=math.sqrt(s2) * scale)


def unloaded_q(fit, s11_f0=None, s22_f0=None):
    ql, t = fit["QL"], fit["S21_peak"]
    out = dict(Qu_symmetric=ql / (1 - t))
    if s11_f0 is not None and s22_f0 is not None:
        g1, g2 = abs(s11_f0), abs(s22_f0)
        b1, b2 = (1 - g1) / (g1 + g2), (1 - g2) / (g1 + g2)
        out.update(beta1=b1, beta2=b2, Qu_general=ql * (1 + b1 + b2))
    return out


def bootstrap(fit, n=200, rng=None):
    rng = rng or np.random.default_rng(0)
    f, s21 = fit["f"], fit["s21"]
    res = s21 - model(fit["params"], f)
    qus, f0s = [], []
    for _ in range(n):
        idx = rng.integers(0, f.size, f.size)
        fb = fit_s21(f, model(fit["params"], f) + res[idx])
        qus.append(unloaded_q(fb)["Qu_symmetric"])
        f0s.append(fb["f0"])
    return np.std(qus), np.std(f0s)


def synth(f, f0, qu, beta, noise_db, leak=1e-4, rng=None):
    """Synthetic symmetric transmission resonator with leakage and complex noise."""
    rng = rng or np.random.default_rng(1)
    ql = qu / (1 + 2 * beta)
    t0 = 2 * beta / (1 + 2 * beta)
    x = 2 * ql * (f - f0) / f0
    s21 = t0 * np.exp(1j * 0.7) / (1 + 1j * x) + leak * np.exp(1j * 2.1)
    r = (1 - 2 * beta / (1 + 2 * beta) * 1 / (1 + 1j * x))  # reflection at each port (approx.)
    sig = 10 ** (noise_db / 20)
    nz = lambda: sig * (rng.standard_normal(f.size) + 1j * rng.standard_normal(f.size)) / math.sqrt(2)
    return s21 + nz(), r + nz(), r + nz()


def selftest():
    f0, qu = 145.33e6, 3930
    f = np.linspace(f0 - 400e3, f0 + 400e3, 1601)
    print("| beta (each port) | |S21(f0)| | noise | Qu true | Qu (symmetric) | Qu (general) | f0 error (Hz) |")
    print("|---|---|---|---|---|---|---|")
    for beta in (0.005, 0.016, 0.05):
        for noise_db in (-100, -80, -70):
            s21, s11, s22 = synth(f, f0, qu, beta, noise_db)
            fit = fit_s21(f, s21)
            i0 = int(np.argmin(np.abs(f - fit["f0"])))
            q = unloaded_q(fit, s11[i0], s22[i0])
            t0 = 2 * beta / (1 + 2 * beta)
            print(f"| {beta} | {20*math.log10(t0):.1f} dB | {noise_db} dB | {qu} | {q['Qu_symmetric']:.0f} | "
                  f"{q['Qu_general']:.0f} | {fit['f0']-f0:+.1f} |")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--csv")
    ap.add_argument("--span-lw", type=float, default=None)
    ap.add_argument("--bootstrap", type=int, default=0)
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    d = np.genfromtxt(a.csv, delimiter=",", names=True)
    f = d["f_Hz"]
    s21 = d["S21_re"] + 1j * d["S21_im"]
    fit = fit_s21(f, s21, a.span_lw)
    i0 = int(np.argmin(np.abs(f - fit["f0"])))
    s11 = d["S11_re"][i0] + 1j * d["S11_im"][i0] if "S11_re" in d.dtype.names else None
    s22 = d["S22_re"][i0] + 1j * d["S22_im"][i0] if "S22_re" in d.dtype.names else None
    q = unloaded_q(fit, s11, s22)
    print(f"f0 = {fit['f0']/1e6:.6f} MHz ± {fit['sig_f0']:.1f} Hz")
    print(f"Q_L = {fit['QL']:.1f} ± {fit['sig_QL']:.1f};  |S21(f0)| = {20*math.log10(fit['S21_peak']):.2f} dB")
    for k, v in q.items():
        print(f"{k} = {v:.4g}")
    if a.bootstrap:
        sq, sf = bootstrap(fit, a.bootstrap)
        print(f"bootstrap ({a.bootstrap}): sigma(Qu_symmetric) = {sq:.1f}, sigma(f0) = {sf:.1f} Hz")


if __name__ == "__main__":
    main()
