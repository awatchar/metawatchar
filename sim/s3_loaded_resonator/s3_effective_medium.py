#!/usr/bin/env python3
"""Effective-medium predictions for the S3 short-circuited line resonators (study 1).

Model: the S2 TEM unit-cell line (normalized impedance 1) holding a homogeneous slab of length
L = N a (wave impedance z, index n) between two empty sections of length g, short-circuited at
x = -g and x = L + g. The transverse-resonance condition is evaluated at the left short:

    Z1 = j tan(k0 g)                                   (right air section + right short)
    Z2 = z (Z1 + j z tan(n k0 L)) / (z + j Z1 tan(n k0 L))   (slab)
    Z3 = (Z2 + j tan(k0 g)) / (1 + j Z2 tan(k0 g))     (left air section)
    resonance:  Z3(f) = 0                              (series resonance at the left short)

Engineering convention exp(+j w t) throughout (n, z are the conjugates of the physics-convention
retrieval output), so lossy roots have Im f > 0 and Q = Re f / (2 Im f).

Two ways of giving the medium:
  A  "retrieved": n(f), z(f) retrieved (Smith/Chen inversion with a Kramers-Kronig branch choice)
     from the S2 slab S-parameters of N_ret = 1, 2, 4 cells (results/S_N*.csv of ../s2_mng_slab),
     on the real 10 kHz grid 130-200 MHz. Roots: X = Im Z3 crosses zero upwards; Q from the slope,
     Q = f0 (dX/df) / (2 R), with R = Re Z3 at the crossing.
  B  "Lorentz": mu(f) = mu_inf - F f^2 / (f^2 - f0^2 - j f f0 / Q_cell), eps constant, continued
     to complex f; complex roots by Newton/secant iteration. Two parameter sets:
       B-fit : f0, F, mu_inf from a least-squares fit to the retrieved mu (N_ret = 4), Q_cell from
               the meta-atom eigenmode;
       B-eig : from the meta-atom eigenmodes only: PMC-box f_m is the k = pi/a point, mu(f_m) =
               (pi c / (2 pi f_m a))^2 / eps; PEC-box f_mp is mu = 0; mu_inf = 1.
     eps = median of the retrieved Re eps over 130-145 and 155-200 MHz.

Usage (needs numpy and scipy, plus the retrieval module `mu_eff_retrieval.py`):
  python s3_effective_medium.py --retrieval-module <dir containing mu_eff_retrieval.py>
writes results/em_predictions.csv, results/em_parameters.json, results/retrieved_nz_N{1,2,4}.csv
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
S2 = HERE.parent / "s2_mng_slab" / "results"
C0 = 299_792_458.0
A = 0.070
# meta-atom eigenmodes, v2 mesh, order 2 (S3 half model, results/L1_g00_*_o2*_eig.csv):
#   "s2": copper entered as in S2 (Rs || Ls = Rs/w in Palace's parallel model = half the copper loss);
#         consistent with the S2 slab S-parameters, hence with model A
#   "cu": correct copper surface impedance Z = Rs (1 + j)
EIG = {"s2": dict(f_m=149.934059e6, Q_m=5035.4, f_mp=179.516422e6, Q_mp=4924.4),
       "cu": dict(f_m=149.921868e6, Q_m=2811.4, f_mp=179.503928e6, Q_mp=2746.0)}


def zin(f, n, z, N, g_m):
    k0 = 2 * np.pi * f / C0
    L = N * A
    t_g = np.tan(k0 * g_m)
    z1 = 1j * t_g
    tn = np.tan(n * k0 * L)
    z2 = z * (z1 + 1j * z * tn) / (z + 1j * z1 * tn)
    return (z2 + 1j * t_g) / (1 + 1j * z2 * t_g)


def roots_real_grid(f, n, z, N, g_m, fmin=None, fmax=None):
    Z = zin(f, n, z, N, g_m)
    X, R = Z.imag, Z.real
    out = []
    for i in range(len(f) - 1):
        if X[i] < 0 <= X[i + 1] and np.isfinite(X[i]) and np.isfinite(X[i + 1]):
            # a pole of Z3 crosses + -> - ; only - -> + crossings are series resonances
            t = -X[i] / (X[i + 1] - X[i])
            f0 = f[i] + t * (f[i + 1] - f[i])
            j0, j1 = max(i - 1, 0), min(i + 2, len(f) - 1)
            dxdf = (X[j1] - X[j0]) / (f[j1] - f[j0])
            r = R[i] + t * (R[i + 1] - R[i])
            if fmin and f0 < fmin or fmax and f0 > fmax:
                continue
            out.append(dict(f_MHz=f0 / 1e6, Q=f0 * dxdf / (2 * r) if r > 0 else float("inf"),
                            R=r, dXdf_perMHz=dxdf * 1e6, Xjump=max(abs(X[i]), abs(X[i + 1]))))
    return out


def lorentz_eng(f, f0, F, q, minf):
    return minf - F * f ** 2 / (f ** 2 - f0 ** 2 - 1j * f * f0 / q)


def nz_lorentz(f, p):
    mu = lorentz_eng(f, p["f0"], p["F"], p["Q"], p["mu_inf"])
    eps = p["eps"]
    n = np.sqrt(mu * eps)
    n = np.where(n.imag > 0, -n, n)          # engineering: Im n <= 0 (decay along +x)
    zz = np.sqrt(mu / eps)
    zz = np.where(zz.real < 0, -zz, zz)
    return n, zz


def complex_roots(p, N, g_m, fgrid):
    """Complex roots of Z3 for the Lorentz medium. Seeds: upward zero crossings of Im Z3 on a real
    grid with the lossless medium, then secant iteration on the lossy function. The result does not
    depend on the branch of n (Z3 is even in n), only on z = sqrt(mu/eps) with mu, eps analytic."""
    pl = dict(p, Q=1e12)
    n, zz = nz_lorentz(fgrid, pl)
    seeds = [r["f_MHz"] * 1e6 for r in roots_real_grid(fgrid, n, zz, N, g_m)]

    def D(fc):
        mu = lorentz_eng(fc, p["f0"], p["F"], p["Q"], p["mu_inf"])
        nn = np.sqrt(mu * p["eps"])
        zz_ = nn / p["eps"]                  # z = mu / n = n / eps (same branch as n: Z3 even in n)
        return zin(fc, nn, zz_, N, g_m)

    out = []
    for s in seeds:
        x0, x1 = complex(s, s / 2e4), complex(s * (1 + 1e-6), s / 2e4)
        d0, d1 = D(x0), D(x1)
        for _ in range(60):
            if d1 == d0:
                break
            x2 = x1 - d1 * (x1 - x0) / (d1 - d0)
            x0, d0, x1, d1 = x1, d1, x2, D(x2)
            if abs(x1 - x0) < 1e-3:
                break
        if not (abs(d1) < 1e-6 and x1.imag > 0):
            continue
        # Bloch-zone filter: a homogeneous Lorentz medium has unbounded k near f0, the lattice
        # does not (k <= pi/a). Keep roots with Re(n) k0 a < 0.97 pi; drop duplicates.
        mu_r = lorentz_eng(x1.real, p["f0"], p["F"], 1e12, p["mu_inf"])
        ka = abs(np.sqrt(complex(mu_r * p["eps"])).real) * 2 * np.pi * x1.real / C0 * A
        if mu_r.real > 0 and ka >= 0.97 * np.pi:
            continue
        if any(abs(o["f_MHz"] - x1.real / 1e6) < 1e-4 for o in out):
            continue
        out.append(dict(f_MHz=x1.real / 1e6, Q=x1.real / (2 * x1.imag), seed_MHz=s / 1e6, ka_over_pi=ka / np.pi))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--retrieval-module", default=str(HERE.parents[2] / "analysis"))
    a = ap.parse_args()
    sys.path.insert(0, a.retrieval_module)
    import mu_eff_retrieval as mr

    out_rows, params = [], {}
    ret = {}
    for Nr in (1, 2, 4):
        d = np.genfromtxt(S2 / f"S_N{Nr}.csv", delimiter=",", names=True)
        f = d["f_Hz"]
        r = mr.retrieve(f, d["S11_re"] + 1j * d["S11_im"], d["S21_re"] + 1j * d["S21_im"], Nr * A,
                        convention="engineering")
        ret[Nr] = (f, r)
        np.savetxt(HERE / "results" / f"retrieved_nz_N{Nr}.csv",
                   np.column_stack([f, r["n"].real, r["n"].imag, r["z"].real, r["z"].imag, r["mu"].real,
                                    r["mu"].imag, r["eps"].real, r["eps"].imag]), delimiter=",",
                   header="f_Hz,n_re,n_im,z_re,z_im,mu_re,mu_im,eps_re,eps_im (physics convention exp(-iwt))",
                   comments="", fmt="%.10g")

    # --- model A: retrieved medium on the real grid
    cases = [(N, g) for g in (0.0, 35.0) for N in (1, 2, 3, 4)]
    for (N, g) in cases:
        for Nr in (1, 2, 4):
            f, r = ret[Nr]
            n_e, z_e = np.conj(r["n"]), np.conj(r["z"])
            for k, rt in enumerate(roots_real_grid(f, n_e, z_e, N, g * 1e-3)):
                ka = abs(np.interp(rt["f_MHz"] * 1e6, f, r["n"].real)) * 2 * np.pi * rt["f_MHz"] * 1e6 / C0 * A
                out_rows.append(dict(N=N, g_mm=g, model=f"A-ret-N{Nr}", f_MHz=rt["f_MHz"], Q=rt["Q"],
                                     ka_over_pi=ka / np.pi))

    # --- model B: Lorentz medium
    f4, r4 = ret[4]
    band = ((f4 > 130e6) & (f4 < 145e6)) | ((f4 > 155e6) & (f4 < 200e6))
    eps = float(np.median(r4["eps"].real[band]))
    fit = mr.fit_lorentz(f4, r4["mu"], fmin=135e6, fmax=200e6)
    from scipy.optimize import brentq
    params = {"eps_eff_median_N4": eps, "lorentz_fit_raw_Q_ill_conditioned": float(fit["Q"]),
              "meta_atom_eigen": EIG}
    sets = []
    for tag, e in EIG.items():
        sfx = "" if tag == "s2" else "-cu"
        pfit = dict(f0=float(fit["f0"]), F=float(fit["F"]), mu_inf=float(fit["mu_inf"]), Q=e["Q_m"], eps=eps)
        # B-eig: mu(f_m) = (c/(2 f_m a))^2 / eps at the PMC-box point (k = pi/a), mu(f_mp) = 0, mu_inf = 1
        mu_m = (C0 / (2 * e["f_m"] * A)) ** 2 / eps

        def resid(f0):
            F = 1 - f0 ** 2 / e["f_mp"] ** 2
            return 1 - F * e["f_m"] ** 2 / (e["f_m"] ** 2 - f0 ** 2) - mu_m
        f0e = brentq(resid, e["f_m"] * (1 + 1e-9), e["f_m"] * 1.05)
        peig = dict(f0=f0e, F=1 - f0e ** 2 / e["f_mp"] ** 2, mu_inf=1.0, Q=e["Q_m"], eps=eps)
        params["B-fit" + sfx] = pfit
        params["B-eig" + sfx] = dict(peig, mu_at_fm=mu_m)
        sets += [("B-fit" + sfx, pfit), ("B-eig" + sfx, peig)]
    fgrid = np.linspace(110e6, 210e6, 200001)
    for (N, g) in cases:
        for name, p in sets:
            for rt in complex_roots(p, N, g * 1e-3, fgrid):
                out_rows.append(dict(N=N, g_mm=g, model=name, f_MHz=rt["f_MHz"], Q=rt["Q"],
                                     ka_over_pi=rt["ka_over_pi"]))

    with open(HERE / "results" / "em_predictions.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["N", "g_mm", "model", "f_MHz", "Q", "ka_over_pi"])
        w.writeheader()
        for r in sorted(out_rows, key=lambda r: (r["g_mm"], r["N"], r["model"], r["f_MHz"])):
            w.writerow({k: (f"{v:.6f}" if k == "f_MHz" else f"{v:.1f}" if k == "Q" else f"{v:.4f}" if k == "ka_over_pi" else v)
                        for k, v in r.items()})
    (HERE / "results" / "em_parameters.json").write_bytes((json.dumps(params, indent=2) + "\n").encode())
    for r in sorted(out_rows, key=lambda r: (r["g_mm"], r["N"], r["model"], r["f_MHz"])):
        if 120 <= r["f_MHz"] <= 200:
            print(f"N={r['N']} g={r['g_mm']:4.0f}  {r['model']:<9} f={r['f_MHz']:9.4f} MHz  Q={r['Q']:9.1f}  ka/pi={r['ka_over_pi']:.3f}")
    print(json.dumps(params, indent=1))


if __name__ == "__main__":
    main()
