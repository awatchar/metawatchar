#!/usr/bin/env python3
"""mu_eff retrieval of the corrected-copper S2 slab runs and comparison with the first runs
(local, light work on pulled results; needs numpy, scipy and the retrieval module
mu_eff_retrieval.py, path given with --retrieval).

The first S2 runs entered copper as Rs, Ls = Rs/w. Palace combines the two IN PARALLEL, so that was
Z = Rs (1 + j)/2, half the copper loss. The "_cu2" runs use Rs_p = 2 Rs, Ls_p = 2 Rs/w
(Z = Rs (1 + j) at 150 MHz).

Writes
  results/S_N<n>_cu2.csv            production S-parameters (order 3), columns as S_N<n>.csv
  results/S_N<n>_cu2.mu_eff.csv     retrieved mu, eps, n, z (physics convention, mu'' > 0 = loss)
  results/mu_eff_cu2_compare.json   Lorentz fits (old and corrected, same fit window) and
                                    mu', mu'' at selected frequencies, old vs corrected

Usage: python retrieve_cu2.py --retrieval <dir containing mu_eff_retrieval.py> [--N 1 2 4]
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
A_X = 0.070
FREQS_MHZ = [140.0, 145.0, 148.0, 149.5, 150.0, 150.5, 151.0, 155.0, 160.0, 170.0]


def load(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    return d["f_Hz"], d["S11_re"] + 1j * d["S11_im"], d["S21_re"] + 1j * d["S21_im"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--retrieval", required=True)
    ap.add_argument("--N", type=int, nargs="+", default=[1, 2, 4])
    ap.add_argument("--fit-window", type=float, nargs=2, default=[None, None], help="Hz (default: whole sweep, as S_N1.lorentz_fit.json)")
    a = ap.parse_args()
    sys.path.insert(0, a.retrieval)
    import mu_eff_retrieval as mr

    out = dict(note="mu in the physics convention (mu'' > 0 = loss); Lorentz fit "
                    "mu = mu_inf - F f^2/(f^2 - f0^2 + i f f0/Q) over the fit window (None = whole sweep), |mu| < 60",
               fit_window_Hz=a.fit_window, cases={})
    for n in a.N:
        src = RES / f"N{n}_z_cu2_o3_S.csv"
        if not src.exists():
            print("missing", src)
            continue
        dst = RES / f"S_N{n}_cu2.csv"
        with open(src) as fi, open(dst, "w", newline="") as fo:
            r = csv.DictReader(fi)
            w = csv.writer(fo)
            w.writerow(["f_Hz", "S11_re", "S11_im", "S21_re", "S21_im"])
            for row in r:
                w.writerow([row["f_Hz"], row["S11_re"], row["S11_im"], row["S21_re"], row["S21_im"]])
        L = n * A_X
        rec = {}
        for tag, path in (("old", RES / f"S_N{n}.csv"), ("cu2", dst)):
            f, s11, s21 = load(path)
            res = mr.retrieve(f, s11, s21, L, convention="engineering")
            if tag == "cu2":
                np.savetxt(RES / f"S_N{n}_cu2.mu_eff.csv",
                           np.column_stack([res["f"], res["mu"].real, res["mu"].imag, res["eps"].real,
                                            res["eps"].imag, res["n"].real, res["n"].imag, res["z"].real,
                                            res["z"].imag]),
                           delimiter=",", header="f_Hz,mu_re,mu_im,eps_re,eps_im,n_re,n_im,z_re,z_im", comments="")
            fit = mr.fit_lorentz(res["f"], res["mu"], *a.fit_window)
            mu = res["mu"]
            pts = {}
            for fm in FREQS_MHZ:
                i = int(np.argmin(np.abs(res["f"] - fm * 1e6)))
                pts[f"{fm:g}"] = dict(mu_re=float(mu[i].real), mu_im=float(mu[i].imag))
            mag = np.abs(s21)
            i0 = int(np.argmin(mag))
            neg = res["f"][np.isfinite(mu) & (mu.real < 0)]
            rec[tag] = dict(lorentz={k: float(v) for k, v in fit.items() if k != "model"},
                            mu_at_MHz=pts, mu_im_max=float(np.nanmax(mu.imag)),
                            f_mu_im_max_Hz=float(res["f"][np.nanargmax(mu.imag)]),
                            S21_min_dB=float(20 * np.log10(mag[i0])), f_S21_min_Hz=float(f[i0]),
                            absorption_max=float(np.max(1 - np.abs(s11) ** 2 - mag ** 2)),
                            mng_band_Hz=[float(neg[0]), float(neg[-1])] if neg.size else None)
        out["cases"][f"N{n}"] = rec
        o, c = rec["old"], rec["cu2"]
        print(f"N={n}: |S21|min {o['S21_min_dB']:.2f} -> {c['S21_min_dB']:.2f} dB "
              f"({o['f_S21_min_Hz']/1e6:.3f} -> {c['f_S21_min_Hz']/1e6:.3f} MHz); "
              f"abs max {o['absorption_max']:.4f} -> {c['absorption_max']:.4f}; "
              f"Lorentz Q {o['lorentz']['Q']:.0f} -> {c['lorentz']['Q']:.0f}, f0 {o['lorentz']['f0']/1e6:.3f} -> "
              f"{c['lorentz']['f0']/1e6:.3f}, F {o['lorentz']['F']:.4f} -> {c['lorentz']['F']:.4f}; "
              f"mu'' max {o['mu_im_max']:.1f} -> {c['mu_im_max']:.1f}")
    (RES / "mu_eff_cu2_compare.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
