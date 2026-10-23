#!/usr/bin/env python3
"""Figures and summary numbers for the S2 fixture case (light work on pulled results; numpy,
matplotlib, scikit-rf; the retrieval module directory is passed with --retrieval).

  python analyze_fixture.py empty   --retrieval <dir>   -> figures/empty_fixture.png, results/empty_summary.json
  python analyze_fixture.py slab    --retrieval <dir>   -> figures/mu_fixture_vs_ideal.png, results/mu_summary.json
  python analyze_fixture.py probes                      -> figures/edge_probes.png, results/edge_summary.json

Conventions: Palace and the CSVs use exp(+j w t); mu_eff is reported in the physics convention
(mu = mu' + i mu'', mu'' > 0 = loss), as in ../s2_mng_slab.
"""
import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
FIG = HERE / "figures"
S2 = HERE.parent / "s2_mng_slab" / "results"
C0 = 299_792_458.0
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a5cd1", "#52514e", "#c2407a", "#7a8a2a"]
INK2, GRID = "#52514e", "#e4e3df"


def plt_setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def load_S(run):
    d = np.genfromtxt(RES / f"{run}_S.csv", delimiter=",", names=True)
    out = {"f": d["f_Hz"], "S11": d["S11_re"] + 1j * d["S11_im"], "S21": d["S21_re"] + 1j * d["S21_im"]}
    if "S22_re" in d.dtype.names:
        out["S22"] = d["S22_re"] + 1j * d["S22_im"]
        out["S12"] = d["S12_re"] + 1j * d["S12_im"]
    return out


def db(x):
    return 20 * np.log10(np.abs(x) + 1e-300)


def read_palace_csv(path):
    with open(path) as fh:
        hdr = [h.strip() for h in fh.readline().split(",")]
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    return hdr, np.atleast_2d(data)


def tdr_impedance(f, s11, z0=50.0):
    """Step-response impedance profile of a 2-port from S11 on a uniform grid f = k df (DC added
    as S11(0) = 0 for a thru; Hann window)."""
    df = f[1] - f[0]
    s = np.concatenate([[0.0], s11 * np.hanning(2 * s11.size + 1)[s11.size + 1:]])
    npad = 16 * 2 * (s.size - 1)                 # zero padding: finer time sampling of the same response
    h = np.fft.irfft(s, n=npad) * npad / (2 * (s.size - 1))
    step = np.cumsum(h) * (2 * (s.size - 1)) / npad
    t = np.arange(step.size) / (npad * df)
    return t, z0 * (1 + step) / (1 - step)


def empty(args):
    plt = plt_setup()
    out = {}
    runs = [("openfar_nx0_al_o2b", "open, Al"), ("openfar_nx0_pec_o2b", "open, PEC"), ("openfar35_nx0_pec_o2b", "open, PEC, box 3.5 m"), ("open_nx0_pec_o2", "open, PEC, small box (0.45 m)"),
            ("open_nx0_box15_pec_o2", "open, PEC, small box x1.5"), ("closed_nx0_al_o2", "closed, w = 344, Al"), ("closedw280_nx0_al_o2", "closed, w = 280, Al")]
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), dpi=140)
    for k, (run, lab) in enumerate(runs):
        if not (RES / f"{run}_S.csv").exists():
            continue
        s = load_S(run)
        f = s["f"]
        band = (f >= 130e6) & (f <= 200e6)
        loss = 1 - np.abs(s["S11"]) ** 2 - np.abs(s["S21"]) ** 2
        t, z = tdr_impedance(f, s["S11"])
        # group delay of S21 -> one-way delay; centre plane at half of it
        ph = np.unwrap(np.angle(s["S21"]))
        tau = -np.polyfit(2 * np.pi * f[(f > 20e6) & (f < 280e6)], ph[(f > 20e6) & (f < 280e6)], 1)[0]
        zc = float(np.interp(tau, t, z))          # round trip to the centre plane = tau
        rec = dict(worst_S11_dB_130_200=float(db(s["S11"][band]).max()),
                   f_worst_Hz=float(f[band][np.argmax(np.abs(s["S11"][band]))]),
                   worst_S11_dB_1_300=float(db(s["S11"][f >= 1e6]).max()),
                   S21_dB_130_200=[float(db(s["S21"][band]).min()), float(db(s["S21"][band]).max())],
                   loss_130_200=[float(loss[band].min()), float(loss[band].max())],
                   loss_at_150=float(np.interp(150e6, f, loss)), loss_at_200=float(np.interp(200e6, f, loss)),
                   loss_at_300=float(loss[-1]), delay_ns=tau * 1e9, tdr_Z_centre_ohm=zc,
                   reciprocity_max=float(np.max(np.abs(s["S21"] - s.get("S12", s["S21"])))),
                   symmetry_S11_S22_max=float(np.max(np.abs(s["S11"] - s.get("S22", s["S11"])))))
        out[run] = rec
        c = PAL[k]
        axs[0, 0].plot(f / 1e6, db(s["S11"]), color=c, lw=1.4, label=lab)
        axs[0, 1].plot(f / 1e6, db(s["S21"]), color=c, lw=1.4, label=lab)
        axs[1, 0].semilogy(f / 1e6, np.clip(loss, 1e-7, None), color=c, lw=1.4, label=lab)
        axs[1, 1].plot(t * 1e9, z, color=c, lw=1.4, label=lab)
        print(run, json.dumps({k_: (round(v, 5) if isinstance(v, float) else v) for k_, v in rec.items()}))
    for ax in axs.flat:
        ax.grid(True, color=GRID, lw=0.6)
        ax.tick_params(colors=INK2, labelsize=8)
    for ax in (axs[0, 0], axs[0, 1], axs[1, 0]):
        ax.axvspan(130, 200, color=PAL[0], alpha=0.06, lw=0)
        ax.set_xlabel("f (MHz)", color=INK2)
    axs[0, 0].axhline(-15, color=INK2, lw=0.8, ls=":")
    axs[0, 0].set_ylabel("|S11| (dB)")
    axs[0, 1].set_ylabel("|S21| (dB)")
    axs[0, 1].set_ylim(-1.0, 0.05)
    axs[1, 0].set_ylabel("1 − |S11|² − |S21|² (loss + radiation)")
    axs[1, 1].set_xlabel("round-trip time (ns)")
    axs[1, 1].set_ylabel("TDR impedance (Ω), 300 MHz bandwidth")
    axs[1, 1].set_xlim(0, 8)
    axs[1, 1].set_ylim(35, 60)
    axs[0, 0].legend(fontsize=7.5)
    fig.suptitle("Empty fixture (the 2x-thru): 50-Ω air microstrip, 300-mm constant-impedance tapers; "
                 "half model, lumped 50-Ω launch at the taper tips", fontsize=9.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "empty_fixture.png")
    (RES / "empty_summary.json").write_text(json.dumps(out, indent=1) + "\n")


def mng_metrics(f, mu):
    """Pole (max mu''), zero crossings of mu' (MNG band), values at selected frequencies."""
    ok = np.isfinite(mu)
    f, mu = f[ok], mu[ok]
    i = int(np.argmax(mu.imag))
    re = mu.real
    cross = [float(f[k] + (f[k + 1] - f[k]) * re[k] / (re[k] - re[k + 1]))
             for k in range(f.size - 1) if re[k] * re[k + 1] < 0 and abs(re[k] - re[k + 1]) < 50]
    neg = f[re < 0]
    pts = {}
    for fm in (140, 145, 155, 160, 165, 170, 175, 185):
        j = int(np.argmin(np.abs(f - fm * 1e6)))
        pts[str(fm)] = [float(mu[j].real), float(mu[j].imag)]
    return dict(f_pole_Hz=float(f[i]), mu_im_max=float(mu.imag[i]), zero_crossings_Hz=cross,
                mng_band_Hz=[float(neg[0]), float(neg[-1])] if neg.size else None, mu_at_MHz=pts)


def fit_series(f, mu_fix, mu_ideal, exclude):
    """1/mu_fix = a / mu_ideal + b, least squares over points outside the `exclude` windows."""
    m = np.isfinite(mu_fix) & np.isfinite(mu_ideal)
    for lo, hi in exclude:
        m &= ~((f >= lo) & (f <= hi))
    A = np.column_stack([1 / mu_ideal[m], np.ones(m.sum())])
    coef, *_ = np.linalg.lstsq(A, 1 / mu_fix[m], rcond=None)
    a, b = coef
    pc = np.sum(np.real(1 / mu_ideal[m] - 1) * np.real(1 / mu_fix[m] - 1)) / np.sum(np.real(1 / mu_ideal[m] - 1) ** 2)
    return dict(a=complex(a), b=complex(b), p_constrained=float(pc), n_points=int(m.sum()))


CASES = [("openfar_nx1_al_o2", 1, "open4", "open, 4 x 1"), ("openfar_nx2_al_o2", 2, "open4", "open, 4 x 2"),
         ("openfar_ny5_nx1_al_o2", 1, "open5", "open, 5 x 1"), ("closed_nx1_al_o2", 1, "closed4_lid300", "closed w 344, 4 x 1"),
         ("closedw280_nx1_al_o2", 1, None, "closed w 280, 4 x 1")]


def slab(args):
    plt = plt_setup()
    qs = json.loads((RES / "quasistatic_design.json").read_text(encoding="utf-8"))
    ideal = {}
    for n in (1, 2):
        d = np.genfromtxt(S2 / f"S_N{n}_cu2.mu_eff.csv", delimiter=",", names=True)
        ideal[n] = (d["f_Hz"], d["mu_re"] + 1j * d["mu_im"])
    out = {"ideal_N1": mng_metrics(*ideal[1]), "ideal_N2": mng_metrics(*ideal[2]), "cases": {}}
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.5), dpi=140)
    fi, mui = ideal[1]
    s2n1 = np.genfromtxt(S2 / "S_N1_cu2.csv", delimiter=",", names=True)
    axs[0, 0].plot(fi / 1e6, mui.real, color="k", lw=1.0, label="ideal unit cell (S2, N = 1)")
    axs[0, 1].plot(fi / 1e6, mui.imag, color="k", lw=1.0, label="ideal unit cell (S2, N = 1)")
    axs[1, 0].plot(s2n1["f_Hz"] / 1e6, db(s2n1["S21_re"] + 1j * s2n1["S21_im"]), color="k", lw=1.0, label="ideal unit cell, N = 1")
    for k, (run, nx, qv, lab) in enumerate(CASES):
        p_mu = RES / f"{run}_slab.mu_eff.csv"
        if not p_mu.exists():
            continue
        d = np.genfromtxt(p_mu, delimiter=",", names=True)
        sl = np.genfromtxt(RES / f"{run}_slab.csv", delimiter=",", names=True, skip_header=1)
        f = d["f_Hz"]
        mu = d["mu_re"] + 1j * d["mu_im"]
        s11 = sl["S11_re"] + 1j * sl["S11_im"]
        s21 = sl["S21_re"] + 1j * sl["S21_im"]
        loss = 1 - np.abs(s11) ** 2 - np.abs(s21) ** 2
        raw = load_S(run)
        fI, muI = ideal[nx]
        mu_id = np.interp(f, fI, muI.real) + 1j * np.interp(f, fI, muI.imag)
        rec = dict(fixture_embedded=mng_metrics(f, mu),
                   raw_S21_min_dB=float(db(raw["S21"]).min()), raw_f_S21_min_Hz=float(raw["f"][np.argmin(np.abs(raw["S21"]))]),
                   deembedded_min_loss=float(loss.min()), f_deembedded_min_loss_Hz=float(f[np.argmin(loss)]),
                   deembedded_min_loss_below_150MHz=float(loss[f < 150e6].min()))
        band = (f >= 135e6) & (f <= 148e6)
        r = mu[band].real / mu_id[band].real
        rec["re_mu_ratio_fixture_over_ideal_135_148"] = [float(r.min()), float(r.max())]
        fitr = fit_series(f, mu, mu_id, exclude=[(146e6, 200e6)])     # fit below the resonance only
        rec["series_fit_below_146MHz"] = dict(a=[fitr["a"].real, fitr["a"].imag], b=[fitr["b"].real, fitr["b"].imag],
                                             p_constrained=fitr["p_constrained"], n=fitr["n_points"])
        p_fw = fitr["p_constrained"]
        pred_fw = 1 / (p_fw / mu_id + 1 - p_fw)
        rec["series_model_with_fitted_p"] = mng_metrics(f, pred_fw)
        if qv:
            p_qs = qs[qv]["p_series"]
            rec["p_quasistatic"] = p_qs
            rec["qs_prediction"] = mng_metrics(f, 1 / (p_qs / mu_id + 1 - p_qs))
        out["cases"][run] = rec
        c = PAL[k]
        mm = np.where(np.abs(mu) < 40, mu, np.nan)
        axs[0, 0].plot(f / 1e6, mm.real, color=c, lw=1.3, label=lab)
        axs[0, 1].plot(f / 1e6, mm.imag, color=c, lw=1.3, label=lab)
        axs[1, 0].plot(raw["f"] / 1e6, db(raw["S21"]), color=c, lw=1.2, label=f"{lab} (fixture ports)")
        axs[1, 1].plot(f / 1e6, loss, color=c, lw=1.2, label=lab)
        print(run, json.dumps({kk: v for kk, v in rec.items() if kk not in ("qs_prediction", "series_model_with_fitted_p")})[:1200])
    for ax, yl in ((axs[0, 0], (-12, 25)), (axs[0, 1], (-0.5, 6)), (axs[1, 0], (-70, 2)), (axs[1, 1], (-0.1, 0.1))):
        ax.set_ylim(*yl)
        ax.set_xlim(135, 195)
        ax.grid(True, color=GRID, lw=0.6)
        ax.axhline(0, color=INK2, lw=0.6)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.set_xlabel("f (MHz)", color=INK2)
        ax.legend(fontsize=6.8)
    axs[0, 0].set_ylabel("Re mu (fixture-embedded, after P370 de-embedding)")
    axs[0, 1].set_ylabel("Im mu")
    axs[1, 0].set_ylabel("|S21| at the fixture ports (dB)")
    axs[1, 1].set_ylabel("1 - |S11|^2 - |S21|^2 of the de-embedded slab")
    fig.suptitle("Fixture-embedded mu_eff (Palace fixture model, IEEE P370 2x-thru, planes moved to the slab faces) "
                 "vs the ideal unit-cell mu_eff of S2", fontsize=9.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "mu_fixture_vs_ideal.png")
    (RES / "mu_summary.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")


def probe_by(run):
    """|B_y| phasors at every probe (excitation 1), sorted by frequency."""
    p = RES / f"{run}_probe-B.csv"
    with open(p, encoding="utf-8") as fh:
        hdr = [h.strip() for h in fh.readline().split(",")]
    d = np.genfromtxt(p, delimiter=",", skip_header=1)
    d = d[np.argsort(d[:, 0])]
    pts = json.loads((HERE / "runs" / f"{run}.json").read_text())["probes"]
    out = []
    for i in range(len(pts)):
        cr = [k for k, h in enumerate(hdr) if h.startswith(f"Re{{B_y[{i + 1}]")][0]
        out.append(d[:, cr] + 1j * d[:, cr + 1])
    return d[:, 0] * 1e9, pts, out


def probes(args):
    plt = plt_setup()
    out = {}
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.2), dpi=140)
    f, pts, B = probe_by("openfar_nx0_al_o2b")
    j = int(np.argmin(np.abs(f - 150e6)))
    prof = [(pt[1], float(abs(b[j]) / abs(B[0][j]))) for pt, b in zip(pts, B) if pt[0] == 0.0 and pt[2] == 35.0]
    out["empty_By_profile_150MHz"] = prof
    axs[0].plot([a for a, _ in prof], [b for _, b in prof], "o-", color=PAL[0], lw=1.2, ms=4)
    for lo in (12.5, 82.5):
        axs[0].axvspan(lo, lo + 45, color=PAL[1], alpha=0.12, lw=0)
    axs[0].axvline(172, color=INK2, lw=0.8, ls=":")
    axs[0].text(174, 0.5, "strip edge", fontsize=7, color=INK2)
    axs[0].set_xlabel("y (mm) at x = 0, z = 35 mm (loop spans shaded)")
    axs[0].set_ylabel("|B_y| / |B_y(y = 0)|, empty open line, 150 MHz")
    axs[0].set_ylim(0, 1.1)
    for k, (run, title, ys) in enumerate((("openfar_nx1_al_o2", "open, 4 × 1: outer / inner loop", (35.0, 105.0)),
                                           ("openfar_ny5_nx1_al_o2", "open, 5 × 1: loops / centre loop", (0.0, 70.0, 140.0)))):
        if not (RES / f"{run}_probe-B.csv").exists():
            continue
        f, pts, B = probe_by(run)
        idx = [next(i for i, pt in enumerate(pts) if pt[0] == 0.0 and pt[1] == y and pt[2] == 35.0) for y in ys]
        ref = B[idx[0]]
        ax = axs[1 + k]
        rec = {}
        for n, i in enumerate(idx[1:]):
            r = B[i] / ref
            ax.plot(f / 1e6, np.abs(r), color=PAL[n], lw=1.3, label=f"|B(y={ys[n + 1]:.0f})/B(y={ys[0]:.0f})|")
            ax.plot(f / 1e6, np.cos(np.angle(r)), color=PAL[n], lw=0.8, ls=":", label="cos(phase difference)")
            rec[f"y{ys[n + 1]:.0f}_over_y{ys[0]:.0f}"] = {str(fm): [float(np.abs(np.interp(fm * 1e6, f, r.real) + 1j * np.interp(fm * 1e6, f, r.imag))),
                                                                   float(np.degrees(np.angle(np.interp(fm * 1e6, f, r.real) + 1j * np.interp(fm * 1e6, f, r.imag))))]
                                                          for fm in (140, 145, 150, 154, 158, 160, 165, 170, 172, 175)}
        for i, y in zip(idx, ys):
            rec[f"peak_y{y:.0f}_Hz"] = float(f[np.argmax(np.abs(B[i]))])
        out[run] = rec
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("f (MHz)")
        ax.set_xlim(140, 185)
        ax.set_ylim(-1.2, 2.5)
        ax.axhline(1, color=INK2, lw=0.6)
        ax.legend(fontsize=7)
    for ax in axs:
        ax.grid(True, color=GRID, lw=0.6)
        ax.tick_params(colors=INK2, labelsize=8)
    fig.tight_layout()
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "edge_probes.png")
    (RES / "edge_summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out)[:2500])


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["empty", "slab", "probes"])
    ap.add_argument("--retrieval", default=str(HERE.parents[2] / "analysis"))
    args = ap.parse_args()
    {"empty": empty, "slab": slab, "probes": probes}[args.cmd](args)


if __name__ == "__main__":
    main()
