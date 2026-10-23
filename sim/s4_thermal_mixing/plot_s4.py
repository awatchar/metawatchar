#!/usr/bin/env python3
"""Figures for S4 from results/s4_*.csv (local; matplotlib). Called by analyze_s4.py."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
RES, FIG = HERE / "results", HERE / "figures"
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"      # categorical slots 1-4
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
F_DUP = 145.2e6


def rows(name):
    p = RES / name
    if not p.exists():
        return []
    with open(p) as fh:
        out = []
        for r in csv.DictReader(fh):
            for k, v in r.items():
                try:
                    r[k] = float(v)
                except (TypeError, ValueError):
                    pass
            out.append(r)
        return out


def style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, color=INK2)
    ax.set_ylabel(ylabel, color=INK2)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2)


def tau(r):
    return r["tau_pf"]          # participation formula (exact first derivative; FD-checked)


def fig_tau():
    box, det, line = rows("s4_box_mixing.csv"), rows("s4_box_detuned.csv"), rows("s4_line.csv")
    if not box:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=150)
    style(ax, "TiO$_2$ share of the gap electric energy, $x_{TiO_2}$ (%)", r"$\tau_f$ of the working mode (ppm/K)")
    ax.axhline(0, color=INK2, lw=1)
    x = np.array([r["x_ti"] for r in box]) * 100
    y = np.array([tau(r) for r in box])
    c = np.polyfit(x, y, 1)
    xx = np.linspace(0, 100, 2)
    ax.plot(xx, np.polyval(c, xx), color=S1, lw=1, ls="--", alpha=0.6)
    ax.plot(x, y, "o", ms=8, color=S1, mec="white", mew=1.5, label="count mixing, k = 0...4 TiO$_2$ cells (2x1x2 box)")
    for r in box:
        if r["k"] in (0, 4) or r["label"].endswith("bottom row") and r["k"] == 1:
            ax.annotate(f"k = {int(r['k'])}", (r["x_ti"] * 100, tau(r)), textcoords="offset points",
                        xytext=(8, -4), color=INK2, fontsize=8)
    dd = [r for r in det if isinstance(r.get("detune"), float)]
    if dd:
        ax.plot([r["x_ti"] * 100 for r in dd], [tau(r) for r in dd], "s-", ms=8, lw=2, color=S2,
                mec="white", mew=1.5, label="k = 1, TiO$_2$ cell detuned up (continuous share)")
    comp = [r for r in det if "composite" in str(r.get("label"))]
    if comp:
        ax.plot([r["x_ti"] * 100 for r in comp], [tau(r) for r in comp], "D", ms=8, color=S4,
                mec="white", mew=1.5, label="composite-gap cells (meta-atom level, reference)")
    if line:
        ax.plot([r["x_ti"] * 100 for r in line], [tau(r) for r in line], "^", ms=8, color=S3,
                mec="white", mew=1.5, label="4-cell line (S3 study-1 geometry)")
    xz = -c[1] / c[0]
    ax.annotate(f"zero at {xz:.1f} %", (xz, 0), textcoords="offset points", xytext=(6, 10), color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "tau_vs_share.png")
    # zoom
    ax.set_xlim(0, 35)
    ax.set_ylim(-80, 80)
    fig.savefig(FIG / "tau_vs_share_zoom.png")
    plt.close(fig)


def fig_q():
    qt = [r for r in rows("s4_q_tand.csv") if "CCCC" not in r["config"]]   # composite: other f
    if not qt:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    style(ax, "TiO$_2$ share of the gap electric energy, $x_{TiO_2}$ (%)", "Q of the working mode")
    cols = {6e-5: S3, 1e-4: S1, 3e-4: S2, 1e-3: S4}
    for tt, col in cols.items():
        rr = sorted([r for r in qt if abs(r["tand_ti"] - tt) < 1e-12], key=lambda r: r["x_ti"])
        ax.plot([r["x_ti"] * 100 for r in rr], [r["Q"] for r in rr], "o", ms=7, color=col, mec="white",
                mew=1.2, label=f"Q, tan $\\delta$(TiO$_2$) = {tt:.0e}")
    rr = sorted([r for r in qt if abs(r["tand_ti"] - 1e-4) < 1e-12], key=lambda r: r["x_ti"])
    ax.plot([r["x_ti"] * 100 for r in rr], [r["Q_c_cu"] for r in rr], "x", ms=7, color=INK2, label="Q$_c$ (copper only)")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIG / "q_vs_share.png")
    plt.close(fig)


def fig_drift():
    """Working-mode frequency against temperature: FD points (permittivity part) plus the linear
    expansion terms, scaled to 145.2 MHz; best linear trim applied (tau_f -> optimum for +-20 K)."""
    import json
    from scipy.optimize import minimize_scalar
    s = json.loads((RES / "s4_summary.json").read_text())
    allrows = {r["config"]: r for k in ("box", "detuned", "composite") for r in s.get(k, [])}
    sys_path = str(HERE)
    import sys
    sys.path.insert(0, sys_path)
    import analyze_s4 as an
    picks = [("all alumina (untrimmed)", "B212_AAAA_o2", S1, False),
             ("1 TiO2 cell, detuned (W3.41), best linear trim", "B212_TAAA_W3p41_o2", S2, True),
             ("1 TiO2 cell, matched (W3.29), linear part removed", "B212_TAAA_o2", S3, True),
             ("composite cells (meta-atom level)", "B212_CCCC_o2", S4, True)]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    style(ax, "temperature change from the design point (K)", "frequency change at 145.2 MHz (kHz)")
    ax.axhspan(-15, 15, color=GRID, alpha=0.6, lw=0, label="±15 kHz window")
    for lab, cfg, col, trim in picks:
        r = an.working(an.load(cfg), "box")
        d = an.fd(cfg, r)
        tp = an.tau_parts(r)
        pts = d["f_pts"]
        T = np.array(sorted(pts))
        y = np.array([pts[t] / pts[0.0] - 1 for t in T]) + (tp["tau_gap"] + tp["tau_geo"]) * T
        p = np.polyfit(T, y, min(4, len(T) - 1))
        TT = np.linspace(-25, 25, 201)
        lin = 0.0
        if trim:
            p0 = p.copy(); p0[-2] = 0
            lin = minimize_scalar(lambda t1: np.ptp(np.polyval(p0, TT[40:161]) + t1 * TT[40:161]),
                                  bounds=(-1e-4, 1e-4), method="bounded").x - p[-2]
        yy = (np.polyval(p, TT) + lin * TT) * F_DUP / 1e3
        ax.plot(TT, yy - np.polyval(p, 0) * F_DUP / 1e3, lw=2, color=col, label=lab)
        ax.plot(T, (y + lin * T) * F_DUP / 1e3, "o", ms=6, color=col, mec="white", mew=1)
    ax.set_ylim(-120, 60)
    ax.legend(frameon=False, fontsize=8, loc="lower center")
    fig.tight_layout()
    fig.savefig(FIG / "drift_vs_T.png")
    plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    fig_tau()
    fig_q()
    fig_drift()


if __name__ == "__main__":
    main()
