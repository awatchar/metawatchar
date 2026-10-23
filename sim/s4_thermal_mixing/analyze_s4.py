#!/usr/bin/env python3
"""S4 analysis: temperature coefficient (tau_f) and Q of the working mode of every configuration,
participation formula against finite differences, effective-medium comparison, tables, figures.

Thermal model (see README):
  tau_f = tau_eps + tau_gap + tau_geo
    tau_eps = -1/2 sum_i p_e,i tau_eps,i      (participation formula; checked against FD runs)
    tau_gap = +1/2 sum_j p_gap,j (alpha_d,j - alpha_Cu)   (gap width follows the gap dielectric)
    tau_geo = -alpha_Cu                        (uniform expansion of all copper)
  FD: runs at T0 +- dT (permittivities only) give the central difference tau_eps,FD and the
  curvature c2 = (f+ + f- - 2 f0) / (f0 dT^2) (1/K^2) of the permittivity part.
Q: Q_c from Palace (copper input evaluated at rs_f) rescaled to the copper Rs at the mode frequency
(Q_c_cu, as in S3), Q_d per material from the participations, Q = (1/Q_c_cu + 1/Q_d,al + 1/Q_d,ti)^-1.

Usage: python analyze_s4.py   (local; numpy + matplotlib) -> results/s4_*.csv, s4_summary.json,
       figures/*.png
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from s4_geometry import ALPHA_CU, MAT  # noqa: E402

RES = HERE / "results"
F_DUP = 145.2e6                       # duplexer channel (for kHz numbers)
LC_BAND = (100e6, 260e6)


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return x


def load(run):
    p = RES / f"{run}_eig.csv"
    if not p.exists():
        return None
    with open(p) as fh:
        rows = [{k: num(v) for k, v in r.items()} for r in csv.DictReader(fh)]
    meta = json.loads((HERE / "runs" / f"{run}.json").read_text())
    lc = [r for r in rows if LC_BAND[0] < r["f_MHz"] * 1e6 < LC_BAND[1]]
    for r in lc:
        r["_meta"] = meta
    return lc


def cell_specs(meta):
    return meta["cells"]


def pe_cells(r, n):
    """Electric energy share of each cell's gap (dielectric + gap air)."""
    out = []
    for j in range(n):
        out.append(sum(r.get(f"pe_c{j}_{m}", 0.0) or 0.0 for m in ("al", "ti", "gap")))
    return out


def tau_parts(r, mat=None, alpha_cu=None):
    """Participation-formula tau_f parts (1/K); `mat`, `alpha_cu` override the material data."""
    mat = mat or MAT
    acu = ALPHA_CU if alpha_cu is None else alpha_cu
    cells = cell_specs(r["_meta"])
    pc = pe_cells(r, len(cells))
    t_eps = 0.0
    for j in range(len(cells)):
        for m in ("al", "ti"):
            t_eps += -0.5 * (r.get(f"pe_c{j}_{m}", 0.0) or 0.0) * mat[m]["tau_eps"]
    t_gap = 0.0
    for j, sp in enumerate(cells):
        a_d = mat["ti"]["alpha"] if sp[0] in "TW" else mat["al"]["alpha"]
        t_gap += 0.5 * pc[j] * (a_d - acu)
    return dict(tau_eps=t_eps, tau_gap=t_gap, tau_geo=-acu, tau_pf=t_eps + t_gap - acu)


def shares(r):
    cells = cell_specs(r["_meta"])
    p_al = sum(r.get(f"pe_c{j}_al", 0.0) or 0.0 for j in range(len(cells)))
    p_ti = sum(r.get(f"pe_c{j}_ti", 0.0) or 0.0 for j in range(len(cells)))
    pc = pe_cells(r, len(cells))
    w_t = sum(p for p, sp in zip(pc, cells) if sp[0] in "TW") / sum(pc)
    return dict(p_al=p_al, p_ti=p_ti, x_ti=p_ti / (p_al + p_ti), w_T=w_t, p_cells=pc,
                p_gap_total=sum(pc), cell_nonuni=(max(pc) - min(pc)) / (sum(pc) / len(pc)))


def q_parts(r, tand_al=None, tand_ti=None):
    meta = r["_meta"]
    ta = MAT["al"]["tand"] if tand_al is None else tand_al
    tt = MAT["ti"]["tand"] if tand_ti is None else tand_ti
    sh = shares(r)
    x = r["f_MHz"] * 1e6 / meta["rs_f_Hz"]
    q_c_cu = r["Q_c"] * (2 * x * x / (1 + x * x)) / math.sqrt(x)
    q_dal = 1 / (sh["p_al"] * ta) if sh["p_al"] > 0 else math.inf
    q_dti = 1 / (sh["p_ti"] * tt) if sh["p_ti"] > 0 else math.inf
    q = 1 / (1 / q_c_cu + 1 / q_dal + 1 / q_dti)
    return dict(Q_palace=r["Q"], Q_c=r["Q_c"], Q_c_cu=q_c_cu, Q_d_al=q_dal, Q_d_ti=q_dti, Q=q)


def match(r0, cands):
    n = len(cell_specs(r0["_meta"]))
    p0 = np.array(pe_cells(r0, n))
    best, bc = None, 1e9
    for c in cands:
        cost = abs(c["f_MHz"] / r0["f_MHz"] - 1) / 0.01 + float(np.abs(np.array(pe_cells(c, n)) - p0).sum())
        if cost < bc:
            best, bc = c, cost
    return best


def fd(base, r0, dts=(5.0, -5.0, 12.5, -12.5, 25.0, -25.0)):
    """Permittivity part from the FD runs base_dT<+-x>: central difference at the smallest
    symmetric step, curvature c2 from the +-25 K pair, and a polynomial (degree <= 4) through all
    temperature points."""
    pts = {0.0: r0["f_MHz"]}
    for dT in dts:
        tag = f"_dT{'p' if dT > 0 else 'm'}{abs(dT):g}".replace(".", "p")
        lc = load(base + tag)
        if lc:
            pts[dT] = match(r0, lc)["f_MHz"]
    steps = sorted({abs(t) for t in pts if t and -t in pts})
    if not steps:
        return None
    f0 = pts[0.0]
    h = steps[0]
    out = dict(tau_eps_fd=(pts[h] - pts[-h]) / (2 * h * f0), fd_step=h, f_pts=pts)
    H = steps[-1]
    out["c2"] = (pts[H] + pts[-H] - 2 * f0) / (f0 * H * H)
    if len(steps) > 1:
        out["tau_eps_fd_25"] = (pts[H] - pts[-H]) / (2 * H * f0)
    T = np.array(sorted(pts))
    y = np.array([pts[t] / f0 - 1 for t in T])
    out["poly"] = np.polyfit(T, y, min(4, len(T) - 1)).tolist()
    return out


def working(lc, kind):
    if not lc:
        return None
    if kind == "line":
        hi = [r for r in lc if r["f_MHz"] > 160]
        return max(hi, key=lambda r: r["f_MHz"]) if hi else None
    return min(lc, key=lambda r: r["f_MHz"])


def drift_pp(tau1, c2, W):
    """Peak-to-peak frequency excursion (relative) of 1 + tau1 dT + c2 dT^2 over dT in [-W, W]."""
    T = np.linspace(-W, W, 2001)
    y = tau1 * T + c2 * T * T
    return float(y.max() - y.min())


def write_csv(path, rows):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.8g}" if isinstance(v, float) else v) for k, v in r.items()})


def summarize(base, kind, label="", extra=None):
    """One row for the working mode of a configuration (nominal run `base`)."""
    lc = load(base)
    r = working(lc, kind)
    if r is None:
        return None
    meta = r["_meta"]
    row = dict(config=base, label=label, cells="/".join(meta["cells"]), order=meta["order"],
               mode=int(r["mode"]), f_MHz=r["f_MHz"])
    sh = shares(r)
    row.update({k: v for k, v in sh.items() if k != "p_cells"})
    row["p_cells"] = " ".join(f"{p:.4f}" for p in sh["p_cells"])
    tp = tau_parts(r)
    row.update({k: v * 1e6 for k, v in tp.items()})           # ppm/K
    d = fd(base, r)
    if d:
        row["tau_eps_fd"] = d["tau_eps_fd"] * 1e6
        row["fd_minus_pf_ppm"] = (d["tau_eps_fd"] - tp["tau_eps"]) * 1e6
        row["c2_ppb_K2"] = d["c2"] * 1e9
        tau1 = d["tau_eps_fd"] + tp["tau_gap"] + tp["tau_geo"]
        row["tau_f_fd"] = tau1 * 1e6
        row["fd_step_K"] = d["fd_step"]
        row["n_T_points"] = len(d["f_pts"])
        lin = tp["tau_gap"] + tp["tau_geo"]
        poly = np.array(d["poly"])
        for W in (10, 20, 25):
            T = np.linspace(-W, W, 2001)
            y = np.polyval(poly, T) + lin * T
            row[f"pp_kHz_pm{W}K_at145"] = float(y.max() - y.min()) * F_DUP / 1e3
            row[f"pp_kHz_pm{W}K_tau1zero"] = drift_pp(0.0, d["c2"], W) * F_DUP / 1e3
    row.update(q_parts(r))
    row["Q_c_cu_at145"] = row["Q_c_cu"] * math.sqrt(145e6 / (r["f_MHz"] * 1e6))
    if extra:
        row.update(extra)
    return row


BOX = [  # (config, k = number of TiO2 cells, label)
    ("B212_AAAA_o2", 0, "all alumina"),
    ("B212_TAAA_o2", 1, "1 TiO2, bottom row"), ("B212_ATAA_o2", 1, "1 TiO2, top row"),
    ("B212_TATA_o2", 2, "2 TiO2, bottom row"), ("B212_ATAT_o2", 2, "2 TiO2, top row"),
    ("B212_TTAA_o2", 2, "2 TiO2, one column"), ("B212_TAAT_o2", 2, "2 TiO2, diagonal"),
    ("B212_ATTT_o2", 3, "3 TiO2, alumina bottom"), ("B212_TATT_o2", 3, "3 TiO2, alumina top"),
    ("B212_TTTT_o2", 4, "all TiO2"),
]
DETUNED = [("B212_TAAA_o2", 1.00), ("B212_TAAA_W3p41_o2", 1.018), ("B212_TAAA_W3p43_o2", 1.020),
           ("B212_TAAA_d030_o2", 1.03), ("B212_TAAA_d060_o2", 1.06), ("B212_TAAA_d100_o2", 1.10)]
LINES = [("L4_AAAA_o2", "all alumina"), ("L4_TAAA_o2", "1 TiO2 at the end"),
         ("L4_ATAA_o2", "1 TiO2 inside")]
CELLS = ["P1_A_o2_pmc", "P1_A_o2_pec", "P1_W3p29_o2_pmc", "P1_W3p29_o2_pec", "P1_W3p0_o2_pmc",
         "P1_W3p6_o2_pmc", "P1_T3p0_o2_pmc", "P1_T3p6_o2_pmc", "P1_T4p2_o2_pmc", "P1_T3p6_o2_pec",
         "P1_C0p6_o2_pmc", "P1_C0p75_o2_pmc", "P1_C0p75_o2_pec", "P1_C1p0_o2_pmc_t",
         "P1_A_d1p01_o2_pmc", "B111_A_o2", "B111_W3p29_o2"]


def cell_row(run):
    lc = load(run)
    if not lc:
        return None
    r = min(lc, key=lambda r: r["f_MHz"])
    sh, tp, q = shares(r), tau_parts(r), q_parts(r)
    row = dict(run=run, cell=r["_meta"]["cells"][0], f_MHz=r["f_MHz"], **q,
               Q_c_cu_at_150=q["Q_c_cu"] * math.sqrt(149.92 / r["f_MHz"]),
               p_al=sh["p_al"], p_ti=sh["p_ti"], p_gap=sh["p_gap_total"], x_ti=sh["x_ti"],
               **{k: v * 1e6 for k, v in tp.items()})
    d = fd(run, r)
    if d:
        row.update(tau_eps_fd=d["tau_eps_fd"] * 1e6, c2_ppb_K2=d["c2"] * 1e9)
    return row


def all_modes(base):
    """tau_f of every LC mode of a box configuration (participation formula + FD where present)."""
    lc = load(base)
    out = []
    for r in lc or []:
        if r["f_MHz"] > 250:
            continue
        tp, sh, q = tau_parts(r), shares(r), q_parts(r)
        row = dict(config=base, mode=int(r["mode"]), f_MHz=r["f_MHz"], Q=q["Q"], x_ti=sh["x_ti"],
                   w_T=sh["w_T"], tau_pf=tp["tau_pf"] * 1e6, tau_gap=tp["tau_gap"] * 1e6,
                   p_cells=" ".join(f"{p:.4f}" for p in sh["p_cells"]))
        d = fd(base, r)
        if d:
            row["tau_f_fd"] = (d["tau_eps_fd"] + tp["tau_gap"] + tp["tau_geo"]) * 1e6
            row["fd_step_K"] = d["fd_step"]
            row["c2_ppb_K2"] = d["c2"] * 1e9
        out.append(row)
    return out


def tau_best(row):
    return row.get("tau_f_fd", row["tau_pf"])


def zero_cross(xs, ys):
    xs, ys = np.asarray(xs), np.asarray(ys)
    o = np.argsort(xs)
    xs, ys = xs[o], ys[o]
    for i in range(len(xs) - 1):
        if ys[i] * ys[i + 1] <= 0:
            return float(xs[i] - ys[i] * (xs[i + 1] - xs[i]) / (ys[i + 1] - ys[i]))
    return float("nan")


def em_compare():
    """Effective-medium (Lorentz, from unit-cell eigenmodes) against full wave: lines and box.
    Compared quantity: the permittivity part tau_eps (the expansion terms are common)."""
    import s4_effective_medium as em
    rows = []
    cp = {}
    for key, run in (("A", "P1_A_o2"), ("T", "P1_W3p29_o2")):
        lm, lp = load(run + "_pmc"), load(run + "_pec")
        if not (lm and lp):
            return rows
        rm, rp = min(lm, key=lambda r: r["f_MHz"]), min(lp, key=lambda r: r["f_MHz"])
        cp[key] = dict(p=em.beig_params(rm["f_MHz"] * 1e6, rp["f_MHz"] * 1e6),
                       tau_m=tau_parts(rm)["tau_eps"], tau_mp=tau_parts(rp)["tau_eps"])

    def at(arr, dT):
        return [em.shift(cp[c]["p"], cp[c]["tau_m"], cp[c]["tau_mp"], dT) for c in arr]
    for base, _ in LINES:
        arr = base.split("_")[1]
        fw = summarize(base, "line")
        if not fw:
            continue
        for model in ("layer", "mix"):
            f0 = max(em.line_modes(at(arr, 0.0), model=model))
            fp = max(em.line_modes(at(arr, 25.0), model=model))
            fm = max(em.line_modes(at(arr, -25.0), model=model))
            rows.append(dict(config=base, model=model, f_fw_MHz=fw["f_MHz"], f_em_MHz=f0 / 1e6,
                             df_pct=(f0 / 1e6 / fw["f_MHz"] - 1) * 100,
                             tau_eps_fw=fw.get("tau_eps_fd", fw["tau_eps"]),
                             tau_eps_em=(fp - fm) / (50 * f0) * 1e6))
    for base, k, _ in BOX:
        arr = base.split("_")[1]
        modes = sorted(all_modes(base), key=lambda r: r["f_MHz"])
        if not modes:
            continue
        f0 = em.mixture_fmp(at(arr, 0.0))
        tau = (em.mixture_fmp(at(arr, 25.0)) - em.mixture_fmp(at(arr, -25.0))) / (50 * f0)
        work, top = modes[0], modes[-1]
        rows.append(dict(config=base, model="mix, mu_eff = 0", f_em_MHz=f0 / 1e6,
                         f_fw_MHz=work["f_MHz"], f_fw_inphase_MHz=top["f_MHz"],
                         tau_eps_em=tau * 1e6,
                         tau_eps_fw=tau_best(work) - work["tau_gap"] + ALPHA_CU * 1e6,
                         tau_eps_fw_inphase=tau_best(top) - top["tau_gap"] + ALPHA_CU * 1e6))
    return rows


def sensitivity():
    """TiO2 share x_ti for tau_f = 0 when the material data change (participations fixed)."""
    out = []
    cases = [("baseline", {}, None)]
    for v in (-700e-6, -800e-6, -1000e-6, -1100e-6):
        cases.append((f"tau_eps_TiO2 = {v * 1e6:+.0f} ppm/K", {"ti": {"tau_eps": v}}, None))
    for v in (80e-6, 130e-6):
        cases.append((f"tau_eps_Al2O3 = {v * 1e6:+.0f} ppm/K", {"al": {"tau_eps": v}}, None))
    for v in (16.5e-6, 17.7e-6):
        cases.append((f"alpha_Cu = {v * 1e6:.1f} ppm/K", {}, v))
    for v in (5e-6, 9e-6):
        cases.append((f"alpha_TiO2 = {v * 1e6:.0f} ppm/K", {"ti": {"alpha": v}}, None))
    for name, ov, acu in cases:
        mat = {k: dict(MAT[k], **ov.get(k, {})) for k in MAT}
        xs, ys = [], []
        for base, k, _ in BOX:
            r = working(load(base), "box")
            if r is None:
                continue
            xs.append(shares(r)["x_ti"])
            ys.append(tau_parts(r, mat, acu)["tau_pf"] * 1e6)
        xd, yd, xdx = [], [], []
        for base, rel in DETUNED:
            r = working(load(base), "box")
            if r is None:
                continue
            xd.append(rel)
            xdx.append(shares(r)["x_ti"])
            yd.append(tau_parts(r, mat, acu)["tau_pf"] * 1e6)
        cfit = np.polyfit(xs, ys, 1) if len(xs) > 1 else [np.nan, np.nan]
        out.append(dict(case=name, x_ti_zero_count_fit=-cfit[1] / cfit[0], slope_ppm_per_share=cfit[0],
                        tau_A=ys[0] if ys else np.nan, tau_allTiO2=ys[-1] if ys else np.nan,
                        detune_rel_zero=zero_cross(xd, yd) if len(xd) > 1 else np.nan,
                        x_ti_zero_detuned=zero_cross(xdx, yd) if len(xd) > 1 else np.nan))
    return out


def main():
    cells = [r for r in (cell_row(c) for c in CELLS) if r]
    write_csv(RES / "s4_cells.csv", cells)
    box = [r for r in (summarize(b, "box", lab, dict(k=k, count_share=k / 4)) for b, k, lab in BOX) if r]
    det = [r for r in (summarize(b, "box", f"k = 1, f_T/f_A = {rel:.2f}", dict(k=1, detune=rel))
                       for b, rel in DETUNED) if r]
    comp = [r for r in (summarize("B212_CCCC_o2", "box", "composite cells (reference)"),) if r]
    line = [r for r in (summarize(b, "line", lab) for b, lab in LINES) if r]
    lev = [r for r in (summarize(b, "box", lab) for b, lab in
                       (("B111_A_o2", "1 cell, tight box"), ("B212_AAAA_o2", "2x1x2, tight box"),
                        ("B141_AAAA_o2", "1x4x1 coaxial stack, tight box"))) if r]
    write_csv(RES / "s4_box_mixing.csv", box)
    write_csv(RES / "s4_box_detuned.csv", det + comp)
    write_csv(RES / "s4_line.csv", line)
    modes = [m for b, _, _ in BOX for m in all_modes(b)] + [m for b, _ in DETUNED[1:] for m in all_modes(b)] \
        + all_modes("B212_CCCC_o2")
    write_csv(RES / "s4_box_all_modes.csv", modes)
    emr = em_compare()
    write_csv(RES / "s4_effective_medium.csv", emr)
    sens = sensitivity()
    write_csv(RES / "s4_sensitivity.csv", sens)
    qs = []
    for row in box + det + comp:
        r = working(load(row["config"]), "box")
        for tt in (6e-5, 1e-4, 3e-4, 1e-3):
            q = q_parts(r, tand_ti=tt)
            qs.append(dict(config=row["config"], x_ti=row["x_ti"], tand_ti=tt, Q=q["Q"],
                           Q_d_ti=q["Q_d_ti"], Q_c_cu=q["Q_c_cu"], Q_d_al=q["Q_d_al"]))
    write_csv(RES / "s4_q_tand.csv", qs)
    qcell = next((c for c in cells if c["run"] == "P1_A_o2_pmc"), None)
    lev5 = []
    for row in lev:
        lev5.append(dict(config=row["config"], label=row["label"], f_MHz=row["f_MHz"], Q=row["Q"],
                         Q_c_cu=row["Q_c_cu"], Q_over_Qcell=row["Q"] / qcell["Q"],
                         Qc_over_Qccell_fscaled=row["Q_c_cu"] / (qcell["Q_c_cu"] * math.sqrt(row["f_MHz"] / qcell["f_MHz"])),
                         p_cells=row["p_cells"], cell_nonuni=row["cell_nonuni"]))
    write_csv(RES / "s4_lever5.csv", lev5)
    summ = dict(box=box, detuned=det, composite=comp, line=line, em=emr, sensitivity=sens,
                lever5=lev5, cells=cells)
    (RES / "s4_summary.json").write_bytes((json.dumps(summ, indent=1, default=float) + "\n").encode())
    for t, rows in (("cells", cells), ("box", box), ("detuned+composite", det + comp), ("line", line),
                    ("em", emr), ("sensitivity", sens), ("lever5", lev5)):
        print(f"--- {t}")
        for r in rows:
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                   if not k.startswith("pp_") and k not in ("mode", "order")})
    if (HERE / "plot_s4.py").exists():
        import plot_s4
        plot_s4.main()


if __name__ == "__main__":
    main()
