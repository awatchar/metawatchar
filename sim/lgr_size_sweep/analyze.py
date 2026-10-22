#!/usr/bin/env python3
"""Local analysis of the summarized LGR-A results (results/*.csv pulled from the cluster).

  python analyze.py tune results/r1.csv [results/r2.csv ...]
      For each (V, a, d, order, h) group, fit 1/f0^2 = A + B*s through the two runs closest to
      the 145.5 MHz target and print the lip height s* that hits it (capacitance is linear in s
      for parallel plates plus a fixed fringing part; L is nearly independent of s).

  python analyze.py analytic
      Print the analytic conductor Q of the coaxial-return LGR model for the nine sweep cases.

Analytic model (shielded LGR, coaxial return flux, no end caps, REALISM = 1):
  rho = r^2 / (R^2 - (r+t)^2)
  Q_c = omega mu0 pi r^2 (1+rho) / (2 pi Rs [r + rho^2 ((r+t) + R)])
The ratio Q_c(Palace) / Q_c(analytic) is the empirical derating ("REALISM") of that model
for this closed-can geometry with end clearance.
"""
import csv
import math
import re
import sys
from collections import defaultdict

F_T = 145.5
MU0 = 4e-7 * math.pi
SIGMA = 5.8e7


def qc_analytic(R, r, t, f=F_T * 1e6):
    rs = math.sqrt(math.pi * f * MU0 / SIGMA)
    w = 2 * math.pi * f
    rho = r * r / (R * R - (r + t) ** 2)
    return w * MU0 * math.pi * r * r * (1 + rho) / (2 * math.pi * rs * (r + rho * rho * ((r + t) + R)))


def load(paths):
    rows = []
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh):
                for k, v in r.items():
                    try:
                        r[k] = float(v)
                    except (TypeError, ValueError):
                        pass
                r["_src"] = p
                rows.append(r)
    return rows


def tune(paths):
    groups = defaultdict(list)
    for r in load(paths):
        groups[(r["V_l"], r["a"], r["d_mm"], r["order"], r["h_scale"])].append(r)
    print("| V (L) | a | d | order | h | runs (s: f0) | s* (mm) | s_max (mm) |")
    print("|---|---|---|---|---|---|---|---|")
    for key in sorted(groups):
        rs = sorted(groups[key], key=lambda r: r["s_mm"])
        # unique s values, keep the latest run for each s
        by_s = {}
        for r in rs:
            by_s[r["s_mm"]] = r
        pts = sorted(by_s.values(), key=lambda r: abs(r["f0_MHz"] - F_T))[:2]
        s_star = float("nan")
        if len(pts) == 2 and pts[0]["s_mm"] != pts[1]["s_mm"]:
            (s1, f1), (s2, f2) = [(p["s_mm"], p["f0_MHz"]) for p in pts]
            y1, y2 = 1 / f1 ** 2, 1 / f2 ** 2
            b = (y2 - y1) / (s2 - s1)
            s_star = s1 + (1 / F_T ** 2 - y1) / b
        s_max = rs[0]["R_mm"] - 3.0 - rs[0]["r_mm"] - 3.0
        runs = ", ".join(f"{p['s_mm']:.2f}: {p['f0_MHz']:.2f}" for p in sorted(by_s.values(), key=lambda r: r["s_mm"]))
        print(f"| {key[0]} | {key[1]} | {key[2]} | {key[3]:.0f} | {key[4]} | {runs} | {s_star:.3f} | {s_max:.2f} |")


def analytic():
    for V in (0.15, 0.25, 0.5):
        R = (V * 1e-3 / (2 * math.pi)) ** (1 / 3)
        for a in (0.4, 0.5, 0.6):
            print(f"V={V} a={a} R={R*1e3:.2f} mm Q_c,analytic={qc_analytic(R, a * R, 0.003):.0f}")


EDGE = [("air_slab_face", "air, slab outer face"), ("air_bore_face", "air, slab bore face"),
        ("air_triple_tip", "air, lip-tip triple line"),
        ("air_lip_corner", "air, lip outer corner"), ("diel_triple_tip", "alumina, lip-tip triple line"),
        ("air_triple_bore", "air, bore triple line"), ("air_end_face", "air, slab end face"),
        ("air_end_tip_corner", "air, tip/end corner")]


def fit_s(rows):
    """Linear 1/f^2(s) through the two runs nearest the target. Returns (s*, s0), with s0 the
    intercept of 1/f^2 = B (s + s0)."""
    by_s = {r["s_mm"]: r for r in rows}
    pts = sorted(by_s.values(), key=lambda r: abs(r["f0_MHz"] - F_T))[:2]
    if len(pts) < 2:
        return float("nan"), float("nan")
    (s1, f1), (s2, f2) = [(p["s_mm"], p["f0_MHz"]) for p in pts]
    b = (1 / f2 ** 2 - 1 / f1 ** 2) / (s2 - s1)
    s_star = s1 + (1 / F_T ** 2 - 1 / f1 ** 2) / b
    s0 = 1 / (f1 ** 2 * b) - s1
    return s_star, s0


REF_BASE = "v025_a05_d0635_s1092"        # case used for the convergence study
PROD = r"v\d{3}_a\d{2}_d\d{4}_s\d{4}_o3_v2"   # production runs (v2 mesh, FE order 3)
V1 = r"v\d{3}_a\d{2}_d\d{4}_s\d{4}_o3"        # first production mesh (v1), kept for comparison
VARIANTS = [   # suffix after REF_BASE -> FE order, description (mesh recipes in lgr_geometry.py)
    ("", 2, "v1 mesh"),
    ("_h070", 2, "v1, all sizes x0.7"),
    ("_o3", 3, "v1 mesh"),
    ("_h070_o3", 3, "v1, all sizes x0.7"),
    ("_mG070_o3", 3, "v1, gap size 0.45 -> 0.32 mm"),
    ("_mE050_o3", 3, "v1, gap-edge size 0.25 -> 0.125 mm"),
    ("_mB060_o3", 3, "v1, bulk size x0.6"),
    ("_mL050_o3", 3, "v1, loop/lip surface size 2 -> 1 mm"),
    ("_mN1_o2", 2, "N1: denser size-field sampling + loop/lip edges 0.4 mm"),
    ("_mN1_o3", 3, "N1"),
    ("_mN1_o4", 4, "N1"),
    ("_mN2_o3", 3, "N1, loop/lip surface size 2 -> 1 mm"),
    ("_mN4_o2", 2, "N4 = N1 with loop/lip edges 0.2 mm (= v2 recipe)"),
    ("_mN4_o3", 3, "N4 (separate run of the v2 recipe)"),
    ("_o3_v2", 3, "v2 production mesh"),
    ("_mN4_o4", 4, "N4 = v2 mesh"),
    ("_mN5_o2", 2, "N5 = N1 with loop/lip edges 0.1 mm"),
]


def air_edge_max(r):
    return max((r.get(f"E_{k}_kV_mm", 0.0), txt) for k, txt in EDGE if k.startswith("air"))


HDR = ("| V (L) | a | R (mm) | r (mm) | z (mm) | s (mm) | clearance (mm) | f0 (MHz) | Q_c | Q_d | Q | p_d | "
       "E gap (kV/mm) | V gap (V) | max air E near edges (kV/mm) | Q_c analytic | Palace/analytic |\n"
       "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")


def tuned_row(r, P, k):
    qa = qc_analytic(r["R_mm"] * 1e-3, r["r_mm"] * 1e-3, 3e-3)
    qc = r["Q_c"] * k
    q = 1 / (1 / qc + 1 / r["Q_d"])
    emax, lab = air_edge_max(r)
    P(f"| {r['V_l']:g} | {r['a']:g} | {r['R_mm']:.2f} | {r['r_mm']:.2f} | {r['z_mm']:.2f} | "
      f"{r['s_mm']:.2f} | {r['clearance_mm']:.2f} | {r['f0_MHz']:.2f} | {qc:,.0f} | "
      f"{r['Q_d']:,.0f} | **{q:,.0f}** | {r['p_d']:.3f} | {r['E_gap_wall_mid_kV_mm']:.2f} | "
      f"{r['V_gap_V']:.0f} | {emax:.2f} ({lab}) | {qa:,.0f} | {qc / qa:.3f} |")


def report(paths):
    byname = {}
    for r in load(paths):          # the latest row per case name wins
        byname[r["name"]] = r
    rows = list(byname.values())
    key = lambda r: (r["d_mm"], r["V_l"], r["a"])
    prod = sorted([r for r in rows if re.fullmatch(PROD, r["name"])], key=key)
    final = [r for r in prod if r["clearance_mm"] >= 2.99]
    relaxed = [r for r in prod if r["clearance_mm"] < 2.99]
    amr = sorted([r for r in rows if r["name"].startswith(REF_BASE + "_mN1_amr_o3")],
                 key=lambda r: r["dofs"] if r["dofs"] != "" else 0)
    ref_prod = byname.get(REF_BASE + "_o3_v2")
    k = 1.0   # Q_c values are reported raw (v2 mesh); see the convergence section for the error estimate
    out = []
    P = out.append
    P("# LGR-A size sweep: results")
    P("")
    P("Palace 0.18.1 eigenmode solutions. Copper is an impedance boundary (Rs = 3.147 mOhm/sq, "
      "Ls = Rs/omega, sigma = 5.8e7 S/m at 145.5 MHz) on the shield, loop and lips; the gap is alumina "
      "(eps_r 9.8, tan_delta 4e-5). Q_c = conductor-only Q (shield + loop + lips), Q_d = 1/(p_d tan_delta) "
      "with p_d the electric-energy participation of the alumina, Q = total = (1/Q_c + 1/Q_d)^-1. Fields are "
      "peak phasor amplitudes for a stored energy of **10 uJ**; V gap = E gap x d. Edge probes sit "
      "0.1 mm from a metal/dielectric/air corner: the field at a sharp metal edge is singular, so these values "
      "depend on the probe distance and, in hardware, on the edge radius. Analytic Q_c: coaxial-return LGR "
      "model without end caps (REALISM = 1, see `analyze.py`), same R, r, t; ratio = Palace / analytic.")
    P("")
    if k != 1.0:
        P(f"**Mesh correction.** Q_c in the tuned-case tables is the v2-mesh value times **{k:.4f}**, the ratio "
          f"of the adaptive-mesh-refinement (AMR) reference to the v2 mesh for the V = 0.25 L, a = 0.5 case (see "
          f"\"Convergence\"). Raw v2 values are in `results/v2_*.csv`. Q is recomputed from the corrected Q_c.")
        P("")
    for d in sorted({r["d_mm"] for r in final}):
        P(f"## Tuned cases, d = {d:g} mm (FE order 3, clearance >= 3 mm)")
        P("")
        P(HDR)
        for r in final:
            if r["d_mm"] == d:
                tuned_row(r, P, k)
        P("")
    if relaxed:
        P("## Near-miss cases tuned with the clearance relaxed below 3 mm (violate the 3 mm rule)")
        P("")
        P(HDR)
        for r in relaxed:
            tuned_row(r, P, k)
        P("")
    pec = [r for r in rows if r["name"].endswith("_pec") and re.fullmatch(PROD, r["name"][:-4])]
    if pec:
        P("## Loss split: lossless (PEC) shield (raw v2-mesh values)")
        P("")
        P("| case | Q_c, all copper lossy | Q_c, PEC shield (loop + lips only) | implied Q_c of the shield alone | "
          "shield share of conductor loss |")
        P("|---|---|---|---|---|")
        for r in pec:
            base = byname.get(r["name"][:-4])
            if base:
                q_sh = 1 / (1 / base["Q_c"] - 1 / r["Q_c"])
                P(f"| {base['name']} | {base['Q_c']:,.0f} | {r['Q_c']:,.0f} | {q_sh:,.0f} | "
                  f"{100 * base['Q_c'] / q_sh:.0f} % |")
        P("")
    P("## Field at every probe (kV/mm at 10 uJ)")
    P("")
    P("| case | gap (wall) | gap (lips) | " + " | ".join(t for _, t in EDGE) + " |")
    P("|---|---|---|" + "---|" * len(EDGE))
    for r in prod:
        P(f"| {r['name']} | {r['E_gap_wall_mid_kV_mm']:.2f} | {r['E_gap_lip_mid_kV_mm']:.2f} | " +
          " | ".join(f"{r.get(f'E_{k_}_kV_mm', float('nan')):.2f}" for k_, _ in EDGE) + " |")
    P("")
    groups = defaultdict(list)
    for r in rows:
        if r["order"] == 2 and r["h_scale"] == 1.0 and re.fullmatch(r"v\d{3}_a\d{2}_d\d{4}_s\d{4}", r["name"]):
            groups[(r["d_mm"], r["V_l"], r["a"])].append(r)
    if "test_o2" in byname:
        groups[(0.635, 0.25, 0.5)].append(byname["test_o2"])
    P("## Tuning runs (FE order 2, v1 mesh) and reachability of 145.5 MHz")
    P("")
    P("s_max = R - 3 mm - (r + t), i.e. lips ending 3 mm from the shield wall. s* comes from a linear fit of "
      "1/f0^2 against s through the two runs nearest the target; s0 is the fit intercept "
      "(1/f0^2 = B (s + s0)), close to the wall depth t = 3 mm plus a small fringing allowance. Q here is "
      "the raw order-2 v1-mesh value (low, see \"Convergence\").")
    P("")
    P("| d (mm) | V (L) | a | runs, s (mm): f0 (MHz), Q | s_max (mm) | s* (mm) | s0 (mm) | 145.5 MHz reachable? |")
    P("|---|---|---|---|---|---|---|---|")
    for gk in sorted(groups):
        g = sorted(groups[gk], key=lambda r: r["s_mm"])
        smax = g[0]["R_mm"] - 6.0 - g[0]["r_mm"]
        s_star, s0 = fit_s(g)
        if s_star <= smax + 0.01:
            verdict = "yes"
        else:
            need = g[0]["R_mm"] - g[0]["r_mm"] - 3.0 - s_star
            verdict = (f"**no**: f0 = {min(r['f0_MHz'] for r in g):.1f} MHz at s_max; would need a clearance of "
                       f"{need:.1f} mm" + (" (lips beyond the shield)" if need < 0 else ""))
        runs = "; ".join(f"{r['s_mm']:.2f}: {r['f0_MHz']:.1f}, {r['Q']:,.0f}" for r in g)
        P(f"| {gk[0]:g} | {gk[1]:g} | {gk[2]:g} | {runs} | {smax:.2f} | {s_star:.2f} | {s0:.2f} | {verdict} |")
    P("")
    P("## Convergence (V = 0.25 L, a = 0.5, d = 0.635 mm, s = 10.92 mm)")
    P("")
    P("Mesh recipes (`lgr_geometry.py`): v1 = 0.45 mm in and around the alumina, 0.25 mm along its edges, "
      "2 mm on loop/lip surfaces, R/5 in the bulk, no extra refinement of the loop rims and lip edges; "
      "v2 = v1 with a denser size-field sampling and 0.2 mm along every loop/lip edge (the magnetic field is "
      "singular at the square tube rims and lip corners, and that dominates the conductor-loss error). "
      "N5 = loop/lip edges 0.1 mm. AMR = Palace adaptive refinement (Doerfler fraction 0.5, 3 iterations) "
      "from N1: its error indicator targets the electric field, refines only ~1 % of the mesh per iteration, "
      "and leaves Q_c unchanged, so it does not help here.")
    P("")
    ref = byname.get(REF_BASE + "_mN4_o4") or ref_prod or byname.get(REF_BASE + "_mN4_o3")
    P(f"Reference: `{ref['name']}`.")
    P("")
    P("| run | FE order | mesh | tets | DOFs | f0 (MHz) | Q_c | Q | E gap (kV/mm) | Q_c vs reference |")
    P("|---|---|---|---|---|---|---|---|---|---|")
    for suf, order, desc in VARIANTS:
        r = byname.get(REF_BASE + suf)
        if r is None:
            continue
        dofs = f"{r['dofs']:,.0f}" if r["dofs"] != "" else ""
        P(f"| `{r['name']}` | {order} | {desc} | {r['n_tet']:,.0f} | {dofs} | {r['f0_MHz']:.3f} | "
          f"{r['Q_c']:,.0f} | {r['Q']:,.0f} | {r['E_gap_wall_mid_kV_mm']:.3f} | "
          f"{100 * (r['Q_c'] / ref['Q_c'] - 1):+.1f} % |")
    for r in amr:
        it = r["name"].split("@")[1] if "@" in r["name"] else "final"
        P(f"| `{r['name']}` | 3 | AMR, {it} | | {r['dofs']:,.0f} | {r['f0_MHz']:.3f} | {r['Q_c']:,.0f} | "
          f"{r['Q']:,.0f} | {r['E_gap_wall_mid_kV_mm']:.3f} | {100 * (r['Q_c'] / ref['Q_c'] - 1):+.1f} % |")
    P("")
    v1 = {r["name"]: r for r in rows if re.fullmatch(V1, r["name"])}
    pairs = [(v1[r["name"][:-3]], r) for r in prod if r["name"][:-3] in v1]
    if pairs:
        P("## v1 vs v2 mesh, all production cases (FE order 3, raw Q_c)")
        P("")
        P("| case | Q_c v1 | Q_c v2 | v2/v1 | f0 v1 (MHz) | f0 v2 (MHz) |")
        P("|---|---|---|---|---|---|")
        for a_, b_ in pairs:
            P(f"| {a_['name']} | {a_['Q_c']:,.0f} | {b_['Q_c']:,.0f} | {b_['Q_c'] / a_['Q_c']:.3f} | "
              f"{a_['f0_MHz']:.2f} | {b_['f0_MHz']:.2f} |")
        P("")
    print("\n".join(out))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "tune":
        tune(sys.argv[2:])
    elif cmd == "analytic":
        analytic()
    elif cmd == "report":
        sys.stdout.reconfigure(encoding="utf-8")
        report(sys.argv[2:])
