#!/usr/bin/env python3
"""Collect the per-case summaries (results/<round>.json written by postprocess.py in the job) into
results/sizing.json (loop choice) and results/summary.json (all study results).

  python analyze.py sizing      # after the sizing round: fit beta(area), choose the loop
  python analyze.py summary     # after all rounds
"""
import glob
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TARGET_DB = -30.0


def records():
    recs = {}
    for p in sorted(glob.glob(str(HERE / "results" / "*.json"))):
        if Path(p).name in ("sizing.json", "summary.json"):
            continue
        for r in json.loads(Path(p).read_text()):
            if "error" not in r:
                recs[r["case"]] = r
    return recs


def beta_target(db=TARGET_DB):
    s = 10 ** (db / 20)
    return s / (2 * (1 - s))      # symmetric: S21 = 2 beta / (1 + 2 beta)


def fit_power(pts):
    """least squares ln(beta) = ln(k) + n ln(A)"""
    xs = [math.log(a) for a, _ in pts]
    ys = [math.log(b) for _, b in pts]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sl = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    return math.exp(my - sl * mx), sl


def sizing(recs):
    sizes = [(3, 3), (5, 5), (8, 6), (10, 8)]
    pts, rows = [], []
    for w, h in sizes:
        r = recs.get(f"e_xy{w}x{h}_o2")
        if not r:
            continue
        b = 0.5 * (r["beta"]["1"] + r["beta"]["2"])
        pts.append((w * h, b))
        rows.append(dict(loop_mm=[w, h], area_mm2=w * h, f0_Hz=r["f0_Hz"], Q_L=r["Q"], Q_u=r["Q_u"],
                         Q_ext=r["Q_ext"], beta=r["beta"], S21_f0_dB=r["S21_f0_dB_pred"]))
    # local fit through the two sizes that bracket the target
    bt = beta_target()
    k, n = fit_power(pts[:2])
    a_t = (bt / k) ** (1 / n)
    chosen = (4, 4)
    b_ch = k * (chosen[0] * chosen[1]) ** n
    base = recs["e_base_o2"]
    f0_est = {f"xy{w}x{h}": r["f0_Hz"] / 1e9 for (w, h), r in
              ((s, recs.get(f"e_xy{s[0]}x{s[1]}_o2")) for s in sizes) if r}
    f0_est["xy4x4"] = 0.5 * (f0_est["xy3x3"] + f0_est["xy5x5"])
    out = dict(target_S21_dB=TARGET_DB, beta_target=bt, fit=dict(k=k, n=n, sizes="3x3, 5x5"),
               area_for_target_mm2=a_t, side_for_target_mm=math.sqrt(a_t), chosen_loop_mm=list(chosen),
               beta_chosen_pred=b_ch, S21_chosen_pred_dB=20 * math.log10(2 * b_ch / (1 + 2 * b_ch)),
               f0_base_GHz=base["f0_Hz"] / 1e9, f0_est_GHz=f0_est, sizes=rows,
               rz5x5=recs.get("e_rz5x5_o2"))
    (HERE / "results" / "sizing.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k_: v for k_, v in out.items() if k_ not in ("sizes", "rz5x5")}, indent=1))
    for r in rows:
        print(r["loop_mm"], round(r["f0_Hz"] / 1e6, 4), round(r["Q_L"]), round(r["Q_u"]), {k_: round(v) for k_, v in r["Q_ext"].items()},
              {k_: round(v, 5) for k_, v in r["beta"].items()}, round(r["S21_f0_dB"], 2))
    if out["rz5x5"]:
        r = out["rz5x5"]
        print("rz5x5", round(r["f0_Hz"] / 1e6, 4), round(r["Q"]), r.get("Q_ext"), r.get("S21_f0_dB_pred"))


def prom_check():
    """Max relative difference between the reduced-order sweep and the full-order solves."""
    import cmath
    import csv
    R = list(csv.DictReader(open(HERE / "results" / "sparams_d_xy4x4_o2.csv")))
    F = [float(r["f_Hz"]) for r in R]
    full = list(csv.reader(open(HERE / "postpro" / "d_xy4x4_o2_full" / "port-S.csv")))[1:]
    worst = 0.0
    for r in full:
        f = float(r[0]) * 1e9
        i = int((f - F[0]) // (F[1] - F[0]))
        t = (f - F[i]) / (F[i + 1] - F[i])
        a = complex(float(R[i]["S21_re"]), float(R[i]["S21_im"]))
        b = complex(float(R[i + 1]["S21_re"]), float(R[i + 1]["S21_im"]))
        s21p = 1 / ((1 - t) / a + t / b)      # 1/S21 is linear in f near the pole
        s21f = 10 ** (float(r[3]) / 20) * cmath.exp(1j * math.radians(float(r[4])))
        worst = max(worst, abs(s21p - s21f) / abs(s21f))
    return dict(full_order_freqs_GHz=[float(r[0]) for r in full], max_rel_diff_S21=worst)


def summary(recs):
    sz = json.loads((HERE / "results" / "sizing.json").read_text())
    W = H = 4
    mn = f"xy{W}x{H}"
    sel = lambda r, ks: {k: r[k] for k in ks if k in r}
    dk = ["f0_Hz", "S21_f0_dB", "Q_L_3dB", "beta1", "beta2", "Q_u", "S11_f0_dB", "S22_f0_dB", "S21_background_dB", "csv"]
    ek = ["f0_Hz", "Q", "Q_u", "Q_c", "Q_d", "Q_ext", "beta", "S21_f0_dB_pred", "p_alumina", "p_ptfe"]
    loops = []
    for w, h in [(3, 3), (4, 4), (5, 5), (8, 6), (10, 8)]:
        n = f"xy{w}x{h}"
        loops.append(dict(loop_depth_mm=w, loop_width_mm=h, area_mm2=w * h,
                          eigen_o2=sel(recs.get(f"e_{n}_o2", {}), ek), driven_o2=sel(recs.get(f"d_{n}_o2", {}), dk)))
    d2, d3, dr = recs[f"d_{mn}_o2"], recs[f"d_{mn}_o3"], recs[f"d_{mn}_r_o2"]
    conv = dict(
        p_refinement=dict(dS21_dB=d3["S21_f0_dB"] - d2["S21_f0_dB"], df0_rel=(d3["f0_Hz"] - d2["f0_Hz"]) / d2["f0_Hz"],
                          dQu_rel=(d3["Q_u"] - d2["Q_u"]) / d2["Q_u"]),
        h_refinement=dict(dS21_dB=dr["S21_f0_dB"] - d2["S21_f0_dB"], df0_rel=(dr["f0_Hz"] - d2["f0_Hz"]) / d2["f0_Hz"],
                          dQu_rel=(dr["Q_u"] - d2["Q_u"]) / d2["Q_u"]),
        criteria="|dS21| < 0.5 dB and |df0| < 0.05 %")
    base3 = recs["e_base_o3"]
    tune = [dict(penetration_mm=0, f0_Hz=base3["f0_Hz"], Q=base3["Q"], df0_Hz=0.0, dQ_rel=0.0,
                 gap_to_lip_tip_mm=sz.get("clearance_mm", 8.086))]
    for p in (2, 4, 6):
        r = recs.get(f"e_screw{p}_o3")
        if r:
            tune.append(dict(penetration_mm=p, f0_Hz=r["f0_Hz"], Q=r["Q"], df0_Hz=r["f0_Hz"] - base3["f0_Hz"],
                             dQ_rel=(r["Q"] - base3["Q"]) / base3["Q"], gap_to_lip_tip_mm=8.086 - p))
    fil = recs.get("e_fillet_o3")
    fillet = dict(radius_mm=0.4, f0_sharp_Hz=base3["f0_Hz"], Q_sharp=base3["Q"], Qc_sharp=base3["Q_c"]) if fil else None
    if fil:
        fillet.update(f0_fillet_Hz=fil["f0_Hz"], Q_fillet=fil["Q"], Qc_fillet=fil["Q_c"], df0_Hz=fil["f0_Hz"] - base3["f0_Hz"],
                      df0_rel=(fil["f0_Hz"] - base3["f0_Hz"]) / base3["f0_Hz"], dQ_rel=(fil["Q"] - base3["Q"]) / base3["Q"])
    b2 = recs["e_base_o2"]
    axial = [dict(rod_diameter_mm=6, penetration_mm=pz, f0_Hz=r["f0_Hz"], Q=r["Q"], df0_Hz=r["f0_Hz"] - b2["f0_Hz"],
                  dQ_rel=(r["Q"] - b2["Q"]) / b2["Q"]) for pz in (20, 30) if (r := recs.get(f"e_axial6_{pz}_o2"))]
    nosp = recs["e_base_nosp_o2"]
    # coupling-wire geometry for assembly (centreline)
    Rcan = 43.0127
    sag = Rcan - math.sqrt(Rcan ** 2 - H ** 2)
    out = dict(
        cell=dict(V_l=0.5, a=0.5, d_mm=1.0, s_mm=10.42, R_mm=Rcan, H_mm=2 * Rcan, r_mm=21.5064, z_mm=64.519),
        model=dict(ptfe_spacers="included (eps_r 2.1, tan_delta 2e-4)", port="50 ohm coaxial lumped port, 1 mm behind the inner wall, 50 ohm air-coax recess",
                   copper_sigma=5.8e7, brass_sigma=1.5e7, alumina="eps_r 9.8, tan_delta 4e-5"),
        recommended_loop=dict(orientation="mid-height plane (normal = can axis)", depth_mm=W, width_mm=H, wire_diameter_mm=1.0,
                              centreline_length_in_can_mm=W + H + W + sag, ground_point_offset_mm=H,
                              ground_point_angle_offset_deg=math.degrees(math.atan2(H, Rcan - sag)),
                              ports_deg=[120, 240], tangential_leg_towards_deg=180,
                              S21_f0_dB=d3["S21_f0_dB"], beta1=d3["beta1"], beta2=d3["beta2"], Q_L=d3["Q_L_3dB"], Q_u=d3["Q_u"], f0_Hz=d3["f0_Hz"],
                              sensitivity=dict(beta_vs_area_exponent=sz["fit"]["n"], area_for_minus30dB_mm2=sz["area_for_target_mm2"],
                                               dS21_dB_per_percent_area=20 * sz["fit"]["n"] * math.log10(1.01))),
        loops=loops,
        rz_orientation_5x5=sel(recs["e_rz5x5_o2"], ek),
        chosen_runs=dict(driven_o2=sel(d2, dk), driven_o3=sel(d3, dk), driven_o2_refined=sel(dr, dk),
                         eigen_o2=sel(recs[f"e_{mn}_o2"], ek), eigen_o3=sel(recs[f"e_{mn}_o3"], ek), eigen_o2_refined=sel(recs[f"e_{mn}_r_o2"], ek)),
        convergence=conv, prom_check=prom_check(),
        spacers=dict(f0_without_Hz=nosp["f0_Hz"], f0_with_Hz=b2["f0_Hz"], df0_Hz=b2["f0_Hz"] - nosp["f0_Hz"],
                     Q_without=nosp["Q"], Q_with=b2["Q"], p_ptfe=b2.get("p_ptfe"), order=2),
        tuning_screw_M4=dict(order=3, runs=tune), axial_rod_tuner=dict(order=2, runs=axial), fillet=fillet,
        as_built=dict(description="4x4 mm loops + R0.4 mm lip-tip fillets; prediction to compare with the VNA measurement",
                      driven_o3=sel(recs["d_xy4x4_fil_o3"], dk), driven_o2=sel(recs["d_xy4x4_fil_o2"], dk)),
        slit_air_film_estimate=dict(note="analytic, not simulated: a total air film delta in series with the 1.0 mm alumina gap "
                                    "scales the gap capacitance by 1/(1 + eps_r delta/d)",
                                    df0_per_um_Hz=0.5 * 9.8 / 1000.0 * 145.2e6,
                                    f0_ratio_25um=math.sqrt(1 + 9.8 * 0.025 / 1.0)),
        shr_spent=0.0)
    (HERE / "results" / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    {"sizing": sizing, "summary": summary}[sys.argv[1]](records())
