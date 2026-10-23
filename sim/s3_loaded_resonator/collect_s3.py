#!/usr/bin/env python3
"""Collect the S3 eigenmode results into the deliverable tables (local, after the runs are copied
back; numpy + matplotlib).

Inputs : results/<run>_eig.csv (written by postprocess_s3.py in the job), runs/*.json, cases/*.json,
         results/em_predictions.csv (s3_effective_medium.py).
Outputs: results/eigen_all.csv            every physical mode of every run
         results/study1_line.csv          study 1: full-wave modes vs effective-medium predictions
         results/study2_cavity.csv        study 2/3: loaded cavities vs empty cavity
         results/s3_summary.json          everything above plus convergence and parameters
         results/study1_line.png          figure
"""
import csv
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
C0 = 299_792_458.0
MU0 = 4e-7 * math.pi
SIGMA = 5.8e7
F_M = 149.934059e6            # meta-atom PMC-box (mu pole) eigenmode, order 2
# models compared with each full-wave set: "s2" = copper entered as in S2 (half the copper loss,
# consistent with the S2 slab S-parameters and so with model A); "exact" = correct copper
MODELS = {"s2": ("A-ret-N1", "A-ret-N2", "A-ret-N4", "B-eig", "B-fit"),
          "exact": ("A-ret-N4", "B-eig-cu", "B-fit-cu")}
MAIN = {"s2": ("A-ret-N4", "B-eig"), "exact": ("A-ret-N4", "B-eig-cu")}


def rcsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def empty_box_modes(X, Y, Z, n=3):
    """All TE/TM modes of a closed rectangular box up to index n, sorted by f (f only)."""
    out = []
    for m in range(n + 1):
        for k in range(n + 1):
            for p in range(n + 1):
                if (m > 0) + (k > 0) + (p > 0) < 2:
                    continue
                f = C0 / 2 * math.sqrt((m / X) ** 2 + (k / Y) ** 2 + (p / Z) ** 2)
                out.append((f, (m, k, p)))
    return sorted(out)


def q_te101(a, b, d, f):
    """Copper-wall Q of the TE101-type mode (E along b, sinusoidal in a and d), Pozar eq. 6.46."""
    rs = math.sqrt(math.pi * f * MU0 / SIGMA)
    k = 2 * math.pi * f / C0
    eta = 376.730313
    return (k * a * d) ** 3 * b * eta / (2 * math.pi ** 2 * rs) / (2 * a ** 3 * b + 2 * b * d ** 3 + a ** 3 * d + a * d ** 3)


def plot_study1(s1, pred):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pal = {"fw": "#2a78d6", "A-ret-N4": "#eb6834", "B-eig-cu": "#1baf7a", "B-eig": "#1baf7a"}
    fig, axs = plt.subplots(3, 2, figsize=(10, 10.5), dpi=140, sharex=True)
    panels = ((None, ("A-ret-N4", "B-eig-cu"), "f (MHz)"),
              ("s2", ("A-ret-N4", "B-eig"), "Q (S2 copper input: half the loss)"),
              ("exact", ("B-eig-cu",), "Q (correct copper)"))
    for j, g in enumerate((35.0, 0.0)):
        for row_i, (zs, models, ylabel) in enumerate(panels):
            ax = axs[row_i, j]
            want = "exact" if zs is None else zs
            fw = [r for r in s1 if abs(r["g_mm"] - g) < 1e-6 and int(r["order"]) == 2 and r["zs"] == want]
            key = "f_MHz" if row_i == 0 else "Q"
            for model in models:
                mk = "x" if model.startswith("A") else "+"
                pts = [(int(p["N"]), num(p["f_MHz"]), num(p["Q"])) for p in pred
                       if p["model"] == model and abs(num(p["g_mm"]) - g) < 1e-6 and 120 <= num(p["f_MHz"]) <= 200]
                dx = -0.08 if model.startswith("A") else 0.08
                ax.scatter([x[0] + dx for x in pts], [x[1] if row_i == 0 else x[2] for x in pts], marker=mk, s=40,
                           color=pal[model], label=f"effective medium {model}")
            ax.scatter([r["N"] for r in fw], [r[key] for r in fw], s=60, facecolors="none", edgecolors=pal["fw"],
                       lw=1.5, label="full wave (Palace)")
            ax.grid(True, color="#e4e3df", lw=0.6)
            ax.set_ylabel(ylabel, fontsize=8.5)
            if row_i == 0:
                ax.axhspan(F_M / 1e6, 179.516, color=pal["fw"], alpha=0.06, lw=0)
                ax.set_title(f"shorted line, end gaps g = {g:.0f} mm", fontsize=10)
            elif fw:
                ax.axhline(fw[0]["Q_cell"], color="#52514e", lw=0.8, ls="--")
                ax.text(0.75, fw[0]["Q_cell"] * 1.005, "meta-atom Q", fontsize=7, color="#52514e")
            if row_i == 2:
                ax.set_xlabel("N (cells)")
                ax.set_xticks([1, 2, 3, 4])
            ax.legend(fontsize=6.5, loc="best")
    fig.suptitle("S3 study 1: shorted-line resonator loaded with N meta-atoms (shaded: MNG band of the retrieved mu_eff)",
                 fontsize=9.5)
    fig.tight_layout()
    fig.savefig(RES / "study1_line.png")

def main():
    runs = {p.stem: json.loads(p.read_text()) for p in (HERE / "runs").glob("*.json")}
    cases = {p.stem: json.loads(p.read_text()) for p in (HERE / "cases").glob("*.json")}
    rows = []
    for p in sorted(RES.glob("*_eig.csv")):
        rows += rcsv(p)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(RES / "eigen_all.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    by_run = {}
    for r in rows:
        by_run.setdefault(r["run"], []).append(r)
    summary = dict(runs={k: dict(runs[k], modes=v) for k, v in by_run.items() if k in runs})

    # ---------------- study 1 ----------------
    pred = rcsv(RES / "em_predictions.csv") if (RES / "em_predictions.csv").exists() else []
    s1 = []
    for run, modes in sorted(by_run.items()):
        meta = runs.get(run)
        if not meta:
            continue
        g = cases[meta["mesh"]]["geometry"]
        if g["kind"] != "line" or meta["xwall"] != "pec" or abs(g["cell"]["scale"] - 1) > 1e-9:
            continue
        zs = meta.get("zs", "s2")
        ref = by_run.get(f"L1_g00_pmc_o{meta['order']}" + ("_zx" if zs == "exact" else ""))
        q_cell = num(ref[0]["Q"]) if ref else float("nan")
        fw = sorted([m for m in modes if 120 <= num(m["f_MHz"]) <= 200], key=lambda m: num(m["f_MHz"]))
        for i, m in enumerate(fw):
            f = num(m["f_MHz"])
            row = dict(run=run, zs=zs, N=g["N"], g_mm=g["g"], order=meta["order"], mode_index=i + 1, f_MHz=f,
                       Q=num(m["Q"]), Q_c=num(m.get("Q_c")), Q_cu=num(m.get("Q_cu")), Q_cell=q_cell,
                       Q_over_Qcell=num(m["Q"]) / q_cell,
                       pm_bore=num(m.get("pm_bore")), pm_air=num(m.get("pm_air")), pm_end=num(m.get("pm_end")),
                       p_d=num(m.get("p_d")))
            for model in MODELS[zs]:
                cand = [num(p["f_MHz"]) for p in pred if int(p["N"]) == g["N"] and abs(num(p["g_mm"]) - g["g"]) < 1e-6
                        and p["model"] == model]
                cq = {num(p["f_MHz"]): num(p["Q"]) for p in pred if int(p["N"]) == g["N"]
                      and abs(num(p["g_mm"]) - g["g"]) < 1e-6 and p["model"] == model}
                if not cand:
                    continue
                fp = min(cand, key=lambda x: abs(x - f))
                row[f"f_{model}"] = fp
                row[f"df_pct_{model}"] = (fp - f) / f * 100
                if not (zs == "exact" and model.startswith("A")):   # model A carries the S2 (half) copper loss
                    row[f"Q_{model}"] = cq[fp]
                row[f"n_pred_{model}"] = len([x for x in cand if 120 <= x <= 200])
            row["shift_from_fm_MHz"] = f - F_M / 1e6
            s1.append(row)
    if s1:
        k1 = []
        for r in s1:
            k1 += [k for k in r if k not in k1]
        with open(RES / "study1_line.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=k1)
            w.writeheader()
            for r in s1:
                w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    summary["study1"] = s1

    # ---------------- study 2 / 3 ----------------
    s2 = []
    for run, modes in sorted(by_run.items()):
        meta = runs.get(run)
        if not meta:
            continue
        g = cases[meta["mesh"]]["geometry"]
        if g["kind"] != "cavity":
            continue
        X, Y, Z = (v * 1e-3 for v in g["size_mm"])
        em = empty_box_modes(X, Y, Z)
        f_e, idx = em[0]
        dimsorted = sorted([X, Y, Z])
        q_e = q_te101(dimsorted[1], dimsorted[0], dimsorted[2], f_e)
        fw = sorted(modes, key=lambda m: num(m["f_MHz"]))
        for i, m in enumerate(fw):
            f = num(m["f_MHz"]) * 1e6
            lam = C0 / f
            s2.append(dict(run=run, array=f"{g['nx']}x{g['ny']}x{g['nz']}", N=g["N"], ext_mm=g["ext"],
                           scale=g["cell"]["scale"], walls=meta["walls"], order=meta["order"],
                           box_mm="x".join(f"{v:.1f}" for v in g["size_mm"]), V_L=X * Y * Z * 1e3,
                           zs=meta.get("zs", "s2"), mode_index=i + 1, f_MHz=f / 1e6, Q=num(m["Q"]),
                           Q_c=num(m.get("Q_c")), Q_cu=num(m.get("Q_cu")), Q_c_cu=num(m.get("Q_c_cu")),
                           p_d=num(m.get("p_d")),
                           pm_bore=num(m.get("pm_bore")), pm_air=num(m.get("pm_air")), pe_air=num(m.get("pe_air")),
                           f_empty_MHz=f_e / 1e6, empty_mode=str(idx), Q_empty_analytic=q_e,
                           miniaturization=f_e / f, max_dim_over_lambda=max(X, Y, Z) / lam,
                           V_over_lambda3=X * Y * Z / lam ** 3))
    if s2:
        k2 = []
        for r in s2:
            k2 += [k for k in r if k not in k2]
        with open(RES / "study2_cavity.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=k2)
            w.writeheader()
            for r in s2:
                w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    summary["study2"] = s2
    # empty-cavity Palace check
    for run, modes in by_run.items():
        meta = runs.get(run)
        if meta and cases[meta["mesh"]]["geometry"]["kind"] == "empty":
            g = cases[meta["mesh"]]["geometry"]
            X, Y, Z = (v * 1e-3 for v in g["size_mm"])
            if run.endswith("allcu_o2"):      # y = 0 plane is copper too: the box is 70 x 35 x 140 mm
                Y = Y / 2
            f_e = empty_box_modes(X, Y, Z)[0][0]
            ds = sorted([X, Y, Z])
            summary.setdefault("empty_check", []).append(dict(run=run, zs=meta.get("zs", "s2"), f_palace_MHz=num(modes[0]["f_MHz"]), Q_palace=num(modes[0]["Q"]),
                                          f_analytic_MHz=f_e / 1e6, Q_analytic=q_te101(ds[1], ds[0], ds[2], f_e),
                                          Q_cu=num(modes[0].get("Q_cu")), box_mm=[X * 1e3, Y * 1e3, Z * 1e3],
                                          note="Rs at 2.4 GHz in the run; analytic Q with Rs at the analytic f"))
    if s1:
        plot_study1(s1, pred)
    (RES / "s3_summary.json").write_bytes((json.dumps(summary, indent=1, default=str) + "\n").encode())

    # ---------------- printout ----------------
    print("STUDY 1")
    for r in s1:
        print(f"{r['run']:<14} N={r['N']} g={r['g_mm']:4.0f} #{r['mode_index']} f={r['f_MHz']:9.4f} Q={r['Q']:7.0f} "
              f"Q/Qcell={r['Q_over_Qcell']:.3f} pm_bore={r['pm_bore']:.3f} pm_end={r['pm_end']:.3f} | "
              + " ".join(f"{m}:{r.get('f_' + m, float('nan')):.4f}({r.get('df_pct_' + m, float('nan')):+.3f}%,Q{r.get('Q_' + m, float('nan')):.0f})"
                         for m in MAIN[r["zs"]]))
    print("STUDY 2/3")
    for r in s2:
        print(f"{r['run']:<18} {r['array']} k={r['scale']} ext={r['ext_mm']} walls={r['walls']} #{r['mode_index']} "
              f"f={r['f_MHz']:9.3f} Q={r['Q']:7.0f} Qc={r['Q_c']:7.0f} Qcu={r['Q_cu']:7.0f} pm_bore={r['pm_bore']:.3f} "
              f"f_empty={r['f_empty_MHz']:.1f} {r['empty_mode']} mini={r['miniaturization']:.2f} "
              f"Dmax/lam={r['max_dim_over_lambda']:.4f} V={r['V_L']:.3f} L")
    if "empty_check" in summary:
        print(summary["empty_check"])


if __name__ == "__main__":
    main()
