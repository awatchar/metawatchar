#!/usr/bin/env python3
"""Collect the S2 deliverables from the per-run outputs (local, light file work on pulled results).

Writes
  results/S_N1.csv, S_N2.csv, S_N4.csv   production S-parameters, columns f_Hz,S11_re,S11_im,S21_re,S21_im
                                         (de-embedded to the slab faces x = 0 and x = L)
  results/S_N1_orient_x.csv              opposite capacitor orientation (also S12, S22)
  (N = 1 mesh-convergence levels: results/N1_z_o2_S.csv, N1_z_o3_S.csv, N1_z_v3_o2_S.csv)
  results/s2_slab.json                   lattice, cell dimensions, slab length L per N, reference planes,
                                         materials, conventions, eigenmode results, run list
  results/driven_summary.csv, eigen_summary.csv   all runs, one row each (mode)
  results/s21_overview.png               |S21|, |S11| and absorption for N = 1, 2, 4 (if matplotlib)

Usage: python collect_s2.py
"""
import csv
import json
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / "results"

# production runs (order 3 on the v2 mesh where available; see README "Convergence")
PRODUCTION = {1: "N1_z_o3", 2: "N2_z_o3", 4: "N4_z_o3"}
FALLBACK = {1: "N1_z_o2", 2: "N2_z_o2", 4: "N4_z_o2"}   # also written as S_N<n>_o2.csv (consistent order-2 set)
ORIENT_X = "N1_x_o2"
EIGEN = ["eig_z_s0620_pmc_o3", "eig_z_s0620_pec_o3", "eig_z_s0620_pmc_o2", "eig_z_s0620_pec_o2"]


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def merge(pattern, out):
    allr = []
    for p in sorted(RES.glob(pattern)):
        if p.name in ("driven_summary.csv", "eigen_summary.csv"):
            continue
        allr += rows(p)
    if not allr:
        return []
    keys = []
    for r in allr:
        keys += [k for k in r if k not in keys]
    with open(RES / out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(allr)
    return allr


def main():
    drv = {r["run"]: r for r in merge("*_summary.csv", "driven_summary.csv")}
    eig = merge("*_eig.csv", "eigen_summary.csv")
    used = {}
    for n in (1, 2, 4):
        run = PRODUCTION[n] if (RES / f"{PRODUCTION[n]}_S.csv").exists() else FALLBACK[n]
        if (RES / f"{run}_S.csv").exists():
            shutil.copyfile(RES / f"{run}_S.csv", RES / f"S_N{n}.csv")
            used[n] = run
    for n, run in FALLBACK.items():
        if (RES / f"{run}_S.csv").exists():
            shutil.copyfile(RES / f"{run}_S.csv", RES / f"S_N{n}_o2.csv")
    if (RES / f"{ORIENT_X}_S.csv").exists():
        shutil.copyfile(RES / f"{ORIENT_X}_S.csv", RES / "S_N1_orient_x.csv")

    cases = {p.stem: json.loads(p.read_text()) for p in (HERE / "cases").glob("*.json")}
    runs = {p.stem: json.loads(p.read_text()) for p in (HERE / "runs").glob("*.json")}
    some = cases[runs[used[1]]["mesh"]]
    g = some["geometry"]
    meta = {
        "description": "S2: slab of N unshielded loop-gap (MNG) cells in a TEM unit-cell waveguide "
                       "(PEC z = +-a_z/2, PMC y = +-a_y/2 = infinite 2D array in y-z), propagation +x",
        "units": "lengths in mm, frequency in Hz",
        "time_convention": "exp(+j omega t): an empty line of length l has S21 = exp(-j k0 l) "
                           "(checked with the empty-line run)",
        "lattice_mm": {"a_x": g["a_x"], "a_y": g["a_y"], "a_z": g["a_z"]},
        "cell_mm": {"loop_inner_radius_r": g["r"], "loop_wall_t": g["t"], "loop_length_along_y": g["ell"],
                    "lip_thickness": g["t_lip"], "lip_height_s": g["s"], "gap_d_alumina": g["d"],
                    "lip_tip_clearance_to_plate": g["clearance"],
                    "loop_axis": "y (parallel to the TEM H-field)",
                    "orientation_production": "slit and lips towards +z: lip plates parallel to E (z), "
                                              "gap field along x (perpendicular to E)",
                    "cell_centre": "x = (i + 1/2) a_x, y = 0, z = 0, i = 0..N-1"},
        "materials": {"alumina": {"eps_r": some["eps_alumina"], "tan_delta": some["tand_alumina"]},
                      "copper": {"sigma_S_per_m": some["sigma_cu"], "Rs_ohm_sq": some["Rs"],
                                 "Ls_H_sq": some["Ls"], "Rs_evaluated_at_Hz": some["f_Rs"],
                                 "model": "surface impedance (frequency-independent Rs, Ls)"},
                      "plates_and_walls": "ideal PEC (z) and PMC (y), lossless"},
        "ports": "numeric wave ports (TEM mode, unit incident power), reference impedance = TEM line "
                 "impedance eta0 a_z / a_y; de-embedded by Palace (port Offset) to the slab faces",
        "slabs": {},
        "eigenmodes_single_cell": [],
        "eigenmode_note": "single cell in the same PEC/PMC box closed at x = 0, a_x: PMC x-walls = alternating "
                          "images (k = pi/a_x), just below the mu_eff pole f_m; PEC x-walls = in-phase images "
                          "(k = 0, zero macroscopic E), the mu_eff = 0 frequency f_mp. Q includes copper + alumina loss",
        "orientation_x_comparison": {"file": "S_N1_orient_x.csv", "run": ORIENT_X,
                                     "description": "N = 1, slit and lips towards +x (gap field parallel to E); "
                                                    "both ports excited (S11, S21, S12, S22); bianisotropic"},
        "frequency_grid": "130-200 MHz, 10 kHz steps (7001 points), Palace adaptive fast sweep (AdaptiveTol 1e-4)",
        "runs": {},
    }
    for n, run in used.items():
        gg = cases[runs[run]["mesh"]]["geometry"]
        meta["slabs"][f"N{n}"] = {
            "file": f"S_N{n}.csv", "run": run, "N": n, "slab_length_L_mm": gg["L"],
            "slab_faces_x_mm": gg["slab_faces"], "reference_planes_x_mm": gg["slab_faces"],
            "port_planes_x_mm": [gg["x_min"], gg["x_max"]], "empty_line_each_side_mm": gg["l_empty"],
            "deembedding_offset_mm": runs[run]["offset_mm"], "fe_order": runs[run]["order"],
            "mesh_recipe": cases[runs[run]["mesh"]]["mesh"]["recipe"],
            "n_tet": cases[runs[run]["mesh"]]["mesh"]["n_tet"], "dofs": drv.get(run, {}).get("dofs")}
    for r in eig:
        if r["run"] in EIGEN and int(r["mode"]) == min(int(x["mode"]) for x in eig if x["run"] == r["run"]):
            meta["eigenmodes_single_cell"].append({
                "run": r["run"], "x_walls": r["xwall"], "fe_order": int(r["order"]), "f_Hz": float(r["f_MHz"]) * 1e6,
                "Q": float(r["Q"]), "Q_c": float(r["Q_c"]), "Q_d": float(r["Q_d"]), "p_d": float(r["p_d"])})
    for run, rm in runs.items():
        if run in drv:
            meta["runs"][run] = {k: drv[run][k] for k in ("N", "orient", "order", "dofs", "f_S21min_Hz",
                                                          "S21min_dB", "absorption_max", "f_absmax_Hz")}
    (RES / "s2_slab.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({k: v["run"] for k, v in meta["slabs"].items()}))

    try:
        import math
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(3, 1, figsize=(7, 8), sharex=True)
    for n, run in used.items():
        d = rows(RES / f"S_N{n}.csv")
        f = [float(r["f_Hz"]) / 1e6 for r in d]
        s11 = [complex(float(r["S11_re"]), float(r["S11_im"])) for r in d]
        s21 = [complex(float(r["S21_re"]), float(r["S21_im"])) for r in d]
        ax[0].plot(f, [20 * math.log10(max(abs(v), 1e-12)) for v in s21], label=f"N = {n}")
        ax[1].plot(f, [20 * math.log10(max(abs(v), 1e-12)) for v in s11], label=f"N = {n}")
        ax[2].plot(f, [1 - abs(a) ** 2 - abs(b) ** 2 for a, b in zip(s11, s21)], label=f"N = {n}")
    ax[0].set_ylabel("|S21| (dB)")
    ax[1].set_ylabel("|S11| (dB)")
    ax[2].set_ylabel("1 - |S11|^2 - |S21|^2")
    ax[2].set_xlabel("f (MHz)")
    for a in ax:
        a.grid(alpha=0.3)
        a.legend()
    fig.tight_layout()
    fig.savefig(RES / "s21_overview.png", dpi=120)


if __name__ == "__main__":
    main()
