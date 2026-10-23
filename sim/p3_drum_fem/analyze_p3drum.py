"""P3 drum-inductor FEM: analysis of the Palace inductance matrices (light work, local).

Reads results/<run>_M.csv (pulled from the cluster) and runs_meta.json and writes
results/summary.csv plus a printed report:
  val_*    thin-ring mutual inductance vs the Neumann double integral (filaments on the
           ring centre-lines), the validation required before trusting the setup
  single_* series inductance of the K representative turns, L_rep; the turn count implied by
           the datasheet inductance, N_L = K sqrt(L_ds / L_rep); an independent estimate N_DCR
           from the datasheet DCR (wire length vs window fill); apparent core gain L_rep/L_air
  pair_*   coupling k = M_ab / sqrt(L_aa L_bb) between two cells (independent of N),
           m^2/L implied by k (k = mu0 m^2 / (4 pi d^3 L) for coplanar cells), and the ratio
           to the point-dipole law, which shows where the lattice model breaks down
Run: python analyze_p3drum.py
"""
import csv
import glob
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MU0 = 4e-7 * math.pi
RHO_CU = 1.72e-8


def read_M(run):
    p = os.path.join(HERE, "results", f"{run}_M.csv")
    if not os.path.exists(p):
        return None
    rows = list(csv.reader(open(p)))
    return np.array([[float(x) for x in r[1:]] for r in rows[1:]])


def neumann(r1, c1, r2, c2, n=2000):
    """Mutual inductance of two coaxial-axis (both along z) circular filaments, centres c1, c2."""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    p1 = np.stack([c1[0] + r1 * np.cos(t), c1[1] + r1 * np.sin(t), np.full(n, c1[2])], 1)
    dl1 = np.stack([-r1 * np.sin(t), r1 * np.cos(t), np.zeros(n)], 1) * (2 * np.pi / n)
    p2 = np.stack([c2[0] + r2 * np.cos(t), c2[1] + r2 * np.sin(t), np.full(n, c2[2])], 1)
    dl2 = np.stack([-r2 * np.sin(t), r2 * np.cos(t), np.zeros(n)], 1) * (2 * np.pi / n)
    M = 0.0
    for i in range(0, n, 250):
        d = np.linalg.norm(p1[i:i + 250, None, :] - p2[None, :, :], axis=2)
        M += np.sum((dl1[i:i + 250] @ dl2.T) / d)
    return MU0 / (4 * np.pi) * M * 1e-3      # coordinates in mm -> m


def n_from_dcr(m):
    """Turn count from DCR with a square-packed winding filling 60 % of the window."""
    fill = 0.6
    a_win = ((m["D"] / 2 - 0.35) - (m["DP"] / 2 + 0.25)) * (m["H"] - 2 * m["TF"] - 0.6) * 1e-6
    r_mean = ((m["D"] / 2 - 0.35) + (m["DP"] / 2 + 0.25)) / 2 * 1e-3
    return math.sqrt(m["DCR"] * fill * a_win / (2 * math.pi * RHO_CU * r_mean))


def main():
    meta = json.load(open(os.path.join(HERE, "runs_meta.json")))
    out = []
    print("== validation: thin rings vs Neumann (filament) ==")
    for run, m in meta.items():
        if m["kind"] != "val":
            continue
        M = read_M(run)
        if M is None:
            print(run, "no result")
            continue
        mn = neumann(m["r"], (-m["dx"] / 2, 0, -m["dz"] / 2), m["r"], (m["dx"] / 2, 0, m["dz"] / 2))
        err = M[0, 1] / mn - 1
        # effective radius of the PEC strip that reproduces the FEM mutual exactly
        from scipy.optimize import brentq
        f = lambda r: neumann(r, (-m["dx"] / 2, 0, -m["dz"] / 2), r, (m["dx"] / 2, 0, m["dz"] / 2), n=1200) - M[0, 1]
        try:
            r_eff = brentq(f, m["r"] - m["w"] / 2, m["r"] + m["w"] / 2, xtol=1e-4)
        except ValueError:
            r_eff = float("nan")
        print(f"{run}: M_fem = {M[0, 1] * 1e9:.4f} nH, M_neumann(centre-line) = {mn * 1e9:.4f} nH, "
              f"diff {err * 100:+.2f} %, effective radius {r_eff:.3f} mm, L11 = {M[0, 0] * 1e9:.3f} nH")
        out.append(dict(run=run, quantity="M_val", fem=M[0, 1], ref=mn, rel_err=err, r_eff=r_eff))

    print("\n== single cells ==")
    l_air = {}
    for run, m in meta.items():
        if m["kind"] == "single" and run.endswith("_air"):
            M = read_M(run)
            if M is not None:
                l_air[(m["type"], m["kr"], m["kz"])] = M.sum()
    for run, m in meta.items():
        if m["kind"] != "single":
            continue
        M = read_M(run)
        if M is None:
            print(run, "no result")
            continue
        K = M.shape[0]
        L_rep = M.sum()
        N_L = K * math.sqrt(m["L_uH"] * 1e-6 / L_rep)
        gain = L_rep / l_air.get((m["type"], m["kr"], m["kz"]), float("nan"))
        print(f"{run}: K={K} L_rep={L_rep * 1e9:.2f} nH  N_L={N_L:.0f}  N_DCR={n_from_dcr(m):.0f}  core gain {gain:.2f}")
        out.append(dict(run=run, quantity="L_rep", fem=L_rep, N_L=N_L, N_DCR=n_from_dcr(m), gain=gain))

    print("\n== pairs: coupling k (independent of N) ==")
    for run, m in sorted(meta.items(), key=lambda kv: (kv[1].get("type", ""), kv[1].get("geom", ""), kv[1].get("s", 0))):
        if m["kind"] != "pair":
            continue
        M = read_M(run)
        if M is None:
            print(run, "no result")
            continue
        K = M.shape[0] // 2
        La, Lb, Mab = M[:K, :K].sum(), M[K:, K:].sum(), M[:K, K:].sum()
        k = Mab / math.sqrt(La * Lb)
        line = f"{run}: k = {k:+.5f}"
        rec = dict(run=run, quantity="k", fem=k, type=m["type"], geom=m["geom"], s=m["s"])
        if m["geom"] == "copl":
            d = m["s"] * 1e-3
            m2_over_L = -k * 4 * math.pi * d ** 3 / MU0          # m^2/L in m^4/H (k < 0 coplanar)
            line += f"   m^2/L = {m2_over_L:.3e} m^4/H"
            rec["m2_over_L"] = m2_over_L
        print(line)
        out.append(rec)
    # dipole-law check: m^2/L should be constant with d if the point-dipole lattice model holds
    for typ in "ABC":
        rows = [r for r in out if r.get("quantity") == "k" and r.get("type") == typ and r.get("geom") == "copl"]
        if len(rows) >= 2:
            far = max(rows, key=lambda r: r["s"])["m2_over_L"]
            print(f"type {typ}: m^2/L(d) / m^2/L(far) = " +
                  ", ".join(f"{r['s']:g} mm: {r['m2_over_L'] / far:.3f}" for r in sorted(rows, key=lambda r: r["s"])))
    # 4 x 4 arrays: cell-level coupling matrix -> collective (bright) mode vs point-dipole lattice
    print("\n== 4 x 4 arrays: collective mode (f / f0 of the isolated cell) ==")
    for run, m in meta.items():
        if m["kind"] != "array":
            continue
        M = read_M(run)
        if M is None:
            print(run, "no result")
            continue
        K = m["kr"] * m["kz"]
        nc = M.shape[0] // K
        Lc = np.array([[M[i * K:(i + 1) * K, j * K:(j + 1) * K].sum() for j in range(nc)] for i in range(nc)])
        Kn = Lc / np.sqrt(np.outer(np.diag(Lc), np.diag(Lc)))          # normalised: 1 on the diagonal
        lam, vec = np.linalg.eigh(Kn)
        bright = (vec.sum(axis=0) ** 2) / nc
        ib = int(np.argmax(bright))
        f_ratio = 1 / math.sqrt(lam[ib])                                  # f_mode / f0 for identical cells
        p = m["pitch"] * 1e-3
        k_nn = Kn[0, 1] if abs(Kn[0, 1]) > 0 else float("nan")
        # point-dipole lattice with m^2/L taken from the far (50 mm) coplanar pair of the same type
        far = [r for r in out if r.get("quantity") == "k" and r.get("type") == m["type"]
               and r.get("geom") == "copl" and r.get("s") == 50.0]
        xs = [(i - 1.5) * p for i in range(4)]
        pos = np.array([(x, y) for x in xs for y in xs])   # same ordering as drum_geometry (i outer, j inner)
        line = (f"{run}: bright mode f/f0 = {f_ratio:.4f} (weight {bright[ib]:.2f}), "
                f"nearest-neighbour k = {k_nn:+.4f}")
        if far:
            m2L = far[0]["m2_over_L"]
            d = np.linalg.norm(pos[:, None] - pos[None], axis=2)
            np.fill_diagonal(d, np.inf)
            Kd = np.eye(nc) - MU0 * m2L / (4 * np.pi * d ** 3)
            lam_d, vec_d = np.linalg.eigh(Kd)
            ib_d = int(np.argmax((vec_d.sum(axis=0) ** 2) / nc))
            line += f"; point-dipole model f/f0 = {1 / math.sqrt(lam_d[ib_d]):.4f}"
        print(line)
        out.append(dict(run=run, quantity="array_bright_f_ratio", fem=f_ratio, type=m["type"], s=m["pitch"]))
    keys = sorted({k for r in out for k in r})
    with open(os.path.join(HERE, "results", "summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    main()
