#!/usr/bin/env python3
"""Summarize Palace eigenmode results for LGR-A cases (runs inside the Slurm job; numpy only).

For each case name: reads cases/<name>.json (geometry, probe labels), postpro/<name>/eig.csv,
domain-E.csv and probe-E.csv, takes the lowest mode (the loop-gap mode), and writes one row:

  f0 (MHz), Q (total), p_d (electric-energy participation of the alumina),
  Q_d = 1/(p_d tan_delta), Q_c = (1/Q - 1/Q_d)^-1 (conductor loss: shield + loop + lips),
  W = E_elec (J, Palace's 1/2 Re int D*.E dV of the mode = total stored energy),
  |E| at every probe scaled to a stored energy of W_REF = 10 uJ (kV/mm, peak phasor amplitude),
  V_gap = |E(gap_wall_mid)| d at W_REF and the effective gap capacitance C = 2 W / V^2.

Usage: python3 postprocess.py OUT.csv name1 [name2 ...]
"""
import csv
import json
import math
import re
import sys
from pathlib import Path

W_REF = 10e-6


def read_csv(path):
    with open(path) as fh:
        rows = [[c.strip() for c in line.split(",")] for line in fh if line.strip()]
    hdr = rows[0]
    data = [[float(x) for x in r] for r in rows[1:]]
    return hdr, data


def col(hdr, pattern):
    for i, h in enumerate(hdr):
        if re.fullmatch(pattern, h):
            return i
    raise KeyError(pattern)


def summarize(name, sub=None):
    """sub: AMR iteration subdirectory (e.g. "iteration1"); the row is then named <name>@<sub>."""
    meta = json.loads(Path(f"cases/{name}.json").read_text())
    g = meta["geometry"]
    out = Path(f"postpro/{name}") / (sub or "")
    hdr, eig = read_csv(out / "eig.csv")
    i_m, i_f, i_q = col(hdr, r"m"), col(hdr, r"Re\{f\} \(GHz\)"), col(hdr, r"Q")
    # drop spurious, critically damped eigenvalues (Q < 10) that the quadratic eigenproblem with
    # impedance boundaries can return; the loop-gap mode has Q ~ 1e3
    eig = [r for r in eig if r[i_q] > 10.0]
    eig.sort(key=lambda r: r[i_f])
    m = eig[0][i_m]
    f0 = eig[0][i_f] * 1e3
    q = eig[0][i_q]
    f_next = eig[1][i_f] * 1e3 if len(eig) > 1 else float("nan")

    hdr, dom = read_csv(out / "domain-E.csv")
    row = next(r for r in dom if r[0] == m)
    w = row[col(hdr, r"E_elec \(J\)")]
    w_mag = row[col(hdr, r"E_mag \(J\)")]
    p_d = row[col(hdr, r"p_elec\[1\]")]
    q_d = 1.0 / (p_d * meta["tand_alumina"])
    q_c = 1.0 / (1.0 / q - 1.0 / q_d)

    try:
        dofs = json.loads((out / "palace.json").read_text())["Problem"]["DegreesOfFreedom"]
    except (FileNotFoundError, KeyError, ValueError):
        dofs = ""
    res = dict(name=name + (f"@{sub}" if sub else ""), V_l=g["V_l"], a=g["a"], R_mm=g["R"], r_mm=g["r"], z_mm=g["z"], s_mm=g["s"],
               d_mm=g["d"], clearance_mm=g["clearance"], order=meta["order"],
               shield_pec=int(meta.get("shield_pec", False)),
               n_tet=meta["mesh"].get("n_tet", ""), h_scale=meta["mesh"].get("h_scale", ""), dofs=dofs,
               f0_MHz=f0, f_next_MHz=f_next, Q=q, Q_c=q_c, Q_d=q_d, p_d=p_d, W_J=w, W_mag_J=w_mag)

    scale = math.sqrt(W_REF / w)
    try:
        hdr, prb = read_csv(out / "probe-E.csv")
        prow = next(r for r in prb if r[0] == m)
        comps = {}
        for i, h in enumerate(hdr):
            mm = re.match(r"(Re|Im)\{E_([xyz])\[(\d+)\]\}", h)   # e.g. "Re{E_x[1]} (V/m)"
            if mm:
                comps.setdefault(int(mm.group(3)), 0.0)
                comps[int(mm.group(3))] += prow[i] ** 2
        for p in meta["probes"]:
            e = math.sqrt(comps.get(p["index"], float("nan"))) * scale
            res[f"E_{p['label']}_kV_mm"] = e / 1e6
        v_gap = res["E_gap_wall_mid_kV_mm"] * 1e6 * g["d"] * 1e-3
        res["V_gap_V"] = v_gap
        res["C_eff_pF"] = 2 * W_REF / v_gap ** 2 * 1e12
    except (FileNotFoundError, StopIteration, KeyError) as exc:
        print(f"{name}: probes unavailable ({exc})", file=sys.stderr)
    return res


def main():
    out_csv = sys.argv[1]
    rows = []
    for name in sys.argv[2:]:
        for it in sorted(Path(f"postpro/{name}").glob("iteration*"), key=lambda q: int(q.name[9:])):
            try:
                rows.append(summarize(name, it.name))      # AMR history
            except Exception as exc:
                print(f"{name}/{it.name}: FAILED ({exc!r})", file=sys.stderr)
        try:
            rows.append(summarize(name))
        except Exception as exc:  # keep going; report the failed case
            print(f"{name}: FAILED ({exc!r})", file=sys.stderr)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    for r in rows:
        print(f"{r['name']:<24} f0={r['f0_MHz']:8.3f} MHz  Q={r['Q']:7.0f}  Qc={r['Q_c']:7.0f}  "
              f"Qd={r['Q_d']:8.0f}  p_d={r['p_d']:.3f}  E_gap={r.get('E_gap_wall_mid_kV_mm', float('nan')):.3f} kV/mm")


if __name__ == "__main__":
    main()
