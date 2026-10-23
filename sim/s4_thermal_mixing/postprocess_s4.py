#!/usr/bin/env python3
"""Summarize Palace eigenmode results of the S4 runs (runs inside the Slurm job; standard library).

For every run (runs/<run>.json) -> results/<run>_eig.csv, one row per physical mode (Q > 10; the
impedance-boundary eigenproblem also returns critically damped spurious values) with
  f (MHz), Q, the electric (pe_<domain>) and magnetic (pm_<domain>) energy shares of every volume
  domain (air, bore, end, c<j>_al, c<j>_ti, c<j>_gap), the dielectric Q of each material
  Q_d_al = 1 / sum_j(pe_cj_al) tan_d_al, Q_d_ti likewise, and Q_c = (1/Q - 1/Q_d_al - 1/Q_d_ti)^-1
  (copper: loops, plus the box walls where they are copper).
Usage: python3 postprocess_s4.py run1 [run2 ...]
"""
import csv
import json
import re
import sys
from pathlib import Path


def read_csv(path):
    with open(path) as fh:
        rows = [[c.strip() for c in line.split(",")] for line in fh if line.strip()]
    return rows[0], [[float(x) for x in r] for r in rows[1:]]


def col(hdr, pattern):
    for i, h in enumerate(hdr):
        if re.fullmatch(pattern, h):
            return i
    return None


def one(run):
    meta = json.loads((Path("runs") / f"{run}.json").read_text())
    out = Path("postpro") / run
    hdr, eig = read_csv(out / "eig.csv")
    i_m, i_f, i_q = col(hdr, r"m"), col(hdr, r"Re\{f\} \(GHz\)"), col(hdr, r"Q")
    i_fi = col(hdr, r"Im\{f\} \(GHz\)")
    try:
        dh, dom = read_csv(out / "domain-E.csv")
    except FileNotFoundError:
        dh, dom = None, []
    try:
        dofs = json.loads((out / "palace.json").read_text())["Problem"]["DegreesOfFreedom"]
    except (FileNotFoundError, KeyError, ValueError):
        dofs = ""
    tand = {"al": meta["mat"]["al"]["tand"], "ti": meta["mat"]["ti"]["tand"]}
    rows = []
    for r in sorted(eig, key=lambda r: r[i_f]):
        q = r[i_q]
        if not q > 10.0:
            continue
        res = dict(run=run, mesh=meta["mesh"], order=meta["order"], dT=meta["dT"], dofs=dofs,
                   mode=int(r[i_m]), f_MHz=r[i_f] * 1e3, f_im_MHz=r[i_fi] * 1e3, Q=q)
        drow = next((x for x in dom if x[0] == r[i_m]), None)
        inv_qd = {"al": 0.0, "ti": 0.0}
        if drow is not None:
            for idx, lab in meta["energy_index"].items():
                pe, pm = col(dh, rf"p_elec\[{idx}\]"), col(dh, rf"p_mag\[{idx}\]")
                res[f"pe_{lab}"] = drow[pe]
                res[f"pm_{lab}"] = drow[pm]
                mat = lab.split("_")[-1]
                if mat in inv_qd:
                    inv_qd[mat] += drow[pe] * tand[mat]
            res["E_elec_J"] = drow[col(dh, r"E_elec \(J\)")]
            res["E_mag_J"] = drow[col(dh, r"E_mag \(J\)")]
        res["Q_d_al"] = 1 / inv_qd["al"] if inv_qd["al"] > 0 else float("inf")
        res["Q_d_ti"] = 1 / inv_qd["ti"] if inv_qd["ti"] > 0 else float("inf")
        res["Q_c"] = 1.0 / (1.0 / q - inv_qd["al"] - inv_qd["ti"])
        rows.append(res)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(f"results/{run}_eig.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.10g}" if isinstance(v, float) else v) for k, v in r.items()})
    for res in rows:
        print(f"{run:<34} mode {res['mode']}: f={res['f_MHz']:10.5f} MHz  Q={res['Q']:8.1f}  "
              f"Qc={res['Q_c']:8.1f}")


def main():
    Path("results").mkdir(exist_ok=True)
    for run in sys.argv[1:]:
        try:
            one(run)
        except Exception as exc:
            print(f"{run}: FAILED ({exc!r})", file=sys.stderr)


if __name__ == "__main__":
    main()
