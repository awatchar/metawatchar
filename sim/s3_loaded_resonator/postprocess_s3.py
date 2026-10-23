#!/usr/bin/env python3
"""Summarize Palace eigenmode results of the S3 runs (runs inside the Slurm job; standard library).

For every run (runs/<run>.json + cases/<mesh>.json) -> results/<run>_eig.csv, one row per physical
mode (Q > 10; the impedance-boundary eigenproblem also returns critically damped spurious values):
  f (MHz), Q, p_d (alumina share of the electric energy), Q_d = 1/(p_d tan_delta),
  Q_c = (1/Q - 1/Q_d)^-1 (all copper: loops, plus cavity walls where they are copper),
  Q_c_cu and Q_cu: Q_c and Q rescaled to the true copper surface resistance Rs(f) = Rs(f_Rs) sqrt(f/f_Rs)
  at the mode frequency. Palace combines the impedance-boundary Rs and Ls in parallel, so the
  surface resistance seen by the run is Re Z(f) = k Rs(f_Rs) x^2 / (1 + x^2), x = f / f_Rs, with
  k = 1 for zs = "s2" (Rs, Ls = Rs/w: half the copper loss at f_Rs) and k = 2 for zs = "exact"
  (2 Rs, 2 Rs/w: Z = Rs (1 + j) at f_Rs). Q_c_cu = Q_c Re Z(f) / Rs(f),
  and the magnetic-energy shares in the air outside the loops (pm_air), inside the loop bores
  (pm_bore) and in the empty line sections (pm_end), plus pe_air / pe_end for the electric energy.
Usage: python3 postprocess_s3.py run1 [run2 ...]
"""
import csv
import json
import math
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


def dofs_of(out):
    try:
        return json.loads((out / "palace.json").read_text())["Problem"]["DegreesOfFreedom"]
    except (FileNotFoundError, KeyError, ValueError):
        return ""


def eigen(run, meta, case):
    out = Path("postpro") / run
    hdr, eig = read_csv(out / "eig.csv")
    i_m, i_f, i_q = col(hdr, r"m"), col(hdr, r"Re\{f\} \(GHz\)"), col(hdr, r"Q")
    try:
        dh, dom = read_csv(out / "domain-E.csv")
    except FileNotFoundError:
        dh, dom = None, []
    g = case["geometry"]
    tand = case["tand_alumina"]
    f_rs = meta.get("rs_wall_f_Hz", 150e6) if g["kind"] == "empty" else meta.get("rs_f_Hz", 150e6)
    k_zs = 1.0 if meta.get("zs", "s2") == "s2" else 2.0
    rows = []
    for r in sorted(eig, key=lambda r: r[i_f]):
        q = r[i_q]
        if not q > 10.0:
            continue
        f = r[i_f] * 1e9
        res = dict(run=run, mesh=meta["mesh"], kind=g["kind"], N=g["N"], order=meta["order"],
                   xwall=meta["xwall"] if g["kind"] == "line" else "", walls=meta["walls"] if g["kind"] != "line" else "",
                   n_tet=case["mesh"]["n_tet"], dofs=dofs_of(out), mode=int(r[i_m]), f_MHz=f / 1e6, Q=q)
        drow = next((x for x in dom if x[0] == r[i_m]), None)
        if drow is not None and g["N"] > 0:
            p = lambda name: drow[col(dh, name)] if col(dh, name) is not None else float("nan")
            p_d = p(r"p_elec\[1\]")
            q_d = 1.0 / (p_d * tand) if p_d > 0 else float("inf")
            q_c = 1.0 / (1.0 / q - 1.0 / q_d)
            x = f / f_rs
            q_cc = q_c * (k_zs * x * x / (1 + x * x)) / math.sqrt(x)
            res.update(p_d=p_d, Q_d=q_d, Q_c=q_c, Q_c_cu=q_cc, Q_cu=1.0 / (1.0 / q_cc + 1.0 / q_d),
                       pe_air=p(r"p_elec\[2\]"), pm_air=p(r"p_mag\[2\]"), pm_bore=p(r"p_mag\[3\]"),
                       pm_alumina=p(r"p_mag\[1\]"), pe_bore=p(r"p_elec\[3\]"),
                       pe_end=p(r"p_elec\[4\]"), pm_end=p(r"p_mag\[4\]"),
                       E_elec_J=p(r"E_elec \(J\)"), E_mag_J=p(r"E_mag \(J\)"))
        if g["N"] == 0:
            x = f / f_rs
            res.update(Q_c=q, Q_c_cu=q * (k_zs * x * x / (1 + x * x)) / math.sqrt(x))
            res["Q_cu"] = res["Q_c_cu"]
        res["zs"] = meta.get("zs", "s2")
        rows.append(res)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(f"results/{run}_eig.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.9g}" if isinstance(v, float) else v) for k, v in r.items()})
    for res in rows:
        print(f"{run:<28} mode {res['mode']}: f={res['f_MHz']:9.4f} MHz  Q={res['Q']:8.1f}  "
              f"Qc={res.get('Q_c', float('nan')):8.1f}  pm_bore={res.get('pm_bore', float('nan')):.3f}")


def main():
    Path("results").mkdir(exist_ok=True)
    for run in sys.argv[1:]:
        try:
            meta = json.loads((Path("runs") / f"{run}.json").read_text())
            case = json.loads((Path("cases") / f"{meta['mesh']}.json").read_text())
            eigen(run, meta, case)
        except Exception as exc:
            print(f"{run}: FAILED ({exc!r})", file=sys.stderr)


if __name__ == "__main__":
    main()
