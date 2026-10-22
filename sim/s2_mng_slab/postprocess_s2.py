#!/usr/bin/env python3
"""Summarize Palace results of the S2 runs (runs inside the Slurm job; standard library only).

For every run name (runs/<run>.json + cases/<mesh>.json):
  * eigen  -> results/<run>_eig.csv, one row per physical mode (Q > 10): f (MHz), Q,
              p_d (alumina electric-energy participation), Q_d = 1/(p_d tan_delta),
              Q_c = (1/Q - 1/Q_d)^-1 (copper loop + lips; the PEC/PMC walls are lossless).
  * driven -> results/<run>_S.csv with f_Hz, S11_re, S11_im, S21_re, S21_im (+ S12, S22 if the
              second port was also excited), de-embedded by Palace to the slab faces (port Offset),
              and results/<run>_summary.csv: frequency and depth of min |S21|,
              max absorption 1 - |S11|^2 - |S21|^2, and the passivity check.
Palace writes |S| in dB and arg(S) in degrees; they are converted to re/im here.

Usage: python3 postprocess_s2.py run1 [run2 ...]
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
    raise KeyError(pattern)


def append_rows(path, rows):
    path = Path(path)
    old = []
    if path.exists():
        with open(path) as fh:
            old = list(csv.DictReader(fh))
    names = {r["run"] for r in rows}
    old = [r for r in old if r["run"] not in names]
    allr = old + [{k: (f"{v:.9g}" if isinstance(v, float) else v) for k, v in r.items()} for r in rows]
    keys = []
    for r in allr:
        keys += [k for k in r if k not in keys]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in allr:
            w.writerow(r)


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
    tand = case["tand_alumina"]
    rows = []
    for r in sorted(eig, key=lambda r: r[i_f]):
        q = r[i_q]
        if not q > 10.0:          # spurious, critically damped eigenvalues of the impedance PEP
            continue
        res = dict(run=run, mesh=meta["mesh"], xwall=meta["xwall"], order=meta["order"],
                   orient=case["geometry"]["orient"], s_mm=case["geometry"]["s"],
                   n_tet=case["mesh"]["n_tet"], dofs=dofs_of(out), mode=int(r[i_m]),
                   f_MHz=r[i_f] * 1e3, Q=q)
        drow = next((x for x in dom if x[0] == r[i_m]), None)
        if drow is not None:
            p_d = drow[col(dh, r"p_elec\[1\]")]
            q_d = 1.0 / (p_d * tand) if p_d > 0 else float("inf")
            res.update(p_d=p_d, Q_d=q_d, Q_c=1.0 / (1.0 / q - 1.0 / q_d),
                       E_elec_J=drow[col(dh, r"E_elec \(J\)")], E_mag_J=drow[col(dh, r"E_mag \(J\)")])
        rows.append(res)
    append_rows(f"results/{run}_eig.csv", rows)
    for res in rows:
        print(f"{run:<28} mode {res['mode']}: f={res['f_MHz']:9.4f} MHz  Q={res['Q']:8.1f}  "
              f"Qc={res.get('Q_c', float('nan')):8.1f}  p_d={res.get('p_d', float('nan')):.4f}")


def driven(run, meta, case):
    out = Path("postpro") / run
    hdr, data = read_csv(out / "port-S.csv")
    i_f = col(hdr, r"f \(GHz\)")
    sp = {}
    for i, h in enumerate(hdr):
        m = re.fullmatch(r"\|S\[(\d+)\]\[(\d+)\]\| \(dB\)", h)
        if m:
            j = col(hdr, rf"arg\(S\[{m.group(1)}\]\[{m.group(2)}\]\) \(deg\.\)")
            sp[(int(m.group(1)), int(m.group(2)))] = (i, j)
    keys = [k for k in [(1, 1), (2, 1), (1, 2), (2, 2)] if k in sp]
    f_hz = [r[i_f] * 1e9 for r in data]
    s = {}
    for k in keys:
        i, j = sp[k]
        s[k] = [10 ** (r[i] / 20) * complex(math.cos(math.radians(r[j])), math.sin(math.radians(r[j])))
                for r in data]
    path = Path("results") / f"{run}_S.csv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["f_Hz"] + [f"S{a}{b}_{p}" for (a, b) in keys for p in ("re", "im")])
        for n, fr in enumerate(f_hz):
            w.writerow([f"{fr:.1f}"] + [f"{v:.10e}" for k in keys for v in (s[k][n].real, s[k][n].imag)])
    s21 = [abs(v) for v in s[(2, 1)]]
    n_min = min(range(len(s21)), key=lambda n: s21[n])
    absorb = [1 - abs(a) ** 2 - abs(b) ** 2 for a, b in zip(s[(1, 1)], s[(2, 1)])]
    n_abs = max(range(len(absorb)), key=lambda n: absorb[n])
    row = dict(run=run, mesh=meta["mesh"], order=meta["order"], N=case["geometry"]["N"],
               orient=case["geometry"]["orient"], s_mm=case["geometry"]["s"], n_tet=case["mesh"]["n_tet"],
               dofs=dofs_of(out), n_freq=len(f_hz), f_lo_Hz=f_hz[0], f_hi_Hz=f_hz[-1],
               f_S21min_Hz=f_hz[n_min], S21min_dB=20 * math.log10(max(s21[n_min], 1e-300)),
               S11_at_S21min_dB=20 * math.log10(max(abs(s[(1, 1)][n_min]), 1e-300)),
               f_absmax_Hz=f_hz[n_abs], absorption_max=absorb[n_abs], absorption_min=min(absorb),
               S21_lo_dB=20 * math.log10(s21[0]), S21_hi_dB=20 * math.log10(s21[-1]),
               argS21_lo_deg=math.degrees(math.atan2(s[(2, 1)][0].imag, s[(2, 1)][0].real)),
               offset_mm=meta["offset_mm"])
    append_rows(f"results/{run}_summary.csv", [row])
    print(f"{run:<28} min|S21| {row['S21min_dB']:8.2f} dB at {row['f_S21min_Hz']/1e6:9.4f} MHz; "
          f"max absorption {row['absorption_max']:.4f} at {row['f_absmax_Hz']/1e6:9.4f} MHz; "
          f"min(1-|S11|^2-|S21|^2) {row['absorption_min']:.2e}")


def main():
    Path("results").mkdir(exist_ok=True)
    for run in sys.argv[1:]:
        try:
            meta = json.loads((Path("runs") / f"{run}.json").read_text())
            case = json.loads((Path("cases") / f"{meta['mesh']}.json").read_text())
            (eigen if meta["type"] == "eigen" else driven)(run, meta, case)
        except Exception as exc:     # keep going; report the failed run
            print(f"{run}: FAILED ({exc!r})", file=sys.stderr)


if __name__ == "__main__":
    main()
