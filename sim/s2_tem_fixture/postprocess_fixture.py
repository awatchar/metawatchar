#!/usr/bin/env python3
"""Summarize Palace results of the fixture runs (runs inside the Slurm job; standard library only).

For every run name (runs/<run>.json):
  * driven       -> results/<run>_S.csv: f_Hz, S11_re, S11_im, S21_re, S21_im (+ S12, S22 when both
                    ports were excited), reference planes at the lumped ports (taper tips), 50 ohm;
                    results/<run>_probe-E.csv, _probe-B.csv, _surface-F.csv, _domain-E.csv (copies);
                    results/<run>_summary.csv: worst |S11| over 130-200 MHz, |S21| and
                    1 - |S11|^2 - |S21|^2 at 150 MHz, min |S21| and its frequency.
  * boundarymode -> results/<run>_*.csv: copies of every CSV Palace wrote.
Palace writes |S| in dB and arg(S) in degrees; they are converted to re/im here.

Usage: python3 postprocess_fixture.py run1 [run2 ...]
"""
import csv
import json
import math
import re
import shutil
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


def driven(run, out):
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
    with open(Path("results") / f"{run}_S.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["f_Hz"] + [f"S{a}{b}_{p}" for (a, b) in keys for p in ("re", "im")])
        for n, fr in enumerate(f_hz):
            w.writerow([f"{fr:.1f}"] + [f"{v:.10e}" for k in keys for v in (s[k][n].real, s[k][n].imag)])
    band = [n for n, fr in enumerate(f_hz) if 130e6 <= fr <= 200e6]
    db = lambda v: 20 * math.log10(max(abs(v), 1e-300))
    n150 = min(range(len(f_hz)), key=lambda n: abs(f_hz[n] - 150e6))
    s21 = [abs(v) for v in s[(2, 1)]]
    nmin = min(band, key=lambda n: s21[n]) if band else 0
    nworst = max(band, key=lambda n: abs(s[(1, 1)][n])) if band else 0
    row = dict(run=run, n_freq=len(f_hz), worst_S11_dB_130_200=db(s[(1, 1)][nworst]),
               f_worst_S11_Hz=f_hz[nworst], S21_dB_150=db(s[(2, 1)][n150]),
               loss_150=1 - abs(s[(1, 1)][n150]) ** 2 - abs(s[(2, 1)][n150]) ** 2,
               min_S21_dB=db(s[(2, 1)][nmin]), f_min_S21_Hz=f_hz[nmin])
    with open(Path("results") / f"{run}_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row))
        w.writeheader()
        w.writerow({k: (f"{v:.9g}" if isinstance(v, float) else v) for k, v in row.items()})
    print(f"{run:<26} worst |S11| (130-200 MHz) {row['worst_S11_dB_130_200']:7.2f} dB at "
          f"{row['f_worst_S11_Hz']/1e6:7.2f} MHz; |S21|(150) {row['S21_dB_150']:7.3f} dB; "
          f"loss(150) {row['loss_150']:.2e}; min |S21| {row['min_S21_dB']:7.2f} dB at {row['f_min_S21_Hz']/1e6:8.3f} MHz")


def main():
    Path("results").mkdir(exist_ok=True)
    for run in sys.argv[1:]:
        try:
            meta = json.loads((Path("runs") / f"{run}.json").read_text())
            out = Path("postpro") / run
            for p in sorted(out.glob("*.csv")):
                if p.name != "port-S.csv" or meta["type"] != "driven":
                    shutil.copy(p, Path("results") / f"{run}_{p.name}")
            if meta["type"] == "driven":
                driven(run, out)
        except Exception as exc:     # keep going; report the failed run
            print(f"{run}: FAILED ({exc!r})", file=sys.stderr)


if __name__ == "__main__":
    main()
