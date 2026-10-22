#!/usr/bin/env python3
"""Summarize LGR-A Palace runs (standard library only, so it runs inside the Slurm job).

  python3 postprocess.py <round> <case> [<case> ...]

Eigenmode cases (name starts with "e_"): f0, loaded Q, port Q_ext (port-Q.csv), electric energy
participation of alumina and PTFE (domain-E.csv), and from them
  Q_u = (1/Q - sum 1/Q_ext,i)^-1,  beta_i = Q_u / Q_ext,i,
  |S21(f0)| = 2 sqrt(beta1 beta2) / (1 + beta1 + beta2),  Q_d = 1/(p tan_delta),  Q_c.
Driven cases ("d_"): port-S.csv -> results/sparams_<case>.csv (f_Hz, S11_re, S11_im, S21_re,
S21_im, S22_re, S22_im, S12_re, S12_im), plus a resonance fit:
  S21(f) = C / (1 + 2j Q_L (f - f0)/f0) + background; S11, S22 = A + B / (1 + 2j Q_L (f - f0)/f0)
  beta_1 = d1 / (2 - d1 - d2), beta_2 = d2 / (2 - d1 - d2) with d_i = |B_i| (circle diameters),
  Q_u = Q_L (1 + beta_1 + beta_2).
Writes results/<round>.json (one record per case) and prints a table.
"""
import cmath
import csv
import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TAND = {"alumina": 4.0e-5, "ptfe": 2.0e-4}


def read_csv(path):
    with open(path) as fh:
        rows = [r for r in csv.reader(fh) if r]
    head = [h.strip() for h in rows[0]]
    data = [[float(x) for x in r] for r in rows[1:]]
    return head, data


def col(head, pattern):
    for i, h in enumerate(head):
        if re.search(pattern, h):
            return i
    return None


def eigen(case):
    pp = HERE / "postpro" / case
    head, data = read_csv(pp / "eig.csv")
    i_f, i_q = col(head, r"Re\{f\}"), col(head, r"^Q$")
    modes = [dict(m=int(r[0]), f=r[i_f] * 1e9, Q=r[i_q]) for r in data]
    good = [md for md in modes if md["Q"] > 10 and 0.08e9 < md["f"] < 0.3e9]
    if not good:
        return dict(case=case, error="no LGR mode", modes=modes)
    md = good[0]
    rec = dict(case=case, f0_Hz=md["f"], Q=md["Q"], mode=md["m"])
    # port coupling Q
    pq = pp / "port-Q.csv"
    qext = {}
    if pq.exists():
        h2, d2 = read_csv(pq)
        row = [r for r in d2 if int(r[0]) == md["m"]][0]
        for i, hh in enumerate(h2):
            mm = re.match(r"Q_ext\[(\d+)\]", hh) or re.match(r"Q\[(\d+)\]", hh)
            if mm:
                qext[int(mm.group(1))] = row[i]
        rec["port_Q_header"] = h2
    # energy participation
    de = pp / "domain-E.csv"
    if de.exists():
        h3, d3 = read_csv(de)
        row = [r for r in d3 if int(r[0]) == md["m"]][0]
        e_tot = row[col(h3, r"^E_elec \(J\)")]
        for idx, mat in ((1, "alumina"), (2, "ptfe")):
            j = col(h3, rf"^E_elec\[{idx}\]")
            if j is not None:
                rec[f"p_{mat}"] = row[j] / e_tot
    inv_qd = sum(rec.get(f"p_{k}", 0.0) * v for k, v in TAND.items())
    inv_qu = 1.0 / rec["Q"] - sum(1.0 / q for q in qext.values() if q > 0)
    rec["Q_u"] = 1.0 / inv_qu
    rec["Q_d"] = 1.0 / inv_qd if inv_qd > 0 else None
    rec["Q_c"] = 1.0 / (inv_qu - inv_qd)
    if qext:
        rec["Q_ext"] = qext
        b = {k: rec["Q_u"] / q for k, q in qext.items()}
        rec["beta"] = b
        if 1 in b and 2 in b:
            s21 = 2 * math.sqrt(b[1] * b[2]) / (1 + b[1] + b[2])
            rec["S21_f0_dB_pred"] = 20 * math.log10(s21)
            rec["Q_L"] = rec["Q"]
    return rec


def solve_lin(rows, rhs):
    """Complex least squares for a small dense system via normal equations (Gauss elimination)."""
    n = len(rows[0])
    a = [[sum(r[i].conjugate() * r[j] for r in rows) for j in range(n)] for i in range(n)]
    b = [sum(r[i].conjugate() * y for r, y in zip(rows, rhs)) for i in range(n)]
    for i in range(n):
        p = max(range(i, n), key=lambda k: abs(a[k][i]))
        a[i], a[p], b[i], b[p] = a[p], a[i], b[p], b[i]
        for k in range(i + 1, n):
            fct = a[k][i] / a[i][i]
            for j in range(i, n):
                a[k][j] -= fct * a[i][j]
            b[k] -= fct * b[i]
    x = [0j] * n
    for i in reversed(range(n)):
        x[i] = (b[i] - sum(a[i][j] * x[j] for j in range(i + 1, n))) / a[i][i]
    return x


def fit_resonance(f, s, f0, ql):
    """Fit s(f) = A + A1 (f - f0) + B / (1 + 2j QL (f - f0)/f0); return A, B."""
    rows = [[1 + 0j, complex((fi - f0) / f0), 1 / (1 + 2j * ql * (fi - f0) / f0)] for fi in f]
    a, a1, b = solve_lin(rows, s)
    res = math.sqrt(sum(abs(y - (a + a1 * r[1] + b * r[2])) ** 2 for r, y in zip(rows, s)) / len(s))
    return a, b, res


def driven(case):
    pp = HERE / "postpro" / case
    head, data = read_csv(pp / "port-S.csv")
    f = [r[0] * 1e9 for r in data]
    S = {}
    for i, hh in enumerate(head):
        mm = re.match(r"\|S\[(\d+)\]\[(\d+)\]\| \(dB\)", hh)
        if mm:
            k = (int(mm.group(1)), int(mm.group(2)))
            ph = col(head, re.escape(f"arg(S[{k[0]}][{k[1]}])"))
            S[k] = [10 ** (r[i] / 20) * cmath.exp(1j * math.radians(r[ph])) for r in data]
    out = HERE / "results" / f"sparams_{case}.csv"
    keys = [(1, 1), (2, 1), (2, 2), (1, 2)]
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["f_Hz", "S11_re", "S11_im", "S21_re", "S21_im", "S22_re", "S22_im", "S12_re", "S12_im"])
        for n, fi in enumerate(f):
            row = [f"{fi:.1f}"]
            for k in keys:
                v = S[k][n] if k in S else complex("nan")
                row += [f"{v.real:.9e}", f"{v.imag:.9e}"]
            w.writerow(row)
    s21 = S[(2, 1)]
    mag2 = [abs(x) ** 2 for x in s21]
    n0 = max(range(len(f)), key=lambda n: mag2[n])
    if 0 < n0 < len(f) - 1:           # parabolic refinement of the peak
        y0, y1, y2 = mag2[n0 - 1], mag2[n0], mag2[n0 + 1]
        dn = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    else:
        dn = 0.0
    df = f[1] - f[0]
    f0 = f[n0] + dn * df
    half = mag2[n0] / 2

    def cross(rng):
        for n in rng:
            if mag2[n] < half:
                n2 = n + (1 if rng.step < 0 else -1)
                return f[n] + (half - mag2[n]) * (f[n2] - f[n]) / (mag2[n2] - mag2[n])
        return None
    fl, fh_ = cross(range(n0, -1, -1)), cross(range(n0, len(f)))
    ql = f0 / (fh_ - fl) if fl and fh_ else None
    rec = dict(case=case, f0_Hz=f0, S21_f0_dB=10 * math.log10(mag2[n0]), Q_L_3dB=ql, n_freq=len(f),
               f_min_Hz=f[0], f_max_Hz=f[-1], f_step_Hz=df, csv=out.name)
    if ql:
        sel = [n for n in range(len(f)) if abs(f[n] - f0) < 6 * f0 / ql]
        fs = [f[n] for n in sel]
        c21, b21, r21 = fit_resonance(fs, [s21[n] for n in sel], f0, ql)
        _, b11, r11 = fit_resonance(fs, [S[(1, 1)][n] for n in sel], f0, ql)
        _, b22, r22 = fit_resonance(fs, [S[(2, 2)][n] for n in sel], f0, ql) if (2, 2) in S else (0, 0, 0)
        d1, d2 = abs(b11), abs(b22)
        beta1, beta2 = d1 / (2 - d1 - d2), d2 / (2 - d1 - d2)
        rec.update(S21_resonant_amp_dB=20 * math.log10(abs(b21)), S21_background_dB=20 * math.log10(abs(c21) + 1e-30),
                   d11=d1, d22=d2, beta1=beta1, beta2=beta2, Q_u=ql * (1 + beta1 + beta2),
                   S11_f0_dB=20 * math.log10(abs(S[(1, 1)][n0])), S22_f0_dB=20 * math.log10(abs(S[(2, 2)][n0])) if (2, 2) in S else None,
                   fit_rms=dict(S21=r21, S11=r11, S22=r22),
                   S21_f0_dB_from_beta=20 * math.log10(2 * math.sqrt(beta1 * beta2) / (1 + beta1 + beta2)))
    pj = pp / "palace.json"
    if pj.exists():
        try:
            j = json.loads(pj.read_text())
            rec["palace_total_s"] = j["ElapsedTime"]["Durations"]["Total"]
            rec["linear_its"] = j.get("LinearSolver", {})
        except Exception:
            pass
    return rec


def main():
    rnd, names = sys.argv[1], sys.argv[2:]
    recs = []
    for n in names:
        try:
            recs.append(eigen(n) if n.startswith("e_") else driven(n))
        except Exception as exc:            # keep going; report the failure
            recs.append(dict(case=n, error=repr(exc)))
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results" / f"{rnd}.json").write_text(json.dumps(recs, indent=1) + "\n")
    for r in recs:
        keep = {k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()
                if k not in ("port_Q_header", "linear_its", "fit_rms")}
        print(json.dumps(keep))


if __name__ == "__main__":
    main()
