#!/usr/bin/env python3
"""Compare Palace eigenmode results for the cylindrical cavity with analytic values.

Geometry and materials follow Palace's examples/cylinder cavity: a circular cylindrical
cavity of radius a = 2.74 cm and height d = 2a, filled with Teflon (eps_r = 2.08,
tan_delta = 4e-4). Walls are either PEC or copper modelled by a surface resistance
Rs = 0.0184 ohm/sq.

Analytic reference:
  * Resonant frequencies of TE_nml / TM_nml modes (Pozar, Microwave Engineering, Sec. 6.4).
  * Q_d = 1 / tan_delta (dielectric loss).
  * Q_c (wall loss) from the closed-form mode fields: Q_c = omega * W / P_c, with
    W = (eps/2) * int |E|^2 dV and P_c = (Rs/2) * surface-int |H_t|^2 dS over all walls.
    The radial integrals are evaluated numerically; for TE_011 this reproduces Pozar's
    Example 6.4 (Q_c ~ 2.94e4), and for TM_010 the textbook 1.2025*eta/(Rs*(1 + a/d)).
  * Total Q = 1 / (1/Q_d + 1/Q_c).

Usage:
  python3 compare_cavity.py postpro/cavity_pec postpro/cavity_impedance > results.md
Each directory must contain Palace's eig.csv. Directories whose name contains "impedance"
are compared against the lossy-wall Q; others against Q_d only (PEC walls).
"""
import csv
import itertools
import math
import sys
from pathlib import Path

import numpy as np
from scipy import integrate, special

C0 = 299792458.0
MU0 = 4e-7 * math.pi
EPS0 = 1.0 / (MU0 * C0**2)

A = 2.74e-2          # radius (m)
D = 2 * A            # height (m)
EPS_R = 2.08
TAN_D = 4e-4
RS = 0.0184          # ohm/sq
N_MODES = 15         # modes requested in the configs ("N": 15)


def bessel_root(n, m, derivative):
    return (special.jnp_zeros if derivative else special.jn_zeros)(n, m)[-1]


def mode(kind, n, m, l, rs=RS):
    """Return (f_Hz, Q_c) for TE_nml or TM_nml."""
    p = bessel_root(n, m, derivative=(kind == "TE"))
    kc, beta = p / A, l * math.pi / D
    k = math.hypot(kc, beta)
    eps, mu = EPS0 * EPS_R, MU0
    omega = k / math.sqrt(mu * eps)
    f = omega / (2 * math.pi)

    J = lambda r: special.jv(n, kc * r)
    Jp = lambda r: special.jvp(n, kc * r)
    R1 = integrate.quad(lambda r: J(r) ** 2 * r, 0, A, limit=200)[0]
    R2 = integrate.quad(lambda r: Jp(r) ** 2 * r, 0, A, limit=200)[0]
    R3 = 0.0 if n == 0 else integrate.quad(lambda r: n**2 * J(r) ** 2 / (kc**2 * r), 0, A, limit=200)[0]
    Ic, Is = (2 * math.pi, 0.0) if n == 0 else (math.pi, math.pi)      # int cos^2, sin^2 (n phi)
    Zc, Zs = (D, 0.0) if l == 0 else (D / 2, D / 2)                   # int cos^2, sin^2 (beta z)

    if kind == "TM":   # Ez = J cos(n phi) cos(beta z)
        E2 = R1 * Ic * Zc + (beta / kc) ** 2 * (R2 * Ic + R3 * Is) * Zs
        h = omega * eps / kc
        side = h**2 * Jp(A) ** 2 * Ic * Zc * A                        # H_phi at r = a
        caps = 2 * h**2 * (R2 * Ic + R3 * Is)                         # H_r, H_phi at z = 0, d
    else:              # Hz = J cos(n phi) sin(beta z)
        E2 = (omega * mu / kc) ** 2 * (R2 * Ic + R3 * Is) * Zs
        side = (J(A) ** 2 * Ic * Zs + (beta * n / (kc**2 * A)) ** 2 * J(A) ** 2 * Is * Zc) * A
        caps = 2 * (beta / kc) ** 2 * (R2 * Ic + R3 * Is)
    W = 0.5 * eps * E2
    P = 0.5 * rs * (side + caps)
    return f, omega * W / P


def analytic_modes(fmin=2.0e9, count=N_MODES):
    """Lowest analytic modes above fmin, with degenerate (cos/sin) pairs listed twice."""
    out = []
    for kind, n, m, l in itertools.product(("TE", "TM"), range(4), range(1, 3), range(3)):
        if kind == "TE" and l == 0:
            continue
        f, qc = mode(kind, n, m, l)
        if f < fmin:
            continue
        for _ in range(1 if n == 0 else 2):
            out.append(dict(name=f"{kind}{n}{m}{l}", f=f, qc=qc))
    out.sort(key=lambda x: x["f"])
    return out[:count]


def read_eig(path):
    rows = []
    with open(path) as fh:
        rd = csv.reader(fh)
        next(rd)
        for r in rd:
            if r:
                rows.append(dict(f=float(r[1]) * 1e9, q=float(r[3])))
    return rows


def match(sim, ref, lossy):
    """Pair sorted simulated and analytic modes; inside groups of (near-)degenerate analytic
    frequencies, choose the permutation that best matches the expected Q."""
    sim = sorted(sim, key=lambda x: x["f"])
    pairs, i = [], 0
    while i < len(ref):
        j = i
        while j + 1 < len(ref) and abs(ref[j + 1]["f"] - ref[i]["f"]) < 1e-6 * ref[i]["f"]:
            j += 1
        group_ref, group_sim = ref[i:j + 1], sim[i:j + 1]
        if lossy and len(group_ref) > 1:
            qexp = lambda r: 1 / (TAN_D + 1 / r["qc"])
            best = min(itertools.permutations(group_ref),
                       key=lambda p: sum(abs(s["q"] - qexp(r)) for s, r in zip(group_sim, p)))
            group_ref = list(best)
        pairs += list(zip(group_sim, group_ref))
        i = j + 1
    return pairs


def report(case_dir):
    lossy = "impedance" in Path(case_dir).name
    sim = read_eig(Path(case_dir) / "eig.csv")
    ref = analytic_modes(count=len(sim))
    qd = 1 / TAN_D
    print(f"\n### {Path(case_dir).name} ({'Rs = %.4f ohm/sq walls' % RS if lossy else 'PEC walls'})\n")
    if lossy:
        print("| # | Mode | f analytic (GHz) | f Palace (GHz) | rel. err | Q analytic | Q Palace | "
              "Q_c analytic | Q_c Palace | Q_c rel. err |")
        print("|---|---|---|---|---|---|---|---|---|---|")
    else:
        print("| # | Mode | f analytic (GHz) | f Palace (GHz) | rel. err | Q analytic (=Q_d) | Q Palace |")
        print("|---|---|---|---|---|---|---|")
    worst_f = worst_q = 0.0
    for k, (s, r) in enumerate(match(sim, ref, lossy), 1):
        ef = (s["f"] - r["f"]) / r["f"]
        worst_f = max(worst_f, abs(ef))
        if lossy:
            q_an = 1 / (1 / qd + 1 / r["qc"])
            qc_sim = 1 / (1 / s["q"] - 1 / qd)
            eq = (qc_sim - r["qc"]) / r["qc"]
            worst_q = max(worst_q, abs(eq))
            print(f"| {k} | {r['name']} | {r['f']/1e9:.6f} | {s['f']/1e9:.6f} | {ef:+.2e} | {q_an:.1f} | "
                  f"{s['q']:.1f} | {r['qc']:.4g} | {qc_sim:.4g} | {eq:+.2e} |")
        else:
            eq = (s["q"] - qd) / qd
            worst_q = max(worst_q, abs(eq))
            print(f"| {k} | {r['name']} | {r['f']/1e9:.6f} | {s['f']/1e9:.6f} | {ef:+.2e} | {qd:.1f} | {s['q']:.1f} |")
    label = "Q_c" if lossy else "Q"
    print(f"\nMax |rel. err|: f = {worst_f:.2e}, {label} = {worst_q:.2e}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for d in sys.argv[1:]:
        report(d)
