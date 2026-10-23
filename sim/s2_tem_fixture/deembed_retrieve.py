#!/usr/bin/env python3
"""2x-thru de-embedding (IEEE P370, scikit-rf IEEEP370_SE_NZC_2xThru) of the simulated fixture
S-parameters, reference-plane shift to the slab faces, and mu_eff retrieval (light work on pulled
results; needs numpy, scipy, scikit-rf >= 1.0 and the retrieval module mu_eff_retrieval.py, whose
directory is passed with --retrieval).

Chain (identical to what will be done with VNA data):
  1. 2x-thru = empty fixture (FIX-FIX), simulated on f = k * 0.5 MHz up to 300 MHz (P370 needs a
     uniform grid starting at df). It is split at its centre plane (x = 0) into FIX-1 and FIX-2.
  2. The error boxes are interpolated (magnitude / unwrapped phase) onto the slab sweep grid and
     removed: D = FIX-1^-1 . (FIX-slab-FIX) . FIX-2^-1. D is the slab with both reference planes
     at x = 0, i.e. the slab with L/2 of air line removed on each side.
  3. Reference planes moved outward by L/2 of air line on each side (the uniform section is air:
     beta = k0): S_slab = D * exp(-j k0 L) (all four entries).
  4. Optional renormalization from the 50-ohm port reference to the uniform line's own impedance
     Z_line (see --zline; the self-test shows which reference the NZC error boxes deliver).
  5. Retrieval (Smith/Chen, Kramers-Kronig branch), engineering-convention input.

Usage
  python deembed_retrieve.py selftest --retrieval <dir>
  python deembed_retrieve.py run --retrieval <dir> --thru <run> --dut <run> --nx 1 [--zline 49.5]
"""
import argparse
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
C0 = 299_792_458.0
A_CELL = 0.070


def load_net(path, symmetric=True):
    import skrf as rf
    d = np.genfromtxt(path, delimiter=",", names=True)
    f = d["f_Hz"]
    s11 = d["S11_re"] + 1j * d["S11_im"]
    s21 = d["S21_re"] + 1j * d["S21_im"]
    names = d.dtype.names
    s12 = d["S12_re"] + 1j * d["S12_im"] if "S12_re" in names else s21
    s22 = d["S22_re"] + 1j * d["S22_im"] if "S22_re" in names else s11
    s = np.stack([np.stack([s11, s12], -1), np.stack([s21, s22], -1)], -2)
    return rf.Network(frequency=rf.Frequency.from_f(f, unit="Hz"), s=s, z0=50)


def interp_polar(net, f_new):
    """Interpolate a smooth 2-port (an error box) in magnitude and unwrapped phase."""
    import skrf as rf
    f = net.f
    s = np.empty((f_new.size, 2, 2), complex)
    for i in range(2):
        for j in range(2):
            v = net.s[:, i, j]
            mag = np.interp(f_new, f, np.abs(v))
            ph = np.interp(f_new, f, np.unwrap(np.angle(v)))
            s[:, i, j] = mag * np.exp(1j * ph)
    return rf.Network(frequency=rf.Frequency.from_f(f_new, unit="Hz"), s=s, z0=net.z0[0, 0])


def error_boxes(thru, z0=50.0):
    from skrf.calibration.deembedding import IEEEP370_SE_NZC_2xThru
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dm = IEEEP370_SE_NZC_2xThru(dummy_2xthru=renorm(thru, z0), z0=z0, name="2xthru")
    return dm.s_side1, dm.s_side2


def deembed(thru, fdf, z0=50.0):
    """P370 NZC with every network first renormalized to z0 (the uniform line's impedance): the
    self-test shows that this puts the centre-plane reference at z0, whereas splitting at 50 ohm
    leaves an error of order (Z_line - 50)/100 in S."""
    s1, s2_ = error_boxes(thru, z0)
    f = fdf.f
    e1, e2 = interp_polar(s1, f), interp_polar(s2_, f)
    return e1.inv ** renorm(fdf, z0) ** e2.flipped().inv, (s1, s2_)


def shift_out(net, L):
    """Move both reference planes outward by L/2 of air line (beta = k0), exp(+j w t)."""
    import skrf as rf
    k0 = 2 * math.pi * net.f / C0
    ph = np.exp(-1j * k0 * L)[:, None, None]
    return rf.Network(frequency=net.frequency, s=net.s * ph, z0=net.z0)


def renorm(net, z):
    out = net.copy()
    out.renormalize(z)
    return out


def net_to_csv(net, path, extra_header=""):
    s = net.s
    np.savetxt(path, np.column_stack([net.f, s[:, 0, 0].real, s[:, 0, 0].imag, s[:, 1, 0].real, s[:, 1, 0].imag,
                                      s[:, 0, 1].real, s[:, 0, 1].imag, s[:, 1, 1].real, s[:, 1, 1].imag]),
               delimiter=",", header="f_Hz,S11_re,S11_im,S21_re,S21_im,S12_re,S12_im,S22_re,S22_im",
               comments=extra_header)


def retrieve_tolerant(mr, f, s11, s21, L, tol=0.02):
    """mu_eff_retrieval.retrieve, but without its passivity z-flip where Im n is only slightly
    negative. De-embedded data can be marginally active (|S21| up to ~0.01 dB above 0 dB from
    de-embedding error); the module then flips the sign of z, which flips the sign of mu and eps.
    Here the flip is undone wherever |Im n| < tol * |Re n| (the passive branch Re z > 0 is kept)."""
    r = mr.retrieve(f, s11, s21, L, convention="engineering")
    z, n = r["z"].copy(), r["n"].copy()
    flipped = z.real < 0
    k0 = 2 * math.pi * np.asarray(f) / C0
    s11p, s21p = np.conj(s11), np.conj(s21)                       # physics convention, as inside
    undo = flipped & (np.abs(n.imag) < tol * np.maximum(np.abs(n.real), 1e-3))
    if np.any(undo):
        zu = -z[undo]
        X = s21p[undo] / (1 - s11p[undo] * (zu - 1) / (zu + 1))
        n[undo] = n.real[undo] - 1j * np.log(np.abs(X)) / (k0[undo] * L)
        z[undo] = zu
    r = dict(r, z=z, n=n, mu=n * z, eps=n / z, unflipped=undo)
    return r


# ---------------- self-test: synthetic fixture, known slab ----------------
def synthetic(f, zline, slab_mu=None, slab_eps=None, L=0.0, launch_c=0.3e-12):
    """Half fixture (port -> 300 mm taper 50 -> zline in 30 steps -> uniform line to the centre)
    and, if a slab is given, the loaded fixture. Air (gamma = j k0), all ports referenced to 50."""
    import skrf as rf
    from skrf.media import DefinedGammaZ0
    freq = rf.Frequency.from_f(f, unit="Hz")
    k0 = 2 * math.pi * f / C0

    def line(d, z):
        m = DefinedGammaZ0(frequency=freq, z0_port=50, z0=z, gamma=1j * k0)
        return m.line(d, unit="m")

    m50 = DefinedGammaZ0(frequency=freq, z0_port=50, z0=50, gamma=1j * k0)
    half = m50.shunt_capacitor(launch_c)
    for i in range(30):
        half = half ** line(0.010, 50 + (zline - 50) * (i + 0.5) / 30)
    full_half = half ** line(0.210, zline)
    thru = full_half ** full_half.flipped()
    if slab_mu is None:
        return thru, None, None
    # slab as a line section: gamma = j k0 conj(n_phys), Z = conj(z_phys) * zline (engineering)
    n = np.sqrt(slab_mu * slab_eps)
    n = np.where(n.imag < 0, -n, n)
    z = np.sqrt(slab_mu / slab_eps)
    z = np.where(z.real < 0, -z, z)
    ms = DefinedGammaZ0(frequency=freq, z0_port=50, z0=np.conj(z) * zline, gamma=1j * k0 * np.conj(n))
    slab = ms.line(L, unit="m")
    side = half ** line(0.210 - L / 2, zline)
    fdf = side ** slab ** side.flipped()
    # the slab alone, referenced to the line impedance (what the retrieval should see)
    slab_ref = renorm(slab, zline)
    return thru, fdf, slab_ref


def selftest(args):
    sys.path.insert(0, args.retrieval)
    import mu_eff_retrieval as mr
    ideal = np.genfromtxt(HERE.parent / "s2_mng_slab" / "results" / "S_N1_cu2.mu_eff.csv", delimiter=",", names=True)
    f_d = ideal["f_Hz"][::10]
    mu = ideal["mu_re"][::10] + 1j * ideal["mu_im"][::10]
    eps = ideal["eps_re"][::10] + 1j * ideal["eps_im"][::10]
    f_t = np.arange(1, 601) * 0.5e6
    L = A_CELL
    out = {}
    for zl in (50.0, 49.5, 47.0):
        thru, _, _ = synthetic(f_t, zl)
        _, fdf, slab_ref = synthetic(f_d, zl, mu, eps, L)
        rec = {}
        # (a) split at 50 ohm, renormalize afterwards; (b) renormalize to Z_line first (used)
        for tag, z_split in (("split50_then_renorm", 50.0), ("renorm_first", zl)):
            D, _ = deembed(thru, fdf, z_split)
            net = renorm(shift_out(D, L), zl)
            r = mr.retrieve(net.f, net.s[:, 0, 0], net.s[:, 1, 0], L, convention="engineering")
            err = np.abs(r["mu"] - mu)
            ok = np.abs(mu) < 30
            rec[tag] = dict(max_abs_S_err=float(np.max(np.abs(net.s - slab_ref.s))),
                            median_mu_err=float(np.median(err[ok])), max_mu_err_absmu_lt30=float(np.max(err[ok])))
        out[f"zline_{zl:g}"] = rec
        print(f"Z_line = {zl:g} ohm:", json.dumps(rec))
    (RES / "deembed_selftest.json").write_text(json.dumps(out, indent=1) + "\n")
    print("wrote", RES / "deembed_selftest.json")


# ---------------- simulated fixture data ----------------
def run(args):
    sys.path.insert(0, args.retrieval)
    import mu_eff_retrieval as mr
    fdf = load_net(RES / f"{args.dut}_S.csv")
    L = args.nx * A_CELL
    z0 = args.zline or 50.0
    tag = args.tag or args.dut
    if args.thru == "none":
        # reference model with wave ports: already de-embedded by Palace to the slab faces,
        # power-normalized to the TEM mode of the uniform line (no 2x-thru step)
        S = fdf
        net_to_csv(S, RES / f"{tag}_slab.csv", "# wave ports de-embedded by Palace to the slab faces (TEM-mode reference)\n")
        r = retrieve_tolerant(mr, S.f, S.s[:, 0, 0], S.s[:, 1, 0], L)
        np.savetxt(RES / f"{tag}_slab.mu_eff.csv",
                   np.column_stack([r["f"], r["mu"].real, r["mu"].imag, r["eps"].real, r["eps"].imag,
                                    r["n"].real, r["n"].imag, r["z"].real, r["z"].imag]),
                   delimiter=",", header="f_Hz,mu_re,mu_im,eps_re,eps_im,n_re,n_im,z_re,z_im", comments="")
        print("wrote", tag, "(direct)")
        return
    thru = load_net(RES / f"{args.thru}_S.csv")
    D, (e1, e2) = deembed(thru, fdf, z0)
    S = shift_out(D, L)
    net_to_csv(D, RES / f"{tag}_deemb_centre.csv", f"# de-embedded to the centre plane x = 0 (reference {z0} ohm)\n")
    net_to_csv(S, RES / f"{tag}_slab.csv", f"# reference planes at the slab faces, reference {z0} ohm\n")
    r = retrieve_tolerant(mr, S.f, S.s[:, 0, 0], S.s[:, 1, 0], L)
    np.savetxt(RES / f"{tag}_slab.mu_eff.csv",
               np.column_stack([r["f"], r["mu"].real, r["mu"].imag, r["eps"].real, r["eps"].imag,
                                r["n"].real, r["n"].imag, r["z"].real, r["z"].imag]),
               delimiter=",", header="f_Hz,mu_re,mu_im,eps_re,eps_im,n_re,n_im,z_re,z_im", comments="")
    # error boxes on the thru grid, for the record
    net_to_csv(e1, RES / f"{args.thru}_fix1.csv", "# P370 NZC error box FIX-1 (port 1 -> centre plane)\n")
    print("wrote", tag)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("selftest")
    a.add_argument("--retrieval", required=True)
    b = sub.add_parser("run")
    b.add_argument("--retrieval", required=True)
    b.add_argument("--thru", required=True, help='2x-thru run, or "none" for wave-port reference runs')
    b.add_argument("--dut", required=True)
    b.add_argument("--nx", type=int, required=True)
    b.add_argument("--zline", type=float, default=None)
    b.add_argument("--tag", default=None)
    args = ap.parse_args()
    selftest(args) if args.cmd == "selftest" else run(args)


if __name__ == "__main__":
    main()
