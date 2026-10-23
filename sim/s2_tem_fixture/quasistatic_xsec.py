#!/usr/bin/env python3
"""Quasi-static (TEM) cross-section model of the S2 measurement fixture: characteristic impedance,
transverse field uniformity, and the "fixture-embedded" permeability of a finite-width slab of
magnetic meta-atoms (demagnetization).

Why a quasi-static model is exact enough here
  * The empty fixture is all air, so its fundamental mode is TEM and its transverse fields are
    exactly the 2D electrostatic (E) and magnetostatic (H) fields of the cross-section.
  * The cross-section (~0.5 m) is small against lambda (2 m at 150 MHz), so with a slab of
    homogenized cells the mode stays quasi-TEM: L' and C' per unit length follow from 2D statics
    and the line behaves as if filled with mu_fix = L'(mu) / L'_0 and eps_fix = C'(eps) / C'_0.
    This is what a TEM-line retrieval measures ("fixture-embedded" mu_eff).

Cross-section (mm), half plane y >= 0 (mirror symmetry about y = 0, a magnetic wall for the
fundamental mode):
  * ground plane z = 0; strip 0 <= y <= w/2, h <= z <= h + t (PEC);
  * open fixture: unbounded (Dirichlet at a far boundary, 4 m away);
  * closed fixture: PEC side wall at y = W_box/2 and lid at z = H_lid (rectangular coax);
  * slab of meta-atoms (homogenized): |y| <= N_y a / 2, 0 <= z <= h (below the strip), optionally
    also h + t <= z <= h + t + a (above the strip). The loop axes are along y, so the slab is
    anisotropic: mu_yy = mu_eff(f) (from the ideal S2 unit-cell retrieval), mu_zz = 1.

Magnetostatics with the flux function A (B = curl(A x_hat)), PEC: A = const (ground 0, strip 1):
  -d/dy( (1/mu_zz) dA/dy ) - d/dz( (1/mu_yy) dA/dz ) = 0,
  1/L' = (1/mu0) * integral[ (1/mu_yy)(dA/dz)^2 + (1/mu_zz)(dA/dy)^2 ]   (x 2 for the half plane)
(bilinear, not Hermitian, so it holds for complex mu). Electrostatics is the same operator with eps.
Finite volumes on a tensor grid, materials per cell, sparse LU.

The 1D "series reluctance" picture (useful to read the results): all flux lines that link the strip
pass under it, so the slab (|y| < N_y a/2) is in SERIES with the air under the strip outside the slab
and with the return path around the strip edges:
  1/mu_fix = p / mu + (1 - p),  p = reluctance share of the slab region in the empty line,
i.e. a demagnetization factor N_d = 1 - p. The apparent pole of mu_fix sits where
mu = -p / (1 - p), i.e. at f0 / sqrt(1 - F N_d) for a Lorentz mu (F = filling factor): it moves from
the unit-cell pole f_m towards the mu = 0 frequency f_mp.

Usage
  python quasistatic_xsec.py design        # Z0, uniformity, p for the open and closed variants
  python quasistatic_xsec.py map --mu-csv <S2 mu_eff csv> --variant open4 [--fmin 130e6 --fmax 200e6]
      writes results/qs_mufix_<variant>.csv: f, ideal mu, predicted fixture-embedded mu_fix
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
ETA0 = 376.730313668
C0 = 299_792_458.0
MU0 = 4e-7 * math.pi
EPS0 = 1 / (MU0 * C0 ** 2)

# fixture (from analysis of note 23 / tem_fixture_cad.py): strip underside at h, thickness t
H, T, W = 70.0, 2.0, 344.0
A_CELL = 70.0

VARIANTS = {
    # name: dict(closed=None or (W_box, H_lid), ny=cells across, above=cells also above the strip)
    "open4": dict(closed=None, ny=4, above=False),
    "open5": dict(closed=None, ny=5, above=False),
    "open6": dict(closed=None, ny=6, above=False),
    "closed4_lid210": dict(closed=(504.0, 210.0), ny=4, above=False),
    "closed4_lid300": dict(closed=(504.0, 300.0), ny=4, above=False),
    "closed4_lid144": dict(closed=(504.0, 144.0), ny=4, above=False),
    "closed4_sym144": dict(closed=(504.0, 144.0), ny=4, above=True),
    "closed6_sym144": dict(closed=(504.0, 144.0), ny=6, above=True),
}


def axis(fine_end, d_fine, far, growth=1.12, extra=()):
    """Tensor-grid axis: uniform d_fine on [0, fine_end] (with the points in `extra` inserted),
    then geometric growth up to `far`."""
    x = list(np.arange(0.0, fine_end + 1e-9, d_fine))
    for e in extra:
        if e <= fine_end:
            x.append(e)
    x = sorted(set(round(v, 6) for v in x))
    d = d_fine
    while x[-1] < far - 1e-9:
        d *= growth
        x.append(min(far, x[-1] + d))
    return np.array(x)


class XSec:
    def __init__(self, closed=None, d=1.0, far=4000.0):
        self.closed = closed
        if closed is None:
            self.y = axis(320.0, d, far, extra=(W / 2,))
            self.z = axis(160.0, d / 2 if d > 0.5 else d, far, extra=(H, H + T))
        else:
            wb, hl = closed
            self.y = axis(wb / 2, d, wb / 2, extra=(W / 2,))
            self.z = axis(hl, d / 2 if d > 0.5 else d, hl, extra=(H, H + T))
        self.ny, self.nz = self.y.size, self.z.size
        yc = 0.5 * (self.y[1:] + self.y[:-1])
        zc = 0.5 * (self.z[1:] + self.z[:-1])
        self.Yc, self.Zc = np.meshgrid(yc, zc, indexing="ij")          # cell centres
        Y, Z = np.meshgrid(self.y, self.z, indexing="ij")
        strip = (Y <= W / 2 + 1e-9) & (Z >= H - 1e-9) & (Z <= H + T + 1e-9)
        self.val = np.full(Y.shape, np.nan)                            # Dirichlet values
        self.val[strip] = 1.0
        self.val[:, 0] = 0.0                                            # ground
        self.val[-1, :] = 0.0                                           # far / side wall
        self.val[:, -1] = 0.0                                           # far / lid
        self.free = np.isnan(self.val)
        # cells inside the strip are excluded (coefficient 0) - they are inside a conductor
        self.in_strip = (self.Yc < W / 2) & (self.Zc > H) & (self.Zc < H + T)

    def region(self, ny, above=False):
        half = ny * A_CELL / 2
        m = (self.Yc < half) & (self.Zc < H)
        if above:
            m |= (self.Yc < half) & (self.Zc > H + T) & (self.Zc < H + T + A_CELL)
        return m

    def solve(self, ky, kz):
        """ky, kz: per-cell coefficients multiplying (dA/dy)^2 and (dA/dz)^2 (complex allowed).
        Returns (A, energy) with energy = sum over the half plane of ky Ay^2 + kz Az^2 (bilinear)."""
        y, z = self.y, self.z
        dy, dz = np.diff(y), np.diff(z)
        ky = np.where(self.in_strip, 0.0, ky)
        kz = np.where(self.in_strip, 0.0, kz)
        n = self.ny * self.nz
        idx = np.arange(n).reshape(self.ny, self.nz)
        rows, cols, vals = [], [], []
        # y-edges between nodes (i, j) and (i+1, j): coupling = sum over adjacent cells of ky * (dz/2) / dy
        cy = np.zeros((self.ny - 1, self.nz), complex)
        cy[:, :-1] += ky * (dz[None, :] / 2) / dy[:, None]
        cy[:, 1:] += ky * (dz[None, :] / 2) / dy[:, None]
        cz = np.zeros((self.ny, self.nz - 1), complex)
        cz[:-1, :] += kz * (dy[:, None] / 2) / dz[None, :]
        cz[1:, :] += kz * (dy[:, None] / 2) / dz[None, :]
        a, b = idx[:-1, :].ravel(), idx[1:, :].ravel()
        c = cy.ravel()
        a2, b2 = idx[:, :-1].ravel(), idx[:, 1:].ravel()
        c2 = cz.ravel()
        I = np.concatenate([a, b, a, b, a2, b2, a2, b2])
        J = np.concatenate([a, b, b, a, a2, b2, b2, a2])
        V = np.concatenate([c, c, -c, -c, c2, c2, -c2, -c2])
        K = sp.csr_matrix((V, (I, J)), shape=(n, n))
        free = self.free.ravel()
        x = np.nan_to_num(self.val.ravel()).astype(complex)
        Kff = K[free][:, free].tocsc()
        rhs = -K[free][:, ~free] @ x[~free]
        x[free] = spla.spsolve(Kff, rhs)
        energy = x @ (K @ x)
        return x.reshape(self.ny, self.nz), energy

    def vacuum(self):
        one = np.ones(self.Yc.shape)
        return self.solve(one, one)


def z0_and_fields(xs):
    A, e = xs.vacuum()
    c_full = 2 * EPS0 * e.real * 1e0          # per metre (grid in mm cancels in the ratio dA/dx * dx)
    z0 = 1 / (C0 * c_full)
    # E_z (or H_y for the TEM mode) along y at mid-height, normalized to y = 0
    jz = np.argmin(np.abs(xs.z - H / 2))
    ez = -(A[:, jz + 1] - A[:, jz - 1]).real / (xs.z[jz + 1] - xs.z[jz - 1])
    ez /= ez[0]
    return z0, A.real, ez


def cell_flux_weights(xs, A, centres, r=21.5, half_len=22.5):
    """Mean TEM H_y (prop. to E_z = -dA/dz) over each loop bore (disc of radius r in x-z,
    centred at mid-height; projected onto z) and over the loop length in y; normalized to the
    first entry."""
    out = []
    zc = H / 2
    jz = (xs.z > zc - r) & (xs.z < zc + r)
    gz = -np.gradient(A.real, xs.z, axis=1)
    for yc in centres:
        iy = (xs.y >= yc - half_len) & (xs.y <= yc + half_len)
        zz = xs.z[jz]
        wz = 2 * np.sqrt(np.clip(r ** 2 - (zz - zc) ** 2, 0, None))    # chord length of the disc
        sub = gz[np.ix_(iy, jz)]
        out.append(float(np.mean(sub @ wz) / wz.sum()))
    out = np.array(out)
    return out / out[0]


def mufix(xs, mu, ny, above=False, base=None):
    """Fixture-embedded mu for slab permeability mu (scalar complex, mu_yy; mu_zz = 1)."""
    reg = xs.region(ny, above)
    one = np.ones(xs.Yc.shape)
    if base is None:
        base = xs.vacuum()[1]
    kz = np.where(reg, 1.0 / mu, 1.0)
    _, e = xs.solve(one, kz)
    return base / e                     # L'(mu)/L'0 = energy0 / energy(mu)  (energy prop. to 1/L')


def series_p(xs, ny, above=False):
    """Reluctance share p of the slab region (from d mu_fix / d(1/mu) at mu = 1)."""
    base = xs.vacuum()[1]
    m1 = mufix(xs, 1 / 0.999, ny, above, base)       # 1/mu = 0.999
    m2 = mufix(xs, 1 / 1.001, ny, above, base)
    # 1/mu_fix = p/mu + (1 - p) -> d(1/mu_fix)/d(1/mu) = p
    return float(((1 / m2 - 1 / m1) / (0.002)).real)


def design(args):
    out = {}
    centres = [A_CELL / 2 + k * A_CELL for k in range(3)]          # 35, 105, 175 mm
    for name, v in VARIANTS.items():
        xs = XSec(v["closed"], d=args.d)
        z0, A, ez = z0_and_fields(xs)
        w = cell_flux_weights(xs, A, centres)
        p = series_p(xs, v["ny"], v["above"])
        # where does the apparent pole move? 1/mu_fix = 0  ->  mu = -p/(1-p)
        mu_pole = -p / (1 - p)
        # Lorentz cell: F = 0.302, f_m = 149.918 MHz (S2 eigenmode, corrected copper)
        F, fm = 0.302, 149.918
        nd = 1 - p
        f_pole = fm / math.sqrt(1 - F * nd)
        iy = {f"{y:.0f}": float(np.interp(y, xs.y, ez)) for y in (35, 70, 105, 127.5, 140, 172)}
        out[name] = dict(variant=v, Z0_ohm=z0, p_series=p, N_d=nd, mu_at_apparent_pole=mu_pole,
                         f_apparent_pole_MHz_lorentz=f_pole, Ez_rel_midheight=iy,
                         loop_flux_rel_35_105_175=w.tolist(), grid=[int(xs.ny), int(xs.nz)])
        print(f"{name:16s} Z0 = {z0:6.2f} ohm  p = {p:.3f}  N_d = {nd:.3f}  pole at mu = {mu_pole:6.2f}"
              f" -> ~{f_pole:6.2f} MHz  Ez(105)/Ez(0) = {iy['105']:.3f}  loop flux 105/35 = {w[1]:.3f}")
    path = HERE / "results" / "quasistatic_design.json"
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(out, indent=1) + "\n")
    print("wrote", path)


def _map_chunk(variant, d, mus):
    v = VARIANTS[variant]
    xs = XSec(v["closed"], d=d)
    base = xs.vacuum()[1]
    return np.array([mufix(xs, m, v["ny"], v["above"], base) for m in mus])


def map_mu(args):
    d = np.genfromtxt(args.mu_csv, delimiter=",", names=True)
    f = d["f_Hz"]
    mu = d["mu_re"] + 1j * d["mu_im"]                 # physics convention (mu'' > 0)
    m = (f >= args.fmin) & (f <= args.fmax)
    f, mu = f[m], mu[m]
    if args.step:
        keep = np.zeros(f.size, bool)
        last = -1e30
        for i, fi in enumerate(f):
            fine = args.fine[0] <= fi <= args.fine[1]
            if fine or fi - last >= args.step - 1:
                keep[i] = True
                last = fi
        f, mu = f[keep], mu[keep]
    v = VARIANTS[args.variant]
    # physics-convention mu'' > 0; the static problem is convention-free (bilinear)
    if args.procs > 1:
        from multiprocessing import Pool
        chunks = np.array_split(mu, args.procs)
        with Pool(args.procs) as pool:
            out = np.concatenate(pool.starmap(_map_chunk, [(args.variant, args.d, c) for c in chunks]))
    else:
        out = _map_chunk(args.variant, args.d, mu)
    path = Path(args.out) if args.out else HERE / "results" / f"qs_mufix_{args.variant}.csv"
    np.savetxt(path, np.column_stack([f, mu.real, mu.imag, out.real, out.imag]), delimiter=",",
               header="f_Hz,mu_ideal_re,mu_ideal_im,mu_fix_re,mu_fix_im", comments="")
    print("wrote", path, f.size, "points")


def cutoff_modes(closed=(504.0, 300.0), d=2.0, n=4):
    """Cut-off frequencies (MHz) of the lowest TE and TM modes of the closed cross-section
    (rectangular coax), for the even (magnetic wall at y = 0, the family of the fundamental mode) and
    the odd family. FD on the same tensor grid: TE uses H_x with Neumann on metal, TM uses E_x with
    Dirichlet on metal; generalized eigenproblem K u = kc^2 M u (lumped M)."""
    import scipy.sparse.linalg as sla
    xs = XSec(closed, d=d)
    one = np.ones(xs.Yc.shape)
    y, z = xs.y, xs.z
    dy, dz = np.diff(y), np.diff(z)
    ky = np.where(xs.in_strip, 0.0, one)
    nn = xs.ny * xs.nz
    idx = np.arange(nn).reshape(xs.ny, xs.nz)
    cy = np.zeros((xs.ny - 1, xs.nz))
    cy[:, :-1] += ky * (dz[None, :] / 2) / dy[:, None]
    cy[:, 1:] += ky * (dz[None, :] / 2) / dy[:, None]
    cz = np.zeros((xs.ny, xs.nz - 1))
    cz[:-1, :] += ky * (dy[:, None] / 2) / dz[None, :]
    cz[1:, :] += ky * (dy[:, None] / 2) / dz[None, :]
    a, b, c = idx[:-1, :].ravel(), idx[1:, :].ravel(), cy.ravel()
    a2, b2, c2 = idx[:, :-1].ravel(), idx[:, 1:].ravel(), cz.ravel()
    I = np.concatenate([a, b, a, b, a2, b2, a2, b2])
    J = np.concatenate([a, b, b, a, a2, b2, b2, a2])
    V = np.concatenate([c, c, -c, -c, c2, c2, -c2, -c2])
    K = sp.csr_matrix((V, (I, J)), shape=(nn, nn))
    area = np.zeros((xs.ny, xs.nz))
    ca = np.where(xs.in_strip, 0.0, dy[:, None] * dz[None, :]) / 4
    for sy in (0, 1):
        for sz in (0, 1):
            area[sy:xs.ny - 1 + sy, sz:xs.nz - 1 + sz] += ca
    Y, Z = np.meshgrid(y, z, indexing="ij")
    metal = ~np.isnan(xs.val)                     # strip, ground, walls, lid
    inside_strip_only = (Y < W / 2 - 1e-9) & (Z > H + 1e-9) & (Z < H + T - 1e-9)
    out = {}
    for fam in ("even", "odd"):
        for kind in ("TE", "TM"):
            if kind == "TE":
                dirich = (Y == 0) if fam == "even" else np.zeros_like(metal)
                dirich = dirich | inside_strip_only
            else:
                dirich = metal | ((Y == 0) if fam == "odd" else False)
            keep = (~dirich.ravel()) & (area.ravel() > 0)
            Kk = K[keep][:, keep].tocsc()
            Mk = sp.diags(area.ravel()[keep])
            vals = sla.eigsh(Kk, k=n + 1, M=Mk, sigma=-1e-9, which="LM", return_eigenvectors=False)
            kc = np.sqrt(np.clip(np.sort(vals.real), 0, None)) * 1e3      # 1/m (grid in mm)
            fc = [float(C0 * k / (2 * math.pi) / 1e6) for k in kc if k > 1e-3]
            out[f"{fam}_{kind}"] = fc[:n]
    return out


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("modes")
    a = sub.add_parser("design")
    a.add_argument("--d", type=float, default=1.0, help="fine grid step (mm)")
    b = sub.add_parser("map")
    b.add_argument("--mu-csv", required=True)
    b.add_argument("--variant", default="open4", choices=sorted(VARIANTS))
    b.add_argument("--fmin", type=float, default=130e6)
    b.add_argument("--fmax", type=float, default=200e6)
    b.add_argument("--step", type=float, default=1e5, help="decimation step outside the fine window (Hz)")
    b.add_argument("--fine", type=float, nargs=2, default=[147e6, 156e6])
    b.add_argument("--d", type=float, default=1.0)
    b.add_argument("--procs", type=int, default=1)
    b.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.cmd == "modes":
        res = {"closed_504x300": cutoff_modes((504.0, 300.0))}
        print(json.dumps(res, indent=1))
        (HERE / "results" / "quasistatic_cutoffs.json").write_bytes((json.dumps(res, indent=1) + "\n").encode())
        return
    design(args) if args.cmd == "design" else map_mu(args)


if __name__ == "__main__":
    main()
