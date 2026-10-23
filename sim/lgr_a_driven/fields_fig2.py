#!/usr/bin/env python3
"""Field maps of the LGR-A fundamental mode from a Palace eigenmode run with saved fields
(paper P1, Fig. 2).

  python fields_fig2.py [<case>] [--step 0.1] [--zoom-step 0.025]

<case> defaults to e_fil_s1055_cu2_o2_fields (build_cases.py group cu2_fields). Reads the ParaView
output in postpro/<case>/paraview/ and needs numpy, matplotlib, pyvista/vtk (cmocean optional).
Run it where the field dump is, i.e. as part of the Slurm job; only the figures and a small JSON
summary are kept.

What it does
  1. Reads the fundamental mode (complex E and B phasors, peak amplitude) and the mesh.
  2. Integrates the stored energy from the fields, W_e = 1/2 int eps |E|^2 dV and
     W_m = 1/2 int |B|^2 / mu0 dV (equal at resonance), and the electric participation of the
     alumina and the PTFE (check against Palace's domain-E.csv).
  3. Scales the mode to P_loss = 1 W dissipated in the resonator: W = Q P / omega, using f0 and Q
     of the same run (eig.csv). Fields then are peak phasor amplitudes at 1 W.
  4. Samples three planes with a VTK probe filter on regular grids. Copper is not meshed; it is
     drawn from the analytic geometry (lgr_a_model.build_mesh, with the lip fillets). Grid points
     in air or dielectric that the probe misses (inside curved boundary elements) take the value
     of the nearest sampled point; their count is in the JSON summary (`sampling_fill`).
       z = 0 (mid-height), full section and a zoom on the gap;  y = 0 (vertical, through the gap).
  5. Writes figures/fig2a_E_midplane, fig2b_Hz_midplane, fig2c_H_vertical, fig2_fields (all
     panels) as PNG (600 dpi) and PDF, and results/fields_<case>.json.
Coordinates: can axis = z, mid-height z = 0, gap and lips on +x (lgr_a_model.py).
"""
import argparse
import csv
import json
import math
import re
import time
from pathlib import Path

import numpy as np

import lgr_a_model as m

HERE = Path(__file__).resolve().parent
EPS0, MU0 = 8.8541878128e-12, 4e-7 * math.pi
EPS_R = {1: 1.0, 2: m.EPS_AL, 3: m.EPS_PTFE, 4: m.EPS_PEEK}
METAL = "#9b9b9b"      # copper (not meshed): neutral grey, greyscale-safe
O3_REF = dict(case="e_fil_s1055_cu2_o3", f0_MHz=145.2048, Q=2151.76)   # results/cu2_eig.json


# ----------------------------------------------------------------------------- input

def read_csv(path):
    with open(path) as fh:
        rows = [[c.strip() for c in r] for r in csv.reader(fh) if r]
    head = rows[0]
    data = [[float(c) for c in r] for r in rows[1:]]
    return head, data


def col(head, pat):
    return next(i for i, h in enumerate(head) if re.search(pat, h))


def load_mode(case):
    import pyvista as pv
    pp = HERE / "postpro" / case / "paraview"
    pvds = sorted(p for p in pp.rglob("*.pvd") if "boundary" not in p.name and "boundary" not in p.parent.name)
    if not pvds:
        raise SystemExit(f"no ParaView output under {pp}")
    rd = pv.get_reader(str(pvds[0]))
    tv = list(getattr(rd, "time_values", []) or [])
    if tv:
        rd.set_active_time_value(tv[0])
    data = rd.read()
    if isinstance(data, pv.MultiBlock):
        data = data.combine()
    print(f"read {pvds[0].relative_to(HERE)} (time values {tv[:4]}): {data.n_cells} cells, "
          f"{data.n_points} points, arrays {list(data.point_data.keys())} / {list(data.cell_data.keys())}",
          flush=True)
    return data


def pick(names, *pats):
    for p in pats:
        for n in names:
            if re.fullmatch(p, n):
                return n
    return None


def complex_field(grid, base):
    names = list(grid.point_data.keys())
    re_n = pick(names, rf"{base}_real", rf"{base}_re", rf"Re\({base}\)", rf"{base}")
    im_n = pick(names, rf"{base}_imag", rf"{base}_im", rf"Im\({base}\)")
    if re_n is None:
        raise SystemExit(f"no {base} field in {names}")
    f = np.asarray(grid.point_data[re_n], dtype=float).astype(complex)
    if im_n is not None:
        f = f + 1j * np.asarray(grid.point_data[im_n], dtype=float)
    return f, re_n, im_n


def to_quadratic_tets(grid):
    """VTK Lagrange tetrahedra of degree 2 (type 71, 10 nodes) use the same node order as
    VTK_QUADRATIC_TETRA (24); recast them so every VTK filter treats them natively."""
    import pyvista as pv
    ct = np.asarray(grid.celltypes)
    if not np.any(ct == 71):
        return grid
    conn = np.asarray(grid.cell_connectivity)
    off = np.asarray(grid.cell_offsets)
    npts = np.diff(off)
    if not np.all(npts[ct == 71] == 10):
        return grid           # higher degree: leave as Lagrange cells
    ct2 = np.where(ct == 71, 24, ct).astype(np.uint8)
    cells = np.insert(conn, off[:-1], npts)   # [n, ids..., n, ids...]
    g2 = pv.UnstructuredGrid(cells, ct2, np.asarray(grid.points))
    for k in grid.point_data.keys():
        g2.point_data[k] = grid.point_data[k]
    for k in grid.cell_data.keys():
        g2.cell_data[k] = grid.cell_data[k]
    print("recast Lagrange order-2 tetrahedra as quadratic tetrahedra", flush=True)
    return g2


# ----------------------------------------------------------------------------- integrals

def cell_integrals(grid, unit_m, fields, g):
    """Integrate |field|^2 per material with the cell volumes of the corner tetrahedra and the
    quadratic-tet quadrature weights (vertices -1/20, edge nodes 1/5; linear tets 1/4)."""
    ct = np.asarray(grid.celltypes)
    conn = np.asarray(grid.cell_connectivity)
    off = np.asarray(grid.cell_offsets)
    npts = np.diff(off)
    ok = np.isin(ct, (10, 24, 71))
    if not np.all(ok):
        print(f"warning: {np.sum(~ok)} non-tetrahedral cells ignored in the integrals", flush=True)
    P = np.asarray(grid.points) * unit_m
    res = {}
    attr = None
    for k in ("attribute", "Attribute", "material", "element_attribute"):
        if k in grid.cell_data:
            attr = np.asarray(grid.cell_data[k]).astype(int)
            break
    for kind, nn in ((10, 4), (24, 10), (71, 10)):
        sel = np.where((ct == kind) & (npts == nn))[0]
        if sel.size == 0:
            continue
        ids = conn[off[sel][:, None] + np.arange(nn)[None, :]]
        a, b, c, d = (P[ids[:, i]] for i in range(4))
        vol = np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))) / 6.0
        w = np.full(nn, 0.25) if nn == 4 else np.r_[np.full(4, -0.05), np.full(6, 0.2)]
        if attr is not None:
            at = attr[sel]
        else:      # material from the cell centroid (mm)
            cen = P[ids[:, :4]].mean(axis=1) / 1e-3
            rr = np.hypot(cen[:, 0], cen[:, 1])
            at = np.ones(sel.size, int)
            at[(np.abs(cen[:, 1]) < g["d"] / 2) & (cen[:, 0] > g["r"] - 1.0) & (np.abs(cen[:, 2]) < g["z"] / 2)] = 2
            at[(np.abs(cen[:, 2]) > g["z"] / 2) & (rr > g["r"]) & (rr < g["r"] + g["t"])] = 3
        for name, f2 in fields.items():
            vals = (f2[ids] * w[None, :]).sum(axis=1) * vol
            for A in np.unique(at):
                key = (name, int(A))
                res[key] = res.get(key, 0.0) + float(vals[at == A].sum())
        res[("volume", 0)] = res.get(("volume", 0), 0.0) + float(vol.sum())
    return res, attr is not None


# ----------------------------------------------------------------------------- sampling

def sample_plane(grid, axes, fixed, lo, hi, step, names, tol):
    """Probe the mesh on a regular grid in a coordinate plane. axes = (i, j) of the in-plane
    coordinates, fixed = (k, value). Grid points are cell-centred (offset by step/2)."""
    import pyvista as pv
    u = np.arange(lo[0] + step / 2, hi[0], step)
    v = np.arange(lo[1] + step / 2, hi[1], step)
    U, V = np.meshgrid(u, v)
    pts = np.zeros((U.size, 3))
    pts[:, axes[0]] = U.ravel()
    pts[:, axes[1]] = V.ravel()
    pts[:, fixed[0]] = fixed[1]
    t0 = time.time()
    out = pv.PolyData(pts).sample(grid, tolerance=tol)
    mask = np.asarray(out.point_data["vtkValidPointMask"]).astype(bool).reshape(U.shape)
    vals = {n: np.asarray(out.point_data[n]).reshape(U.shape + (-1,)) for n in names}
    print(f"sampled {U.size} points in {time.time() - t0:.1f} s ({mask.mean() * 100:.1f} % inside the mesh)",
          flush=True)
    return u, v, mask, vals


# ----------------------------------------------------------------------------- plotting

def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 7, "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"], "mathtext.fontset": "dejavusans",
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600, "figure.dpi": 150})
    return plt


def cmap_div():
    try:
        import cmocean
        return cmocean.cm.balance
    except Exception:
        import matplotlib
        return matplotlib.colormaps["RdBu_r"]


def metal_layer(ax, u, v, mask, inside):
    """Grey where the plane is inside the can but not meshed (copper), white outside the can."""
    import matplotlib.colors as mc
    rgba = np.zeros(mask.shape + (4,))
    rgba[(~mask) & inside] = mc.to_rgba(METAL)
    du, dv = u[1] - u[0], v[1] - v[0]
    ax.imshow(rgba, origin="lower", extent=(u[0] - du / 2, u[-1] + du / 2, v[0] - dv / 2, v[-1] + dv / 2),
              interpolation="nearest", zorder=2)


def field_image(ax, u, v, data, **kw):
    du, dv = u[1] - u[0], v[1] - v[0]
    return ax.imshow(data, origin="lower", extent=(u[0] - du / 2, u[-1] + du / 2, v[0] - dv / 2, v[-1] + dv / 2),
                     interpolation="nearest", zorder=1, **kw)


def label(ax, text, xy, xytext, color="k", **kw):
    bb = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8)
    ax.annotate(text, xy=xy, xytext=xytext, fontsize=6.5, color=color, bbox=bb, zorder=5,
                ha=kw.pop("ha", "center"), va=kw.pop("va", "center"),
                arrowprops=dict(arrowstyle="-", lw=0.5, color="k", shrinkA=0, shrinkB=0) if xy != xytext else None,
                **kw)


def can_outline(ax, R, plane):
    import matplotlib.patches as mp
    if plane == "xy":
        ax.add_patch(mp.Circle((0, 0), R, fill=False, lw=0.6, color="k", zorder=4))
    else:
        pass


def colorbar(fig, ax, im, text, extend="neither"):
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, extend=extend)
    cb.set_label(text)
    cb.ax.tick_params(width=0.5, length=2)
    return cb


# ----------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", nargs="?", default="e_fil_s1055_cu2_o2_fields")
    ap.add_argument("--step", type=float, default=0.1, help="sampling step of the full planes (mm)")
    ap.add_argument("--zoom-step", type=float, default=0.025, help="sampling step of the gap zoom (mm)")
    ap.add_argument("--power", type=float, default=1.0, help="dissipated power for the scaling (W)")
    a = ap.parse_args()
    t_start = time.time()
    case = a.case
    meta = json.loads((HERE / "cases" / f"{case}.json").read_text())
    g = meta["geometry"]
    R, H, r, t, z, s, d, c = (g[k] for k in ("R", "H", "r", "t", "z", "s", "d", "c"))
    x_lip0, x_tip, t_lip = r + t / 2, r + t + s, g["t_lip"]

    # eigenvalue of this run
    head, rows = read_csv(HERE / "postpro" / case / "eig.csv")
    i_f, i_q = col(head, r"Re\{f\}"), col(head, r"^Q$")
    mode1 = min((rw for rw in rows if rw[i_q] > 10), key=lambda rw: rw[i_f])
    f0, Q = mode1[i_f] * 1e9, mode1[i_q]
    w0 = 2 * math.pi * f0
    print(f"f0 = {f0 / 1e6:.4f} MHz, Q = {Q:.1f}", flush=True)

    grid = to_quadratic_tets(load_mode(case))
    E, e_re, e_im = complex_field(grid, "E")
    B, b_re, b_im = complex_field(grid, "B")
    rmax = float(np.max(np.hypot(grid.points[:, 0], grid.points[:, 1])))
    unit_mm = R / rmax                    # file length unit in mm (mesh in mm -> 1)
    unit_mm = 1.0 if abs(unit_mm - 1) < 0.01 else unit_mm
    print(f"max radius in file {rmax:.4g} -> length unit {unit_mm:.6g} mm", flush=True)

    # energies (peak phasors): W_e = 1/2 int eps |E|^2 = 2 <W_e>; at resonance W = <W_e> + <W_m>
    E2 = np.sum(np.abs(E) ** 2, axis=1)
    B2 = np.sum(np.abs(B) ** 2, axis=1)
    ints, have_attr = cell_integrals(grid, unit_mm * 1e-3, dict(E2=E2, B2=B2), g)
    We_by = {A: 0.5 * EPS0 * EPS_R.get(A, 1.0) * v for (n, A), v in ints.items() if n == "E2"}
    We = sum(We_by.values())
    Wm = sum(0.5 * v / MU0 for (n, A), v in ints.items() if n == "B2")
    p_al, p_ptfe = We_by.get(2, 0.0) / We, We_by.get(3, 0.0) / We
    vol_can = math.pi * (R * 1e-3) ** 2 * H * 1e-3
    print(f"W_e = {We:.4e}, W_m = {Wm:.4e} (ratio {We / Wm:.4f}); p_alumina {p_al:.4f}, p_PTFE {p_ptfe:.5f}; "
          f"meshed volume {ints[('volume', 0)]:.4e} m^3 vs can {vol_can:.4e}", flush=True)
    pal_E = None
    try:
        hd, dr = read_csv(HERE / "postpro" / case / "domain-E.csv")
        r1 = dr[0]
        pal_E = dict(E_elec=r1[col(hd, r"^E_elec \(J\)")], E_mag=r1[col(hd, r"^E_mag \(J\)")],
                     p_alumina=r1[col(hd, r"^p_elec\[1\]")], p_ptfe=r1[col(hd, r"^p_elec\[2\]")])
    except Exception as ex:      # noqa: BLE001
        print("domain-E.csv not read:", ex)

    # scale to P_loss: W = Q P / omega
    W1 = Q * a.power / w0
    kE, kB = math.sqrt(W1 / We), math.sqrt(W1 / Wm)

    # phase: make H real (in-phase part) with positive H_z in the bore; E separately real
    def align(F):
        s2 = np.sum(F * F)
        return np.exp(-0.5j * np.angle(s2))
    pB = align(B[:, 2])
    pE = align(E[:, 1])
    names = ["Er", "Ei", "Hr", "Hi"]
    H_ = B * pB * kB / MU0
    E_ = E * pE * kE
    grid.point_data.clear()
    grid.point_data["Er"], grid.point_data["Ei"] = E_.real, E_.imag
    grid.point_data["Hr"], grid.point_data["Hi"] = H_.real, H_.imag
    sc = 1.0 / unit_mm                     # mm -> file units for the probe points
    del E, B, E_, H_, E2, B2

    fil = float(meta.get("mesh_kw", {}).get("fillet", 0.0))
    fill_stats = {}

    def metal_xy(X, Y):
        """Copper in the mid-height plane (analytic, as lgr_a_model.build_mesh): loop wall with
        the slit, and the two lips with the R `fil` fillets on their four tip edges."""
        rho = np.hypot(X, Y)
        wall = (rho > r) & (rho < r + t) & ~((np.abs(Y) < d / 2) & (X > 0))
        ay = np.abs(Y)
        lips = (X > x_lip0) & (X < x_tip) & (ay > d / 2) & (ay < d / 2 + t_lip)
        if fil > 0:
            for yc in (d / 2 + fil, d / 2 + t_lip - fil):     # fillet centres (|y|), x = x_tip - fil
                corner = (X > x_tip - fil) & (np.abs(ay - yc) < fil) & (
                    (ay - yc) * (1 if yc < d / 2 + t_lip / 2 else -1) < 0)
                lips &= ~(corner & (np.hypot(X - (x_tip - fil), ay - yc) > fil))
        return wall | lips

    def metal_xz(X, Z):
        """Copper in the plane y = 0: the loop wall on the -x side (the +x side is the slit)."""
        return (X < -r) & (X > -r - t) & (np.abs(Z) < z / 2)

    def plane(tag, axes, fixed, lo, hi, step, inside, metal):
        from scipy import ndimage
        u, v, mask, vals = sample_plane(grid, axes, (fixed[0], fixed[1] * sc),
                                        [x * sc for x in lo], [x * sc for x in hi], step * sc, names,
                                        tol=1e-4 * sc)
        u, v = u / sc, v / sc
        UU, VV = np.meshgrid(u, v)
        field = inside(UU, VV) & ~metal(UU, VV)     # air and dielectrics (analytic geometry)
        good = mask & field
        # probe misses inside curved boundary elements: take the nearest sampled value
        miss = field & ~good
        idx = ndimage.distance_transform_edt(~good, return_distances=False, return_indices=True)
        for n_ in names:
            vals[n_][miss] = vals[n_][idx[0][miss], idx[1][miss]]
        fill_stats[tag] = dict(points_in_field_region=int(field.sum()), filled_from_neighbour=int(miss.sum()),
                               probe_hits_in_analytic_metal=int((mask & ~field & inside(UU, VV)).sum()))
        Ec = vals["Er"] + 1j * vals["Ei"]
        Hc = vals["Hr"] + 1j * vals["Hi"]
        Ec[~field] = np.nan
        Hc[~field] = np.nan
        return u, v, field, Ec, Hc

    pad = 0.0
    in_can_xy = lambda X, Y: np.hypot(X, Y) < R
    in_can_xz = lambda X, Z: (np.abs(X) < R) & (np.abs(Z) < H / 2)
    xy = plane("z0_full", (0, 1), (2, 0.0), (-R - pad, -R - pad), (R + pad, R + pad), a.step, in_can_xy, metal_xy)
    zm_lo, zm_hi = (r - 3.0, -6.0), (x_tip + 4.0, 6.0)
    zm = plane("z0_gap_zoom", (0, 1), (2, 0.0), zm_lo, zm_hi, a.zoom_step, in_can_xy, metal_xy)
    xz = plane("y0_vertical", (0, 2), (1, 0.0), (-R - pad, -H / 2 - pad), (R + pad, H / 2 + pad), a.step,
               in_can_xz, metal_xz)
    print("sampling:", json.dumps(fill_stats), flush=True)

    # sign: H_z > 0 in the bore
    u, v, mask, Ec, Hc = xy
    UU, VV = np.meshgrid(u, v)
    rr = np.hypot(UU, VV)
    bore = mask & (rr < r)
    sgn = np.sign(np.nanmean(Hc[..., 2][bore].real))
    for P_ in (xy, zm, xz):
        P_[4][...] *= sgn

    # numbers for the summary
    Hz = xy[4][..., 2].real
    dA = (u[1] - u[0]) * (v[1] - v[0]) * 1e-6
    ann = mask & (rr > r + t + 0.2) & (rr < R - 0.05)
    phi_bore = float(np.nansum(MU0 * Hz[bore]) * dA)
    phi_ret = float(np.nansum(MU0 * Hz[mask & (rr > r + t)]) * dA)
    uz, vz, mz, Ez, Hz_ = zm
    Uz, Vz = np.meshgrid(uz, vz)
    in_al = mz & (np.abs(Vz) < d / 2) & (Uz > r) & (Uz < x_tip)
    Emag_z = np.linalg.norm(np.nan_to_num(Ez), axis=-1)
    xm = 0.5 * (x_lip0 + x_tip)
    colm = np.argmin(np.abs(uz - xm))
    ey = Ez[:, colm, 1].real
    sl = np.abs(vz) < d / 2
    v_gap = float(np.abs(np.nansum(ey[sl]) * (vz[1] - vz[0]) * 1e-3))
    Hmag_xy = np.linalg.norm(np.nan_to_num(np.abs(xy[4])), axis=-1)
    i0 = np.argmin(np.abs(v)); j0 = np.argmin(np.abs(u))
    summary = dict(
        case=case, order=2, f0_MHz=f0 / 1e6, Q=Q,
        order3_reference=dict(O3_REF, df0_ppm=(f0 / 1e6 / O3_REF["f0_MHz"] - 1) * 1e6,
                              dQ_percent=(Q / O3_REF["Q"] - 1) * 100),
        field_arrays=dict(E=[e_re, e_im], B=[b_re, b_im]), material_attribute_in_file=have_attr,
        energy_check=dict(W_e_over_W_m=We / Wm, p_alumina=p_al, p_ptfe=p_ptfe,
                          meshed_volume_m3=ints[("volume", 0)], palace_domain_E=pal_E),
        scaling=dict(P_loss_W=a.power, W_stored_J=W1,
                     convention="peak phasor amplitudes; W = 1/2 int eps|E|^2 = 1/2 int |B|^2/mu0 = Q P / omega"),
        at_P_loss=dict(
            E_max_alumina_kV_per_m=float(np.nanmax(Emag_z[in_al])) / 1e3,
            E_alumina_mid_lip_kV_per_m=float(np.nanmedian(Emag_z[:, colm][sl & mz[:, colm]])) / 1e3,
            E_max_air_zoom_kV_per_m=float(np.nanmax(np.where(mz & ~in_al, Emag_z, np.nan))) / 1e3,
            E_max_air_zoom_at_xy_mm=[float(Uz.flat[k]), float(Vz.flat[k])] if (
                k := int(np.nanargmax(np.where(mz & ~in_al, Emag_z, np.nan)))) >= 0 else None,
            E_max_alumina_at_xy_mm=[float(Uz.flat[k2]), float(Vz.flat[k2])] if (
                k2 := int(np.nanargmax(np.where(in_al, Emag_z, np.nan)))) >= 0 else None,
            V_gap_mid_lip_V=v_gap, x_mid_lip_mm=float(uz[colm]),
            H_bore_centre_A_per_m=float(Hmag_xy[i0, j0]),
            H_z_bore_mean_A_per_m=float(np.nanmean(Hz[bore])),
            H_z_annulus_mean_A_per_m=float(np.nanmean(Hz[ann])),
            flux_bore_Wb=phi_bore, flux_outside_loop_Wb=phi_ret,
            flux_return_ratio=-phi_ret / phi_bore,
            sampled_area_bore_mm2=float(bore.sum() * dA * 1e6),
            sampled_area_outside_loop_mm2=float((mask & (rr > r + t)).sum() * dA * 1e6),
            geometric_area_bore_mm2=math.pi * r ** 2,
            geometric_area_outside_loop_mm2=math.pi * (R ** 2 - (r + t) ** 2) - 2 * t_lip * (x_tip - (r + t))),
        sampling_mm=dict(full=a.step, zoom=a.zoom_step), sampling_fill=fill_stats)
    print(json.dumps(summary["at_P_loss"], indent=1), flush=True)

    # ------------------------------------------------------------------ figures
    plt = style()
    import matplotlib.colors as mc
    from matplotlib import colormaps
    figdir = HERE / "figures"
    figdir.mkdir(exist_ok=True)
    cm_E, cm_H, cm_D = colormaps["magma"], colormaps["viridis"], cmap_div()
    for cmx in (cm_E, cm_H):
        cmx.set_bad(alpha=0)
    pw = f"{a.power:g} W"

    def panel_E_full(fig, ax):
        u, v, mask, Ec, Hc = xy
        Em = np.linalg.norm(np.abs(Ec), axis=-1) / 1e3
        top = float(np.nanpercentile(Em, 99.95))
        im = field_image(ax, u, v, Em, cmap=cm_E, norm=mc.LogNorm(vmin=top * 1e-3, vmax=top))
        UU, VV = np.meshgrid(u, v)
        metal_layer(ax, u, v, mask, np.hypot(UU, VV) < R)
        can_outline(ax, R, "xy")
        ax.add_patch(plt.Rectangle((zm_lo[0], zm_lo[1]), zm_hi[0] - zm_lo[0], zm_hi[1] - zm_lo[1],
                                   fill=False, lw=0.6, ls="--", color="w", zorder=4))
        label(ax, "bore", (0, 0), (0, 0))
        label(ax, "Cu loop\n(t = 3 mm)", (-(r + t / 2) * 0.707, (r + t / 2) * 0.707), (-30, 30))
        label(ax, "Cu can", (-R * 0.94, -R * 0.34), (-34, -40))
        label(ax, "air", (-0.5 * (r + t + R), 0), (-0.5 * (r + t + R), 0))
        label(ax, "lips +\nalumina", (x_tip - 2, 3.5), (30, 22))
        ax.set_xlim(-R - pad, R + pad); ax.set_ylim(-R - pad, R + pad)
        ax.set_aspect("equal"); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
        colorbar(fig, ax, im, f"|E| (kV/m), {pw} dissipated", extend="max")
        ax.set_title("|E|, mid-height plane z = 0 (log scale)")

    def panel_E_zoom(fig, ax):
        u, v, mask, Ec, Hc = zm
        Em = np.linalg.norm(np.abs(Ec), axis=-1) / 1e3
        UU, VV = np.meshgrid(u, v)
        gap = mask & (np.abs(VV) < d / 2) & (UU > x_lip0) & (UU < x_tip - 1.0)
        im = field_image(ax, u, v, Em, cmap=cm_E, vmin=0, vmax=2.0 * float(np.nanmedian(Em[gap])))
        UU, VV = np.meshgrid(u, v)
        metal_layer(ax, u, v, mask, np.hypot(UU, VV) < R)
        label(ax, f"alumina, {d:g} mm", (0.5 * (x_lip0 + x_tip), 0), (0.5 * (x_lip0 + x_tip), -4.8))
        label(ax, "Cu lip", (x_lip0 + 7, d / 2 + t_lip / 2), (x_lip0 + 7, 4.8))
        label(ax, "Cu lip", (x_lip0 + 3, -d / 2 - t_lip / 2), (x_lip0 + 1, -2.8))
        label(ax, "fringing\n(air)", (x_tip + 0.6, 0.9), (x_tip + 2.6, 3.6))
        label(ax, "loop wall", (r + t / 2, 4.5), (r + t / 2, 4.5))
        label(ax, "bore", (r - 1.6, -4.5), (r - 1.6, -4.5))
        ax.set_xlim(zm_lo[0], zm_hi[0]); ax.set_ylim(zm_lo[1], zm_hi[1])
        ax.set_aspect("equal"); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
        colorbar(fig, ax, im, "|E| (kV/m)", extend="max")
        ax.set_title("|E| at the gap, z = 0 (linear scale)")

    def panel_Hz(fig, ax):
        u, v, mask, Ec, Hc = xy
        Hz = Hc[..., 2].real
        lim = float(np.nanpercentile(np.abs(Hz), 99.9))
        im = field_image(ax, u, v, Hz, cmap=cm_D, vmin=-lim, vmax=lim)
        UU, VV = np.meshgrid(u, v)
        metal_layer(ax, u, v, mask, np.hypot(UU, VV) < R)
        can_outline(ax, R, "xy")
        label(ax, "bore flux\n($H_z$ > 0)", (0, 0), (0, 0))
        label(ax, "return flux\n($H_z$ < 0)", (-0.5 * (r + t + R), 0), (-0.5 * (r + t + R), 0))
        label(ax, "return flux", (0, -0.5 * (r + t + R)), (0, -0.5 * (r + t + R)))
        label(ax, "Cu loop", ((r + t / 2) * 0.5, -(r + t / 2) * 0.866), (12, -16))
        label(ax, "gap", (x_tip - 1, 0), (30, -22))
        ax.plot([-R, 0], [0, 0], ls=":", lw=0.7, color="k", zorder=4)
        ax.set_xlim(-R - pad, R + pad); ax.set_ylim(-R - pad, R + pad)
        ax.set_aspect("equal"); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
        colorbar(fig, ax, im, f"$H_z$ (A/m), {pw} dissipated")
        ax.set_title("Axial magnetic field $H_z$, z = 0")

    def panel_Hz_line(fig, ax):
        u, v, mask, Ec, Hc = xy
        i = np.argmin(np.abs(v))
        neg = u < 0
        rad = -u[neg]
        hz = Hc[i, neg, 2].real
        o = np.argsort(rad)
        ax.plot(rad[o], hz[o], lw=1.0, color="k")
        ax.axhline(0, lw=0.5, color="0.5")
        ax.axvspan(r, r + t, color=METAL, alpha=0.6, lw=0)
        ax.axvline(R, lw=0.6, color="k")
        ok = np.isfinite(hz)
        hb = float(np.nanmean(hz[ok & (rad < r)]))
        ha_ = float(np.nanmean(hz[ok & (rad > r + t)]))
        hi_, lo_ = float(np.nanmax(hz)), float(np.nanmin(hz))
        span = hi_ - lo_
        ax.set_ylim(lo_ - 0.3 * span, hi_ + 0.3 * span)     # room above the bore and below the annulus
        ax.text(r / 2, hb + 0.08 * span, "bore", ha="center", va="bottom", fontsize=6.5)
        ax.text(r + t / 2, 0.5, "Cu\nloop", transform=ax.get_xaxis_transform(), ha="center", va="center", fontsize=6.5)
        ax.text(0.5 * (r + t + R), ha_ - 0.08 * span, "annulus", ha="center", va="top", fontsize=6.5)
        ax.set_xlim(0, R + 1)
        ax.set_xlabel("radius along φ = 180° (mm)")
        ax.set_ylabel("$H_z$ (A/m)")
        ax.set_title("$H_z$ along the dotted line")

    def panel_H_vertical(fig, ax):
        u, v, mask, Ec, Hc = xz
        Hm = np.linalg.norm(np.abs(Hc), axis=-1)
        im = field_image(ax, u, v, Hm, cmap=cm_H, vmin=0, vmax=float(np.nanpercentile(Hm, 99.9)))
        UU, VV = np.meshgrid(u, v)
        inside = (np.abs(UU) < R) & (np.abs(VV) < H / 2)
        metal_layer(ax, u, v, mask, inside)
        for x0_, z0_, w_, h_ in ((r, -z / 2, x_tip - r, z), (-r - t, z / 2, t, c), (-r - t, -H / 2, t, c),
                                 (r, z / 2, t, c), (r, -H / 2, t, c)):
            ax.add_patch(plt.Rectangle((x0_, z0_), w_, h_, fill=False, lw=0.6, ls="--", color="w", zorder=4))
        # in-plane direction of the (real, in-phase) H field
        k = max(1, int(round(3.0 / (u[1] - u[0]))))
        Hx, Hz = Hc[::k, ::k, 0].real, Hc[::k, ::k, 2].real
        nrm = np.hypot(Hx, Hz)
        ok = np.isfinite(nrm) & (nrm > 0.03 * np.nanmax(nrm))
        Uq, Vq = UU[::k, ::k], VV[::k, ::k]
        ax.quiver(Uq[ok], Vq[ok], (Hx / nrm)[ok], (Hz / nrm)[ok], color="w", scale=55, width=0.0022,
                  headwidth=3.5, headlength=4, zorder=3, pivot="mid")
        label(ax, "bore", (0, 0), (0, 0))
        label(ax, "Cu loop", (-(r + t / 2), 0), (-(r + t / 2), 0))
        label(ax, "alumina slab\n(gap plane)", (0.5 * (r + x_tip), 10), (0.5 * (r + x_tip), 10))
        label(ax, "PTFE spacer", (-(r + t / 2), z / 2 + c / 2), (-(r + t / 2) - 9, z / 2 + c / 2 + 3))
        label(ax, "PTFE spacer", ((r + t / 2), -z / 2 - c / 2), ((r + t / 2) + 9, -z / 2 - c / 2 - 3))
        label(ax, "return\nflux", (-0.5 * (r + t + R), -12), (-0.5 * (r + t + R), -12))
        label(ax, "Cu can", (R, 30), (R - 4, 38))
        ax.set_xlim(-R - pad, R + pad); ax.set_ylim(-H / 2 - pad, H / 2 + pad)
        ax.set_aspect("equal"); ax.set_xlabel("x (mm)"); ax.set_ylabel("z (mm)")
        colorbar(fig, ax, im, f"|H| (A/m), {pw} dissipated")
        ax.set_title("|H| and field direction, vertical plane y = 0")

    def save(fig, stem):
        for ext in ("png", "pdf"):
            fig.savefig(figdir / f"{stem}.{ext}", bbox_inches="tight", pad_inches=0.02)
        plt.close(fig)
        print("wrote", f"figures/{stem}.png|pdf", flush=True)

    fig, axs = plt.subplots(1, 2, figsize=(7.16, 3.1), gridspec_kw=dict(width_ratios=[1, 1.05]))
    panel_E_full(fig, axs[0]); panel_E_zoom(fig, axs[1]); fig.tight_layout(); save(fig, "fig2a_E_midplane")
    fig, axs = plt.subplots(1, 2, figsize=(7.16, 3.1), gridspec_kw=dict(width_ratios=[1, 0.9]))
    panel_Hz(fig, axs[0]); panel_Hz_line(fig, axs[1]); fig.tight_layout(); save(fig, "fig2b_Hz_midplane")
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.6))
    panel_H_vertical(fig, ax); fig.tight_layout(); save(fig, "fig2c_H_vertical")

    fig = plt.figure(figsize=(7.16, 8.3))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.45], width_ratios=[1, 1])
    panels = [(panel_E_full, gs[0, 0], "(a)"), (panel_E_zoom, gs[0, 1], "(b)"), (panel_Hz, gs[1, 0], "(c)"),
              (panel_H_vertical, gs[1, 1], "(d)"), (panel_Hz_line, gs[2, :], "(e)")]
    for fn, gsi, tag in panels:
        ax = fig.add_subplot(gsi)
        fn(fig, ax)
        ttl = ax.get_title()
        ax.set_title("")
        ax.set_title(f"{tag} {ttl}", loc="left")
    fig.tight_layout()
    save(fig, "fig2_fields")

    summary["figures"] = sorted(p.name for p in figdir.glob("fig2*"))
    summary["wall_s"] = round(time.time() - t_start, 1)
    (HERE / "results" / f"fields_{case}.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("wrote", f"results/fields_{case}.json", flush=True)


if __name__ == "__main__":
    main()
