#!/usr/bin/env python3
"""Geometry, mesh and Palace configuration for the shielded loop-gap resonator "LGR-A".

LGR-A (all lengths in mm in the mesh; Palace L0 = 1e-3):
  * Shield: closed copper can, inner radius R, inner height H = 2R, so V = 2*pi*R^3.
  * Loop: coaxial copper tube, inner radius r = a*R, wall t = 3 mm, axial length
    z = H - 2c with end clearance c = 0.25*R (centred), so z = 1.5*R.
  * Slit: one axial slit through the loop wall at +x, width d (default 0.635 mm), full length z.
  * Lips: two flat copper plates, thickness 3 mm, extend radially outward from the loop's outer
    surface by s, full length z, separated by d (they continue the slit faces).
  * Dielectric: alumina (eps_r 9.8, tan_delta 4e-5) fills the d-wide gap through the wall and
    between the lips (radial extent t + s, length z).
  * Conductors (shield, loop, lips) are not meshed: they are impedance boundaries.
    Palace combines the Rs, Ls (and Cs) of an "Impedance" boundary IN PARALLEL,
    1/Z = 1/Rs + 1/(j w Ls), so the good-conductor impedance Z = Rs (1 + j) needs
    Rs_param = 2 Rs, Ls_param = 2 Rs / w (--zs exact, default since 2026-10-23; case names
    ending in "_cu2"). The sweep runs before that used Rs, Ls = Rs / w (--zs s2), which is
    Z = Rs (1 + j) / 2: half the copper loss, so their Q_c and Q are about 2x too high
    (frequencies and fields are unaffected).

Mesh attributes:
  volumes  1 = air (inside the can, outside the copper), 2 = alumina gap
  surfaces 1 = shield inner walls, 2 = loop + lip surfaces

Usage (local, gmsh >= 4.11 with OpenCASCADE):
  python lgr_geometry.py --name v025_a05 --V 0.25 --a 0.5 --s 4.0 [--d 0.635] [--h 1.0]
      [--order 2] [--outdir .]
writes mesh/<name>.msh (Gmsh 2.2 binary, 2nd-order curved tets), configs/<name>.json (Palace
config from palace_template.json) and cases/<name>.json (geometry + probe metadata).
"""
import argparse
import copy
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

F0 = 145.5e6                     # design frequency (Hz)
SIGMA_CU = 5.8e7                 # S/m
MU0 = 4e-7 * math.pi
RS = math.sqrt(math.pi * F0 * MU0 / SIGMA_CU)   # 3.147e-3 ohm/sq at 145.5 MHz
LS = RS / (2 * math.pi * F0)                    # surface inductance of a good conductor, H/sq
T_WALL = 3.0                     # loop wall (mm)
T_LIP = 3.0                      # lip plate thickness (mm)
CLEAR_MIN = 3.0                  # required clearance R - (r + t + s) (mm)
EPS_AL, TAND_AL = 9.8, 4.0e-5


def dims(V_l, a, s, d):
    R = (V_l * 1e-3 / (2 * math.pi)) ** (1 / 3) * 1e3   # mm
    H = 2 * R
    r = a * R
    c = 0.25 * R
    z = H - 2 * c
    return dict(V_l=V_l, a=a, R=R, H=H, r=r, t=T_WALL, z=z, c=c, s=s, d=d, t_lip=T_LIP,
                clearance=R - (r + T_WALL + s))


def probes(g):
    """Probe points (mm). Edges are probed at 0.1 mm from the metal/dielectric corner."""
    r, t, s, d, z, tl = g["r"], g["t"], g["s"], g["d"], g["z"], g["t_lip"]
    ro = r + t + s                       # radial position of the lip tips / slab outer face
    e = 0.1
    p = [
        ("gap_wall_mid", [r + t / 2, 0.0, 0.0]),             # inside alumina, in the wall slit
        ("gap_lip_mid", [r + t + max(s, 0) / 2, 0.0, 0.0]),  # inside alumina, between the lips
        ("air_slab_face", [ro + e, 0.0, 0.0]),               # air, 0.1 mm off the slab's outer face
        ("air_triple_tip", [ro + e, d / 2 + e, 0.0]),        # air, 0.1 mm off Cu/alumina/air line at lip tip
        ("air_lip_corner", [ro + e, d / 2 + tl + e, 0.0]),   # air, 0.1 mm off the lip's outer corner
        ("diel_triple_tip", [ro - e, d / 2 - e / 2, 0.0]),   # alumina, near the tip triple line
        ("air_bore_face", [r - e, 0.0, 0.0]),                # air in the bore, off the slab's inner face
        ("air_triple_bore", [math.sqrt(r * r - (d / 2 + e) ** 2) - e, d / 2 + e, 0.0]),
        ("air_end_face", [r + t + max(s, 0) / 2, 0.0, z / 2 + e]),   # air beyond the slab end
        ("air_end_tip_corner", [ro + e, d / 2 + e, z / 2 + e]),      # air near tip/end corner (3 surfaces)
    ]
    return p


def build_mesh(g, h=1.0, path=None, gui=False, mult=(1.0, 1.0, 1.0, 1.0, 1.0)):
    """mult = extra size multipliers for (gap, gap edges, loop surfaces, bulk, loop/lip edges),
    for convergence tests."""
    import gmsh

    R, H, r, t, z, s, d, tl = g["R"], g["H"], g["r"], g["t"], g["z"], g["s"], g["d"], g["t_lip"]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("lgr")
    occ = gmsh.model.occ

    can = occ.addCylinder(0, 0, -H / 2, 0, 0, H, R)
    tube_o = occ.addCylinder(0, 0, -z / 2, 0, 0, z, r + t)
    tube_i = occ.addCylinder(0, 0, -z / 2 - 1, 0, 0, z + 2, r)
    tube, _ = occ.cut([(3, tube_o)], [(3, tube_i)])
    slit = occ.addBox(r - 2, -d / 2, -z / 2 - 1, t + 4, d, z + 2)
    tube, _ = occ.cut(tube, [(3, slit)])
    cond = tube
    if s > 0:
        x0 = r + t / 2
        lip1 = occ.addBox(x0, d / 2, -z / 2, r + t + s - x0, tl, z)
        lip2 = occ.addBox(x0, -d / 2 - tl, -z / 2, r + t + s - x0, tl, z)
        cond, _ = occ.fuse(cond, [(3, lip1), (3, lip2)])
    # dielectric slab: from the bore surface to the lip tips (or to the loop's outer surface)
    x_out = r + t + max(s, 0.0)
    slab = occ.addBox(r - 1, -d / 2, -z / 2, x_out - (r - 1), d, z)
    bore = occ.addCylinder(0, 0, -z / 2 - 1, 0, 0, z + 2, r)
    slab_c, _ = occ.cut([(3, slab)], [(3, bore)])
    if s <= 0:
        # no lips: clip slab to the loop's outer cylinder
        outer = occ.addCylinder(0, 0, -z / 2 - 1, 0, 0, z + 2, r + t)
        slab_c, _ = occ.intersect(slab_c, [(3, outer)])
    dom, _ = occ.cut([(3, can)], cond)
    out, outmap = occ.fragment(dom, slab_c)
    occ.synchronize()

    vols = [v for (dim, v) in gmsh.model.getEntities(3)]
    diel = []
    air = []
    for v in vols:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(3, v)
        if (ymax - ymin) < d * 1.01 and xmin > r - 1.5:
            diel.append(v)
        else:
            air.append(v)
    assert len(diel) == 1 and len(air) == 1, (diel, air)

    b_air = {abs(tag) for (_, tag) in gmsh.model.getBoundary([(3, air[0])], oriented=False)}
    b_die = {abs(tag) for (_, tag) in gmsh.model.getBoundary([(3, diel[0])], oriented=False)}
    interface = b_air & b_die
    exterior = (b_air | b_die) - interface
    shield, metal = [], []
    tol = 1e-3 * R
    for sf in exterior:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, sf)
        on_can = (max(abs(xmin), abs(xmax), abs(ymin), abs(ymax)) > R - tol and zmax - zmin > H - tol) \
            or (abs(zmin - zmax) < tol and abs(abs(zmin) - H / 2) < tol)
        (shield if on_can else metal).append(sf)

    gmsh.model.addPhysicalGroup(3, air, 1, name="air")
    gmsh.model.addPhysicalGroup(3, diel, 2, name="alumina")
    gmsh.model.addPhysicalGroup(2, shield, 1, name="shield")
    gmsh.model.addPhysicalGroup(2, metal, 2, name="loop_lips")

    # ---- mesh size fields ----
    h_gap = 0.45 * h * mult[0] * max(1.0, d / 0.635) ** 0.5     # in/around the gap
    h_edge = 0.25 * h * mult[1]                         # triple lines / lip edges
    h_loop = 2.0 * h * mult[2]                          # loop and lip surfaces
    h_bulk = max(R / 5.0, 4.0) * h * mult[3]            # far field (can)
    h_medge = 0.2 * h * mult[4]                         # loop rims, lip edges (singular H, sets Q_c)
    gap_surfs = list(b_die)
    gap_curves = sorted({abs(c) for (_, c) in gmsh.model.getBoundary([(2, sfc) for sfc in gap_surfs],
                                                                      combined=False, oriented=False)})
    f = gmsh.model.mesh.field
    f1 = f.add("Distance")
    f.setNumbers(f1, "SurfacesList", gap_surfs)
    f.setNumber(f1, "Sampling", 150)
    f2 = f.add("Threshold")
    f.setNumber(f2, "InField", f1)
    f.setNumber(f2, "SizeMin", h_gap)
    f.setNumber(f2, "SizeMax", h_bulk)
    f.setNumber(f2, "DistMin", 0.5)
    f.setNumber(f2, "DistMax", 12.0)
    f3 = f.add("Distance")
    f.setNumbers(f3, "CurvesList", gap_curves)
    f.setNumber(f3, "Sampling", 200)
    f4 = f.add("Threshold")
    f.setNumber(f4, "InField", f3)
    f.setNumber(f4, "SizeMin", h_edge)
    f.setNumber(f4, "SizeMax", h_bulk)
    f.setNumber(f4, "DistMin", 0.2)
    f.setNumber(f4, "DistMax", 4.0)
    f5 = f.add("Distance")
    f.setNumbers(f5, "SurfacesList", metal)
    f.setNumber(f5, "Sampling", 150)
    f6 = f.add("Threshold")
    f.setNumber(f6, "InField", f5)
    f.setNumber(f6, "SizeMin", h_loop)
    f.setNumber(f6, "SizeMax", h_bulk)
    f.setNumber(f6, "DistMin", 1.0)
    f.setNumber(f6, "DistMax", 15.0)
    metal_curves = sorted({abs(c) for (_, c) in gmsh.model.getBoundary([(2, sfc) for sfc in metal],
                                                                        combined=False, oriented=False)})
    f7 = f.add("Distance")
    f.setNumbers(f7, "CurvesList", metal_curves)
    f.setNumber(f7, "Sampling", 400)
    f8 = f.add("Threshold")
    f.setNumber(f8, "InField", f7)
    f.setNumber(f8, "SizeMin", h_medge)
    f.setNumber(f8, "SizeMax", h_bulk)
    f.setNumber(f8, "DistMin", 0.3)
    f.setNumber(f8, "DistMax", 6.0)
    fmin = f.add("Min")
    f.setNumbers(fmin, "FieldsList", [f2, f4, f6, f8])
    f.setAsBackgroundMesh(fmin)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 16)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_bulk)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)       # HXT
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.model.mesh.generate(3)

    etypes, etags, _ = gmsh.model.mesh.getElements(3)
    n_tet = sum(len(x) for x in etags)
    info = dict(n_tet=int(n_tet), h_scale=h, h_mult=list(mult), h_gap=h_gap, h_edge=h_edge, h_loop=h_loop, h_bulk=h_bulk, h_medge=h_medge,
                n_shield_surf=len(shield), n_metal_surf=len(metal))
    if path is not None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(path))
    if gui:
        gmsh.fltk.run()
    gmsh.finalize()
    return info


def zs_entry_values(f=F0, zs="exact"):
    """(Rs_param, Ls_param) for copper at f; parallel Rs || Ls model (see the docstring)."""
    rs = math.sqrt(math.pi * f * MU0 / SIGMA_CU)
    k = 2.0 if zs == "exact" else 1.0
    return k * rs, k * rs / (2 * math.pi * f)


def make_config(name, g, order, mesh_name=None, template=HERE / "palace_template.json", n_modes=2,
                shield_pec=False, zs="exact", f_rs=F0):
    cfg = json.loads(Path(template).read_text())
    cfg = copy.deepcopy(cfg)
    cfg["Problem"]["Output"] = f"postpro/{name}"
    cfg["Model"]["Mesh"] = f"mesh/{mesh_name or name}.msh"
    cfg["Solver"]["Order"] = order
    cfg["Solver"]["Eigenmode"]["N"] = n_modes
    for m in cfg["Domains"]["Materials"]:
        if m["Attributes"] == [2]:
            m["Permittivity"], m["LossTan"] = EPS_AL, TAND_AL
    for b in cfg["Boundaries"]["Impedance"]:
        b["Rs"], b["Ls"] = zs_entry_values(f_rs, zs)
    if shield_pec:   # diagnostic: lossless shield, so Q_c is the loop + lip loss alone
        cfg["Boundaries"]["Impedance"][0]["Attributes"] = [2]
        cfg["Boundaries"]["PEC"] = {"Attributes": [1]}
    cfg["Domains"]["Postprocessing"]["Probe"] = [
        {"Index": i + 1, "Center": xyz} for i, (_, xyz) in enumerate(probes(g))]
    return cfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--V", type=float, required=True, help="shield inner volume (L)")
    ap.add_argument("--a", type=float, required=True, help="r / R")
    ap.add_argument("--s", type=float, required=True, help="lip height (mm)")
    ap.add_argument("--d", type=float, default=0.635, help="gap width = alumina thickness (mm)")
    ap.add_argument("--h", type=float, default=1.0, help="global mesh-size scale")
    ap.add_argument("--order", type=int, default=2, help="Palace FE order")
    ap.add_argument("--mult", type=float, nargs=5, default=[1.0, 1.0, 1.0, 1.0, 1.0],
                    metavar=("GAP", "EDGE", "LOOP", "BULK", "MEDGE"), help="extra mesh-size multipliers (convergence tests)")
    ap.add_argument("--outdir", default=str(HERE))
    ap.add_argument("--no-mesh", action="store_true", help="reuse an existing mesh (see --mesh-name)")
    ap.add_argument("--mesh-name", default=None, help="mesh file stem (default: --name)")
    ap.add_argument("--shield-pec", action="store_true", help="diagnostic: PEC shield (loop/lip loss only)")
    ap.add_argument("--zs", choices=["exact", "s2"], default="exact",
                    help="copper input: exact = 2Rs || 2Rs/w (correct); s2 = Rs || Rs/w (first runs, half the loss)")
    ap.add_argument("--f-rs", type=float, default=F0, help="frequency (Hz) at which Rs is evaluated")
    ap.add_argument("--gui", action="store_true")
    args = ap.parse_args()

    out = Path(args.outdir)
    g = dims(args.V, args.a, args.s, args.d)
    info = {}
    if not args.no_mesh:
        info = build_mesh(g, args.h, out / "mesh" / f"{args.mesh_name or args.name}.msh", gui=args.gui,
                          mult=tuple(args.mult))
    else:
        prev = out / "cases" / f"{args.mesh_name}.json"
        if prev.exists():
            info = json.loads(prev.read_text())["mesh"]
    cfg = make_config(args.name, g, args.order, args.mesh_name, shield_pec=args.shield_pec, zs=args.zs,
                      f_rs=args.f_rs)
    (out / "configs").mkdir(parents=True, exist_ok=True)
    (out / "configs" / f"{args.name}.json").write_text(json.dumps(cfg, indent=2) + "\n")
    meta = dict(name=args.name, geometry=g, mesh=info, order=args.order, shield_pec=args.shield_pec,
                Rs=math.sqrt(math.pi * args.f_rs * MU0 / SIGMA_CU), copper_input=args.zs, f_Rs=args.f_rs,
                impedance=cfg["Boundaries"]["Impedance"], f_design=F0,
                eps_alumina=EPS_AL, tand_alumina=TAND_AL,
                probes=[dict(index=i + 1, label=lab, xyz_mm=xyz) for i, (lab, xyz) in enumerate(probes(g))])
    (out / "cases").mkdir(parents=True, exist_ok=True)
    (out / "cases" / f"{args.name}.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(dict(name=args.name, R=round(g["R"], 3), r=round(g["r"], 3), z=round(g["z"], 3),
                          s=g["s"], d=g["d"], clearance=round(g["clearance"], 2), **info)))


if __name__ == "__main__":
    main()
