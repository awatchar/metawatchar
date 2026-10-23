#!/usr/bin/env python3
"""Geometry, mesh and Palace configs for stage S2: a slab of unshielded loop-gap (MNG) cells in a
TEM-equivalent unit-cell waveguide (infinite transverse array).

Coordinates (mm in the mesh, Palace L0 = 1e-3):
  * x = propagation direction. The slab occupies 0 <= x <= L = N a_x; cell i is centred at
    x = (i + 1/2) a_x, y = 0, z = 0.
  * PEC plates at z = +-a_z/2 (E-field of the TEM mode along z), PMC walls at y = +-a_y/2 (H along y).
    With these mirror walls the cell row is an infinite 2D array in y-z.
  * Driven model: an empty line of length l_e on each side of the slab, numeric wave ports at
    x = -l_e and x = L + l_e, de-embedded (Palace "Offset" = l_e) to the slab faces x = 0 and x = L.
  * Eigen model: one cell, closed at x = +-a_x/2 by PEC or PMC (chosen in the config, same mesh).

Unit cell (LGR-A loop without its can, loop axis along y, parallel to the TEM H-field):
  * copper tube, inner radius r, wall t = 3 mm, length ell along y, centred;
  * one axial slit of width d through the wall, with two copper lips (thickness 3 mm, full length
    ell) extending radially outward by s from the outer wall surface, separated by d;
  * alumina (eps_r 9.8, tan_delta 4e-5) fills the slit and the space between the lips.
  * orient = "z": slit and lips point to +z (towards the top PEC plate). The lip plates are y-z planes,
    i.e. parallel to the incident E (z); the gap field is along x, perpendicular to E. The cell is
    mirror-symmetric in x, so the uniform incident E_z does not drive the LC mode (no first-order
    electric coupling / bianisotropy).
  * orient = "x": the same cell rotated by +90 deg about y: slit and lips point to +x (downstream).
    The gap field is along z, parallel to the incident E, so E drives the capacitor directly
    (electric coupling, magneto-electric/bianisotropic response; S11 != S22).
  * Conductors are not meshed: their surfaces are impedance boundaries (copper, Rs at F_RS).
    Palace combines the Rs, Ls (and Cs) of an "Impedance" boundary IN PARALLEL,
    1/Z = 1/Rs + 1/(j w Ls). The good-conductor impedance Z = Rs (1 + j) therefore needs
    Rs_param = 2 Rs, Ls_param = 2 Rs / w (config --zs exact, the default since 2026-10-23; run names
    with "_cu2"). The first S2 runs used Rs, Ls = Rs / w (--zs s2), i.e. Z = Rs (1 + j) / 2: half
    the copper loss (Q_c about 2x too high, frequencies and Re mu_eff unaffected). With a fixed
    parallel Rs_p, Ls_p the surface resistance seen at f is Re Z = 2 Rs(f_Rs) x^2 / (1 + x^2),
    x = f / f_Rs, against the true Rs(f_Rs) sqrt(x): within +-8 % over 130-200 MHz for f_Rs = 150 MHz.

Mesh attributes:
  volumes  1 = air, 2 = alumina
  surfaces 1 = top plate (z = +a_z/2), 2 = bottom plate (z = -a_z/2), 3 = PMC walls (y = +-a_y/2),
           4 = x = x_min face (port 1 / eigen x-wall), 5 = x = x_max face (port 2 / eigen x-wall),
           6 = copper (loops + lips of all cells)

Usage (local, gmsh >= 4.11 with OpenCASCADE):
  python s2_geometry.py mesh --kind driven --N 1 --orient z --s 6.0 [--recipe v2] [--name ...]
  python s2_geometry.py mesh --kind eigen  --orient z --s 6.0
  python s2_geometry.py mesh --kind empty
writes mesh/<name>.msh (Gmsh 2.2 binary, 2nd-order curved tets) and cases/<name>.json (geometry).
  python s2_geometry.py config --mesh <mesh name> --name <run> --type eigen --xwall pmc|pec --order 2
  python s2_geometry.py config --mesh <mesh name> --name <run> --type driven --order 2
      [--fmin 0.13 --fmax 0.20 --fstep 1e-5] [--both-ports] [--tol 1e-4]
writes configs/<run>.json.
"""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

# ---- fixed parameters ----
A_X = A_Y = A_Z = 70.0           # lattice constants (mm), about lambda/29 at 150 MHz
R_IN = 21.5                      # loop inner radius (mm), as LGR-A 0.5 L / a = 0.5
T_WALL = 3.0                     # loop wall (mm)
T_LIP = 3.0                      # lip plate thickness (mm)
D_GAP = 1.0                      # slit width = alumina thickness (mm)
ELL = 45.0                       # loop length along y (mm)
L_EMPTY = 2 * A_X                # empty line on each side of the slab (mm)
CLEAR_MIN = 3.0                  # minimum clearance lip tip -> plate / cell boundary (mm)
EPS_AL, TAND_AL = 9.8, 4.0e-5
SIGMA_CU = 5.8e7
MU0 = 4e-7 * math.pi
F_RS = 150e6                     # frequency at which the copper surface resistance is evaluated
RS = math.sqrt(math.pi * F_RS * MU0 / SIGMA_CU)   # 3.196e-3 ohm/sq
LS = RS / (2 * math.pi * F_RS)                     # internal (surface) inductance, H/sq


def zs_entry(attrs, f=F_RS, zs="exact"):
    """Palace "Impedance" entry for copper at frequency f (parallel Rs || Ls, see the docstring)."""
    rs = math.sqrt(math.pi * f * MU0 / SIGMA_CU)
    k = 2.0 if zs == "exact" else 1.0
    return {"Attributes": list(attrs), "Rs": k * rs, "Ls": k * rs / (2 * math.pi * f)}

# mesh recipes: sizes in mm. v2 = LGR-A production recipe (note: 0.2 mm on every loop/lip edge).
RECIPES = {
    "v2": dict(gap=0.565, gap_edge=0.25, loop=2.0, medge=0.2, bulk=8.0),
    "v3": dict(gap=0.40, gap_edge=0.18, loop=1.4, medge=0.12, bulk=6.0),   # refined (convergence)
    "c1": dict(gap=0.8, gap_edge=0.4, loop=3.0, medge=0.4, bulk=10.0),     # coarse (quick tuning)
}


def s_max():
    """Largest lip height with CLEAR_MIN to the plate / cell boundary."""
    return A_Z / 2 - CLEAR_MIN - (R_IN + T_WALL)


def dims(kind, N, orient, s):
    L = N * A_X if kind != "eigen" else A_X
    if kind == "eigen":
        x0, x1 = 0.0, A_X                      # single cell, centred at x = a_x/2
        N = 1
    elif kind == "empty":
        x0, x1 = -L_EMPTY, L_EMPTY
        N, L = 0, 0.0
    else:
        x0, x1 = -L_EMPTY, L + L_EMPTY
    return dict(kind=kind, N=N, orient=orient, a_x=A_X, a_y=A_Y, a_z=A_Z, r=R_IN, t=T_WALL,
                t_lip=T_LIP, d=D_GAP, ell=ELL, s=s, L=L, l_empty=(L_EMPTY if kind != "eigen" else 0.0),
                x_min=x0, x_max=x1, slab_faces=[0.0, L], s_max=s_max(),
                clearance=A_Z / 2 - (R_IN + T_WALL + s),
                cell_centres=[(i + 0.5) * A_X for i in range(N)])


def build_cell(occ, xc, g):
    """Return (copper dimtags, alumina dimtags) of one cell centred at (xc, 0, 0), slit at +z,
    rotated to +x for orient = "x"."""
    r, t, d, ell, s, tl = g["r"], g["t"], g["d"], g["ell"], g["s"], g["t_lip"]
    y0 = -ell / 2
    tube_o = occ.addCylinder(xc, y0, 0, 0, ell, 0, r + t)
    tube_i = occ.addCylinder(xc, y0 - 1, 0, 0, ell + 2, 0, r)
    tube, _ = occ.cut([(3, tube_o)], [(3, tube_i)])
    slit = occ.addBox(xc - d / 2, y0 - 1, r - 2, d, ell + 2, t + 4)
    tube, _ = occ.cut(tube, [(3, slit)])
    cond = tube
    if s > 0:
        z0 = r + t / 2
        lip1 = occ.addBox(xc + d / 2, y0, z0, tl, ell, r + t + s - z0)
        lip2 = occ.addBox(xc - d / 2 - tl, y0, z0, tl, ell, r + t + s - z0)
        cond, _ = occ.fuse(cond, [(3, lip1), (3, lip2)])
    z_out = r + t + max(s, 0.0)
    slab = occ.addBox(xc - d / 2, y0, r - 1, d, ell, z_out - (r - 1))
    bore = occ.addCylinder(xc, y0 - 1, 0, 0, ell + 2, 0, r)
    slab_c, _ = occ.cut([(3, slab)], [(3, bore)])
    if s <= 0:
        outer = occ.addCylinder(xc, y0 - 1, 0, 0, ell + 2, 0, r + t)
        slab_c, _ = occ.intersect(slab_c, [(3, outer)])
    if g["orient"] == "x":
        occ.rotate(cond + slab_c, xc, 0, 0, 0, 1, 0, math.pi / 2)    # +z -> +x
    return cond, slab_c


def build_mesh(g, recipe="v2", path=None, gui=False):
    import gmsh

    hs = RECIPES[recipe]
    ax, ay, az = g["a_x"], g["a_y"], g["a_z"]
    x0, x1 = g["x_min"], g["x_max"]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("s2")
    occ = gmsh.model.occ
    box = occ.addBox(x0, -ay / 2, -az / 2, x1 - x0, ay, az)
    conds, slabs = [], []
    for xc in g["cell_centres"]:
        c, sl = build_cell(occ, xc, g)
        conds += c
        slabs += sl
    if conds:
        dom, _ = occ.cut([(3, box)], conds)
        out, _ = occ.fragment(dom, slabs)
    occ.synchronize()

    vols = [v for (_, v) in gmsh.model.getEntities(3)]
    diel, air = [], []
    for v in vols:
        bb = gmsh.model.getBoundingBox(3, v)
        thin = min(bb[3] - bb[0], bb[5] - bb[2]) < g["d"] * 1.01      # slit slab: d thick in x or z
        (diel if thin else air).append(v)
    assert len(air) == 1 and len(diel) == len(g["cell_centres"]), (air, diel)

    b = {}
    for v in vols:
        b[v] = {abs(tag) for (_, tag) in gmsh.model.getBoundary([(3, v)], oriented=False)}
    allb = set().union(*b.values())
    interface = set()
    for dv in diel:
        interface |= b[air[0]] & b[dv]
    exterior = allb - interface
    tol = 1e-6 * ax
    groups = {1: [], 2: [], 3: [], 4: [], 5: [], 6: []}
    for sf in exterior:
        xmn, ymn, zmn, xmx, ymx, zmx = gmsh.model.getBoundingBox(2, sf)
        if abs(zmn - az / 2) < tol and abs(zmx - az / 2) < tol:
            groups[1].append(sf)
        elif abs(zmn + az / 2) < tol and abs(zmx + az / 2) < tol:
            groups[2].append(sf)
        elif (abs(abs(ymn) - ay / 2) < tol and abs(ymx - ymn) < tol):
            groups[3].append(sf)
        elif abs(xmn - x0) < tol and abs(xmx - x0) < tol:
            groups[4].append(sf)
        elif abs(xmn - x1) < tol and abs(xmx - x1) < tol:
            groups[5].append(sf)
        else:
            groups[6].append(sf)
    assert len(groups[1]) == 1 and len(groups[2]) == 1 and len(groups[3]) == 2, groups
    assert len(groups[4]) == 1 and len(groups[5]) == 1, groups

    gmsh.model.addPhysicalGroup(3, air, 1, name="air")
    if diel:
        gmsh.model.addPhysicalGroup(3, diel, 2, name="alumina")
    names = {1: "pec_top", 2: "pec_bottom", 3: "pmc_sides", 4: "xmin", 5: "xmax", 6: "copper"}
    for k, v in groups.items():
        if v:
            gmsh.model.addPhysicalGroup(2, v, k, name=names[k])

    f = gmsh.model.mesh.field
    fields = []

    def threshold(entities, dim, smin, dmin, dmax, sampling):
        fd = f.add("Distance")
        f.setNumbers(fd, "SurfacesList" if dim == 2 else "CurvesList", entities)
        f.setNumber(fd, "Sampling", sampling)
        ft = f.add("Threshold")
        f.setNumber(ft, "InField", fd)
        f.setNumber(ft, "SizeMin", smin)
        f.setNumber(ft, "SizeMax", hs["bulk"])
        f.setNumber(ft, "DistMin", dmin)
        f.setNumber(ft, "DistMax", dmax)
        fields.append(ft)

    if diel:
        gap_surfs = sorted(set().union(*[b[dv] for dv in diel]))
        gap_curves = sorted({abs(c) for (_, c) in gmsh.model.getBoundary(
            [(2, sfc) for sfc in gap_surfs], combined=False, oriented=False)})
        metal = groups[6]
        metal_curves = sorted({abs(c) for (_, c) in gmsh.model.getBoundary(
            [(2, sfc) for sfc in metal], combined=False, oriented=False)})
        threshold(gap_surfs, 2, hs["gap"], 0.5, 12.0, 150)
        threshold(gap_curves, 1, hs["gap_edge"], 0.2, 4.0, 200)
        threshold(metal, 2, hs["loop"], 1.0, 15.0, 150)
        threshold(metal_curves, 1, hs["medge"], 0.3, 6.0, 400)
        fmin = f.add("Min")
        f.setNumbers(fmin, "FieldsList", fields)
        f.setAsBackgroundMesh(fmin)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 16)
    gmsh.option.setNumber("Mesh.MeshSizeMax", hs["bulk"])
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)       # HXT
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.model.mesh.generate(3)
    _, etags, _ = gmsh.model.mesh.getElements(3)
    info = dict(recipe=recipe, sizes_mm=hs, n_tet=int(sum(len(x) for x in etags)),
                n_surf={names[k]: len(v) for k, v in groups.items()})
    if path is not None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(path))
    if gui:
        gmsh.fltk.run()
    gmsh.finalize()
    return info


def make_config(run, mesh, case, ptype, order, xwall="pmc", fmin=0.13, fmax=0.20, fstep=1e-5,
                tol=1e-4, max_samples=80, both_ports=False, offset=None, target=0.12, n_modes=3,
                zs="exact", f_rs=F_RS):
    g = case["geometry"]
    has_diel = g["N"] > 0
    mats = [{"Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0}]
    if has_diel:
        mats.append({"Attributes": [2], "Permeability": 1.0, "Permittivity": EPS_AL, "LossTan": TAND_AL})
    cfg = {
        "Problem": {"Type": "Eigenmode" if ptype == "eigen" else "Driven", "Verbose": 2,
                    "Output": f"postpro/{run}",
                    "OutputFormats": {"Paraview": False, "GridFunction": False}},
        "Model": {"Mesh": f"mesh/{mesh}.msh", "L0": 1.0e-3},
        "Domains": {"Materials": mats},
        "Boundaries": {"PEC": {"Attributes": [1, 2]}, "PMC": {"Attributes": [3]}},
        "Solver": {"Order": order, "Device": "CPU",
                   "Linear": {"Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-8, "MaxIts": 400}},
    }
    if has_diel:
        cfg["Domains"]["Postprocessing"] = {"Energy": [{"Index": 1, "Attributes": [2]},
                                                       {"Index": 2, "Attributes": [1]}]}
        cfg["Boundaries"]["Impedance"] = [zs_entry([6], f_rs, zs)]
    if ptype == "eigen":
        wall = "PEC" if xwall == "pec" else "PMC"
        cfg["Boundaries"][wall]["Attributes"] += [4, 5]
        cfg["Solver"]["Eigenmode"] = {"N": n_modes, "Tol": 1.0e-8, "Target": target, "Save": 0}
    else:
        off = g["l_empty"] if offset is None else offset
        cfg["Boundaries"]["WavePort"] = [
            {"Index": 1, "Attributes": [4], "Mode": 1, "Offset": off, "Excitation": 1,
             "PolarityAttributes": [1, 2]},
            {"Index": 2, "Attributes": [5], "Mode": 1, "Offset": off, "Excitation": 2 if both_ports else 0,
             "PolarityAttributes": [1, 2]},
        ]
        cfg["Solver"]["Driven"] = {
            "Samples": [{"Type": "Linear", "MinFreq": fmin, "MaxFreq": fmax, "FreqStep": fstep}],
            "Save": [], "AdaptiveTol": tol, "AdaptiveMaxSamples": max_samples}
        if tol == 0.0:
            del cfg["Solver"]["Driven"]["AdaptiveMaxSamples"]
    return cfg


def mesh_name(kind, N, orient, s, recipe):
    if kind == "empty":
        return "empty"
    n = "eig" if kind == "eigen" else f"N{N}"
    return f"{n}_{orient}_s{round(s * 100):04d}_{recipe}"


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("mesh")
    m.add_argument("--kind", choices=["driven", "eigen", "empty"], required=True)
    m.add_argument("--N", type=int, default=1)
    m.add_argument("--orient", choices=["z", "x"], default="z")
    m.add_argument("--s", type=float, default=0.0, help="lip height (mm)")
    m.add_argument("--recipe", default="v2", choices=sorted(RECIPES))
    m.add_argument("--name", default=None)
    m.add_argument("--gui", action="store_true")
    c = sub.add_parser("config")
    c.add_argument("--mesh", required=True)
    c.add_argument("--name", required=True)
    c.add_argument("--type", choices=["eigen", "driven"], required=True)
    c.add_argument("--order", type=int, default=2)
    c.add_argument("--xwall", choices=["pec", "pmc"], default="pmc")
    c.add_argument("--fmin", type=float, default=0.13)
    c.add_argument("--fmax", type=float, default=0.20)
    c.add_argument("--fstep", type=float, default=1e-5)
    c.add_argument("--tol", type=float, default=1e-4)
    c.add_argument("--max-samples", type=int, default=80)
    c.add_argument("--both-ports", action="store_true")
    c.add_argument("--offset", type=float, default=None, help="de-embedding offset (mm); default l_empty")
    c.add_argument("--target", type=float, default=0.12)
    c.add_argument("--nmodes", type=int, default=3)
    c.add_argument("--zs", choices=["exact", "s2"], default="exact",
                   help="copper input: exact = 2Rs || 2Rs/w (correct); s2 = Rs || Rs/w (first runs, half the loss)")
    c.add_argument("--f-rs", type=float, default=F_RS, help="frequency (Hz) at which Rs is evaluated")
    args = ap.parse_args()

    if args.cmd == "mesh":
        g = dims(args.kind, args.N, args.orient, args.s)
        name = args.name or mesh_name(args.kind, args.N, args.orient, args.s, args.recipe)
        info = build_mesh(g, args.recipe, HERE / "mesh" / f"{name}.msh", gui=args.gui)
        meta = dict(name=name, geometry=g, mesh=info, Rs=RS, Ls=LS, f_Rs=F_RS, sigma_cu=SIGMA_CU,
                    eps_alumina=EPS_AL, tand_alumina=TAND_AL)
        (HERE / "cases").mkdir(exist_ok=True)
        (HERE / "cases" / f"{name}.json").write_text(json.dumps(meta, indent=2) + "\n")
        print(json.dumps(dict(name=name, s=args.s, clearance=g["clearance"], n_tet=info["n_tet"])))
    else:
        case = json.loads((HERE / "cases" / f"{args.mesh}.json").read_text())
        cfg = make_config(args.name, args.mesh, case, args.type, args.order, args.xwall, args.fmin,
                          args.fmax, args.fstep, args.tol, args.max_samples, args.both_ports,
                          args.offset, args.target, args.nmodes, args.zs, args.f_rs)
        (HERE / "configs").mkdir(exist_ok=True)
        (HERE / "configs" / f"{args.name}.json").write_bytes((json.dumps(cfg, indent=2) + "\n").encode())
        run_meta = dict(run=args.name, mesh=args.mesh, type=args.type, order=args.order,
                        xwall=args.xwall if args.type == "eigen" else None,
                        offset_mm=cfg["Boundaries"]["WavePort"][0]["Offset"] if args.type == "driven" else None,
                        both_ports=args.both_ports, copper_input=args.zs, f_Rs=args.f_rs,
                        impedance=cfg["Boundaries"].get("Impedance"))
        (HERE / "runs").mkdir(exist_ok=True)
        (HERE / "runs" / f"{args.name}.json").write_bytes((json.dumps(run_meta, indent=2) + "\n").encode())
        print(args.name)


if __name__ == "__main__":
    main()
