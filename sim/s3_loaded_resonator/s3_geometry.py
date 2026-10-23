#!/usr/bin/env python3
"""Geometry, mesh and Palace eigenmode configs for stage S3: resonators loaded with the S2 MNG
meta-atom array (loop-gap cells, s2_mng_slab).

The unit cell is built by `build_cell()` of ../s2_mng_slab/s2_geometry.py (imported, unchanged):
copper loop r = 21.5 mm, t = 3 mm, length 45 mm along y, lips 3 mm x s, alumina gap d, slit and
lips towards +z ("z" orientation of S2). Cells can be scaled (study 3) and translated.

All models use the y = 0 mirror plane of the cell: only y >= 0 is meshed and y = 0 is a PMC
symmetry plane. This keeps every mode whose loop flux (H_y) is even in y, which includes all
loop-gap (LC) modes of the array; modes odd in y are not computed.

Case kinds
  line   : S2 TEM unit-cell waveguide (PEC z = +-a/2, PMC y = a/2) holding N cells along x,
           0 <= x <= N a, with empty line sections of length g at both ends, closed by walls at
           x = -g and x = N a + g (PEC: short-circuited line; PMC possible). Study 1.
  cavity : rectangular copper box around an n_x x n_y x n_z array (lattice a, loop axes along y),
           cells at ((i+1/2) a, (j+1/2-n_y/2) a, (k+1/2) a). The box is the lattice boundary
           (walls a/2 from the cell centres) enlarged by `ext` on every side. Study 2 / 3.
  empty  : the same box without cells (study 2 reference).

Mesh attributes
  volumes  1 = air (outside the loops), 2 = alumina, 3 = air inside the loop bores (r, length ell),
           4 = air in the empty end sections of a line (x < 0 or x > N a)
  surfaces 1 = z max wall, 2 = z min wall, 3 = y max wall, 4 = x min wall, 5 = x max wall,
           6 = copper of the cells, 7 = y = 0 symmetry plane (PMC)

Usage (local, gmsh >= 4.11 with OpenCASCADE):
  python s3_geometry.py mesh --kind line --N 4 --g 35 [--scale 1 --s 6.2 --d 1.0] [--recipe v2]
  python s3_geometry.py mesh --kind cavity --nx 2 --ny 1 --nz 2 [--ext 0]
  python s3_geometry.py mesh --kind empty  --nx 1 --ny 1 --nz 2
  python s3_geometry.py config --mesh <mesh> --name <run> --order 2 --xwall pec|pmc
         --walls cu|pec --target 0.12 --nmodes 6 [--rs-f 150e6]
"""
import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "s2_mng_slab"))
import s2_geometry as s2  # noqa: E402  (read-only reuse of the S2 unit cell)

A0 = s2.A_X                       # 70 mm lattice of S2
EPS_AL, TAND_AL = s2.EPS_AL, s2.TAND_AL
SIGMA_CU, MU0 = s2.SIGMA_CU, s2.MU0
RECIPES = s2.RECIPES


def rs_of(f):
    return math.sqrt(math.pi * f * MU0 / SIGMA_CU)


def zs_entry(attrs, f, zs):
    """Palace "Impedance" entry for a good conductor at frequency f.

    Palace combines Rs, Ls (and Cs) of an impedance boundary IN PARALLEL: 1/Z = 1/Rs + 1/(j w Ls).
      zs = "s2"   : Rs, Ls = Rs/w as in ../s2_mng_slab. In Palace's parallel model this gives
                    Z = Rs (1 + j)/2, i.e. HALF the copper loss of a real conductor (kept only to
                    reproduce S2 consistently).
      zs = "exact": parallel Rs_p = 2 Rs, Ls_p = 2 Rs / w, which gives Z = Rs (1 + j) at f, the
                    good-conductor surface impedance.
    """
    rs = rs_of(f)
    w = 2 * math.pi * f
    if zs == "s2":
        return {"Attributes": attrs, "Rs": rs, "Ls": rs / w}
    return {"Attributes": attrs, "Rs": 2 * rs, "Ls": 2 * rs / w}


def cell_params(scale=1.0, s=6.2, d=None):
    """S2 cell scaled by `scale` (all lengths), with lip height s and gap d given explicitly."""
    return dict(orient="z", r=s2.R_IN * scale, t=s2.T_WALL * scale, t_lip=s2.T_LIP * scale,
                ell=s2.ELL * scale, d=(s2.D_GAP * scale if d is None else d), s=s, scale=scale,
                a=A0 * scale)


def dims(kind, cp, N=1, g=0.0, nx=1, ny=1, nz=1, ext=0.0):
    a = cp["a"]
    if kind == "line":
        cells = [((i + 0.5) * a, 0.0, 0.0) for i in range(N)]
        box = dict(x0=-g, x1=N * a + g, y1=a / 2, z0=-a / 2, z1=a / 2)
        return dict(kind=kind, N=N, g=g, L=N * a, cells=cells, box=box, cell=cp, a=a)
    cells = []
    if kind == "cavity":
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    cells.append(((i + 0.5) * a, (j + 0.5 - ny / 2) * a, (k + 0.5) * a))
    box = dict(x0=-ext, x1=nx * a + ext, y1=ny * a / 2 + ext, z0=-ext, z1=nz * a + ext)
    return dict(kind=kind, N=len(cells), nx=nx, ny=ny, nz=nz, ext=ext, cells=cells, box=box,
                cell=cp, a=a, size_mm=[box["x1"] - box["x0"], 2 * box["y1"], box["z1"] - box["z0"]])


def build_mesh(gd, recipe="v2", path=None):
    import gmsh

    hs = RECIPES[recipe]
    b = gd["box"]
    cp = gd["cell"]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("s3")
    occ = gmsh.model.occ
    box = occ.addBox(b["x0"], 0.0, b["z0"], b["x1"] - b["x0"], b["y1"], b["z1"] - b["z0"])
    conds, slabs, bores = [], [], []
    for (xc, yc, zc) in gd["cells"]:
        c, sl = s2.build_cell(occ, xc, cp)
        if yc or zc:
            occ.translate(c + sl, 0.0, yc, zc)
        conds += c
        slabs += sl
        bores.append((3, occ.addCylinder(xc, yc - cp["ell"] / 2, zc, 0, cp["ell"], 0, cp["r"])))
    tools = []
    if conds:
        dom, _ = occ.cut([(3, box)], conds)
        # keep only the y >= 0 halves of alumina slabs and bores
        for obj in slabs + bores:
            hb = occ.addBox(b["x0"] - 1, 0.0, b["z0"] - 1, b["x1"] - b["x0"] + 2, b["y1"] + 1,
                            b["z1"] - b["z0"] + 2)
            part, _ = occ.intersect([obj], [(3, hb)])
            tools += part
    else:
        dom = [(3, box)]
    if gd["kind"] == "line" and gd["g"] > 0:
        tools.append((3, occ.addBox(b["x0"], 0.0, b["z0"], gd["g"], b["y1"], b["z1"] - b["z0"])))
        tools.append((3, occ.addBox(gd["L"], 0.0, b["z0"], gd["g"], b["y1"], b["z1"] - b["z0"])))
    if tools:
        occ.fragment(dom, tools)
    occ.synchronize()

    vols = [v for (_, v) in gmsh.model.getEntities(3)]
    groups_v = {1: [], 2: [], 3: [], 4: []}
    for v in vols:
        bb = gmsh.model.getBoundingBox(3, v)
        cx, cy, cz = occ.getCenterOfMass(3, v)
        thin = min(bb[3] - bb[0], bb[5] - bb[2]) < cp["d"] * 1.01
        in_bore = any((cx - xc) ** 2 + (cz - zc) ** 2 < cp["r"] ** 2 and abs(cy - yc) < cp["ell"] / 2
                      and bb[3] - bb[0] < 2 * cp["r"] + 1e-6 for (xc, yc, zc) in gd["cells"])
        if thin and gd["N"] > 0:
            groups_v[2].append(v)
        elif in_bore:
            groups_v[3].append(v)
        elif gd["kind"] == "line" and (cx < 0 or cx > gd["L"]):
            groups_v[4].append(v)
        else:
            groups_v[1].append(v)
    n_half = sum(1 for (xc, yc, zc) in gd["cells"] if yc >= 0)   # cells cut by / lying in y >= 0
    assert len(groups_v[1]) == 1, groups_v
    assert len(groups_v[2]) == len(groups_v[3]) == n_half, groups_v

    bnd = {}
    for v in vols:
        bnd[v] = [abs(t) for (_, t) in gmsh.model.getBoundary([(3, v)], oriented=False)]
    count = {}
    for v in vols:
        for sf in bnd[v]:
            count[sf] = count.get(sf, 0) + 1
    exterior = [sf for sf, n in count.items() if n == 1]
    tol = 1e-6 * A0
    gs = {k: [] for k in range(1, 8)}
    for sf in exterior:
        xmn, ymn, zmn, xmx, ymx, zmx = gmsh.model.getBoundingBox(2, sf)
        if abs(zmn - b["z1"]) < tol and abs(zmx - b["z1"]) < tol:
            gs[1].append(sf)
        elif abs(zmn - b["z0"]) < tol and abs(zmx - b["z0"]) < tol:
            gs[2].append(sf)
        elif abs(ymn - b["y1"]) < tol and abs(ymx - b["y1"]) < tol:
            gs[3].append(sf)
        elif abs(ymn) < tol and abs(ymx) < tol:
            gs[7].append(sf)
        elif abs(xmn - b["x0"]) < tol and abs(xmx - b["x0"]) < tol:
            gs[4].append(sf)
        elif abs(xmn - b["x1"]) < tol and abs(xmx - b["x1"]) < tol:
            gs[5].append(sf)
        else:
            gs[6].append(sf)
    vnames = {1: "air", 2: "alumina", 3: "bore_air", 4: "end_air"}
    snames = {1: "zmax", 2: "zmin", 3: "ymax", 4: "xmin", 5: "xmax", 6: "copper", 7: "sym_y0"}
    for k, v in groups_v.items():
        if v:
            gmsh.model.addPhysicalGroup(3, v, k, name=vnames[k])
    for k, v in gs.items():
        if v:
            gmsh.model.addPhysicalGroup(2, v, k, name=snames[k])

    def on_sym(dim, tag):
        bb = gmsh.model.getBoundingBox(dim, tag)
        return abs(bb[1]) < tol and abs(bb[4]) < tol

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

    if groups_v[2]:
        gap_surfs = sorted({sf for dv in groups_v[2] for sf in bnd[dv] if not on_sym(2, sf)})
        gap_curves = sorted({abs(c) for (_, c) in gmsh.model.getBoundary(
            [(2, sf) for sf in gap_surfs], combined=False, oriented=False) if not on_sym(1, abs(c))})
        metal = gs[6]
        metal_curves = sorted({abs(c) for (_, c) in gmsh.model.getBoundary(
            [(2, sf) for sf in metal], combined=False, oriented=False) if not on_sym(1, abs(c))})
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
    vol = {vnames[k]: sum(occ.getMass(3, v) for v in vs) * 2 for k, vs in groups_v.items() if vs}
    info = dict(recipe=recipe, sizes_mm=hs, n_tet=int(sum(len(x) for x in etags)),
                n_surf={snames[k]: len(v) for k, v in gs.items()}, volumes_full_mm3=vol)
    if path is not None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(path))
    gmsh.finalize()
    return info


def make_config(run, mesh, case, order, xwall="pec", walls="cu", target=0.12, n_modes=4, rs_f=150e6,
                rs_wall_f=None, zs="exact"):
    g = case["geometry"]
    has_cells = g["N"] > 0
    present = {"air": 1, "alumina": 2, "bore_air": 3, "end_air": 4}
    vattr = sorted(present[k] for k in case["mesh"]["volumes_full_mm3"])
    mats = [{"Attributes": [v for v in vattr if v != 2], "Permeability": 1.0, "Permittivity": 1.0,
             "LossTan": 0.0}]
    if has_cells:
        mats.append({"Attributes": [2], "Permeability": 1.0, "Permittivity": EPS_AL, "LossTan": TAND_AL})
    bnd = {"PMC": {"Attributes": [7]}}
    imp = []
    if has_cells:
        imp.append(zs_entry([6], rs_f, zs))
    if g["kind"] == "line":
        bnd["PEC"] = {"Attributes": [1, 2]}
        bnd["PMC"]["Attributes"] += [3]
        (bnd["PEC"] if xwall == "pec" else bnd["PMC"])["Attributes"] += [4, 5]
    elif walls == "pec":
        bnd["PEC"] = {"Attributes": [1, 2, 3, 4, 5]}
    else:
        fw = rs_wall_f or rs_f
        imp.append(zs_entry([1, 2, 3, 4, 5], fw, zs))
    if imp:
        bnd["Impedance"] = imp
    cfg = {
        "Problem": {"Type": "Eigenmode", "Verbose": 2, "Output": f"postpro/{run}",
                    "OutputFormats": {"Paraview": False, "GridFunction": False}},
        "Model": {"Mesh": f"mesh/{mesh}.msh", "L0": 1.0e-3},
        "Domains": {"Materials": mats},
        "Boundaries": bnd,
        "Solver": {"Order": order, "Device": "CPU",
                   "Linear": {"Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-8, "MaxIts": 400},
                   "Eigenmode": {"N": n_modes, "Tol": 1.0e-8, "Target": target, "Save": 0}},
    }
    if has_cells:
        en = [{"Index": 1, "Attributes": [2]}, {"Index": 2, "Attributes": [1]},
              {"Index": 3, "Attributes": [3]}]
        if g["kind"] == "line" and g["g"] > 0:
            en.append({"Index": 4, "Attributes": [4]})
        cfg["Domains"]["Postprocessing"] = {"Energy": en}
    return cfg


def mesh_name(kind, a, recipe):
    sc = "" if abs(a.scale - 1) < 1e-9 else f"_k{round(a.scale * 100):03d}"
    cell = f"s{round(a.s * 100):04d}" + ("" if a.d is None else f"_d{round(a.d * 100):03d}")
    if kind == "line":
        base = f"L{a.N}_g{round(a.g):02d}"
    elif kind == "cavity":
        base = f"C{a.nx}{a.ny}{a.nz}_e{round(a.ext):02d}"
    else:
        return f"E{a.nx}{a.ny}{a.nz}_e{round(a.ext):02d}{sc}"
    return f"{base}{sc}_{cell}_{recipe}"


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("mesh")
    m.add_argument("--kind", choices=["line", "cavity", "empty"], required=True)
    m.add_argument("--N", type=int, default=1)
    m.add_argument("--g", type=float, default=0.0, help="empty end section of the line (mm)")
    m.add_argument("--nx", type=int, default=1)
    m.add_argument("--ny", type=int, default=1)
    m.add_argument("--nz", type=int, default=1)
    m.add_argument("--ext", type=float, default=0.0, help="cavity enlargement beyond the lattice box (mm)")
    m.add_argument("--scale", type=float, default=1.0)
    m.add_argument("--s", type=float, default=6.2)
    m.add_argument("--d", type=float, default=None)
    m.add_argument("--recipe", default="v2", choices=sorted(RECIPES))
    m.add_argument("--name", default=None)
    c = sub.add_parser("config")
    c.add_argument("--mesh", required=True)
    c.add_argument("--name", required=True)
    c.add_argument("--order", type=int, default=2)
    c.add_argument("--xwall", choices=["pec", "pmc"], default="pec")
    c.add_argument("--walls", choices=["cu", "pec"], default="cu")
    c.add_argument("--target", type=float, default=0.12)
    c.add_argument("--nmodes", type=int, default=4)
    c.add_argument("--rs-f", type=float, default=150e6)
    c.add_argument("--rs-wall-f", type=float, default=None)
    c.add_argument("--zs", choices=["exact", "s2"], default="exact",
                   help="copper surface impedance input (see zs_entry); s2 = as in S2 (half loss)")
    args = ap.parse_args()

    if args.cmd == "mesh":
        cp = cell_params(args.scale, args.s, args.d)
        gd = dims(args.kind, cp, args.N, args.g, args.nx, args.ny, args.nz, args.ext)
        name = args.name or mesh_name(args.kind, args, args.recipe)
        info = build_mesh(gd, args.recipe, HERE / "mesh" / f"{name}.msh")
        meta = dict(name=name, geometry=gd, mesh=info, eps_alumina=EPS_AL, tand_alumina=TAND_AL,
                    sigma_cu=SIGMA_CU)
        (HERE / "cases").mkdir(exist_ok=True)
        (HERE / "cases" / f"{name}.json").write_bytes((json.dumps(meta, indent=2) + "\n").encode())
        print(json.dumps(dict(name=name, n_tet=info["n_tet"], volumes=info["volumes_full_mm3"])))
    else:
        case = json.loads((HERE / "cases" / f"{args.mesh}.json").read_text())
        cfg = make_config(args.name, args.mesh, case, args.order, args.xwall, args.walls, args.target,
                          args.nmodes, args.rs_f, args.rs_wall_f, args.zs)
        (HERE / "configs").mkdir(exist_ok=True)
        (HERE / "configs" / f"{args.name}.json").write_bytes((json.dumps(cfg, indent=2) + "\n").encode())
        meta = dict(run=args.name, mesh=args.mesh, order=args.order, xwall=args.xwall, walls=args.walls,
                    target_GHz=args.target, n_modes=args.nmodes, rs_f_Hz=args.rs_f,
                    rs_wall_f_Hz=args.rs_wall_f or args.rs_f, zs=args.zs)
        (HERE / "runs").mkdir(exist_ok=True)
        (HERE / "runs" / f"{args.name}.json").write_bytes((json.dumps(meta, indent=2) + "\n").encode())
        print(args.name)


if __name__ == "__main__":
    main()
