#!/usr/bin/env python3
"""Geometry, mesh and Palace eigenmode configs for stage S4: thermally self-compensating MNG
media made by mixing loop-gap meta-atoms whose gap dielectrics have opposite temperature
coefficients (alumina-gap cells, tau_f < 0, and rutile-TiO2-gap cells, tau_f > 0).

The copper part of every cell is the S2/S3 meta-atom, built by `build_cell()` of
../s2_mng_slab/s2_geometry.py (imported, unchanged): copper loop r = 21.5 mm, wall 3 mm, length
45 mm along y, two 3 mm lips of height s, slit width d = 1.0 mm, slit and lips towards +z. The
enclosures (TEM unit-cell line, tight copper box) are those of ../s3_loaded_resonator
(`s3_geometry.dims`, imported, unchanged). Only the gap filling differs between cell types:

  "A"        alumina fills the whole slit and lip gap (the S2/S3 cell);
  "T<l>"     a rutile TiO2 block of length l (mm, along y, centred) fills the gap; the rest of
             the slit and lip gap (45 - l mm) is air. The TiO2 cell keeps the copper (inductance,
             coupling) of the alumina cell and is tuned only through l;
  "W<d>[s<s>]" alternative TiO2 cell: TiO2 fills the whole slit (full length), slit width d and
             lip height s (default 0, no lips) of its own (same loop radius, wall and length);
  "C<l>"     composite (meta-atom-level reference): alumina gap with a centred TiO2 block of
             length l (mm).

All models mesh only y >= 0 and use the y = 0 symmetry plane (PMC), as in S3; every gap filling
is mirror-symmetric in y, so every cell keeps its own y-mirror symmetry.

Mesh attributes
  volumes  1 = air outside the loops, 3 = air in the loop bores, 4 = air in the empty end sections
           of a line; per cell j (0-based, order of s3_geometry.dims): 100 + 10 j + 1 = alumina,
           + 2 = TiO2, + 3 = air inside the slit / lip gap (TiO2 cells only)
  surfaces 1 = z max wall, 2 = z min wall, 3 = y max wall, 4 = x min wall, 5 = x max wall,
           6 = copper of the cells, 7 = y = 0 symmetry plane (PMC)

Thermal model (applied in analyze_s4.py; Palace has no temperature input):
  * permittivity: eps_i(T) = eps_i (1 + tau_eps_i dT), entered directly in the FD configs (--dT);
  * thermal expansion: the copper (loops and box) scales uniformly with alpha_Cu, which shifts
    every frequency by -alpha_Cu dT exactly (Maxwell scaling); the gap width follows the
    expansion of the gap dielectric instead of copper, which adds +1/2 sum_j p_gap,j
    (alpha_d,j - alpha_Cu) dT (first-order capacitance perturbation, p_gap,j = electric energy
    share in the gap of cell j, dielectric + gap air).

Usage (local, gmsh >= 4.11 with OpenCASCADE):
  python s4_geometry.py mesh --kind line --N 1 --g 0 --cells T4.0                  # single cell box
  python s4_geometry.py mesh --kind cavity --nx 2 --ny 1 --nz 2 --cells A,T4.0,A,A
  python s4_geometry.py mesh --kind line --N 4 --g 35 --cells T4.0,A,A,A
  python s4_geometry.py config --mesh <mesh> --name <run> --order 2 [--xwall pmc|pec]
         [--walls cu|pec] --target 0.165 --nmodes 5 --rs-f 172.5e6 [--dT 25]
"""
import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "s2_mng_slab"))
sys.path.insert(0, str(HERE.parent / "s3_loaded_resonator"))
import s2_geometry as s2  # noqa: E402  (read-only reuse of the S2 unit cell)
import s3_geometry as s3  # noqa: E402  (read-only reuse of the S3 enclosures and copper input)

A0 = s2.A_X
RECIPES = s2.RECIPES

# ---- material data (see README "Material data"; [VERIFY] = not a verified VHF value) ----
#   tau_eps = temperature coefficient of eps_r (1/K); alpha = linear thermal expansion (1/K).
#   alumina 99.5-99.6 %: eps 9.8, tan d 4e-5 (measured 99.5 % grades near 450 MHz), tau_f(DR)
#   about -60 ppm/K -> tau_eps = -2 (tau_f + alpha) = +106 ppm/K [VERIFY]
#   rutile TiO2 (dense, doped polycrystalline): eps 100 [VERIFY], tan d 1e-4 [VERIFY at VHF],
#   tau_f(DR) about +450 ppm/K -> tau_eps = -916 ppm/K [VERIFY], alpha 8 ppm/K [VERIFY]
MAT = {
    "al": dict(eps=9.8, tand=4.0e-5, tau_eps=+106e-6, alpha=7.0e-6),
    "ti": dict(eps=100.0, tand=1.0e-4, tau_eps=-916e-6, alpha=8.0e-6),
}
ALPHA_CU = 17.0e-6            # copper, textbook 16.5-17 ppm/K [VERIFY]
VOL_AIR, VOL_BORE, VOL_END = 1, 3, 4


def cell_attr(j, what):
    return 100 + 10 * j + {"al": 1, "ti": 2, "gap": 3}[what]


def parse_spec(spec):
    spec = spec.strip()
    if spec == "A":
        return dict(kind="A", l=0.0)
    if spec[0] in "TC":
        return dict(kind=spec[0], l=float(spec[1:]))
    if spec[0] == "W":                      # full-length TiO2 slab, own slit width d and lip height s
        d, _, s = spec[1:].partition("s")
        return dict(kind="W", l=0.0, d=float(d), s=float(s) if s else 0.0)
    raise ValueError(spec)


def build_mesh(gd, specs, recipe="v2", path=None):
    import gmsh

    hs = RECIPES[recipe]
    b = gd["box"]
    cp = gd["cell"]
    assert len(specs) == len(gd["cells"]), (specs, gd["cells"])
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("s4")
    occ = gmsh.model.occ
    box = occ.addBox(b["x0"], 0.0, b["z0"], b["x1"] - b["x0"], b["y1"], b["z1"] - b["z0"])
    big = dict(x=b["x0"] - 5, dx=b["x1"] - b["x0"] + 10, z=b["z0"] - 5, dz=b["z1"] - b["z0"] + 10)

    def half(dimtags):
        hb = occ.addBox(b["x0"] - 1, 0.0, b["z0"] - 1, b["x1"] - b["x0"] + 2, b["y1"] + 1,
                        b["z1"] - b["z0"] + 2)
        part, _ = occ.intersect(dimtags, [(3, hb)])
        return part

    conds, tools = [], []          # tools: (dimtags, label)
    for j, ((xc, yc, zc), sp) in enumerate(zip(gd["cells"], specs)):
        cpj = dict(cp, d=sp["d"], s=sp["s"]) if sp["kind"] == "W" else cp
        c, sl = s2.build_cell(occ, xc, cpj)
        if yc or zc:
            occ.translate(c + sl, 0.0, yc, zc)
        conds += c
        if sp["kind"] == "A":
            tools.append((half(sl), ("al", j)))
        elif sp["kind"] == "W":
            tools.append((half(sl), ("ti", j)))
        else:
            l = sp["l"]
            yb = occ.addBox(big["x"], yc - l / 2, big["z"], big["dx"], l, big["dz"])
            cen, _ = occ.intersect(sl, [(3, yb)], removeObject=False, removeTool=True)
            out, _ = occ.cut(sl, cen, removeObject=True, removeTool=False)
            tools.append((half(cen), ("ti", j)))
            tools.append((half(out), ("al", j) if sp["kind"] == "C" else ("gap", j)))
        bore = occ.addCylinder(xc, yc - cp["ell"] / 2, zc, 0, cp["ell"], 0, cp["r"])
        tools.append((half([(3, bore)]), ("bore", j)))
    dom = occ.cut([(3, box)], conds)[0] if conds else [(3, box)]
    if gd["kind"] == "line" and gd["g"] > 0:
        for x0 in (b["x0"], gd["L"]):
            tools.append(([(3, occ.addBox(x0, 0.0, b["z0"], gd["g"], b["y1"], b["z1"] - b["z0"]))],
                          ("end", -1)))
    flat = [dt for (dts, _) in tools for dt in dts]
    owner = [lab for (dts, lab) in tools for _ in dts]
    outmap = occ.fragment(dom, flat)[1] if flat else [dom]
    occ.synchronize()

    label_of = {}
    for k, lab in enumerate(owner):
        for (_, v) in outmap[len(dom) + k]:
            label_of[v] = lab
    vols = [v for (_, v) in gmsh.model.getEntities(3)]
    groups_v = {}
    for v in vols:
        lab = label_of.get(v, ("air", -1))
        if lab[0] == "air":
            a = VOL_AIR
        elif lab[0] == "bore":
            a = VOL_BORE
        elif lab[0] == "end":
            a = VOL_END
        else:
            a = cell_attr(lab[1], lab[0])
        groups_v.setdefault(a, []).append(v)
    assert len(groups_v[VOL_AIR]) == 1, groups_v

    bnd = {v: [abs(t) for (_, t) in gmsh.model.getBoundary([(3, v)], oriented=False)] for v in vols}
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

    vname = {VOL_AIR: "air", VOL_BORE: "bore_air", VOL_END: "end_air"}
    for a, v in sorted(groups_v.items()):
        nm = vname.get(a) or f"c{(a - 100) // 10}_{ {1: 'al', 2: 'ti', 3: 'gap'}[a % 10] }"
        gmsh.model.addPhysicalGroup(3, v, a, name=nm.replace(" ", ""))
    snames = {1: "zmax", 2: "zmin", 3: "ymax", 4: "xmin", 5: "xmax", 6: "copper", 7: "sym_y0"}
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

    gap_vols = [v for a, vs in groups_v.items() if a >= 100 for v in vs]
    gap_surfs = sorted({sf for dv in gap_vols for sf in bnd[dv] if not on_sym(2, sf)})
    gap_curves = sorted({abs(c) for (_, c) in gmsh.model.getBoundary(
        [(2, sf) for sf in gap_surfs], combined=False, oriented=False) if not on_sym(1, abs(c))})
    metal = gs[6]
    metal_curves = sorted({abs(c) for (_, c) in gmsh.model.getBoundary(
        [(2, sf) for sf in metal], combined=False, oriented=False) if not on_sym(1, abs(c))})
    if gap_vols:
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
    vol = {str(a): sum(occ.getMass(3, v) for v in vs) * 2 for a, vs in groups_v.items()}
    info = dict(recipe=recipe, sizes_mm=hs, n_tet=int(sum(len(x) for x in etags)),
                n_surf={snames[k]: len(v) for k, v in gs.items()}, volumes_full_mm3=vol)
    if path is not None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(path))
    gmsh.finalize()
    return info


def vol_label(a):
    a = int(a)
    if a < 100:
        return {VOL_AIR: "air", VOL_BORE: "bore", VOL_END: "end"}[a]
    return f"c{(a - 100) // 10}_{ {1: 'al', 2: 'ti', 3: 'gap'}[a % 10] }".replace(" ", "")


def make_config(run, mesh, case, order, xwall="pec", walls="cu", target=0.165, n_modes=5,
                rs_f=172.5e6, dT=0.0, dT_mats=("al", "ti"), eps_override=None, zs="exact"):
    g = case["geometry"]
    attrs = sorted(int(a) for a in case["mesh"]["volumes_full_mm3"])
    mats = []
    air = [a for a in attrs if a < 100 or a % 10 == 3]
    mats.append({"Attributes": air, "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0})
    eps_used = {}
    for key, code in (("al", 1), ("ti", 2)):
        aa = [a for a in attrs if a >= 100 and a % 10 == code]
        if not aa:
            continue
        m = dict(MAT[key])
        if eps_override and key in eps_override:
            m["eps"] = eps_override[key]
        e = m["eps"] * (1 + (m["tau_eps"] * dT if key in dT_mats else 0.0))
        eps_used[key] = e
        mats.append({"Attributes": aa, "Permeability": 1.0, "Permittivity": e, "LossTan": m["tand"]})
    energy = [{"Index": i + 1, "Attributes": [a]} for i, a in enumerate(attrs)]
    bnd = {"PMC": {"Attributes": [7]}}
    imp = [s3.zs_entry([6], rs_f, zs)] if g["N"] > 0 else []
    if g["kind"] == "line":
        bnd["PEC"] = {"Attributes": [1, 2]}
        bnd["PMC"]["Attributes"] += [3]
        (bnd["PEC"] if xwall == "pec" else bnd["PMC"])["Attributes"] += [4, 5]
    elif walls == "pec":
        bnd["PEC"] = {"Attributes": [1, 2, 3, 4, 5]}
    else:
        imp.append(s3.zs_entry([1, 2, 3, 4, 5], rs_f, zs))
    bnd["Impedance"] = imp
    cfg = {
        "Problem": {"Type": "Eigenmode", "Verbose": 2, "Output": f"postpro/{run}",
                    "OutputFormats": {"Paraview": False, "GridFunction": False}},
        "Model": {"Mesh": f"mesh/{mesh}.msh", "L0": 1.0e-3},
        "Domains": {"Materials": mats, "Postprocessing": {"Energy": energy}},
        "Boundaries": bnd,
        "Solver": {"Order": order, "Device": "CPU",
                   "Linear": {"Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-8, "MaxIts": 400},
                   "Eigenmode": {"N": n_modes, "Tol": 1.0e-8, "Target": target, "Save": 0}},
    }
    meta = dict(run=run, mesh=mesh, order=order, xwall=xwall, walls=walls, target_GHz=target,
                n_modes=n_modes, rs_f_Hz=rs_f, zs=zs, dT=dT, dT_mats=list(dT_mats),
                eps_used=eps_used, mat=MAT, alpha_cu=ALPHA_CU,
                energy_index={str(i + 1): vol_label(a) for i, a in enumerate(attrs)},
                tand={vol_label(a): (MAT["al"]["tand"] if a % 10 == 1 else MAT["ti"]["tand"])
                      for a in attrs if a >= 100 and a % 10 in (1, 2)})
    return cfg, meta


def spec_tag(specs):
    return "".join(s.replace(".", "p") for s in specs)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("mesh")
    m.add_argument("--kind", choices=["line", "cavity", "empty"], required=True)
    m.add_argument("--N", type=int, default=1)
    m.add_argument("--g", type=float, default=0.0)
    m.add_argument("--nx", type=int, default=1)
    m.add_argument("--ny", type=int, default=1)
    m.add_argument("--nz", type=int, default=1)
    m.add_argument("--ext", type=float, default=0.0)
    m.add_argument("--s", type=float, default=6.2)
    m.add_argument("--d", type=float, default=None)
    m.add_argument("--cells", default="", help="comma list of cell specs: A, T<l>, C<l>")
    m.add_argument("--recipe", default="v2", choices=sorted(RECIPES))
    m.add_argument("--name", required=True)
    c = sub.add_parser("config")
    c.add_argument("--mesh", required=True)
    c.add_argument("--name", required=True)
    c.add_argument("--order", type=int, default=2)
    c.add_argument("--xwall", choices=["pec", "pmc"], default="pec")
    c.add_argument("--walls", choices=["cu", "pec"], default="cu")
    c.add_argument("--target", type=float, default=0.165)
    c.add_argument("--nmodes", type=int, default=5)
    c.add_argument("--rs-f", type=float, default=172.5e6)
    c.add_argument("--dT", type=float, default=0.0)
    c.add_argument("--dT-mats", default="al,ti")
    c.add_argument("--eps-ti", type=float, default=None)
    args = ap.parse_args()

    if args.cmd == "mesh":
        specs_s = [s.strip() for s in args.cells.split(",") if s.strip()]
        specs = [parse_spec(s) for s in specs_s]
        cp = s3.cell_params(1.0, args.s, args.d)
        gd = s3.dims(args.kind, cp, args.N, args.g, args.nx, args.ny, args.nz, args.ext)
        info = build_mesh(gd, specs, args.recipe, HERE / "mesh" / f"{args.name}.msh")
        meta = dict(name=args.name, geometry=gd, cells=specs_s, mesh=info, mat=MAT)
        (HERE / "cases").mkdir(exist_ok=True)
        (HERE / "cases" / f"{args.name}.json").write_bytes((json.dumps(meta, indent=2) + "\n").encode())
        print(json.dumps(dict(name=args.name, n_tet=info["n_tet"], volumes=info["volumes_full_mm3"])))
    else:
        case = json.loads((HERE / "cases" / f"{args.mesh}.json").read_text())
        eo = {"ti": args.eps_ti} if args.eps_ti else None
        cfg, meta = make_config(args.name, args.mesh, case, args.order, args.xwall, args.walls,
                                args.target, args.nmodes, args.rs_f, args.dT,
                                tuple(args.dT_mats.split(",")), eo)
        meta["cells"] = case["cells"]
        (HERE / "configs").mkdir(exist_ok=True)
        (HERE / "configs" / f"{args.name}.json").write_bytes((json.dumps(cfg, indent=2) + "\n").encode())
        (HERE / "runs").mkdir(exist_ok=True)
        (HERE / "runs" / f"{args.name}.json").write_bytes((json.dumps(meta, indent=2) + "\n").encode())
        print(args.name)


if __name__ == "__main__":
    main()
