#!/usr/bin/env python3
"""Full-wave model of the S2 measurement fixture (50-ohm air microstrip TEM line with constant-
impedance tapers) with and without a slab of LGR meta-atoms: Gmsh geometry/mesh and Palace configs.

Fixture (mm; hardware/tem_fixture/tem_fixture_cad.py, design note in README.md):
  * ground plane z = 0 (the domain floor);
  * top strip, 2 mm thick, underside at h = 70 mm, width w = 344 mm, uniform for |x| <= 210 mm;
  * constant-impedance tapers 210 <= |x| <= 510 mm: underside height 70 -> 3 mm and width
    344 -> 14.7 mm, both linear in x (w/h = 4.91 held), a planar trapezoid;
  * ports: a 50-ohm lumped port across the 3 mm gap at each tip (x = +-510 mm), standing in for
    the N-type jack launch; the 2x-thru de-embedding removes it together with the taper;
  * 4 PTFE support posts (diameter 12 mm, 70 mm tall) at x = +-195, y = +-162 mm;
  * variant "open": air box with first-order absorbing faces. Default |x| <= 700, y, z <= 450 mm is
    TOO SMALL (the absorbing faces take up the line's near field); production runs use
    --box 2510 2500 2500 --bulk 150 (converged against 3510 3500 3500);
  * variant "closed": rectangular-coax TEM cell made from the same parts plus side walls at
    |y| = 252 mm (the ground-plate edges), a lid at z = 300 mm and end walls at |x| = 520 mm
    (the jack flanges); no absorbing boundary. --w 280 narrows the strip (50 ohm in the box).
  * variant "closedline": the closed uniform section only (|x| <= 150 mm), wave ports, no tapers,
    posts or launches (reference model; the Palace port-mode solve stalled on it, see results.md);
  * --ny 5: 5 cells across (y = 0, +-70, +-140 mm; the centre cell is cut by the symmetry plane).
Half model (y >= 0): the plane y = 0 is a magnetic wall (PMC) for the fundamental, even mode.
The lumped port of the half model therefore has R = 100 ohm (half the current at the same
voltage), which reproduces the S-parameters of the full 50-ohm structure.

Meta-atom slab: N_y = 4 cells across (centres y = +-35, +-105 mm; the half model holds the two at
y = 35, 105) by N_x = 1 (x = 0) or 2 (x = +-35 mm) along the line. The cell is exactly the S2
unshielded LGR cell (../s2_mng_slab/s2_geometry.py: r = 21.5 mm, wall 3 mm, length 45 mm along y,
lips s = 6.2 mm towards +z, alumina gap 1.0 mm), with its axis at mid-height z = h/2 = 35 mm as in
the S2 unit cell (lip tips 4.3 mm below the strip). Copper uses the corrected Palace input
(parallel Rs || Ls: Rs_p = 2 Rs, Ls_p = 2 Rs / w at 150 MHz, sigma = 5.8e7 S/m).

Fixture metal: aluminium 6061-T6, sigma = 2.5e7 S/m (about 43 % IACS), same parallel input at
150 MHz, on the strip, the ground and (closed variant) the walls; or PEC (--metal pec) to isolate
radiation (open variant: 1 - |S11|^2 - |S21|^2 is then the radiated fraction).

Mesh attributes
  volumes  1 air, 2 alumina, 3 PTFE
  surfaces 1 ground (z = 0), 2 strip (all faces), 3 symmetry plane y = 0, 4 outer box faces
           (absorbing for "open", aluminium walls for "closed"), 5 port 1 (x = -510),
           6 port 2 (x = +510), 7 copper (meta-atoms)

Usage (local, gmsh >= 4.11 with OpenCASCADE):
  python fixture_geometry.py mesh --variant open --nx 0          # empty fixture (the 2x-thru)
  python fixture_geometry.py mesh --variant open --nx 1          # 4 x 1 slab
  python fixture_geometry.py config --mesh open_nx0 --name open_nx0_al_o2 --order 2 --metal al --band thru
  python fixture_geometry.py xsec --variant open                 # 2D cross-section (boundary mode)
  python fixture_geometry.py xconfig --variant open --sym even --freq 0.5 --name xs_open_even_f500
"""
import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "s2_mng_slab"))
import s2_geometry as s2  # noqa: E402  (cell geometry, mesh recipe, copper input)

# ---- fixture parameters (mm) ----
H, T_S, W = 70.0, 2.0, 344.0          # strip underside height, thickness, width
L_U, L_T = 420.0, 300.0               # uniform section, taper length
H_C, W_C = 3.0, 14.7                  # strip height and width at the connector end
X_TIP = L_U / 2 + L_T                 # 510
POSTS = [(-195.0, 162.0), (195.0, 162.0)]  # half model (y > 0), y = w/2 - 10; mirrored at y < 0
R_POST = 6.0
EPS_PTFE, TAND_PTFE = 2.1, 2.0e-4
SIGMA_AL = 2.5e7
A_CELL = s2.A_X
Z_CELL = H / 2                        # loop axis at mid-height (as in the S2 unit cell)
Y_CELLS = [A_CELL / 2, 3 * A_CELL / 2]  # 35, 105 (half model of the 4-wide row)
Y_CELLS_NY = {4: Y_CELLS, 5: [0.0, A_CELL, 2 * A_CELL]}   # 5 wide: the centre cell is cut by y = 0
BOX = {"open": dict(x=700.0, y=450.0, z=450.0), "closed": dict(x=520.0, y=252.0, z=300.0),
       # reference model: the closed uniform section only (no tapers, posts or launches), wave ports
       # on the cross-sections x = +-150 mm, de-embedded by Palace to the slab faces
       "closedline": dict(x=150.0, y=252.0, z=300.0)}
LINE = "closedline"
R_PORT_HALF = 100.0                   # 2 x 50 ohm in the half model


def zs_metal(attrs, sigma, f=150e6):
    rs = math.sqrt(math.pi * f * s2.MU0 / sigma)
    return {"Attributes": list(attrs), "Rs": 2 * rs, "Ls": 2 * rs / (2 * math.pi * f)}


def h_loc(x):
    s = min(max((abs(x) - L_U / 2) / L_T, 0.0), 1.0)
    return H + (H_C - H) * s


def halfw_loc(x):
    s = min(max((abs(x) - L_U / 2) / L_T, 0.0), 1.0)
    return W / 2 + (W_C / 2 - W / 2) * s


def cell_xs(nx):
    return [] if nx == 0 else ([0.0] if nx == 1 else [-A_CELL / 2, A_CELL / 2])


def build_strip(occ):
    """Strip solid (half, y >= 0): uniform plate plus two sloped trapezoidal plates, 2 mm thick
    (extruded vertically)."""
    parts = [(3, occ.addBox(-L_U / 2, 0, H, L_U, W / 2, T_S))]
    for sg in (+1, -1):
        x0, x1 = sg * L_U / 2, sg * X_TIP
        pts = [occ.addPoint(x0, 0, H), occ.addPoint(x0, W / 2, H),
               occ.addPoint(x1, W_C / 2, H_C), occ.addPoint(x1, 0, H_C)]
        lines = [occ.addLine(pts[i], pts[(i + 1) % 4]) for i in range(4)]
        loop = occ.addCurveLoop(lines)
        face = occ.addPlaneSurface([loop])
        ext = occ.extrude([(2, face)], 0, 0, T_S)
        parts += [e for e in ext if e[0] == 3]
    out, _ = occ.fuse([parts[0]], parts[1:])
    return out


def port_face(occ, x):
    pts = [occ.addPoint(x, 0, 0), occ.addPoint(x, W_C / 2, 0), occ.addPoint(x, W_C / 2, H_C),
           occ.addPoint(x, 0, H_C)]
    lines = [occ.addLine(pts[i], pts[(i + 1) % 4]) for i in range(4)]
    return (2, occ.addPlaneSurface([occ.addCurveLoop(lines)]))


def on_strip(c, tol=0.05):
    """Is the point c = (x, y, z) on the strip surface?"""
    x, y, z = c
    if abs(x) > X_TIP + tol or y > W / 2 + tol:
        return False
    zb = h_loc(x)
    hw = halfw_loc(x)
    if y <= hw + tol and (abs(z - zb) < tol or abs(z - zb - T_S) < tol):
        return True
    if abs(y - hw) < tol * 4 and zb - tol <= z <= zb + T_S + tol:
        return True
    if abs(abs(x) - X_TIP) < tol and H_C - tol <= z <= H_C + T_S + tol:
        return True
    if abs(abs(x) - L_U / 2) < tol and H - tol <= z <= H + T_S + tol:
        return True
    return False


def build_mesh(variant, nx, recipe="v2", path=None, bulk=40.0, strip_div=10.0, box_scale=1.0, box=None,
               grade=40.0):
    import gmsh

    hs = s2.RECIPES[recipe]
    bx = {k: (v * box_scale if variant == "open" else v) for k, v in BOX[variant].items()}
    if box is not None:
        bx = dict(x=box[0], y=box[1], z=box[2])
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("fixture")
    occ = gmsh.model.occ
    box = occ.addBox(-bx["x"], 0, 0, 2 * bx["x"], bx["y"], bx["z"])
    strip = build_strip(occ)
    g = s2.dims("driven", 1, "z", 6.2)       # S2 cell dictionary (s = 6.2 mm, orient z)
    conds, slabs = [], []
    for xc in cell_xs(nx):
        for yc in Y_CELLS:
            c, sl = s2.build_cell(occ, xc, g)
            occ.translate(c + sl, 0, yc, Z_CELL)
            if abs(yc) < 1e-9:     # keep the y >= 0 half of the alumina (the copper is cut with the box)
                clip = occ.addBox(xc - 100, 0, 0, 200, 100, 100)
                sl, _ = occ.intersect(sl, [(3, clip)])
            conds += c
            slabs += sl
    posts = [(3, occ.addCylinder(px, W / 2 - 10.0, 0, 0, 0, H, R_POST)) for px, _ in POSTS]
    ports = [port_face(occ, -X_TIP), port_face(occ, X_TIP)]
    if variant == LINE:
        occ.remove(posts + ports, recursive=True)
        posts, ports = [], []
    dom, _ = occ.cut([(3, box)], strip + conds)
    out, _ = occ.fragment(dom, slabs + posts + ports)
    occ.synchronize()

    vols = [v for (_, v) in gmsh.model.getEntities(3)]
    air, diel, ptfe = [], [], []
    for v in vols:
        bb = gmsh.model.getBoundingBox(3, v)
        if min(bb[3] - bb[0], bb[5] - bb[2]) < g["d"] * 1.01:
            diel.append(v)
        elif (bb[3] - bb[0]) < 2 * R_POST + 0.1 and (bb[4] - bb[1]) < 2 * R_POST + 0.1:
            ptfe.append(v)
        else:
            air.append(v)
    assert len(air) == 1 and len(diel) == len(slabs) and len(ptfe) == len(posts), (air, diel, ptfe)

    b = {v: {abs(t) for (_, t) in gmsh.model.getBoundary([(3, v)], oriented=False)} for v in vols}
    groups = {k: [] for k in range(1, 8)}
    tol = 1e-3
    surfs = [s for (_, s) in gmsh.model.getEntities(2)]
    interface = set()
    for v in diel + ptfe:
        interface |= b[air[0]] & b[v]
    for sf in surfs:
        xmn, ymn, zmn, xmx, ymx, zmx = gmsh.model.getBoundingBox(2, sf)
        com = gmsh.model.occ.getCenterOfMass(2, sf)
        if variant == LINE and abs(xmx - xmn) < tol and abs(abs(xmn) - bx["x"]) < tol:
            groups[5 if xmn < 0 else 6].append(sf)
        elif abs(xmx - xmn) < tol and abs(abs(xmn) - X_TIP) < tol and zmx <= H_C + tol and ymx <= W_C / 2 + tol:
            groups[5 if xmn < 0 else 6].append(sf)
        elif sf in interface:
            continue
        elif abs(zmn) < tol and abs(zmx) < tol:
            groups[1].append(sf)
        elif abs(ymn) < tol and abs(ymx) < tol:
            groups[3].append(sf)
        elif (abs(abs(xmn) - bx["x"]) < tol and abs(xmx - xmn) < tol) or abs(ymn - bx["y"]) < tol or abs(zmn - bx["z"]) < tol:
            groups[4].append(sf)
        elif on_strip(com):
            groups[2].append(sf)
        elif any(abs(sf) in b[v] for v in ptfe):    # post faces against the strip / ground
            if abs(zmn - H) < tol and abs(zmx - H) < tol:
                groups[2].append(sf)
            elif abs(zmn) < tol and abs(zmx) < tol:
                groups[1].append(sf)
        else:
            groups[7].append(sf)
    assert groups[5] and groups[6] and groups[1] and groups[2] and groups[3] and groups[4], groups
    assert bool(groups[7]) == bool(conds), groups[7]

    gmsh.model.addPhysicalGroup(3, air, 1, name="air")
    if diel:
        gmsh.model.addPhysicalGroup(3, diel, 2, name="alumina")
    gmsh.model.addPhysicalGroup(3, ptfe, 3, name="ptfe")
    names = {1: "ground", 2: "strip", 3: "symmetry_y0", 4: "outer", 5: "port1", 6: "port2", 7: "copper"}
    for k, v in groups.items():
        if v:
            gmsh.model.addPhysicalGroup(2, v, k, name=names[k])

    f = gmsh.model.mesh.field
    fields = []

    def threshold(entities, dim, smin, dmin, dmax, sampling, smax=None):
        fd = f.add("Distance")
        f.setNumbers(fd, "SurfacesList" if dim == 2 else "CurvesList", entities)
        f.setNumber(fd, "Sampling", sampling)
        ft = f.add("Threshold")
        f.setNumber(ft, "InField", fd)
        f.setNumber(ft, "SizeMin", smin)
        f.setNumber(ft, "SizeMax", smax or bulk)
        f.setNumber(ft, "DistMin", dmin)
        f.setNumber(ft, "DistMax", dmax)
        if smax is not None:          # S2 near-cell recipe: impose nothing beyond DistMax
            f.setNumber(ft, "StopAtDistMax", 1)
        fields.append(ft)

    # strip: size follows the local strip height (h/strip_div on the faces, half that on the edges),
    # growing at 0.25 per unit distance
    hexpr = f"Max({H_C}, Min({H}, {H} + ({H_C} - {H}) * (Fabs(x) - {L_U / 2}) / {L_T}))"
    strip_curves = sorted({abs(c) for (_, c) in gmsh.model.getBoundary(
        [(2, s) for s in groups[2]], combined=False, oriented=False)})
    for ents, dim, div, minsz in ((groups[2], 2, strip_div, 0.6), (strip_curves, 1, 2 * strip_div, 0.4)):
        fd = f.add("Distance")
        f.setNumbers(fd, "SurfacesList" if dim == 2 else "CurvesList", ents)
        f.setNumber(fd, "Sampling", 400 if dim == 2 else 1500)
        fm = f.add("MathEval")
        f.setString(fm, "F", f"Min({bulk}, Max({minsz}, {hexpr} / {div}) + 0.25 * F{fd})")
        fields.append(fm)
    port_curves = sorted({abs(c) for (_, c) in gmsh.model.getBoundary(
        [(2, s) for s in groups[5] + groups[6]], combined=False, oriented=False)})
    threshold(port_curves, 1, 0.5, 1.0, 40.0, 50)
    post_surfs = sorted(set().union(*[b[v] for v in ptfe]))
    threshold(post_surfs, 2, 3.0, 2.0, 60.0, 60)
    if diel:
        gap_surfs = sorted(set().union(*[b[dv] for dv in diel]))
        gap_curves = sorted({abs(c) for (_, c) in gmsh.model.getBoundary(
            [(2, sfc) for sfc in gap_surfs], combined=False, oriented=False)})
        metal = groups[7]
        metal_curves = sorted({abs(c) for (_, c) in gmsh.model.getBoundary(
            [(2, sfc) for sfc in metal], combined=False, oriented=False)})
        threshold(gap_surfs, 2, hs["gap"], 0.5, 12.0, 150, smax=hs["bulk"])
        threshold(gap_curves, 1, hs["gap_edge"], 0.2, 4.0, 200, smax=hs["bulk"])
        threshold(metal, 2, hs["loop"], 1.0, 15.0, 150, smax=hs["bulk"])
        threshold(metal_curves, 1, hs["medge"], 0.3, 6.0, 400, smax=hs["bulk"])
        # the S2 bulk size (8 mm) around the cells, growing to the fixture bulk
        threshold(metal, 2, hs["bulk"], 15.0, 120.0, 60, smax=grade)
    fmin = f.add("Min")
    f.setNumbers(fmin, "FieldsList", fields)
    f.setAsBackgroundMesh(fmin)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 16)
    gmsh.option.setNumber("Mesh.MeshSizeMax", bulk)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.model.mesh.generate(3)
    _, etags, _ = gmsh.model.mesh.getElements(3)
    info = dict(recipe=recipe, bulk=bulk, grade=grade, strip_div=strip_div, box=bx,
                n_tet=int(sum(len(x) for x in etags)),
                n_surf={names[k]: len(v) for k, v in groups.items()})
    if path is not None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(path))
    gmsh.finalize()
    return info


def probes(nx):
    """Probe points (mm). Empty line: transverse profile at x = 0 and mid-height, plus heights.
    Slab: the centre of every loop bore and the mid-height points between/around cells."""
    pts = []
    if nx == 0:
        for y in [0, 17.5, 35, 52.5, 70, 87.5, 105, 122.5, 140, 157.5, 172, 190, 210, 250]:
            pts.append([0.0, y, Z_CELL])
        for z in [10.0, 60.0]:
            for y in [0.0, 105.0]:
                pts.append([0.0, y, z])
        pts.append([-150.0, 0.0, Z_CELL])
        pts.append([150.0, 0.0, Z_CELL])
    else:
        for xc in cell_xs(nx):
            for yc in Y_CELLS:
                pts.append([xc, yc, Z_CELL])          # loop bore centre
        for y in [0.0, 70.0, 140.0, 172.0]:
            pts.append([cell_xs(nx)[0], y, Z_CELL])
        pts.append([-150.0, 0.0, Z_CELL])
        pts.append([150.0, 0.0, Z_CELL])
    return pts


BANDS = {
    # the 2x-thru: uniform grid from df (P370 needs f = k df) to 300 MHz
    "thru": [{"Type": "Linear", "MinFreq": 0.0005, "MaxFreq": 0.300, "FreqStep": 0.0005}],
    # slab: 50 kHz over 130-200 MHz plus 10 kHz over 146-168 MHz
    "slab": [{"Type": "Linear", "MinFreq": 0.130, "MaxFreq": 0.200, "FreqStep": 0.00005},
             {"Type": "Linear", "MinFreq": 0.146, "MaxFreq": 0.168, "FreqStep": 0.00001}],
    # coarser 2x-thru grid for the large open air box (low-frequency solves are slow there):
    # P370 only needs f = k df; df = 2.5 MHz still gives a 400-ns alias-free time window
    "thru25": [{"Type": "Linear", "MinFreq": 0.0025, "MaxFreq": 0.300, "FreqStep": 0.0025}],
    "test": [{"Type": "Linear", "MinFreq": 0.010, "MaxFreq": 0.300, "FreqStep": 0.010}],
}


def make_config(run, mesh, case, order, metal="al", band="slab", both_ports=False, tol=1e-4,
                max_samples=100, flux=True):
    variant, nx = case["variant"], case["nx"]
    mats = [{"Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0},
            {"Attributes": [3], "Permeability": 1.0, "Permittivity": EPS_PTFE, "LossTan": TAND_PTFE}]
    if nx:
        mats.append({"Attributes": [2], "Permeability": 1.0, "Permittivity": s2.EPS_AL, "LossTan": s2.TAND_AL})
    fixture = [1, 2] + ([4] if variant in ("closed", LINE) else [])
    bnd = {"PMC": {"Attributes": [3]}}
    imp = []
    if metal == "pec":
        bnd["PEC"] = {"Attributes": fixture}
    else:
        imp.append(zs_metal(fixture, SIGMA_AL))
    if nx:
        imp.append(s2.zs_entry([7], s2.F_RS, "exact"))
    if imp:
        bnd["Impedance"] = imp
    if variant == "open":
        bnd["Absorbing"] = {"Attributes": [4], "Order": 1}
        if flux:
            bnd["Postprocessing"] = {"SurfaceFlux": [{"Index": 1, "Attributes": [4], "Type": "Power",
                                                      "Center": [0.0, 0.0, 0.0]}]}
    if variant == LINE:
        off = BOX[LINE]["x"] - nx * A_CELL / 2         # port plane -> slab face
        bnd["WavePort"] = [
            {"Index": 1, "Attributes": [5], "Mode": 1, "Offset": off, "Excitation": 1},
            {"Index": 2, "Attributes": [6], "Mode": 1, "Offset": off, "Excitation": 2 if both_ports else 0}]
    else:
        bnd["LumpedPort"] = [
            {"Index": 1, "R": R_PORT_HALF, "Excitation": 1, "Attributes": [5], "Direction": "+Z"},
            {"Index": 2, "R": R_PORT_HALF, "Excitation": 2 if both_ports else 0, "Attributes": [6], "Direction": "+Z"},
        ]
    dom_pp = {"Probe": [{"Index": i + 1, "Center": p} for i, p in enumerate(case["probes"])]}
    if nx:
        dom_pp["Energy"] = [{"Index": 1, "Attributes": [2]}, {"Index": 2, "Attributes": [1]}]
    cfg = {
        "Problem": {"Type": "Driven", "Verbose": 2, "Output": f"postpro/{run}",
                    "OutputFormats": {"Paraview": False, "GridFunction": False}},
        "Model": {"Mesh": f"mesh/{mesh}.msh", "L0": 1.0e-3},
        "Domains": {"Materials": mats, "Postprocessing": dom_pp},
        "Boundaries": bnd,
        "Solver": {"Order": order, "Device": "CPU",
                   "Linear": {"Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-8, "MaxIts": 500},
                   "Driven": {"Samples": BANDS[band], "Save": [], "AdaptiveTol": tol,
                              "AdaptiveMaxSamples": max_samples}},
    }
    return cfg


# ---- 2D cross-section (boundary-mode analysis: Z0, higher-order modes) ----
def build_xsec(variant, path, d_strip=1.0, bulk=25.0):
    """2D mesh of the uniform-section cross-section (mesh coordinates: X = y, Y = z), half plane.
    Attributes: domain 1 air; curves 1 ground, 2 strip, 3 symmetry line y = 0, 4 outer."""
    import gmsh
    bx = BOX[variant]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("xsec")
    occ = gmsh.model.occ
    dom = occ.addRectangle(0, 0, 0, bx["y"], bx["z"])
    st = occ.addRectangle(0, H, 0, W / 2, T_S)
    out, _ = occ.cut([(2, dom)], [(2, st)])
    occ.synchronize()
    tol = 1e-6
    grp = {1: [], 2: [], 3: [], 4: []}
    for (_, c) in gmsh.model.getEntities(1):
        xmn, ymn, _, xmx, ymx, _ = gmsh.model.getBoundingBox(1, c)
        if abs(ymn) < tol and abs(ymx) < tol:
            grp[1].append(c)
        elif abs(xmn) < tol and abs(xmx) < tol:
            grp[3].append(c)
        elif abs(xmn - bx["y"]) < tol or abs(ymn - bx["z"]) < tol:
            grp[4].append(c)
        else:
            grp[2].append(c)
    gmsh.model.addPhysicalGroup(2, [out[0][1]], 1, name="air")
    for k, v in grp.items():
        gmsh.model.addPhysicalGroup(1, v, k)
    f = gmsh.model.mesh.field
    fd = f.add("Distance")
    f.setNumbers(fd, "CurvesList", grp[2])
    f.setNumber(fd, "Sampling", 2000)
    ft = f.add("Threshold")
    f.setNumber(ft, "InField", fd)
    f.setNumber(ft, "SizeMin", d_strip)
    f.setNumber(ft, "SizeMax", bulk)
    f.setNumber(ft, "DistMin", 1.0)
    f.setNumber(ft, "DistMax", 150.0)
    f.setAsBackgroundMesh(ft)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.model.mesh.generate(2)
    _, et, _ = gmsh.model.mesh.getElements(2)
    n = int(sum(len(x) for x in et))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    gmsh.write(str(path))
    gmsh.finalize()
    return n


def make_xconfig(run, variant, sym, freq, n_modes=8, order=2, metal="pec"):
    bx = BOX[variant]
    bnd = {}
    pec = [1, 2] + ([4] if variant == "closed" else [])
    bnd["PEC"] = {"Attributes": pec + ([3] if sym == "odd" else [])}
    if sym == "even":
        bnd["PMC"] = {"Attributes": [3]}
    if variant == "open":
        bnd["Absorbing"] = {"Attributes": [4], "Order": 1}
    yv = H / 2
    bnd["Postprocessing"] = {"Impedance": [{
        "Index": 1, "VoltagePath": [[0.0, 0.0], [0.0, H]], "NSamples": 200,
        "CurrentPath": [[0.0, H - 10], [W / 2 + 15, H - 10], [W / 2 + 15, H + T_S + 10], [0.0, H + T_S + 10]]}]}
    cfg = {
        "Problem": {"Type": "BoundaryMode", "Verbose": 2, "Output": f"postpro/{run}",
                    "OutputFormats": {"Paraview": False, "GridFunction": False}},
        "Model": {"Mesh": f"mesh/xs_{variant}.msh", "L0": 1.0e-3},
        "Domains": {"Materials": [{"Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0}]},
        "Boundaries": bnd,
        "Solver": {"Order": order, "Device": "CPU",
                   "BoundaryMode": {"Freq": freq, "N": n_modes, "Save": 0, "Tol": 1e-8, "Target": 0.999},
                   "Linear": {"Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-10, "MaxIts": 400}},
    }
    _ = yv, bx
    return cfg


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("mesh")
    m.add_argument("--variant", choices=["open", "closed", LINE], required=True)
    m.add_argument("--nx", type=int, choices=[0, 1, 2], required=True)
    m.add_argument("--recipe", default="v2", choices=sorted(s2.RECIPES))
    m.add_argument("--bulk", type=float, default=40.0)
    m.add_argument("--strip-div", type=float, default=10.0)
    m.add_argument("--box-scale", type=float, default=1.0, help="open variant: scale the air box")
    m.add_argument("--ny", type=int, choices=[4, 5], default=4, help="cells across the strip")
    m.add_argument("--w", type=float, default=344.0, help="strip width in the uniform section (mm)")
    m.add_argument("--box", type=float, nargs=3, default=None, help="air box half-length x, y, z (mm)")
    m.add_argument("--name", default=None)
    c = sub.add_parser("config")
    c.add_argument("--mesh", required=True)
    c.add_argument("--name", required=True)
    c.add_argument("--order", type=int, default=2)
    c.add_argument("--metal", choices=["al", "pec"], default="al")
    c.add_argument("--band", choices=sorted(BANDS), default="slab")
    c.add_argument("--both-ports", action="store_true")
    c.add_argument("--tol", type=float, default=1e-4)
    x = sub.add_parser("xsec")
    x.add_argument("--variant", choices=["open", "closed"], required=True)
    xc = sub.add_parser("xconfig")
    xc.add_argument("--variant", choices=["open", "closed"], required=True)
    xc.add_argument("--sym", choices=["even", "odd"], required=True)
    xc.add_argument("--freq", type=float, required=True, help="GHz")
    xc.add_argument("--name", required=True)
    xc.add_argument("--nmodes", type=int, default=8)
    args = ap.parse_args()

    if args.cmd == "mesh":
        global W, Y_CELLS
        W = args.w
        Y_CELLS = Y_CELLS_NY[args.ny]
        name = args.name or f"{args.variant}_nx{args.nx}"
        info = build_mesh(args.variant, args.nx, args.recipe, HERE / "mesh" / f"{name}.msh",
                          args.bulk, args.strip_div, args.box_scale, args.box)
        meta = dict(name=name, variant=args.variant, nx=args.nx, ny=args.ny,
                    cell_centres_half=[[xc_, yc, Z_CELL] for xc_ in cell_xs(args.nx) for yc in Y_CELLS],
                    slab_length_mm=args.nx * A_CELL, slab_faces_x_mm=[-args.nx * A_CELL / 2, args.nx * A_CELL / 2],
                    fixture=dict(h=H, t_strip=T_S, w=W, L_uniform=L_U, L_taper=L_T, h_tip=H_C, w_tip=W_C,
                                 x_tip=X_TIP, posts=[(px, W / 2 - 10.0) for px, _ in POSTS], r_post=R_POST, sigma_al=SIGMA_AL,
                                 eps_ptfe=EPS_PTFE, tand_ptfe=TAND_PTFE, port_R_half=R_PORT_HALF),
                    cell=dict(s2_dims=s2.dims("driven", 1, "z", 6.2), z_axis=Z_CELL),
                    probes=probes(args.nx), mesh=info)
        (HERE / "cases").mkdir(exist_ok=True)
        (HERE / "cases" / f"{name}.json").write_text(json.dumps(meta, indent=2) + "\n")
        print(json.dumps(dict(name=name, n_tet=info["n_tet"], n_surf=info["n_surf"])))
    elif args.cmd == "config":
        case = json.loads((HERE / "cases" / f"{args.mesh}.json").read_text())
        cfg = make_config(args.name, args.mesh, case, args.order, args.metal, args.band, args.both_ports, args.tol)
        (HERE / "configs").mkdir(exist_ok=True)
        (HERE / "configs" / f"{args.name}.json").write_bytes((json.dumps(cfg, indent=2) + "\n").encode())
        meta = dict(run=args.name, mesh=args.mesh, type="driven", order=args.order, metal=args.metal,
                    band=args.band, both_ports=args.both_ports, variant=case["variant"], nx=case["nx"],
                    probes=case["probes"], impedance=cfg["Boundaries"].get("Impedance"))
        (HERE / "runs").mkdir(exist_ok=True)
        (HERE / "runs" / f"{args.name}.json").write_bytes((json.dumps(meta, indent=2) + "\n").encode())
        print(args.name)
    elif args.cmd == "xsec":
        n = build_xsec(args.variant, HERE / "mesh" / f"xs_{args.variant}.msh")
        print(f"xs_{args.variant}: {n} triangles")
    else:
        cfg = make_xconfig(args.name, args.variant, args.sym, args.freq, args.nmodes)
        (HERE / "configs").mkdir(exist_ok=True)
        (HERE / "configs" / f"{args.name}.json").write_bytes((json.dumps(cfg, indent=2) + "\n").encode())
        meta = dict(run=args.name, type="boundarymode", variant=args.variant, sym=args.sym, freq_GHz=args.freq)
        (HERE / "runs").mkdir(exist_ok=True)
        (HERE / "runs" / f"{args.name}.json").write_bytes((json.dumps(meta, indent=2) + "\n").encode())
        print(args.name)


if __name__ == "__main__":
    main()
