"""P3 drum-inductor FEM (Palace magnetostatic): meshes and configs.

Goal: predict the coupling between lumped LC cells built from unshielded ferrite drum
inductors, before the coupon boards are measured, and test where the point-dipole
lattice model breaks down at short spacing.

Model of one cell (units mm):
- Ferrite drum core: bottom flange, centre post, top flange (one volume, relative
  permeability MU swept because the vendors do not publish it).
- Winding: K representative turns, thin PEC annuli (width W_ANN) in the winding window,
  each closed by a small surface-current port across a gap. Each turn is its own terminal;
  the series winding is the sum over a cell's block of the inductance matrix. An N-turn
  winding scales every entry by (N/K)^2, so the coupling k = M_ab / L_aa between two cells
  does not depend on N.
- Far boundary: PEC sphere.

Geometry (diameter of flange D, height H, flange thickness TF, post diameter DP) is an
ESTIMATE from datasheet outlines; the post diameter is not published. Re-run with measured
core dimensions when the parts arrive.

Runs written here (configs/<run>.json, mesh/<mesh>.msh, lists/<round>.txt):
  val_*    validation: two thin rings in air, coaxial and coplanar, vs Neumann (postprocess)
  single_* one cell, per type and permeability (self block, K convergence)
  pair_*   two cells coplanar at d = 11.5..50 mm, or coaxial (stacked) at several gaps

Run with the project venv:  D:/PhD/.venv/Scripts/python.exe drum_geometry.py
"""
import json
import math
import os

import gmsh

HERE = os.path.dirname(os.path.abspath(__file__))
W_ANN = 0.4      # annulus (representative turn) radial width, mm
GAP = 0.3        # port gap length along the circumference, mm
R_AIR = 200.0    # far-boundary sphere radius, mm

# Estimated core geometry (mm). DP (post) is a guess; datasheet outlines give D and H only.
TYPES = {
    "A": dict(D=9.6, H=5.4, TF=0.9, DP=5.0, L_uH=47.0, DCR=0.170, name="Sumida CD105T125NP-470KC"),
    "B": dict(D=7.4, H=5.0, TF=0.8, DP=3.6, L_uH=100.0, DCR=0.273, name="SXN SM7850-101MT"),
    "C": dict(D=7.6, H=5.3, TF=0.8, DP=3.6, L_uH=100.0, DCR=0.450, name="Bourns SDR0805-101KL"),
}


def turns_layout(t, kr, kz):
    """Representative turn radii and heights inside the winding window."""
    r0 = t["DP"] / 2 + 0.25 + W_ANN / 2
    r1 = t["D"] / 2 - 0.35 - W_ANN / 2
    z0 = t["TF"] + 0.3
    z1 = t["H"] - t["TF"] - 0.3
    rs = [r0 + (r1 - r0) * (i + 0.5) / kr for i in range(kr)] if kr > 1 else [(r0 + r1) / 2]
    zs = [z0 + (z1 - z0) * (j + 0.5) / kz for j in range(kz)] if kz > 1 else [(z0 + z1) / 2]
    return [(r, z) for z in zs for r in rs]


class Model:
    def __init__(self, name):
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(name)
        self.name = name
        self.occ = gmsh.model.occ
        self.ferrite = []         # volume tags
        self.rings = []           # list of (pec_surface_tag, port_surface_tag, direction)
        self.fine_curves = []

    def drum(self, t, xc, yc, zb):
        o = self.occ
        f1 = o.addCylinder(xc, yc, zb, 0, 0, t["TF"], t["D"] / 2)
        p = o.addCylinder(xc, yc, zb + t["TF"], 0, 0, t["H"] - 2 * t["TF"], t["DP"] / 2)
        f2 = o.addCylinder(xc, yc, zb + t["H"] - t["TF"], 0, 0, t["TF"], t["D"] / 2)
        v, _ = o.fuse([(3, f1)], [(3, p), (3, f2)])
        self.ferrite += [tag for dim, tag in v]

    def ring(self, xc, yc, z, r):
        """Annulus of radius r (width W_ANN) in the plane z, split into a PEC part and a
        small port piece at angle 0 (the +x side). Current direction at the port: +y (CCW from +z)."""
        o = self.occ
        outer = o.addDisk(xc, yc, z, r + W_ANN / 2, r + W_ANN / 2)
        inner = o.addDisk(xc, yc, z, r - W_ANN / 2, r - W_ANN / 2)
        ann, _ = o.cut([(2, outer)], [(2, inner)])
        cutter = o.addRectangle(xc + r - W_ANN, yc - GAP / 2, z, 2 * W_ANN, GAP)
        port, _ = o.intersect(ann, [(2, cutter)], removeObject=False, removeTool=True)
        pec, _ = o.cut(ann, port, removeObject=True, removeTool=False)
        self.rings.append(dict(pec=[t for _, t in pec], port=[t for _, t in port], dir=[0.0, 1.0, 0.0]))

    def cell(self, t, xc, yc, zb, kr, kz, with_core=True):
        if with_core:
            self.drum(t, xc, yc, zb)
        for r, z in turns_layout(t, kr, kz):
            self.ring(xc, yc, zb + z, r)

    def build(self, h_fine=0.15, h_core=0.5, h_far=20.0):
        o = self.occ
        air = o.addSphere(0, 0, 0, R_AIR)
        objs = [(3, air)] + [(3, v) for v in self.ferrite]
        tools = []
        for rg in self.rings:
            tools += [(2, s) for s in rg["pec"] + rg["port"]]
        out, outmap = o.fragment(objs, tools)
        o.synchronize()
        # map original entities to fragment results
        k = 0
        new_ferrite = []
        for i, (dim, tag) in enumerate(objs):
            if i > 0:
                new_ferrite += [t for d, t in outmap[i] if d == 3]
        ferrite_set = set(new_ferrite)
        idx = len(objs)
        for rg in self.rings:
            npec, nport = [], []
            for s in rg["pec"]:
                npec += [t for d, t in outmap[idx] if d == 2]
                idx += 1
            for s in rg["port"]:
                nport += [t for d, t in outmap[idx] if d == 2]
                idx += 1
            rg["pec"], rg["port"] = npec, nport
        vols = [t for d, t in gmsh.model.getEntities(3)]
        air_vols = [v for v in vols if v not in ferrite_set]
        gmsh.model.addPhysicalGroup(3, air_vols, 1, "air")
        if new_ferrite:
            gmsh.model.addPhysicalGroup(3, new_ferrite, 2, "ferrite")
        # far boundary: surfaces of the air sphere not shared with anything else
        bnd = gmsh.model.getBoundary([(3, v) for v in vols], oriented=False, combined=True)
        far = [t for d, t in bnd if d == 2]
        gmsh.model.addPhysicalGroup(2, far, 3, "farfield")
        pec_all = [t for rg in self.rings for t in rg["pec"]]
        gmsh.model.addPhysicalGroup(2, pec_all, 4, "rings_pec")
        for i, rg in enumerate(self.rings):
            gmsh.model.addPhysicalGroup(2, rg["port"], 10 + i, f"port{i + 1}")
        # mesh size: fine on ring edges, medium in ferrite, coarse far away
        ring_curves = []
        for rg in self.rings:
            for s in rg["pec"] + rg["port"]:
                ring_curves += [t for d, t in gmsh.model.getBoundary([(2, s)], oriented=False)]
        core_surfs = [t for d, t in gmsh.model.getBoundary([(3, v) for v in new_ferrite], oriented=False)] if new_ferrite else []
        f = gmsh.model.mesh.field
        d1 = f.add("Distance")
        f.setNumbers(d1, "CurvesList", list(set(ring_curves)))
        f.setNumber(d1, "Sampling", 60)
        t1 = f.add("Threshold")
        f.setNumber(t1, "InField", d1)
        f.setNumber(t1, "SizeMin", h_fine)
        f.setNumber(t1, "SizeMax", h_far)
        f.setNumber(t1, "DistMin", 0.3)
        f.setNumber(t1, "DistMax", 60.0)
        fields = [t1]
        if core_surfs:
            d2 = f.add("Distance")
            f.setNumbers(d2, "SurfacesList", core_surfs)
            t2 = f.add("Threshold")
            f.setNumber(t2, "InField", d2)
            f.setNumber(t2, "SizeMin", h_core)
            f.setNumber(t2, "SizeMax", h_far)
            f.setNumber(t2, "DistMin", 0.5)
            f.setNumber(t2, "DistMax", 60.0)
            fields.append(t2)
        mn = f.add("Min")
        f.setNumbers(mn, "FieldsList", fields)
        f.setAsBackgroundMesh(mn)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
        gmsh.model.mesh.generate(3)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        os.makedirs(os.path.join(HERE, "mesh"), exist_ok=True)
        path = os.path.join(HERE, "mesh", self.name + ".msh")
        gmsh.write(path)
        ntet = sum(len(t) for t in gmsh.model.mesh.getElementsByType(4)[0:1])
        gmsh.finalize()
        return path, len(self.rings), ntet


def config(run, mesh, nports, mu):
    cfg = {
        "Problem": {"Type": "Magnetostatic", "Verbose": 2, "Output": f"postpro/{run}"},
        "Model": {"Mesh": f"mesh/{mesh}.msh", "L0": 1.0e-3},
        "Domains": {"Materials": [{"Attributes": [1], "Permeability": 1.0}]
                    + ([{"Attributes": [2], "Permeability": mu}] if mu is not None else [])},
        "Boundaries": {
            "PEC": {"Attributes": [3, 4]},
            "SurfaceCurrent": [{"Index": i + 1, "Attributes": [10 + i], "Direction": [0.0, 1.0, 0.0]}
                               for i in range(nports)],
        },
        "Solver": {"Order": 2, "Device": "CPU",
                   "Magnetostatic": {"Save": 0},
                   "Linear": {"Type": "AMS", "KSPType": "CG", "Tol": 1.0e-8, "MaxIts": 300}},
    }
    os.makedirs(os.path.join(HERE, "configs"), exist_ok=True)
    with open(os.path.join(HERE, "configs", run + ".json"), "w") as f:
        json.dump(cfg, f, indent=1)


def write_list(name, runs):
    os.makedirs(os.path.join(HERE, "lists"), exist_ok=True)
    # LF line endings: bash on the cluster reads these (CRLF breaks the run names)
    with open(os.path.join(HERE, "lists", name + ".txt"), "w", newline="\n") as f:
        f.write("\n".join(runs) + "\n")


def main(rounds=("val", "single", "pair")):
    meta = {}
    runs_by_round = {r: [] for r in rounds}
    if "val" in rounds:
        # two thin rings (r = 5 mm) in air: coaxial gap 5 mm, coplanar centre distance 12 mm
        for kind, dx, dz in (("coax5", 0.0, 5.0), ("copl12", 12.0, 0.0)):
            name = f"val_{kind}"
            m = Model(name)
            m.ring(-dx / 2, 0, -dz / 2, 5.0)
            m.ring(dx / 2, 0, dz / 2, 5.0)
            path, n, ntet = m.build()
            config(name, name, n, None)
            meta[name] = dict(kind="val", r=5.0, w=W_ANN, dx=dx, dz=dz, ports=n, tets=ntet)
            runs_by_round["val"].append(name)
    if "single" in rounds:
        for typ, t in TYPES.items():
            for kr, kz in ((1, 2), (2, 4)):
                mesh = f"single_{typ}_k{kr}x{kz}"
                m = Model(mesh)
                m.cell(t, 0, 0, -t["H"] / 2, kr, kz)
                path, n, ntet = m.build()
                for mu in (50, 200, 1000):
                    run = f"{mesh}_mu{mu}"
                    config(run, mesh, n, mu)
                    meta[run] = dict(kind="single", type=typ, kr=kr, kz=kz, mu=mu, ports=n, tets=ntet, **t)
                    runs_by_round["single"].append(run)
            # coreless reference (air core), K = 2x4
            mesh = f"single_{typ}_k2x4"
            run = f"{mesh}_air"
            config(run, mesh, 8, 1.0)
            meta[run] = dict(kind="single", type=typ, kr=2, kz=4, mu=1, ports=8, **t)
            runs_by_round["single"].append(run)
    if "pair" in rounds:
        kr, kz = 2, 4
        for typ, t in TYPES.items():
            geoms = [("copl", d) for d in (11.5, 14.0, 18.0, 25.0, 35.0, 50.0)] + \
                    [("coax", g) for g in (1.0, 5.0, 15.0)]
            for kind, s in geoms:
                mesh = f"pair_{typ}_{kind}{s:g}"
                m = Model(mesh)
                if kind == "copl":
                    m.cell(t, -s / 2, 0, -t["H"] / 2, kr, kz)
                    m.cell(t, s / 2, 0, -t["H"] / 2, kr, kz)
                else:   # stacked: gap s between the top of the lower drum and the bottom of the upper one
                    m.cell(t, 0, 0, -t["H"] - s / 2, kr, kz)
                    m.cell(t, 0, 0, s / 2, kr, kz)
                path, n, ntet = m.build()
                for mu in (200,):
                    run = f"{mesh}_mu{mu}"
                    config(run, mesh, n, mu)
                    meta[run] = dict(kind="pair", geom=kind, s=s, type=typ, kr=kr, kz=kz, mu=mu,
                                     ports=n, tets=ntet, **t)
                    runs_by_round["pair"].append(run)
    if "array" in rounds:
        # the 4 x 4 mini-arrays of board 2 (A at 15 mm, B/C at 12.5 mm), K = 1 x 2 turns per cell
        # (32 terminals); the full 32 x 32 matrix feeds the coupled-circuit model of the array
        kr, kz = 1, 2
        for typ, p in (("A", 15.0), ("B", 12.5), ("C", 12.5)):
            t = TYPES[typ]
            mesh = f"array_{typ}_4x4_p{p:g}"
            m = Model(mesh)
            for i in range(4):
                for j in range(4):
                    m.cell(t, (i - 1.5) * p, (j - 1.5) * p, -t["H"] / 2, kr, kz)
            path, n, ntet = m.build(h_fine=0.2, h_core=0.6)
            run = f"{mesh}_mu200"
            config(run, mesh, n, 200)
            meta[run] = dict(kind="array", type=typ, pitch=p, kr=kr, kz=kz, mu=200, ports=n, tets=ntet, **t)
            runs_by_round["array"].append(run)
    for r, runs in runs_by_round.items():
        write_list(r, runs)
    old = {}
    p = os.path.join(HERE, "runs_meta.json")
    if os.path.exists(p):
        old = json.load(open(p))
    old.update(meta)
    json.dump(old, open(p, "w"), indent=1)
    for r, runs in runs_by_round.items():
        print(r, len(runs), "runs")


if __name__ == "__main__":
    import sys
    main(tuple(sys.argv[1:]) or ("val", "single", "pair"))
