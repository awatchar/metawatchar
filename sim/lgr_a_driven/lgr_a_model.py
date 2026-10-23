#!/usr/bin/env python3
"""LGR-A meta-atom in its copper can, with SMA coupling loops, brass tuning screw and lip fillets:
geometry, Gmsh mesh and Palace configuration (eigenmode with lumped ports, or 2-port driven).

Cell (fixed, from the size sweep in ../lgr_size_sweep): V = 0.5 L, a = r/R = 0.5, alumina gap
d = 1.0 mm, lip height s = 10.42 mm, loop wall t = 3 mm, lip thickness 3 mm, can H = 2R,
end clearance c = 0.25 R. All lengths in mm (Palace L0 = 1e-3). Coordinates: can axis = z,
mid-height z = 0, gap and lips on +x (0 deg).

Additions in this model:
  * PTFE spacer rings (optional, eps_r 2.1, tan_delta 2e-4): radius r..r+t, between each loop end
    and the can end (height c), as in the CAD (hardware/lgr_a).
  * Two SMA ports at 120 and 240 deg, mid-height. Each is a short 50 ohm air coax recess in the can
    wall (inner conductor = the 1.0 mm loop wire, outer diameter 2.3 mm, 60 ln(2.3) = 49.97 ohm)
    ending in a coaxial lumped port (annulus, 50 ohm, "Direction": "+R") 1 mm behind the inner
    wall surface. The wire runs from the port into the can and forms a rectangular loop that is
    grounded on the can wall. Orientation:
      "xy": loop in the mid-height plane (normal along the can axis), so it links the axial
            return flux in the annulus between the LGR and the can. Radial leg of depth w,
            tangential leg of width h (pointing away from the gap, towards 180 deg), radial
            return leg grounded on the wall. Ports are mirror images about the x-z plane.
      "rz": loop in the plane containing the radial direction and the can axis (radial leg w,
            axial leg h upwards, radial return leg to the wall). Normal is azimuthal.
    w and h are wire-centreline dimensions; wire diameter 1.0 mm.
  * Brass M4 tuning screw at 0 deg, mid-height, facing the lip tips: a 4.0 mm diameter cylinder
    whose flat end is `screw` mm inside the can inner wall (0 = no screw, plain wall).
  * Exploratory alternative tuner (`axial`): a brass rod on the can axis through the lid, reaching
    into the loop bore, where the magnetic field is strongest (inductive tuning).
  * Lip-tip fillets: radius `fillet` on the four axial edges at the lip tips (as in the CAD).
  * PEEK slab clips (optional, `clips=True`; eps_r 3.2, tan_delta 3e-3): the U-clips of
    hardware/lgr_a/lgr_a_assembly.py `slab_clip`, one at each axial end of the lips. Each covers
    the outer 6 mm of both lips (2 mm side walls outside the lips, 2 mm back wall beyond the lip
    tips) over 4 mm along z, flush with the loop end. The 0.05 mm interference is not modelled
    (the opening equals the lip pack, 2 t_lip + d).
Conductors are not meshed; they are impedance boundaries (copper sigma 5.8e7, brass 1.5e7 S/m).

Copper surface impedance (corrected 2026-10-23). Palace's "Impedance" boundary combines Rs, Ls
(and Cs) IN PARALLEL: 1/Z = 1/Rs + 1/(j w Ls). A good conductor has Z = Rs (1 + j), so the
correct input is Rs_param = 2 Rs, Ls_param = 2 Rs / w, with Rs = sqrt(pi f mu0 / sigma) at the
operating frequency (zs = "exact", case names ending in "_cu2"). The first runs of this study used
Rs, Ls = Rs / w (zs = "s2"), which is Z = Rs (1 + j) / 2, i.e. half the copper loss: their Q_c and
Q values are about 2x too high (frequencies are unaffected). The empty-can case `e_emptycan_cu2_o2`
checks the input against the analytic TM010 copper Q.

Mesh attributes
  volumes   1 air, 2 alumina, 3 PTFE, 4 PEEK clips
  surfaces  1 can (shield) inner walls, 2 LGR loop + lips, 3 coupling wires + coax recess walls,
            4 brass screw, 11 port 1 (120 deg), 12 port 2 (240 deg)

Usage: see build_cases.py (it defines every case in this study).
"""
import copy
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

# ---- cell (note: identical to lgr_size_sweep case v050_a05_d1000_s1042) ----
V_L, A_RATIO, D_GAP, S_LIP = 0.5, 0.5, 1.0, 10.42
T_WALL, T_LIP = 3.0, 3.0
F0 = 145.33e6
MU0 = 4e-7 * math.pi
SIGMA_CU, SIGMA_BRASS = 5.8e7, 1.5e7


def rs(sigma, f=F0):
    return math.sqrt(math.pi * f * MU0 / sigma)


RS_CU, RS_BR = rs(SIGMA_CU), rs(SIGMA_BRASS)
LS_CU, LS_BR = RS_CU / (2 * math.pi * F0), RS_BR / (2 * math.pi * F0)
EPS_AL, TAND_AL = 9.8, 4.0e-5
EPS_PTFE, TAND_PTFE = 2.1, 2.0e-4
EPS_PEEK, TAND_PEEK = 3.2, 3.0e-3
CLIP = dict(wall=2.0, depth=6.0, length=4.0)   # mm, as hardware/lgr_a/lgr_a_assembly.py CLIP


def zs_entry(attrs, sigma, f, zs="exact"):
    """Palace "Impedance" entry for a good conductor of conductivity sigma at frequency f.
    Palace combines Rs and Ls IN PARALLEL (1/Z = 1/Rs + 1/(j w Ls)).
      zs = "exact": Rs_p = 2 Rs, Ls_p = 2 Rs / w  ->  Z = Rs (1 + j) at f (correct)
      zs = "s2"   : Rs, Ls = Rs / w (input of the first runs)  ->  Z = Rs (1 + j) / 2 (half the loss)
    """
    r, w = rs(sigma, f), 2 * math.pi * f
    k = 2.0 if zs == "exact" else 1.0
    return {"Attributes": list(attrs), "Rs": k * r, "Ls": k * r / w}

# ---- ports ----
PORT_ANG = (120.0, 240.0)
R_WIRE = 0.5            # 1.0 mm wire
B_COAX = 1.15           # outer radius of the 50 ohm air coax recess (60 ln(1.15/0.5) = 49.97 ohm)
L_COAX = 1.0            # recess depth behind the inner wall; port annulus sits here
Z_PORT = 50.0
SCREW_R = 2.0           # M4


def dims(s=S_LIP):
    R = (V_L * 1e-3 / (2 * math.pi)) ** (1 / 3) * 1e3
    H = 2 * R
    r = A_RATIO * R
    c = 0.25 * R
    z = H - 2 * c
    return dict(V_l=V_L, a=A_RATIO, R=R, H=H, r=r, t=T_WALL, z=z, c=c, s=s, d=D_GAP,
                t_lip=T_LIP, clearance=R - (r + T_WALL + s))


def port_frame(ang_deg, sign):
    a = math.radians(ang_deg)
    u = (math.cos(a), math.sin(a), 0.0)               # radial, outward
    v = (-sign * math.sin(a), sign * math.cos(a), 0.0)  # tangential (sign +1: increasing phi)
    return u, v


def add(p, q, k=1.0):
    return tuple(pi + k * qi for pi, qi in zip(p, q))


def scale(p, k):
    return tuple(k * pi for pi in p)


def loop_points(g, ang, sign, w, h, orient):
    """Centreline points of the coupling wire: port end (beyond the recess end), corner 1,
    corner 2, ground end (beyond the inner wall)."""
    R = g["R"]
    u, v = port_frame(ang, sign)
    e2 = v if orient == "xy" else (0.0, 0.0, 1.0)
    p0 = scale(u, R + L_COAX + 0.6)
    p1 = scale(u, R - w)
    p2 = add(p1, e2, h)
    p3 = add(scale(u, R + 1.5), e2, h)
    return [p0, p1, p2, p3]


def build_empty_can(g, path, h=1.0):
    """Empty copper can (no LGR): check of the impedance input against the analytic TM010 Q."""
    import gmsh
    R, H = g["R"], g["H"]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("can")
    occ = gmsh.model.occ
    v = occ.addCylinder(0, 0, -H / 2, 0, 0, H, R)
    occ.synchronize()
    gmsh.model.addPhysicalGroup(3, [v], 1, name="air")
    gmsh.model.addPhysicalGroup(2, [abs(t) for (_, t) in gmsh.model.getBoundary([(3, v)], oriented=False)],
                                1, name="can")
    gmsh.option.setNumber("Mesh.MeshSizeMax", 6.0 * h)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.model.mesh.generate(3)
    n_tet = sum(len(x) for x in gmsh.model.mesh.getElements(3)[1])
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    gmsh.write(str(path))
    gmsh.finalize()
    return dict(n_tet=int(n_tet), h=h), []


def empty_can_tm010(g, sigma=SIGMA_CU):
    """Analytic TM010 of the empty can: f = c p01 / (2 pi R), Q_c = R H / (delta (R + H))."""
    c0, p01 = 299792458.0, 2.404825557695773
    R, H = g["R"] * 1e-3, g["H"] * 1e-3
    f = c0 * p01 / (2 * math.pi * R)
    delta = 1.0 / math.sqrt(math.pi * f * MU0 * sigma)
    return dict(f_Hz=f, Q_c=R * H / (delta * (R + H)))


def build_mesh(g, path, loops=None, spacers=True, screw=0.0, fillet=0.0, h=1.0, hl=1.0,
               axial=None, clips=False, gui=False):
    """loops: None or dict(w=, h=, orient='xy'|'rz'). h = global size scale (as in the sweep),
    hl = extra scale on the coupling-loop / port region. clips: PEEK slab clips (see module doc)."""
    import gmsh

    R, H, r, t, z, s, d, tl = g["R"], g["H"], g["r"], g["t"], g["z"], g["s"], g["d"], g["t_lip"]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("lgr_a")
    occ = gmsh.model.occ

    # ---- LGR conductor: tube with slit + lips (+ fillets) ----
    tube_o = occ.addCylinder(0, 0, -z / 2, 0, 0, z, r + t)
    tube_i = occ.addCylinder(0, 0, -z / 2 - 1, 0, 0, z + 2, r)
    tube, _ = occ.cut([(3, tube_o)], [(3, tube_i)])
    slit = occ.addBox(r - 2, -d / 2, -z / 2 - 1, t + 4, d, z + 2)
    tube, _ = occ.cut(tube, [(3, slit)])
    x0 = r + t / 2
    lip1 = occ.addBox(x0, d / 2, -z / 2, r + t + s - x0, tl, z)
    lip2 = occ.addBox(x0, -d / 2 - tl, -z / 2, r + t + s - x0, tl, z)
    cond, _ = occ.fuse(tube, [(3, lip1), (3, lip2)])
    occ.synchronize()
    xt = r + t + s
    if fillet > 0:
        cv = []
        curves = {abs(cc) for (dim, cc) in gmsh.model.getEntities(1)}
        for cc in curves:
            bb = gmsh.model.getBoundingBox(1, cc)
            if (abs(bb[0] - xt) < 1e-4 and abs(bb[3] - xt) < 1e-4 and bb[5] - bb[2] > z - 1e-3
                    and abs(bb[1] - bb[4]) < 1e-6
                    and min(abs(abs(bb[1]) - d / 2), abs(abs(bb[1]) - (d / 2 + tl))) < 1e-4):
                cv.append(cc)
        assert len(cv) == 4, cv
        cond = occ.fillet([cond[0][1]], cv, [fillet], removeVolume=True)
        cond = [(3, cond[0][1])] if cond and cond[0][0] == 3 else cond
        occ.synchronize()

    # ---- air domain: can + coax recesses ----
    can = occ.addCylinder(0, 0, -H / 2, 0, 0, H, R)
    tools = []
    port_info = []
    wire_vols = []
    if loops:
        for k, ang in enumerate(PORT_ANG):
            sign = +1 if k == 0 else -1          # mirror images about the x-z plane
            u, _ = port_frame(ang, sign)
            a0 = scale(u, R - 0.5)
            rec = occ.addCylinder(*a0, *scale(u, L_COAX + 0.5), B_COAX)
            tools.append((3, rec))
            pts = loop_points(g, ang, sign, loops["w"], loops["h"], loops["orient"])
            segs = []
            for pa, pb in zip(pts[:-1], pts[1:]):
                dv = tuple(b - a for a, b in zip(pa, pb))
                segs.append((3, occ.addCylinder(*pa, *dv, R_WIRE)))
            for pc in pts[1:3]:
                segs.append((3, occ.addSphere(*pc, R_WIRE)))
            wv, _ = occ.fuse([segs[0]], segs[1:])
            wire_vols += wv
            port_info.append(dict(index=k + 1, angle_deg=ang, center=scale(u, R + L_COAX), u=u,
                                  points=pts))
        dom, _ = occ.fuse([(3, can)], tools)
    else:
        dom = [(3, can)]
    metal = cond + wire_vols
    if screw > 0:
        sc = occ.addCylinder(R - screw, 0, 0, screw + 3.0, 0, 0, SCREW_R)
        metal = metal + [(3, sc)]
    if axial:   # exploratory tuner: brass rod of radius axial[0] through the lid on the can axis,
        # penetrating axial[1] mm below the inner lid surface (into the loop bore)
        ar, ap = axial
        metal = metal + [(3, occ.addCylinder(0, 0, H / 2 - ap, 0, 0, ap + 3.0, ar))]
    dom, _ = occ.cut(dom, metal)

    # ---- dielectrics ----
    x_out = r + t + s
    slab = occ.addBox(r - 1, -d / 2, -z / 2, x_out - (r - 1), d, z)
    bore = occ.addCylinder(0, 0, -z / 2 - 1, 0, 0, z + 2, r)
    slab_c, _ = occ.cut([(3, slab)], [(3, bore)])
    diels = list(slab_c)
    if spacers:
        for z0 in (z / 2, -H / 2):
            so = occ.addCylinder(0, 0, z0, 0, 0, g["c"], r + t)
            si = occ.addCylinder(0, 0, z0 - 1, 0, 0, g["c"] + 2, r)
            ring, _ = occ.cut([(3, so)], [(3, si)])
            diels += ring
    if clips:
        # PEEK U-clips (hardware/lgr_a/lgr_a_assembly.py slab_clip): body = (depth + wall) x
        # (lip pack + 2 wall) x length, minus the slot that takes the lips and the slab
        cw, cdp, cln = CLIP["wall"], CLIP["depth"], CLIP["length"]
        pack = 2 * tl + d
        for z0 in (z / 2 - cln, -z / 2):
            body = occ.addBox(x_out - cdp, -pack / 2 - cw, z0, cdp + cw, pack + 2 * cw, cln)
            slot = occ.addBox(x_out - cdp - 1.0, -pack / 2, z0 - 1.0, cdp + 1.0, pack, cln + 2.0)
            clip, _ = occ.cut([(3, body)], [(3, slot)])
            diels += clip
    out, _ = occ.fragment(dom, diels)
    occ.synchronize()

    vols = [v for (dim, v) in gmsh.model.getEntities(3)]
    al, ptfe, peek, air = [], [], [], []
    for v in vols:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(3, v)
        if (ymax - ymin) < d * 1.01 and xmin > r - 1.5:
            al.append(v)
        elif clips and xmin > x_out - CLIP["depth"] - 1e-3 and xmax < x_out + CLIP["wall"] + 1e-3 \
                and zmax - zmin < CLIP["length"] + 1e-3:
            peek.append(v)
        elif zmin > z / 2 - 1e-3 or zmax < -z / 2 + 1e-3:
            ptfe.append(v)
        else:
            air.append(v)
    assert len(al) == 1 and len(air) in (1, 2) and len(ptfe) == (2 if spacers else 0), (al, ptfe, air)
    assert len(peek) == (2 if clips else 0), peek

    def bnd(vs):
        return {abs(tg) for (_, tg) in gmsh.model.getBoundary([(3, v) for v in vs], combined=False,
                                                              oriented=False)}
    b_all = [bnd([v]) for v in vols]
    from collections import Counter
    cnt = Counter(sf for b in b_all for sf in b)
    exterior = {sf for sf, n in cnt.items() if n == 1}
    shield, lgr, wires, brass, ports = [], [], [], [], {1: [], 2: []}
    tol = 1e-3 * R
    for sf in exterior:
        bb = gmsh.model.getBoundingBox(2, sf)
        xmin, ymin, zmin, xmax, ymax, zmax = bb
        cm = gmsh.model.occ.getCenterOfMass(2, sf)
        typ = gmsh.model.getType(2, sf)
        # ports: plane annulus at the recess end
        hit = None
        for p in port_info:
            if typ == "Plane" and math.dist(cm, p["center"]) < 0.05 and max(xmax - xmin, ymax - ymin, zmax - zmin) < 2 * B_COAX + 0.1:
                hit = p["index"]
        if hit:
            ports[hit].append(sf)
            continue
        on_can = (max(abs(xmin), abs(xmax), abs(ymin), abs(ymax)) > R - tol and zmax - zmin > H - tol) \
            or (abs(zmin - zmax) < tol and abs(abs(zmin) - H / 2) < tol)
        if on_can:
            shield.append(sf)
            continue
        if screw > 0 and xmin > R - screw - 1e-3 and max(abs(ymin), abs(ymax), abs(zmin), abs(zmax)) < SCREW_R + 1e-3:
            brass.append(sf)
            continue
        if axial and max(abs(xmin), abs(xmax), abs(ymin), abs(ymax)) < axial[0] + 1e-3 and zmin > H / 2 - axial[1] - 1e-3:
            brass.append(sf)
            continue
        near_port = any(math.dist(cm, p["points"][1]) < max(loops["w"], loops["h"]) + 4.0 or
                        math.dist(cm, p["center"]) < 4.0 for p in port_info) if loops else False
        if near_port and math.hypot(cm[0], cm[1]) > r + t + 0.5:
            wires.append(sf)
        else:
            lgr.append(sf)
    if loops:
        assert len(ports[1]) == 1 and len(ports[2]) == 1, ports
    gmsh.model.addPhysicalGroup(3, air, 1, name="air")
    gmsh.model.addPhysicalGroup(3, al, 2, name="alumina")
    if ptfe:
        gmsh.model.addPhysicalGroup(3, ptfe, 3, name="ptfe")
    if peek:
        gmsh.model.addPhysicalGroup(3, peek, 4, name="peek")
    gmsh.model.addPhysicalGroup(2, shield, 1, name="can")
    gmsh.model.addPhysicalGroup(2, lgr, 2, name="loop_lips")
    if wires:
        gmsh.model.addPhysicalGroup(2, wires, 3, name="coupling_wires")
    if brass:
        gmsh.model.addPhysicalGroup(2, brass, 4, name="screw")
    for k in (1, 2):
        if ports[k]:
            gmsh.model.addPhysicalGroup(2, ports[k], 10 + k, name=f"port{k}")

    # ---- mesh size fields (sweep recipe v2 + wire/port/screw refinement) ----
    h_gap = 0.45 * h * max(1.0, d / 0.635) ** 0.5
    h_edge = 0.25 * h
    h_loop = 2.0 * h
    h_bulk = max(R / 5.0, 4.0) * h
    h_medge = 0.2 * h
    h_wire = 0.3 * h * hl
    h_port = 0.15 * h * hl
    al_surfs = list(bnd(al))
    f = gmsh.model.mesh.field
    fl = []

    def thr(kind, ents, smin, dmin, dmax, sampling=150):
        a = f.add("Distance")
        f.setNumbers(a, kind, ents)
        f.setNumber(a, "Sampling", sampling)
        b = f.add("Threshold")
        f.setNumber(b, "InField", a)
        f.setNumber(b, "SizeMin", smin)
        f.setNumber(b, "SizeMax", h_bulk)
        f.setNumber(b, "DistMin", dmin)
        f.setNumber(b, "DistMax", dmax)
        fl.append(b)

    def curves_of(surfs):
        return sorted({abs(c) for (_, c) in gmsh.model.getBoundary([(2, x) for x in surfs], combined=False,
                                                                    oriented=False)})
    thr("SurfacesList", al_surfs, h_gap, 0.5, 12.0)
    thr("CurvesList", curves_of(al_surfs), h_edge, 0.2, 4.0, 200)
    thr("SurfacesList", lgr, h_loop, 1.0, 15.0)
    thr("CurvesList", curves_of(lgr), h_medge, 0.3, 6.0, 400)
    if wires:
        thr("SurfacesList", wires, h_wire, 0.4, 10.0, 200)
        thr("SurfacesList", ports[1] + ports[2], h_port, 0.2, 3.0, 100)
    if brass:
        thr("SurfacesList", brass, 0.5 * h, 0.5, 8.0, 150)
        thr("CurvesList", curves_of(brass), h_medge * 1.5, 0.3, 4.0, 200)
    fmin = f.add("Min")
    f.setNumbers(fmin, "FieldsList", fl)
    f.setAsBackgroundMesh(fmin)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 16)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_bulk)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.model.mesh.generate(3)
    etypes, etags, _ = gmsh.model.mesh.getElements(3)
    n_tet = sum(len(x) for x in etags)
    info = dict(n_tet=int(n_tet), h=h, hl=hl, h_gap=h_gap, h_edge=h_edge, h_loop=h_loop, h_bulk=h_bulk,
                h_medge=h_medge, h_wire=h_wire, h_port=h_port,
                n_surf=dict(can=len(shield), loop_lips=len(lgr), wires=len(wires), screw=len(brass)))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    gmsh.write(str(path))
    if gui:
        gmsh.fltk.run()
    gmsh.finalize()
    return info, port_info


def make_config(name, mesh_name, order, problem="Eigenmode", loops=True, spacers=True, screw=False,
                f_lo=None, f_hi=None, f_step=1e-6, adaptive_tol=1e-4, target=0.12, n_modes=2,
                linear=None, points=None, zs="exact", f_rs=F0, clips=False, empty=False, save=0):
    """zs: copper impedance input, "exact" (correct, default) or "s2" (the first runs, half the
    copper loss; see the module docstring). f_rs: frequency (Hz) at which Rs is evaluated, i.e. the
    expected f0 of the case. clips: PEEK slab clips present (volume attribute 4). empty: empty can.
    save: number of eigenmodes whose fields are written (ParaView output); 0 = none."""
    cfg = {
        "Problem": {"Type": problem, "Verbose": 2, "Output": f"postpro/{name}",
                    "OutputFormats": {"Paraview": save > 0, "GridFunction": False}},
        "Model": {"Mesh": f"mesh/{mesh_name}.msh", "L0": 1.0e-3},
        "Domains": {
            "Materials": [
                {"Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0},
                {"Attributes": [2], "Permeability": 1.0, "Permittivity": EPS_AL, "LossTan": TAND_AL}],
            "Postprocessing": {"Energy": [{"Index": 1, "Attributes": [2]}]}},
        "Boundaries": {"Impedance": [zs_entry([1] if empty else [1, 2] + ([3] if loops else []),
                                              SIGMA_CU, f_rs, zs)]},
        "Solver": {"Order": order, "Device": "CPU",
                   "Linear": {"Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-8, "MaxIts": 400}},
    }
    if spacers:
        cfg["Domains"]["Materials"].append(
            {"Attributes": [3], "Permeability": 1.0, "Permittivity": EPS_PTFE, "LossTan": TAND_PTFE})
        cfg["Domains"]["Postprocessing"]["Energy"].append({"Index": 2, "Attributes": [3]})
    if clips:
        cfg["Domains"]["Materials"].append(
            {"Attributes": [4], "Permeability": 1.0, "Permittivity": EPS_PEEK, "LossTan": TAND_PEEK})
        cfg["Domains"]["Postprocessing"]["Energy"].append({"Index": 3, "Attributes": [4]})
    if empty:
        del cfg["Domains"]["Materials"][1:]
        cfg["Domains"]["Postprocessing"]["Energy"] = [{"Index": 1, "Attributes": [1]}]
    if screw:
        cfg["Boundaries"]["Impedance"].append(zs_entry([4], SIGMA_BRASS, f_rs, zs))
    if loops:
        cfg["Boundaries"]["LumpedPort"] = [
            {"Index": k, "Attributes": [10 + k], "Direction": "+R", "R": Z_PORT} for k in (1, 2)]
    if problem == "Eigenmode":
        cfg["Solver"]["Eigenmode"] = {"N": n_modes, "Tol": 1.0e-8, "Target": target, "Save": save}
        cfg["Solver"]["Linear"]["MaxIts"] = 200
    else:
        for k, p in enumerate(cfg["Boundaries"]["LumpedPort"]):
            p["Excitation"] = k + 1          # two excitations -> full 2x2 S-matrix
        del cfg["Domains"]["Postprocessing"]
        if points:   # full-order solves at explicit frequencies (check of the reduced-order model)
            cfg["Solver"]["Driven"] = {"Samples": [{"Type": "Point", "Freq": points}], "AdaptiveTol": 0.0}
        else:
            cfg["Solver"]["Driven"] = {
                "Samples": [{"Type": "Linear", "MinFreq": f_lo, "MaxFreq": f_hi, "FreqStep": f_step}],
                "AdaptiveTol": adaptive_tol, "AdaptiveMaxSamples": 40}
    if linear:
        cfg["Solver"]["Linear"].update(linear)
    return cfg
