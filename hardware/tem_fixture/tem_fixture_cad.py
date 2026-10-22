"""Parametric CAD of the S2 measurement fixture (air microstrip TEM line) with a meta-atom slab.

Geometry (mm), from analysis/tem_fixture_design.py:
  - Ground plate: aluminium 3 mm, length L_tot + 2*margin, width W_g.
  - Top strip: aluminium 2 mm. Uniform section at height h = 70 mm (= S2 lattice a_z), width
    w = 344 mm (50 ohm), length L_u = 420 mm (room for up to 6 cells along x).
    Each end: a constant-impedance taper (w/h kept at 4.91) sloping down to h_c = 3 mm and
    w_c = 14.7 mm over L_t = 300 mm, ending on an N-type panel jack centre pin.
  - Supports: PTFE posts under the strip edges, outside the cell region.
  - Meta-atom slab: N_y x N_x LGR loops (from S2: r = 21.5, t = 3, length 45 mm along y, lips
    s = 6.2 mm towards +z), period 70 mm, centred under the strip; loop axes along y (across).
Visual/manufacturing draft only; full-wave check of the fixture (edge effects with 4 cells
across, taper match) is a separate Palace task.
"""
import math
from pathlib import Path

import cadquery as cq

OUT = Path(__file__).parent / "out"

P = dict(h=70.0, w=344.0, L_u=420.0, L_t=300.0, h_c=3.0, w_c=14.7, t_g=3.0, t_s=2.0,
         margin=50.0, W_g=504.0, a=70.0, r=21.5, t_loop=3.0, len_loop=45.0, s=6.2, t_lip=3.0,
         slit=1.02, d_slab=1.0, n_x=2, n_y=4)


def ground(p):
    L = p["L_u"] + 2 * p["L_t"] + 2 * p["margin"]
    return cq.Workplane("XY").box(L, p["W_g"], p["t_g"]).translate((0, 0, -p["t_g"] / 2))


def strip(p):
    hu, w = p["h"], p["w"]
    mid = cq.Workplane("XY").box(p["L_u"], w, p["t_s"]).translate((0, 0, hu + p["t_s"] / 2))
    parts = [mid]
    for sgn in (+1, -1):
        x0 = sgn * p["L_u"] / 2
        x1 = sgn * (p["L_u"] / 2 + p["L_t"])
        # tapered plate: trapezoid in plan, sloping from hu to h_c (lofted between two thin rects)
        r0 = cq.Workplane("YZ", origin=(x0, 0, hu)).rect(w, p["t_s"])
        prof = (cq.Workplane("YZ", origin=(x0, 0, hu + p["t_s"] / 2)).rect(w, p["t_s"])
                .workplane(offset=sgn * p["L_t"]).center(0, (p["h_c"] - hu)).rect(p["w_c"], p["t_s"]).loft())
        parts.append(prof)
    out = parts[0]
    for q in parts[1:]:
        out = out.union(q)
    return out


def n_jacks(p):
    jacks = []
    for sgn in (+1, -1):
        x = sgn * (p["L_u"] / 2 + p["L_t"] + 8)
        body = (cq.Workplane("YZ", origin=(x, 0, p["h_c"] + 1)).circle(8.0).extrude(sgn * 18))
        flange = cq.Workplane("XY").box(4, 25.4, 25.4).translate((x, 0, 12.7 - p["t_g"]))
        jacks.append(body.union(flange))
    return jacks


def posts(p):
    out = None
    ys = (p["w"] / 2 - 10, -(p["w"] / 2 - 10))
    xs = (-p["L_u"] / 2 + 15, p["L_u"] / 2 - 15)
    for x in xs:
        for y in ys:
            c = cq.Workplane("XY").circle(6).extrude(p["h"]).translate((x, y, 0))
            out = c if out is None else out.union(c)
    return out


def cell(p):
    """LGR loop with lips (axis along y), lips pointing +z, plus its alumina slab."""
    r, t, L, s, tl, w = p["r"], p["t_loop"], p["len_loop"], p["s"], p["t_lip"], p["slit"]
    tube = cq.Workplane("XZ").circle(r + t).circle(r).extrude(L / 2, both=True)
    lip_len = t + s
    for sgn in (+1, -1):
        lip = cq.Workplane("XY").box(tl, L, lip_len).translate((sgn * (w / 2 + tl / 2), 0, r + lip_len / 2))
        tube = tube.union(lip)
    slot = cq.Workplane("XY").box(w, L + 2, lip_len + 2).translate((0, 0, r + lip_len / 2 - 0.5))
    loop = tube.cut(slot)
    slab = cq.Workplane("XY").box(p["d_slab"], L, lip_len).translate((0, 0, r + lip_len / 2))
    return loop, slab


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p = P
    asm = cq.Assembly(name="TEM fixture")
    asm.add(ground(p), name="Ground plate (Al 3 mm)", color=cq.Color(0.72, 0.74, 0.77))
    asm.add(strip(p), name="Top strip with tapers (Al 2 mm)", color=cq.Color(0.80, 0.82, 0.85))
    for i, j in enumerate(n_jacks(p)):
        asm.add(j, name=f"N-type panel jack {'AB'[i]}", color=cq.Color(0.9, 0.78, 0.3))
    asm.add(posts(p), name="PTFE support posts", color=cq.Color("white"))
    loop, slab = cell(p)
    zc = p["h"] / 2
    for ix in range(p["n_x"]):
        for iy in range(p["n_y"]):
            x = (ix - (p["n_x"] - 1) / 2) * p["a"]
            y = (iy - (p["n_y"] - 1) / 2) * p["a"]
            zc_loop = zc - (p["s"] + p["t_loop"]) / 2  # centre the loop+lips envelope vertically
            asm.add(loop, name=f"Meta-atom loop x{ix}y{iy}", loc=cq.Location((x, y, zc_loop)),
                    color=cq.Color(0.78, 0.45, 0.22))
            asm.add(slab, name=f"Alumina slab x{ix}y{iy}", loc=cq.Location((x, y, zc_loop)),
                    color=cq.Color(0.96, 0.95, 0.88))
    asm.export(str(OUT / "tem_fixture.glb"))
    asm.export(str(OUT / "tem_fixture.step"))
    print("wrote", OUT / "tem_fixture.glb")


if __name__ == "__main__":
    main()
