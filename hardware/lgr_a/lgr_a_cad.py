"""Parametric CAD of the LGR-A loop-gap-resonator unit-cell prototype (CadQuery).

LGR-A geometry (shared with the Palace size sweep, sim/lgr_size_sweep):
  - Shield: closed copper can, inner radius R, inner height H = 2R (inner volume 2*pi*R^3).
  - Loop: coaxial copper tube, inner radius r = a*R, wall t, length z = H - 2c, c = 0.25*R,
    held on axis by two PTFE spacer rings (height c) at the bottom and under the lid.
  - Slit: one axial slit through the loop wall, width d (0.635 mm alumina substrate + clearance).
  - Lips: two radial copper lips (thickness t_lip, radial height s) at the slit; the alumina slab
    (thickness d, radial extent t + s, length z) sits between them and forms the gap capacitor.
  - Two SMA bulkhead ports at 120/240 deg, mid-height, each with a 4.7 x 4.7 mm wire coupling loop
    (radial depth from the inner wall x tangential width, wire centreline) lying in the MID-HEIGHT
    PLANE (normal along the can axis, links the axial return flux). The return leg runs parallel to
    the radial leg, so the ground end is soldered at a 1.2 mm hole asin(4.7/R) towards 180 deg
    (6.27 deg at 0.5 L). Palace with corrected copper: |S21(f0)| = -29.9 dB.
  - No tuning screw (Palace: < 5 kHz range). f0 is set by the lip height s; order s variants.

Outputs STEP files per part and an assembly, plus a parameter JSON. Units: mm.
Run with the project venv:  D:\\PhD\\.venv\\Scripts\\python.exe lgr_a_cad.py [--volume-l 0.25] [--a 0.5] [--s 2.5]
"""
import argparse
import json
import math
from pathlib import Path

import cadquery as cq


def params(volume_l=0.25, a=0.5, s=2.5, t=3.0, t_lip=3.0, d_slab=0.635, slit_clear=0.025,
           wall=4.0, base=4.0, lid=4.0, flange=7.0, n_bolts=8, loop_w=4.7, loop_h=4.7):
    R = (volume_l * 1e6 / (2 * math.pi)) ** (1 / 3)  # mm, inner radius (V = 2 pi R^3)
    H = 2 * R
    r = a * R
    c = 0.25 * R
    z = H - 2 * c
    p = dict(volume_l=volume_l, a=a, R=R, H=H, r=r, t=t, c=c, z=z, s=s, t_lip=t_lip,
             d_slab=d_slab, slit=d_slab + slit_clear, wall=wall, base=base, lid=lid,
             flange=flange, n_bolts=n_bolts, bolt_clear=3.4, bolt_tap="M3",
             sma_hole=6.5, loop_w=loop_w, loop_h=loop_h,
             loop_gnd_deg=math.degrees(math.asin(loop_h / R)), loop_gnd_hole=1.2)
    assert r + t + s <= R - 3, "lips do not fit inside the can (r + t + s > R - 3 mm)"
    return p


def can_body(p):
    Ro = p["R"] + p["wall"]
    Rf = Ro + p["flange"]
    body = (cq.Workplane("XY").circle(Ro).extrude(p["base"] + p["H"])
            .faces(">Z").workplane().circle(p["R"]).cutBlind(-p["H"]))
    flange = (cq.Workplane("XY").workplane(offset=p["base"] + p["H"] - 6)
              .circle(Rf).circle(Ro - 0.01).extrude(6))
    body = body.union(flange)
    # tapped bolt holes in the flange (modelled as tap-drill holes; thread called out on drawing)
    bc = (Ro + Rf) / 2
    pts = [(bc * math.cos(2 * math.pi * k / p["n_bolts"]), bc * math.sin(2 * math.pi * k / p["n_bolts"]))
           for k in range(p["n_bolts"])]
    body = body.faces(">Z").workplane().pushPoints(pts).hole(2.5, 6)
    zc = p["base"] + p["H"] / 2
    # SMA ports at 120 and 240 deg; small holes for the coupling-loop ground ends asin(loop_h/R) towards 180 deg
    # (the M4 tuning screw was dropped: Palace showed < 5 kHz of tuning, sim/lgr_a_driven)
    for ang, dia in ((120, p["sma_hole"]), (240, p["sma_hole"]),
                     (120 + p["loop_gnd_deg"], p["loop_gnd_hole"]), (240 - p["loop_gnd_deg"], p["loop_gnd_hole"])):
        hole = (cq.Workplane("XY").workplane(offset=zc).transformed(rotate=(0, 90, ang))
                .circle(dia / 2).extrude(Ro + 5))
        hole = cq.Workplane("XY").add(
            cq.Solid.makeCylinder(dia / 2, Ro + 5, cq.Vector(0, 0, zc),
                                  cq.Vector(math.cos(math.radians(ang)), math.sin(math.radians(ang)), 0)))
        body = body.cut(hole)
    return body


def lid(p):
    Rf = p["R"] + p["wall"] + p["flange"]
    Ro = p["R"] + p["wall"]
    bc = (Ro + Rf) / 2
    pts = [(bc * math.cos(2 * math.pi * k / p["n_bolts"]), bc * math.sin(2 * math.pi * k / p["n_bolts"]))
           for k in range(p["n_bolts"])]
    return (cq.Workplane("XY").circle(Rf).extrude(p["lid"])
            .faces(">Z").workplane().pushPoints(pts).hole(p["bolt_clear"]))


def loop(p):
    r, t, z, s, tl, w = p["r"], p["t"], p["z"], p["s"], p["t_lip"], p["slit"]
    tube = cq.Workplane("XY").circle(r + t).circle(r).extrude(z)
    # lips: two plates on either side of the slit (slit centred on the +X axis)
    lip_len = t + s
    for sign in (+1, -1):
        y0 = sign * (w / 2 + tl / 2)
        lip = (cq.Workplane("XY").center(r + lip_len / 2, y0).rect(lip_len, tl).extrude(z))
        tube = tube.union(lip)
    slit = cq.Workplane("XY").center(r + lip_len / 2 - 0.5, 0).rect(lip_len + 2, w).extrude(z)
    return tube.cut(slit)


def spacer(p):
    return cq.Workplane("XY").circle(p["r"] + p["t"]).circle(p["r"]).extrude(p["c"])


def alumina(p):
    lip_len = p["t"] + p["s"]
    return cq.Workplane("XY").center(p["r"] + lip_len / 2, 0).rect(lip_len, p["d_slab"]).extrude(p["z"])


def build(p, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    parts = {
        "can_body_Cu-T2": can_body(p),
        "lid_Cu-T2": lid(p),
        "loop_Cu-T2": loop(p),
        "spacer_PTFE": spacer(p),
        "slab_alumina-99.6": alumina(p),
    }
    for name, wp in parts.items():
        cq.exporters.export(wp, str(out / f"{name}.step"))
    zb = p["base"]
    asm = (cq.Assembly(name="LGR-A")
           .add(parts["can_body_Cu-T2"], name="can", color=cq.Color(0.80, 0.50, 0.20, 0.35))
           .add(parts["lid_Cu-T2"], name="lid", loc=cq.Location((0, 0, zb + p["H"])),
                color=cq.Color(0.80, 0.50, 0.20, 0.35))
           .add(parts["spacer_PTFE"], name="spacer_bottom", loc=cq.Location((0, 0, zb)), color=cq.Color("white"))
           .add(parts["loop_Cu-T2"], name="loop", loc=cq.Location((0, 0, zb + p["c"])),
                color=cq.Color(0.85, 0.55, 0.25))
           .add(parts["slab_alumina-99.6"], name="slab", loc=cq.Location((0, 0, zb + p["c"])),
                color=cq.Color(0.95, 0.95, 0.90))
           .add(parts["spacer_PTFE"], name="spacer_top", loc=cq.Location((0, 0, zb + p["c"] + p["z"])),
                color=cq.Color("white")))
    asm.save(str(out / "LGR-A_assembly.step"))
    (out / "params.json").write_text(json.dumps(p, indent=2))
    return parts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--volume-l", type=float, default=0.25)
    ap.add_argument("--a", type=float, default=0.5)
    ap.add_argument("--s", type=float, default=10.55, help="lip radial height (mm); 10.55 restores ~145.2 MHz with R0.4 fillets (sim/lgr_a_driven)")
    ap.add_argument("--d", type=float, default=0.635, help="alumina slab thickness (mm): 0.635 or 1.0")
    ap.add_argument("--out", default=str(Path(__file__).parent / "out"))
    args = ap.parse_args()
    p = params(args.volume_l, args.a, args.s, d_slab=args.d)
    build(p, Path(args.out))
    print(json.dumps({k: round(v, 3) if isinstance(v, float) else v for k, v in p.items()}, indent=1))


if __name__ == "__main__":
    main()
