"""Full LGR-A assembly for visualization: machined parts + bought parts, exported as GLB and STEP.

Machined parts come from lgr_a_cad.py (same parameters as the Palace sweep). The loop lips get
edge fillets (FILLET mm) for manufacturing and to soften edge fields; Palace modelled sharp edges,
so f0 may rise slightly (re-tune with the screw). Bought parts (SMA bulkhead jacks with coupling
loops, brass M4 tuning screw + lock nut, M3 bolts) are simplified look-alikes for visualization,
NOT manufacturing models.

Run:  D:\\PhD\\.venv\\Scripts\\python.exe lgr_a_assembly.py
Writes out/LGR-A_full.glb, out/LGR-A_full.step, out/explode.json (per-part explode vectors).
"""
import json
import math
from pathlib import Path

import cadquery as cq

from lgr_a_cad import params, can_body, lid, loop, spacer, alumina, build

OUT = Path(__file__).parent / "out"
FILLET = 0.4  # mm, lip tip edges

# Chosen first prototype (notes 17, 21): 0.5 L, a = 0.5, 1.0 mm alumina; s = 10.55 mm restores ~145.2 MHz
# with the R0.4 lip fillets. Loop variants for frequency trim are exported separately (S_VARIANTS).
P = params(volume_l=0.5, a=0.5, s=10.55, d_slab=1.0)
S_VARIANTS = (10.42, 10.55, 10.70)
CLIP_OPENINGS = (6.95, 7.00, 7.05)  # mm, PEEK clip inner opening

COLORS = {
    "copper": cq.Color(0.78, 0.45, 0.22),
    "copper_dark": cq.Color(0.66, 0.37, 0.18),
    "brass": cq.Color(0.83, 0.69, 0.32),
    "alumina": cq.Color(0.96, 0.95, 0.88),
    "ptfe": cq.Color(0.93, 0.93, 0.95),
    "steel": cq.Color(0.62, 0.64, 0.67),
    "gold": cq.Color(0.95, 0.78, 0.25),
    "sma_dielectric": cq.Color(0.95, 0.95, 0.95),
}


def loop_filleted(p):
    lp = loop(p)
    try:  # round the vertical edges at the lip tips (outermost in +X)
        lp = lp.edges("|Z").edges(cq.selectors.BoxSelector(
            (p["r"] + p["t"] + p["s"] - 0.5, -20, -1), (p["r"] + p["t"] + p["s"] + 0.5, 20, p["z"] + 1))).fillet(FILLET)
    except Exception as exc:  # keep sharp edges rather than fail the build
        print(f"lip fillet skipped: {exc}")
    return lp


def radial(ang_deg):
    a = math.radians(ang_deg)
    return math.cos(a), math.sin(a)


def sma_jack(p, ang):
    """Single-hole SMA bulkhead jack: outside body + hex nut, threaded barrel through the wall,
    PTFE insulator and centre pin inside, with a small wire coupling loop to the can wall."""
    ux, uy = radial(ang)
    Ro = p["R"] + p["wall"]
    zc = p["base"] + p["H"] / 2
    d = cq.Vector(ux, uy, 0)
    parts = {}
    base = cq.Vector(0, 0, zc)
    # barrel through the wall (from R - 1 to Ro + 9)
    barrel = cq.Solid.makeCylinder(3.15, Ro + 9 - (p["R"] - 1), base + d * (p["R"] - 1), d)
    nut = cq.Workplane().add(cq.Solid.makeCylinder(5.0, 2.5, base + d * (Ro + 0.2), d))
    body = cq.Solid.makeCylinder(4.0, 9.0, base + d * (Ro + 3.0), d)
    pin = cq.Solid.makeCylinder(0.64, 1.5, base + d * (p["R"] - 1.0), d)
    insul = cq.Solid.makeCylinder(2.05, 1.5, base + d * (p["R"] - 1.0), d)
    parts["barrel"] = cq.Workplane().add(barrel).union(cq.Workplane().add(body))
    parts["nut"] = nut
    parts["pin"] = cq.Workplane().add(pin)
    parts["insul"] = cq.Workplane().add(insul)
    # coupling loop: 1 mm Cu wire in the MID-HEIGHT plane, loop_w (radial depth from the inner
    # wall) x loop_h (tangential), wire centreline, as the Palace model (sim/lgr_a_driven):
    # leg 1 radially inward from the pin, leg 2 tangential towards 180 deg, leg 3 parallel to
    # leg 1 back out to the wall (ground hole at asin(loop_h/R))
    sgn = 1 if ang < 180 else -1                      # mirror the two loops about 180 deg
    t = cq.Vector(-uy, ux, 0) * sgn                    # tangential unit vector towards 180 deg
    p0 = base + d * (p["R"] - 0.5)
    p1 = base + d * (p["R"] - p["loop_w"])
    p2 = p1 + t * p["loop_h"]
    p3 = base + t * p["loop_h"] + d * math.sqrt(p["R"] ** 2 - p["loop_h"] ** 2)
    loop_pts = [p0, p1, p2, p3]
    segs = []
    for a_, b_ in zip(loop_pts[:-1], loop_pts[1:]):
        segs.append(cq.Solid.makeCylinder(0.5, (b_ - a_).Length, a_, (b_ - a_).normalized()))
    lw = cq.Workplane().add(segs[0])
    for s_ in segs[1:]:
        lw = lw.union(cq.Workplane().add(s_))
    parts["loop_wire"] = lw
    return parts


def tuning_screw(p, ang=0.0, penetration=4.0):
    ux, uy = radial(ang)
    Ro = p["R"] + p["wall"]
    zc = p["base"] + p["H"] / 2
    d = cq.Vector(ux, uy, 0)
    base = cq.Vector(0, 0, zc)
    shank_start = p["R"] - penetration
    shank = cq.Solid.makeCylinder(2.0, Ro + 6 - shank_start, base + d * shank_start, d)
    head = cq.Solid.makeCylinder(3.5, 2.8, base + d * (Ro + 6), d)
    nut = cq.Solid.makeCylinder(3.9, 3.2, base + d * (Ro + 0.3), d)
    screw = cq.Workplane().add(shank).union(cq.Workplane().add(head))
    return {"screw": screw, "locknut": cq.Workplane().add(nut)}


CLIP = dict(wall=2.0, depth=6.0, length=4.0, interference=0.05)  # mm, PEEK


def slab_clip(p, at_top, interference=None):
    """PEEK U-clip straddling both lips at one axial end of the loop, squeezing the lips onto the
    alumina slab so no air film forms (sim/lgr_a_driven: an air film adds ~0.7 MHz/um and drifts with T).
    Inner opening = 2*t_lip + d_slab - interference; it covers the outer `depth` mm of the lips."""
    t_lip, d, s_ = p["t_lip"], p["d_slab"], p["s"]
    opening = 2 * t_lip + d - (CLIP["interference"] if interference is None else interference)
    outer = opening + 2 * CLIP["wall"]
    x_tip = p["r"] + p["t"] + s_
    x0 = x_tip - CLIP["depth"]
    body = cq.Workplane("XY").box(CLIP["depth"] + CLIP["wall"], outer, CLIP["length"]).translate(
        (x0 + (CLIP["depth"] + CLIP["wall"]) / 2, 0, 0))
    slot = cq.Workplane("XY").box(CLIP["depth"] + 0.01, opening, CLIP["length"] + 2).translate(
        (x0 + CLIP["depth"] / 2 - 0.005, 0, 0))
    z = p["base"] + p["c"] + (p["z"] - CLIP["length"] / 2 if at_top else CLIP["length"] / 2)
    return body.cut(slot).translate((0, 0, z))


def bolts(p):
    Ro = p["R"] + p["wall"]
    Rf = Ro + p["flange"]
    bc = (Ro + Rf) / 2
    ztop = p["base"] + p["H"] + p["lid"]
    out = None
    for k in range(p["n_bolts"]):
        x, y = bc * math.cos(2 * math.pi * k / p["n_bolts"]), bc * math.sin(2 * math.pi * k / p["n_bolts"])
        b = (cq.Workplane("XY").workplane(offset=ztop - p["lid"] - 6).center(x, y).circle(1.5).extrude(p["lid"] + 6)
             .faces(">Z").workplane().circle(2.75).extrude(2.0))
        out = b if out is None else out.union(b)
    return out


def main():
    build(P, OUT)  # refresh per-part STEP files and params.json with the chosen dimensions
    zb = P["base"]
    asm = cq.Assembly(name="LGR-A")
    asm.add(can_body(P), name="Copper can (T2)", color=COLORS["copper"])
    asm.add(lid(P), name="Copper lid (T2)", loc=cq.Location((0, 0, zb + P["H"])), color=COLORS["copper_dark"])
    asm.add(spacer(P), name="PTFE spacer (bottom)", loc=cq.Location((0, 0, zb)), color=COLORS["ptfe"])
    asm.add(loop_filleted(P), name="Copper loop with lips (T2)", loc=cq.Location((0, 0, zb + P["c"])),
            color=COLORS["copper"])
    asm.add(alumina(P), name="Alumina slab 99.6% (1.0 mm)", loc=cq.Location((0, 0, zb + P["c"])),
            color=COLORS["alumina"])
    asm.add(spacer(P), name="PTFE spacer (top)", loc=cq.Location((0, 0, zb + P["c"] + P["z"])),
            color=COLORS["ptfe"])
    for ang, tag in ((120, "A"), (240, "B")):
        sp = sma_jack(P, ang)
        asm.add(sp["barrel"], name=f"SMA bulkhead jack {tag}", color=COLORS["gold"])
        asm.add(sp["nut"], name=f"SMA panel nut {tag}", color=COLORS["gold"])
        asm.add(sp["insul"], name=f"SMA PTFE insulator {tag}", color=COLORS["sma_dielectric"])
        asm.add(sp["pin"], name=f"SMA centre pin {tag}", color=COLORS["gold"])
        asm.add(sp["loop_wire"], name=f"Coupling loop {tag} (Cu wire)", color=COLORS["copper"])
    asm.add(bolts(P), name="M3 bolts x8", color=COLORS["steel"])
    peek = cq.Color(0.82, 0.70, 0.52)
    asm.add(slab_clip(P, at_top=False), name="PEEK slab clip (bottom)", color=peek)
    asm.add(slab_clip(P, at_top=True), name="PEEK slab clip (top)", color=peek)

    asm.export(str(OUT / "LGR-A_full.glb"))
    for sv in S_VARIANTS:  # loop variants: each 0.1 mm of lip height ~ 0.5 MHz (notes 17, 21)
        cq.exporters.export(loop_filleted(params(volume_l=0.5, a=0.5, s=sv, d_slab=1.0)),
                            str(OUT / f"loop_Cu-T2_s{str(sv).replace('.', 'p')}.step"))
    asm.export(str(OUT / "LGR-A_full.step"))
    for op in CLIP_OPENINGS:  # clip variants: 6.95 (0.05 mm interference), 7.00, 7.05 mm openings
        clip = slab_clip(P, at_top=False, interference=2 * P["t_lip"] + P["d_slab"] - op)
        cq.exporters.export(clip, str(OUT / f"clip_PEEK_open{str(op).replace('.', 'p')}.step"))

    # explode vectors (mm) for the viewer: lid/top parts up, ports and screw radially out
    zc = zb + P["H"] / 2
    ex = {
        "Copper lid (T2)": [0, 0, 70], "M3 bolts x8": [0, 0, 95], "PTFE spacer (top)": [0, 0, 50],
        "Copper loop with lips (T2)": [0, 0, 25], "Alumina slab 99.6% (1.0 mm)": [40, 0, 25],
        "PTFE spacer (bottom)": [0, 0, 0], "Copper can (T2)": [0, 0, 0],
        "PEEK slab clip (top)": [30, 0, 45], "PEEK slab clip (bottom)": [30, 0, -5],
    }
    for ang, tag in ((120, "A"), (240, "B")):
        ux, uy = radial(ang)
        for nm, k in (("SMA bulkhead jack", 40), ("SMA panel nut", 28), ("SMA PTFE insulator", 40),
                      ("SMA centre pin", 40), ("Coupling loop", 40)):
            label = f"{nm} {tag}" + (" (Cu wire)" if nm == "Coupling loop" else "")
            ex[label] = [round(ux * k, 2), round(uy * k, 2), 0]
    (OUT / "explode.json").write_text(json.dumps({"params": P, "explode_mm": ex, "fillet_mm": FILLET}, indent=1))
    print("wrote", OUT / "LGR-A_full.glb")
    print(json.dumps({k: round(v, 2) if isinstance(v, float) else v for k, v in P.items()}))


if __name__ == "__main__":
    main()
