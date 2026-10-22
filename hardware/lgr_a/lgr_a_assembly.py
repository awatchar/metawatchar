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

# Chosen first prototype (note 17): 0.5 L, a = 0.5, 1.0 mm alumina, s = 10.42 mm
P = params(volume_l=0.5, a=0.5, s=10.42, d_slab=1.0)

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
    pin = cq.Solid.makeCylinder(0.64, 6.0, base + d * (p["R"] - 6.5), d)
    insul = cq.Solid.makeCylinder(2.05, 1.5, base + d * (p["R"] - 1.0), d)
    parts["barrel"] = cq.Workplane().add(barrel).union(cq.Workplane().add(body))
    parts["nut"] = nut
    parts["pin"] = cq.Workplane().add(pin)
    parts["insul"] = cq.Workplane().add(insul)
    # coupling loop: 1 mm wire, ~8 mm loop in the plane containing the radial and z axes
    tip = base + d * (p["R"] - 6.5)
    path = (cq.Workplane("XY").moveTo(0, 0).lineTo(0, 8).lineTo(6.0, 8).lineTo(6.0, 0)).val()
    # build the loop in a local frame then place it: local x = -radial, local y = +z
    loop_pts = [tip, tip + cq.Vector(0, 0, 8), tip + cq.Vector(0, 0, 8) + d * 6.0, tip + d * 6.0]
    wire = cq.Workplane().polyline([(q.x, q.y, q.z) for q in loop_pts])
    segs = []
    for a_, b_ in zip(loop_pts[:-1], loop_pts[1:]):
        L = (b_ - a_).Length
        segs.append(cq.Solid.makeCylinder(0.5, L, a_, (b_ - a_).normalized()))
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
    ts = tuning_screw(P)
    asm.add(ts["screw"], name="Brass M4 tuning screw", color=COLORS["brass"])
    asm.add(ts["locknut"], name="Brass M4 lock nut", color=COLORS["brass"])
    asm.add(bolts(P), name="M3 bolts x8", color=COLORS["steel"])

    asm.export(str(OUT / "LGR-A_full.glb"))
    asm.export(str(OUT / "LGR-A_full.step"))

    # explode vectors (mm) for the viewer: lid/top parts up, ports and screw radially out
    zc = zb + P["H"] / 2
    ex = {
        "Copper lid (T2)": [0, 0, 70], "M3 bolts x8": [0, 0, 95], "PTFE spacer (top)": [0, 0, 50],
        "Copper loop with lips (T2)": [0, 0, 25], "Alumina slab 99.6% (1.0 mm)": [40, 0, 25],
        "PTFE spacer (bottom)": [0, 0, 0], "Copper can (T2)": [0, 0, 0],
        "Brass M4 tuning screw": [45, 0, 0], "Brass M4 lock nut": [30, 0, 0],
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
