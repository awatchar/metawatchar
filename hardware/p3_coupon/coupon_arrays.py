"""P3 board 2: 4 x 4 mini-arrays of each drum-cell type + printed-spiral reference cells.

  "C:/Program Files/KiCad/10.0/bin/python.exe" coupon_arrays.py

- Mini-arrays (one per cell type) show the collective (uniform) mode of a small lattice
  against the coupled-circuit model: type A at 15 mm pitch, B and C at 12.5 mm. The
  capacitor sits in the gap to the right of each inductor; no trim pads (no room).
- Printed reference cells: 38 mm square spiral, 17 turns, 0.6 mm trace / 0.25 mm gap,
  top layer only, crossover on the bottom layer; L ~ 7.7 uH (modified Wheeler), closed by
  2 x 220 nF C0G (86.7 kHz nominal) plus a DNP trim pad. Two closed cells plus one spiral
  with an edge-SMA footprint (no capacitor) for a direct Z measurement of L and R.
- Order with 2 oz copper (raises the printed-cell Q; upper bound ~6 single-layer).
"""
import os

import pcbnew

import p3pcb
from p3pcb import pt

OUT = os.path.join(p3pcb.HERE, "out")
NAME = "p3_array_coupon"

SP = dict(d_out=38.0, n=17, w=0.6, s=0.25)


def printed_cell(B, cid, x, y, caps=True, sma=False):
    """Spiral centred at (x, y); connections leave from its left side."""
    net = B.net(cid)
    (xo, yo), (xi, yi) = B.spiral(net, net, x, y, **SP)
    # crossover: via at the inner end, bottom-layer trace straight left, out of the spiral
    xl = x - SP["d_out"] / 2 - 2.0
    B.via(pt(xi, yi), net)
    B.track(pt(xi, yi), pt(xl, yi), net, width=0.8, layer=pcbnew.B_Cu)
    B.via(pt(xl, yi), net)
    # outer end: top-layer stub left
    B.track(pt(xo, yo), pt(xl, yo), net, width=0.8)
    if caps:
        refs = []
        for k, dx in enumerate((3.0, 7.0, 11.0)):
            dnp = k == 2
            xc = xl - dx
            fc = B.place_part(p3pcb.FP_CAP, p3pcb.CAP_220N["c_fp"], "C",
                              dict(val=p3pcb.CAP_220N["c_val"], mpn=p3pcb.CAP_220N["c_mpn"],
                                   lcsc=p3pcb.CAP_220N["c_lcsc"]),
                              xc, (yo + yi) / 2, 90, dnp=dnp, val=None if not dnp else "DNP trim")
            pc = sorted(fc.Pads(), key=lambda p: p.GetPosition().y)
            for p in pc:
                p.SetNet(net)
            refs.append(fc.GetReference())
            # top pad up to the outer-end rail, bottom pad down to the crossover rail
            B.track(pc[0].GetPosition(), pt(xc, yo), net, width=0.8)
            B.track(pc[1].GetPosition(), pt(xc, yi), net, width=0.8)
        xr = xl - 11.0
        B.track(pt(xl, yo), pt(xr, yo), net, width=0.8)
        B.track(pt(xl, yi), pt(xr, yi), net, width=0.8)
        B.cells.append(dict(id=cid, type="printed", x_mm=x, y_mm=y, has_c=True, group="printed", c=refs))
    if sma:
        # Edge SMA (hand-soldered, DNP) on the bottom edge, pointing out of the board.
        # Footprint: board edge at local x = +2.1, so rotate -90 and put the origin 2.1 mm inside.
        xs = xl - 6.0
        fs = B.add_fp(p3pcb.FP_COAX, "SMA_Samtec_SMA-J-P-H-ST-EM1_EdgeMount", B.ref("J"), "SMA edge (hand solder)",
                      xs, B.h - 2.1, -90, dnp=True)
        for p in fs.Pads():
            p.SetNet(net)
        sig = [p for p in fs.Pads() if p.GetNumber() == "1"][0].GetPosition()
        gnd_f = sorted([p for p in fs.Pads() if p.GetNumber() == "2" and p.IsOnLayer(pcbnew.F_Cu)],
                       key=lambda p: p.GetPosition().x)
        # outer end -> centre pin: left to x = xs, then straight down (stays left of the ground rail)
        B.track(pt(xl, yo), pt(xs, yo), net, width=0.6)
        B.track(pt(xs, yo), sig, net, width=0.6)
        # crossover -> right-hand ground pad: left to that pad's x, then down
        gx = pcbnew.ToMM(gnd_f[1].GetPosition().x)
        B.track(pt(xl, yi), pt(gx, yi), net, width=0.6)
        B.track(pt(gx, yi), gnd_f[1].GetPosition(), net, width=0.6)
        # tie both ground pads front/back with vias and join them on the bottom layer
        for g in gnd_f:
            B.via(g.GetPosition(), net)
        B.track(gnd_f[0].GetPosition(), gnd_f[1].GetPosition(), net, width=0.6, layer=pcbnew.B_Cu)
        B.cells.append(dict(id=cid, type="printed-Z", x_mm=x, y_mm=y, has_c=False, group="printed"))
    B.text(cid, x, y - SP["d_out"] / 2 - 3.0)


def build():
    B = p3pcb.Builder(NAME, 270.0, 185.0)
    B.outline()
    arrays = (("A", 20.0, 15.0, 7.5), ("B", 115.0, 12.5, 6.25), ("C", 200.0, 12.5, 6.25))
    for typ, x0, p, dx in arrays:
        for i in range(4):
            for j in range(4):
                B.cell(f"M{typ}{j + 1}{i + 1}", typ, x0 + p * i, 25.0 + p * j, trim=False, dx1=dx,
                       group=f"array{typ}", label="none")
        B.text(f"Array {typ}: 4 x 4, pitch {p:g} mm (rows 1-4 top-down, cols 1-4 left-right)",
               x0 + 1.5 * p, 25.0 + 3 * p + 10.0, 1.0)
    printed_cell(B, "S1", 62.0, 128.0)
    printed_cell(B, "S2", 152.0, 128.0)
    printed_cell(B, "SZ", 240.0, 135.0, caps=False, sma=True)
    B.text("P3 board 2 v1: mini-arrays + printed reference cells  (metawatchar, 2026)", 120, B.h - 9, 1.4)
    B.text("Order 2 oz copper. No copper pours. J1 (SMA edge) is hand-soldered.", 120, B.h - 5.5, 1.0)
    return B


if __name__ == "__main__":
    B = build()
    B.save(OUT)
    print(NAME, len(B.cells), "cells", len(B.parts), "parts")
    print(p3pcb.export(OUT, NAME))
