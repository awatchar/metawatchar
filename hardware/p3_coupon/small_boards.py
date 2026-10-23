"""P3 small hand-assembled boards: Z-tab, shielded loop probes, Helmholtz coil.

  "C:/Program Files/KiCad/10.0/bin/python.exe" small_boards.py

- p3_ztab (2-layer, 11 x 52 mm): one drum inductor (universal land fits types A, B, C)
  at the tip, an edge SMA at the other end, feed as a broadside pair (top/bottom) to keep
  the lead loop small. Two tabs side by side (down to 11.5 mm centre spacing) or stacked
  give Z21 = j w M vs distance, i.e. the inductor's external dipole moment m_eff. One tab
  alone gives L, ESR and Q of the inductor.
- p3_probe10 / p3_probe20 (4-layer): shielded loop probes, loop OD 10 / 20 mm,
  2 turns on each inner layer (4 in series). Outer layers are a ground shield with a slit
  from the loop centre to the tip, so the shield passes H but blocks E and forms no shorted
  turn. Loop return meets ground through a net tie at the connector.
- p3_helmholtz (2-layer, 110 x 110 mm, order 2 oz): one coil of a Helmholtz pair,
  10-turn circular spiral, mean radius 40 mm; two boards on 40 mm spacers. Known field for
  calibrating the probes (field computed from the actual track geometry).
SMA connectors and the Z-tab inductors are hand-soldered (not JLCPCB-assembled).
"""
import math
import os

import pcbnew

import p3pcb
from p3pcb import pt

OUT = os.path.join(p3pcb.HERE, "out")
FP_NETTIE = os.path.join(p3pcb.FP_ROOT, "NetTie.pretty")
FP_HDR = os.path.join(p3pcb.FP_ROOT, "Connector_PinHeader_2.54mm.pretty")


def edge_sma(B, x, net_sig, net_gnd):
    """Edge SMA on the bottom edge (pointing +y). Returns (signal pad pos, [ground front pads])."""
    fs = B.add_fp(p3pcb.FP_COAX, "SMA_Samtec_SMA-J-P-H-ST-EM1_EdgeMount", B.ref("J"), "SMA edge (hand solder)",
                  x, B.h - 2.1, -90, dnp=True)
    for p in fs.Pads():
        p.SetNet(net_sig if p.GetNumber() == "1" else net_gnd)
    sig = [p for p in fs.Pads() if p.GetNumber() == "1"][0].GetPosition()
    gf = sorted([p for p in fs.Pads() if p.GetNumber() == "2" and p.IsOnLayer(pcbnew.F_Cu)],
                key=lambda p: p.GetPosition().x)
    for g in gf:
        B.via(g.GetPosition(), net_gnd)
    return sig, gf


def ztab():
    B = p3pcb.Builder("p3_ztab", 11.0, 52.0)
    B.outline(holes=False, fiducials=False)
    cx = 5.5
    na, nb = B.net("Z_a"), B.net("Z_b")
    fl = B.add_fp(p3pcb.FP_LOCAL, "L_Drum_Universal_A_B_C", B.ref("L"), "drum L (hand solder)", cx, 7.0, 0, dnp=True)
    pads = sorted(fl.Pads(), key=lambda p: p.GetPosition().y)
    pads[0].SetNet(na)
    pads[1].SetNet(nb)
    sig, gf = edge_sma(B, cx, nb, na)
    # pad 2 -> centre pin on top; pad 1 -> via -> bottom, directly under the top trace
    B.track(pads[1].GetPosition(), sig, nb, width=0.6)
    B.via(pt(cx, 4.0), na)
    B.track(pt(cx, 4.0), pt(cx, 44.0), na, width=0.6, layer=pcbnew.B_Cu)
    g = gf[1].GetPosition()
    B.track(pt(cx, 44.0), pt(pcbnew.ToMM(g.x), 44.0), na, width=0.6, layer=pcbnew.B_Cu)
    B.track(pt(pcbnew.ToMM(g.x), 44.0), g, na, width=0.6, layer=pcbnew.B_Cu)
    B.track(gf[0].GetPosition(), gf[1].GetPosition(), na, width=0.6, layer=pcbnew.B_Cu)
    B.text("Z", cx, 22.0, 1.5)
    B.text("P3", cx, 27.0, 1.0)
    return B


def zone(B, layer, net, poly):
    z = pcbnew.ZONE(B.b)
    z.SetLayer(layer)
    z.SetNet(net)
    o = z.Outline()
    o.NewOutline()
    for x, y in poly:
        o.Append(pcbnew.FromMM(x), pcbnew.FromMM(y))
    z.SetLocalClearance(pcbnew.FromMM(0.3))
    z.SetMinThickness(pcbnew.FromMM(0.25))
    B.b.Add(z)
    return z


def loop_probe(od):
    w_board = od + 2.0
    B = p3pcb.Builder(f"p3_probe{od:g}", w_board, 85.0, layers=4)
    B.outline(holes=False, fiducials=False)
    cx, cy = w_board / 2, od / 2 + 1.5
    gnd, sig = B.net("GND"), B.net("PROBE_SIG")
    tw, pitch, turns = 0.25, 0.5, 2
    r_out = od / 2 - 0.5
    nseg = 48

    def arc(r0, r1, th0, nturn):
        pts = []
        for k in range(nturn * nseg + 1):
            f = k / (nturn * nseg)
            th = th0 + 2 * math.pi * nturn * f
            r = r0 + (r1 - r0) * f
            pts.append((cx + r * math.sin(th), cy + r * math.cos(th)))   # th = 0 points to +y (the stem)
        return pts

    th0 = 0.0
    p1 = arc(r_out, r_out - turns * pitch, th0, turns)                                  # In1, inwards
    p2 = arc(r_out - turns * pitch, r_out, th0 + 2 * math.pi * turns, turns)            # In2, outwards, same sense
    for (a, b) in zip(p1[:-1], p1[1:]):
        B.track(pt(*a), pt(*b), sig, width=tw, layer=pcbnew.In1_Cu)
    for (a, b) in zip(p2[:-1], p2[1:]):
        B.track(pt(*a), pt(*b), sig, width=tw, layer=pcbnew.In2_Cu)
    # The loop is one conductor (net PROBE_SIG): In1 inwards, a via at the inner end, In2 outwards.
    # It meets ground only at the connector, through a net tie (the one place SIG touches GND).
    B.via(pt(*p1[-1]), sig, size=0.5)
    y_sig_via, y_ret_via = B.h - 6.0, B.h - 12.0
    xs = cx
    B.track(pt(*p1[0]), pt(xs, y_sig_via), sig, width=tw, layer=pcbnew.In1_Cu)
    B.track(pt(*p2[-1]), pt(xs, y_ret_via), sig, width=tw, layer=pcbnew.In2_Cu)
    B.track(pt(xs, y_ret_via), pt(xs + 2.0, y_ret_via), sig, width=tw, layer=pcbnew.In2_Cu)
    B.via(pt(xs, y_sig_via), sig, size=0.5)
    B.via(pt(xs + 2.0, y_ret_via), sig, size=0.5)
    sigpad, gf = edge_sma(B, xs, sig, gnd)
    B.track(pt(xs, y_sig_via), sigpad, sig, width=0.4)
    nt = B.add_fp(FP_NETTIE, "NetTie-2_SMD_Pad0.5mm", B.ref("NT"), "loop return to GND", xs + 2.0, y_ret_via - 2.0, 90)
    ntp = sorted(nt.Pads(), key=lambda p: p.GetPosition().y)
    ntp[0].SetNet(gnd)
    ntp[1].SetNet(sig)
    B.track(pt(xs + 2.0, y_ret_via), ntp[1].GetPosition(), sig, width=0.4)
    # shield zones (F and B): whole board with a slit from the tip edge to the loop centre
    m, s = 0.3, 0.5
    poly = [(m, m), (cx - s, m), (cx - s, cy), (cx + s, cy), (cx + s, m), (w_board - m, m),
            (w_board - m, B.h - m), (m, B.h - m)]
    for lay in (pcbnew.F_Cu, pcbnew.B_Cu):
        zone(B, lay, gnd, poly)
    # stitching vias along both long edges, below the loop
    y = od + 4.0
    while y < B.h - 8.0:
        for x in (0.9, w_board - 0.9):
            B.via(pt(x, y), gnd, size=0.5)
        y += 4.0
    pcbnew.ZONE_FILLER(B.b).Fill(B.b.Zones())
    B.text(f"P{od:g}", cx, od + 8.0, 1.0)
    return B


def helmholtz():
    B = p3pcb.Builder("p3_helmholtz", 110.0, 110.0)
    B.outline(holes=False, fiducials=False)
    for i, (x, y) in enumerate(((7, 7), (103, 7), (7, 103), (103, 103))):
        B.add_fp(p3pcb.FP_HOLE, "MountingHole_3.2mm_M3", f"H{i + 1}", "M3 spacer", x, y)
    c = 55.0
    net = B.net("COIL")
    n, pitch, w, r0 = 10, 1.2, 0.9, 46.0
    th_end = 2 * math.pi * n + math.radians(10.0)
    nseg = int(n * 72 + 2)
    pts = []
    for k in range(nseg + 1):
        th = th_end * k / nseg
        r = r0 - pitch * th / (2 * math.pi)
        pts.append((c + r * math.cos(th), c + r * math.sin(th)))
    for a, b in zip(pts[:-1], pts[1:]):
        B.track(pt(*a), pt(*b), net, width=w)
    hdr = B.add_fp(FP_HDR, "PinHeader_1x02_P2.54mm_Vertical", B.ref("J"), "coil leads", 104.0, 55.0, 0, dnp=True)
    for p in hdr.Pads():
        p.SetNet(net)
    pins = sorted(hdr.Pads(), key=lambda p: p.GetPosition().y)
    B.track(pt(*pts[0]), pins[0].GetPosition(), net, width=w)
    # crossover: inner end radially out on the bottom layer, then to pin 2 on top
    th = math.radians(10.0)
    xo, yo = c + 52.0 * math.cos(th), c + 52.0 * math.sin(th)
    B.via(pt(*pts[-1]), net, size=0.8, drill=0.4)
    B.track(pt(*pts[-1]), pt(xo, yo), net, width=w, layer=pcbnew.B_Cu)
    B.via(pt(xo, yo), net, size=0.8, drill=0.4)
    B.track(pt(xo, yo), pins[1].GetPosition(), net, width=w)
    B.text("P3 Helmholtz coil: 10 turns, r 46 -> 34 mm (mean 40). Pair on 40 mm spacers.", c, 103.0, 1.0)
    B.text("Order 2 oz copper.", c, 106.0, 1.0)
    return B


if __name__ == "__main__":
    for B, layers in ((ztab(), 2), (loop_probe(10), 4), (loop_probe(20), 4), (helmholtz(), 2)):
        B.save(OUT)
        print(B.name, p3pcb.export(OUT, B.name, layers=layers))
