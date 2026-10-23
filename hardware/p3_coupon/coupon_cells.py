"""P3 cell coupon v1: board generator (run with KiCad's Python).

Builds `out/p3_cell_coupon.kicad_pcb`: a bare FR-4 board (no copper pours) carrying
isolated LC cells for measuring single-cell f0/Q and their spread, coplanar cell pairs
at three spacings for the inter-cell mutual inductance (m_eff), and L-only cells
(no capacitor) as controls. Every cell is a closed loop: one unshielded ferrite drum
inductor and one C0G capacitor, plus an unpopulated pad for a trim capacitor.

Cell types (parts chosen from JLCPCB stock, 2026-10):
  A  Sumida CD105T125NP-470KC  47 uH, 10.4 mm  + 68 nF C0G 1206  -> ~89 kHz
  B  SXN SM7850-101MT         100 uH,  7.8 mm  + 33 nF C0G 0805  -> ~88 kHz
  C  Bourns SDR0805-101KL     100 uH,  7.8 mm  + 33 nF C0G 0805  -> ~88 kHz

Run (then DRC and exports, from out/):
  "C:/Program Files/KiCad/10.0/bin/python.exe" coupon_cells.py
  kicad-cli pcb drc --format json -o drc.json p3_cell_coupon.kicad_pcb
  kicad-cli pcb export gerbers -l "F.Cu,B.Cu,F.Paste,F.SilkS,B.SilkS,F.Mask,B.Mask,Edge.Cuts" -o jlc/gerbers/ p3_cell_coupon.kicad_pcb
  kicad-cli pcb export drill --format excellon --excellon-separate-th -o jlc/gerbers/ p3_cell_coupon.kicad_pcb
  kicad-cli pcb export pos --format csv --units mm --side front -o pos.csv p3_cell_coupon.kicad_pcb
  python make_jlc.py
"""
import json
import os

import pcbnew

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FP_LOCAL = os.path.join(HERE, "footprints.pretty")
FP_CAP = os.path.join(os.environ.get("KICAD10_FOOTPRINT_DIR",
                                      "C:/Program Files/KiCad/10.0/share/kicad/footprints"),
                      "Capacitor_SMD.pretty")
FP_FID = os.path.join(os.path.dirname(FP_CAP), "Fiducial.pretty")
FP_HOLE = os.path.join(os.path.dirname(FP_CAP), "MountingHole.pretty")

BOARD_W, BOARD_H = 225.0, 210.0
TRACE_W = 0.8

# type -> inductor footprint, rotation, part, LCSC; capacitor footprint, part, LCSC; offsets of C1/C2
TYPES = {
    "A": dict(l_fp="L_Sumida_CD105T125", l_rot=0, l_val="47uH", l_mpn="CD105T125NP-470KC", l_lcsc="C171573",
              c_fp="C_1206_3216Metric", c_val="68nF C0G", c_mpn="GRM3195C1H683JA05D", c_lcsc="C3847777",
              dx1=7.5, dx2=11.0),
    "B": dict(l_fp="L_SXN_SM7850", l_rot=90, l_val="100uH", l_mpn="SM7850-101MT", l_lcsc="C10155",
              c_fp="C_0805_2012Metric", c_val="33nF C0G", c_mpn="C2012C0G1H333JT000N", c_lcsc="C76628",
              dx1=6.0, dx2=9.0),
    "C": dict(l_fp="L_Bourns_SDR0805", l_rot=0, l_val="100uH", l_mpn="SDR0805-101KL", l_lcsc="C139487",
              c_fp="C_0805_2012Metric", c_val="33nF C0G", c_mpn="C2012C0G1H333JT000N", c_lcsc="C76628",
              dx1=6.0, dx2=9.0),
}

mm = pcbnew.FromMM


def pt(x, y):
    return pcbnew.VECTOR2I(mm(x), mm(y))


class Builder:
    def __init__(self):
        self.b = pcbnew.BOARD()
        self.b.SetCopperLayerCount(2)
        self.parts = []      # for BOM/CPL and the README table
        self.cells = []      # measurement map
        self.n = {"L": 0, "C": 0}

    def add_fp(self, lib, name, ref, val, x, y, rot=0, dnp=False):
        fp = pcbnew.FootprintLoad(lib, name)
        fp.SetReference(ref)
        fp.SetValue(val)
        fp.SetPosition(pt(x, y))
        fp.SetOrientationDegrees(rot)
        fp.Reference().SetVisible(False)   # cell IDs are printed instead
        if dnp:
            fp.SetDNP(True)
            fp.SetExcludedFromBOM(True)
            fp.SetExcludedFromPosFiles(True)
        self.b.Add(fp)
        return fp

    def track(self, a, b, net):
        t = pcbnew.PCB_TRACK(self.b)
        t.SetStart(a)
        t.SetEnd(b)
        t.SetWidth(mm(TRACE_W))
        t.SetLayer(pcbnew.F_Cu)
        t.SetNet(net)
        self.b.Add(t)

    def text(self, s, x, y, size=1.2, layer=pcbnew.F_SilkS):
        t = pcbnew.PCB_TEXT(self.b)
        t.SetText(s)
        t.SetPosition(pt(x, y))
        t.SetLayer(layer)
        t.SetTextSize(pcbnew.VECTOR2I(mm(size), mm(size)))
        t.SetTextThickness(mm(0.15 * size))
        self.b.Add(t)

    def cell(self, cid, typ, x, y, with_c=True, group="", label_left=False):
        """One LC cell centred on the inductor; capacitors to the right (+x)."""
        T = TYPES[typ]
        self.n["L"] += 1
        lref = f"L{self.n['L']}"
        fl = self.add_fp(FP_LOCAL, T["l_fp"], lref, T["l_val"], x, y, T["l_rot"])
        self.parts.append(dict(ref=lref, val=T["l_val"], fp=T["l_fp"], mpn=T["l_mpn"], lcsc=T["l_lcsc"]))
        pads_l = sorted(fl.Pads(), key=lambda p: p.GetPosition().y)    # top, bottom
        net_a = pcbnew.NETINFO_ITEM(self.b, f"{cid}_a")
        net_b = pcbnew.NETINFO_ITEM(self.b, f"{cid}_b")
        self.b.Add(net_a)
        self.b.Add(net_b)
        pads_l[0].SetNet(net_a)
        pads_l[1].SetNet(net_b)
        caps = []
        if with_c:
            for k, dx in enumerate((T["dx1"], T["dx2"])):
                self.n["C"] += 1
                cref = f"C{self.n['C']}"
                dnp = k == 1
                fc = self.add_fp(FP_CAP, T["c_fp"], cref, T["c_val"] if not dnp else "DNP trim", x + dx, y, 90, dnp)
                if not dnp:
                    self.parts.append(dict(ref=cref, val=T["c_val"], fp=T["c_fp"], mpn=T["c_mpn"], lcsc=T["c_lcsc"]))
                pc = sorted(fc.Pads(), key=lambda p: p.GetPosition().y)
                pc[0].SetNet(net_a)
                pc[1].SetNet(net_b)
                caps.append(pc)
            # inductor pad centre -> C1 pad, then C1 -> C2 straight across
            self.track(pads_l[0].GetPosition(), caps[0][0].GetPosition(), net_a)
            self.track(pads_l[1].GetPosition(), caps[0][1].GetPosition(), net_b)
            self.track(caps[0][0].GetPosition(), caps[1][0].GetPosition(), net_a)
            self.track(caps[0][1].GetPosition(), caps[1][1].GetPosition(), net_b)
        if label_left:
            self.text(cid, x - 14.0, y)
        else:
            self.text(cid, x, y - 8.0)
        # tick on the inductor axis (y), used to check the probe jig; x and y come from the jig holes
        for (x0, y0, x1, y1) in ((x - 7.5, y, x - 6.0, y),):
            s = pcbnew.PCB_SHAPE(self.b)
            s.SetShape(pcbnew.SHAPE_T_SEGMENT)
            s.SetStart(pt(x0, y0))
            s.SetEnd(pt(x1, y1))
            s.SetLayer(pcbnew.F_SilkS)
            s.SetWidth(mm(0.15))
            self.b.Add(s)
        self.cells.append(dict(id=cid, type=typ, x_mm=x, y_mm=y, has_c=with_c, group=group,
                               inductor=lref, c1=None if not with_c else f"C{self.n['C'] - 1}"))

    def outline(self):
        s = pcbnew.PCB_SHAPE(self.b)
        s.SetShape(pcbnew.SHAPE_T_RECT)
        s.SetStart(pt(0, 0))
        s.SetEnd(pt(BOARD_W, BOARD_H))
        s.SetLayer(pcbnew.Edge_Cuts)
        s.SetWidth(mm(0.1))
        self.b.Add(s)
        for i, (x, y) in enumerate(((5, 5), (BOARD_W - 5, 5), (5, BOARD_H - 5), (BOARD_W - 5, BOARD_H - 5))):
            h = self.add_fp(FP_HOLE, "MountingHole_3.2mm_M3", f"H{i + 1}", "M3 jig hole", x, y)
        for i, (x, y) in enumerate(((14, 5), (BOARD_W - 14, 5), (14, BOARD_H - 5))):
            self.add_fp(FP_FID, "Fiducial_1mm_Mask2mm", f"FID{i + 1}", "Fiducial", x, y)


def main():
    os.makedirs(OUT, exist_ok=True)
    B = Builder()
    B.outline()

    # 1) statistics rows: 7 isolated cells per type at 30 mm pitch
    for r, typ in enumerate("ABC"):
        for k in range(7):
            B.cell(f"{typ}{k + 1}", typ, 20 + 30 * k, 22 + 30 * r, group="stat")

    # 2) coplanar pairs (centre-to-centre spacing d along y), types A and B
    spacings = (11.5, 14.0, 18.0)   # k ~ d^-3: keep the pair splitting above the linewidth (see README)
    y0 = 117.0
    col = 0
    for typ in "AB":
        for d in spacings:
            x = 20 + 35 * col
            B.cell(f"{typ}P{d:g}a", typ, x, y0, group=f"pair{typ}{d:g}", label_left=True)
            B.cell(f"{typ}P{d:g}b", typ, x, y0 + d, group=f"pair{typ}{d:g}", label_left=True)
            col += 1

    # 3) type C pairs and L-only controls
    y1 = 168.0
    for k, d in enumerate(spacings):
        x = 20 + 35 * k
        B.cell(f"CP{d:g}a", "C", x, y1, group=f"pairC{d:g}", label_left=True)
        B.cell(f"CP{d:g}b", "C", x, y1 + d, group=f"pairC{d:g}", label_left=True)
    for k, typ in enumerate("ABC"):
        B.cell(f"{typ}L", typ, 130 + 30 * k, y1 + 8, with_c=False, group="L-only")

    B.text("P3 cell coupon v1  (metawatchar, 2026)", 112, BOARD_H - 9, 1.6)
    B.text("Bare FR-4, no copper pours. Rows A/B/C: isolated cells, 30 mm pitch. "
           "xPdd: pairs at dd mm. xL: L only.", 112, BOARD_H - 5.5, 1.0)

    path = os.path.join(OUT, "p3_cell_coupon.kicad_pcb")
    pcbnew.SaveBoard(path, B.b)
    with open(os.path.join(OUT, "cell_map.json"), "w") as f:
        json.dump(B.cells, f, indent=1)
    with open(os.path.join(OUT, "parts.json"), "w") as f:
        json.dump(B.parts, f, indent=1)
    print(path, len(B.cells), "cells", len(B.parts), "placed parts")


if __name__ == "__main__":
    main()
