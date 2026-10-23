"""Shared helpers for the P3 boards (run inside KiCad's Python).

Board builder, the three drum-inductor cell types, printed square spirals,
and the export chain (DRC, Gerbers, drill, position file, JLCPCB BOM/CPL, renders).
"""
import collections
import csv
import glob
import json
import os
import subprocess
import zipfile

import pcbnew

HERE = os.path.dirname(os.path.abspath(__file__))
FP_LOCAL = os.path.join(HERE, "footprints.pretty")
KICAD = "C:/Program Files/KiCad/10.0"
FP_ROOT = os.environ.get("KICAD10_FOOTPRINT_DIR", KICAD + "/share/kicad/footprints")
FP_CAP = os.path.join(FP_ROOT, "Capacitor_SMD.pretty")
FP_FID = os.path.join(FP_ROOT, "Fiducial.pretty")
FP_HOLE = os.path.join(FP_ROOT, "MountingHole.pretty")
FP_COAX = os.path.join(FP_ROOT, "Connector_Coaxial.pretty")
KICAD_CLI = KICAD + "/bin/kicad-cli.exe"

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
CAP_220N = dict(c_fp="C_1206_3216Metric", c_val="220nF C0G", c_mpn="GRM31C5C1H224JE02L", c_lcsc="C464977")

mm = pcbnew.FromMM


def pt(x, y):
    return pcbnew.VECTOR2I(mm(x), mm(y))


class Builder:
    def __init__(self, name, w, h, layers=2):
        self.name = name
        self.w, self.h = w, h
        self.b = pcbnew.BOARD()
        self.b.SetCopperLayerCount(layers)
        self.parts = []      # assembled parts (BOM/CPL)
        self.cells = []      # measurement map
        self.n = collections.Counter()
        self.trace_w = 0.8

    # --- primitives -------------------------------------------------------
    def ref(self, prefix):
        self.n[prefix] += 1
        return f"{prefix}{self.n[prefix]}"

    def net(self, name):
        n = pcbnew.NETINFO_ITEM(self.b, name)
        self.b.Add(n)
        return n

    def add_fp(self, lib, name, ref, val, x, y, rot=0, dnp=False, layer_back=False):
        fp = pcbnew.FootprintLoad(lib, name)
        fp.SetReference(ref)
        fp.SetValue(val)
        fp.SetPosition(pt(x, y))
        fp.SetOrientationDegrees(rot)
        fp.Reference().SetVisible(False)
        if layer_back:
            fp.Flip(pt(x, y), pcbnew.FLIP_DIRECTION_LEFT_RIGHT)
        if dnp:
            fp.SetDNP(True)
            fp.SetExcludedFromBOM(True)
            fp.SetExcludedFromPosFiles(True)
        self.b.Add(fp)
        return fp

    def place_part(self, lib, fpname, prefix, part, x, y, rot=0, dnp=False, val=None):
        ref = self.ref(prefix)
        fp = self.add_fp(lib, fpname, ref, val or part.get("val", ""), x, y, rot, dnp)
        if not dnp and part.get("lcsc"):
            self.parts.append(dict(ref=ref, val=val or part["val"], fp=fpname, mpn=part["mpn"], lcsc=part["lcsc"]))
        return fp

    def track(self, a, b, net, width=None, layer=pcbnew.F_Cu):
        t = pcbnew.PCB_TRACK(self.b)
        t.SetStart(a)
        t.SetEnd(b)
        t.SetWidth(mm(width or self.trace_w))
        t.SetLayer(layer)
        t.SetNet(net)
        self.b.Add(t)

    def via(self, p, net, drill=0.3, size=0.6):
        v = pcbnew.PCB_VIA(self.b)
        v.SetPosition(p)
        v.SetDrill(mm(drill))
        v.SetWidth(mm(size))
        v.SetNet(net)
        self.b.Add(v)

    def text(self, s, x, y, size=1.2, layer=pcbnew.F_SilkS):
        t = pcbnew.PCB_TEXT(self.b)
        t.SetText(s)
        t.SetPosition(pt(x, y))
        t.SetLayer(layer)
        t.SetTextSize(pcbnew.VECTOR2I(mm(size), mm(size)))
        t.SetTextThickness(mm(0.15 * size))
        self.b.Add(t)

    def seg(self, x0, y0, x1, y1, layer=pcbnew.F_SilkS, width=0.15):
        s = pcbnew.PCB_SHAPE(self.b)
        s.SetShape(pcbnew.SHAPE_T_SEGMENT)
        s.SetStart(pt(x0, y0))
        s.SetEnd(pt(x1, y1))
        s.SetLayer(layer)
        s.SetWidth(mm(width))
        self.b.Add(s)

    def outline(self, holes=True, fiducials=True):
        s = pcbnew.PCB_SHAPE(self.b)
        s.SetShape(pcbnew.SHAPE_T_RECT)
        s.SetStart(pt(0, 0))
        s.SetEnd(pt(self.w, self.h))
        s.SetLayer(pcbnew.Edge_Cuts)
        s.SetWidth(mm(0.1))
        self.b.Add(s)
        if holes:
            for i, (x, y) in enumerate(((5, 5), (self.w - 5, 5), (5, self.h - 5), (self.w - 5, self.h - 5))):
                self.add_fp(FP_HOLE, "MountingHole_3.2mm_M3", f"H{i + 1}", "M3 jig hole", x, y)
        if fiducials:
            for i, (x, y) in enumerate(((14, 5), (self.w - 14, 5), (14, self.h - 5))):
                self.add_fp(FP_FID, "Fiducial_1mm_Mask2mm", f"FID{i + 1}", "Fiducial", x, y)

    # --- drum-inductor cell -------------------------------------------------
    def cell(self, cid, typ, x, y, with_c=True, trim=True, dx1=None, group="", label="top"):
        """Closed LC loop centred on the inductor; capacitor(s) to the right (+x)."""
        T = TYPES[typ]
        fl = self.place_part(FP_LOCAL, T["l_fp"], "L",
                             dict(val=T["l_val"], mpn=T["l_mpn"], lcsc=T["l_lcsc"]), x, y, T["l_rot"])
        pads_l = sorted(fl.Pads(), key=lambda p: p.GetPosition().y)
        na, nb = self.net(f"{cid}_a"), self.net(f"{cid}_b")
        pads_l[0].SetNet(na)
        pads_l[1].SetNet(nb)
        c1_ref = None
        if with_c:
            offs = [dx1 if dx1 is not None else T["dx1"]] + ([T["dx2"]] if trim else [])
            caps = []
            for k, dx in enumerate(offs):
                dnp = k == 1
                fc = self.place_part(FP_CAP, T["c_fp"], "C",
                                     dict(val=T["c_val"], mpn=T["c_mpn"], lcsc=T["c_lcsc"]),
                                     x + dx, y, 90, dnp=dnp, val=None if not dnp else "DNP trim")
                if k == 0:
                    c1_ref = fc.GetReference()
                pc = sorted(fc.Pads(), key=lambda p: p.GetPosition().y)
                pc[0].SetNet(na)
                pc[1].SetNet(nb)
                caps.append(pc)
            self.track(pads_l[0].GetPosition(), caps[0][0].GetPosition(), na)
            self.track(pads_l[1].GetPosition(), caps[0][1].GetPosition(), nb)
            if len(caps) > 1:
                self.track(caps[0][0].GetPosition(), caps[1][0].GetPosition(), na)
                self.track(caps[0][1].GetPosition(), caps[1][1].GetPosition(), nb)
        if label == "top":
            self.text(cid, x, y - 8.0)
        elif label == "left":
            self.text(cid, x - 14.0, y)
        if label != "none":
            self.seg(x - 7.5, y, x - 6.0, y)
        self.cells.append(dict(id=cid, type=typ, x_mm=x, y_mm=y, has_c=with_c, group=group,
                               inductor=fl.GetReference(), c1=c1_ref))

    # --- printed square spiral ---------------------------------------------
    def spiral(self, net_outer, net_inner, x, y, d_out, n, w, s, layer=pcbnew.F_Cu):
        """Square spiral, centre-line, starting at the top-left outer corner and winding inwards.
        Returns (outer end point, inner end point) in mm."""
        p = w + s
        H = (d_out - w) / 2
        pts = [(-H, -H)]
        for k in range(n):
            a = H - k * p
            pts += [(a, -a), (a, a), (-a, a), (-a, -a + p)]
        pts = [(x + u, y + v) for u, v in pts]
        net = net_outer
        for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
            self.track(pt(x0, y0), pt(x1, y1), net, width=w, layer=layer)
        return pts[0], pts[-1]

    # --- outputs -------------------------------------------------------------
    def save(self, outdir):
        os.makedirs(outdir, exist_ok=True)
        path = os.path.join(outdir, self.name + ".kicad_pcb")
        pcbnew.SaveBoard(path, self.b)
        with open(os.path.join(outdir, self.name + "_map.json"), "w") as f:
            json.dump(self.cells, f, indent=1)
        with open(os.path.join(outdir, self.name + "_parts.json"), "w") as f:
            json.dump(self.parts, f, indent=1)
        return path


def run(*args):
    r = subprocess.run([KICAD_CLI, *args], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr or r.stdout)
    return r.stdout


def export(outdir, name, layers=2, render=True):
    """DRC + Gerbers + drill + position file + JLCPCB BOM/CPL + renders. Returns DRC summary."""
    pcb = os.path.join(outdir, name + ".kicad_pcb")
    drc = os.path.join(outdir, name + "_drc.json")
    run("pcb", "drc", "--format", "json", "-o", drc, pcb)
    d = json.load(open(drc))
    summary = dict(violations=collections.Counter((v["severity"], v["type"]) for v in d["violations"]),
                   unconnected=len(d["unconnected_items"]))
    jlc = os.path.join(outdir, "jlc_" + name)
    gdir = os.path.join(jlc, "gerbers")
    os.makedirs(gdir, exist_ok=True)
    cu = "F.Cu,B.Cu" if layers == 2 else "F.Cu,In1.Cu,In2.Cu,B.Cu"
    run("pcb", "export", "gerbers", "-l", cu + ",F.Paste,F.SilkS,B.SilkS,F.Mask,B.Mask,Edge.Cuts", "-o", gdir + "/", pcb)
    run("pcb", "export", "drill", "--format", "excellon", "--excellon-separate-th", "-o", gdir + "/", pcb)
    with zipfile.ZipFile(os.path.join(jlc, f"gerbers_{name}.zip"), "w", zipfile.ZIP_DEFLATED) as z:
        for fn in glob.glob(os.path.join(gdir, "*")):
            z.write(fn, os.path.basename(fn))
    parts = {p["ref"]: p for p in json.load(open(os.path.join(outdir, name + "_parts.json")))}
    if parts:
        pos = os.path.join(outdir, name + "_pos.csv")
        run("pcb", "export", "pos", "--format", "csv", "--units", "mm", "--side", "front", "-o", pos, pcb)
        rows = [r for r in csv.DictReader(open(pos)) if r["Ref"] in parts]
        with open(os.path.join(jlc, f"CPL_{name}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Designator", "Mid X", "Mid Y", "Layer", "Rotation"])
            for r in rows:
                w.writerow([r["Ref"], r["PosX"] + "mm", r["PosY"] + "mm", "Top", r["Rot"]])
        g = collections.OrderedDict()
        for ref, p in parts.items():
            g.setdefault((p["val"], p["fp"], p["lcsc"], p["mpn"]), []).append(ref)
        with open(os.path.join(jlc, f"BOM_{name}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Comment", "Designator", "Footprint", "LCSC Part #", "Manufacturer Part", "Qty"])
            for (val, fp, lcsc, mpn), refs in g.items():
                w.writerow([val, ",".join(refs), fp, lcsc, mpn, len(refs)])
        summary["placements"] = len(rows)
        summary["bom"] = [(k[0], k[2], len(v)) for k, v in g.items()]
    if render:
        run("pcb", "render", "--side", "top", "--width", "1600", "--height", "1300", "--quality", "basic",
            "-o", os.path.join(outdir, name + "_top.png"), pcb)
        run("pcb", "render", "--width", "1800", "--height", "1200", "--rotate", "-50,0,25", "--zoom", "1.1",
            "--quality", "high", "-o", os.path.join(outdir, name + "_iso.png"), pcb)
    return summary
