"""Write JLCPCB BOM/CPL and zip the Gerbers from KiCad exports (run in out/)."""
import collections, csv, glob, json, os, zipfile
parts = {p["ref"]: p for p in json.load(open("parts.json"))}
rows = [r for r in csv.DictReader(open("pos.csv")) if r["Ref"] in parts]
with open("jlc/CPL_p3_cell_coupon.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["Designator", "Mid X", "Mid Y", "Layer", "Rotation"])
    for r in rows: w.writerow([r["Ref"], r["PosX"] + "mm", r["PosY"] + "mm", "Top", r["Rot"]])
g = collections.OrderedDict()
for ref, p in parts.items(): g.setdefault((p["val"], p["fp"], p["lcsc"], p["mpn"]), []).append(ref)
with open("jlc/BOM_p3_cell_coupon.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["Comment", "Designator", "Footprint", "LCSC Part #", "Manufacturer Part", "Qty"])
    for (val, fp, lcsc, mpn), refs in g.items(): w.writerow([val, ",".join(refs), fp, lcsc, mpn, len(refs)])
with zipfile.ZipFile("jlc/gerbers_p3_cell_coupon.zip", "w", zipfile.ZIP_DEFLATED) as z:
    for fn in glob.glob("jlc/gerbers/*"): z.write(fn, os.path.basename(fn))
print(len(rows), "placements;", [(k[0], k[2], len(v)) for k, v in g.items()])
