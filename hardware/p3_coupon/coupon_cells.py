"""P3 cell coupon v1 (board 1): generator. Run with KiCad's Python:

  "C:/Program Files/KiCad/10.0/bin/python.exe" coupon_cells.py

A bare FR-4 board (no copper pours) with isolated LC cells for single-cell f0/Q and
their spread, coplanar cell pairs at three spacings (in-situ coupling check), and
L-only cells (no capacitor) as non-resonant controls. Each cell is a closed loop of one
unshielded ferrite drum inductor and one C0G capacitor, plus a DNP pad for a trim cap.
Cell types A/B/C are defined in p3pcb.TYPES. Writes out/ (board, DRC, JLCPCB files, renders).
"""
import os

import p3pcb

OUT = os.path.join(p3pcb.HERE, "out")
NAME = "p3_cell_coupon"


def build():
    B = p3pcb.Builder(NAME, 225.0, 210.0)
    B.outline()

    # 1) isolated cells: 7 per type at 30 mm pitch (neighbour k <~ 1e-3)
    for r, typ in enumerate("ABC"):
        for k in range(7):
            B.cell(f"{typ}{k + 1}", typ, 20 + 30 * k, 22 + 30 * r, group="stat")

    # 2) coplanar pairs, centre spacing d along y; k ~ d^-3, so keep the splitting above 1/Q
    spacings = (11.5, 14.0, 18.0)
    col = 0
    for typ in "AB":
        for d in spacings:
            x = 20 + 35 * col
            B.cell(f"{typ}P{d:g}a", typ, x, 117.0, group=f"pair{typ}{d:g}", label="left")
            B.cell(f"{typ}P{d:g}b", typ, x, 117.0 + d, group=f"pair{typ}{d:g}", label="left")
            col += 1
    for k, d in enumerate(spacings):
        x = 20 + 35 * k
        B.cell(f"CP{d:g}a", "C", x, 168.0, group=f"pairC{d:g}", label="left")
        B.cell(f"CP{d:g}b", "C", x, 168.0 + d, group=f"pairC{d:g}", label="left")

    # 3) L-only controls
    for k, typ in enumerate("ABC"):
        B.cell(f"{typ}L", typ, 130 + 30 * k, 176.0, with_c=False, group="L-only")

    B.text("P3 cell coupon v1  (metawatchar, 2026)", 112, B.h - 9, 1.6)
    B.text("Bare FR-4, no copper pours. Rows A/B/C: isolated cells, 30 mm pitch. "
           "xPdd: pairs at dd mm. xL: L only.", 112, B.h - 5.5, 1.0)
    return B


if __name__ == "__main__":
    B = build()
    B.save(OUT)
    print(NAME, len(B.cells), "cells", len(B.parts), "parts")
    print(p3pcb.export(OUT, NAME))
