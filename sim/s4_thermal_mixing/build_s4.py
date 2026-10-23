#!/usr/bin/env python3
"""Build the S4 meshes, Palace configs and job lists, round by round (local, gmsh >= 4.11).

  python build_s4.py <round> [--mesh-only | --config-only]

Rounds (each one writes lists/<round>*.txt):
  r1  single cells (S2 unit-cell box, PMC / PEC x-walls), TiO2-cell and composite-cell tuning,
      gap-width check, empty-box copper validation; all-alumina 2x1x2 tight box, 4-cell line and
      1x4x1 coaxial stack, single cell in its tight box (Q references)
  r2  count mixing k = 0..4 in the 2x1x2 tight box (all distinct arrangements), the 4-cell line
      with one TiO2 cell (end / interior), TiO2 cell in its own box
  r3  detuned TiO2 cell (continuous share), composite reference cell, order-3 confirmations
Naming: P1_* single cell in the unit-cell box (line N = 1, g = 0), B<nx><ny><nz>_<cells> tight copper
box, L4_<cells> 4-cell line with g = a/2, E112 empty box; run suffix _o<order>[_pmc|_pec]
[_dTp25 | _dTm25 | ...] (permittivities at T0 + dT, both dielectrics).
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import s4_geometry as g4  # noqa: E402

F_PMC, F_PEC, F_BOX, F_LINE, F_B111 = 149.92e6, 179.50e6, 172.5e6, 171.5e6, 186.6e6
DT_FD = (25.0, -25.0)


def tag(dT):
    return "" if dT == 0 else f"_dT{'p' if dT > 0 else 'm'}{abs(dT):g}".replace(".", "p")


class Builder:
    def __init__(self, mesh_only=False, config_only=False):
        self.mesh_only, self.config_only = mesh_only, config_only
        self.lists = {}

    def mesh(self, name, kind, cells, **kw):
        if self.config_only or (HERE / "mesh" / f"{name}.msh").exists() and (HERE / "cases" / f"{name}.json").exists():
            return name
        specs_s = [c for c in cells.split(",") if c]
        cp = g4.s3.cell_params(1.0, kw.pop("s", 6.2), kw.pop("d", None))
        gd = g4.s3.dims(kind, cp, kw.get("N", 1), kw.get("g", 0.0), kw.get("nx", 1), kw.get("ny", 1),
                        kw.get("nz", 1), kw.get("ext", 0.0))
        info = g4.build_mesh(gd, [g4.parse_spec(s) for s in specs_s], "v2", HERE / "mesh" / f"{name}.msh")
        meta = dict(name=name, geometry=gd, cells=specs_s, mesh=info, mat=g4.MAT)
        (HERE / "cases" / f"{name}.json").write_bytes((json.dumps(meta, indent=2) + "\n").encode())
        print(f"mesh {name}: {info['n_tet']} tets", flush=True)
        return name

    def run(self, lst, mesh, name, order=2, dT=0.0, **kw):
        if not self.mesh_only:
            case = json.loads((HERE / "cases" / f"{mesh}.json").read_text())
            cfg, meta = g4.make_config(name, mesh, case, order, dT=dT, **kw)
            meta["cells"] = case["cells"]
            (HERE / "configs" / f"{name}.json").write_bytes((json.dumps(cfg, indent=2) + "\n").encode())
            (HERE / "runs" / f"{name}.json").write_bytes((json.dumps(meta, indent=2) + "\n").encode())
        self.lists.setdefault(lst, []).append(name)

    def fd(self, lst, mesh, base, dts=(0.0,) + DT_FD, **kw):
        for dT in dts:
            self.run(lst, mesh, base + tag(dT), dT=dT, **kw)

    def write_lists(self):
        for k, v in self.lists.items():
            (HERE / "lists" / f"{k}.txt").write_bytes(("\n".join(v) + "\n").encode())   # LF only
            print(f"lists/{k}.txt: {len(v)} runs")


def r1(b):
    # single cells in the S2 unit-cell box (PMC x-walls: mu pole f_m; PEC: mu_eff = 0, f_mp)
    pmc = dict(xwall="pmc", target=0.14, n_modes=2, rs_f=F_PMC)
    pec = dict(xwall="pec", target=0.17, n_modes=2, rs_f=F_PEC)
    m = b.mesh("P1_A", "line", "A", N=1, g=0.0)
    b.fd("r1a", m, "P1_A_o2_pmc", **pmc)
    b.run("r1a", m, "P1_A_o2_pec", **pec)
    for l in ("3.0", "3.6", "4.2"):
        m = b.mesh(f"P1_T{l.replace('.', 'p')}", "line", f"T{l}", N=1, g=0.0)
        b.run("r1a", m, f"{m}_o2_pmc", **pmc)
    m = "P1_T3p6"
    b.run("r1a", m, f"{m}_o2_pmc{tag(25.0)}", dT=25.0, **pmc)
    b.run("r1a", m, f"{m}_o2_pmc{tag(-25.0)}", dT=-25.0, **pmc)
    b.run("r1a", m, f"{m}_o2_pec", **pec)
    for l in ("0.6", "1.0"):
        m = b.mesh(f"P1_C{l.replace('.', 'p')}", "line", f"C{l}", N=1, g=0.0)
        b.run("r1a", m, f"{m}_o2_pmc", **pmc)
    m = b.mesh("P1_A_d1p01", "line", "A", N=1, g=0.0, d=1.01)       # gap 1 % wider (expansion check)
    b.run("r1a", m, f"{m}_o2_pmc", **pmc)
    m = b.mesh("E112", "empty", "", nx=1, ny=1, nz=2)                 # copper-input validation
    b.run("r1a", m, "E112_o2", target=2.3, n_modes=2, rs_f=2.394e9)
    # all-alumina references
    m = b.mesh("B212_AAAA", "cavity", "A,A,A,A", nx=2, ny=1, nz=2)
    b.fd("r1b", m, "B212_AAAA_o2", target=0.16, n_modes=5, rs_f=F_BOX)
    m = b.mesh("L4_AAAA", "line", "A,A,A,A", N=4, g=35.0)
    b.fd("r1c", m, "L4_AAAA_o2", xwall="pec", target=0.12, n_modes=7, rs_f=F_LINE)
    m = b.mesh("B111_A", "cavity", "A", nx=1, ny=1, nz=1)
    b.run("r1c", m, "B111_A_o2", target=0.17, n_modes=2, rs_f=F_B111)
    m = b.mesh("B141_AAAA", "cavity", "A,A,A,A", nx=1, ny=4, nz=1)    # coaxial stack (lever 5)
    b.run("r1c", m, "B141_AAAA_o2", target=0.15, n_modes=4, rs_f=F_BOX)


def tuning():
    """TiO2 block length for a given cell f_m (fit 1/f^2 = a + b l over the r1 PMC runs) and the
    composite TiO2 length for tau_f = 0 (linear interpolation of the participation tau_f)."""
    import numpy as np
    import analyze_s4 as an
    fa = an.load("P1_A_o2_pmc")[0]["f_MHz"]
    L, y = [], []
    for l in ("3.0", "3.6", "4.2"):
        r = an.load(f"P1_T{l.replace('.', 'p')}_o2_pmc")[0]
        L.append(float(l))
        y.append(1 / r["f_MHz"] ** 2)
    b_, a_ = np.polyfit(L, y, 1)
    l_of = lambda rel: (1 / (fa * rel) ** 2 - a_) / b_      # rel = f_T / f_A
    lc, tc = [], []
    for l in ("0.6", "1.0"):
        r = an.load(f"P1_C{l.replace('.', 'p')}_o2_pmc")[0]
        lc.append(float(l))
        tc.append(an.tau_parts(r)["tau_pf"])
    l0 = lc[0] - tc[0] * (lc[1] - lc[0]) / (tc[1] - tc[0])
    return dict(f_A=fa, fit=(a_, b_), l_T=l_of, l_C0=l0)


def lstr(x):
    return f"{x:.2f}"


def tuning_w():
    """Slit width d of the full-length TiO2 cell (W) for a given f_m: fit 1/f^2 = a + b/d over the
    r1d PMC runs (C ~ 1/d)."""
    import numpy as np
    import analyze_s4 as an
    fa = an.load("P1_A_o2_pmc")[0]["f_MHz"]
    D, y = [], []
    for d in ("3.0", "3.6"):
        r = an.load(f"P1_W{d.replace('.', 'p')}_o2_pmc")[0]
        D.append(1 / float(d))
        y.append(1 / r["f_MHz"] ** 2)
    b_, a_ = np.polyfit(D, y, 1)
    return lambda rel: b_ / (1 / (fa * rel) ** 2 - a_)


def r2(b):
    """Mixing studies with the full-length TiO2 cell W<d> (chosen after r1/r1d, see README)."""
    d_of = tuning_w()
    W = f"W{lstr(d_of(1.0))}"
    C = "C0.75"                                  # composite cell with tau_f = 0 (r1d)
    print(f"TiO2 cell for f_T = f_A: {W}")
    pmc = dict(xwall="pmc", target=0.14, n_modes=3, rs_f=F_PMC)
    pec = dict(xwall="pec", target=0.17, n_modes=2, rs_f=F_PEC)
    box = dict(target=0.16, n_modes=5, rs_f=F_BOX)
    line = dict(xwall="pec", target=0.12, n_modes=7, rs_f=F_LINE)
    nm = lambda x: x.replace(".", "p")
    m = b.mesh(f"P1_{nm(W)}", "line", W, N=1, g=0.0)
    b.fd("r2a", m, f"{m}_o2_pmc", **pmc)
    b.run("r2a", m, f"{m}_o2_pec", **pec)
    m = b.mesh(f"B111_{nm(W)}", "cavity", W, nx=1, ny=1, nz=1)
    b.run("r2a", m, f"{m}_o2", target=0.17, n_modes=2, rs_f=F_B111)
    b.run("r2a", "P1_C0p75", "P1_C0p75_o2_pec", **pec)
    cellstr = lambda arr: ",".join({"T": W, "A": "A", "C": C}[c] for c in arr)
    # 4-cell line, one TiO2 cell at the end / in the interior
    for lst, arr in (("r2a", "TAAA"), ("r2b", "ATAA")):
        m = b.mesh(f"L4_{arr}", "line", cellstr(arr), N=4, g=35.0)
        b.fd(lst, m, f"{m}_o2", **line)
    # 2x1x2 tight box; cell order (x0,z0) bottom-left, (x0,z1) top-left, (x1,z0), (x1,z1)
    for lst, arr, fdq in (("r2b", "TTAA", False), ("r2b", "ATTT", False), ("r2b", "TATT", False),
                          ("r2c", "TAAA", True), ("r2c", "ATAA", True),
                          ("r2d", "TATA", True), ("r2d", "TAAT", True),
                          ("r2e", "TTTT", True), ("r2e", "CCCC", True), ("r2g", "ATAT", False)):
        m = b.mesh(f"B212_{arr}", "cavity", cellstr(arr), nx=2, ny=1, nz=2)
        b.fd(lst, m, f"{m}_o2", **box) if fdq else b.run(lst, m, f"{m}_o2", **box)
    # detuned TiO2 cell (continuous share), k = 1 at the bottom-left position
    for lst, rel in (("r2f", 1.03), ("r2f", 1.06), ("r2g", 1.10)):
        Wd = f"W{lstr(d_of(rel))}"
        m = b.mesh(f"B212_TAAA_d{round((rel - 1) * 1000):03d}", "cavity", f"{Wd},A,A,A", nx=2, ny=1, nz=2)
        b.fd(lst, m, f"{m}_o2", **box)
        print(f"detuned f_T/f_A = {rel}: {Wd}")


def r1d(b):
    """Full-length TiO2 slab cell (no lips) tuning; composite cell near tau_f = 0."""
    pmc = dict(xwall="pmc", target=0.14, n_modes=3, rs_f=F_PMC)
    for d in ("3.0", "3.6"):
        m = b.mesh(f"P1_W{d.replace('.', 'p')}", "line", f"W{d}", N=1, g=0.0)
        b.run("r1d", m, f"{m}_o2_pmc", **pmc)
    m = b.mesh("P1_C0p75", "line", "C0.75", N=1, g=0.0)
    b.fd("r1d", m, "P1_C0p75_o2_pmc", **dict(pmc, target=0.135))
    m = "P1_C1p0"
    b.run("r1d", m, "P1_C1p0_o2_pmc_t", **dict(pmc, target=0.13))


def r3(b):
    """Zero-crossing design (k = 1, TiO2 cell detuned to f_T/f_A ~ 1.018-1.021), fine FD in T, and
    order-3 confirmations."""
    box = dict(target=0.16, n_modes=5, rs_f=F_BOX)
    m = b.mesh("B212_TAAA_W3p41", "cavity", "W3.41,A,A,A", nx=2, ny=1, nz=2)
    b.fd("r3a", m, f"{m}_o2", dts=(0.0, 5.0, -5.0, 12.5, -12.5, 25.0, -25.0), **box)
    m = b.mesh("B212_TAAA_W3p43", "cavity", "W3.43,A,A,A", nx=2, ny=1, nz=2)
    b.fd("r3b", m, f"{m}_o2", dts=(0.0, 5.0, -5.0, 25.0, -25.0), **box)
    b.run("r3c", "B212_AAAA", "B212_AAAA_o3", order=3, **box)
    b.fd("r3c", "B212_TAAA_W3p41", "B212_TAAA_W3p41_o3", dts=(0.0, 5.0, -5.0), order=3, **box)
    pmc = dict(xwall="pmc", target=0.14, n_modes=3, rs_f=F_PMC)
    b.run("r3c", "P1_A", "P1_A_o3_pmc", order=3, **pmc)
    b.run("r3c", "P1_W3p29", "P1_W3p29_o3_pmc", order=3, **pmc)


ROUNDS = {"r1": r1, "r1d": r1d, "r2": r2, "r3": r3}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("round", choices=sorted(ROUNDS))
    ap.add_argument("--mesh-only", action="store_true")
    ap.add_argument("--config-only", action="store_true")
    a = ap.parse_args()
    for d in ("mesh", "cases", "configs", "runs", "lists"):
        (HERE / d).mkdir(exist_ok=True)
    b = Builder(a.mesh_only, a.config_only)
    ROUNDS[a.round](b)
    b.write_lists()


if __name__ == "__main__":
    main()
