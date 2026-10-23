#!/usr/bin/env python3
"""Define, mesh and configure every case of the LGR-A driven study (run locally with gmsh >= 4.11).

  python build_cases.py <group> [--mesh] [--jobs N]

Groups
  sizing   eigenmode, FE order 2: baseline with/without PTFE spacers, coupling loops in the
           mid-height plane ("xy") of 3x3, 5x5, 8x6, 10x8 mm (depth x width), and one 5x5 loop in
           the radial-axial plane ("rz") for comparison of orientations
  tune     eigenmode, FE order 3: no loops, screw penetration 0, 2, 4, 6 mm, plus R0.4 lip fillets
  driven   2-port driven (adaptive fast sweep) for every loop size, order 2; the chosen size
           (4 x 4 mm) also at order 3 and on an h-refined mesh, full-order point solves, and the
           matching eigenmode runs. Needs results/sizing.json (analyze.py sizing).
  axial    exploratory: 6 mm brass rod through the lid on the axis (order 2)
  asbuilt  chosen loops + R0.4 mm lip fillets, driven, orders 2 and 3 (measurement prediction)
  (all groups above use the first, wrong copper input zs = "s2": half the copper loss; they are
   kept only to reproduce the first runs)

Corrected copper input (zs = "exact": Rs_p = 2 Rs, Ls_p = 2 Rs / w at the case's f0; names end
in "_cu2"):
  cu2_eig    eigenmode, order 3, no loops: sharp-edged reference (s = 10.42 mm) and the final
             build geometry (R0.4 lip fillets, PTFE spacers) at s = 10.42 / 10.55 / 10.70 mm,
             s = 10.55 mm with the two PEEK slab clips, and the empty-can check (order 2)
  cu2_loops  s = 10.55 mm as built: eigenmode with 50 ohm ports for loops of 4.1 / 4.5 / 4.7 /
             4.9 mm (order 2, Q_ext), the chosen loop at order 3 (eigenmode and driven)
  cu2_o4     s = 10.55 mm as built, no loops, order 4 (convergence of Q_c)
  cu2_fields s = 10.55 mm as built, no loops, order 2, fields of the fundamental mode saved
             (ParaView output; field maps made by fields_fig2.py)
Writes mesh/<mesh>.msh, configs/<case>.json, cases/<case>.json and lists/<group>.txt.
"""
import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import lgr_a_model as m

HERE = Path(__file__).resolve().parent
SIZES = [(3, 3), (5, 5), (8, 6), (10, 8)]
CHOSEN = None   # (w, h) in mm, set after the sizing round

# Expected f0 (Hz) used to evaluate Rs for the corrected runs. Sharp / R0.4 fillets at s = 10.42 mm:
# the order-3 results of the first runs (frequencies are not affected by the copper input). Other s:
# f0 ~ 1 / sqrt(s + s0 - ds) with s0 = 3.35 mm (size-sweep fit) and ds = 0.13 mm for the fillets.
F_FIL_1042 = 145.9088e6


def _f_fil(sv):
    return F_FIL_1042 * ((10.42 - 0.13 + 3.35) / (sv - 0.13 + 3.35)) ** 0.5


F_EST = {("sharp", 10.42): 145.2117e6, ("fil", 10.42): F_FIL_1042,
         ("fil", 10.55): _f_fil(10.55), ("fil", 10.70): _f_fil(10.70)}
# loop sizes (square, wire centreline, mm) for the corrected coupling: Q_ext does not depend on the
# copper input, so the first runs' beta(A) fit with the corrected Q_u (~2,130) predicts -30 dB at
# about 4.7 x 4.7 mm; 4.1 x 4.1 mm is the loop now drawn in the CAD
LOOPS_CU2 = (4.1, 4.5, 4.7, 4.9)
CHOSEN_CU2 = 4.7


def mesh_spec(mesh, **kw):
    return dict(mesh=mesh, **kw)


def cases(group):
    """Return a list of (case name, mesh name, mesh kwargs, config kwargs)."""
    out = []
    if group == "sizing":
        out.append(("e_base_nosp_o2", "base_nosp", dict(spacers=False), dict(order=2, loops=False, spacers=False)))
        out.append(("e_base_o2", "base", dict(), dict(order=2, loops=False)))
        for w, h in SIZES:
            mn = f"xy{w}x{h}"
            out.append((f"e_{mn}_o2", mn, dict(loops=dict(w=w, h=h, orient="xy")), dict(order=2)))
        out.append(("e_rz5x5_o2", "rz5x5", dict(loops=dict(w=5, h=5, orient="rz")), dict(order=2)))
        # first driven run: solver behaviour, cost, and a wide-band S-parameter check
        out.append(("d_xy5x5_o2_wide", "xy5x5", dict(loops=dict(w=5, h=5, orient="xy")),
                    dict(order=2, problem="Driven", f_lo=0.1440, f_hi=0.1458, f_step=2e-6)))
    elif group == "tune":
        out.append(("e_base_o3", "base", dict(), dict(order=3, loops=False)))
        for p in (2, 4, 6):
            out.append((f"e_screw{p}_o3", f"screw{p}", dict(screw=float(p)), dict(order=3, loops=False, screw=True)))
        out.append(("e_fillet_o3", "fillet04", dict(fillet=0.4), dict(order=3, loops=False)))
    elif group == "axial":
        # exploratory: M6 brass rod (6 mm) through the lid on the axis, 20 and 30 mm deep
        for pz in (20, 30):
            out.append((f"e_axial6_{pz}_o2", f"axial6_{pz}", dict(axial=(3.0, float(pz))),
                        dict(order=2, loops=False, screw=True)))
    elif group == "asbuilt":
        # recommended 4x4 loops + R0.4 lip fillets (geometry as drawn in hardware/lgr_a)
        lp = dict(loops=dict(w=4, h=4, orient="xy"), fillet=0.4)
        f0 = 0.1459105
        band = dict(f_lo=round(f0 - 0.0006, 6), f_hi=round(f0 + 0.0006, 6), f_step=1e-6)
        out.append(("d_xy4x4_fil_o3", "xy4x4_fil", lp, dict(order=3, problem="Driven", **band)))
        out.append(("d_xy4x4_fil_o2", "xy4x4_fil", lp, dict(order=2, problem="Driven", **band)))
    elif group == "driven":
        sz = json.loads((HERE / "results" / "sizing.json").read_text())
        chosen = tuple(sz["chosen_loop_mm"]) if CHOSEN is None else CHOSEN
        sizes = SIZES + ([chosen] if chosen not in SIZES else [])
        for w, h in sizes:
            mn = f"xy{w:g}x{h:g}".replace(".", "p")
            lp = dict(loops=dict(w=w, h=h, orient="xy"))
            f0 = sz["f0_est_GHz"].get(mn, sz["f0_base_GHz"])
            band = dict(f_lo=round(f0 - 0.0006, 6), f_hi=round(f0 + 0.0006, 6), f_step=1e-6)
            out.append((f"d_{mn}_o2", mn, lp, dict(order=2, problem="Driven", **band)))
            if (w, h) == chosen:
                # convergence: p-refinement (order 3, same mesh) and h-refinement (all sizes x0.8,
                # loop/port region x0.6 on top), plus full-order solves to check the reduced model
                mr = mn + "_r"
                lpr = dict(lp, h=0.8, hl=0.6)
                out.append((f"d_{mn}_o3", mn, lp, dict(order=3, problem="Driven", **band)))
                out.append((f"d_{mr}_o2", mr, lpr, dict(order=2, problem="Driven", **band)))
                pts = [round(f0 + df, 7) for df in (-0.0001, -0.000019, 0.0, 0.000019, 0.0001)]
                out.append((f"d_{mn}_o2_full", mn, lp, dict(order=2, problem="Driven", points=pts)))
                out.append((f"e_{mn}_o2", mn, lp, dict(order=2)))
                out.append((f"e_{mn}_o3", mn, lp, dict(order=3)))
                out.append((f"e_{mr}_o2", mr, lpr, dict(order=2)))
    elif group == "cu2_eig":
        out.append(("e_base_cu2_o3", "base", dict(), dict(order=3, loops=False, f_rs=F_EST[("sharp", 10.42)])))
        out.append(("e_fil_s1042_cu2_o3", "fillet04", dict(fillet=0.4),
                    dict(order=3, loops=False, f_rs=F_EST[("fil", 10.42)])))
        for sv in (10.55, 10.70):
            tag = f"s{round(sv * 100):04d}"
            out.append((f"e_fil_{tag}_cu2_o3", f"fil_{tag}", dict(fillet=0.4, s=sv),
                        dict(order=3, loops=False, f_rs=F_EST[("fil", sv)])))
        out.append(("e_fil_s1055_clip_cu2_o3", "fil_s1055_clip", dict(fillet=0.4, s=10.55, clips=True),
                    dict(order=3, loops=False, clips=True, f_rs=F_EST[("fil", 10.55)])))
        tm = m.empty_can_tm010(m.dims())
        out.append(("e_emptycan_cu2_o2", "emptycan", dict(empty=True),
                    dict(order=2, loops=False, spacers=False, empty=True, f_rs=tm["f_Hz"],
                         target=round(tm["f_Hz"] * 0.98e-9, 4), n_modes=2)))
    elif group == "cu2_loops":
        f0 = F_EST[("fil", 10.55)]
        for w in LOOPS_CU2:
            mn = f"fil_s1055_xy{w:g}x{w:g}".replace(".", "p")
            lp = dict(loops=dict(w=w, h=w, orient="xy"), fillet=0.4, s=10.55)
            out.append((f"e_{mn}_cu2_o2", mn, lp, dict(order=2, f_rs=f0)))
            if w == CHOSEN_CU2:
                band = dict(f_lo=round(f0 * 1e-9 - 0.001, 6), f_hi=round(f0 * 1e-9 + 0.001, 6), f_step=1e-6)
                out.append((f"e_{mn}_cu2_o3", mn, lp, dict(order=3, f_rs=f0)))
                out.append((f"d_{mn}_cu2_o3", mn, lp, dict(order=3, problem="Driven", f_rs=f0, **band)))
    elif group == "cu2_o4":
        out.append(("e_fil_s1055_cu2_o4", "fil_s1055", dict(fillet=0.4, s=10.55),
                    dict(order=4, loops=False, f_rs=F_EST[("fil", 10.55)])))
    elif group == "cu2_fields":
        # field maps for the paper (P1 Fig. 2): the default build without loops, order 2, with the
        # fundamental mode's fields written for ParaView (figures: fields_fig2.py)
        out.append(("e_fil_s1055_cu2_o2_fields", "fil_s1055", dict(fillet=0.4, s=10.55),
                    dict(order=2, loops=False, f_rs=F_EST[("fil", 10.55)], save=1)))
    else:
        raise SystemExit(f"unknown group {group}")
    if not group.startswith("cu2"):
        for c in out:
            c[3].setdefault("zs", "s2")      # reproduce the first runs (half the copper loss)
    return out


def do_mesh(args):
    mesh, kw = args
    kw = dict(kw)
    g = m.dims(kw.pop("s", m.S_LIP))
    if kw.pop("empty", False):
        info, ports = m.build_empty_can(g, HERE / "mesh" / f"{mesh}.msh")
    else:
        info, ports = m.build_mesh(g, HERE / "mesh" / f"{mesh}.msh", **kw)
    return mesh, info, ports


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("group")
    ap.add_argument("--mesh", action="store_true", help="(re)build the meshes of this group")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--force", action="store_true", help="rebuild meshes that already exist")
    a = ap.parse_args()
    cs = cases(a.group)
    meshes = {}
    for name, mesh, mkw, ckw in cs:
        meshes.setdefault(mesh, mkw)
    minfo = {}
    if a.mesh:
        todo = [(mn, kw) for mn, kw in meshes.items()
                if a.force or not (HERE / "mesh" / f"{mn}.msh").exists()]
        with ProcessPoolExecutor(a.jobs) as ex:
            for mn, info, ports in ex.map(do_mesh, todo):
                minfo[mn] = dict(info=info, ports=ports)
                print(mn, info["n_tet"], flush=True)
        for mn, v in minfo.items():
            (HERE / "mesh" / f"{mn}.json").write_text(json.dumps(v, indent=1))
    (HERE / "configs").mkdir(exist_ok=True)
    (HERE / "cases").mkdir(exist_ok=True)
    (HERE / "lists").mkdir(exist_ok=True)
    names = []
    for name, mesh, mkw, ckw in cs:
        g = m.dims(mkw.get("s", m.S_LIP))
        ckw = dict(ckw)
        ckw.setdefault("spacers", mkw.get("spacers", True))
        cfg = m.make_config(name, mesh, **ckw)
        (HERE / "configs" / f"{name}.json").write_text(json.dumps(cfg, indent=2) + "\n")
        mj = HERE / "mesh" / f"{mesh}.json"
        meta = dict(name=name, mesh=mesh, geometry=g, mesh_kw=mkw, config_kw=ckw,
                    mesh_info=json.loads(mj.read_text())["info"] if mj.exists() else None,
                    copper_input=ckw.get("zs", "exact"), f_Rs=ckw.get("f_rs", m.F0),
                    impedance=cfg["Boundaries"]["Impedance"],
                    Rs_cu=m.rs(m.SIGMA_CU, ckw.get("f_rs", m.F0)), Rs_brass=m.RS_BR, sigma_cu=m.SIGMA_CU, sigma_brass=m.SIGMA_BR if hasattr(m, "SIGMA_BR") else m.SIGMA_BRASS,
                    port=dict(Z=m.Z_PORT, wire_d=2 * m.R_WIRE, coax_outer_d=2 * m.B_COAX, recess=m.L_COAX))
        (HERE / "cases" / f"{name}.json").write_text(json.dumps(meta, indent=1) + "\n")
        names.append(name)
    (HERE / "lists" / f"{a.group}.txt").write_text("\n".join(names) + "\n")
    print("cases:", " ".join(names))


if __name__ == "__main__":
    sys.exit(main())
