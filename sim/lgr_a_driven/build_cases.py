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
    else:
        raise SystemExit(f"unknown group {group}")
    return out


def do_mesh(args):
    mesh, kw = args
    g = m.dims()
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
    g = m.dims()
    names = []
    for name, mesh, mkw, ckw in cs:
        ckw = dict(ckw)
        ckw.setdefault("spacers", mkw.get("spacers", True))
        cfg = m.make_config(name, mesh, **ckw)
        (HERE / "configs" / f"{name}.json").write_text(json.dumps(cfg, indent=2) + "\n")
        mj = HERE / "mesh" / f"{mesh}.json"
        meta = dict(name=name, mesh=mesh, geometry=g, mesh_kw=mkw, config_kw=ckw,
                    mesh_info=json.loads(mj.read_text())["info"] if mj.exists() else None,
                    Rs_cu=m.RS_CU, Rs_brass=m.RS_BR, sigma_cu=m.SIGMA_CU, sigma_brass=m.SIGMA_BR if hasattr(m, "SIGMA_BR") else m.SIGMA_BRASS,
                    port=dict(Z=m.Z_PORT, wire_d=2 * m.R_WIRE, coax_outer_d=2 * m.B_COAX, recess=m.L_COAX))
        (HERE / "cases" / f"{name}.json").write_text(json.dumps(meta, indent=1) + "\n")
        names.append(name)
    (HERE / "lists" / f"{a.group}.txt").write_text("\n".join(names) + "\n")
    print("cases:", " ".join(names))


if __name__ == "__main__":
    sys.exit(main())
