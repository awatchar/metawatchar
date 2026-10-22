#!/usr/bin/env python3
"""Generate meshes/configs for a list of (V, a, d, s) runs and write lists/<round>_<k>.txt.

  python make_round.py [--tag _v2] <round> <n_lists> V,a,d,s[,h[,order[,pec]]] ...
      (s in mm; "max" = s_max). Example (production case): --tag _v2 f 1 0.25,0.5,0.635,10.92,1.0,3
Case names: v<V*100>_a<a*10>_d<d*1000>_s<s*100>[_h<h*100>][_o<order>]<tag>[_pec]; the mesh is shared
by every order of the same geometry (mesh/v..._s....[_h...]<tag>.msh). The tuning rounds of the sweep
were made without a tag (v1 mesh recipe); the production runs with --tag _v2 (current recipe).
"""
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def s_max(V, a):
    R = (V * 1e-3 / (2 * math.pi)) ** (1 / 3) * 1e3
    return math.floor((R - a * R - 3.0 - 3.0) * 100) / 100


def case_name(V, a, d, s, h=1.0, order=2):
    n = f"v{round(V*100):03d}_a{round(a*10):02d}_d{round(d*1000):04d}_s{round(s*100):04d}"
    if h != 1.0:
        n += f"_h{round(h*100):03d}"
    if order != 2:
        n += f"_o{order}"
    return n


def main():
    argv = sys.argv[1:]
    tag = ""
    if argv[0] == "--tag":
        tag, argv = argv[1], argv[2:]
    rnd, nlists = argv[0], int(argv[1])
    names = []
    for spec in argv[2:]:
        parts = spec.split(",")
        V, a, d = map(float, parts[:3])
        s = s_max(V, a) if parts[3] == "max" else float(parts[3])
        h = float(parts[4]) if len(parts) > 4 else 1.0
        order = int(parts[5]) if len(parts) > 5 else 2
        pec = len(parts) > 6 and parts[6] == "pec"
        name = case_name(V, a, d, s, h, order) + tag + ("_pec" if pec else "")
        mesh = case_name(V, a, d, s, h) + tag
        args = [sys.executable, str(HERE / "lgr_geometry.py"), "--name", name, "--V", str(V), "--a", str(a),
                "--s", f"{s:.2f}", "--d", str(d), "--h", str(h), "--order", str(order), "--outdir", str(HERE),
                "--mesh-name", mesh]
        if pec:
            args += ["--shield-pec"]
        if (HERE / "mesh" / f"{mesh}.msh").exists():
            args += ["--no-mesh"]
        subprocess.run(args, check=True)
        names.append(name)
    (HERE / "lists").mkdir(exist_ok=True)
    for k in range(nlists):
        chunk = names[k::nlists]
        if chunk:
            (HERE / "lists" / f"{rnd}_{k+1}.txt").write_bytes(("\n".join(chunk) + "\n").encode())


if __name__ == "__main__":
    main()
