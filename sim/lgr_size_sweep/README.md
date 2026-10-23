# LGR-A size sweep: shielded loop-gap resonator at 145.5 MHz

Eigenmode simulations (open-source FEM solver [Palace](https://github.com/awslabs/palace) 0.18.1) of
a single shielded loop-gap-resonator (LGR) unit cell, swept over shield volume and loop radius, to
choose the size of a first hardware prototype for a 2 m band (145.2/145.8 MHz) duplexer resonator.
Reported per case: resonant frequency, conductor Q, dielectric Q, total Q, alumina energy
participation, and peak electric fields normalized to a stored energy of 10 uJ.
Results: [`results.md`](results.md).

## Geometry "LGR-A"

All parts are coaxial with the z axis (`lgr_geometry.py`):

| Part | Definition |
|---|---|
| Shield | closed copper can, inner radius R, inner height H = 2R, volume V = 2 pi R^3; V = 0.15, 0.25, 0.50 L gives R = 28.79, 34.14, 43.01 mm |
| Loop | copper tube, inner radius r = a R (a = 0.4, 0.5, 0.6), wall t = 3 mm, length z = H - 2c, end clearance c = 0.25 R (so z = 1.5 R), centred |
| Slit | one axial slit through the loop wall, width d = 0.635 mm (standard alumina substrate thickness), full length z |
| Lips | two flat copper plates, 3 mm thick, extending radially outward from the loop's outer surface by s, full length z, separated by d |
| Gap dielectric | alumina, eps_r = 9.8, tan_delta = 4e-5, filling the d-wide gap through the wall and between the lips (radial extent t + s, length z) |

s is tuned so the loop-gap mode sits at 145.5 +- 1 MHz, subject to r + t + s <= R - 3 mm.
Mechanical supports are ignored. A variant with d = 1.0 mm (thicker alumina substrate, easier to
machine) is also run for the 0.5 L cases.

## Model

- **Conductors** (shield, loop, lips) are not meshed. They are surface-impedance boundaries with
  Rs = sqrt(pi f mu0 / sigma) = 3.147 mOhm/sq for copper (sigma = 5.8e7 S/m) at 145.5 MHz, plus the
  matching surface inductance Ls = Rs / omega.
  **Correction:** Palace combines an `Impedance` boundary's Rs and Ls in parallel, so this input
  (Rs, Ls = Rs/omega) is Z = Rs (1 + j)/2, half the copper loss, and the sweep's Q values are too
  high. The correct input is Rs_p = 2 Rs, Ls_p = 2 Rs/omega (`lgr_geometry.py --zs exact`, now the
  default; `--zs s2` reproduces the sweep). A corrected subset (`*_cu2`) is in `results.md`.
- **Dielectric loss** comes from the alumina loss tangent. Q_d = 1/(p_d tan_delta), where p_d is the
  share of electric energy in the alumina (Palace `domain-E.csv`). The conductor-only Q follows from
  the eigenvalue's total Q: Q_c = (1/Q - 1/Q_d)^-1.
- **Fields** are sampled with Palace probes: in the gap (wall slit and between the lips), and in the
  air or alumina 0.1 mm from the copper/alumina/air triple lines at the lip tips, the bore side and the
  slab ends. The eigenvector is rescaled so that the stored energy W (Palace's E_elec = 1/2 Re
  int D*.E dV, which equals the total stored energy of the mode) is 10 uJ. Fields at sharp metal edges
  are singular, so edge values depend on the probe offset and, in hardware, on the edge radius.
- **Mesh**: Gmsh 4.15 (OpenCASCADE), second-order curved tetrahedra (`lgr_geometry.py`): 0.45 mm in
  and around the alumina slab, 0.25 mm along its edges, 2 mm on the loop and lip surfaces, 0.2 mm
  along every loop and lip edge, and R/5 in the bulk ("v2" recipe, about 0.22 M tetrahedra).
  The magnetic field is singular at the square tube rims and lip corners, and resolving those
  edges is what sets the conductor-Q accuracy: a first recipe without the edge refinement ("v1")
  gave Q_c about 7–8 % low in every case. Frequency tuning used order-2 Nedelec elements on v1
  meshes; reported cases use order 3 on v2 meshes (about 4.5 M unknowns). The convergence study
  (order 2–4, edge size 0.4/0.2/0.1 mm, surface, bulk and gap refinement, adaptive refinement) is
  in `results.md`.
- **Analytic reference** (`analyze.py`): coaxial-return LGR model without end caps,
  Q_c = omega mu0 pi r^2 (1 + rho) / (2 pi Rs [r + rho^2 ((r + t) + R)]), rho = r^2/(R^2 - (r + t)^2).
  The Palace/analytic ratio calibrates that model for this closed-can geometry.

## Files

| File | Purpose |
|---|---|
| `lgr_geometry.py` | Builds the geometry and mesh for one case (Gmsh Python API) and writes its Palace config and metadata |
| `palace_template.json` | Palace eigenmode configuration template (materials, impedance boundaries, energy and probe postprocessing, solver) |
| `make_round.py` | Builds a batch of cases (`V,a,d,s[,h[,order[,pec]]]`) and the case-list files for the job |
| `run_sweep.sbatch` | Slurm job: runs Palace for every case in a list file, then `postprocess.py` |
| `postprocess.py` | Summarizes each case's `eig.csv`, `domain-E.csv` and `probe-E.csv` into one CSV row (numpy-free) |
| `analyze.py` | Local analysis: lip-height tuning fit, analytic Q_c, and the `results.md` tables |
| `results/*.csv` | Summarized results of every run (one row per case) and job timings |
| `postpro/<case>/` | Palace's small CSV outputs per case (field dumps are disabled) |
| `cases/`, `configs/`, `lists/` | Case metadata, generated Palace configs, and job case lists |

Case names are `v<V*100>_a<a*10>_d<d*1000>_s<s*100>[_h<scale*100>|_m<variant>][_o<order>][_v2][_pec]`
(s in units of 0.01 mm). `_v2` marks the production mesh recipe, `_h`/`_m` mark mesh-convergence
variants, and `_pec` is a diagnostic run with a lossless shield (loop and lip losses only).
Production runs are the `*_o3_v2` cases.

## How to run

Meshes are not stored (they are large); regenerate them locally:

```bash
python make_round.py --tag _v2 <round> <n_lists> 0.25,0.5,0.635,10.92,1.0,3 ...
```

This needs Python 3 with `gmsh` (>= 4.11, OpenCASCADE). Copy the folder to the cluster, then
submit one job per list file:

```bash
sbatch -A <ACCOUNT> --export=ALL,PALACE_ENV=<PALACE_INSTALL>/env.sh run_sweep.sbatch lists/<round>_1.txt
```

`env.sh` loads the module set Palace was built with and puts `palace-x86_64.bin` on `PATH`.
Each job uses 32 MPI ranks on one node. On a v1 mesh one case takes about 2 min at order 2
(~0.6 M unknowns) and 8–11 min at order 3 (~1.7 M unknowns, ~15 GB). On a v2 mesh at order 3
(~4.5 M unknowns, ~40 GB) it takes 35–90 min, depending on node load. Copy back `results/`,
`postpro/*/*.csv` and `postpro/*/*.json`, then run `python analyze.py report results/*.csv > results.md`.
