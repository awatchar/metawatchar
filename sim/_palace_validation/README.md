# Palace validation: cylindrical cavity eigenmodes

Checks a source build of the open-source FEM solver [Palace](https://github.com/awslabs/palace)
0.18.1 (Apache-2.0) on the LANTA HPC system against closed-form results for a cylindrical cavity
resonator: resonant frequencies and quality factors (Q).

## Problem

Palace's `examples/cylinder` cavity: a circular cylinder of radius a = 2.74 cm and height
d = 2a, filled with Teflon (eps_r = 2.08, tan_delta = 4e-4). The lowest modes above 2 GHz are
computed with 4th-order Nédélec elements. The eigensolver is SLEPc (shift-and-invert at 2 GHz),
and the linear solver is Palace's default: GMRES with p-multigrid (p = 1, 2, 4) and a
SuperLU_DIST direct solve on the coarse level.

| Config | Walls | Mesh | ND unknowns | Expected Q |
|---|---|---|---|---|
| `cavity_pec.json` | PEC | 80 curved hexahedra (`cylinder_hex.msh`) | 16,544 | Q = Q_d = 1/tan_delta = 2500 |
| `cavity_impedance.json` | copper, Rs = 0.0184 ohm/sq (5 GHz) | 96 curved prisms (`cylinder_prism.msh`) | 11,224 | 1/Q = 1/Q_d + 1/Q_c |

The two JSON files are Palace's example configurations, unchanged (Apache-2.0, © Amazon.com,
Inc. or its affiliates). The job copies the meshes from the Palace source tree at run time.

## Analytic reference (`compare_cavity.py`)

- Frequencies of the TE_nml / TM_nml modes from the Bessel-function roots (Pozar, *Microwave
  Engineering*, Sec. 6.4), with c0 = 299 792 458 m/s.
- Q_d = 1/tan_delta.
- Q_c = omega W / P_c, computed from the closed-form mode fields for every mode (TE and TM):
  W = (eps/2) ∫|E|² dV and P_c = (Rs/2) ∮|H_tan|² dS over the side wall and both end caps,
  with the radial integrals evaluated numerically. This reproduces Pozar's Example 6.4
  (TE_011, Q_c ≈ 2.94e4) and the textbook TM_010 result 1.2025 eta / (Rs (1 + a/d)).
- Simulated Q_c is recovered as (1/Q − 1/Q_d)^−1.
- Degenerate modes (cos/sin pairs, and TE_011 with TM_111 at 5.000 GHz) are matched by
  frequency order, then by the Q that fits best inside each degenerate group.

## How to run (LANTA)

Copy this folder to `<PROJECT_DIR>/sim/_palace_validation` on LANTA, then:

```bash
cd <PROJECT_DIR>/sim/_palace_validation
sbatch -A <ACCOUNT> --export=ALL,PALACE_ENV=<PROJECT_DIR>/software/palace/env.sh run_cavity.sbatch
```

`env.sh` loads the module set Palace was built with and puts `palace-x86_64.bin` on `PATH`.
The job runs both configs with 16 MPI ranks on one node, then writes `results.md` (the
comparison tables) and `timings.txt`. Copy back only `postpro/*/{eig,domain-E,error-indicators}.csv`,
`postpro/*/*.json`, `results.md`, and `timings.txt`; skip `postpro/*/paraview` (field dumps,
~35 MB), `mesh/`, `logs/`, and `slurm-*.out`.

## Results (22 Oct 2026)

Run: 1 node, 16 MPI ranks (AMD EPYC 7543), Palace 0.18.1. Wall time per config ≈ 10 s
(Palace total 8.7 s each); peak memory ≈ 1.2 GB in total.

| Case | Modes compared | Max \|f error\| | Q check | Max \|Q error\| |
|---|---|---|---|---|
| PEC walls (hex mesh) | 18 (TM010 … TE311) | 4.0e-4 | Q vs Q_d = 2500 | 7e-8 |
| Copper walls, Rs = 0.0184 ohm/sq (prism mesh) | 15 (TM010 … TE212) | 7.9e-5 | Q_c vs analytic Q_c | 2.5e-5 |

Selected modes, copper walls:

| Mode | f analytic (GHz) | f Palace (GHz) | Q analytic | Q Palace | Q_c analytic | Q_c Palace |
|---|---|---|---|---|---|---|
| TM010 | 2.903636 | 2.903862 | 2049.7 | 2049.7 | 1.138e4 | 1.138e4 |
| TE111 | 2.922197 | 2.922328 | 2072.3 | 2072.3 | 1.211e4 | 1.211e4 |
| TM011 | 3.468175 | 3.468364 | 2007.7 | 2007.7 | 1.019e4 | 1.019e4 |
| TE011 | 5.000146 | 5.000483 | 2304.0 | 2304.0 | 2.940e4 | 2.939e4 |
| TM111 | 5.000146 | 5.000486 | 2136.6 | 2136.6 | 1.470e4 | 1.470e4 |

All 15 lossy-wall modes agree with the analytic Q_c to within 2.5e-5 and in frequency to within
8e-5. The remaining frequency error is discretization error of these small example meshes
(hex: ≤ 4e-4, prism: ≤ 8e-5). Palace's documentation shows it decreasing under mesh/order
refinement. The eigenvalues also match Palace's own regression reference outputs for these
two configs (`test/data/regression/ref/cylinder/`): relative differences ≤ 6e-10 in f and
≤ 1e-5 in Q, so this build reproduces the upstream results. Full per-mode tables are in
[`results.md`](results.md), and the raw Palace outputs are in `postpro/`.
