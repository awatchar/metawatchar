# S3: resonators loaded with the MNG loop-gap meta-atom array

Eigenmode simulations, made with the open-source FEM solver [Palace](https://github.com/awslabs/palace) 0.18.1, of resonators loaded with the loop-gap (LGR) magnetic meta-atoms of stage S2 ([`../s2_mng_slab`](../s2_mng_slab)). They test three things:

1. whether the loaded resonator is predicted by the effective-medium parameters (n, z, µ_eff) retrieved in S2;
2. the loaded-resonator Q compared with the meta-atom Q;
3. the miniaturization against an empty cavity of the same size, and the Q cost of finer granularity (more, smaller cells in the same volume).

Results: [`results.md`](results.md). Files are listed under [Files](#files).

## Meta-atom (from S2)

Copper loop r = 21.5 mm, wall 3 mm, length 45 mm along y; two 3 mm lips of height s = 6.2 mm; 1.0 mm alumina gap (εr 9.8, tan δ 4e-5); lattice a = 70 mm; slit and lips towards +z ("z" orientation, no first-order magneto-electric coupling). In the S2 unit-cell box (PEC z-plates, PMC y-walls) the meta-atom has

| x-walls | meaning | f (MHz) order 2 / 3 | Q order 2 / 3 |
|---|---|---|---|
| PMC | k = π/a, just below the µ pole f_m | 149.934 / 149.930 | 5,016 / 5,083 |
| PEC | k = 0, µ_eff = 0 (f_mp) | 179.516 / 179.503 | 4,905 / 4,982 |

**These S2 Q values carry half the copper loss** (copper input Rs with Ls = Rs/ω, which Palace combines in parallel; see [Common settings](#common-settings)). With the correct copper impedance this study finds f_m = 149.922 / 149.918 MHz with **Q = 2,811 / 2,847** (Q_c 3,148 / 3,193), and f_mp = 179.504 / 179.491 MHz with Q = 2,746 / 2,789 (order 2 / 3, half model). Frequencies change by < 0.01 %.

## Models

All models mesh only y ≥ 0 and put a PMC symmetry plane at y = 0 (the cell is mirror-symmetric in y). This keeps every mode whose loop flux H_y is even in y, which includes all loop-gap (LC) modes of the arrays. The S2 single-cell eigenmodes are reproduced by the half model to 2e-6 in f and 0.4 % in Q (see results).

### Study 1: short-circuited TEM line (1D effective-medium resonator)

The S2 TEM unit-cell line (PEC z = ±a/2, PMC y = a/2, i.e. an infinite transverse array) holds N = 1–4 cells at x = (i + ½) a, 0 ≤ x ≤ L = N a. Empty line sections of length g sit at both ends, and PEC walls at x = −g and x = L + g short-circuit the line.

- g = a/2 = 35 mm: the end walls sit half a lattice period beyond the slab faces (the S2 reference planes), i.e. at the mirror plane of an empty cell.
- g = 0: the end walls are at the slab faces (for N = 1 this is the S2 PEC box).

Effective-medium prediction (`s3_effective_medium.py`): the slab is replaced by a homogeneous medium of index n(f) and wave impedance z(f), and the transverse-resonance condition is solved at the left short,

Z1 = j tan(k0 g), Z2 = z (Z1 + j z tan(n k0 L)) / (z + j Z1 tan(n k0 L)), Z3 = (Z2 + j tan(k0 g)) / (1 + j Z2 tan(k0 g)), resonance: Z3 = 0,

with impedances normalized to the empty line. Two descriptions of the medium:

- **A (retrieved)**: n(f), z(f) retrieved from the S2 slab S-parameters of N_ret = 1, 2 and 4 cells (Smith/Chen inversion with a Kramers–Kronig branch choice; tables in `results/retrieved_nz_N*.csv`). Roots on the real 10 kHz grid where Im Z3 crosses zero upwards; Q = f (dX/df)/(2R).
- **B (Lorentz)**: µ(f) = µ∞ − F f²/(f² − f0² − j f f0/Q_cell) with constant ε = 2.385 (median retrieved ε_eff), continued to complex f. **B-eig** takes f0 and F from the two meta-atom eigenmodes only (PMC box = the k = π/a point, PEC box = µ_eff = 0, µ∞ = 1) and the meta-atom Q (`B-eig`: 5,035 with the S2 copper input; `B-eig-cu`: 2,811 with correct copper); **B-fit** (`B-fit`, `B-fit-cu`) takes f0, F, µ∞ from a least-squares fit to the retrieved µ_eff (N_ret = 4). Roots outside the first Brillouin zone (Re n k0 a ≥ 0.97 π), which a continuous Lorentz medium has near f0 but the lattice cannot support, are dropped.

### Study 2: compact copper cavity

A closed rectangular box with copper walls (surface impedance, σ = 5.8e7 S/m) around an n_x × n_y × n_z array of the same cells (loop axes along y, all lips towards +z, cells at ((i+½) a, (j+½−n_y/2) a, (k+½) a)). The tight box ("e00") is the lattice boundary itself, i.e. walls a/2 from the cell centres: 4.3 mm from the lip tips, 10.5 mm from the loop sides and 12.5 mm from the loop ends. The roomy box ("e24") is 24.5 mm larger on every side, which puts the walls about a/2 = 35 mm from the cell surfaces. Cases: 1×1×1, 1×1×2 and 2×1×2 (tight), 1×1×1 and 1×1×2 (roomy), and 1×1×2 with lossless PEC walls (to separate wall and loop loss). The empty-box reference is analytic (all modes; copper Q of the lowest mode from Pozar's TE10ℓ formula) and checked with Palace (`E112_e00_cu_o2_zx`, `E112_e00_cu_o2`, `E112_e00_allcu_o2`).

### Study 3: granularity (scaled cells)

The cell is scaled by 0.8 (a = 56 mm; r, wall, lips, length all × 0.8). Pure scaling would raise f by 1.25, so the cell is re-tuned to the same µ pole: the alumina gap is 0.5 mm (a standard substrate thickness) and the lip height is s = 4.91 mm, from two tuning runs at s = 4.0 and 5.5 mm and a linear fit of 1/f² against s. Two pairs of boxes of nearly the same volume are compared: 1×1×2 at a = 70 mm (0.686 L, 2 cells) against 2×1×2 at a = 56 mm (0.702 L, 4 cells), and 2×1×2 at a = 70 mm (1.372 L, 4 cells) against 2×2×2 at a = 56 mm (1.405 L, 8 cells).

### Common settings

- Copper loops and copper walls: impedance boundary with Rs = √(π f µ0/σ) = 3.196 mΩ/sq (150 MHz). **Palace combines the Rs, Ls (and Cs) of an `Impedance` boundary in parallel**, 1/Z = 1/Rs + 1/(jωLs) (Palace boundary guide; `palace/models/surfaceimpedanceoperator.cpp` adds 1/Ls to the stiffness and 1/Rs to the damping). Two inputs are used:
  - `zs = "exact"` (run names ending in `_zx`; the production set): Rs_p = 2 Rs and Ls_p = 2 Rs/ω, which gives the good-conductor impedance Z = Rs (1 + j) at 150 MHz. Checked against the analytic copper Q of the empty 70 × 70 × 140 mm box (TE101-type mode, 2.394 GHz, Rs at 2.4 GHz): Palace 18,531, analytic 18,510.
  - `zs = "s2"` (no suffix): Rs and Ls = Rs/ω, the input used in `../s2_mng_slab`. In the parallel model this is Z = Rs (1 + j)/2, i.e. **half the copper loss**. Empty-box checks: Q 37,062 against 18,510 analytic, and 27,309 against 13,639 for a 70 × 35 × 140 mm all-copper box. This set is kept because the S2 slab S-parameters were computed with the same input, so it is the like-for-like test of the retrieved-medium model A, Q included.
  - Rs is fixed in each run. `postprocess_s3.py` also reports `Q_c_cu` and `Q_cu`, which rescale the copper loss to the true Rs(f) = Rs(150 MHz) √(f/150 MHz) at the mode frequency, given the surface resistance Re Z(f) = k Rs x²/(1 + x²) (x = f/150 MHz; k = 1 for `s2`, 2 for `exact`) that the parallel model applies. This is first-order perturbation; applied to the `s2` runs it reproduces the `_zx` runs to about 1 %. Q_cu is the best estimate of a real copper build.
- Alumina loss: Q_d = 1/(p_d tan δ) with p_d the alumina share of the electric energy; Q_c = (1/Q − 1/Q_d)⁻¹.
- Energy postprocessing per domain: alumina, air inside the loop bores (r, length 45 mm; "pm_bore"), air outside the loops ("pm_air") and, in study 1, the empty end sections ("pm_end"); pm = magnetic energy share, pe = electric energy share.
- Mesh: Gmsh 4.15 (OpenCASCADE), second-order curved tetrahedra, S2 recipe `v2` (0.2 mm on every loop rim and lip edge, 0.565 mm in and around the alumina, 2 mm on copper, 8 mm bulk), about 0.1 M tetrahedra per half cell. Order-2 elements in production; order 3 for the convergence checks.
- Solver: Palace eigenmode (SLEPc Krylov–Schur, shift-and-invert at the target), GMRES with p-multigrid and a sparse direct coarse solve, tol 1e-8. Spurious, critically damped eigenvalues of the impedance eigenproblem (Q < 10) are discarded.

## Files

| File | Purpose |
|---|---|
| `s3_geometry.py` | Builds the geometry and mesh (line, cavity, empty box) from the S2 cell (`../s2_mng_slab/s2_geometry.py`, imported unchanged) and writes the Palace configs |
| `run_s3.sbatch` | Slurm job: runs Palace for each run in a list file, then `postprocess_s3.py` |
| `postprocess_s3.py` | Inside the job: eigenmode table per run (f, Q, Q_c, Q_d, energy shares). Standard library only |
| `s3_effective_medium.py` | Effective-medium (transverse-resonance) predictions for study 1, models A and B |
| `collect_s3.py` | Builds the deliverable tables and the figure from the per-run outputs |
| **`results/study1_line.csv`** | Study 1: every full-wave mode 120–200 MHz with f, Q, energy shares and the matched predictions (f, Δf %, Q) of each model |
| **`results/em_predictions.csv`**, `em_parameters.json` | Effective-medium predictions (all models, all N, g) and the model parameters |
| **`results/study2_cavity.csv`** | Studies 2 and 3: loaded-cavity modes with f, Q, energy shares, the empty-box lowest mode, miniaturization factor and size in wavelengths |
| `results/eigen_all.csv`, `results/<run>_eig.csv` | Every physical mode of every run |
| `results/s3_summary.json` | All of the above in one file, plus the empty-box check |
| `results/retrieved_nz_N*.csv` | n, z, µ_eff, ε_eff retrieved from the S2 slabs (physics convention, exp(−iωt)) |
| `results/study1_line.png` | Figure: study 1 frequencies and Q against N |
| `results/*-timings.txt`, `*-log-excerpt.txt` | Job timings and solver summaries |
| `cases/`, `configs/`, `runs/`, `lists/` | Mesh/geometry metadata, Palace configs, run metadata, job lists |

Run names (suffix `_zx` = correct copper impedance; no suffix = S2 copper input): `L<N>_g<g>_o<order>` (study 1, PEC ends), `L1_g00_<pec|pmc>_o2` (S2 box check), `C<nx><ny><nz>_e<ext>_<cu|pec>_o<order>` (cavities), `E112_e00_cu_o2` (empty box), `K080_*` (cells scaled by 0.8; `K080_s<s·100>_<pmc|pec>_o2` = scaled meta-atom).

## How to run

Meshes are not stored (they are large); regenerate them locally with Python 3 and `gmsh` ≥ 4.11, e.g.

```bash
python s3_geometry.py mesh --kind line --N 4 --g 35                 # -> mesh/L4_g35_s0620_v2.msh
python s3_geometry.py config --mesh L4_g35_s0620_v2 --name L4_g35_o2_zx --target 0.12 --nmodes 7   # --zs exact (default)
python s3_geometry.py mesh --kind cavity --nx 2 --ny 1 --nz 2       # -> mesh/C212_e00_s0620_v2.msh
python s3_geometry.py config --mesh C212_e00_s0620_v2 --name C212_e00_cu_o2_zx --target 0.13 --nmodes 6
python s3_geometry.py mesh --kind cavity --nx 2 --ny 2 --nz 2 --scale 0.8 --d 0.5 --s 4.91
```

Copy the folder to the cluster and submit a list of run names:

```bash
sbatch -A <ACCOUNT> -p <PARTITION> --ntasks=64 -t 06:00:00 \
       --export=ALL,PALACE_ENV=<PALACE_INSTALL>/env.sh run_s3.sbatch lists/<round>.txt
```

Copy back `results/` and `postpro/*/{eig.csv,domain-E.csv,palace.json}`, then run `python s3_effective_medium.py --retrieval-module <dir of mu_eff_retrieval.py>` and `python collect_s3.py`.
