# S2: MNG slab of unshielded loop-gap cells in a TEM unit-cell waveguide

Frequency-domain simulations, made with the open-source FEM solver [Palace](https://github.com/awslabs/palace) 0.18.1, of a slab of N = 1, 2 and 4 loop-gap (LGR) magnetic meta-atoms along the propagation direction. The slab sits in an ideal TEM unit-cell waveguide, so the model describes an infinite 2D array in the transverse plane.

The outputs are:

- de-embedded complex S11 and S21 for retrieving the effective permeability mu_eff(f);
- a check that mu_eff does not depend on the slab thickness (homogenization);
- the single-cell eigenmodes that bracket the MNG band.

Results: [`results.md`](results.md). Deliverable files are listed under [Files](#files).

## Model

| Item | Definition |
|---|---|
| Propagation | +x. The slab occupies 0 ≤ x ≤ L = N a_x, and cell i is centred at x = (i + ½) a_x, y = z = 0 |
| Lattice | a_x = a_y = a_z = 70 mm (about λ/29 at 150 MHz) |
| Transverse walls | PEC at z = ±a_z/2 (TEM E along z) and PMC at y = ±a_y/2 (TEM H along y). The mirror images make an infinite 2D array in y–z (neighbouring cells are mirror images) |
| Ports | Numeric wave ports (TEM mode) at x = −140 mm and x = L + 140 mm, so each side has 2 a_x of empty line. Higher-order modes (cutoff about 2.1 GHz) decay by more than 55 dB before reaching a port. Palace's port `Offset` (140 mm) de-embeds the reference planes to the slab faces x = 0 and x = L |
| Conventions | Time dependence exp(+jωt): an empty line of length l has S21 = exp(−j k0 l). Ports are power-normalized, so the reference impedance is the TEM line impedance η0 a_z/a_y = η0 |
| Loop | Copper tube with its axis along **y** (parallel to the TEM H-field, so H threads the loop): inner radius r = 21.5 mm, wall t = 3 mm, length 45 mm (12.5 mm from each PMC wall), centred in the cell |
| Capacitor | One axial slit of width d = 1.0 mm through the wall. Two copper lips, 3 mm thick and full loop length, extend radially outward by s = 6.2 mm. Alumina (εr 9.8, tan δ 4e-5) fills the slit and the space between the lips. This is the LGR-A loop of `../lgr_size_sweep` (0.5 L, a = 0.5, d = 1.0 mm) without its can |
| Orientation (production, `z`) | Slit and lips point to **+z**, towards the top plate. The lip plates are y–z planes, **parallel to E**, and the gap field is along x, perpendicular to E. The cell is mirror-symmetric in x, so the uniform incident E_z does not drive the LC mode, there is no first-order electric (magneto-electric) coupling, and S22 = S11. Lip tips end 4.3 mm from the plate |
| Orientation (comparison, `x`) | The same cell rotated by +90° about y: slit and lips point to +x (downstream), and the gap field is along z, parallel to E. E drives the capacitor directly (bianisotropic response, S11 ≠ S22). Run for N = 1 only |
| Copper | Not meshed. Copper surfaces are impedance boundaries with Rs = √(π f µ0/σ) = 3.196 mΩ/sq (σ = 5.8e7 S/m) evaluated at 150 MHz, plus Ls = Rs/ω. Rs is held constant over the sweep (it varies as √f, ±8 % over 130–200 MHz) |
| Plates and walls | Ideal PEC/PMC (lossless) |

Tuning: s was set from two order-2 eigenmode runs (s = 4.0 and 7.5 mm) with a linear fit of 1/f² against s. The fit is 1/f² ∝ (s + 3.33 mm). The target put the PMC-x-wall eigenmode, which approximates the µ pole, at 150 MHz, giving s = 6.2 mm.

### Single-cell eigenmodes (reference)

The same PEC/PMC box is closed at x = 0 and x = a_x by either:

- **PMC x-walls:** the images alternate along x, which is the k = π/a_x lattice mode. In an effective medium this needs µ_eff ≈ (π c/(ω a_x))² ≈ 200, so it sits just below the **pole f_m** of µ_eff.
- **PEC x-walls:** the images are in phase (k = 0) and the macroscopic E is zero. A uniform H then needs B = µ_eff H = 0, so this mode is the **zero µ_eff = 0** (magnetic plasma frequency f_mp, the upper edge of the MNG band).

The ratio gives the Lorentz filling factor F = 1 − (f_m/f_mp)².

## Solver

- Palace driven solver with the adaptive fast frequency sweep (PROM, `AdaptiveTol` = 1e-4). The output grid is 130–200 MHz in 10 kHz steps (7001 points).
- The adaptive sweep needed only 6–9 full solves per case. It was validated against 7 uniform (full) solves across the N = 1 resonance (149–152 MHz): |ΔS| ≤ 3e-5.
- The empty-line reference reproduces S21 = exp(−j k0 · 280 mm) without de-embedding, and S21 = 1 (phase error < 1e-6°) with |S11| < −178 dB with it.
- Linear solver: GMRES with p-multigrid and a sparse direct coarse solve, tol 1e-8.

### Mesh

Gmsh 4.15 (OpenCASCADE), second-order curved tetrahedra (`s2_geometry.py`), recipe `v2`. This is the LGR-A production recipe:

- 0.565 mm in and around the alumina;
- 0.25 mm along its edges;
- 2 mm on the copper surfaces;
- **0.2 mm along every loop rim and lip edge** (the magnetic field is singular there, and those edges control the conductor loss);
- 8 mm in the bulk.

That gives about 0.19 M tetrahedra per cell. Recipe `v3` refines every one of these sizes by about 0.6–0.7× (0.12 mm on the metal edges).

## Files

| File | Purpose |
|---|---|
| `s2_geometry.py` | Builds the geometry and mesh (driven slab, single-cell eigen box, empty line) and writes the Palace configs |
| `run_s2.sbatch` | Slurm job: runs Palace for each run in a list file, then `postprocess_s2.py` |
| `postprocess_s2.py` | Inside the job: converts `port-S.csv` (dB/deg) to re/im, summarizes the S21 dip and absorption, and summarizes the eigenmodes (Q, Q_c, Q_d, p_d). Standard library only |
| `collect_s2.py` | Gathers the deliverables below from the per-run outputs |
| **`results/S_N1.csv`, `S_N2.csv`, `S_N4.csv`** | **Production S-parameters**: `f_Hz, S11_re, S11_im, S21_re, S21_im`, de-embedded to the slab faces |
| **`results/s2_slab.json`** | Lattice, cell dimensions, slab length L per N, reference-plane and port positions, materials, conventions, eigenmodes, run list |
| `results/S_N1_orient_x.csv` | N = 1 with the opposite capacitor orientation (also `S12`, `S22`) |
| `results/S_N1_o2.csv`, `S_N2_o2.csv`, `S_N4_o2.csv` | The same slabs at order 2 (a consistent lower level). N = 1 convergence levels: `N1_z_o2_S.csv`, `N1_z_o3_S.csv`, `N1_z_v3_o2_S.csv` |
| `results/<run>_S.csv`, `<run>_summary.csv`, `<run>_eig.csv` | Per-run outputs; `driven_summary.csv` and `eigen_summary.csv` merge them |
| `results/*-timings.txt`, `*-log-excerpt.txt` | Job timings and solver summaries |
| `postpro/<run>/` | Palace's own small outputs (`port-S.csv`, `eig.csv`, `palace.json`) |
| `cases/`, `configs/`, `runs/`, `lists/` | Mesh/geometry metadata, Palace configs, run metadata, job lists |

Run names: `N<N>_<orient>_o<order>` (driven, v2 mesh), `N1_z_v3_o<order>` (refined mesh), `eig_z_s<s·100>_<pmc|pec>_o<order>` (eigenmode), `empty_off<offset>_o2` (empty line).

## How to run

Meshes are not stored (they are large); regenerate them locally with Python 3 and `gmsh` ≥ 4.11:

```bash
python s2_geometry.py mesh --kind driven --N 4 --orient z --s 6.2            # -> mesh/N4_z_s0620_v2.msh
python s2_geometry.py config --mesh N4_z_s0620_v2 --name N4_z_o3 --type driven --order 3
python s2_geometry.py mesh --kind eigen --orient z --s 6.2
python s2_geometry.py config --mesh eig_z_s0620_v2 --name eig_z_s0620_pmc_o3 --type eigen --xwall pmc --order 3 --target 0.13 --nmodes 2
```

Copy the folder to the cluster and submit a list of run names:

```bash
sbatch -A <ACCOUNT> -p <PARTITION> --ntasks=64 -t 06:00:00 \
       --export=ALL,PALACE_ENV=<PALACE_INSTALL>/env.sh run_s2.sbatch lists/<round>.txt
```

Copy back `results/` and `postpro/*/port-S.csv`, `eig.csv`, `*.json`, then run `python collect_s2.py`.
