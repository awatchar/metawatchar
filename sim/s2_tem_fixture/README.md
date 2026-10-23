# S2 fixture: full-wave model of the TEM-line μ_eff measurement fixture

Frequency-domain simulations, made with the open-source FEM solver [Palace](https://github.com/awslabs/palace) 0.18.1, of the fixture planned for measuring the effective permeability μ_eff of a slab of loop-gap (LGR) magnetic meta-atoms at VHF. The fixture is a 50-Ω air microstrip line ("open") or the same line closed by side walls and a lid into a rectangular-coax TEM cell ("closed"). The slab is 4 cells across the line by N_x = 1 or 2 cells along it (also 5 across as a check). The cells are the unshielded LGR cells of [`../s2_mng_slab`](../s2_mng_slab), the ideal unit-cell model (infinite 2D array) this fixture is meant to reproduce.

The question: after IEEE P370 2x-thru de-embedding and S-parameter retrieval, does the fixture return the ideal unit-cell μ_eff, and which variant (open or closed) should be built?

Results and the answer: [`results.md`](results.md).

## Model

| Item | Definition |
|---|---|
| Line | Ground plane z = 0. Top strip 2 mm thick, underside at h = 70 mm (= the S2 lattice period), width w = 344 mm (w/h = 4.91, 50 Ω), uniform for \|x\| ≤ 210 mm |
| Tapers | 210 ≤ \|x\| ≤ 510 mm: underside height 70 → 3 mm and width 344 → 14.7 mm, both linear in x, so w/h stays 4.91 (constant-impedance taper; a planar trapezoid) |
| Launch | A 50-Ω lumped port across the 3-mm gap at each taper tip (x = ±510 mm), standing in for the N-type jack; removed with the taper by the 2x-thru de-embedding |
| Supports | 4 PTFE posts (Ø12 mm × 70 mm, εr 2.1, tan δ 2e-4) under the strip at x = ±195 mm, y = ±(w/2 − 10 mm) |
| Open variant | Air box \|x\| ≤ 2510 mm, y ≤ 2500 mm, z ≤ 2500 mm with first-order absorbing faces (runs `openfar_*`). The ground plane is the box floor (infinite ground). A first set with a small box (\|x\| ≤ 700, y, z ≤ 450 mm, runs `open_*`) is kept as a documented failure: a first-order absorbing boundary that close to the line absorbs its near field (see results) |
| Closed variant | The same line plus PEC-like aluminium side walls at y = ±252 mm (the ground-plate edges), a lid at z = 300 mm and end walls at x = ±520 mm (the jack flanges). `closed_*`: strip as above (w = 344 mm); `closedw280_*`: strip narrowed to w = 280 mm (= 4 cell periods) in the uniform section, which brings the closed line back to 50 Ω; its taper runs linearly from 280 × 70 mm to 14.7 × 3 mm |
| Symmetry | Half model, y ≥ 0. The plane y = 0 is a magnetic wall for the fundamental (even) mode. The half-model lumped port has R = 100 Ω, which gives the S-parameters of the full 50-Ω structure |
| Cells | S2 LGR cell: copper tube r = 21.5 mm, wall 3 mm, length 45 mm along y (the line's H-field direction), axial slit 1.0 mm with alumina (εr 9.8, tan δ 4e-5), lips s = 6.2 mm towards +z (towards the strip). Loop axis at mid-height z = 35 mm, lip tips 4.3 mm below the strip, as in the S2 unit cell. Period 70 mm in x and y. 4 across: y = ±35, ±105 mm; 5 across: y = 0, ±70, ±140 mm. N_x = 1: x = 0; N_x = 2: x = ±35 mm |
| Copper (cells) | Impedance boundary, σ = 5.8e7 S/m, Rs at 150 MHz. Palace's `Impedance` boundary is a parallel Rs ‖ Ls model, so the input is Rs_p = 2 Rs, Ls_p = 2 Rs/ω (the corrected input of `../s2_mng_slab`, runs `_cu2`) |
| Fixture metal | Aluminium 6061-T6, σ = 2.5e7 S/m (about 43 % IACS), same parallel input at 150 MHz, on the strip, the ground and the closed-variant walls. `_pec` runs use PEC instead, to separate radiation from conductor loss |
| Conventions | exp(+jωt) in Palace and in every CSV; μ_eff is reported in the physics convention (μ = μ′ + iμ″, μ″ > 0 = loss), as in `../s2_mng_slab` |

### Mesh

Gmsh 4.15 (OpenCASCADE), second-order curved tetrahedra (`fixture_geometry.py`):

- strip faces: h_loc/10 (7 mm in the uniform section, 0.6 mm minimum at the tips), strip edges h_loc/20, growing by 0.25 per unit distance;
- cells: exactly the S2 "v2" recipe (0.2 mm on every loop rim and lip edge, 0.565 mm in and around the alumina, 2 mm on copper, 8 mm out to 15 mm), graded to 40 mm at 120 mm;
- bulk: 40 mm (small boxes) or growing to 150 mm (large open box, λ/10 at 200 MHz at order 2).

About 0.3 M tetrahedra for the empty fixture and 0.19 M per meta-atom (0.65–1.1 M with cells).

## Solver

- Palace driven solver, order 2, adaptive fast sweep (PROM, `AdaptiveTol` 1e-4), GMRES with p-multigrid and a sparse direct coarse solve, tol 1e-8.
- 2x-thru (empty fixture): uniform grid f = k·df up to 300 MHz (df = 0.5 MHz for the small boxes, 2.5 MHz for the large open box, where low-frequency solves are slow); IEEE P370 needs the grid to start at df.
- Slab runs: 130–200 MHz in 50 kHz steps plus 146–168 MHz in 10 kHz steps (3161 points), port 1 excited; S22 = S11 and S12 = S21 by the x-mirror symmetry of the structure (checked on the empty fixture: \|S11 − S22\| ≤ 7e-5).
- Postprocessing in Palace: port S-parameters, field probes (B and E) at the centre of every loop bore and along the line cross-section, and the power flux through the absorbing boundary (radiation).

### De-embedding and retrieval (`deembed_retrieve.py`)

1. The empty fixture is the 2x-thru. All networks are first renormalized from 50 Ω to the uniform line's impedance Z_line (from `quasistatic_xsec.py`; the band-limited TDR of the 2x-thru agrees to 0.2 % for the open line), then split with scikit-rf's `IEEEP370_SE_NZC_2xThru`. A synthetic self-test (taper + line + known slab built from the S2 μ_eff) shows that splitting at 50 Ω leaves an error of order (Z_line − 50)/100 in S, while renormalizing first does not (`results/deembed_selftest.json`).
2. The error boxes are interpolated (magnitude, unwrapped phase) onto the slab grid and removed. The result has both reference planes at the line centre x = 0; the planes are then moved outward by L/2 of air line each (β = k0): S_slab = D·exp(−j k0 L).
3. μ_eff retrieval: Smith/Chen inversion with a Kramers–Kronig branch choice (`mu_eff_retrieval.py`, passed with `--retrieval`). The z-sign flip that the module applies where Im n < 0 is undone where Im n is only slightly negative (marginally active de-embedded data), which otherwise flips the sign of μ and ε.

### Cross-section model (`quasistatic_xsec.py`)

The empty fixture is all air, so its mode is TEM and its transverse fields are the 2D static fields of the cross-section. The script solves 2D electro- and magnetostatics on a tensor grid for Z0, the transverse field profile, the TE/TM cut-off frequencies of the closed cross-section, and the "series reluctance" share p of a homogenized slab (1/μ_fix = p/μ + 1 − p, a demagnetization factor N_d = 1 − p). It runs in seconds.

## Files

| File | Purpose |
|---|---|
| `fixture_geometry.py` | Geometry, mesh and Palace configs (driven fixture, 2D cross-section). Imports the cell from `../s2_mng_slab/s2_geometry.py` |
| `run_fixture.sbatch` | Slurm job: runs Palace for each run in a list file (optional per-run rank count and time cap), then `postprocess_fixture.py`. Account, partition and install path are given at submit time |
| `postprocess_fixture.py` | Inside the job: port-S (dB/deg) to re/im, copies of the probe, flux and energy CSVs, a one-line summary per run. Standard library only |
| `quasistatic_xsec.py` | Cross-section statics: Z0, field uniformity, demagnetization share p, closed-cell cut-offs |
| `deembed_retrieve.py` | P370 2x-thru de-embedding, reference-plane shift, μ_eff retrieval, synthetic self-test |
| `analyze_fixture.py` | Figures and summary JSON files |
| `results/<run>_S.csv` | S-parameters at the fixture ports (taper tips, 50 Ω): `f_Hz, S11_re, S11_im, S21_re, S21_im` (+ `S12`, `S22` when both ports were excited) |
| `results/<run>_slab.csv`, `<run>_slab.mu_eff.csv` | De-embedded slab S-parameters (reference planes at the slab faces, reference Z_line) and the retrieved μ, ε, n, z |
| `results/<run>_probe-B.csv`, `_probe-E.csv`, `_surface-F.csv` | Palace probe fields (1 W incident) and radiated power through the absorbing boundary |
| `results/*_summary.json`, `quasistatic_*.json`, `deembed_selftest.json` | Summaries used in `results.md` |
| `figures/*.png` | Figures in `results.md` |
| `cases/`, `configs/`, `runs/`, `lists/` | Geometry records, Palace configs, run records, job lists |

Meshes (`mesh/`), solver logs and the raw Palace output directories are not kept (they are regenerated by `fixture_geometry.py` and the job script).

## Reproduce

```bash
python fixture_geometry.py mesh --variant open --nx 0 --bulk 150 --box 2510 2500 2500 --name openfar_nx0
python fixture_geometry.py config --mesh openfar_nx0 --name openfar_nx0_al_o2b --order 2 --metal al --band thru25
# ... as listed in lists/*.txt; then on the cluster:
sbatch -A <ACCOUNT> -p <PARTITION> --ntasks=128 -t 01:55:00 \
       --export=ALL,PALACE_ENV=<PALACE_INSTALL>/env.sh run_fixture.sbatch lists/p9.txt
# locally, on the pulled results:
python deembed_retrieve.py run --retrieval <dir of mu_eff_retrieval.py> --thru openfar_nx0_al_o2b --dut openfar_nx1_al_o2 --nx 1 --zline 49.55
python analyze_fixture.py slab --retrieval <dir of mu_eff_retrieval.py>
```
