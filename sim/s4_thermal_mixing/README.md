# S4: thermally self-compensating MNG medium by cell mixing

Eigenmode simulations, made with the open-source FEM solver [Palace](https://github.com/awslabs/palace) 0.18.1, of small MNG media built from the loop-gap (LGR) meta-atoms of stages S2/S3 ([`../s2_mng_slab`](../s2_mng_slab), [`../s3_loaded_resonator`](../s3_loaded_resonator)). Two kinds of meta-atom with opposite temperature coefficients of resonant frequency (τf) are mixed inside one 4-cell medium:

- **A cells:** the S2/S3 meta-atom, with an alumina-filled gap (τf ≈ −72 ppm/K);
- **TiO₂ cells:** the same copper loop with a rutile TiO₂ gap (τf ≈ +420 ppm/K).

Questions:

1. How does τf of the working mode of a 4-cell resonator change with the TiO₂ share?
2. What is its Q?

Results: [`results.md`](results.md). Files are listed under [Files](#files).

## Meta-atoms

All cells share the S2/S3 copper loop: inner radius 21.5 mm, wall 3 mm, length 45 mm along y, one axial slit facing +z, lattice a = 70 mm. The copper is built by `build_cell()` of [`../s2_mng_slab/s2_geometry.py`](../s2_mng_slab/s2_geometry.py), imported unchanged. Only the capacitor differs:

| Spec | Gap | Use |
|---|---|---|
| `A` | 1.0 mm slit + two 3 mm × 6.2 mm lips, alumina in the whole slit and lip gap (S2/S3 cell) | alumina cell |
| `W<d>` | slit of width d, **no lips**, TiO₂ filling the whole slit (3 mm × 45 mm × d) | **TiO₂ cell used in the mixing studies**; d = 3.29 mm matches the alumina cell's resonance |
| `T<l>` | A-cell copper (1.0 mm slit + lips), a TiO₂ block of length l (along y, centred) in the gap, air in the rest | first TiO₂-cell design; dropped (see results: current crowding costs 14 % of Q_c) |
| `C<l>` | A cell with a centred TiO₂ block of length l replacing alumina | composite gap: meta-atom-level compensation, reference only |

Rutile has εr ≈ 10× that of alumina, so a TiO₂ cell with the alumina cell's copper would need a 10 mm gap. Two ways to keep the resonance were tried: a short TiO₂ block in the existing gap (`T`), and a full-length TiO₂ slab in a wider slit without lips (`W`).

## Resonators and mode naming

- **Unit-cell box** (`P1_*`): one cell in the S2 TEM unit cell (PEC z-plates, PMC y-walls), closed at x = 0, a by PMC walls (k = π/a, just below the µ pole f_m) or PEC walls (µ_eff = 0, f_mp). These two modes give each meta-atom's Lorentz parameters and its own τf.
- **Tight copper box** (`B212_*`, the resonator): a 2 × 1 × 2 array (x × y × z) of cells in a copper box at the lattice boundary, 140 × 70 × 140 mm = 1.372 L (S3 study 2 "tight" box, the most compact S3 resonator). Cell order in the configuration names: (x0, z0) bottom left, (x0, z1) top left, (x1, z0) bottom right, (x1, z1) top right; "T" marks a TiO₂ cell. All lips point to +z, so the top-row cells face the box lid (4.3 mm away) and the bottom-row cells face the top-row loops; the two rows are not equivalent.
- **Working mode:** the lowest of the four loop-gap modes of the box (172.5 MHz, adjacent cells in antiphase along x and z), 5 MHz below the next mode. The other three modes are also reported (`s4_box_all_modes.csv`).
- **4-cell line** (`L4_*`): the S3 study-1 resonator (TEM unit-cell line, PEC shorts a/2 beyond the slab faces). Its working mode is the MNG–air mode (µ_eff L + 2g ≈ 0), which S3 showed to be predicted by the retrieved µ_eff. It is used here to test the effective-medium model of the mixed medium.
- **Coaxial stack** (`B141_AAAA`): 1 × 4 × 1 cells stacked along their axes in a tight box of the same volume, for the bulk-lattice Q question.
- **Single cell in its tight box** (`B111_*`).

All models mesh only y ≥ 0 and use y = 0 as a PMC symmetry plane, as in S3. Every gap filling is y-symmetric, so each cell keeps its mirror symmetry.

## Temperature model

Palace has no temperature input. The frequency change with temperature T is built from three parts:

τf = τ_ε + τ_gap + τ_geo

- **Permittivity (τ_ε).** Each dielectric's εr is perturbed by its τε. To first order, τ_ε = −½ Σ_i p_e,i τε,i, where p_e,i is the share of the mode's electric energy in dielectric i (Palace domain-energy postprocessing). The formula was checked against finite differences, i.e. eigenmode runs with εr,i(T0 + ΔT) = εr,i (1 + τε,i ΔT) on the same mesh, at ΔT = ±25 K for every configuration and at ±5 and ±12.5 K for the zero-τf design.
- **Copper expansion (τ_geo = −α_Cu).** All copper (loops and box) expands uniformly. Maxwell scaling shifts every frequency by exactly −α_Cu ΔT.
- **Gap width (τ_gap).** The gap width follows the gap dielectric (bonded slab) rather than the copper, which adds τ_gap = +½ Σ_j p_gap,j (α_d,j − α_Cu). Here p_gap,j is the electric-energy share in the gap of cell j, dielectric plus gap air. This first-order term was checked by a run with a 1 % wider gap: Δf/f = +4.764×10⁻³ against +4.753×10⁻³ predicted.

The FD runs also give the curvature c₂ = (f₊ + f₋ − 2f₀)/(f₀ ΔT²) of the permittivity part, which sets the drift over a temperature window once the linear term is cancelled. Nonlinearity of εr(T) itself is not modelled: εr is linear in T, and its true curvature is unknown.

## Material data

| Quantity | Alumina 99.5–99.6 % | Rutile TiO₂ (dense, doped) | Copper |
|---|---|---|---|
| εr | 9.8 | 100 **[VERIFY]** (polycrystalline; single crystal 86 ⊥c, ≈ 170 ∥c) | — |
| tan δ at VHF | 4×10⁻⁵ (measured on 99.5 % grades near 450 MHz) | 1×10⁻⁴ **[VERIFY]** (vendor < 1×10⁻⁴ at 1 MHz; 6×10⁻⁵ at 3 GHz for Mg/Al-doped; undoped dense TiO₂ is much lossier) | σ = 5.8×10⁷ S/m |
| τε | +106 ppm/K **[VERIFY]**, from the dielectric-resonator τf ≈ −60 ppm/K via τε = −2(τf + α) | −916 ppm/K **[VERIFY]**, from τf ≈ +450 ppm/K | — |
| α | 7 ppm/K | 8 ppm/K **[VERIFY]** | 17 ppm/K **[VERIFY]** |

Q_d and τf are linear in tan δ and τε at fixed participations, so `analyze_s4.py` re-evaluates the results for other values without new runs (tables `s4_sensitivity.csv` and `s4_q_tand.csv`). εr of TiO₂ changes the tuned slit width (d ∝ εr at fixed C), not the participation-based results.

## Common numerical settings

- **Copper:** impedance boundary on all copper (loops and box walls), entered as the parallel Palace input Rs_p = 2Rs, Ls_p = 2Rs/ω, with Rs evaluated at the expected mode frequency of each model (149.92 MHz for unit-cell PMC, 179.5 MHz for PEC, 172.5 MHz for the box, 171.5 MHz for the line, 186.6 MHz for the single cell in its box). `Q_c_cu` rescales the copper loss to Rs at the actual mode frequency, as in S3.
- **Copper validation:** the empty 70 × 70 × 140 mm copper box (`E112_o2`) gives Q = 18,508 against 18,510 analytic (TE10ℓ-type mode, Rs at 2.394 GHz).
- **Q split:** Q_d,i = 1/(p_e,i tan δ_i) per dielectric, Q_c = (1/Q − Σ 1/Q_d,i)⁻¹, Q = (1/Q_c,cu + Σ 1/Q_d,i)⁻¹.
- **Mesh:** Gmsh 4.15, second-order tetrahedra, S2/S3 recipe `v2` (0.2 mm on every loop and lip edge, 0.565 mm in and around every gap, 2 mm on copper, 8 mm bulk). The all-alumina meshes are identical to S3's (same element counts).
- **Solver:** FE order 2 in production; order 3 for the convergence check of the key cases. SLEPc Krylov–Schur with shift-and-invert, tolerance 1e-8. Modes are matched between temperature points by frequency and per-cell energy pattern; for the near-degenerate box modes 2 and 3 (0.26 MHz apart) this matching can swap modes, so their FD values are not used.

## Files

| File | Purpose |
|---|---|
| `s4_geometry.py` | Geometry, mesh (cells with per-cell gap fillings, enclosures from S3) and Palace eigenmode configs, including the temperature-perturbed permittivities |
| `build_s4.py` | Builds the meshes, configs and job lists, round by round (r1: cells and references; r1d: TiO₂-cell redesign; r2: mixing studies; r3: zero-τf design, order 3) |
| `run_s4.sbatch` | Slurm job: runs Palace for each run in a list, then `postprocess_s4.py` |
| `postprocess_s4.py` | Inside the job: per-mode f, Q and per-domain energy shares → `results/<run>_eig.csv` |
| `analyze_s4.py` | τf (participation formula, FD), Q split, mixing rules, effective-medium comparison, sensitivity → `results/s4_*.csv`, `s4_summary.json` |
| `s4_effective_medium.py` | Lorentz µ_eff of each meta-atom from its two unit-cell eigenmodes; homogeneous and cell-layered line models; mixture µ_eff = 0 frequency |
| `plot_s4.py` | Figures in `figures/` |
| **`results/s4_box_mixing.csv`** | Count mixing k = 0…4 in the box, all distinct arrangements: f, Q parts, shares, τf parts, FD, drift |
| **`results/s4_box_detuned.csv`** | One TiO₂ cell detuned (continuous share), zero-τf designs, composite-cell box |
| `results/s4_box_all_modes.csv` | τf and shares of all four loop-gap modes of every box configuration |
| `results/s4_line.csv`, `s4_effective_medium.csv` | 4-cell line and the µ_eff model comparison |
| `results/s4_cells.csv` | Single meta-atoms: f, Q parts, participations, τf |
| `results/s4_sensitivity.csv`, `s4_q_tand.csv` | TiO₂ share for τf = 0 under other material data; Q against TiO₂ tan δ |
| `results/s4_lever5.csv` | Q_loaded/Q_cell for 1, 4 (2×1×2) and 4 (coaxial) cells |
| `results/<run>_eig.csv`, `*-timings.txt`, `*-log-excerpt.txt` | Per-run modes, job timings, solver summaries |
| `figures/tau_vs_share.png` (+ `_zoom`), `q_vs_share.png`, `drift_vs_T.png` | Main figures |
| `cases/`, `configs/`, `runs/`, `lists/` | Mesh metadata, Palace configs, run metadata, job lists |

Run names: `<mesh>_o<order>[_pmc|_pec][_dT{p|m}<K>]`, where `_dTp25` means all permittivities at T0 + 25 K.

## How to run

Meshes are not stored (they are large). Regenerate them locally with Python 3, `gmsh` ≥ 4.11 and numpy:

```bash
python build_s4.py r1        # and r1d, r2, r3 (r2/r3 read the r1/r1d results for the tuning)
```

Copy the folder to the cluster and submit each list:

```bash
sbatch -A <ACCOUNT> -p <PARTITION> --ntasks=64 -t 02:00:00 \
       --export=ALL,PALACE_ENV=<PALACE_INSTALL>/env.sh run_s4.sbatch lists/<round>.txt
```

Copy back `results/` and `postpro/*/{eig.csv,domain-E.csv,palace.json}`, then run `python analyze_s4.py` (numpy, matplotlib).
