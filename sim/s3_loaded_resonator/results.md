# S3 results: resonators loaded with the MNG loop-gap meta-atom array

Palace 0.18.1 eigenmode solutions; the models are described in [`README.md`](README.md). All tables are order-2 elements on the `v2` mesh unless marked. "Correct copper" means the `_zx` runs (copper surface impedance Z = Rs (1 + j)); "S2 copper input" means the runs that use the input of `../s2_mng_slab`, which in Palace's parallel Rs‖Ls model carries half the copper loss (README, Common settings).

## 0. Copper surface impedance check (affects the S2 Q values)

| case | copper input | Palace Q | analytic Q | ratio |
|---|---|---|---|---|
| empty box 70 × 70 × 140 mm, lowest mode 2.394 GHz | exact (2Rs ‖ 2Rs/ω) | 18,531 | 18,510 | 1.001 |
| same | S2 (Rs ‖ Rs/ω) | 37,062 | 18,510 | 2.002 |
| empty box 70 × 35 × 140 mm (all copper) | S2 | 27,309 | 13,639 | 2.002 |

Palace adds Rs, Ls and Cs of an `Impedance` boundary in parallel, so the S2 input Rs, Ls = Rs/ω gives Z = Rs (1 + j)/2: half the copper loss. The meta-atom of S2 therefore has a lower Q than S2 reported:

| meta-atom (S2 unit-cell box) | f (MHz) | Q | Q_c | Q_d |
|---|---|---|---|---|
| PMC x-walls (≈ µ pole f_m), S2 input, order 3 (S2 value) | 149.930 | 5,083 | 6,301 | 26,292 |
| PMC x-walls, **correct copper**, order 2 / 3 | 149.922 / **149.918** | 2,811 / **2,847** | 3,148 / 3,193 | 26,300 |
| PEC x-walls (µ_eff = 0, f_mp), **correct copper**, order 2 / 3 | 179.504 / 179.491 | 2,746 / 2,789 | 3,065 / 3,119 | 26,350 |

The frequencies change by < 0.01 %, so the S2 frequencies, band edges and Re µ_eff are unaffected; Im µ_eff (the copper part of the loss) is too small by about 2 in the S2 retrieval.

## 1. Study 1: short-circuited TEM line loaded with N cells

Cases run: N = 1–4 with g = a/2 = 35 mm, and N = 1, 2, 4 with g = 0 (N = 3, g = 0 has predictions only). Every case has exactly N loop-gap modes between 120 and 200 MHz, and the effective-medium model predicts the same number. With the end walls at the slab faces (g = 0) they are N − 1 modes of the lower polariton branch, k = mπ/(N a) (m = 1 … N − 1), plus the µ_eff = 0 mode at f_mp, whose frequency does not depend on N. With empty end sections g = a/2 there are N − 1 lower-branch modes plus one **MNG–air mode inside the MNG band** (162.2 → 171.5 MHz for N = 1 → 4).

### Frequencies: full wave against effective medium (correct copper)

Δ = (f_model − f_fullwave)/f_fullwave. A = retrieved n, z of the N_ret = 4 slab; B-eig-cu = Lorentz medium built from the two meta-atom eigenmodes only.

| N | g (mm) | mode | f full wave (MHz) | f A (MHz) | Δ A | f B-eig-cu (MHz) | Δ B |
|---|---|---|---|---|---|---|---|
| 1 | 35 | MNG–air | 162.196 | 162.101 | −0.058 % | 162.783 | +0.362 % |
| 2 | 35 | lower, m = 1 | 149.233 | 149.211 | −0.015 % | 149.053 | −0.121 % |
| 2 | 35 | MNG–air | 167.103 | 167.070 | −0.020 % | 167.706 | +0.361 % |
| 3 | 35 | lower, m = 1 | 147.777 | 147.724 | −0.036 % | 147.495 | −0.191 % |
| 3 | 35 | lower, m = 2 | 149.687 | 149.679 | −0.006 % | 149.575 | −0.075 % |
| 3 | 35 | MNG–air | 169.779 | 169.772 | −0.004 % | 170.330 | +0.324 % |
| 4 | 35 | lower, m = 1 | 145.563 | 145.581 | +0.012 % | 145.091 | −0.325 % |
| 4 | 35 | lower, m = 2 | 149.222 | 149.233 | +0.008 % | 149.081 | −0.094 % |
| 4 | 35 | lower, m = 3 | 149.799 | 149.813 | +0.009 % | 149.708 | −0.061 % |
| 4 | 35 | MNG–air | 171.466 | 171.470 | +0.002 % | 171.960 | +0.288 % |
| 1 | 0 | µ_eff = 0 | 179.504 | 179.531 | +0.015 % | 179.504 | 0 (by construction) |
| 2 | 0 | lower, k = π/2a | 149.213 | 149.253 | +0.026 % | 149.107 | −0.071 % |
| 2 | 0 | µ_eff = 0 | 179.505 | 179.531 | +0.014 % | 179.504 | −0.001 % |
| 4 | 0 | lower, k = π/4a | 145.913 | 146.007 | +0.065 % | 145.585 | −0.225 % |
| 4 | 0 | lower, k = π/2a | 149.213 | 149.253 | +0.026 % | 149.107 | −0.071 % |
| 4 | 0 | lower, k = 3π/4a | 149.797 | 149.812 | +0.010 % | 149.713 | −0.056 % |
| 4 | 0 | µ_eff = 0 | 179.505 | 179.531 | +0.014 % | 179.504 | −0.001 % |

- **The retrieved medium predicts every mode to ≤ 0.065 %** (≤ 0.10 MHz). Measured against each mode's distance from the µ pole (|f − f_m| = 0.12 … 29.6 MHz), the error is ≤ 0.8 % for the MNG–air and µ_eff = 0 modes, 0.4–5.6 % for the lower-branch modes, and 11–12 % for the two modes nearest the zone edge (0.12 MHz below f_m, error 0.015 MHz). The N_ret = 1, 2 and 4 retrievals agree with each other to ≤ 0.07 % (`results/em_predictions.csv`), which is the homogenization result of S2 seen through a resonator. The full-wave runs with the S2 copper input (same input as the retrieval) give the same picture (`results/study1_line.csv`).
- The Lorentz medium built from the meta-atom alone (no slab retrieval) is within 0.36 %.
- **MNG–air mode = µ_eff design rule.** The retrieved µ_eff at the full-wave frequency is −1.003, −0.512, −0.348, −0.265 for N = 1 … 4, against the quasi-static condition µ_eff L + 2g = 0, i.e. µ_eff = −2g/L = −1, −0.5, −0.333, −0.25. The resonance is fixed by the ratio of air length to MNG length, not by the total length; the resonator is 0.08–0.20 λ long.
- The highest lower-branch mode of each N sits close to the Brillouin-zone edge (k a → π), where the homogeneous Lorentz medium has infinitely many spurious roots; these are removed with a k a < 0.97 π filter. The lattice limits the lower branch to N − 1 modes, as the full-wave result shows.

### Q (study 1)

| N | g (mm) | mode | f (MHz) | Q correct copper | Q / Q_cell | Q S2 copper input | Q model A (S2 input) | Q B-eig-cu |
|---|---|---|---|---|---|---|---|---|
| 1 | 35 | MNG–air | 162.20 | 2,891 | 1.028 | 5,161 | 5,230 | 2,598 |
| 2 | 35 | MNG–air | 167.10 | 2,884 | 1.026 | 5,150 | 5,210 | 2,524 |
| 3 | 35 | MNG–air | 169.78 | 2,869 | 1.020 | 5,126 | 5,180 | 2,487 |
| 4 | 35 | MNG–air | 171.47 | 2,857 | 1.016 | 5,107 | 5,158 | 2,465 |
| 2 | 35 | lower m = 1 | 149.23 | 2,839 | 1.010 | 5,087 | 5,212 | 2,880 |
| 3 | 35 | lower m = 1 | 147.78 | 2,893 | 1.029 | 5,186 | 5,247 | 2,986 |
| 4 | 35 | lower m = 1 | 145.56 | 2,978 | 1.059 | 5,346 | 5,410 | 3,175 |
| 4 | 35 | lower m = 3 | 149.80 | 2,819 | 1.003 | 5,050 | 4,351 | 2,839 |
| 1 | 0 | µ_eff = 0 | 179.50 | 2,746 | 0.977 | 4,924 | 4,971 | 2,352 |
| 4 | 0 | lower k = π/4a | 145.91 | 2,951 | 1.050 | 5,295 | 5,327 | 3,113 |

Q_cell = 2,811 (correct copper) or 5,035 (S2 input), the PMC-box meta-atom at the same order and mesh recipe.

- **Q_loaded ≈ Q_cell in every case: 0.98–1.06 Q_cell.** Q does not fall with N. The modes furthest below the pole keep a little energy in the line and are slightly above Q_cell; the µ_eff = 0 mode is 2 % below.
- **The retrieved medium predicts Q** (like-for-like, S2 copper input) to within +0.6 … +2.5 % for the MNG–air, µ_eff = 0 and m = 1 modes; it fails only for the mode nearest the zone edge (N = 4, m = 3: 4,351 against 5,050), where the retrieval is ill-conditioned.
- A Lorentz medium with constant damping (B) gets Q wrong by −10 … −14 % for the µ_eff ≤ 0 modes and +1 … +7 % for the lower-branch modes. In that model Q(µ = 0) = Q_cell f0/f_mp, while the full wave keeps Q ≈ Q_cell; the copper loss of a loop-gap cell is not a fixed Lorentz damping rate.
- Energy (`pm_*` columns): 52–74 % of the magnetic energy is inside the loop bores; the empty end sections hold 7–9 % of it in the MNG–air modes and < 1 % in the lower-branch modes. 95 % of the electric energy is in the alumina.

## 2. Study 2: compact copper cavity

Fundamental (lowest) loop-gap mode, correct copper. Q_cu rescales the copper loss to Rs(f) at the mode frequency (README); Q is the raw run value.

| box | cells | box (mm) | V (L) | f_loaded (MHz) | Q | **Q_cu** | Q_cu / Q_cell | f_empty (MHz) | **miniaturization f_empty/f** | largest side / λ | magnetic energy in loop bores |
|---|---|---|---|---|---|---|---|---|---|---|---|
| tight 1×1×1 | 1 | 70 × 70 × 70 | 0.343 | 186.62 | 2,018 | 2,183 | 0.78 | 3,028 | **16.2** | 0.044 (λ/23) | 49 % |
| tight 1×1×2 | 2 | 70 × 70 × 140 | 0.686 | 177.46 | 2,037 | 2,172 | 0.77 | 2,394 | **13.5** | 0.083 (λ/12) | 54 % |
| tight 2×1×2 | 4 | 140 × 70 × 140 | 1.372 | 172.49 | 2,060 | 2,177 | 0.77 | 1,514 | **8.8** | 0.081 (λ/12) | 57 % |
| roomy 1×1×1 | 1 | 119 × 119 × 119 | 1.685 | 165.86 | 2,569 | 2,676 | 0.95 | 1,781 | **10.7** | 0.066 (λ/15) | 61 % |
| roomy 1×1×2 | 2 | 119 × 119 × 189 | 2.676 | 161.83 | 2,559 | 2,641 | 0.94 | 1,489 | **9.2** | 0.102 (λ/9.8) | 64 % |

- f_empty is the lowest mode of the same box without cells (analytic; TE101-type). Its copper Q is 18,500–25,300.
- Second loaded modes: 1×1×2 at 186.79 MHz (Q_cu 2,247); 2×1×2 at 177.46, 177.73 and 186.79 MHz; roomy 1×1×2 at 168.54 MHz.
- **Where the loaded modes sit.** In the tight boxes the loaded modes lie at 172–187 MHz, around f_mp (179.5 MHz): a closed copper box enforces zero net flux, so the loops return their flux through the 4–12 mm of air between cells and walls, which is the µ_eff ≈ 0 condition. More cells let neighbours carry each other's return flux, and the fundamental moves down (186.6 → 177.5 → 172.5 MHz). With room around the array (walls ≈ a/2 from the cells) the fundamental drops into the MNG band (162–166 MHz). This is the same MNG–air mode as in study 1 (162.2 MHz for one cell with 35 mm of air on each side).
- **Wall loss.** The 1×1×2 tight box with lossless walls has Q_cu 2,495, so the copper walls take 13 % of the loss (Q_walls ≈ 16,800). The loops themselves lose more than in the open lattice: 2,495 against the meta-atom's 2,811–2,978, because the 12.5 mm end gaps squeeze the return flux.
- **Q against the meta-atom:** Q_cu = 0.77–0.78 Q_cell in the tight boxes and 0.94–0.95 Q_cell with walls a/2 from the cells.
- **Stored energy:** 95 % of the electric energy is in the alumina gaps in every case. 49–64 % of the magnetic energy is inside the loop bores; the rest is the return flux in the air between loops and walls. In the tight boxes that air is part of the unit cells' own lattice volume, so effectively all of the energy is in the cells; nothing is stored in a cavity field of the empty-box type.

## 3. Study 3: granularity (cells scaled by 0.8)

Scaled cell (a = 56 mm, r = 17.2 mm, length 36 mm, 0.5 mm alumina, s = 4.91 mm), correct copper:

| meta-atom | a (mm) | f_m (MHz) | Q | Q_c | f_mp (MHz) | F = 1 − (f_m/f_mp)² |
|---|---|---|---|---|---|---|
| S2 cell | 70 | 149.922 | 2,811 | 3,148 | 179.504 | 0.302 |
| scaled cell | 56 | 149.943 | 2,281 | 2,502 | 179.544 | 0.303 |
| ratio | 0.8 | | **0.811** | **0.795** | | |

Same-volume cavities (tight boxes, correct copper, fundamental mode):

| pair | box | cells | V (L) | f (MHz) | Q | Q_cu | Q ratio (fine/coarse) | N^(−1/3) |
|---|---|---|---|---|---|---|---|---|
| 1 | 1×1×2, a = 70 | 2 | 0.686 | 177.46 | 2,037 | 2,172 | | |
| 1 | 2×1×2, a = 56 | 4 | 0.702 | 172.51 | 1,662 | 1,758 | **0.816 (Q), 0.809 (Q_cu)** | 0.794 |
| 2 | 2×1×2, a = 70 | 4 | 1.372 | 172.49 | 2,060 | 2,177 | | |
| 2 | 2×2×2, a = 56 | 8 | 1.405 | 162.72 | 1,856 | 1,921 | 0.90 (Q), 0.88 (Q_cu) | 0.794 |

- **Pair 1 follows the granularity law:** doubling N at fixed volume costs a factor 0.81 in Q, against N^(−1/3) = 0.794. The ratio equals the meta-atom ratio (0.811): Q_loaded ≈ Q_cell and Q_cell ∝ cell size. The copper part alone scales as 0.795 ≈ 0.8.
- Pair 2 is not like-for-like: the 2×2×2 box stacks two loops coaxially along y, which relieves the end-wall flux squeeze and moves the fundamental to 162.7 MHz (Q_cu/Q_cell = 0.84 against 0.77 for 2×1×2). Its ratio (0.88–0.90) mixes granularity with the arrangement.

## 4. Convergence (Q_c, order 2 → 3, same mesh)

| case | copper | f change | Q_c order 2 → 3 | change |
|---|---|---|---|---|
| meta-atom PMC box | correct | −0.002 % | 3,148 → 3,193 | +1.4 % |
| meta-atom PEC box | correct | −0.007 % | 3,065 → 3,119 | +1.8 % |
| study 1, N = 2, g = 35, mode 1 (lower) | correct | −0.002 % | 3,179 → 3,223 | +1.4 % |
| study 1, N = 2, g = 35, mode 2 (MNG–air) | correct | −0.005 % | 3,237 → 3,290 | +1.6 % |
| study 2, tight 1×1×2, mode 1 | correct | −0.015 % | 2,207 → 2,248 | +1.8 % |
| study 2, tight 1×1×2, mode 2 | correct | −0.015 % | 2,255 → 2,296 | +1.8 % |
| study 1, N = 2, g = 35, modes 1 / 2 | S2 input | −0.002 / −0.006 % | 6,291 → 6,378 / 6,399 → 6,501 | +1.4 / +1.6 % |
| study 2, tight 1×1×2, mode 1 | S2 input | −0.015 % | 4,359 → 4,438 | +1.8 % |

All ≤ 2 %, the criterion. Order-2 Q values are therefore 1–2 % low; ratios between cases (Q/Q_cell, granularity) are taken at the same order and cancel this.

Half-model check: the y = 0 symmetry plane reproduces the S2 full-model eigenmodes (same S2 copper input, order 2) to 2e-6 in f and 0.4 % in Q (149.9341 against 149.9338 MHz, Q 5,035 against 5,016; 179.5164 against 179.5164 MHz, Q 4,924 against 4,905).

## Cost

All runs on LANTA's free `compute-limited` partition (CPU billing weight 0.0): **0.000 SHr** charged. 11 jobs with 64 MPI ranks each (one 8-rank check), 47 runs, ≈ 500 core-hours in total (≈ 3.9 node-hours). Single runs took 4 s – 40 min at order 2 and 7–34 min at order 3 (`results/*-timings.txt`).

## Not modelled

Cell supports and the slab holder, Rs(f) inside a run (corrected afterwards to first order), surface roughness, the TiO₂ compensation layer, and coupling ports (all Q values are unloaded).
