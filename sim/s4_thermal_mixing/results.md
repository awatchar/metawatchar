# S4 results: thermal self-compensation by mixing meta-atoms

Setup, material data and the temperature model are in [README.md](README.md). Unless stated otherwise, values are FE order 2. τf includes the copper expansion (−17 ppm/K) and the gap-width term (≈ −4.6 ppm/K). The resonator is the 2 × 1 × 2 tight copper box (1.372 L) and its working mode is the lowest loop-gap mode, at about 172.5 MHz. Values marked **[VERIFY]** in the README (TiO₂ εr, tan δ, τε; alumina τε; α) carry through as shown in [Sensitivity](#sensitivity-to-the-verify-material-data).

## Headline

1. **One mixing rule holds for every configuration: τf = −72.8 + 491.5 · x_TiO₂ ppm/K.** Here x_TiO₂ is the TiO₂ share of the electric energy stored in the gaps. The rule covers 20 configurations: count mixing k = 0…4 in all distinct arrangements, the detuned series, the composite-gap cells, and the 4-cell line. It fits all of them to within 1.6 ppm/K. **τf = 0 at x_TiO₂ = 14.8 %.** This agrees with the 14.7 % capacitance share of the analytic cell model and with the composite cell, whose 14.8 % share gives τf = −0.002 ppm/K.
2. **Count mixing overshoots.** One TiO₂ cell out of four takes 27–29 % of the gap energy, more than the count fraction of 25 %, and gives τf = +60 to +70 ppm/K depending on its position. Two cells take 52–57 % (+182 to +205 ppm/K) and four cells give +420 ppm/K.
3. **The share is tuned continuously by detuning the TiO₂ cell's own resonance upward.** With one TiO₂ cell detuned by f_T/f_A = 1.018 (slit 3.41 mm instead of 3.29 mm), the share is 13.5 % and τf = −6.5 ppm/K (order 3: −6.3). The zero lies at f_T/f_A ≈ 1.016 (slit ≈ 3.39 mm). The knob is steep: about −3.3 ppm/K per 10 µm of TiO₂ slit width. Between the matched and the zero design the resonator frequency moves by only +0.6 MHz (about 1/5 of the TiO₂-cell shift), so the TiO₂ cell's tuning acts as a τf trim that is largely independent of f.
4. **The price of medium-level mixing is curvature.** The TiO₂ cell drifts against the alumina cells at about 490 ppm/K, so its energy share changes with temperature. With the best linear trim, the frequency excursion (scaled to 145.2 MHz) is **8.7 kHz peak to peak over ±10 K, 19.6 kHz over ±15 K, 35 kHz over ±20 K and 56 kHz over ±25 K**. A ±15 kHz window therefore holds for about **±18.5 K** around the design temperature. The meta-atom-level reference (composite gap in every cell) is flat to < 2 kHz over ±25 K, because all cells drift together.
5. **Q at τf ≈ 0 is 2,162 against 2,180 for all-alumina (−0.8 %; order 3: 2,198 against 2,217).** That assumes TiO₂ tan δ = 1×10⁻⁴. The TiO₂ cell (full-length TiO₂ slab, no lips) has a slightly higher copper Q than the alumina cell, so the cost is almost entirely the TiO₂ dielectric loss: Q_d,TiO₂ = 77,900 at tan δ 1×10⁻⁴. At 3×10⁻⁴ the Q is 2,048 (−6 %); at 1×10⁻³ it is 1,730 (−21 %).
6. **The effective-medium model works where the mode is a homogeneous-medium mode.** In the 4-cell line, a cell-layered µ_eff model built only from each meta-atom's two unit-cell eigenmodes predicts τ_ε within 0.4 ppm/K, including the position of the TiO₂ cell (end: 77.7 against 77.4 ppm/K; inside: 67.8 against 67.5). It predicts the mode frequency within 0.29 %. The homogeneous volume-average µ_eff misses the position effect by about ±5 ppm/K. The box's working mode is a zone-corner (antiphase) mode, which a homogeneous µ_eff cannot describe: the µ_eff = 0 mixture rule under-predicts its τ_ε by 8–32 ppm/K. The energy-share rule (1) is exact for both.
7. **Lever 5 (bulk-lattice Q) is not confirmed for 4 cells.** Q_loaded/Q_cell is 0.778 for one cell in its tight box and 0.775 for the 2 × 1 × 2 box. A coaxial 1 × 4 × 1 stack of the same volume reaches 0.858 (+10.6 % Q, +9 % copper Q after frequency scaling), far from the ×1.47 upper bound. In the stack, 81 % of the energy sits in the inner pair of cells.

## 1. Meta-atoms

Unit-cell box (S2 eigen box), PMC walls (µ pole f_m). Q values use the copper loss at the mode frequency. Q_c is scaled to 149.9 MHz (Q_c ∝ √f) for comparison.

| Cell | f_m (MHz) | Q | Q_c (at 149.9 MHz) | Q_d | p_e in dielectric | τ_ε (ppm/K) | τf (ppm/K) |
|---|---|---|---|---|---|---|---|
| A (alumina, S2/S3) | 149.922 | 2,811 (order 3: 2,846) | 3,147 | 26,300 | 0.951 | −50.4 (FD −50.4) | **−72.1** |
| **W3.29 (TiO₂ slab 3.29 mm, no lips)** | 150.012 | 2,494 (order 3: 2,520) | **3,283** | 10,370 (TiO₂) | 0.965 | +441.8 (FD +442.0) | **+420.5** |
| T3.6 (TiO₂ block 3.6 mm in the A gap) | 155.943 | 2,243 | 2,719 | 11,740 (TiO₂) | 0.852 (+0.096 gap air) | +390.1 (FD +390.2) | +368.8 |
| C0.75 (alumina gap, 0.75 mm TiO₂ block) | 140.017 | 2,645 | 3,124 | 30,670 / 70,480 | 0.815 / 0.142 | +21.8 (FD +21.8) | **−0.002** |

- **TiO₂-cell design.** The first design (`T`, a short TiO₂ block in the alumina cell's 1 mm gap) loses 14 % of the copper Q. The loop current has to converge on the 3–4 mm long capacitor (current crowding along the lips). A TiO₂ slab filling the whole 45 mm slit, with the slit widened to 3.29 mm and the lips removed (`W`), has **4 % more** copper Q than the alumina cell. It also matches the alumina cell's f_m to 0.06 %. The `W` cell was used for all mixing studies.
- **TiO₂-cell geometry (W3.29).** It uses the same loop (r 21.5 mm, wall 3 mm, length 45 mm) with a 3.29 mm slit through the wall and no lips. The TiO₂ is a 45 × 3.0 × 3.29 mm bar. Tuning: 1/f² = a + b/d (from d = 3.0 and 3.6 mm: 143.45 and 156.69 MHz). In its own tight box it resonates at 186.02 MHz (alumina cell: 186.61 MHz), with Q 2,003 against 2,187.
- **PEC box (µ_eff = 0):** A cell 179.501 MHz, W cell 179.106 MHz. The Lorentz filling factor is the same for both cells (same loop).
- **Checks.**
  - The participation formula agrees with FD (±25 K) to ≤ 0.2 ppm/K for every single cell.
  - The gap-width term was checked with a 1 % wider gap: Δf/f = +4.764×10⁻³ against +4.753×10⁻³ predicted.
  - The copper input was validated on an empty box: Q 18,508 against 18,510 analytic.

## 2. Count mixing in the 2 × 1 × 2 box (TiO₂ cells matched, W3.29)

| k | Arrangement (TiO₂ cells) | f (MHz) | x_TiO₂ | τf (ppm/K) | Q (tan δ_TiO₂ 1e-4) | Q_c | Q_d,Al | Q_d,TiO₂ | per-cell gap energy (BL TL BR TR) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | — | 172.484 | 0 | −72.1 | 2,180 | 2,377 | 26,300 | — | 0.222 0.253 0.222 0.253 |
| 1 | bottom row | 172.329 | 0.272 | +60.4 | 2,131 | 2,407 | 35,990 | 38,550 | **0.259** 0.252 0.214 0.229 |
| 1 | top row | 172.405 | 0.291 | +69.8 | 2,122 | 2,401 | 36,940 | 35,980 | 0.219 **0.278** 0.210 0.248 |
| 2 | bottom row | 172.181 | 0.522 | +182.8 | 2,088 | 2,435 | 54,590 | 20,030 | **0.250** 0.229 **0.250** 0.229 |
| 2 | top row | 172.326 | 0.567 | +205.4 | 2,068 | 2,424 | 60,230 | 18,400 | 0.208 **0.272** 0.208 **0.272** |
| 2 | one column | 172.249 | 0.555 | +199.2 | 2,076 | 2,430 | 58,580 | 18,820 | **0.255 0.276** 0.203 0.224 |
| 2 | diagonal | 172.258 | 0.521 | +182.4 | 2,083 | 2,428 | 54,450 | 20,060 | **0.246** 0.247 0.212 **0.252** |
| 3 | alumina bottom-left | 172.180 | 0.791 | +316.3 | 2,032 | 2,451 | 124,400 | 13,150 | 0.201 **0.247 0.243 0.271** |
| 3 | alumina top-left | 172.109 | 0.766 | +303.6 | 2,042 | 2,456 | 111,100 | 13,590 | **0.237** 0.225 **0.247 0.252** |
| 4 | all | 172.040 | 1 | +420.2 | 1,999 | 2,476 | — | 10,370 | 0.235 0.248 0.235 0.248 |

- **Energy non-uniformity.** In the all-alumina box, the top-row cells (lips facing the lid) hold 14 % more gap energy than the bottom row. A matched TiO₂ cell takes more than its count share: its tight-box resonance is 0.3 % lower. The share then depends on position: 27.2 % at the bottom, 29.1 % at the top.
- **Other modes.** For each configuration, the other three loop-gap modes (177.0–197.7 MHz) have their own shares and τf, from −72 to +422 ppm/K (`s4_box_all_modes.csv`). Medium-level mixing compensates one chosen mode, not the whole spectrum. The composite box compensates all four modes (+0.23 to +0.29 ppm/K).
- **FD check.** At ±25 K the central difference differs from the participation formula by up to 4 ppm/K in the matched k = 1 boxes. This comes from the strong temperature nonlinearity (cubic term). At ±5 K the difference is 0.10–0.11 ppm/K (zero-τf designs), so the formula is the exact first derivative.

## 3. Continuous share: one TiO₂ cell detuned upward (bottom-left position)

| TiO₂ slit (mm) | f_T/f_A (cell) | f (MHz) | x_TiO₂ | τf (ppm/K) | c₂ ±25 K (ppb/K²) | Q (1e-4) | per-cell gap energy |
|---|---|---|---|---|---|---|---|
| 3.29 | 1.001 | 172.329 | 0.272 | +60.4 | −2,525 | 2,131 | 0.259 0.252 0.214 0.229 |
| **3.41** | **1.018** | **172.902** | **0.135** | **−6.5** (FD ±5 K −6.4; order 3 −6.3) | −1,215 | **2,162** (order 3 2,198) | **0.128** 0.266 0.228 0.329 |
| 3.43 | 1.021 | 172.965 | 0.121 | −13.1 | −1,063 | 2,164 | 0.115 0.266 0.229 0.341 |
| 3.49 | 1.030 | 173.113 | 0.090 | −28.5 | −714 | 2,171 | 0.085 0.266 0.229 0.370 |
| 3.71 | 1.061 | 173.409 | 0.038 | −53.6 | −206 | 2,184 | 0.036 0.258 0.225 0.431 |
| 4.00 | 1.100 | 173.571 | 0.018 | −63.3 | −62 | 2,189 | 0.017 0.250 0.220 0.463 |

- **Zero.** τf = 0 is reached at d ≈ 3.39 mm (f_T/f_A ≈ 1.016), x_TiO₂ = 14.8 %. Detuning makes the field strongly non-uniform: the diagonally opposite alumina cell holds 33 % of the gap energy, against 12.8 % in the TiO₂ cell. At 2 W this is harmless (gap voltage +15 % against uniform). It matters for the 100 W scaling path.
- **Temperature curve of the d = 3.41 mm design** (7 FD points, −25 to +25 K, quartic fit; best linear trim):

| Window | ±5 K | ±10 K | ±13 K | ±15 K | ±20 K | ±25 K |
|---|---|---|---|---|---|---|
| Peak-to-peak at 145.2 MHz | 2.2 kHz | 8.7 kHz | 14.7 kHz | 19.6 kHz | 35.2 kHz | 55.6 kHz |

  The d = 3.43 mm design gives 7.6 / 17.1 / 30.8 / 48.6 kHz over ±10 / ±15 / ±20 / ±25 K. The curvature near T0 is −0.59 ppm/K² (quartic fit). Its magnitude grows towards larger |ΔT|.
- **Composite reference (meta-atom level).** The 2 × 1 × 2 box of C0.75 cells has x_TiO₂ = 14.9 %, τf = +0.26 ppm/K, c₂ = +1.2 ppb/K² and Q 2,057 at 161.0 MHz. Its copper Q, scaled to 145 MHz, is 2,159 (all-alumina: 2,179). Peak-to-peak drift is < 2 kHz over ±25 K.

## 4. Effective-medium (µ_eff) model

Each meta-atom gets a Lorentz µ from its own PMC- and PEC-box eigenmodes (S3 model B-eig). The cell τf enters through those two eigenfrequencies. The comparison is for the permittivity part τ_ε, because the expansion terms are common.

| Resonator | Full wave f (MHz) / τ_ε | Cell-layered µ_eff | Homogeneous ⟨µ⟩ |
|---|---|---|---|
| Line, all alumina | 171.464 / −50.11 | 171.958 (+0.29 %) / −50.17 | same |
| Line, TiO₂ cell at the end | 171.436 / +77.75 | 171.893 (+0.27 %) / +77.39 | 171.896 / +72.31 |
| Line, TiO₂ cell inside | 171.420 / +67.80 | 171.898 (+0.28 %) / +67.54 | 171.896 / +72.31 |
| Box, k = 1 bottom (working mode / in-phase mode) | 172.33 / +85.9 ; 186.59 / +68.6 | — | µ_eff = 0 rule: 179.40 / +72.1 |
| Box, k = 2 (four arrangements, working mode) | 172.2–172.3 / +202.7…+227.0 | — | +194.6 |
| Box, all TiO₂ | 172.04 / +441.8 | — | +441.2 |

- The ≈ 0.28 % frequency offset is the B-eig model's own error (S3: ≤ 0.36 %). It is common to all line configurations, so the mixing-induced change is predicted to 0.01 %.
- **"Designed through µ_eff" holds at the cell-resolved level.** The TiO₂ cell enters the medium as a frequency-shifted Lorentz term with the same filling factor, and the resonator's τf follows to < 0.5 ppm/K. A fully homogenized mixture is accurate only to the extent that the TiO₂ cells sit where the field is average.

## 5. Q against TiO₂ share

Q of the working mode in the box, by TiO₂ tan δ [VERIFY] (alumina tan δ 4×10⁻⁵):

| x_TiO₂ | Q_c | Q (6e-5) | Q (1e-4) | Q (3e-4) | Q (1e-3) |
|---|---|---|---|---|---|
| 0 (all alumina) | 2,377 | 2,180 | 2,180 | 2,180 | 2,180 |
| **0.135 (τf ≈ 0, d 3.41)** | 2,399 | 2,186 | **2,162** | 2,048 | 1,730 |
| 0.272 (k = 1 matched) | 2,407 | 2,175 | 2,131 | 1,917 | 1,421 |
| 0.522 (k = 2) | 2,435 | 2,179 | 2,088 | 1,726 | 1,074 |
| 1 (all TiO₂) | 2,476 | 2,166 | 1,999 | 1,443 | 731 |

Q_c rises slightly with the TiO₂ share because the W cell has a higher copper Q than the alumina cell. Q at τf = 0 stays within 1 % of all-alumina for tan δ_TiO₂ ≤ 1×10⁻⁴. At 1×10⁻³ the compensated resonator loses 21 %. **TiO₂ tan δ at VHF remains the critical unknown for Q.** The order-3 values are 1.7 % higher (2,217 and 2,198).

## 6. Lever 5: does a 4-cell medium approach lattice Q?

| Resonator (same cell) | V (L) | f (MHz) | Q | Q/Q_cell | Q_c/Q_c,cell (√f scaled) | Energy distribution |
|---|---|---|---|---|---|---|
| Unit cell (lattice, PMC box) | 0.343 | 149.92 | 2,811 | 1 | 1 | — |
| 1 cell, tight box | 0.343 | 186.61 | 2,187 | 0.778 | 0.679 | — |
| 2 × 1 × 2, tight box | 1.372 | 172.48 | 2,180 | 0.775 | 0.704 | 0.22–0.25 per cell |
| 1 × 4 × 1 coaxial stack, tight box | 1.372 | 180.56 | 2,411 | **0.858** | 0.769 | inner pair 81 %, outer pair 14 % |

A 4-cell medium in a tight box does not approach the lattice value. All cells are boundary cells. Coaxial stacking, where the loops share flux along their axes, recovers about +10 % (S3 had already hinted at this). The ×1.47 upper bound would need interior cells, i.e. larger arrays.

## 7. Sensitivity to the [VERIFY] material data

Re-evaluated from the same participations (no new runs):

| Change | x_TiO₂ for τf = 0 | τf, all alumina | τf, all TiO₂ |
|---|---|---|---|
| Baseline (τε TiO₂ −916, Al₂O₃ +106 ppm/K; α 17 / 7 / 8) | **14.9 %** | −72.1 | +420.2 |
| τε(TiO₂) −700 / −800 / −1000 / −1100 ppm/K | 18.8 / 16.8 / 13.8 / 12.7 % | — | +316 / +364 / +461 / +509 |
| τε(Al₂O₃) +80 / +130 ppm/K | 12.7 / 16.8 % | −59.8 / −83.5 | — |
| α(Cu) 16.5 / 17.7 ppm/K | 14.7 / 15.1 % | −71.4 / −73.2 | |
| α(TiO₂) 5 / 9 ppm/K | 14.9 / 14.9 % | | |

- **Scaling.** The required share is x* ≈ |τf,A| / (τf,T − τf,A). The alumina τε matters as much as the TiO₂ τε. A ±10 % error in either moves x* by about ±1.5–2 points.
- **Trim.** Detuning the TiO₂ cell trims τf at about 3.3 ppm/K per 10 µm of slit width, or about 2.3 ppm/K per 0.1 % of TiO₂-cell frequency near the zero. Any of these errors can be trimmed after assembly with the TiO₂ cell's own tuning. Measuring τf of one A cell and one TiO₂ cell first is the efficient way to set it.
- **TiO₂ εr** (100 assumed) changes the TiO₂ slab width needed (d ∝ εr), not the share rule.
- **TiO₂ tan δ:** see section 5.

## 8. Convergence (order 2 → 3)

| Case | f (MHz) | Q | τf (ppm/K) |
|---|---|---|---|
| A cell, PMC box | 149.9220 → 149.9184 | 2,811 → 2,846 (+1.3 %) | −72.13 → −72.15 |
| W3.29 cell, PMC box | 150.0116 → 149.9947 | 2,494 → 2,520 (+1.0 %) | +420.51 → +420.54 |
| Box, all alumina | 172.4842 → 172.4580 | 2,180 → 2,217 (+1.7 %) | −72.14 → −72.15 |
| Box, k = 1, d 3.41 (τf ≈ 0 design) | 172.9025 → 172.8742 | 2,162 → 2,198 (+1.7 %) | −6.52 → −6.28 (FD ±5 K: −6.41 → −6.17) |

τf converges to 0.25 ppm/K, and f to 0.02 %. Q is 1–2 % low at order 2, as in S3.

## Cost

- 15 productive Slurm jobs (64 or 128 cores; plus three that stopped within seconds on a list-file line-ending error), 89 Palace eigenmode runs, all on the billed `compute-devel` partition. The free partition was estimated at a 22–32 h wait each time it was checked.
- **7.34 SHr** (elapsed × cores / 128).
- Typical run times at 64 ranks: 3–5 min per single cell, 9–16 min per 4-cell box (order 2, ≈ 3 M DOF), 10 min per 4-cell line. Order-3 box runs took 19–32 min at 128 ranks.
