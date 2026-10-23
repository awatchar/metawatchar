# LGR-A size sweep: results

> **Correction (copper surface impedance).** Palace's `Impedance` boundary combines Rs and Ls in
> parallel (1/Z = 1/Rs + 1/(jωLs)). The sweep entered Rs and Ls = Rs/ω, which is
> Z = Rs(1 + j)/2: **half the copper loss**. Every Q_c below is about 2× too high and every Q about
> 1.8× too high; frequencies, tuning, C, fields at a given stored energy and Q_d are unaffected.
> A subset was re-run with the correct input (Rs_p = 2 Rs, Ls_p = 2 Rs/ω), see the next section.
> The Palace/analytic Q_c ratio ("REALISM") is about 0.63, not 1.27.

## Corrected copper input (subset, `_cu2` runs, v2 mesh, order 3)

| V | a | d (mm) | s (mm) | f0 (MHz) | Q_c | Q_d | **Q** | First runs: Q_c / Q | Q_c ratio | Palace/analytic Q_c |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 L | 0.5 | 0.635 | 10.92 | 145.269 | 1,839 | 25,621 | **1,716** | 3,635 / 3,184 | 0.506 | 0.63 (analytic 2,942) |
| 0.5 L | 0.5 | 0.635 | 5.51 | 145.332 | 2,353 | 25,862 | **2,157** | 4,647 / 3,939 | 0.506 | 0.63 (analytic 3,761) |
| 0.5 L | 0.5 | 1.0 | 10.42 | 145.311 | 2,344 | 25,858 | **2,149** | 4,630 / 3,927 | 0.506 | — |

- Rs is evaluated at each case's f0 (145.29–145.35 MHz). f0 is 15–21 kHz lower than in the first
  runs (the surface reactance doubled).
- The other cases of the sweep were not re-run; their Q_c scale by about 0.506 as well (the same
  ratio holds for all three re-run cases and for the LGR-A can model with loops and fillets in
  `../lgr_a_driven`), so Q ≈ (1/(0.506 Q_c,old) + 1/Q_d)⁻¹ is a good estimate for them.
- Order-3 values are still 1–2 % low (convergence section below; order 3 → 4 gives +0.6–0.8 %).
- Files: `results/cu2.csv` (same columns as the other round files), configs `configs/*_cu2.json`.


## Summary

- **Reaching 145.5 MHz.** The cell needs C ≈ 65–135 pF across the gap. The d = 0.635 mm alumina
  gap gives C ≈ eps0 eps_r z (s + s0) / d with s0 = 3.2–3.4 mm (wall depth t = 3 mm plus about
  0.3 mm of fringing), so the lip height s sets the frequency.
  - **V = 0.15 L: no** case reaches 145.5 MHz. Even with the lips ending 3 mm from the shield, f0 is 192–204 MHz. The lips would have to extend past the shield wall.
  - **V = 0.25 L: only a = 0.5** reaches it (s = 10.92 mm, 3.15 mm clearance). a = 0.4 and a = 0.6 stop at 150 and 152 MHz. They would need about 1.9 mm clearance; these cases are tabulated separately.
  - **V = 0.5 L: all a** reach it, with 10–14 mm clearance at d = 0.635 mm and 6–8 mm at d = 1.0 mm.
  - **d = 1.0 mm** (thicker alumina) fits only at 0.5 L. At 0.25 L, a = 0.5, it would need s ≈ 18.9 mm against s_max = 11.07 mm.
- **Q.** Q_total is 3,180 at 0.25 L (a = 0.5) and 3,650–3,940 at 0.5 L (best a = 0.5).
  - The d = 1.0 mm gap costs about 1 % of Q or less.
  - The alumina holds 96–98 % of the electric energy, so Q_d ≈ 25,600–26,100 ≈ 1/tan_delta. Conductor loss dominates: Q ≈ (0.85–0.88) Q_c.
  - The shield carries about 15 % of the conductor loss (0.25 L, a = 0.5). The rest is in the loop and lips.
- **Fields at 10 uJ stored energy.**
  - The gap field is uniform: 0.69–0.88 kV/mm at d = 0.635 mm and 0.44–0.56 kV/mm at d = 1.0 mm. The gap voltage is 390–560 V.
  - In air, 0.1 mm from the slab faces and the copper/alumina/air triple lines, the field is ≤ 0.69 kV/mm.
  - Fields scale as sqrt(W).
- **Analytic model.** The coaxial-return LGR formula (no end caps) underestimates Q_c for this closed-can geometry. Palace/analytic = 1.23–1.32 (mean 1.27).
  - The model squeezes all return flux through the loop–shield annulus. In LGR-A (z = 1.5 R, end clearance 0.25 R), the return flux spreads through the end regions.
  - The effective inductance is L ≈ 0.75 mu0 pi r^2 / (z (1 + rho)), with rho = r^2/(R^2 - (r + t)^2).
- **Accuracy.** Conductor Q converges slowly because the magnetic field is singular at the square tube rims and lip corners. The production mesh (v2) refines those edges to 0.2 mm.
  - On the convergence case, raising the FE order from 3 to 4 on the v2 mesh changes Q_c by +0.8 %.
  - Halving the edge size to 0.1 mm changes Q_c by +1.6 % at order 2.
  - The reported v2 order-3 Q_c values are therefore estimated to be 1–2 % low. Without the edge refinement (v1), Q_c was 7–9 % low in all nine cases.
  - f0 converged to < 0.1 % on every mesh.

## Details

Palace 0.18.1 eigenmode solutions. Copper is an impedance boundary (Rs = 3.147 mOhm/sq, Ls = Rs/omega, sigma = 5.8e7 S/m at 145.5 MHz) on the shield, loop and lips; the gap is alumina (eps_r 9.8, tan_delta 4e-5). Q_c = conductor-only Q (shield + loop + lips), Q_d = 1/(p_d tan_delta) with p_d the electric-energy participation of the alumina, Q = total = (1/Q_c + 1/Q_d)^-1. Fields are peak phasor amplitudes for a stored energy of **10 uJ**; V gap = E gap x d. Edge probes sit 0.1 mm from a metal/dielectric/air corner: the field at a sharp metal edge is singular, so these values depend on the probe distance and, in hardware, on the edge radius. Analytic Q_c: coaxial-return LGR model without end caps (REALISM = 1, see `analyze.py`), same R, r, t; ratio = Palace / analytic.

### Tuned cases, d = 0.635 mm (FE order 3, clearance >= 3 mm)

| V (L) | a | R (mm) | r (mm) | z (mm) | s (mm) | clearance (mm) | f0 (MHz) | Q_c | Q_d | Q | p_d | E gap (kV/mm) | V gap (V) | max air E near edges (kV/mm) | Q_c analytic | Palace/analytic |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 0.5 | 34.14 | 17.07 | 51.21 | 10.92 | 3.15 | 145.29 | 3,635 | 25,620 | **3,184** | 0.976 | 0.71 | 448 | 0.56 (air, slab bore face) | 2,942 | 1.236 |
| 0.5 | 0.4 | 43.01 | 17.21 | 64.52 | 8.53 | 14.28 | 145.34 | 4,316 | 25,605 | **3,693** | 0.976 | 0.69 | 439 | 0.55 (air, slab bore face) | 3,278 | 1.317 |
| 0.5 | 0.5 | 43.01 | 21.51 | 64.52 | 5.51 | 13.00 | 145.35 | 4,647 | 25,862 | **3,939** | 0.967 | 0.80 | 509 | 0.64 (air, slab bore face) | 3,761 | 1.236 |
| 0.5 | 0.6 | 43.01 | 25.81 | 64.52 | 4.06 | 10.15 | 145.34 | 4,551 | 26,136 | **3,876** | 0.957 | 0.88 | 556 | 0.69 (air, slab bore face) | 3,562 | 1.278 |

### Tuned cases, d = 1 mm (FE order 3, clearance >= 3 mm)

| V (L) | a | R (mm) | r (mm) | z (mm) | s (mm) | clearance (mm) | f0 (MHz) | Q_c | Q_d | Q | p_d | E gap (kV/mm) | V gap (V) | max air E near edges (kV/mm) | Q_c analytic | Palace/analytic |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 0.4 | 43.01 | 17.21 | 64.52 | 15.12 | 7.69 | 145.34 | 4,261 | 25,620 | **3,653** | 0.976 | 0.44 | 438 | 0.38 (air, slab bore face) | 3,278 | 1.300 |
| 0.5 | 0.5 | 43.01 | 21.51 | 64.52 | 10.42 | 8.09 | 145.33 | 4,630 | 25,858 | **3,927** | 0.967 | 0.51 | 508 | 0.44 (air, slab bore face) | 3,761 | 1.231 |
| 0.5 | 0.6 | 43.01 | 25.81 | 64.52 | 8.15 | 6.06 | 145.33 | 4,534 | 26,131 | **3,864** | 0.957 | 0.55 | 555 | 0.48 (air, slab bore face) | 3,562 | 1.273 |

### Near-miss cases tuned with the clearance relaxed below 3 mm (violate the 3 mm rule)

| V (L) | a | R (mm) | r (mm) | z (mm) | s (mm) | clearance (mm) | f0 (MHz) | Q_c | Q_d | Q | p_d | E gap (kV/mm) | V gap (V) | max air E near edges (kV/mm) | Q_c analytic | Palace/analytic |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 0.4 | 34.14 | 13.66 | 51.21 | 15.63 | 1.85 | 145.21 | 3,330 | 25,495 | **2,945** | 0.981 | 0.61 | 387 | 0.49 (air, slab bore face) | 2,591 | 1.285 |
| 0.25 | 0.6 | 34.14 | 20.48 | 51.21 | 8.71 | 1.95 | 145.25 | 3,519 | 25,835 | **3,097** | 0.968 | 0.77 | 487 | 0.61 (air, slab bore face) | 2,720 | 1.293 |

### Loss split: lossless (PEC) shield (raw v2-mesh values)

| case | Q_c, all copper lossy | Q_c, PEC shield (loop + lips only) | implied Q_c of the shield alone | shield share of conductor loss |
|---|---|---|---|---|
| v025_a05_d0635_s1092_o3_v2 | 3,635 | 4,293 | 23,742 | 15 % |

### Field at every probe (kV/mm at 10 uJ)

| case | gap (wall) | gap (lips) | air, slab outer face | air, slab bore face | air, lip-tip triple line | air, lip outer corner | alumina, lip-tip triple line | air, bore triple line | air, slab end face | air, tip/end corner |
|---|---|---|---|---|---|---|---|---|---|---|
| v025_a04_d0635_s1563_o3_v2 | 0.61 | 0.61 | 0.48 | 0.49 | 0.36 | 0.16 | 0.62 | 0.34 | 0.48 | 0.33 |
| v025_a05_d0635_s1092_o3_v2 | 0.71 | 0.71 | 0.56 | 0.56 | 0.41 | 0.12 | 0.72 | 0.40 | 0.56 | 0.36 |
| v025_a06_d0635_s0871_o3_v2 | 0.77 | 0.77 | 0.60 | 0.61 | 0.44 | 0.19 | 0.78 | 0.43 | 0.60 | 0.41 |
| v050_a04_d0635_s0853_o3_v2 | 0.69 | 0.69 | 0.55 | 0.55 | 0.40 | 0.07 | 0.71 | 0.40 | 0.54 | 0.38 |
| v050_a05_d0635_s0551_o3_v2 | 0.80 | 0.80 | 0.63 | 0.64 | 0.46 | 0.09 | 0.82 | 0.46 | 0.63 | 0.44 |
| v050_a06_d0635_s0406_o3_v2 | 0.88 | 0.88 | 0.69 | 0.69 | 0.48 | 0.09 | 0.89 | 0.49 | 0.69 | 0.46 |
| v050_a04_d1000_s1512_o3_v2 | 0.44 | 0.44 | 0.38 | 0.38 | 0.31 | 0.09 | 0.46 | 0.30 | 0.37 | 0.35 |
| v050_a05_d1000_s1042_o3_v2 | 0.51 | 0.51 | 0.44 | 0.44 | 0.36 | 0.10 | 0.53 | 0.34 | 0.43 | 0.37 |
| v050_a06_d1000_s0815_o3_v2 | 0.55 | 0.56 | 0.48 | 0.48 | 0.38 | 0.12 | 0.57 | 0.39 | 0.47 | 0.40 |

### Tuning runs (FE order 2, v1 mesh) and reachability of 145.5 MHz

s_max = R - 3 mm - (r + t), i.e. lips ending 3 mm from the shield wall. s* comes from a linear fit of 1/f0^2 against s through the two runs nearest the target; s0 is the fit intercept (1/f0^2 = B (s + s0)), close to the wall depth t = 3 mm plus a small fringing allowance. Q here is the raw order-2 v1-mesh value (low, see "Convergence").

| d (mm) | V (L) | a | runs, s (mm): f0 (MHz), Q | s_max (mm) | s* (mm) | s0 (mm) | 145.5 MHz reachable? |
|---|---|---|---|---|---|---|---|
| 0.635 | 0.15 | 0.4 | 7.89: 225.8, 2,602; 11.27: 197.6, 2,458 | 11.28 | 23.50 | 3.20 | **no**: f0 = 197.6 MHz at s_max; would need a clearance of -9.2 mm (lips beyond the shield) |
| 0.635 | 0.15 | 0.5 | 5.87: 216.4, 2,711; 8.39: 191.7, 2,602 | 8.40 | 16.98 | 3.29 | **no**: f0 = 191.7 MHz at s_max; would need a clearance of -5.6 mm (lips beyond the shield) |
| 0.635 | 0.15 | 0.6 | 3.86: 225.7, 2,669; 5.51: 203.6, 2,565 | 5.52 | 14.01 | 3.36 | **no**: f0 = 203.6 MHz at s_max; would need a clearance of -5.5 mm (lips beyond the shield) |
| 0.635 | 0.25 | 0.4 | 10.14: 172.9, 2,788; 14.48: 150.1, 2,725 | 14.48 | 15.62 | 3.18 | **no**: f0 = 150.1 MHz at s_max; would need a clearance of 1.9 mm |
| 0.635 | 0.25 | 0.5 | 4.00: 203.1, 3,102; 7.74: 165.2, 2,953; 10.92: 145.5, 2,910; 11.06: 144.8, 2,913 | 11.07 | 10.92 | 3.16 | yes |
| 0.635 | 0.25 | 0.6 | 5.35: 171.3, 2,901; 7.65: 152.4, 2,856 | 7.66 | 8.71 | 3.35 | **no**: f0 = 152.4 MHz at s_max; would need a clearance of 1.9 mm |
| 0.635 | 0.5 | 0.4 | 9.04: 142.5, 3,350; 19.80: 104.0, 3,424 | 19.81 | 8.53 | 3.22 | yes |
| 0.635 | 0.5 | 0.5 | 5.76: 143.5, 3,580; 15.50: 99.6, 3,727 | 15.51 | 5.51 | 3.30 | yes |
| 0.635 | 0.5 | 0.6 | 4.34: 142.8, 3,528; 11.20: 103.8, 3,647 | 11.21 | 4.06 | 3.35 | yes |
| 1 | 0.25 | 0.5 | 7.00: 213.0, 3,150; 11.06: 180.5, 2,985 | 11.07 | 18.85 | 3.37 | **no**: f0 = 180.5 MHz at s_max; would need a clearance of -4.8 mm (lips beyond the shield) |
| 1 | 0.5 | 0.4 | 15.96: 142.3, 3,334; 19.80: 129.8, 3,319 | 19.81 | 15.13 | 3.14 | yes |
| 1 | 0.5 | 0.5 | 10.80: 143.5, 3,565; 15.50: 124.4, 3,572 | 15.51 | 10.42 | 3.35 | yes |
| 1 | 0.5 | 0.6 | 8.56: 143.0, 3,509; 11.20: 129.5, 3,514 | 11.21 | 8.15 | 3.43 | yes |

### Convergence (V = 0.25 L, a = 0.5, d = 0.635 mm, s = 10.92 mm)

Mesh recipes (`lgr_geometry.py`): v1 = 0.45 mm in and around the alumina, 0.25 mm along its edges, 2 mm on loop/lip surfaces, R/5 in the bulk, no extra refinement of the loop rims and lip edges; v2 = v1 with a denser size-field sampling and 0.2 mm along every loop/lip edge (the magnetic field is singular at the square tube rims and lip corners, and that dominates the conductor-loss error). N5 = loop/lip edges 0.1 mm. AMR = Palace adaptive refinement (Doerfler fraction 0.5, 3 iterations) from N1: its error indicator targets the electric field, refines only ~1 % of the mesh per iteration, and leaves Q_c unchanged, so it does not help here.

Reference: `v025_a05_d0635_s1092_mN4_o4`.

| run | FE order | mesh | tets | DOFs | f0 (MHz) | Q_c | Q | E gap (kV/mm) | Q_c vs reference |
|---|---|---|---|---|---|---|---|---|---|
| `v025_a05_d0635_s1092` | 2 | v1 mesh | 80,416 | 596,858 | 145.497 | 3,283 | 2,910 | 0.705 | -10.4 % |
| `v025_a05_d0635_s1092_h070` | 2 | v1, all sizes x0.7 | 221,167 | 1,572,010 | 145.404 | 3,390 | 2,994 | 0.705 | -7.5 % |
| `v025_a05_d0635_s1092_o3` | 3 | v1 mesh | 80,416 | 1,671,534 | 145.373 | 3,380 | 2,986 | 0.705 | -7.8 % |
| `v025_a05_d0635_s1092_h070_o3` | 3 | v1, all sizes x0.7 | 221,167 | 4,451,931 | 145.334 | 3,466 | 3,053 | 0.705 | -5.4 % |
| `v025_a05_d0635_s1092_mG070_o3` | 3 | v1, gap size 0.45 -> 0.32 mm | 140,368 | 2,871,972 | 145.374 | 3,380 | 2,986 | 0.705 | -7.8 % |
| `v025_a05_d0635_s1092_mE050_o3` | 3 | v1, gap-edge size 0.25 -> 0.125 mm | 143,315 | 2,913,483 | 145.376 | 3,374 | 2,982 | 0.705 | -7.9 % |
| `v025_a05_d0635_s1092_mB060_o3` | 3 | v1, bulk size x0.6 | 114,912 | 2,334,789 | 145.360 | 3,439 | 3,032 | 0.705 | -6.2 % |
| `v025_a05_d0635_s1092_mL050_o3` | 3 | v1, loop/lip surface size 2 -> 1 mm | 173,148 | 3,509,313 | 145.313 | 3,503 | 3,081 | 0.705 | -4.4 % |
| `v025_a05_d0635_s1092_mN1_o2` | 2 | N1: denser size-field sampling + loop/lip edges 0.4 mm | 103,765 | 779,940 | 145.347 | 3,288 | 2,914 | 0.705 | -10.3 % |
| `v025_a05_d0635_s1092_mN1_o3` | 3 | N1 | 103,765 | 2,176,803 | 145.305 | 3,430 | 3,025 | 0.705 | -6.4 % |
| `v025_a05_d0635_s1092_mN1_o4` | 4 | N1 | 103,765 | 4,659,988 | 145.295 | 3,511 | 3,088 | 0.705 | -4.2 % |
| `v025_a05_d0635_s1092_mN2_o3` | 3 | N1, loop/lip surface size 2 -> 1 mm | 197,361 | 4,013,922 | 145.299 | 3,467 | 3,053 | 0.705 | -5.4 % |
| `v025_a05_d0635_s1092_mN4_o2` | 2 | N4 = N1 with loop/lip edges 0.2 mm (= v2 recipe) | 218,042 | 1,592,826 | 145.310 | 3,577 | 3,139 | 0.705 | -2.4 % |
| `v025_a05_d0635_s1092_mN4_o3` | 3 | N4 (separate run of the v2 recipe) | 218,042 | 4,477,131 | 145.290 | 3,635 | 3,184 | 0.705 | -0.8 % |
| `v025_a05_d0635_s1092_o3_v2` | 3 | v2 production mesh | 218,042 | 4,477,131 | 145.290 | 3,635 | 3,184 | 0.705 | -0.8 % |
| `v025_a05_d0635_s1092_mN4_o4` | 4 | N4 = v2 mesh | 218,042 | 9,625,532 | 145.288 | 3,665 | 3,206 | 0.705 | +0.0 % |
| `v025_a05_d0635_s1092_mN5_o2` | 2 | N5 = N1 with loop/lip edges 0.1 mm | 989,157 | 6,796,506 | 145.307 | 3,634 | 3,183 | 0.705 | -0.8 % |
| `v025_a05_d0635_s1092_mN1_amr_o3@iteration1` | 3 | AMR, iteration1 | | 2,176,803 | 145.305 | 3,430 | 3,025 | 0.705 | -6.4 % |
| `v025_a05_d0635_s1092_mN1_amr_o3@iteration2` | 3 | AMR, iteration2 | | 2,189,058 | 145.305 | 3,431 | 3,026 | 0.705 | -6.4 % |
| `v025_a05_d0635_s1092_mN1_amr_o3@iteration3` | 3 | AMR, iteration3 | | 2,201,373 | 145.304 | 3,432 | 3,027 | 0.705 | -6.4 % |
| `v025_a05_d0635_s1092_mN1_amr_o3` | 3 | AMR, final | | 2,230,425 | 145.304 | 3,432 | 3,027 | 0.705 | -6.3 % |

### v1 vs v2 mesh, all production cases (FE order 3, raw Q_c)

| case | Q_c v1 | Q_c v2 | v2/v1 | f0 v1 (MHz) | f0 v2 (MHz) |
|---|---|---|---|---|---|
| v025_a04_d0635_s1563_o3 | 3,118 | 3,330 | 1.068 | 145.28 | 145.21 |
| v025_a05_d0635_s1092_o3 | 3,380 | 3,635 | 1.076 | 145.37 | 145.29 |
| v025_a06_d0635_s0871_o3 | 3,282 | 3,519 | 1.072 | 145.33 | 145.25 |
| v050_a04_d0635_s0853_o3 | 4,004 | 4,316 | 1.078 | 145.41 | 145.34 |
| v050_a05_d0635_s0551_o3 | 4,281 | 4,647 | 1.085 | 145.42 | 145.35 |
| v050_a06_d0635_s0406_o3 | 4,212 | 4,551 | 1.081 | 145.40 | 145.34 |
| v050_a04_d1000_s1512_o3 | 3,959 | 4,261 | 1.076 | 145.40 | 145.34 |
| v050_a05_d1000_s1042_o3 | 4,253 | 4,630 | 1.089 | 145.40 | 145.33 |
| v050_a06_d1000_s0815_o3 | 4,192 | 4,534 | 1.081 | 145.40 | 145.33 |

