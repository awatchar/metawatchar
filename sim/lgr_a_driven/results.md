# LGR-A driven study: results

Palace 0.18.1. Model and method: [`README.md`](README.md). All numbers are in
[`results/summary.json`](results/summary.json); per-run records are in `results/<round>.json`.
Copper σ = 5.8e7 S/m, alumina εr 9.8 / tan δ 4e-5, PTFE spacers included (εr 2.1, tan δ 2e-4),
brass σ = 1.5e7 S/m (assumed). "Sharp" means the lip edges as in the size sweep; "as built" means
with the R0.4 mm lip-tip fillets drawn in `hardware/lgr_a`.

> **Correction (copper surface impedance).** Palace's `Impedance` boundary combines Rs and Ls in
> parallel (1/Z = 1/Rs + 1/(jωLs)). The first runs of this study (sections 1–6 below) entered
> Rs and Ls = Rs/ω, which is Z = Rs(1 + j)/2: **half the copper loss**, so every Q, Q_u and Q_c
> in sections 1–6 is too high (Q_c by 2×, Q by about 1.8×) and the loop sizes chosen there give
> too weak a coupling. Frequencies, Q_ext and field shapes are not affected (f0 moves by
> −16 kHz). The section **"Corrected copper input"** directly below replaces them; it uses
> Rs_p = 2 Rs and Ls_p = 2 Rs/ω (Z = Rs(1 + j), Rs at each case's f0), checked against the
> analytic copper Q of the empty can. Run names end in `_cu2`; numbers are in
> [`results/summary_cu2.json`](results/summary_cu2.json).

## Corrected copper input (`_cu2` runs)

Geometry of the build ([`../../hardware/lgr_a`](../../hardware/lgr_a)): V = 0.5 L, a = 0.5,
1.0 mm alumina, R0.4 mm lip-tip fillets, PTFE spacers, lip height s = 10.42 / 10.55 / 10.70 mm
(the three loop variants; s = 10.55 mm is the default).

### Check of the impedance input: empty copper can

| Mode | f analytic (GHz) | f Palace (GHz) | Q_c analytic, R H/(δ (R + H)) | Q Palace (order 2) | Ratio |
|---|---|---|---|---|---|
| TM010 | 2.667643 | 2.667593 | 22,411 | 22,410 | 0.99995 |

(`e_emptycan_cu2_o2`, Rs evaluated at the TM010 frequency.)

### Eigenmodes without coupling loops (order 3 unless noted)

Q_d = 1/(Σ p_i tan δ_i) over alumina, PTFE (and PEEK); Q_c = (1/Q − 1/Q_d)⁻¹ (all copper).

| Case | f0 (MHz) | **Q_u** | Q_c | Q_d | First runs: f0 / Q / Q_c |
|---|---|---|---|---|---|
| Sharp lip edges, s = 10.42 mm | 145.1960 | **2,151** | 2,349 | 25,509 | 145.2117 / 3,931 / 4,647 |
| R0.4 fillets, s = 10.42 mm | 145.8930 | **2,157** | 2,356 | 25,565 | 145.9088 / 3,934 / 4,650 |
| **R0.4 fillets, s = 10.55 mm (default build)** | **145.2048** | **2,152** | 2,350 | 25,561 | — |
| R0.4 fillets, s = 10.70 mm | 144.4247 | **2,146** | 2,342 | 25,557 | — |
| R0.4 fillets, s = 10.55 mm, **order 4** | 145.2043 | **2,165** | 2,365 | 25,561 | — |
| R0.4 fillets, s = 10.55 mm, **with PEEK slab clips** | 145.0880 | **2,124** | 2,350 | 22,079 | — |

- **Q_u of the build is about 2,150 (order 3), not 3,930.** Q_c is 0.506 of the first runs' value
  and Q is 0.547–0.548 of it (the alumina and PTFE loss, Q_d ≈ 25,500, is unchanged).
- Order 3 → 4: Q +0.61 %, Q_c +0.67 %, f0 −0.0003 %. With the size-sweep convergence study
  (order 3 about 1–2 % low), the best estimate for smooth copper with σ = 5.8e7 S/m is
  **Q_u ≈ 2,170–2,190**.
- Lip-height trim with fillets: df0/ds = −5.24 MHz/mm (10.42 → 10.70 mm); 145.200 MHz needs
  s = 10.551 mm (1/f² linear in s). Q_u changes by < 0.6 % over the three variants.
- The doubled surface reactance lowers f0 by 16 kHz against the first runs (−0.011 %).

### PEEK slab clips (s = 10.55 mm, R0.4 fillets, order 3)

Two U-clips as in `hardware/lgr_a/lgr_a_assembly.py` (`slab_clip`): they cover the outer 6 mm of
both lips over 4 mm along z at each loop end, with 2 mm walls; PEEK εr = 3.2, tan δ = 3e-3
(assumed typical values; to be checked against the supplier's VHF data). The 0.05 mm interference
is not modelled.

| | f0 (MHz) | Q_u | Q_c | p (PEEK) | Q_d,PEEK = 1/(p tan δ) |
|---|---|---|---|---|---|
| Without clips | 145.2048 | 2,152 | 2,350 | — | — |
| With clips | 145.0880 | 2,124 | 2,350 | 0.208 % | 160,500 |
| Change | **−117 kHz (−0.080 %)** | **−1.29 %** | +0.02 % | | |

The clips sit in the fringing field at the mouth of the gap; they add capacitance (f0 down by
about 1.7 linewidths) and a small dielectric loss. Q_u ∝ 1/tan δ_PEEK for this term: at
tan δ = 5e-3 the Q_u cost would be about 2.1 %.

### Coupling loops for |S21(f0)| ≈ −30 dB (s = 10.55 mm, R0.4 fillets, PTFE spacers)

Q_ext depends only on the loop (it is the same as in the first runs to < 0.1 %), but the lower
Q_u weakens the coupling β = Q_u/Q_ext. Eigenmode with 50 Ω ports, order 2 (Q_ext), and
|S21(f0)| = 2β/(1 + 2β) with the order-3 Q_u (2,152):

| Loop, depth × width (mm, wire centreline) | Area (mm²) | Q_ext (each port) | β | \|S21(f0)\| (dB) |
|---|---|---|---|---|
| 4.1 × 4.1 | 16.8 | 231,460 | 0.0093 | −34.8 |
| 4.5 × 4.5 | 20.3 | 156,530 | 0.0137 | −31.5 |
| **4.7 × 4.7** | **22.1** | **130,560** | **0.0165** | **−29.9** |
| 4.9 × 4.9 | 24.0 | 109,830 | 0.0196 | −28.5 |

Fit: 1/Q_ext = 1.18e-8 · A^2.09 (A in mm²); −30.0 dB needs A = 22.0 mm² (4.69 × 4.69 mm). One
side +0.1 mm changes |S21| by about +0.38 dB.

**Recommended loop: 4.7 × 4.7 mm** (1.0 mm Cu wire, mid-height plane, radial leg from the inner
wall to 4.7 mm inside it, tangential leg 4.7 mm towards 180°, return leg to the wall). The ground
point on the inner wall is 4.7 mm (tangentially) from the port axis, i.e. asin(4.7/R) = 6.27°
from the SMA hole (5.3° was for a 4.0 mm loop); the return leg is 4.44 mm long because the wall
curves towards it (sagitta 0.26 mm), so the centreline inside the can is ≈ 13.8 mm.

**Prediction for the VNA measurement** (4.7 × 4.7 mm loops, s = 10.55 mm, R0.4 fillets, PTFE
spacers, no clips; driven, order 3, 1 kHz steps over f0 ± 1 MHz;
[`results/sparams_d_fil_s1055_xy4p7x4p7_cu2_o3.csv`](results/sparams_d_fil_s1055_xy4p7x4p7_cu2_o3.csv)):

| f0 (MHz) | \|S21(f0)\| (dB) | \|S11(f0)\| = \|S22(f0)\| (dB) | Q_L (3 dB) | 3 dB bandwidth (kHz) | β1 = β2 | Q_u (circle fit) | S21 off resonance (dB) |
|---|---|---|---|---|---|---|---|
| 145.2043 | −29.93 | −0.288 | 2,081 | 69.8 | 0.01647 | 2,150 | −86 |

The eigenmode solution with resistive ports on the same mesh gives f0 = 145.2042 MHz,
Q_L = 2,082, β = 0.01647, |S21(f0)| = −29.93 dB and Q_u = 2,150. The loops cost < 0.1 % of
Q_u. With the PEEK clips fitted, expect f0 ≈ 145.09 MHz and Q_u ≈ 2,124 (and |S21| about
0.1 dB weaker).

### Mode fields of the default build (field maps, paper P1 Fig. 2)

Case `e_fil_s1055_cu2_o2_fields` (group `cu2_fields`): the default build (s = 10.55 mm, R0.4
lip fillets, PTFE spacers, corrected copper input), no coupling loops (they cost < 0.1 % of Q_u
and do not change the field picture), FE order 2, fields of the fundamental mode saved. Field
maps and numbers: `fields_fig2.py`, run on the cluster after the solve; only the figures and
[`results/fields_e_fil_s1055_cu2_o2_fields.json`](results/fields_e_fil_s1055_cu2_o2_fields.json)
are kept (the ParaView dump, about 360 MB, is not).

| Run | f0 (MHz) | Q_u | Q_c | Q_d |
|---|---|---|---|---|
| This run, order 2 | 145.2126 | 2,121 | 2,313 | 25,563 |
| `e_fil_s1055_cu2_o3` (order 3) | 145.2048 | 2,152 | 2,350 | 25,561 |
| Difference | +54 ppm | −1.4 % | −1.6 % | |

The order-2 values match the order-2 loop runs (Q_u 2,120, f0 145.2128 MHz with 4.7 × 4.7 mm
loops); the order-3 values remain the reference for Q_u.

Figures in [`figures/`](figures) (PNG at 600 dpi and PDF):

| File | Content |
|---|---|
| `fig2a_E_midplane` | \|E\| in the mid-height plane z = 0: whole section (log scale) and the gap (linear scale) |
| `fig2b_Hz_midplane` | H_z in z = 0 (diverging scale): axial bore flux and the opposite return flux in the annulus; H_z along the radius at φ = 180° |
| `fig2c_H_vertical` | \|H\| and the in-plane field direction in the vertical plane y = 0, through the gap |
| `fig2_fields` | All five panels in one figure (double-column width) |

Fields are peak phasor amplitudes, scaled to 1 W dissipated in the resonator (stored energy
W = Q P/ω = 2.33 µJ). They scale as √P: multiply by 1.414 for 2 W. Copper is drawn from the
analytic geometry, since conductors are boundaries and are not meshed. Colormaps: magma (\|E\|),
viridis (\|H\|) and cmocean balance (H_z), all perceptually uniform.

| Quantity at 1 W dissipated | Value |
|---|---|
| \|E\| in the alumina, mid-lip (x = 29.0 mm) | 245 kV/m, uniform across the 1.0 mm slab |
| Gap voltage at mid-lip | 245 V (347 V at 2 W) |
| Largest \|E\| in the air, at the lip-tip fillet next to the slab (x = 34.77 mm, y = 0.51 mm) | 1.57 MV/m (indicative: this is the air wedge that the fillet opens at the alumina face; the value depends on the mesh there) |
| Largest \|E\| in the alumina (same corner) | 383 kV/m |
| H_z in the bore (mean) / centre | 146.9 / 146.2 A/m |
| H_z in the annulus outside the loop (mean) | −55.4 A/m |
| Return flux / bore flux through z = 0 | 1.002 |

Checks: the electric and magnetic energies integrated from the field dump agree to 0.26 %, and
the alumina participation from the field dump (0.9629) matches Palace's own value (0.9632).
The electric energy is 96 % in the alumina, so the field picture is that of a lumped LC: the
capacitance is the alumina-filled slit, and the inductance is the bore (uniform axial H), whose
flux returns through the annulus at about 0.38 of the bore field, over 2.66 × the bore area.

Cost: one 32-core job on a paid queue (the free queue was estimated at about 22 h), 7.6 min
including the figure step: 0.032 SHr. Three later figure-only re-runs (4 cores, < 30 s each) cost
< 0.001 SHr each.

## Summary

- **Coupling-loop orientation matters by more than 100 dB.** The return flux between the LGR and
  the can is axial, so the loop must lie in the mid-height plane (normal along the can axis). A
  5 × 5 mm loop in that plane gives |S21(f0)| = −22.9 dB; the same loop in the radial–axial plane
  (normal azimuthal) gives −126.6 dB, i.e. no usable signal.
- **Recommended loop: 4.0 mm deep × 4.0 mm wide (wire centreline), 1.0 mm copper wire,** in the
  mid-height plane, tangential leg pointing away from the gap, grounded on the can wall 4.0 mm from
  the port axis. Order-3 result (sharp edges): **|S21(f0)| = −30.57 dB, β1 = β2 = 0.0153,
  Q_L = 3,808, Q_u = 3,925, f0 = 145.214 MHz** (3 dB bandwidth 38.1 kHz).
  - Coupling scales as β ∝ A^2.13 (A = loop area), so |S21| changes by 0.18 dB per 1 % of area,
    or about 0.46 dB per 0.1 mm on one side. −30.0 dB exactly needs A = 16.7 mm² (4.1 × 4.1 mm).
  - The loops barely load the resonator apart from the intended coupling: f0 moves by +2 kHz and
    the wire loss lowers Q_u by < 0.1 % (4 × 4 mm).
- **Prediction for the measurement (as built: 4 × 4 mm loops and R0.4 mm lip fillets, order 3):
  f0 = 145.906 MHz, |S21(f0)| = −30.52 dB, |S11(f0)| = |S22(f0)| = −0.265 dB, Q_L = 3,812,
  β1 = β2 = 0.0153, Q_u = 3,929** (`results/sparams_d_xy4x4_fil_o3.csv`). The symmetric
  reduction Q_u = Q_L/(1 − |S21(f0)|) gives the same Q_u. Ideal copper (σ = 5.8e7 S/m, smooth)
  and an ideal gap are assumed; see "Not modelled".
- **Convergence (4 × 4 mm loop):** order 2 → 3 changes |S21(f0)| by +0.11 dB and f0 by −0.012 %;
  an h-refined mesh (all sizes × 0.8, loop/port region × 0.48) changes them by +0.07 dB and
  −0.007 %. Both meet the criteria (0.5 dB, 0.05 %). Q_u changes by +1.3 % (order 3) and +0.8 %
  (h-refined): as in the size sweep, conductor Q converges more slowly than f0 and S21.
- **Driven and eigenmode agree.** The resonance fit of the driven S-parameters gives the same β,
  Q_L and Q_u as the eigenmode solution with resistive ports (to < 0.1 %). The adaptive
  reduced-order sweep matches full-order solves at five frequencies across the resonance to
  6e-7 relative error in S21.
- **M4 tuning screw at 0°: no useful tuning range.** Penetration 2 / 4 / 6 mm moves f0 by
  +2.4 / +4.4 / +2.4 kHz (max ≈ 0.1 linewidth) and lowers Q by ≤ 0.13 %. The screw sits on the
  mirror plane of the gap, where the potential of the two lips (±V/2) cancels, so it hardly
  couples to the electric field; the small upward shift is the magnetic (inductive) perturbation.
- **Exploratory axial tuner:** a 6 mm brass rod through the lid on the axis, reaching into the
  loop bore, raises f0 by +100 kHz (20 mm deep, tip 9 mm inside the bore) and +211 kHz (30 mm
  deep), at a Q cost of 1.3 % and 2.5 % (order 2).
- **R0.4 mm lip-tip fillets: f0 +0.697 MHz (+0.48 %), Q +0.09 %** (order 3, no loops). The fillets
  on the slit-side edges open an air wedge between the copper and the alumina at the lip tips,
  which removes gap capacitance. The shift is 18 linewidths and 160 times the screw range.
- **PTFE spacers:** f0 −122 kHz (−0.08 %); 0.30 % of the electric energy is in the PTFE, so their
  dielectric loss limits Q only at 1.7 × 10⁶. Q changes by +0.5 % (within mesh noise).

## 1. Loop sizing (eigenmode with 50 Ω resistive ports, FE order 2, sharp edges)

Q_u = (1/Q_L − 1/Q_ext,1 − 1/Q_ext,2)⁻¹, β_i = Q_u/Q_ext,i,
|S21(f0)| = 2√(β1β2)/(1 + β1 + β2). Loop depth × width are wire-centreline dimensions.

| Loop (mm) | Area (mm²) | f0 (MHz) | Q_L | Q_u | Q_ext,1 | Q_ext,2 | β1 | β2 | \|S21(f0)\| (dB) | driven \|S21(f0)\| (dB) |
|---|---|---|---|---|---|---|---|---|---|---|
| none | — | 145.2286 | 3,879 | 3,879 | — | — | — | — | — | — |
| 3 × 3 | 9 | 145.2299 | 3,844 | 3,877 | 891,500 | 891,900 | 0.00435 | 0.00435 | −41.29 | −41.29 |
| **4 × 4** | **16** | **145.2311** | **3,762** | **3,876** | **257,300** | **257,300** | **0.01507** | **0.01506** | **−30.68** | **−30.68** |
| 5 × 5 | 25 | 145.2319 | 3,600 | 3,876 | 101,100 | 101,100 | 0.03833 | 0.03833 | −22.95 | −22.95 |
| 8 × 6 | 48 | 145.2319 | 2,998 | 3,873 | 26,510 | 26,510 | 0.1461 | 0.1461 | −12.91 | −12.91 |
| 10 × 8 | 80 | 145.2363 | 2,150 | 3,871 | 9,668 | 9,668 | 0.4004 | 0.4004 | −7.04 | −7.04 |
| 5 × 5, radial–axial plane | 25 | 145.2311 | 3,877 | 3,877 | 1.68e10 | 1.63e10 | 2.3e-7 | 2.4e-7 | −126.6 | — |

Fit through 3 × 3 and 5 × 5: β = 4.03e-5 · A^2.130 (A in mm²). The 4 × 4 mm loop was then
simulated directly (β = 0.0151 vs 0.0148 predicted by the fit).

## 2. S-parameters (driven, two excitations, adaptive fast sweep)

Each CSV has columns `f_Hz, S11_re, S11_im, S21_re, S21_im, S22_re, S22_im, S12_re, S12_im`,
1 kHz steps over f0 ± 0.6 MHz (the `_wide` file: 144.0–145.8 MHz in 2 kHz steps). Reference
planes: the 50 Ω port 1 mm behind the inner can wall (the SMA connector adds only phase).

Resonance fit: S21 = C/(1 + 2jQ_L(f − f0)/f0) + background; S11, S22 = A + B/(1 + 2jQ_L(f − f0)/f0);
β_i = d_i/(2 − d1 − d2) with d_i = |B_i| (circle diameter); Q_u = Q_L(1 + β1 + β2).

| Case | Order | f0 (MHz) | \|S21(f0)\| (dB) | Q_L (3 dB) | β1 | β2 | Q_u | \|S11(f0)\| (dB) | S21 off resonance (dB) | CSV |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 × 3 | 2 | 145.2299 | −41.29 | 3,843 | 0.00435 | 0.00435 | 3,877 | −0.077 | −103 | `sparams_d_xy3x3_o2.csv` |
| 4 × 4 | 2 | 145.2311 | −30.68 | 3,761 | 0.01507 | 0.01506 | 3,875 | −0.261 | −92 | `sparams_d_xy4x4_o2.csv` |
| 4 × 4, h-refined | 2 | 145.2203 | −30.61 | 3,790 | 0.01519 | 0.01519 | 3,905 | −0.263 | −92 | `sparams_d_xy4x4_r_o2.csv` |
| **4 × 4** | **3** | **145.2136** | **−30.57** | **3,808** | **0.01526** | **0.01526** | **3,925** | **−0.264** | **−92** | **`sparams_d_xy4x4_o3.csv`** |
| 5 × 5 | 2 | 145.2319 | −22.95 | 3,599 | 0.03834 | 0.03834 | 3,875 | −0.645 | −84 | `sparams_d_xy5x5_o2.csv` |
| 8 × 6 | 2 | 145.2319 | −12.91 | 2,997 | 0.1461 | 0.1461 | 3,873 | −2.23 | −72 | `sparams_d_xy8x6_o2.csv` |
| 10 × 8 | 2 | 145.2363 | −7.04 | 2,149 | 0.4005 | 0.4005 | 3,871 | −5.12 | −63 | `sparams_d_xy10x8_o2.csv` |
| 4 × 4, as built (R0.4 fillets) | 2 | 145.9135 | −30.63 | 3,764 | 0.01515 | 0.01515 | 3,878 | −0.262 | −92 | `sparams_d_xy4x4_fil_o2.csv` |
| **4 × 4, as built (R0.4 fillets)** | **3** | **145.9056** | **−30.52** | **3,812** | **0.01535** | **0.01535** | **3,929** | **−0.265** | **−92** | **`sparams_d_xy4x4_fil_o3.csv`** |

The "S21 off resonance" column is the fitted non-resonant (direct port-to-port) term; it is
60 dB below the resonant peak for every loop size, so it does not affect the Q extraction.

Reduced-order check (`d_xy4x4_o2_full`): full-order solves at f0 − 100, −19, 0, +19, +100 kHz
agree with the adaptive sweep to ≤ 6.3e-7 relative error in S21.

## 3. Convergence (4 × 4 mm loop)

| Run | f0 (MHz) | \|S21(f0)\| (dB) | Q_L | Q_u | Δ\|S21\| vs order 2 (dB) | Δf0 vs order 2 |
|---|---|---|---|---|---|---|
| Order 2, production mesh (273k tets) | 145.2311 | −30.68 | 3,761 | 3,875 | — | — |
| Order 2, h-refined mesh (593k tets) | 145.2203 | −30.61 | 3,790 | 3,905 | +0.07 | −0.007 % |
| Order 3, production mesh | 145.2136 | −30.57 | 3,808 | 3,925 | +0.11 | −0.012 % |

Eigenmode on the same meshes gives the same numbers (|S21| −30.68 / −30.61 / −30.57 dB). Both
refinements meet the criteria (|Δ\|S21\|| < 0.5 dB, |Δf0| < 0.05 %). Q_u still rises by 1.3 %
from order 2 to 3; from the size-sweep convergence study the order-3 value is itself 1–2 % low,
so the best estimate for the sharp-edged cell is Q_u ≈ 3,960–4,000.

## 4. M4 tuning screw (eigenmode, order 3, no loops, sharp edges)

Brass Ø4.0 mm cylinder at 0°, mid-height, facing the lip tips (8.09 mm from the inner wall).

| Penetration (mm) | Gap to lip tips (mm) | f0 (MHz) | Δf0 (kHz) | Q | ΔQ |
|---|---|---|---|---|---|
| 0 (plain wall) | 8.09 | 145.2117 | 0 | 3,931 | — |
| 2 | 6.09 | 145.2141 | +2.4 | 3,929 | −0.04 % |
| 4 | 4.09 | 145.2161 | +4.4 | 3,928 | −0.06 % |
| 6 | 2.09 | 145.2141 | +2.4 | 3,926 | −0.13 % |

The tuning range is about +4 kHz, 0.1 of the 3 dB linewidth: the screw is not a usable tuner. It
lies on the mirror plane between the lips, where the gap potential is zero.

Exploratory alternative (order 2, reference f0 = 145.2286 MHz, Q = 3,879): a brass Ø6 mm rod
through the lid on the can axis.

| Rod depth below the lid (mm) | Tip inside the bore (mm) | Δf0 (kHz) | ΔQ |
|---|---|---|---|
| 20 | 9.2 | +100 | −1.3 % |
| 30 | 19.2 | +211 | −2.5 % |

## 5. Lip-tip fillets (eigenmode, order 3, no loops)

| Lip-tip edges | f0 (MHz) | Q | Q_c | p (alumina) |
|---|---|---|---|---|
| Sharp | 145.2117 | 3,931 | 4,647 | 0.9653 |
| R0.4 mm (4 axial tip edges) | 145.9088 | 3,934 | 4,650 | 0.9630 |
| Change | +0.697 MHz (+0.48 %) | +0.09 % | +0.07 % | |

Q is unchanged within the mesh error. The f0 shift is consistent with losing about 0.13 mm of
effective lip height: from the size-sweep fit (C ∝ s + s0, s0 ≈ 3.35 mm) a lip height of about
10.55 mm instead of 10.42 mm restores the sharp-edge f0 (estimate, not simulated).

## 6. PTFE spacers (eigenmode, order 2, no loops)

| Spacers | f0 (MHz) | Q | p (PTFE) |
|---|---|---|---|
| Omitted | 145.3510 | 3,861 | — |
| Included | 145.2286 | 3,879 | 0.0030 |

## Not modelled (affects the comparison with measurement)

- **Air films in the gap.** The alumina is assumed to fill the slit and touch both lips. Any air
  film between the slab and the copper (total thickness δ, summed over both faces) is in series
  with the 1.0 mm alumina and scales the gap capacitance by 1/(1 + εr δ/d). Estimate (analytic,
  not simulated): f0 rises by about **0.71 MHz per µm** of total film, and by about 11.6 % (to
  about 162 MHz) for the 25 µm clearance of a 1.025 mm slit with a loose 1.0 mm slab. This is far
  beyond the range of any tuner studied here, so the slab-to-lip interface and the slab
  thickness set f0 in hardware.
- **Lip height tolerance:** df0/ds ≈ −5.3 MHz/mm (from the size-sweep fit), i.e. −53 kHz per
  10 µm.
- **Real copper:** machined-surface conductivity and roughness lower Q_c (a 5–15 % effect was
  estimated in the size sweep; to be verified against the measurement).
- **SMA connector:** modelled as a 50 Ω air coax ending 1 mm behind the inner wall; the real
  jack's PTFE insulator and pin protrusion change the loop's self-inductance slightly and add a
  phase shift. Coupling is set mostly by the loop area.

## Cost

All runs used a free partition: 0.000 SHr. Wall time on 32 cores per run: driven order 2
4.5–5.5 min, driven order 3 19–24 min, eigenmode order 2 7–11 min (32 min on the 593k-tet mesh), eigenmode
order 3 36–46 min.

Corrected runs (`_cu2`): also on a free partition, 0.000 SHr. Wall time on 32 cores (shared
nodes): eigenmode order 3 37–70 min, driven order 3 (2001 points) 29 min, eigenmode order 2 with
ports 12–15 min; order 4 on 128 cores 38 min.
