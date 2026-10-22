# LGR-A driven study: results

Palace 0.18.1. Model and method: [`README.md`](README.md). All numbers are in
[`results/summary.json`](results/summary.json); per-run records are in `results/<round>.json`.
Copper σ = 5.8e7 S/m, alumina εr 9.8 / tan δ 4e-5, PTFE spacers included (εr 2.1, tan δ 2e-4),
brass σ = 1.5e7 S/m (assumed). "Sharp" means the lip edges as in the size sweep; "as built" means
with the R0.4 mm lip-tip fillets drawn in `hardware/lgr_a`.

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
