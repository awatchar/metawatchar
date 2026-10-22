# S2 results: loop-gap MNG slab in a TEM unit-cell waveguide

Palace 0.18.1, frequency-domain (adaptive sweep) and eigenmode solutions. The model is described in [`README.md`](README.md).

- Lattice: a_x = a_y = a_z = 70 mm.
- Cell: loop with r = 21.5 mm, t = 3 mm and length 45 mm along y (the H direction); lips with s = 6.2 mm and d = 1.0 mm alumina.
- Orientation `z`: lip plates parallel to E.
- S-parameters are de-embedded to the slab faces x = 0 and x = L = N a_x, with time convention exp(+jωt).

## Lip-height tuning (single-cell eigenmode, order 2)

| s (mm) | f, PMC x-walls (≈ µ pole f_m) (MHz) | f, PEC x-walls (µ_eff = 0, f_mp) (MHz) | f_mp/f_m | F = 1 − (f_m/f_mp)² |
|---|---|---|---|---|
| 4.0 | 170.888 | 204.366 | 1.196 | 0.301 |
| 6.2 | 149.934 | 179.516 | 1.197 | 0.302 |
| 7.5 | 140.599 | 168.447 | 1.198 | 0.303 |

The fit through s = 4.0 and 7.5 mm is 1/f² ∝ (s + 3.33 mm), the same form and offset as LGR-A in its can. It predicted 149.9 MHz at s = 6.2 mm, and the run gave 149.934 MHz. The filling factor F ≈ 0.30 does not depend on s, so the MNG band spans about f_m … 1.20 f_m (150 … 180 MHz).

## Single-cell eigenmode at s = 6.2 mm (Q)

| x-walls | order | f (MHz) | Q | Q_c (copper) | Q_d (alumina) | p_d |
|---|---|---|---|---|---|---|
| PMC | 2 | 149.934 | 5,016 | 6,198 | 26,300 | 0.951 |
| PEC | 2 | 179.516 | 4,905 | 6,027 | 26,351 | 0.949 |
| PMC | 3 | **149.930** | **5,083** | 6,301 | 26,292 | 0.951 |
| PEC | 3 | **179.503** | 4,982 | 6,145 | 26,343 | 0.949 |

- The plates and walls are ideal, so Q_c is the loop and lip loss only.
- Without the can, Q_c is about 30 % higher than for LGR-A in its 0.5 L can (Q_c 4,630). The can took about 15 % of the conductor loss, and the higher frequency helps further (Q_c ∝ √f).
- Order-2 Q_c is a few per cent low (LGR-A: order 2 → 3 gave +1.6 %, and order 3 → 4 a further +0.8 %).

## Slab S-parameters (orientation z)

| N | L (mm) | run | DOFs | min \|S21\| (dB) at (MHz) | \|S11\| there (dB) | max absorption 1 − \|S11\|² − \|S21\|² at (MHz) | \|S11\| null (MHz) |
|---|---|---|---|---|---|---|---|
| 1 | 70 | N1_z_o3 | 4.15 M | −44.07 at 150.49 | −0.05 | 1.25 % at 150.06 | 137.18 (−74.8 dB) |
| 2 | 140 | N2_z_o3 | 8.08 M | −32.43 at 150.52 | −0.06 | 10.7 % at 149.42 | 136.97 (−69.9 dB), 149.29 (−27.9 dB) |
| 4 | 280 | N4_z_o3 | 15.9 M | −63.02 at 150.41 | −0.06 | 50.6 % at 149.81 | 136.76 (−68.4 dB), 146.07 (−51.9 dB) |

All three are order 3 on the v2 mesh. The order-2 runs of the same meshes (`S_N*_o2.csv`) differ by ≤ 10 kHz in the dip frequency and ≤ 0.09 dB in the dip depth.

|S21| above the dip (dB):

| N | 155 MHz | 160 MHz | 170 MHz | 200 MHz |
|---|---|---|---|---|
| 1 | −1.92 | −0.85 | −0.49 | −0.39 |
| 2 | −6.27 | −3.05 | −1.76 | −1.31 |
| 4 | −17.08 | −9.42 | −5.32 | −3.22 |

Features:

- **Deep stopband just above the pole.** The minimum of |S21| sits at 150.4–150.5 MHz for every N, 0.5 MHz above the eigenmode estimate of the pole (149.93 MHz). Above it, transmission recovers slowly: the slab is only 0.03–0.14 λ thick, so a moderate negative µ attenuates only a few dB. The attenuation grows roughly in proportion to N, as an effective medium would.
- **Reflection null near 137 MHz.** |S11| drops to −68 … −75 dB at 136.8–137.2 MHz for all N. A frequency-independent null of this kind suggests a wave-impedance match (µ_eff ≈ ε_eff), with ε_eff > 1 from the metal and alumina loading. This needs to be confirmed by the retrieval.
- **Absorption is small and narrow.** It peaks at 1.25 % for N = 1, 10.8 % for N = 2 and 51 % for N = 4, in narrow Fabry–Perot-type peaks just below the pole, where |µ_eff| is large. The narrowest peak (N = 4, 149.82 MHz) has an 80 kHz FWHM, which the 10 kHz grid resolves with 8 points.
- **Passivity and reciprocity.** min(1 − |S11|² − |S21|²) > 0 (≥ 5.7e-5) in every run, and S12 = S21 to 1e-11 where both ports were excited.

## Capacitor orientation and electric coupling (N = 1, order 2)

| orientation | gap field | min \|S21\| at (MHz) | depth (dB) | max \|S11 − S22\| | arg S11 − arg S22 at 140 / 145 / 150 MHz |
|---|---|---|---|---|---|
| `z` (production) | along x, ⟂ E | 150.50 | −44.11 | 0 (mirror symmetry in x) | 0° |
| `x` (comparison) | along z, ∥ E | 150.79 | −44.50 | 0.90 (at 150.48 MHz) | −118° / −77° / −56° |

- With the gap field parallel to E, the incident electric field drives the capacitor directly.
- |S11| = |S22| to 0.01 dB, but the reflection phases differ by up to 118°. That is the signature of magneto-electric (omega-type) bianisotropy, which a µ-ε retrieval can't represent.
- The `z` orientation removes it by symmetry, which is why it is used for the µ_eff slabs.

## Mesh and solver convergence (N = 1, orientation z)

| level | mesh | order | DOFs | min \|S21\| at (MHz) | depth (dB) | max absorption | max \|ΔS21\| vs v2/o3 | max \|ΔS11\| vs v2/o3 |
|---|---|---|---|---|---|---|---|---|
| A | v2 (0.2 mm edges) | 2 | 1.48 M | 150.500 | −44.11 | 1.264 % | 1.6e-3 | 1.1e-3 |
| **B (production)** | v2 | 3 | 4.15 M | 150.490 | −44.07 | 1.248 % | — | — |
| C | v3 (0.12 mm edges, all sizes ×0.6–0.7) | 2 | 5.24 M | 150.490 | −43.91 | 1.249 % | 4.0e-4 | 2.9e-4 |

Both refinements change the dip frequency by ≤ 0.007 % and the dip depth by ≤ 0.2 dB. The criteria were < 1 % in frequency and < 1 dB in depth.

The adaptive sweep (6 full solves for N = 1, 9 for N = 4) agrees with 7 uniform full solves at 149.0–152.0 MHz to |ΔS| ≤ 3e-5.

Empty-line check (no cells, 280 mm of line):

- without de-embedding, arg S21 = −43.71° at 130 MHz, against −k0 · 280 mm = −43.72°;
- with Offset = 140 mm on both ports, arg S21 = 3e-7° and |S11| = −178 dB.

## Cost

Everything ran on LANTA's free `compute-limited` partition (CPU billing weight 0.0), so **0.000 SHr** was charged. That is about 710 core-hours in total (≈ 5.6 node-hours, which would have cost ≈ 5.6 SHr on `compute`).

| run | ranks | wall time |
|---|---|---|
| N1_z_o3 | 64 | 60 min |
| N2_z_o3 | 64 | 95 min |
| N4_z_o3 | 128 (2 nodes) | 105 min |
| single-cell eigenmode, order 3 | 64 | 15–21 min |

Almost all of the driven time is the PROM online phase at 7001 output frequencies. Each run needed only 6–9 full solves.
