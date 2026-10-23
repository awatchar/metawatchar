# LGR-A first build: bill of materials and order list

This is the order list for the first LGR-A meta-atom prototype. Dimensions are in mm and come from `out/params.json`. The 3D assembly is in `viewer/lgr-a-viewer.html`.

Simulated target: f₀ ≈ 145.1 MHz with the PEEK clips fitted, and Qᵤ ≈ 2,120 for ideal copper. Real machined copper will give less. See `README.md`.

The build has no tuner, so f₀ is set by the choice of lip height. The loop is ordered in three lip heights, each with its own alumina slab.

## Machined parts (JLCCNC or similar)

| # | Part | File | Material | Qty | Key dimensions and callouts |
|---|---|---|---|---|---|
| 1 | Can body | `out/can_body_Cu-T2.step` | Copper T2 | 1 | Bore Ø86.03 × 86.03 deep, wall 4, floor 4. Flange Ø108.03 × 6 with 8 × M3 tapped holes (6 deep). Ø6.5 SMA holes at 120° and 240°, mid-height. Ø1.2 loop-ground holes 6.27° from each SMA hole towards 180°. Ra ≤ 0.8 µm on the bore and floor. No coating. |
| 2 | Lid | `out/lid_Cu-T2.step` | Copper T2 | 1 | Ø108.03 × 4, 8 × Ø3.4 clearance holes. Contact face flat to ≤ 0.02 mm. Ra ≤ 0.8 µm on the inner face. |
| 3a | Loop, s = 10.42 | `out/loop_Cu-T2_s10p42.step` | Copper T2 | 1 | Tube ID 43.01, wall 3, length 64.52. Two 3 mm lips. **Slit 1.025 (+0.01/−0) by wire EDM.** R0.4 lip-tip fillets. Ra ≤ 0.8 µm. Deburr the slit. |
| 3b | Loop, s = 10.55 (nominal) | `out/loop_Cu-T2_s10p55.step` | Copper T2 | 1 | As 3a; this is the nominal build (≈ 145.1 MHz with clips) |
| 3c | Loop, s = 10.70 | `out/loop_Cu-T2_s10p7.step` | Copper T2 | 1 | As 3a; f₀ about 0.8 MHz lower |
| 4a | Slab clip, opening 6.95 | `out/clip_PEEK_open6p95.step` | PEEK | 4 | U-clip 8 × 10.95 × 4, 6 deep. 0.05 mm interference on the lip pack. Two per build, plus spares. |
| 4b | Slab clip, opening 7.00 | `out/clip_PEEK_open7p0.step` | PEEK | 4 | Zero interference (fallback if 6.95 does not fit) |
| 4c | Slab clip, opening 7.05 | `out/clip_PEEK_open7p05.step` | PEEK | 2 | Loose fit, used for the "unclamped" control measurement |
| 5 | Spacer ring | `out/spacer_PTFE.step` | PTFE (or Rexolite) | 4 | ID 43.01, OD 49.01, height 10.75. Two per build, plus spares. Not POM or nylon. |

## Bought parts

| # | Part | Qty | Specification |
|---|---|---|---|
| 6a | Alumina slab, for loop 3a | 2 | 99.6 % Al₂O₃, 1.00 ± 0.01 mm thick, 13.42 × 64.52, laser-cut, as-fired surface |
| 6b | Alumina slab, for loop 3b | 2 | As 6a, 13.55 × 64.52 |
| 6c | Alumina slab, for loop 3c | 2 | As 6a, 13.70 × 64.52 |
| 7 | SMA female bulkhead jack | 4 | Single-hole panel mount, 1/4-36 thread, fits a Ø6.5 hole in a 4 mm wall, with a solder-cup or extended centre pin. Two for the build, two spare. |
| 8 | Copper wire | 0.5 m | Solid Cu, Ø1.0 mm. Two coupling loops of 4.7 × 4.7 mm (centreline), as drawn in `README.md`. |
| 9 | M3 screws | 8 (+4) | M3 × 8, brass or A2 stainless |

Each slab's radial length equals t + s (loop wall 3 mm plus the lip height), so the slab and loop variants must be paired: 6a with 3a, 6b with 3b, 6c with 3c.

## Possible addition: cells for the slab fixture (later paper)

To get every meta-atom from one machining batch, the unshielded array cells (lip height s = 6.2 mm, 8 cells for a 4 × 2 slab) could go into the same order. They are on hold until the full-wave fixture simulation has fixed the cell count and fixture type.

## Notes for the order

- Order all three loop variants at once. With no tuner, a second order would cost weeks.
- The slit width (1.025 mm for a 1.00 mm slab) is the tightest tolerance in the build. Ask for EDM and for a measured slit width on each loop.
- If the alumina supplier's thickness tolerance is wider than ±0.01 mm, ask for the measured thickness of each slab. Each 1 µm of air film moves f₀ by ≈ 0.7 MHz.
- Copper T2 (C11000) without plating is what the simulation assumes. Silver plating raises Q but changes the baseline; leave it for a later build.

## Before ordering

1. Review the 3D viewer (`viewer/lgr-a-viewer.html`): part placement, the loop orientation (flat in the mid-height plane) and the clip fit.
2. Check JLCCNC's minimum wall and slit capabilities, and get their quote for EDM on copper.
3. Confirm the SMA jack's panel-thickness range covers the 4 mm wall.
