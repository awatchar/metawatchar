# LGR-A: loop-gap-resonator unit-cell prototype (≈ 145.5 MHz)

This is the first hardware prototype of the MNG unit cell (**the meta-atom**: one cell in a test can, not yet a metamaterial; arrays and μ_eff retrieval follow) for the compact same-band VHF duplexer. It is a shielded copper loop-gap resonator (LGR). Its capacitance comes from a 99.6% alumina slab clamped between two radial lips at the loop slit.

The prototype measures:

1. **Unloaded Q** of a VHF unit cell. Simulation: [`../../sim/lgr_size_sweep`](../../sim/lgr_size_sweep).
2. **Resonant-frequency drift against temperature** in a heated chamber. This is the baseline for thermal self-compensation.

> **Status: dimensions fixed from the Palace size sweep** ([`../../sim/lgr_size_sweep`](../../sim/lgr_size_sweep)). Volume 0.5 L, a = r/R = 0.5, alumina slab d = 1.0 mm, lip height s = 10.42 mm, lip-tip fillet R 0.4 mm. Simulated f₀ = 145.33 MHz and Q ≈ 3,930 with sharp edges; expect f₀ slightly high after filleting, then tune with the screw.

## 3D assembly viewer

`viewer/lgr-a-viewer.html` is an interactive three.js view of the full assembly, with explode, cutaway, X-ray and a label for every part. To run it locally, serve the folder, for example `python -m http.server` in `viewer/`, and open the page. `lgr-a-model.js` is the GLB model, base64-encoded. Regenerate it with `lgr_a_assembly.py`.

## Files

| File | What it is |
|---|---|
| `lgr_a_assembly.py` | Full assembly, with bought-part look-alikes, exported to `out/LGR-A_full.glb/.step` |
| `lgr_a_cad.py` | Parametric CadQuery model. Regenerate with `python lgr_a_cad.py --volume-l 0.25 --a 0.5 --s <mm> --d <mm>` |
| `out/*.step` | One STEP file per part, plus `LGR-A_assembly.step` |
| `out/params.json` | Dimensions used for the STEP files |

## Parts

| Part | Material | Process | Notes |
|---|---|---|---|
| `can_body` | Copper T2 | CNC turning and milling | Inner bore and inner floor need Ra ≤ 0.8 µm (RF surface). Eight M3 tapped holes in the flange. One M4 tapped hole for the tuning screw (at 0°, facing the lips). Two Ø6.5 mm holes for SMA bulkhead jacks (at 120° and 240°). **No anodizing or coating.** |
| `lid` | Copper T2 | CNC | Contact face flat to ≤ 0.02 mm so it makes RF contact all round. Eight Ø3.4 mm clearance holes. |
| `loop` | Copper T2 | CNC; the slit needs a fine slot cutter or **wire EDM** | The slit width must equal the alumina thickness plus 0.02–0.03 mm. Inner and outer surfaces need Ra ≤ 0.8 µm. Deburr the slit edges with no burrs left inside the slit, because burrs concentrate the field and change C. |
| `clip` ×2 (`out/clip_PEEK.step`) | PEEK | CNC | U-clip at each axial end of the lips; squeezes both lips onto the alumina slab (opening 6.95 mm, i.e. 0.05 mm interference). Prevents the air film that would shift f₀ ≈ 0.7 MHz/µm and add ≈ +50 ppm/K of spurious drift (notes/21). |
| `spacer` ×2 | PTFE (or Rexolite) | CNC turning | Holds the loop on axis. If PTFE isn't available, use another low-loss plastic. **Avoid POM and nylon**, which are lossy. |
| `slab` | 99.6% alumina substrate, 0.635 or 1.0 mm | **Bought, laser-cut to size** (not CNC) | Size is (t + s) × z. Use high-purity ≥ 99.5%. Lower grades are about 10× lossier at VHF. |

## Bought parts

- 2 × SMA female bulkhead jacks (single-hole, 1/4-36 thread). Solder a small wire coupling loop to each. Keep coupling weak (S21 ≤ −30 dB) so the measured Q is close to Qu.
- 8 × M3 screws, stainless or brass, length to suit the 6 mm flange.
- No tuning screw: simulation showed it tunes less than 5 kHz. Order the loop in three lip heights (`out/loop_Cu-T2_s10p42.step`, `_s10p55`, `_s10p7`) and pick the one closest to 145 MHz.

## Ordering at JLCCNC

1. Upload `can_body_Cu-T2.step`, `lid_Cu-T2.step` and `loop_Cu-T2.step`, and choose **Copper (T2)**.
2. Add a 2D note covering:
   - threads: M3 ×8 in the flange, M4 ×1 for the tuning screw;
   - the slit width tolerance (+0.02/+0.03 mm over the slab thickness);
   - Ra ≤ 0.8 µm on the RF surfaces;
   - no anodize, no coating.
3. Order `spacer_PTFE.step` as a PTFE part if PTFE is offered; otherwise source it separately.
4. Order the alumina slab from a ceramic-substrate supplier, laser-cut to (t + s) × z, in the thickness chosen.

## Measuring Qu (VNA)

- Use weak, symmetric loop coupling with |S21| at resonance ≤ −30 dB.
- Measure the 3 dB bandwidth to get Q_L, then Qu = Q_L / (1 − |S21|) for symmetric coupling.
- Calibrate at the SMA reference planes.
- Record f0 and Qu against temperature in the heated chamber. This drift is the baseline for thermal self-compensation.

## License

CAD and code: MIT. Measurement data and figures: CC BY 4.0. See the repository root.
