# LGR-A: loop-gap-resonator unit-cell prototype (≈ 145.2 MHz)

This is the first hardware prototype of the MNG unit cell (**the meta-atom**: one cell in a test can, not yet a metamaterial; arrays and μ_eff retrieval follow) for the compact same-band VHF duplexer. It is a shielded copper loop-gap resonator (LGR). Its capacitance comes from a 99.6% alumina slab clamped between two radial lips at the loop slit.

The prototype measures:

1. **Unloaded Q** of a VHF unit cell. Simulation: [`../../sim/lgr_size_sweep`](../../sim/lgr_size_sweep).
2. **Resonant-frequency drift against temperature** in a heated chamber. This is the baseline for thermal self-compensation.

> **Status: dimensions fixed from Palace** ([`../../sim/lgr_size_sweep`](../../sim/lgr_size_sweep), [`../../sim/lgr_a_driven`](../../sim/lgr_a_driven)). Volume 0.5 L, a = r/R = 0.5, alumina slab d = 1.0 mm, lip height s = 10.55 mm, lip-tip fillet R 0.4 mm, coupling loops 4.7 × 4.7 mm.
>
> Simulated with the corrected copper surface impedance (see the correction note in [`../../sim/lgr_a_driven/results.md`](../../sim/lgr_a_driven/results.md)):
>
> | Quantity | No clips | With PEEK clips |
> |---|---|---|
> | f₀ | 145.205 MHz | ≈ 145.09 MHz |
> | Qᵤ (ideal copper, σ = 5.8 × 10⁷ S/m) | 2,152 | 2,124 |
> | \|S21(f₀)\| with the 4.7 mm loops | −29.9 dB | ≈ −30 dB |
>
> f₀ falls 5.24 MHz per mm of lip height. Real machined copper will give a lower Qᵤ than ideal copper (roughness, lower conductivity), so treat 2,150 as an upper bound.
>
> Earlier versions of this page quoted Qᵤ ≈ 3,930. That value came from a surface-impedance input that halved the copper loss and is superseded.

## 3D assembly viewer

`viewer/lgr-a-viewer.html` is an interactive three.js view of the full assembly, with explode, cutaway, X-ray and a label for every part. To run it locally, serve the folder, for example `python -m http.server` in `viewer/`, and open the page. `lgr-a-model.js` is the GLB model, base64-encoded. Regenerate it with `lgr_a_assembly.py`.

## Files

| File | What it is |
|---|---|
| `lgr_a_assembly.py` | Full assembly, with bought-part look-alikes, exported to `out/LGR-A_full.glb/.step` |
| `lgr_a_cad.py` | Parametric CadQuery model. Regenerate with `python lgr_a_cad.py --volume-l 0.25 --a 0.5 --s <mm> --d <mm>` |
| `out/*.step` | One STEP file per part, plus `LGR-A_assembly.step` |
| `out/params.json` | Dimensions used for the STEP files |
| `BOM.md` | Order list for the first build: quantities, variants, tolerances |

## Parts

| Part | Material | Process | Notes |
|---|---|---|---|
| `can_body` | Copper T2 | CNC turning and milling | Inner bore and inner floor need Ra ≤ 0.8 µm (RF surface). Eight M3 tapped holes in the flange. Two Ø6.5 mm holes for SMA bulkhead jacks (at 120° and 240°), plus two Ø1.2 mm holes for the coupling-loop ground ends, 6.27° from each SMA hole towards 180°. No tuning-screw hole. **No anodizing or coating.** |
| `lid` | Copper T2 | CNC | Contact face flat to ≤ 0.02 mm so it makes RF contact all round. Eight Ø3.4 mm clearance holes. |
| `loop` | Copper T2 | CNC; the slit needs a fine slot cutter or **wire EDM** | The slit width must equal the alumina thickness plus 0.02–0.03 mm. Inner and outer surfaces need Ra ≤ 0.8 µm. Deburr the slit edges with no burrs left inside the slit, because burrs concentrate the field and change C. |
| `clip` ×2 (`out/clip_PEEK.step`) | PEEK | CNC | U-clip at each axial end of the lips; squeezes both lips onto the alumina slab (opening 6.95 mm, i.e. 0.05 mm interference; also order 7.00 and 7.05 mm openings in case the fit is too tight). Simulated cost: f₀ −0.12 MHz, Qᵤ −1.3 %. Prevents the air film that would shift f₀ ≈ 0.7 MHz/µm and add ≈ +50 ppm/K of spurious drift (sim/lgr_a_driven). |
| `spacer` ×2 | PTFE (or Rexolite) | CNC turning | Holds the loop on axis. If PTFE isn't available, use another low-loss plastic. **Avoid POM and nylon**, which are lossy. |
| `slab` | 99.6% alumina substrate, 0.635 or 1.0 mm | **Bought, laser-cut to size** (not CNC) | Size is (t + s) × z. Use high-purity ≥ 99.5%. Lower grades are about 10× lossier at VHF. |

## Bought parts

- 2 × SMA female bulkhead jacks (single-hole, 1/4-36 thread). Solder a 1.0 mm copper-wire coupling loop to each centre pin: 4.7 mm radially in from the inner wall, 4.7 mm towards 180° in the mid-height plane, then parallel to the first leg back to the wall, soldered at the Ø1.2 mm hole. Simulated |S21(f₀)| = −29.9 dB (β ≈ 0.017 per port). 4.5 mm gives −31.5 dB and 4.9 mm −28.5 dB, so the loop size can be trimmed after the first measurement.
- 8 × M3 screws, stainless or brass, length to suit the 6 mm flange.
- No tuning screw: simulation showed it tunes less than 5 kHz. Order the loop in three lip heights (`out/loop_Cu-T2_s10p42.step`, `_s10p55`, `_s10p7`) and pick the one closest to 145 MHz.

## Ordering at JLCCNC

1. Upload `can_body_Cu-T2.step`, `lid_Cu-T2.step` and `loop_Cu-T2.step`, and choose **Copper (T2)**.
2. Add a 2D note covering:
   - threads: M3 ×8 in the flange;
   - the slit width tolerance (+0.02/+0.03 mm over the slab thickness);
   - Ra ≤ 0.8 µm on the RF surfaces;
   - no anodize, no coating.
3. Order `spacer_PTFE.step` as a PTFE part if PTFE is offered; otherwise source it separately.
4. Order the alumina slab from a ceramic-substrate supplier, laser-cut to (t + s) × z, in the thickness chosen.

## Measuring Qu (VNA)

- Use weak, symmetric loop coupling with |S21| at resonance ≈ −30 dB. The simulated prediction to compare against is [`../../sim/lgr_a_driven/results/sparams_d_fil_s1055_xy4p7x4p7_cu2_o3.csv`](../../sim/lgr_a_driven/results/sparams_d_fil_s1055_xy4p7x4p7_cu2_o3.csv): Q_L ≈ 2,080 (3 dB bandwidth ≈ 70 kHz).
- Measure the 3 dB bandwidth to get Q_L, then Qu = Q_L / (1 − |S21|) for symmetric coupling.
- Calibrate at the SMA reference planes.
- Record f0 and Qu against temperature in the heated chamber. This drift is the baseline for thermal self-compensation.

## License

CAD and code: MIT. Measurement data and figures: CC BY 4.0. See the repository root.
