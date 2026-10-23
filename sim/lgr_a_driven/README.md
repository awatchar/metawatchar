# LGR-A in its copper can: coupling loops, S-parameters, tuning screw, lip fillets (Palace)

Driven 2-port and eigenmode Palace models of the LGR-A meta-atom prototype
([`../../hardware/lgr_a`](../../hardware/lgr_a)) with its two SMA coupling loops, the brass M4
tuning screw and the R0.4 mm lip-tip fillets. The cell is the one chosen in the size sweep
([`../lgr_size_sweep`](../lgr_size_sweep)): V = 0.5 L, a = r/R = 0.5, 1.0 mm alumina gap, lip
height s = 10.42 mm (eigenmode there: f0 = 145.33 MHz, Q = 3,927 with sharp edges, no spacers; first-run copper input, corrected value 2,149, see below).

Goals:

1. Size the coupling loops for weak, symmetric coupling (|S21(f0)| ≈ −30 dB).
2. Predict the S11/S21/S22 curves that the VNA measurement will be compared with.
3. Tuning range of the M4 screw (penetration 0–6 mm).
4. Effect of the R0.4 mm lip-tip fillets on f0 and Q.

Results: [`results.md`](results.md), [`results/summary_cu2.json`](results/summary_cu2.json)
(corrected copper input, current) and [`results/summary.json`](results/summary.json) (first runs).

> **Copper input correction.** Palace's `Impedance` boundary combines Rs and Ls **in parallel**.
> The first runs used Rs with Ls = Rs/ω, i.e. half the copper loss, so their Q values are too high
> (Q_c 2×). The `_cu2` runs use Rs_p = 2 Rs, Ls_p = 2 Rs/ω (Z = Rs(1 + j) at the case's f0) and
> supersede them: as built, Q_u ≈ 2,150 (not 3,930), and the loops for −30 dB are 4.7 × 4.7 mm
> (not 4 × 4 mm). See [`results.md`](results.md#corrected-copper-input-_cu2-runs).

## Model

| Item | Model |
|---|---|
| Can, loop, lips | Copper impedance boundary, σ = 5.8e7 S/m. Palace input (parallel Rs ‖ Ls): `_cu2` runs Rs_p = 2 Rs, Ls_p = 2 Rs/ω with Rs = √(πfµ0/σ) at the case's f0 (≈ 6.29 mΩ/sq entered for Rs = 3.144 mΩ/sq at 145.2 MHz); first runs Rs = 3.145 mΩ/sq, Ls = Rs/ω (half the loss) |
| Gap | Alumina, εr 9.8, tan δ 4e-5, fills the 1.0 mm slit and the space between the lips |
| PTFE spacers | **Included**: rings r…r+t, height c, between each loop end and the can end; εr 2.1, tan δ 2e-4 |
| PEEK slab clips (`_clip` case) | Two U-clips as in `hardware/lgr_a` (`slab_clip`): outer 6 mm of both lips, 4 mm along z at each loop end, 2 mm walls; εr 3.2, tan δ 3e-3 (assumed); interference not modelled |
| Lip height | s = 10.42 mm (first runs); `_cu2` runs also s = 10.55 mm (default build) and 10.70 mm |
| SMA ports (120°, 240°, mid-height) | 50 Ω coaxial **lumped port** (annulus, `"Direction": "+R"`) at the end of a 1 mm deep, 50 Ω air-coax recess in the can wall (inner Ø 1.0 mm = the loop wire, outer Ø 2.3 mm). The reference plane is therefore 1 mm behind the inner wall; the SMA connector adds only a phase shift. |
| Coupling loops | 1.0 mm diameter copper wire (impedance boundary). From the port the wire runs radially into the can by the depth *w*, turns and runs tangentially by the width *h* (towards 180°, i.e. away from the gap), and returns radially to the wall, where it is grounded. The loop lies in the mid-height plane, so its normal is the can axis and it links the axial return flux between the LGR and the can. Dimensions are wire-centreline values. The two loops are mirror images about the x–z plane. |
| Loop orientation check | One 5 × 5 mm loop in the radial–axial plane (normal azimuthal), `e_rz5x5_o2` |
| Tuning screw (0°, mid-height, facing the lip tips) | Brass impedance boundary, σ = 1.5e7 S/m (typical free-cutting brass; to be verified for the actual screw), Ø 4.0 mm cylinder with a flat end, penetration p = 0 (plain wall), 2, 4, 6 mm past the inner wall. Threads are not modelled. |
| Lip fillets | R0.4 mm on the four axial edges at the lip tips, as in `hardware/lgr_a/lgr_a_assembly.py` |
| Exploratory tuner | Brass Ø6 mm rod through the lid on the can axis, 20 or 30 mm deep (reaches into the loop bore) |

Mesh: Gmsh 4.15, second-order curved tetrahedra, the size-sweep "v2" recipe (0.2 mm on every
loop/lip edge, which controls the conductor Q) plus 0.3 mm on the coupling wires and 0.15 mm at
the port annuli. Palace 0.18.1, FE order 2 (sizing, driven sweeps) and 3 (tuning, fillets,
convergence).

- **Eigenmode with lumped ports.** A 50 Ω port is a resistive boundary, so the eigenvalue
  includes the port loading. Palace's `port-Q.csv` gives the external Q of each port.
  Q_u = (1/Q_L − 1/Q_ext,1 − 1/Q_ext,2)⁻¹, β_i = Q_u/Q_ext,i, and
  |S21(f0)| = 2√(β1β2)/(1 + β1 + β2).
- **Driven.** Two excitations (port 1, then port 2) give the full 2 × 2 S-matrix. The adaptive
  fast frequency sweep (`AdaptiveTol` 1e-4) builds a reduced-order model; the output step is
  1 kHz against a 3 dB linewidth of about 38 kHz. `postprocess.py` fits the resonance
  (S21 peak, 3 dB Q_L, S11/S22 circle diameters d_i) and gives
  β_i = d_i/(2 − d1 − d2) and Q_u = Q_L(1 + β1 + β2), the same reduction that is applied to the
  measured data.

## Files

| File | What it is |
|---|---|
| `lgr_a_model.py` | Geometry (Gmsh/OpenCASCADE), mesh and Palace configuration |
| `build_cases.py` | Defines every case; meshes locally and writes `configs/`, `cases/`, `lists/` |
| `run_cases.sbatch` | Slurm job: runs the cases of a list, several at a time on one node, then `postprocess.py` |
| `postprocess.py` | Eigenmode and driven summaries (standard library only), S-parameter CSVs |
| `analyze.py` | Loop-size fit (`results/sizing.json`) and the study summary (`results/summary.json`) |
| `fields_fig2.py` | Field maps of the fundamental mode from a run with saved fields (paper P1 Fig. 2): `figures/fig2*.png/pdf`, `results/fields_<case>.json`; runs on the cluster after the solve (numpy, matplotlib, pyvista/vtk, cmocean) |
| `figures/` | Field-map figures (PNG 600 dpi and PDF) |
| `configs/*.json` | Palace configurations |
| `results/sparams_<case>.csv` | S-parameters: `f_Hz, S11_re, S11_im, S21_re, S21_im, S22_re, S22_im, S12_re, S12_im` |
| `results/<round>.json` | Per-case summaries; `results/summary.json` collects the first runs' numbers, `results/summary_cu2.json` the corrected ones (`python analyze.py cu2`) |

Case names: `e_` eigenmode, `d_` driven; `xy<w>x<h>` loop of depth *w* × width *h* (mm, `p` =
decimal point); `fil` R0.4 fillets; `s1055` lip height 10.55 mm; `clip` PEEK clips; `_cu2`
corrected copper input; `_o2`/`_o3`/`_o4` FE order.
Groups in `build_cases.py`: `sizing`, `tune`, `driven`, `axial`, `asbuilt` (first runs, copper
input reproduced exactly with `zs="s2"`), `cu2_eig`, `cu2_loops`, `cu2_o4` (corrected), `cu2_fields`
(default build, order 2, fields saved for the field maps).

## Running

```bash
python build_cases.py cu2_eig --mesh     # corrected copper input; also cu2_loops, cu2_o4
# first runs: sizing, tune, driven (after `analyze.py sizing`), axial, asbuilt
# copy the folder to the cluster, then from the case directory:
sbatch -A <ACCOUNT> --export=ALL,PALACE_ENV=<PALACE_INSTALL>/env.sh run_cases.sbatch lists/<round>.txt
```

Meshes (`mesh/*.msh`, 30–40 MB each) are not committed; `build_cases.py --mesh` regenerates them.

## License

Code: MIT. Results: CC BY 4.0. See the repository root.
