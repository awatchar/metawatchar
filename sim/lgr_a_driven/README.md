# LGR-A in its copper can: coupling loops, S-parameters, tuning screw, lip fillets (Palace)

Driven 2-port and eigenmode Palace models of the LGR-A meta-atom prototype
([`../../hardware/lgr_a`](../../hardware/lgr_a)) with its two SMA coupling loops, the brass M4
tuning screw and the R0.4 mm lip-tip fillets. The cell is the one chosen in the size sweep
([`../lgr_size_sweep`](../lgr_size_sweep)): V = 0.5 L, a = r/R = 0.5, 1.0 mm alumina gap, lip
height s = 10.42 mm (eigenmode there: f0 = 145.33 MHz, Q = 3,927 with sharp edges, no spacers).

Goals:

1. Size the coupling loops for weak, symmetric coupling (|S21(f0)| ≈ −30 dB).
2. Predict the S11/S21/S22 curves that the VNA measurement will be compared with.
3. Tuning range of the M4 screw (penetration 0–6 mm).
4. Effect of the R0.4 mm lip-tip fillets on f0 and Q.

Results: [`results.md`](results.md) and [`results/summary.json`](results/summary.json).

## Model

| Item | Model |
|---|---|
| Can, loop, lips | Copper impedance boundary, σ = 5.8e7 S/m (Rs = 3.145 mΩ/sq at 145.33 MHz, Ls = Rs/ω) |
| Gap | Alumina, εr 9.8, tan δ 4e-5, fills the 1.0 mm slit and the space between the lips |
| PTFE spacers | **Included**: rings r…r+t, height c, between each loop end and the can end; εr 2.1, tan δ 2e-4 |
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
| `configs/*.json` | Palace configurations |
| `results/sparams_<case>.csv` | S-parameters: `f_Hz, S11_re, S11_im, S21_re, S21_im, S22_re, S22_im, S12_re, S12_im` |
| `results/<round>.json` | Per-case summaries; `results/summary.json` collects the final numbers |

Case names: `e_` eigenmode, `d_` driven; `xy<w>x<h>` loop of depth *w* × width *h* (mm);
`_o2`/`_o3` FE order.

## Running

```bash
python build_cases.py sizing --mesh      # then: tune, driven (after `analyze.py sizing`), axial, asbuilt
# copy the folder to the cluster, then from the case directory:
sbatch -A <ACCOUNT> --export=ALL,PALACE_ENV=<PALACE_INSTALL>/env.sh run_cases.sbatch lists/<round>.txt
```

Meshes (`mesh/*.msh`, 30–40 MB each) are not committed; `build_cases.py --mesh` regenerates them.

## License

Code: MIT. Results: CC BY 4.0. See the repository root.
