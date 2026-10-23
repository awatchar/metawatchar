# Measurement chain: SV4401A VNA → unloaded Q and thermal drift

These scripts measure the LGR-A meta-atom (`../hardware/lgr_a`) with a SYSJOINT SV4401A vector network analyzer, controlled from a PC over its USB serial console. They cover:

1. calibration;
2. the noise floor and trace noise;
3. the unloaded Q;
4. f₀(T) and Qᵤ(T) in a heated chamber.

The same chain runs against a built-in simulator (`--sim`), so every step can be tested before the hardware arrives.

| File | Contents |
|---|---|
| `sv4401a.py` | Console driver (`SV4401A`), sweep container (`Sweep`: CSV/Touchstone I/O), simulator (`SimulatedSerial`, `ResonatorDUT`, `StandardDUT`) |
| `vna_cli.py` | Command-line tool: `info`, `calibrate`, `noise`, `q`, `fit`, `drift`, `drift-analyze` |
| `q_extraction.py` | Complex Lorentzian + background fit of S21, Qᵤ (symmetric and general coupling), bootstrap uncertainty |

Requirements: Python 3.10+ with numpy, scipy, matplotlib, pyserial and scikit-rf.

## Why a separate driver?

SYSJOINT ships interactive PC programs for the SV4401A: modified NanoVNA-Saver and NanoVNA-App builds, and a screen-casting tool. They are good for looking at traces and saving Touchstone files by hand. They don't log a chamber temperature, re-centre the sweep as f₀ drifts, wait for thermal settling, or fit Qᵤ automatically.

This driver is a thin layer over the console commands documented in the user manual (chapter 8), and everything above it is the measurement method itself. Files saved by the vendor software, or by the instrument's own SAVE S2P, can still be analysed with `vna_cli.py fit file.s2p`.

## Instrument facts that shape the method

Values are from the SV4401A User Manual V2.0/V4.0.

- **Console:** CR-terminated commands, echoed, with the output ending in a `ch>` prompt. The commands used here are `scan`, `frequencies`, `data 0|1`, `cal open|short|load|thru|done|on|off`, `save`/`recall` and `help`/`info`/`version`/`sn`.
- **Touch-screen-only settings:** IF bandwidth, output power and on-screen averaging can't be set from the console. Set them on the screen and record them with `--setup "IFBW 1 kHz, -12 dBm, avg 1"`. `info` reports whether the firmware's `help` lists any extra commands.
- **S11 and S21 only:** port 1 always drives. S22 needs the DUT reversed (`q --reverse`).
- **Sweep and dynamic range:**
  - 101–1001 points, IFBW 12.5 kHz down to 30 Hz, power −42 to −12 dBm;
  - S21 dynamic range 75 dB at 12.5 kHz IFBW and no averaging, better at narrower IFBW.
  - The LGR-A resonance is at |S21| ≈ −30 dB with a 70 kHz linewidth. The fit recovers Qᵤ to better than 0.3 % with noise at −70 dB, so 1 kHz IFBW is ample.
- **Frequency accuracy is ±5 kHz above 140 MHz.**
  - This is an absolute offset of up to 34 ppm at 145 MHz. It shifts the reported f₀ but affects Q only at the ppm level.
  - For the drift run, what matters is that the VNA's frequency is **stable** during the run. Keep the VNA outside the chamber at a steady room temperature, and let it warm up for 30 min.
  - Check its CW output against a spectrum analyzer or frequency counter before and after a run. The expected drift of an alumina-gap cell (about −75 ppm/K, or −11 kHz/K) is far above any plausible VNA reference drift.

## Procedure

Run every command from this folder. Add `--port COMx` if the VNA is not found automatically, and `--sim` to rehearse.

### 0. Identify

```
python vna_cli.py info
```

### 1. Calibrate at the SMA cable ends

Use the SMA OPEN/SHORT/LOAD standards and the SMA-KK thru adapter supplied with the instrument, on the two RG316 cables that will connect the DUT.

First set the IFBW, power and averaging on the screen. The tool then sets the sweep, asks for each standard in turn, saves the state in the chosen slot, and checks the result with the thru still connected. The check passes if |S21| is flat to 0.1 dB and |S11| is below −30 dB.

```
python vna_cli.py --setup "IFBW 1 kHz, -12 dBm, avg 1, 2 x 50 cm RG316" calibrate --start 144.2e6 --stop 146.2e6 --slot 1
```

The calibration covers the band in which every lip-height variant, both clip states and the full drift range will fall: about 144.2–146.2 MHz, with f₀ dropping ≈ 0.5 MHz over 25–65 °C. For a Q measurement the calibration matters little. Q comes from the linewidth, and the coupling correction uses |S21(f₀)|, where a 0.1 dB error changes Qᵤ by only 0.04 %.

### 2. Noise floor and trace noise

```
python vna_cli.py noise --start 144.2e6 --stop 146.2e6 --n 10
```

The tool first asks for both cable ends on LOADs, which gives the isolation floor, and then for the thru, which gives the trace noise. Both numbers go into the uncertainty budget.

### 3. Unloaded Q

```
python vna_cli.py --setup "..." --label s1055_clip695_a1 q --find 140e6 150e6 --dut "loop s=10.55, clip 6.95, assembly 1" --reverse
```

The tool does a coarse search, then sweeps 12 linewidths around f₀ (≈ 840 kHz; 1001 points, 4 averages). It fits S21 and reports f₀, Q_L, |S21(f₀)| and Qᵤ with its bootstrap uncertainty. With `--reverse`, it asks for the DUT to be reversed so it can measure S22 and use the general coupling formula.

Outputs in `runs/<time>_q_<label>/`:
- `forward.csv` / `forward.s2p` (and `reverse.csv` with `--reverse`);
- `q.json`;
- `q_fit.png`, showing |S21| with the fit and the Q circle.

Repeat this for each lip-height variant, each clip opening, and one re-assembly per variant, to get the assembly scatter.

### 4. Thermal drift

```
python vna_cli.py --label s1055_clip695 drift --center 145.1e6 --temp file:chamber_log.csv --interval 60 --duration 600
```

- **Loop:** every `--interval` seconds the tool reads the chamber temperature, sweeps 800 kHz around the last f₀, fits it, and appends a row to `drift.csv`. All sweeps are kept in `sweeps/`.
- **Settled rows:** a row is marked `settled` when f₀ has moved less than `--settle-hz` (300 Hz, about 0.03 K) over the last `--settle-n` (3) rows. Hold each temperature step until it settles.
- **Temperature sources:**
  - `file:<csv>` reads the last number on the last line of a log the chamber controller writes;
  - `cmd:<command>` runs a command that prints the temperature;
  - `manual` asks for it at each point.
- **Analysis:** `drift-analyze` fits f₀(T) and Qᵤ(T) through the settled rows. It gives τf (ppm/K) with its standard error, dQᵤ/dT, and separate up and down slopes for hysteresis, and writes `drift_summary.json` and `drift.png`.
- **Control run:** repeat the run with the loose 7.05 mm clips. Comparing τf with the clamped run measures the air-film contribution.

```
python vna_cli.py drift-analyze runs/<run>/drift.csv
```

## Rehearsal with the simulator

```
python vna_cli.py --sim --label test q --find 140e6 150e6
python vna_cli.py --sim --yes calibrate --start 144.2e6 --stop 146.2e6 --slot 1
python vna_cli.py --sim --yes noise --start 144.2e6 --stop 146.2e6 --n 5
python vna_cli.py --sim --label ramp drift --center 145.2e6 --temp sim:25,65,9
```

The simulated DUT is the predicted LGR-A: f₀ 145.2 MHz, Qᵤ 2,150, β = 0.0165 per port, τf −75 ppm/K and dQᵤ/dT −0.18 %/K, with −75 dB complex noise. The chain recovers:

| Quantity | Recovered |
|---|---|
| Qᵤ | 2,151 |
| Q_L | 2,082.5 ± 0.8 |
| τf | −75.06 ± 0.04 ppm/K |
| dQᵤ/dT | −0.179 %/K |

The same fit applied to the Palace driven prediction (`../sim/lgr_a_driven/results/sparams_d_fil_s1055_xy4p7x4p7_cu2_o3.csv`) gives Q_L 2,081.6 and Qᵤ 2,150.2.

## Status

The driver follows the documented console protocol and has been tested against the simulator only. On first contact with a real SV4401A, check:

- the `help` list;
- the line format returned by `data 0|1` and `frequencies`;
- how long each `cal` step takes before the prompt returns. The timeout is 120 s; raise it if a narrow IFBW makes the steps slower.

Measured data will be added under `data/` (CC BY 4.0). Working run folders (`runs/`) are not tracked.
