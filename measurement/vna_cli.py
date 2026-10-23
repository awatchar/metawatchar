"""Measurement chain for the LGR-A meta-atom with the SV4401A VNA: calibration, noise floor,
unloaded Q, and f0/Qu drift against temperature.

Run from this folder (Python with numpy, scipy, matplotlib, pyserial, scikit-rf):

    python vna_cli.py info                                    # identify, list console commands
    python vna_cli.py calibrate --start 144.2e6 --stop 146.2e6 --slot 1
    python vna_cli.py noise --start 144.2e6 --stop 146.2e6 --n 10
    python vna_cli.py q --find 140e6 150e6 --label s1055_clip695 [--reverse]
    python vna_cli.py drift --center 145.1e6 --temp file:chamber_log.csv --interval 60
    python vna_cli.py drift-analyze runs/<run>/drift.csv
    python vna_cli.py fit some_sweep.s2p                      # file saved by the VNA or NanoVNA-Saver

Add `--sim` to any command to run against the built-in simulator instead of hardware, and
`--port COM5` to choose the serial port (default: found by USB ID). Each command writes a run
folder `runs/<UTC time>_<command>[_label]/` with the raw sweeps (CSV), Touchstone files, a JSON
summary and figures. See README.md for the full procedure.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

import q_extraction as qx
from sv4401a import SV4401A, ResonatorDUT, SimulatedSerial, StandardDUT, Sweep

HERE = Path(__file__).resolve().parent


def db(x):
    return 20 * np.log10(np.maximum(np.abs(x), 1e-15))


def run_dir(cmd, label=None, root=None):
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    d = Path(root or HERE / "runs") / (f"{stamp}_{cmd}" + (f"_{label}" if label else ""))
    d.mkdir(parents=True, exist_ok=True)
    return d


def sim_dut(vna, a, kind):
    """With --sim, swap the simulated DUT to match what the operator was asked to connect."""
    if a.sim:
        vna.ser.dut = StandardDUT(kind)


def wait_user(msg, a):
    if a.sim or a.yes:
        print(f"[auto] {msg}")
        return
    input(f"\n>>> {msg}\n    Press Enter when ready (Ctrl-C to abort)... ")


def connect(a, dut=None):
    if a.sim:
        return SV4401A(SimulatedSerial(dut or ResonatorDUT(), noise_db=a.sim_noise_db, f_err_Hz=a.sim_f_err))
    return SV4401A.open(a.port)


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


# ---------------------------------------------------------------------------------------------
def cmd_info(a):
    vna = connect(a)
    ident = vna.identify()
    for k, v in ident.items():
        print(f"{k:9s}: {v}")
    extra = [c for c in ("bandwidth", "power", "avg", "average", "ifbw") if c in vna.commands()]
    print("extra console controls:", extra or "none (set IFBW, power and averaging on the touch screen)")


def cmd_calibrate(a):
    """Guided OSLT calibration at the cable ends (SMA), then a thru verification sweep."""
    vna = connect(a)
    out = run_dir("cal", a.label)
    log = dict(identify=vna.identify(), start_Hz=a.start, stop_Hz=a.stop, points=a.points, slot=a.slot,
               note="Set IFBW, power and averaging on the touch screen first; record them with --setup.",
               setup=a.setup, steps=[])
    vna.command(f"scan {int(a.start)} {int(a.stop)} {a.points}")   # the calibration uses the current sweep
    vna.cal("reset")
    for std, where in (("open", "OPEN on the PORT 1 cable end"), ("short", "SHORT on the PORT 1 cable end"),
                       ("load", "LOAD on the PORT 1 cable end"),
                       ("thru", "PORT 1 and PORT 2 cable ends joined with the SMA-KK (female-female) adapter")):
        wait_user(f"Connect {where}.", a)
        t0 = time.time()
        vna.cal(std)
        log["steps"].append(dict(step=std, seconds=round(time.time() - t0, 1)))
        print(f"  cal {std}: done")
    vna.cal("done")
    vna.cal("on")
    sim_dut(vna, a, "thru")
    if a.slot is not None:
        vna.save_state(a.slot)
        print(f"  calibration saved in slot {a.slot}")
    # verification with the thru still connected
    m, sd, _ = vna.sweep_avg(a.start, a.stop, a.points, n=3)
    m.meta.update(setup=a.setup)
    m.save_csv(out / "verify_thru.csv")
    v = dict(thru_S21_dB_mean=float(np.mean(db(m.s21))), thru_S21_dB_maxdev=float(np.max(np.abs(db(m.s21)))),
             thru_S21_phase_deg_maxdev=float(np.max(np.abs(np.degrees(np.angle(m.s21))))),
             thru_S11_dB_max=float(np.max(db(m.s11))))
    v["pass"] = bool(v["thru_S21_dB_maxdev"] < 0.1 and v["thru_S11_dB_max"] < -30)
    log["verify"] = v
    write_json(out / "calibration.json", log)
    print(json.dumps(v, indent=1))
    print("PASS" if v["pass"] else "CHECK: thru not flat to 0.1 dB or S11 above -30 dB; redo the calibration",
          f"\nwritten: {out}")


def cmd_noise(a):
    """Isolation floor (both cable ends on LOADs) and trace noise (thru, n repeated sweeps)."""
    vna = connect(a)
    out = run_dir("noise", a.label)
    res = dict(setup=a.setup, start_Hz=a.start, stop_Hz=a.stop, points=a.points)
    wait_user("Terminate BOTH cable ends with 50-ohm loads (isolation floor).", a)
    sim_dut(vna, a, "load")
    iso, _, _ = vna.sweep_avg(a.start, a.stop, a.points, n=a.n)
    iso.save_csv(out / "isolation.csv")
    res.update(isolation_S21_dB_mean=float(np.mean(db(iso.s21))), isolation_S21_dB_max=float(np.max(db(iso.s21))))
    wait_user("Join the cable ends with the thru adapter (trace noise).", a)
    sim_dut(vna, a, "thru")
    m, sd, raw = vna.sweep_avg(a.start, a.stop, a.points, n=a.n, keep=True)
    m.save_csv(out / "thru_mean.csv")
    np.savetxt(out / "thru_std_S21.csv", np.column_stack([m.f, sd]), delimiter=",",
               header="f_Hz,std_S21_linear (single sweep)", comments="")
    res.update(n=a.n, trace_noise_dB=float(db(np.median(sd))),
               dynamic_range_dB=float(-np.mean(db(iso.s21))))
    write_json(out / "noise.json", res)
    print(json.dumps(res, indent=1), f"\nwritten: {out}")


def fit_sweep(m: Sweep, s22_f0=None, span_lw=6.0, boot=0):
    fit = qx.fit_s21(m.f, m.s21, span_lw=span_lw)
    i0 = int(np.argmin(np.abs(m.f - fit["f0"])))
    q = qx.unloaded_q(fit, m.s11[i0], s22_f0)
    r = dict(f0_Hz=fit["f0"], sig_f0_Hz=fit["sig_f0"], QL=fit["QL"], sig_QL=fit["sig_QL"],
             S21_f0_dB=20 * math.log10(fit["S21_peak"]), S11_f0_dB=float(db(m.s11[i0])),
             fit_rms=fit["rms"], n_fit_points=int(fit["f"].size), **q)
    r["Qu"] = r.get("Qu_general", r["Qu_symmetric"])
    if boot:
        r["boot_sig_Qu_symmetric"], r["boot_sig_f0_Hz"] = qx.bootstrap(fit, boot)
    return r, fit


def plot_fit(m, fit, r, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fm = np.linspace(m.f.min(), m.f.max(), 4000)
    model = qx.model(fit["params"], fm)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot((m.f - r["f0_Hz"]) / 1e3, db(m.s21), ".", ms=2, label="measured")
    ax[0].plot((fm - r["f0_Hz"]) / 1e3, db(model), lw=1.5, label="fit")
    ax[0].set(xlabel="f − f₀ (kHz)", ylabel="|S21| (dB)",
              title=f"f₀ = {r['f0_Hz']/1e6:.6f} MHz, Q_L = {r['QL']:.0f}, Qᵤ = {r['Qu']:.0f}")
    ax[0].legend()
    ax[0].grid(alpha=0.3)
    ax[1].plot(m.s21.real, m.s21.imag, ".", ms=2)
    ax[1].plot(model.real, model.imag, lw=1.5)
    ax[1].set(xlabel="Re S21", ylabel="Im S21", title="Q circle")
    ax[1].set_aspect("equal", "datalim")
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def find_resonance(vna, lo, hi):
    """Coarse 1001-point search, then a 2 MHz zoom; returns (f0 estimate, 3 dB width estimate)."""
    m = vna.sweep(lo, hi, 1001)
    fpk = m.f[int(np.argmax(np.abs(m.s21)))]
    z = vna.sweep(fpk - 1e6, fpk + 1e6, 1001)
    mag = np.abs(z.s21)
    i = int(np.argmax(mag))
    above = np.flatnonzero(mag >= mag[i] / math.sqrt(2))
    bw = max(z.f[above[-1]] - z.f[above[0]], 3 * (z.f[1] - z.f[0]))
    return z.f[i], bw


def cmd_q(a):
    vna = connect(a)
    out = run_dir("q", a.label)
    if a.find:
        f0, bw = find_resonance(vna, *a.find)
        print(f"resonance near {f0/1e6:.4f} MHz, 3 dB width ≈ {bw/1e3:.1f} kHz")
    else:
        f0, bw = a.center, a.center / 2000
    span = a.span or 12 * bw
    m, sd, _ = vna.sweep_avg(f0 - span / 2, f0 + span / 2, a.points, n=a.avg)
    m.meta.update(label=a.label, setup=a.setup, dut=a.dut)
    m.save_csv(out / "forward.csv")
    m.to_touchstone(out / "forward")
    s22 = None
    if a.reverse:
        wait_user("Swap the two cables at the DUT (port 2 cable to DUT port 1 and vice versa).", a)
        rv, _, _ = vna.sweep_avg(f0 - span / 2, f0 + span / 2, a.points, n=a.avg)
        rv.save_csv(out / "reverse.csv")
        r_rev, _ = fit_sweep(rv)
        s22 = 10 ** (r_rev["S11_f0_dB"] / 20)
    r, fit = fit_sweep(m, s22, boot=a.boot)
    r.update(label=a.label, span_Hz=span, points=a.points, n_avg=a.avg, setup=a.setup, dut=a.dut)
    write_json(out / "q.json", r)
    plot_fit(m, fit, r, out / "q_fit.png")
    print(f"f0 = {r['f0_Hz']/1e6:.6f} MHz ± {r['sig_f0_Hz']:.0f} Hz   Q_L = {r['QL']:.1f} ± {r['sig_QL']:.1f}")
    print(f"|S21(f0)| = {r['S21_f0_dB']:.2f} dB   Qu = {r['Qu']:.1f}"
          + (f" (general; symmetric {r['Qu_symmetric']:.1f})" if "Qu_general" in r else " (symmetric)"))
    print(f"written: {out}")
    return r


def cmd_fit(a):
    """Fit a sweep saved by any software: our CSV, or a Touchstone .s2p (e.g. saved on the
    SV4401A's memory card with SAVE S2P, or exported by NanoVNA-Saver / NanoVNA-App)."""
    p = Path(a.file)
    if p.suffix.lower() == ".s2p":
        import skrf as rf
        nw = rf.Network(str(p))
        m = Sweep(nw.f, nw.s[:, 0, 0], nw.s[:, 1, 0], dict(source=p.name))
        s22 = None
        if a.use_s22:
            fit0 = qx.fit_s21(m.f, m.s21)
            s22 = abs(nw.s[int(np.argmin(np.abs(nw.f - fit0["f0"]))), 1, 1])
    else:
        m, s22 = Sweep.load_csv(p), None
    r, fit = fit_sweep(m, s22, boot=a.boot)
    write_json(p.with_suffix(".q.json"), r)
    plot_fit(m, fit, r, p.with_suffix(".q_fit.png"))
    print(json.dumps({k: v for k, v in r.items() if not isinstance(v, (list, dict))}, indent=1))


class TempSource:
    """Temperature reader: 'manual' (typed), 'file:<csv>' (last number on the last line of a log
    the chamber writes), 'cmd:<shell command>' (prints a number), 'sim:<T0>,<T1>,<steps>' (ramp
    for the simulator; also sets the simulated DUT temperature)."""

    def __init__(self, spec: str, dut=None):
        self.kind, _, self.arg = spec.partition(":")
        self.dut = dut
        if self.kind == "sim":
            t0, t1, n = (float(x) for x in (self.arg or "25,65,9").split(","))
            up = np.linspace(t0, t1, int(n))
            self.ramp = list(np.concatenate([up, up[-2::-1]]))   # up, then back down (hysteresis)

    def read(self):
        if self.kind == "manual":
            return float(input("    chamber temperature (°C): "))
        if self.kind == "file":
            last = Path(self.arg).read_text(encoding="utf-8").strip().splitlines()[-1]
            nums = [x for x in last.replace(";", ",").replace("\t", ",").split(",") if _isnum(x)]
            return float(nums[-1])
        if self.kind == "cmd":
            return float(subprocess.run(self.arg, shell=True, capture_output=True, text=True).stdout.strip().split()[-1])
        if self.kind == "sim":
            if not self.ramp:
                raise StopIteration
            t = float(self.ramp.pop(0))
            if self.dut is not None:
                self.dut.temperature = t
            return t
        raise ValueError(self.kind)


def _isnum(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def cmd_drift(a):
    """Repeated Q measurements while the chamber temperature changes; re-centres on f0 each time.
    A row is flagged settled when f0 moved less than --settle-hz over the last --settle-n rows."""
    dut = ResonatorDUT() if a.sim else None
    vna = connect(a, dut)
    temp = TempSource(a.temp if not a.sim else (a.temp if a.temp.startswith("sim") else "sim:25,65,9"), dut)
    out = run_dir("drift", a.label)
    (out / "sweeps").mkdir()
    rows, center, span = [], a.center, a.span
    t_end = time.time() + a.duration * 60 if a.duration else math.inf
    hdr = "utc,elapsed_s,T_C,f0_Hz,sig_f0_Hz,QL,Qu,S21_f0_dB,fit_rms,settled"
    (out / "drift.csv").write_text(hdr + "\n")
    t_start = time.time()
    print(hdr)
    try:
        while time.time() < t_end:
            T_before = temp.read()
            m, _, _ = vna.sweep_avg(center - span / 2, center + span / 2, a.points, n=a.avg)
            T_after = T_before if temp.kind in ("sim", "manual") else temp.read()
            r, _ = fit_sweep(m)
            k = len(rows)
            m.meta.update(T_C=0.5 * (T_before + T_after), label=a.label)
            m.save_csv(out / "sweeps" / f"sweep_{k:04d}.csv")
            f0s = [x["f0_Hz"] for x in rows[-(a.settle_n - 1):]] + [r["f0_Hz"]] if a.settle_n > 1 else [r["f0_Hz"]]
            settled = len(f0s) >= a.settle_n and (max(f0s) - min(f0s)) < a.settle_hz
            row = dict(utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), elapsed_s=round(time.time() - t_start, 1),
                       T_C=0.5 * (T_before + T_after), f0_Hz=r["f0_Hz"], sig_f0_Hz=r["sig_f0_Hz"], QL=r["QL"],
                       Qu=r["Qu"], S21_f0_dB=r["S21_f0_dB"], fit_rms=r["fit_rms"], settled=int(settled))
            rows.append(row)
            line = ",".join(f"{row[c]:.6g}" if isinstance(row[c], float) else str(row[c]) for c in hdr.split(","))
            with open(out / "drift.csv", "a") as fh:
                fh.write(line + "\n")
            print(line)
            if abs(r["f0_Hz"] - center) < span / 2:
                center = r["f0_Hz"]          # keep the resonance centred as it drifts
            if not a.sim:
                time.sleep(a.interval)
    except (KeyboardInterrupt, StopIteration):
        pass
    print(f"written: {out / 'drift.csv'}  ({len(rows)} points)")
    if rows:
        a.csv, a.all_rows = str(out / "drift.csv"), a.sim   # the simulator ramp has no dwell; use all rows
        cmd_drift_analyze(a)


def cmd_drift_analyze(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = np.genfromtxt(a.csv, delimiter=",", names=True, dtype=None, encoding="utf-8")
    sel = np.ones(d.size, bool) if getattr(a, "all_rows", False) else d["settled"] == 1
    if sel.sum() < 3:
        print(f"only {sel.sum()} settled points; using all {d.size} points")
        sel[:] = True
    T, f0, qu = d["T_C"][sel], d["f0_Hz"][sel], d["Qu"][sel]
    A = np.column_stack([np.ones_like(T), T - 25.0])
    (c0, c1), res, *_ = np.linalg.lstsq(A, f0, rcond=None)
    dof = max(1, T.size - 2)
    s2 = float(np.sum((f0 - A @ [c0, c1]) ** 2) / dof)
    sig_c1 = math.sqrt(s2 * np.linalg.inv(A.T @ A)[1, 1])
    (q0, q1), *_ = np.linalg.lstsq(A, qu, rcond=None)
    out = dict(n_points=int(T.size), T_range_C=[float(T.min()), float(T.max())],
               f0_25C_Hz=float(c0), df0_dT_Hz_per_K=float(c1), sig_df0_dT_Hz_per_K=sig_c1,
               tau_f_ppm_per_K=float(c1 / c0 * 1e6), sig_tau_f_ppm_per_K=float(sig_c1 / c0 * 1e6),
               Qu_25C=float(q0), dQu_dT_per_K=float(q1), dQu_dT_pct_per_K=float(q1 / q0 * 100),
               f0_resid_rms_Hz=math.sqrt(s2))
    rising = np.gradient(d["T_C"])[sel] >= 0
    if rising.any() and (~rising).any():
        for tag, msk in (("up", rising), ("down", ~rising)):
            if msk.sum() >= 2:
                out[f"df0_dT_{tag}_Hz_per_K"] = float(np.polyfit(T[msk] - 25, f0[msk], 1)[0])
    p = Path(a.csv).with_name("drift_summary.json")
    write_json(p, out)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(T, (f0 - c0) / 1e3, "o", ms=4)
    tt = np.linspace(T.min(), T.max(), 50)
    ax[0].plot(tt, c1 * (tt - 25) / 1e3, "-")
    ax[0].set(xlabel="T (°C)", ylabel="f₀ − f₀(25 °C) (kHz)", title=f"τf = {out['tau_f_ppm_per_K']:.1f} ± "
              f"{out['sig_tau_f_ppm_per_K']:.1f} ppm/K")
    ax[1].plot(T, qu, "o", ms=4)
    ax[1].set(xlabel="T (°C)", ylabel="Qᵤ", title=f"dQᵤ/dT = {out['dQu_dT_pct_per_K']:.3f} %/K")
    for x in ax:
        x.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Path(a.csv).with_name("drift.png"), dpi=150)
    plt.close(fig)
    print(json.dumps(out, indent=1))


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port")
    ap.add_argument("--sim", action="store_true", help="use the built-in simulator")
    ap.add_argument("--sim-noise-db", type=float, default=-75.0)
    ap.add_argument("--sim-f-err", type=float, default=0.0, help="simulated VNA frequency error (Hz)")
    ap.add_argument("--yes", action="store_true", help="do not wait for Enter between steps")
    ap.add_argument("--setup", default="", help="free text recorded with the data: IFBW, power, averaging, cables")
    ap.add_argument("--label")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("info")
    c = sub.add_parser("calibrate")
    c.add_argument("--start", type=float, required=True)
    c.add_argument("--stop", type=float, required=True)
    c.add_argument("--points", type=int, default=1001)
    c.add_argument("--slot", type=int)
    n = sub.add_parser("noise")
    n.add_argument("--start", type=float, required=True)
    n.add_argument("--stop", type=float, required=True)
    n.add_argument("--points", type=int, default=1001)
    n.add_argument("--n", type=int, default=10)
    q = sub.add_parser("q")
    g = q.add_mutually_exclusive_group(required=True)
    g.add_argument("--find", type=float, nargs=2, metavar=("LO", "HI"))
    g.add_argument("--center", type=float)
    q.add_argument("--span", type=float, help="sweep span (Hz); default 12 x the 3 dB width")
    q.add_argument("--points", type=int, default=1001)
    q.add_argument("--avg", type=int, default=4, help="software averages of complex sweeps")
    q.add_argument("--reverse", action="store_true", help="also measure with the DUT reversed (gives S22)")
    q.add_argument("--boot", type=int, default=200, help="bootstrap resamples for the Qu uncertainty")
    q.add_argument("--dut", default="", help="free text: loop variant, clip opening, assembly number")
    dr = sub.add_parser("drift")
    dr.add_argument("--center", type=float, required=True)
    dr.add_argument("--span", type=float, default=800e3)
    dr.add_argument("--points", type=int, default=1001)
    dr.add_argument("--avg", type=int, default=2)
    dr.add_argument("--temp", default="manual", help="manual | file:<csv> | cmd:<command> | sim:T0,T1,steps")
    dr.add_argument("--interval", type=float, default=60.0, help="seconds between points")
    dr.add_argument("--duration", type=float, help="minutes (default: until Ctrl-C)")
    dr.add_argument("--settle-hz", type=float, default=300.0)
    dr.add_argument("--settle-n", type=int, default=3)
    fi = sub.add_parser("fit", help="fit a saved CSV or .s2p file")
    fi.add_argument("file")
    fi.add_argument("--boot", type=int, default=200)
    fi.add_argument("--use-s22", action="store_true", help="use the file's S22 (only if really measured)")
    da = sub.add_parser("drift-analyze")
    da.add_argument("csv")
    da.add_argument("--all-rows", action="store_true")
    a = ap.parse_args()
    {"info": cmd_info, "calibrate": cmd_calibrate, "noise": cmd_noise, "q": cmd_q, "drift": cmd_drift,
     "drift-analyze": cmd_drift_analyze, "fit": cmd_fit}[a.cmd](a)


if __name__ == "__main__":
    main()
