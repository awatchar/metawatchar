"""Driver for the SYSJOINT SV4401A handheld VNA over its USB serial console, plus a simulator.

Protocol (SV4401A User Manual V2.0/V4.0, ch. 8): a text console on a USB CDC serial port (baud rate
auto-detected; 115200 is used). Each command line ends with CR. The device echoes the command,
prints its output lines, and ends with the prompt "ch>". Documented commands used here:

    scan {start_Hz} {stop_Hz} [points 101..1001] [outmask]   set and run a sweep
    frequencies                                             frequency list of the sweep (Hz)
    data 0 | data 1                                         S11 | S21 as "re im" per line
    cal open|short|load|thru|done|reset|on|off             calibration steps
    save [0..17] / recall [0..17]                           calibration-state slots
    help, info, version, sn

The console has no documented commands for IF bandwidth, output power or averaging; those are set
on the touch screen (the driver reports whether `help` lists any extra commands, e.g. `bandwidth`).
The instrument measures S11 and S21 only (port 1 drives); S22 needs the DUT reversed.

The simulator (`SimulatedSerial`) answers the same console commands for a model DUT, so the whole
measurement chain can be exercised without hardware: `SV4401A.simulated(dut)`.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np

PROMPT = "ch>"
USB_IDS = {(0x04B4, 0x0008)}  # SV4401A USB VID:PID as listed by the vendor's PC software


@dataclass
class Sweep:
    f: np.ndarray                  # Hz
    s11: np.ndarray                # complex
    s21: np.ndarray                # complex
    meta: dict = field(default_factory=dict)

    def save_csv(self, path, extra_header: str = ""):
        cols = np.column_stack([self.f, self.s11.real, self.s11.imag, self.s21.real, self.s21.imag])
        hdr = "".join(f"# {k}: {v}\n" for k, v in self.meta.items()) + extra_header
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(hdr)
            fh.write("f_Hz,S11_re,S11_im,S21_re,S21_im\n")
            np.savetxt(fh, cols, delimiter=",", fmt="%.10g")

    @staticmethod
    def load_csv(path) -> "Sweep":
        meta = {}
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                k, _, v = line[1:].partition(":")
                meta[k.strip()] = v.strip()
        d = np.genfromtxt(path, delimiter=",", names=True, comments="#")
        return Sweep(d["f_Hz"], d["S11_re"] + 1j * d["S11_im"], d["S21_re"] + 1j * d["S21_im"], meta)

    def to_touchstone(self, path_stem):
        """Write S11 (.s1p) and a 2-port file with S12 = S21 and S22 = S11 marked as assumed."""
        import skrf as rf
        s = np.zeros((self.f.size, 2, 2), complex)
        s[:, 0, 0], s[:, 1, 0], s[:, 0, 1], s[:, 1, 1] = self.s11, self.s21, self.s21, self.s11
        nw = rf.Network(frequency=rf.Frequency.from_f(self.f, unit="hz"), s=s, z0=50,
                        comments="SV4401A: S11, S21 measured; S12 := S21, S22 := S11 (assumed, one-path VNA)")
        nw.write_touchstone(str(path_stem))


class SV4401A:
    def __init__(self, ser, settle_s: float = 0.0):
        self.ser = ser
        self.settle_s = settle_s
        self._commands = None

    # ---- connection -------------------------------------------------------------------------
    @classmethod
    def open(cls, port: str | None = None, timeout: float = 2.0, **kw) -> "SV4401A":
        import serial
        import serial.tools.list_ports as lp
        if port is None:
            cands = [p.device for p in lp.comports() if (p.vid, p.pid) in USB_IDS]
            if not cands:
                names = ", ".join(f"{p.device} ({p.description})" for p in lp.comports()) or "none"
                raise IOError(f"SV4401A not found by USB ID; pass --port. Serial ports: {names}")
            port = cands[0]
        ser = serial.Serial(port, 115200, timeout=timeout)
        vna = cls(ser, **kw)
        vna._flush()
        return vna

    @classmethod
    def simulated(cls, dut=None, **kw) -> "SV4401A":
        return cls(SimulatedSerial(dut), **kw)

    def close(self):
        self.ser.close()

    def _flush(self):
        self.ser.write(b"\r\n")
        time.sleep(0.1)
        self.ser.reset_input_buffer()

    # ---- console ----------------------------------------------------------------------------
    def command(self, cmd: str, timeout: float = 30.0) -> list[str]:
        """Send one command line; return the output lines (echo and prompt removed)."""
        self.ser.reset_input_buffer()
        self.ser.write((cmd + "\r").encode("ascii"))
        out, t0 = [], time.monotonic()
        while True:
            raw = self.ser.readline()
            if not raw:
                if time.monotonic() - t0 > timeout:
                    raise TimeoutError(f"no prompt after '{cmd}' within {timeout} s; got {out[-3:]}")
                continue
            line = raw.decode("ascii", "replace").strip()
            if line.startswith(PROMPT):
                return out
            if line and line != cmd:
                out.append(line)

    def commands(self) -> list[str]:
        if self._commands is None:
            self._commands = " ".join(self.command("help")).replace("Commands:", " ").split()
        return self._commands

    def identify(self) -> dict:
        return dict(version=" ".join(self.command("version")), info=" | ".join(self.command("info")),
                    sn=" ".join(self.command("sn")), commands=" ".join(self.commands()))

    # ---- sweeps -----------------------------------------------------------------------------
    def sweep(self, start: float, stop: float, points: int = 1001, timeout: float = 120.0) -> Sweep:
        if not 101 <= points <= 1001:
            raise ValueError("points must be 101..1001")
        self.command(f"scan {int(round(start))} {int(round(stop))} {points}", timeout=timeout)
        if self.settle_s:
            time.sleep(self.settle_s)
        f = np.array([float(x.split()[0]) for x in self.command("frequencies")])
        s11 = _parse_complex(self.command("data 0"))
        s21 = _parse_complex(self.command("data 1"))
        if not (f.size == s11.size == s21.size == points):
            raise IOError(f"sweep size mismatch: f {f.size}, s11 {s11.size}, s21 {s21.size}, expected {points}")
        return Sweep(f, s11, s21, dict(start_Hz=start, stop_Hz=stop, points=points,
                                       time=time.strftime("%Y-%m-%dT%H:%M:%S%z")))

    def sweep_avg(self, start, stop, points=1001, n=1, keep=False):
        """Mean of n complex sweeps (on top of any on-screen averaging). Returns (mean, std_s21, raw)."""
        raw = [self.sweep(start, stop, points) for _ in range(n)]
        if any(not np.array_equal(r.f, raw[0].f) for r in raw):
            raise IOError("frequency grid changed between sweeps")
        s11 = np.mean([r.s11 for r in raw], axis=0)
        s21 = np.mean([r.s21 for r in raw], axis=0)
        sd = np.std([r.s21 for r in raw], axis=0) if n > 1 else np.zeros(points)
        m = Sweep(raw[0].f, s11, s21, dict(raw[0].meta, n_avg=n))
        return m, sd, (raw if keep else None)

    # ---- calibration ------------------------------------------------------------------------
    def cal(self, step: str, timeout: float = 120.0):
        if step not in ("open", "short", "load", "thru", "done", "reset", "on", "off"):
            raise ValueError(step)
        return self.command(f"cal {step}", timeout=timeout)

    def save_state(self, slot: int):
        return self.command(f"save {int(slot)}")

    def recall_state(self, slot: int):
        return self.command(f"recall {int(slot)}")


def _parse_complex(lines) -> np.ndarray:
    v = [tuple(map(float, ln.split()[:2])) for ln in lines if ln and ln[0] in "+-.0123456789"]
    return np.array([a + 1j * b for a, b in v])


# ---------------------------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------------------------
class ResonatorDUT:
    """Two-port weakly coupled resonator with leakage; f0 and Qu may depend on temperature.

    f0(T) = f0 * (1 + tf * (T - 25)), Qu(T) = qu * (1 + tq * (T - 25)); `temperature` is set by
    the caller (e.g. the drift test)."""

    def __init__(self, f0=145.2e6, qu=2150.0, beta=0.0165, tf=-75e-6, tq=-1.8e-3, leak=3e-5):
        self.f0, self.qu, self.beta, self.tf, self.tq, self.leak = f0, qu, beta, tf, tq, leak
        self.temperature = 25.0

    def __call__(self, f):
        dT = self.temperature - 25.0
        f0 = self.f0 * (1 + self.tf * dT)
        qu = self.qu * (1 + self.tq * dT)
        ql = qu / (1 + 2 * self.beta)
        x = 2 * ql * (f - f0) / f0
        t0 = 2 * self.beta / (1 + 2 * self.beta)
        s21 = t0 * np.exp(1j * 0.4) / (1 + 1j * x) + self.leak * np.exp(1j * 1.9)
        s11 = -(1 - t0 / (1 + 1j * x)) * np.exp(-1j * 0.3)
        return s11, s21


class StandardDUT:
    """Ideal 'thru' (cable ends joined) or 'load' (both ends terminated) for the simulator."""

    def __init__(self, kind):
        self.kind = kind

    def __call__(self, f):
        one, zero = np.ones(f.size, complex), np.zeros(f.size, complex)
        return (zero, one) if self.kind == "thru" else (zero, zero)


class SimulatedSerial:
    """Minimal stand-in for the SV4401A console: echo, output lines, 'ch>' prompt.
    Noise: complex Gaussian at `noise_db` (dB re 1) on S11 and S21; frequency offset `f_err_Hz`."""

    def __init__(self, dut=None, noise_db=-75.0, f_err_Hz=0.0, seed=0):
        self.dut = dut or ResonatorDUT()
        self.noise = 10 ** (noise_db / 20)
        self.f_err = f_err_Hz
        self.rng = np.random.default_rng(seed)
        self.buf: list[bytes] = []
        self.sweep = (144.7e6, 145.7e6, 1001)
        self.data = None
        self.cal_steps = []

    def write(self, b: bytes):
        cmd = b.decode("ascii").strip()
        if not cmd:
            return
        self.buf.append((cmd + "\r\n").encode())
        for line in self._exec(cmd):
            self.buf.append((line + "\r\n").encode())
        self.buf.append(b"ch> ")

    def readline(self):
        return self.buf.pop(0) if self.buf else b""

    def reset_input_buffer(self):
        self.buf.clear()

    def close(self):
        pass

    def _measure(self):
        a, b, n = self.sweep
        f = np.linspace(a, b, n)
        s11, s21 = self.dut(f + self.f_err)   # instrument frequency error: true f differs from reported
        nz = lambda: self.noise * (self.rng.standard_normal(n) + 1j * self.rng.standard_normal(n)) / math.sqrt(2)
        self.data = (f, s11 + nz(), s21 + nz())

    def _exec(self, cmd):
        w = cmd.split()
        if w[0] == "help":
            return ["Commands: help reset data frequencies scan sweep pause resume cal save recall marker "
                    "version info sn resolution"]
        if w[0] == "version":
            return ["SIMULATOR 0.0"]
        if w[0] == "info":
            return ["Simulated SV4401A console (measurement/sv4401a.py)"]
        if w[0] == "sn":
            return ["SIM-0000"]
        if w[0] == "scan":
            self.sweep = (float(w[1]), float(w[2]), int(w[3]) if len(w) > 3 else self.sweep[2])
            self._measure()
            return []
        if w[0] == "frequencies":
            return [f"{int(round(x))}" for x in self.data[0]]
        if w[0] == "data":
            arr = self.data[2] if (len(w) > 1 and w[1] == "1") else self.data[1]
            return [f"{z.real:.9e} {z.imag:.9e}" for z in arr]
        if w[0] == "cal":
            self.cal_steps.append(w[1] if len(w) > 1 else "")
            return []
        if w[0] in ("save", "recall"):
            return []
        return [f"{w[0]}?"]
