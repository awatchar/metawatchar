#!/usr/bin/env python3
"""Effective-medium (mu_eff) predictions for the mixed S4 media.

Each meta-atom type j is described by a Lorentz permeability built from its own two unit-cell
eigenmodes only (the "B-eig" model of ../s3_loaded_resonator/s3_effective_medium.py):

    mu_j(f) = 1 - F_j f^2 / (f^2 - f0_j^2),

with mu_j(f_mp,j) = 0 (PEC-box mode, mu_eff = 0) and mu_j(f_m,j) = (c / (2 f_m,j a))^2 / eps
(PMC-box mode, k = pi / a), eps = 2.385 (median retrieved eps_eff of the S2 slab). Temperature
enters only through the cells' own resonances: f_m,j(T) = f_m,j (1 + tau_j dT) and
f_mp,j(T) = f_mp,j (1 + tau_mp,j dT), with tau_j, tau_mp,j the permittivity part of the cell's
tau_f from its own unit-cell eigenmode (participation formula). The common copper expansion
(-alpha_Cu) and the gap-width term are added outside this model, exactly as for the full-wave
results, so the comparison below is for the permittivity part only.

Two descriptions of a mixed medium of N cells:
  mix     homogeneous medium with the volume-average permeability <mu>(f) = (1/N) sum_j mu_j(f)
          (exact homogenization for H tangential to layers stacked along the line);
  layer   every cell is a layer of thickness a with its own mu_j (cell-resolved effective medium;
          sees the position of a TiO2 cell, the homogeneous mixture does not).
Resonator: the S3 study-1 line (TEM unit-cell line, PEC shorts, empty end sections g = a/2), the
MNG-air mode (mu_eff L + 2 g ~ 0 for a homogeneous slab, S3). For the closed copper box the
medium-level prediction is the mixture rule for the mu_eff = 0 frequency f_mp of <mu>.

Usage: python s4_effective_medium.py   (reads results/*_eig.csv via analyze_s4 helpers)
"""
import math

import numpy as np

C0 = 299_792_458.0
A = 0.070
EPS_EFF = 2.38463026265212          # S3 em_parameters.json, median retrieved eps_eff (N = 4)


def beig_params(f_m, f_mp, eps=EPS_EFF):
    """Lorentz (f0, F) from the PMC-box (k = pi/a) and PEC-box (mu = 0) unit-cell eigenmodes."""
    mu_m = (C0 / (2 * f_m * A)) ** 2 / eps
    # 1 - F f_m^2/(f_m^2 - f0^2) = mu_m with F = 1 - f0^2/f_mp^2 -> solve for f0 by bisection
    lo, hi = f_m * (1 + 1e-9), f_m * 1.05

    def resid(f0):
        F = 1 - f0 ** 2 / f_mp ** 2
        return 1 - F * f_m ** 2 / (f_m ** 2 - f0 ** 2) - mu_m
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if resid(lo) * resid(mid) <= 0:
            hi = mid
        else:
            lo = mid
    f0 = 0.5 * (lo + hi)
    return dict(f0=f0, F=1 - f0 ** 2 / f_mp ** 2, f_m=f_m, f_mp=f_mp, mu_m=mu_m)


def mu_lorentz(f, p):
    return 1 - p["F"] * f ** 2 / (f ** 2 - p["f0"] ** 2)


def zin_layers(f, mus, g_m, eps=EPS_EFF):
    """Normalized impedance at the left short of [air g | layers (each thickness a) | air g | short]."""
    k0 = 2 * np.pi * f / C0
    tg = np.tan(k0 * g_m)
    Z = 1j * tg                                   # right air section in front of the right short
    for mu in reversed(mus):                       # right-most layer first
        n = np.sqrt(mu * eps + 0j)
        z = n / eps
        t = np.tan(n * k0 * A)
        Z = z * (Z + 1j * z * t) / (z + 1j * Z * t)
    return (Z + 1j * tg) / (1 + 1j * Z * tg)


def zin_homog(f, mu, N, g_m, eps=EPS_EFF):
    k0 = 2 * np.pi * f / C0
    tg = np.tan(k0 * g_m)
    Z1 = 1j * tg
    n = np.sqrt(mu * eps + 0j)
    z = n / eps
    t = np.tan(n * k0 * N * A)
    Z2 = z * (Z1 + 1j * z * t) / (z + 1j * Z1 * t)
    return (Z2 + 1j * tg) / (1 + 1j * Z2 * tg)


def roots(fgrid, X):
    """Upward zero crossings of X = Im Z (series resonances at the short), linear interpolation."""
    out = []
    for i in range(len(fgrid) - 1):
        if np.isfinite(X[i]) and np.isfinite(X[i + 1]) and X[i] < 0 <= X[i + 1]:
            out.append(fgrid[i] - X[i] * (fgrid[i + 1] - fgrid[i]) / (X[i + 1] - X[i]))
    return out


def refine(fun, f_lo, f_hi, it=80):
    a, b = f_lo, f_hi
    fa = fun(a)
    for _ in range(it):
        m = 0.5 * (a + b)
        fm = fun(m)
        if fa * fm <= 0:
            b = m
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


def line_modes(cell_params, g_m=0.035, model="layer", fmin=120e6, fmax=200e6, n=80001):
    """Lossless resonances of the 4-cell line for a list of per-cell Lorentz parameter dicts."""
    f = np.linspace(fmin, fmax, n)
    mus = [mu_lorentz(f, p) for p in cell_params]
    if model == "layer":
        X = zin_layers(f, mus, g_m).imag
        fun = lambda x: zin_layers(np.array([x]), [mu_lorentz(np.array([x]), p) for p in cell_params], g_m).imag[0]
    else:
        mu = sum(mus) / len(mus)
        X = zin_homog(f, mu, len(cell_params), g_m).imag
        fun = lambda x: zin_homog(np.array([x]), sum(mu_lorentz(np.array([x]), p) for p in cell_params)
                                  / len(cell_params), len(cell_params), g_m).imag[0]
    out = []
    df = f[1] - f[0]
    for r in roots(f, X):
        # reject crossings at poles of X (Lorentz poles of a cell): require small |X| nearby
        fr = refine(fun, r - df, r + df)
        if abs(fun(fr)) < 1e-3:
            out.append(fr)
    # drop near-pole artefacts within 20 kHz of any cell's f0
    poles = [p["f0"] for p in cell_params]
    return [x for x in out if min(abs(x - q) for q in poles) > 2e4]


def mixture_fmp(cell_params):
    """Frequency where the volume-average Lorentz permeability crosses zero (above all f0)."""
    fun = lambda x: sum(mu_lorentz(x, p) for p in cell_params) / len(cell_params)
    lo = max(p["f0"] for p in cell_params) * (1 + 1e-6)
    return refine(fun, lo, 2.5e8)


def shift(p, tau_m, tau_mp, dT):
    return beig_params(p["f_m"] * (1 + tau_m * dT), p["f_mp"] * (1 + tau_mp * dT))
