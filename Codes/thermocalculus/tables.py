"""Tabulated input: P-V-T data on a rectangular (T, P) grid and Cp0(T) tables.

The compressibility factor Z(T, P) is fitted with a bicubic tensor-product
spline.  Derivatives come from the spline analytically, and residual functions
from the pressure-explicit integrals

    HR(T,P) = -R T^2 Int_0^P (dZ/dT)_P dP/P
    SR(T,P) = -R Int_0^P [ (Z - 1) + T (dZ/dT)_P ] dP/P

A table never reaches P = 0.  Below the lowest tabulated pressure P_min the
fluid is assumed to follow a second-virial form, Z - 1 proportional to P, so
that each integrand is constant on [0, P_min].  The size of that correction is
reported so its importance can be judged.

There is NO extrapolation: a request outside the grid raises OutOfRangeError.
"""
import csv
import math
import re

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import RectBivariateSpline

from .constants import P_REF, R, T_REF
from .errors import (InconsistentInputError, InputError, OutOfRangeError,
                     UnitError)
from .units import VAR_DIMS, dims_str, parse_unit

_HEADER = re.compile(r"^\s*([A-Za-z0-9]+)\s*\[(.*)\]\s*$")
_TABLE_VARS = {"T": "T", "P": "P", "V": "V", "Cp0": "Cp"}


def read_table_csv(path):
    """Read a CSV whose header cells look like ``T [K]`` or ``P [MPa]``.

    Lines starting with '#' are comments.  Returns {variable: SI array}.
    """
    with open(path, newline="") as fh:
        rows = list(csv.reader(line for line in fh if not line.lstrip().startswith("#")))
    rows = [r for r in rows if r]
    if len(rows) < 2:
        raise InputError(f"table '{path}' has no data rows")
    factors, names = [], []
    for cell in rows[0]:
        m = _HEADER.match(cell)
        if not m:
            raise InputError(
                f"table '{path}': column header '{cell}' must be written as "
                "'name [unit]', for example 'P [MPa]'")
        name, unit = m.group(1), m.group(2)
        if name not in _TABLE_VARS:
            raise InputError(f"table '{path}': unknown column '{name}' "
                             f"(allowed: {', '.join(_TABLE_VARS)})")
        factor, dims = parse_unit(unit)
        want = VAR_DIMS[_TABLE_VARS[name]]
        if tuple(dims) != tuple(want):
            raise UnitError(
                f"table '{path}': column {name} is given in '{unit}' "
                f"({dims_str(dims)}); expected {dims_str(want)}")
        names.append(name)
        factors.append(factor)
    try:
        data = np.array([[float(x) for x in r] for r in rows[1:]], float)
    except ValueError as exc:
        raise InputError(f"table '{path}': non-numeric or ragged entry ({exc})") from exc
    return {n: data[:, k] * f for k, (n, f) in enumerate(zip(names, factors))}


class PVTTable:
    """Z(T, P) spline through tabulated molar volumes on a full (T, P) grid."""

    def __init__(self, T, P, V, source=""):
        T, P, V = (np.asarray(x, float) for x in (T, P, V))
        self.Tg, self.Pg = np.unique(T), np.unique(P)
        nT, nP = len(self.Tg), len(self.Pg)
        if nT < 4 or nP < 4:
            raise InputError("a P-V-T table needs at least four temperatures "
                             f"and four pressures (got {nT} and {nP})")
        if len(T) != nT * nP:
            raise InputError(
                f"the P-V-T table is not a full rectangular grid: {nT} "
                f"temperatures x {nP} pressures needs {nT * nP} rows, got {len(T)}")
        if np.any(V <= 0):
            raise InputError("the P-V-T table contains non-positive volumes")
        Zg = np.full((nT, nP), np.nan)
        iT, iP = np.searchsorted(self.Tg, T), np.searchsorted(self.Pg, P)
        Zg[iT, iP] = P * V / (R * T)
        if np.isnan(Zg).any():
            raise InputError("the P-V-T table has duplicate or missing grid points")
        self.spl = RectBivariateSpline(self.Tg, self.Pg, Zg, kx=3, ky=3, s=0)
        self.source = source

    @classmethod
    def from_csv(cls, path):
        cols = read_table_csv(path)
        lack = [c for c in ("T", "P", "V") if c not in cols]
        if lack:
            raise InputError(f"table '{path}' lacks column(s) {', '.join(lack)}")
        return cls(cols["T"], cols["P"], cols["V"], source=path)

    def check(self, T, P):
        if not (self.Tg[0] <= T <= self.Tg[-1]) or not (self.Pg[0] <= P <= self.Pg[-1]):
            raise OutOfRangeError(
                f"state T = {T:g} K, P = {P:g} Pa is outside the table "
                f"(T from {self.Tg[0]:g} to {self.Tg[-1]:g} K, P from "
                f"{self.Pg[0]:g} to {self.Pg[-1]:g} Pa); the engine does not "
                "extrapolate tabulated data")

    def Z(self, T, P, dT=0, dP=0):
        return float(self.spl(T, P, dx=dT, dy=dP)[0, 0])

    def _integrals(self, T, P):
        """I1 = Int_0^P Z_T dP/P,  I2 = Int_0^P (Z-1) dP/P,  I1_T = dI1/dT,
        each split into (low-pressure correction, tabulated part)."""
        p0 = self.Pg[0]
        pts = [p for p in self.Pg if p0 < p < P]
        kw = dict(epsrel=1e-10, epsabs=0.0, limit=200)
        if pts:
            kw["points"] = pts

        def integ(f):
            return quad(f, p0, P, **kw)[0] if P > p0 else 0.0

        low = (self.Z(T, p0, dT=1), self.Z(T, p0) - 1.0, self.Z(T, p0, dT=2))
        tab = (integ(lambda p: self.Z(T, p, dT=1) / p),
               integ(lambda p: (self.Z(T, p) - 1.0) / p),
               integ(lambda p: self.Z(T, p, dT=2) / p))
        return low, tab

    def state(self, T, P, cp0=None, M=None):
        self.check(T, P)
        Z, ZT, ZP = self.Z(T, P), self.Z(T, P, dT=1), self.Z(T, P, dP=1)
        V = Z * R * T / P
        alpha = 1.0 / T + ZT / Z
        kappaT = 1.0 / P - ZP / Z
        if kappaT <= 0:
            raise InconsistentInputError(
                f"the table implies kappaT <= 0 at T = {T:g} K, P = {P:g} Pa "
                "(volume does not decrease with pressure)")
        low, tab = self._integrals(T, P)
        I1, I2, I1T = (low[k] + tab[k] for k in range(3))
        HR = -R * T * T * I1
        SR = -R * (I2 + T * I1)
        out = {"T": T, "P": P, "V": V, "Z": Z, "alpha": alpha, "kappaT": kappaT,
               "gammaV": alpha / kappaT, "HR": HR, "SR": SR, "GR": HR - T * SR,
               "HR_lowP_correction": -R * T * T * low[0],
               "SR_lowP_correction": -R * (low[1] + T * low[0])}
        out["lnphi"] = out["GR"] / (R * T)
        out["missing"] = []
        if M:
            out["rho"] = M / V
        if cp0 is None:
            out["missing"] = ["Cp0(T)"]
            return out
        c0 = cp0.cp(T)
        Cp = c0 - 2.0 * R * T * I1 - R * T * T * I1T
        Cv = Cp - T * V * alpha ** 2 / kappaT
        H = cp0.dh(cp0.tref, T) + HR
        S = cp0.ds(cp0.tref, T) - R * math.log(P / P_REF) + SR
        out.update({"Cp0": c0, "Cp": Cp, "Cv": Cv, "gamma": Cp / Cv,
                    "kappaS": kappaT * Cv / Cp,
                    "muJT": V * (T * alpha - 1.0) / Cp,
                    "H": H, "S": S, "U": H - P * V,
                    "A": H - P * V - T * S, "G": H - T * S})
        out["warnings"] = []
        if Cv <= 0 or out["kappaS"] <= 0:
            out["warnings"].append(
                f"the spline through the table gives Cv = {Cv:.4g} J/(mol K) at "
                f"T = {T:g} K, P = {P:g} Pa, which is unphysical: the grid is too "
                "coarse for second-derivative properties here (near-critical region)")
        elif M:
            out["w"] = math.sqrt(V / (M * out["kappaS"]))
        return out
