"""Fluids, ideal-gas heat capacity models and pressure-explicit equations of state.

Cubic equations are written in the generic two-parameter form

    P = R T / (V - b)  -  a(T) / ((V + eps b)(V + sig b))

which covers van der Waals (vdW), Redlich-Kwong (RK), Soave-Redlich-Kwong
(SRK) and Peng-Robinson (PR).  All derivatives and residual (departure)
functions below are closed-form.
"""
import csv
import math
import os
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .constants import R, T_REF
from .errors import (ConvergenceError, InputError, MissingDataError,
                     OutOfRangeError)

_DATA = os.path.join(os.path.dirname(__file__), "data")


# --------------------------------------------------------------------------
# ideal-gas heat capacity
# --------------------------------------------------------------------------
class IdealGasCp:
    """Ideal-gas isobaric heat capacity Cp0(T) with a stated validity range.

    kind = 'shomate' : Cp0 = A + B t + C t^2 + D t^3 + E/t^2, t = T/1000 K
                       (J/(mol K); the form used by the NIST WebBook)
    kind = 'constant': Cp0 = A
    kind = 'table'   : cubic spline through tabulated (T, Cp0) points
    """

    def __init__(self, kind, coeffs=None, tmin=0.0, tmax=math.inf,
                 table=None, source="", extrapolate=False):
        self.kind, self.source = kind, source
        self.tmin, self.tmax = tmin, tmax
        self.extrapolate = extrapolate
        self.warnings = []
        if kind == "shomate":
            self.c = tuple(float(x) for x in coeffs)
        elif kind == "constant":
            self.c = (float(coeffs[0]),)
        elif kind == "table":
            from scipy.interpolate import CubicSpline
            T, cp = (np.asarray(v, float) for v in table)
            if len(T) < 4 or np.any(np.diff(T) <= 0):
                raise InputError("a Cp0 table needs at least four points with "
                                 "strictly increasing temperature")
            self.spline = CubicSpline(T, cp)
            self._hint = self.spline.antiderivative()
            self._sint = CubicSpline(T, cp / T).antiderivative()
            self.tmin, self.tmax = float(T[0]), float(T[-1])
        else:
            raise InputError(f"unknown Cp0 model '{kind}'")

    @property
    def tref(self):
        """Reference temperature for H and S: T_REF when the model covers it,
        otherwise the lowest temperature of the model."""
        return T_REF if self.tmin <= T_REF <= self.tmax else self.tmin

    def _check(self, T):
        if self.tmin <= T <= self.tmax:
            return
        msg = (f"ideal-gas Cp model ({self.kind}) is valid for "
               f"{self.tmin:g} K to {self.tmax:g} K; requested {T:g} K")
        if self.kind == "table" or not self.extrapolate:
            raise OutOfRangeError(msg)
        if msg not in self.warnings:
            self.warnings.append("EXTRAPOLATED: " + msg)

    def cp(self, T):
        self._check(T)
        if self.kind == "shomate":
            A, B, C, D, E = self.c
            t = T / 1000.0
            return A + B * t + C * t * t + D * t ** 3 + E / (t * t)
        if self.kind == "constant":
            return self.c[0]
        return float(self.spline(T))

    def dh(self, T0, T):
        """Integral of Cp0 dT from T0 to T, J/mol."""
        self._check(T0), self._check(T)
        if self.kind == "shomate":
            A, B, C, D, E = self.c

            def F(T):
                t = T / 1000.0
                return 1000.0 * (A * t + B * t ** 2 / 2 + C * t ** 3 / 3
                                 + D * t ** 4 / 4 - E / t)
            return F(T) - F(T0)
        if self.kind == "constant":
            return self.c[0] * (T - T0)
        return float(self._hint(T) - self._hint(T0))

    def ds(self, T0, T):
        """Integral of Cp0/T dT from T0 to T, J/(mol K)."""
        self._check(T0), self._check(T)
        if self.kind == "shomate":
            A, B, C, D, E = self.c

            def F(T):
                t = T / 1000.0
                return (A * math.log(t) + B * t + C * t ** 2 / 2
                        + D * t ** 3 / 3 - E / (2 * t * t))
            return F(T) - F(T0)
        if self.kind == "constant":
            return self.c[0] * math.log(T / T0)
        return float(self._sint(T) - self._sint(T0))


# --------------------------------------------------------------------------
# fluids
# --------------------------------------------------------------------------
@dataclass
class Fluid:
    name: str
    formula: str
    M: float              # kg/mol
    Tc: float             # K
    Pc: float             # Pa
    omega: float
    cp0: IdealGasCp = None
    source_critical: str = ""
    source_cp: str = ""

    def require_cp0(self):
        if self.cp0 is None:
            raise MissingDataError(
                f"no ideal-gas heat capacity model is available for {self.name}; "
                "caloric properties (H, S, U, A, G, Cp, Cv, w, muJT) cannot be "
                "computed. Supply cp0_shomate, cp0_constant or cp0_table.",
                ["Cp0(T)"])
        return self.cp0


def load_fluids(path=None):
    """Read the fluid database (CSV) and return {name: Fluid}."""
    path = path or os.path.join(_DATA, "fluids.csv")
    out = {}
    with open(path, newline="") as fh:
        rows = csv.DictReader(line for line in fh if not line.startswith("#"))
        for r in rows:
            cp0 = None
            if r["cpA"].strip():
                cp0 = IdealGasCp(
                    "shomate", [r[k] for k in ("cpA", "cpB", "cpC", "cpD", "cpE")],
                    float(r["cp_Tmin_K"]), float(r["cp_Tmax_K"]),
                    source=r["source_cp"])
            out[r["name"]] = Fluid(
                r["name"], r["formula"], float(r["M_kg_per_mol"]),
                float(r["Tc_K"]), float(r["Pc_Pa"]), float(r["omega"]),
                cp0, r["source_critical"], r["source_cp"])
    return out


def get_fluid(name, path=None):
    fluids = load_fluids(path)
    for key, f in fluids.items():
        if name.lower() in (key.lower(), f.formula.lower()):
            return f
    raise MissingDataError(
        f"fluid '{name}' is not in the database (available: "
        + ", ".join(sorted(fluids)) + "). Supply Tc, Pc, omega and M directly.",
        ["Tc", "Pc", "omega"])


# --------------------------------------------------------------------------
# equations of state
# --------------------------------------------------------------------------
_SQ2 = math.sqrt(2.0)
_CBRT2 = 2.0 ** (1.0 / 3.0)
_KINDS = {
    #        eps        sig        Omega_a                     Omega_b
    "vdW": (0.0,       0.0,       27.0 / 64.0,                1.0 / 8.0),
    "RK":  (0.0,       1.0,       1.0 / (9.0 * (_CBRT2 - 1)), (_CBRT2 - 1) / 3.0),
    "SRK": (0.0,       1.0,       1.0 / (9.0 * (_CBRT2 - 1)), (_CBRT2 - 1) / 3.0),
    "PR":  (1 - _SQ2,  1 + _SQ2,  0.45724,                    0.07780),
}
EOS_NAMES = ("IG",) + tuple(_KINDS)


class EoS:
    """Ideal gas ('IG') or cubic ('vdW', 'RK', 'SRK', 'PR') equation of state."""

    def __init__(self, kind, Tc=None, Pc=None, omega=None):
        if kind not in EOS_NAMES:
            raise InputError(f"unknown equation of state '{kind}'; choose one "
                             f"of {', '.join(EOS_NAMES)}")
        self.kind = kind
        self.ideal = kind == "IG"
        if self.ideal:
            self.b = 0.0
            return
        missing = [n for n, v in (("Tc", Tc), ("Pc", Pc)) if v is None]
        if kind in ("SRK", "PR") and omega is None:
            missing.append("omega")
        if missing:
            raise MissingDataError(
                f"the {kind} equation needs " + ", ".join(missing), missing)
        self.Tc, self.Pc, self.omega = Tc, Pc, omega
        self.eps, self.sig, Oa, Ob = _KINDS[kind]
        self.ac = Oa * R * R * Tc * Tc / Pc
        self.b = Ob * R * Tc / Pc
        if kind == "SRK":
            self.m = 0.480 + 1.574 * omega - 0.176 * omega ** 2
        elif kind == "PR":
            self.m = 0.37464 + 1.54226 * omega - 0.26992 * omega ** 2

    # ---- a(T) and its temperature derivatives --------------------------
    def a(self, T):
        """Return (a, da/dT, d2a/dT2)."""
        if self.ideal:
            return 0.0, 0.0, 0.0
        if self.kind == "vdW":
            return self.ac, 0.0, 0.0
        Tr = T / self.Tc
        if self.kind == "RK":
            return (self.ac * Tr ** -0.5,
                    -0.5 * self.ac * Tr ** -1.5 / self.Tc,
                    0.75 * self.ac * Tr ** -2.5 / self.Tc ** 2)
        s, m = math.sqrt(Tr), self.m
        g = 1.0 + m * (1.0 - s)
        return (self.ac * g * g,
                -self.ac * m * g / (s * self.Tc),
                self.ac * m * (1.0 + m) / (2.0 * self.Tc ** 2 * s ** 3))

    # ---- pressure and its derivatives ----------------------------------
    def _den(self, V):
        return (V + self.eps * self.b) * (V + self.sig * self.b)

    def pressure(self, T, V):
        if self.ideal:
            return R * T / V
        return R * T / (V - self.b) - self.a(T)[0] / self._den(V)

    def dPdT(self, T, V):
        if self.ideal:
            return R / V
        return R / (V - self.b) - self.a(T)[1] / self._den(V)

    def dPdV(self, T, V):
        if self.ideal:
            return -R * T / (V * V)
        b = self.b
        return (-R * T / (V - b) ** 2
                + self.a(T)[0] * (2 * V + (self.eps + self.sig) * b) / self._den(V) ** 2)

    # ---- volume roots ---------------------------------------------------
    def volumes(self, T, P):
        """All physically admissible molar volumes (V > b) at (T, P), ascending."""
        if T <= 0 or P <= 0:
            raise InputError("temperature and pressure must be positive "
                             f"(got T = {T:g} K, P = {P:g} Pa)")
        if self.ideal:
            return [R * T / P]
        A = self.a(T)[0] * P / (R * T) ** 2
        B = self.b * P / (R * T)
        e, s = self.eps, self.sig
        coeff = [1.0,
                 (e + s - 1.0) * B - 1.0,
                 e * s * B * B - (e + s) * (B * B + B) + A,
                 -(e * s * (B ** 3 + B * B) + A * B)]
        roots = np.roots(coeff)
        out = []
        for z in roots:
            if abs(z.imag) > 1e-8 * max(1.0, abs(z.real)) or z.real <= B:
                continue
            V = z.real * R * T / P
            for _ in range(50):                      # Newton polish on P(T,V) = P
                dV = (self.pressure(T, V) - P) / self.dPdV(T, V)
                if not (V - dV > self.b):
                    break
                V -= dV
                if abs(dV) < 1e-15 * V:
                    break
            if V > self.b and abs(self.pressure(T, V) - P) < 1e-8 * P:
                out.append(V)
        out.sort()
        uniq = []
        for V in out:
            if not uniq or abs(V - uniq[-1]) > 1e-9 * V:
                uniq.append(V)
        if not uniq:
            raise ConvergenceError(f"no volume root found at T = {T:g} K, P = {P:g} Pa")
        return uniq

    def volume(self, T, P, phase="stable"):
        """Molar volume for the requested root: 'stable', 'liquid' or 'vapor'."""
        vols = self.volumes(T, P)
        if phase == "liquid":
            return vols[0]
        if phase == "vapor":
            return vols[-1]
        if phase != "stable":
            raise InputError(f"phase must be stable, liquid or vapor (got '{phase}')")
        return min(vols, key=lambda V: self.lnphi(T, V))

    # ---- residual functions --------------------------------------------
    def _q(self, V):
        """Integral of dV / ((V + eps b)(V + sig b)) from V to infinity."""
        if self.kind == "vdW":
            return 1.0 / V
        b, d = self.b, self.sig - self.eps
        return math.log((V + self.sig * b) / (V + self.eps * b)) / (b * d)

    def residual(self, T, V):
        """Residual properties at the same T and V as the ideal gas.

        Returns AR, SR, UR (J/mol or J/(mol K)) and CvR.
        """
        if self.ideal:
            return {"AR": 0.0, "SR": 0.0, "UR": 0.0, "CvR": 0.0}
        a, da, d2a = self.a(T)
        q = self._q(V)
        ln = math.log(1.0 - self.b / V)
        return {"AR": -R * T * ln - a * q,
                "SR": R * ln + da * q,
                "UR": -(a - T * da) * q,
                "CvR": T * d2a * q}

    def lnphi(self, T, V):
        """ln of the fugacity coefficient (= GR/RT at the same T and P)."""
        if self.ideal:
            return 0.0
        Z = self.pressure(T, V) * V / (R * T)
        return self.residual(T, V)["AR"] / (R * T) + Z - 1.0 - math.log(Z)

    # ---- saturation -----------------------------------------------------
    def spinodal(self, T):
        """Volumes where (dP/dV)_T = 0 on a subcritical isotherm (liquid, vapor)."""
        grid = self.b * (1.0 + np.logspace(-4, 4, 4000))
        vals = np.array([self.dPdV(T, V) for V in grid])
        idx = np.nonzero(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
        if len(idx) < 2:
            return None
        f = lambda V: self.dPdV(T, V)
        return (brentq(f, grid[idx[0]], grid[idx[0] + 1]),
                brentq(f, grid[idx[-1]], grid[idx[-1] + 1]))

    def psat(self, T):
        """Saturation pressure from equality of fugacity; returns (P, VL, VV)."""
        if self.ideal:
            raise MissingDataError("the ideal-gas equation has no vapour-liquid "
                                   "transition", ["a real-fluid equation of state"])
        if T >= self.Tc:
            raise OutOfRangeError(
                f"no saturation state above the critical temperature "
                f"({T:g} K >= {self.Tc:g} K)")
        sp = self.spinodal(T)
        if sp is None:
            raise ConvergenceError(f"no van der Waals loop found at {T:g} K")
        p_hi = self.pressure(T, sp[1])
        p_lo = max(self.pressure(T, sp[0]), 1e-10 * p_hi)
        vmin = self.b * (1.0 + 1e-13)

        def roots(P):
            # liquid root left of the liquid spinodal, vapour root right of the
            # vapour spinodal; bracketing keeps both even at very low pressure
            f = lambda V: self.pressure(T, V) - P
            return (brentq(f, vmin, sp[0], xtol=1e-18, rtol=1e-14),
                    brentq(f, sp[1], 10.0 * R * T / P + sp[1], xtol=1e-18, rtol=1e-14))

        def gap(P):
            VL, VV = roots(P)
            return self.lnphi(T, VL) - self.lnphi(T, VV)

        lo = p_lo + 1e-9 * (p_hi - p_lo)
        hi = p_hi - 1e-9 * (p_hi - p_lo)
        try:
            P = brentq(gap, lo, hi, xtol=1e-14 * hi, rtol=1e-13)
        except ValueError as exc:
            raise ConvergenceError(f"saturation pressure solve failed at {T:g} K") from exc
        VL, VV = roots(P)
        return P, VL, VV

    def describe(self):
        if self.ideal:
            return "ideal gas, P = RT/V"
        txt = (f"{self.kind}: a_c = {self.ac:.6g} Pa m6/mol2, "
               f"b = {self.b:.6g} m3/mol")
        if hasattr(self, "m"):
            txt += f", m = {self.m:.6g}"
        return txt
