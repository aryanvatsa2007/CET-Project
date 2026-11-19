"""Generic numerical scheme for any pressure-explicit model P(T, V).

Nothing here knows the functional form of the equation of state.  Derivatives
are central differences with one Richardson extrapolation step (error O(h^4));
residual functions are adaptive Gauss-Kronrod quadratures in molar density:

    AR(T,V)/RT = Int_0^rho (Z - 1)/rho' drho'            (same T and V)
    SR(T,V)    = -(dAR/dT)_V
    CvR(T,V)   = -T (d2AR/dT2)_V

The closed-form cubic results in :mod:`eos` are verified against this scheme,
and the scheme is what makes a user-supplied P(T, V) function usable.
"""
import math

from scipy.integrate import quad

from .constants import P_REF, R, T_REF
from .errors import ConvergenceError

REL_STEP_1 = 1e-3    # relative step, first derivatives
REL_STEP_2 = 1e-2    # relative step, second derivatives
QUAD_TOL = 1e-12     # relative tolerance of the quadrature


def d1(f, x, rel=REL_STEP_1):
    """First derivative, central difference + Richardson, O(h^4)."""
    h = rel * abs(x)
    c = lambda s: (f(x + s) - f(x - s)) / (2.0 * s)
    return (4.0 * c(h / 2.0) - c(h)) / 3.0


def d2(f, x, rel=REL_STEP_2):
    """Second derivative, central difference + Richardson, O(h^4)."""
    h = rel * abs(x)
    f0 = f(x)
    c = lambda s: (f(x + s) - 2.0 * f0 + f(x - s)) / (s * s)
    return (4.0 * c(h / 2.0) - c(h)) / 3.0


class NumericalModel:
    """Properties of a fluid described only by a callable ``pressure(T, V)``."""

    def __init__(self, pressure, cp0=None, M=None):
        self.pressure, self.cp0, self.M = pressure, cp0, M

    def dPdT(self, T, V):
        return d1(lambda t: self.pressure(t, V), T)

    def dPdV(self, T, V):
        return d1(lambda v: self.pressure(T, v), V)

    def volume(self, T, P, guess=None):
        """Solve P(T, V) = P by a damped Newton iteration from ``guess``."""
        V = guess if guess else R * T / P
        for _ in range(200):
            f = self.pressure(T, V) - P
            step = f / self.dPdV(T, V)
            while V - step <= 0:
                step *= 0.5
            V -= step
            if abs(step) < 1e-13 * V:
                return V
        raise ConvergenceError(f"volume solve failed at T = {T:g} K, P = {P:g} Pa")

    def AR(self, T, V):
        rho = 1.0 / V
        g = lambda r: (self.pressure(T, 1.0 / r) / (r * R * T) - 1.0) / r
        val, _ = quad(g, 0.0, rho, epsrel=QUAD_TOL, epsabs=0.0, limit=200)
        return R * T * val

    def residual(self, T, V):
        ar = lambda t: self.AR(t, V)
        AR = ar(T)
        SR = -d1(ar, T)
        return {"AR": AR, "SR": SR, "UR": AR + T * SR, "CvR": -T * d2(ar, T)}

    def state(self, T, P, guess=None):
        V = self.volume(T, P, guess)
        dPdT, dPdV = self.dPdT(T, V), self.dPdV(T, V)
        res = self.residual(T, V)
        Z = P * V / (R * T)
        out = {"T": T, "P": P, "V": V, "Z": Z,
               "alpha": -dPdT / (V * dPdV), "kappaT": -1.0 / (V * dPdV),
               "gammaV": dPdT,
               "HR": res["UR"] + P * V - R * T,
               "SR": res["SR"] + R * math.log(Z)}
        out["GR"] = out["HR"] - T * out["SR"]
        out["lnphi"] = out["GR"] / (R * T)
        if self.cp0 is not None:
            c0 = self.cp0.cp(T)
            Cv = c0 - R + res["CvR"]
            Cp = Cv - T * dPdT ** 2 / dPdV
            H = self.cp0.dh(self.cp0.tref, T) + out["HR"]
            S = self.cp0.ds(self.cp0.tref, T) - R * math.log(P / P_REF) + out["SR"]
            out.update({"Cp0": c0, "Cv": Cv, "Cp": Cp, "gamma": Cp / Cv,
                        "kappaS": out["kappaT"] * Cv / Cp,
                        "muJT": V * (T * out["alpha"] - 1.0) / Cp,
                        "H": H, "S": S, "U": H - P * V,
                        "A": H - P * V - T * S, "G": H - T * S})
            if self.M:
                out["w"] = math.sqrt(V / (self.M * out["kappaS"]))
        return out
