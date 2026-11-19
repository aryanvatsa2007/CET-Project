"""Thermodynamic state from an equation of state plus an ideal-gas Cp model.

Reference state: ideal gas at T_REF = 298.15 K and P_REF = 1 bar has H = 0 and
S = 0.  A real-fluid property is the ideal-gas value at (T, P) plus the
residual (departure) at the same T and P:

    H(T,P) = Int_{T_REF}^{T} Cp0 dT                       + HR(T,P)
    S(T,P) = Int_{T_REF}^{T} Cp0/T dT - R ln(P/P_REF)     + SR(T,P)
"""
import math

from .constants import P_REF, R, T_REF
from .errors import InconsistentInputError, MissingDataError

# output order and SI units
STATE_UNITS = [
    ("T", "K"), ("P", "Pa"), ("V", "m3/mol"), ("Z", "1"), ("rho", "kg/m3"),
    ("alpha", "1/K"), ("kappaT", "1/Pa"), ("gammaV", "Pa/K"),
    ("lnphi", "1"), ("fugacity", "Pa"),
    ("HR", "J/mol"), ("SR", "J/(mol*K)"), ("GR", "J/mol"),
    ("Cp0", "J/(mol*K)"), ("Cp", "J/(mol*K)"), ("Cv", "J/(mol*K)"),
    ("gamma", "1"), ("kappaS", "1/Pa"), ("w", "m/s"), ("muJT", "K/Pa"),
    ("H", "J/mol"), ("S", "J/(mol*K)"), ("U", "J/mol"),
    ("A", "J/mol"), ("G", "J/mol"),
]
UNIT_OF = dict(STATE_UNITS)


def phase_label(eos, T, P, V):
    if eos.ideal:
        return "ideal gas"
    if T >= eos.Tc:
        return "supercritical" if P >= eos.Pc else "gas"
    try:
        ps = eos.psat(T)[0]
    except Exception:
        return "undetermined"
    return "liquid" if P > ps else "vapor"


def state(fluid, eos, T, P, phase="stable"):
    """Return a dict of properties at (T, P).

    Volumetric quantities need only the equation of state.  Caloric quantities
    also need the fluid's ideal-gas Cp model; if it is absent they are left out
    and ``out['missing']`` says what would be needed.
    """
    V = eos.volume(T, P, phase)
    dPdT, dPdV = eos.dPdT(T, V), eos.dPdV(T, V)
    if dPdV >= 0:
        raise InconsistentInputError(
            f"mechanically unstable root at T = {T:g} K, P = {P:g} Pa "
            "((dP/dV)_T >= 0)")
    Z = P * V / (R * T)
    res = eos.residual(T, V)
    out = {
        "T": T, "P": P, "V": V, "Z": Z,
        "alpha": -dPdT / (V * dPdV),
        "kappaT": -1.0 / (V * dPdV),
        "gammaV": dPdT,
        "HR": res["UR"] + P * V - R * T,
        "SR": res["SR"] + R * math.log(Z),
    }
    out["GR"] = out["HR"] - T * out["SR"]
    out["lnphi"] = out["GR"] / (R * T)
    out["fugacity"] = P * math.exp(out["lnphi"])
    out["phase"] = phase_label(eos, T, P, V)
    out["missing"] = []
    M = getattr(fluid, "M", None)
    if M:
        out["rho"] = M / V
    cp0 = getattr(fluid, "cp0", None)
    if cp0 is None:
        out["missing"] = ["Cp0(T)"]
        return out
    c0 = cp0.cp(T)
    Cv = c0 - R + res["CvR"]
    Cp = Cv - T * dPdT ** 2 / dPdV
    H = cp0.dh(cp0.tref, T) + out["HR"]
    S = cp0.ds(cp0.tref, T) - R * math.log(P / P_REF) + out["SR"]
    out.update({
        "Cp0": c0, "Cv": Cv, "Cp": Cp, "gamma": Cp / Cv,
        "kappaS": out["kappaT"] * Cv / Cp,
        "muJT": V * (T * out["alpha"] - 1.0) / Cp,
        "H": H, "S": S, "U": H - P * V,
        "A": H - P * V - T * S, "G": H - T * S,
    })
    if M:
        out["w"] = math.sqrt(V / (M * out["kappaS"]))
    else:
        out["missing"].append("M")
    return out


def require(st, names):
    """Raise MissingDataError if any of ``names`` is absent from a state dict."""
    lack = [n for n in names if n not in st]
    if lack:
        raise MissingDataError(
            "cannot compute " + ", ".join(lack) + "; missing input: "
            + ", ".join(st.get("missing", []) or ["unknown"]),
            st.get("missing", []))


def process(fluid, eos, T1, P1, T2, P2, phase1="stable", phase2="stable"):
    """Property changes for state 1 -> state 2 (state functions, path independent)."""
    s1 = state(fluid, eos, T1, P1, phase1)
    s2 = state(fluid, eos, T2, P2, phase2)
    require(s1, ["H", "S"])
    out = {"d" + k: s2[k] - s1[k] for k in ("H", "S", "U", "A", "G", "V")}
    # ideal-gas and residual contributions, reported separately
    out["dH_ideal"] = fluid.cp0.dh(T1, T2)
    out["dS_ideal"] = fluid.cp0.ds(T1, T2) - R * math.log(P2 / P1)
    out["dH_residual"] = s2["HR"] - s1["HR"]
    out["dS_residual"] = s2["SR"] - s1["SR"]
    return out, s1, s2
