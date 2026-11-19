"""Computability: what can be obtained from the data actually supplied.

A state point is described by a set of quantities.  Each one known to the
engine carries a provenance:

    stated         given directly by the user
    tabulated      read from a user table
    parameterised  evaluated from a model with parameters (EoS, Cp0 correlation)
    inferred       derived from other quantities through an exact relation

Exact thermodynamic identities are stored as rules (output <- inputs).  Forward
chaining applies every rule whose inputs are known until nothing new appears.
A rule that reproduces an already-known quantity is a consistency check: a
disagreement beyond the tolerance is reported as an inconsistency.  If a target
stays unknown, a backward search over the same rules lists the smallest sets of
additional quantities that would make it computable.
"""
import math
from dataclasses import dataclass, field

from .constants import R
from .errors import InconsistentInputError, MissingDataError
from .symbolic import SYMS, derivative, parse_derivative

QUANTITIES = ("T", "P", "V", "Z", "M", "rho", "S", "Cp", "Cv", "gamma",
              "alpha", "kappaT", "kappaS", "gammaV", "w", "muJT")


@dataclass
class Rule:
    name: str
    out: str
    ins: tuple
    func: callable
    relation: str


def _rules():
    r = []

    def add(name, out, ins, func, relation):
        r.append(Rule(name, out, tuple(ins), func, relation))

    # compressibility factor
    add("Z-def", "Z", "PVT", lambda P, V, T: P * V / (R * T), "Z = P V / (R T)")
    add("Z-def", "V", "ZTP", lambda Z, T, P: Z * R * T / P, "V = Z R T / P")
    add("Z-def", "P", "ZTV", lambda Z, T, V: Z * R * T / V, "P = Z R T / V")
    add("Z-def", "T", "PVZ", lambda P, V, Z: P * V / (Z * R), "T = P V / (Z R)")
    # mass density
    add("density", "rho", ("M", "V"), lambda M, V: M / V, "rho = M / V")
    add("density", "V", ("M", "rho"), lambda M, rho: M / rho, "V = M / rho")
    # thermal pressure coefficient (cyclic rule)
    add("cyclic", "gammaV", ("alpha", "kappaT"), lambda a, k: a / k,
        "(dP/dT)_V = alpha / kappaT")
    add("cyclic", "alpha", ("gammaV", "kappaT"), lambda g, k: g * k,
        "alpha = (dP/dT)_V kappaT")
    add("cyclic", "kappaT", ("alpha", "gammaV"), lambda a, g: a / g,
        "kappaT = alpha / (dP/dT)_V")
    # heat capacity difference
    add("Cp-Cv", "Cv", ("Cp", "T", "V", "alpha", "kappaT"),
        lambda Cp, T, V, a, k: Cp - T * V * a * a / k,
        "Cv = Cp - T V alpha^2 / kappaT")
    add("Cp-Cv", "Cp", ("Cv", "T", "V", "alpha", "kappaT"),
        lambda Cv, T, V, a, k: Cv + T * V * a * a / k,
        "Cp = Cv + T V alpha^2 / kappaT")
    # heat capacity ratio
    add("gamma-def", "gamma", ("Cp", "Cv"), lambda Cp, Cv: Cp / Cv, "gamma = Cp / Cv")
    add("gamma-def", "Cp", ("gamma", "Cv"), lambda g, Cv: g * Cv, "Cp = gamma Cv")
    add("gamma-def", "Cv", ("Cp", "gamma"), lambda Cp, g: Cp / g, "Cv = Cp / gamma")
    # compressibility ratio
    add("kappa-ratio", "kappaS", ("kappaT", "gamma"), lambda k, g: k / g,
        "kappaS = kappaT / gamma")
    add("kappa-ratio", "kappaT", ("kappaS", "gamma"), lambda k, g: k * g,
        "kappaT = gamma kappaS")
    add("kappa-ratio", "gamma", ("kappaT", "kappaS"), lambda kt, ks: kt / ks,
        "gamma = kappaT / kappaS")
    add("kappaS-direct", "kappaS", ("kappaT", "T", "V", "alpha", "Cp"),
        lambda k, T, V, a, Cp: k - T * V * a * a / Cp,
        "kappaS = kappaT - T V alpha^2 / Cp")
    add("kappaS-direct", "kappaT", ("kappaS", "T", "V", "alpha", "Cp"),
        lambda k, T, V, a, Cp: k + T * V * a * a / Cp,
        "kappaT = kappaS + T V alpha^2 / Cp")
    # speed of sound
    add("sound", "w", ("V", "M", "kappaS"), lambda V, M, k: math.sqrt(V / (M * k)),
        "w = sqrt(V / (M kappaS))")
    add("sound", "kappaS", ("V", "M", "w"), lambda V, M, w: V / (M * w * w),
        "kappaS = V / (M w^2)")
    # Joule-Thomson coefficient
    add("Joule-Thomson", "muJT", ("V", "T", "alpha", "Cp"),
        lambda V, T, a, Cp: V * (T * a - 1.0) / Cp, "muJT = V (T alpha - 1) / Cp")
    add("Joule-Thomson", "alpha", ("muJT", "Cp", "V", "T"),
        lambda mu, Cp, V, T: (mu * Cp / V + 1.0) / T,
        "alpha = (muJT Cp / V + 1) / T")
    return r


RULES = _rules()


@dataclass
class Fact:
    value: float
    provenance: str
    how: str = ""


@dataclass
class Report:
    facts: dict = field(default_factory=dict)
    checks: list = field(default_factory=list)       # (relation, quantity, known, derived, rel diff, ok)
    violations: list = field(default_factory=list)   # text

    @property
    def consistent(self):
        return not self.violations and all(c[5] for c in self.checks)


def infer(known, provenance=None, rtol=1e-2):
    """Forward-chain from ``known`` ({name: SI value}).

    ``rtol`` is the relative disagreement above which two routes to the same
    quantity are declared inconsistent.
    """
    provenance = provenance or {}
    rep = Report()
    for k, v in known.items():
        rep.facts[k] = Fact(float(v), provenance.get(k, "stated"), "input")
    done = set()
    changed = True
    while changed:
        changed = False
        for i, rule in enumerate(RULES):
            if i in done or not all(n in rep.facts for n in rule.ins):
                continue
            done.add(i)
            try:
                val = rule.func(*(rep.facts[n].value for n in rule.ins))
            except (ZeroDivisionError, ValueError):
                rep.violations.append(
                    f"{rule.relation} cannot be evaluated with the supplied values "
                    "(division by zero or square root of a negative number)")
                continue
            if not math.isfinite(val):
                continue
            if rule.out not in rep.facts:
                rep.facts[rule.out] = Fact(
                    val, "inferred", f"{rule.relation}  [from {', '.join(rule.ins)}]")
                changed = True
            else:
                old = rep.facts[rule.out].value
                diff = abs(val - old) / max(abs(old), 1e-300)
                if diff > 1e-9:                      # not a mere round trip
                    rep.checks.append((rule.relation, rule.out, old, val,
                                       diff, diff <= rtol))
    f = rep.facts
    for name in ("T", "P", "V", "Cp", "Cv", "M"):
        if name in f and f[name].value <= 0:
            rep.violations.append(f"{name} = {f[name].value:g} must be positive")
    for name, why in (("kappaT", "mechanical stability"), ("kappaS", "mechanical stability")):
        if name in f and f[name].value <= 0:
            rep.violations.append(f"{name} = {f[name].value:g} violates {why} ({name} > 0)")
    if "Cp" in f and "Cv" in f and f["Cp"].value < f["Cv"].value * (1 - 1e-12):
        rep.violations.append(
            f"Cp = {f['Cp'].value:g} is smaller than Cv = {f['Cv'].value:g}; "
            "Cp - Cv = T V alpha^2 / kappaT cannot be negative")
    return rep


def alternatives(target, known, max_size=3, limit=6):
    """Smallest sets of extra quantities that would make ``target`` computable."""
    known = set(known)

    def need(q, visiting):
        if q in known:
            return {frozenset()}
        opts = {frozenset([q])}                      # supply it directly
        if q in visiting:
            return opts
        for rule in RULES:
            if rule.out != q:
                continue
            combos = {frozenset()}
            for inp in rule.ins:
                sub = need(inp, visiting | {q})
                combos = {a | b for a in combos for b in sub if len(a | b) <= max_size}
                if not combos:
                    break
            opts |= combos
        return {o for o in opts if not any(p < o for p in opts)}

    opts = need(target, frozenset())
    opts.discard(frozenset([target]))
    return [sorted(o) for o in sorted(opts, key=lambda s: (len(s), sorted(s)))][:limit]


def require(rep, target):
    """Value of ``target`` from a Report, or MissingDataError with alternatives."""
    if target in rep.facts:
        return rep.facts[target].value
    alts = alternatives(target, rep.facts)
    msg = f"{target} cannot be computed from the supplied data"
    if alts:
        msg += "; supply any one of: " + " | ".join("{" + ", ".join(a) + "}" for a in alts)
    raise MissingDataError(msg, alts[0] if alts else [target])


def derivative_from_facts(text, rep):
    """Evaluate (dX/dY)_Z from a Report; returns (Ratio, value)."""
    expr = derivative(*parse_derivative(text))
    lack = sorted(s for s in expr.symbols() if s not in rep.facts)
    if lack:
        parts = []
        for s in lack:
            alts = alternatives(s, rep.facts)
            extra = (" (or " + " | ".join("{" + ", ".join(a) + "}" for a in alts[:3]) + ")"
                     if alts else "")
            parts.append(s + extra)
        raise MissingDataError(
            f"{text} = {expr} needs " + ", ".join(sorted(expr.symbols()))
            + "; missing: " + "; ".join(parts), lack)
    return expr, expr.eval({s: rep.facts[s].value for s in SYMS if s in rep.facts})
