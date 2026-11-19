"""Symbolic reduction of thermodynamic partial derivatives.

Any first derivative (dX/dY)_Z with X, Y, Z drawn from
{T, P, V, S, U, H, A, G} is reduced to an expression in the state variables
T, P, V, S and the three measurable response functions

    Cp     = T (dS/dT)_P          isobaric heat capacity
    alpha  = (1/V)(dV/dT)_P       isobaric expansivity
    kappaT = -(1/V)(dV/dP)_T      isothermal compressibility

using the Jacobian identity (dX/dY)_Z = J(X,Z)/J(Y,Z), with
J(X,Z) = X_T Z_P - X_P Z_T taken in the (T, P) plane.

The algebra is exact: expressions are ratios of Laurent polynomials with
rational coefficients, implemented here so the engine needs no computer
algebra package.
"""
import re
from fractions import Fraction
from functools import reduce
from math import gcd

from .errors import InputError, MissingDataError, SingularDerivativeError

SYMS = ("T", "P", "V", "S", "Cp", "alpha", "kappaT")
VARIABLES = ("T", "P", "V", "S", "U", "H", "A", "G")
_N = len(SYMS)
_IDX = {s: i for i, s in enumerate(SYMS)}
_LATEX = {"T": "T", "P": "P", "V": "V", "S": "S", "Cp": "C_P",
          "alpha": r"\alpha", "kappaT": r"\kappa_T"}


class Poly:
    """Laurent polynomial in SYMS with Fraction coefficients."""

    __slots__ = ("t",)

    def __init__(self, terms=None):
        self.t = {e: Fraction(c) for e, c in (terms or {}).items() if c != 0}

    @staticmethod
    def const(c):
        return Poly({(0,) * _N: c})

    @staticmethod
    def mono(c=1, **powers):
        e = [0] * _N
        for name, p in powers.items():
            e[_IDX[name]] = p
        return Poly({tuple(e): c})

    def is_zero(self):
        return not self.t

    def __add__(self, other):
        out = dict(self.t)
        for e, c in other.t.items():
            out[e] = out.get(e, 0) + c
        return Poly(out)

    def __neg__(self):
        return Poly({e: -c for e, c in self.t.items()})

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if not isinstance(other, Poly):
            return Poly({e: c * other for e, c in self.t.items()})
        out = {}
        for e1, c1 in self.t.items():
            for e2, c2 in other.t.items():
                e = tuple(a + b for a, b in zip(e1, e2))
                out[e] = out.get(e, 0) + c1 * c2
        return Poly(out)

    def __eq__(self, other):
        return self.t == other.t

    def shift(self, exps):
        """Divide by the monomial with exponents ``exps``."""
        return Poly({tuple(a - b for a, b in zip(e, exps)): c
                     for e, c in self.t.items()})

    def eval(self, values):
        total = 0.0
        for e, c in self.t.items():
            term = float(c)
            for s, p in zip(SYMS, e):
                if p:
                    term *= values[s] ** p
            total += term
        return total

    def symbols(self):
        return {s for e in self.t for s, p in zip(SYMS, e) if p}


def _divide(a, b):
    """Exact multivariate division a/b; returns the quotient or None."""
    if b.is_zero():
        return None
    lead_b = max(b.t)
    cb = b.t[lead_b]
    rem, quo = dict(a.t), {}
    while rem:
        lead = max(rem)
        mono = tuple(x - y for x, y in zip(lead, lead_b))
        if min(mono) < 0:
            return None
        c = rem[lead] / cb
        quo[mono] = c
        for e, cc in b.t.items():
            ee = tuple(x + y for x, y in zip(e, mono))
            v = rem.get(ee, 0) - c * cc
            if v == 0:
                rem.pop(ee, None)
            else:
                rem[ee] = v
    return Poly(quo)


def _content(poly):
    """Split poly into (rational content, monomial exponents, primitive part)."""
    exps = tuple(min(e[i] for e in poly.t) for i in range(_N))
    num = reduce(gcd, (c.numerator for c in poly.t.values()))
    den = reduce(lambda x, y: x * y // gcd(x, y),
                 (c.denominator for c in poly.t.values()))
    cont = Fraction(num, den)
    prim = Poly({tuple(a - b for a, b in zip(e, exps)): c / cont
                 for e, c in poly.t.items()})
    return cont, exps, prim


class Ratio:
    """Quotient of two Poly objects, kept in a normal form.

    Normalisation cancels the common monomial factor and the rational content,
    and performs an exact division when one polynomial divides the other.  A
    general multivariate polynomial GCD is not implemented; the 336 first
    derivatives of the eight variables do not need one (checked in the tests).
    """

    def __init__(self, num, den=None):
        den = Poly.const(1) if den is None else den
        if den.is_zero():
            raise ZeroDivisionError("zero denominator")
        if num.is_zero():
            self.num, self.den = num, Poly.const(1)
            return
        allexp = list(num.t) + list(den.t)
        low = tuple(min(e[i] for e in allexp) for i in range(_N))
        num, den = num.shift(low), den.shift(low)
        cn, _, _ = _content(num)
        cd, _, _ = _content(den)
        num, den = num * (1 / cn), den * (1 / cd)
        scale = cn / cd
        q = _divide(num, den)
        if q is not None:
            num, den = q, Poly.const(1)
        else:
            q = _divide(den, num)
            if q is not None:
                num, den = Poly.const(1), q
        num = num * scale.numerator
        den = den * scale.denominator
        if den.t[max(den.t)] < 0:
            num, den = -num, -den
        self.num, self.den = num, den

    def __eq__(self, other):
        return (self.num * other.den - other.num * self.den).is_zero()

    def __mul__(self, other):
        if isinstance(other, Ratio):
            return Ratio(self.num * other.num, self.den * other.den)
        return Ratio(self.num * other, self.den)

    def __truediv__(self, other):
        return Ratio(self.num * other.den, self.den * other.num)

    def __sub__(self, other):
        return Ratio(self.num * other.den - other.num * self.den,
                     self.den * other.den)

    def __add__(self, other):
        return Ratio(self.num * other.den + other.num * self.den,
                     self.den * other.den)

    def symbols(self):
        """Symbols the reduced expression actually depends on."""
        return self.num.symbols() | self.den.symbols()

    def eval(self, values):
        missing = sorted(s for s in self.symbols() if s not in values)
        if missing:
            raise MissingDataError(
                "cannot evaluate: missing " + ", ".join(missing), missing)
        den = self.den.eval(values)
        # the denominator is a sum of products; it is zero when the terms cancel
        size = sum(abs(Poly({e: c}).eval(values)) for e, c in self.den.t.items())
        if abs(den) <= 1e-12 * size:
            raise SingularDerivativeError(
                f"the derivative is singular at this state: its denominator "
                f"{self._poly(self.den, False)[0]} vanishes")
        return self.num.eval(values) / den

    # ---- printing -------------------------------------------------------
    @staticmethod
    def _mono(exps, latex):
        parts = []
        for s, p in zip(SYMS, exps):
            if not p:
                continue
            name = _LATEX[s] if latex else s
            if p == 1:
                parts.append(name)
            elif latex:
                parts.append(f"{name}^{{{p}}}")
            else:
                parts.append(f"{name}^{p}")
        return (" " if latex else "*").join(parts)

    @classmethod
    def _poly(cls, poly, latex):
        """Return (text, is_sum) with the monomial content factored out."""
        if poly.is_zero():
            return "0", False
        cont, exps, prim = _content(poly)
        terms = sorted(prim.t.items(), key=lambda kv: kv[0], reverse=True)
        if len(terms) > 1 and terms[0][1] < 0 and any(c > 0 for _, c in terms):
            terms = ([t for t in terms if t[1] > 0]
                     + [t for t in terms if t[1] < 0])
        if len(terms) > 1 and all(c < 0 for _, c in terms):
            cont, terms = -cont, [(e, -c) for e, c in terms]
        pieces = []
        for k, (e, c) in enumerate(terms):
            m = cls._mono(e, latex)
            mag = abs(c)
            body = m if (mag == 1 and m) else (
                f"{mag}" + ((" " if latex else "*") + m if m else ""))
            sign = "-" if c < 0 else "+"
            pieces.append((f"{sign} " if k else ("-" if c < 0 else "")) + body)
        inner = " ".join(pieces)
        outer = cls._mono(exps, latex)
        if len(terms) == 1:
            sep = " " if latex else "*"
            factors = [f for f in (outer, inner) if f and f != "1"]
            text = sep.join(factors) if factors else "1"
            if cont != 1:
                text = ("-" if cont < 0 else "") + (
                    f"{abs(cont)}{sep}" if abs(cont) != 1 else "") + text
            return text, False
        prefix = ""
        if cont != 1:
            prefix = "-" if cont == -1 else f"{cont}" + ("" if latex else "*")
        if outer or prefix:
            if latex:
                return f"{prefix}{outer}\\left({inner}\\right)", False
            return f"{prefix}{outer}{'*' if outer else ''}({inner})", False
        return inner, True

    def __str__(self):
        n, n_sum = self._poly(self.num, False)
        if self.den == Poly.const(1):
            return n
        d, d_sum = self._poly(self.den, False)
        n = f"({n})" if n_sum else n
        d = f"({d})" if (d_sum or "*" in d) else d
        return f"{n}/{d}"

    def latex(self):
        n, _ = self._poly(self.num, True)
        if self.den == Poly.const(1):
            return n
        d, _ = self._poly(self.den, True)
        if n.startswith("-"):
            return f"-\\frac{{{n[1:]}}}{{{d}}}"
        return f"\\frac{{{n}}}{{{d}}}"


def _differentials():
    """(dX/dT)_P and (dX/dP)_T for the eight variables, from the fundamental
    relations dU = T dS - P dV, H = U + PV, A = U - TS, G = H - TS."""
    m, c = Poly.mono, Poly.const
    dT = (c(1), c(0))
    dP = (c(0), c(1))
    dV = (m(V=1, alpha=1), m(-1, V=1, kappaT=1))
    dS = (m(Cp=1, T=-1), m(-1, V=1, alpha=1))          # Maxwell: (dS/dP)_T = -(dV/dT)_P
    T, P, V, S = m(T=1), m(P=1), m(V=1), m(S=1)
    dU = tuple(T * s - P * v for s, v in zip(dS, dV))
    dH = tuple(T * s + V * p for s, p in zip(dS, dP))
    dA = tuple(-(S * t) - P * v for t, v in zip(dT, dV))
    dG = tuple(-(S * t) + V * p for t, p in zip(dT, dP))
    return {"T": dT, "P": dP, "V": dV, "S": dS,
            "U": dU, "H": dH, "A": dA, "G": dG}


DIFF = _differentials()


def jacobian(x, z):
    """J(x, z) = d(x, z)/d(T, P)."""
    return DIFF[x][0] * DIFF[z][1] - DIFF[x][1] * DIFF[z][0]


def derivative(x, y, z):
    """Reduced form of (dx/dy)_z as a :class:`Ratio`."""
    for v in (x, y, z):
        if v not in DIFF:
            raise InputError(f"'{v}' is not one of {', '.join(VARIABLES)}")
    if y == z:
        raise InputError(f"(d{x}/d{y})_{z} is undefined: {y} cannot vary at constant {z}")
    if x == z:
        return Ratio(Poly.const(0))
    return Ratio(jacobian(x, z), jacobian(y, z))


_DERIV = re.compile(r"^\s*\(?\s*d\s*([A-Za-z]+)\s*/\s*d\s*([A-Za-z]+)\s*\)?\s*[_|]\s*([A-Za-z]+)\s*$")


def parse_derivative(text):
    """Read '(dT/dP)_H' or 'dT/dP|H' and return ('T', 'P', 'H')."""
    m = _DERIV.match(text)
    if not m:
        raise InputError(f"cannot read derivative '{text}'; write it as (dX/dY)_Z")
    return m.group(1), m.group(2), m.group(3)


def reduce_derivative(text):
    return derivative(*parse_derivative(text))


def all_derivatives():
    """Every (x, y, z) with three distinct variables: 8*7*6 = 336 triples."""
    return [(x, y, z) for x in VARIABLES for y in VARIABLES for z in VARIABLES
            if len({x, y, z}) == 3]


# The four Maxwell relations, each as (lhs triple, sign, rhs triple):
# lhs = sign * rhs.
MAXWELL = [
    (("T", "V", "S"), -1, ("P", "S", "V"), "dU = T dS - P dV"),
    (("T", "P", "S"), +1, ("V", "S", "P"), "dH = T dS + V dP"),
    (("S", "V", "T"), +1, ("P", "T", "V"), "dA = -S dT - P dV"),
    (("S", "P", "T"), -1, ("V", "T", "P"), "dG = -S dT + V dP"),
]


def check_maxwell():
    """Return a list of (text, holds) for the four Maxwell relations."""
    out = []
    for lhs, sign, rhs, origin in MAXWELL:
        a, b = derivative(*lhs), derivative(*rhs) * sign
        txt = (f"(d{lhs[0]}/d{lhs[1]})_{lhs[2]} = "
               f"{'-' if sign < 0 else ''}(d{rhs[0]}/d{rhs[1]})_{rhs[2]}")
        out.append((txt, origin, str(a), a == b))
    return out
