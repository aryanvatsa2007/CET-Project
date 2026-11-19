"""Unit handling without third-party dependencies.

Every input quantity is written as ``<number> <unit>`` (for example ``350 K``,
``10 MPa``, ``25 degC``, ``0.19 L/mol``, ``37.1 J/(mol*K)``).  The parser converts
it to SI and returns its dimension vector, so that a pressure given in kelvin,
or a number with no unit at all, is rejected instead of silently used.

Compound units are read strictly left to right: ``J/mol/K`` and ``J/(mol*K)``
are both joule per mole per kelvin, whereas ``J/mol K`` means (J/mol)*K.
"""
import re

from .errors import UnitError

BASE = ("kg", "m", "s", "K", "mol")


def _d(kg=0, m=0, s=0, K=0, mol=0):
    return (kg, m, s, K, mol)


_PRESSURE = _d(kg=1, m=-1, s=-2)
_ENERGY = _d(kg=1, m=2, s=-2)
_VOLUME = _d(m=3)

# name -> (factor to SI, dimension vector)
UNITS = {
    "1": (1.0, _d()),
    "%": (0.01, _d()),
    "K": (1.0, _d(K=1)),
    "kg": (1.0, _d(kg=1)),
    "g": (1e-3, _d(kg=1)),
    "m": (1.0, _d(m=1)),
    "dm": (1e-1, _d(m=1)),
    "cm": (1e-2, _d(m=1)),
    "mm": (1e-3, _d(m=1)),
    "s": (1.0, _d(s=1)),
    "mol": (1.0, _d(mol=1)),
    "kmol": (1e3, _d(mol=1)),
    "Pa": (1.0, _PRESSURE),
    "kPa": (1e3, _PRESSURE),
    "MPa": (1e6, _PRESSURE),
    "GPa": (1e9, _PRESSURE),
    "bar": (1e5, _PRESSURE),
    "mbar": (1e2, _PRESSURE),
    "atm": (101325.0, _PRESSURE),
    "psi": (6894.757293168361, _PRESSURE),
    "torr": (101325.0 / 760.0, _PRESSURE),
    "J": (1.0, _ENERGY),
    "kJ": (1e3, _ENERGY),
    "MJ": (1e6, _ENERGY),
    "cal": (4.184, _ENERGY),
    "kcal": (4184.0, _ENERGY),
    "L": (1e-3, _VOLUME),
    "l": (1e-3, _VOLUME),
    "mL": (1e-6, _VOLUME),
    "m3": (1.0, _VOLUME),
    "dm3": (1e-3, _VOLUME),
    "cm3": (1e-6, _VOLUME),
    "N": (1.0, _d(kg=1, m=1, s=-2)),
}

# affine temperature scales: T[K] = scale * (value + shift)
AFFINE = {
    "degC": (1.0, 273.15),
    "°C": (1.0, 273.15),
    "degF": (5.0 / 9.0, 459.67),
    "°F": (5.0 / 9.0, 459.67),
    "degR": (5.0 / 9.0, 0.0),
}

# dimension expected for each engine variable
VAR_DIMS = {
    "T": _d(K=1),
    "P": _PRESSURE,
    "V": _d(m=3, mol=-1),
    "S": _d(kg=1, m=2, s=-2, K=-1, mol=-1),
    "U": _d(kg=1, m=2, s=-2, mol=-1),
    "H": _d(kg=1, m=2, s=-2, mol=-1),
    "A": _d(kg=1, m=2, s=-2, mol=-1),
    "G": _d(kg=1, m=2, s=-2, mol=-1),
    "Cp": _d(kg=1, m=2, s=-2, K=-1, mol=-1),
    "Cv": _d(kg=1, m=2, s=-2, K=-1, mol=-1),
    "alpha": _d(K=-1),
    "kappaT": _d(kg=-1, m=1, s=2),
    "kappaS": _d(kg=-1, m=1, s=2),
    "gammaV": _d(kg=1, m=-1, s=-2, K=-1),
    "muJT": _d(kg=-1, m=1, s=2, K=1),
    "w": _d(m=1, s=-1),
    "M": _d(kg=1, mol=-1),
    "rho": _d(kg=1, m=-3),
    "Z": _d(),
    "gamma": _d(),
    "omega": _d(),
}

# SI unit string used when writing each variable
SI_UNIT = {
    "T": "K", "P": "Pa", "V": "m3/mol", "S": "J/(mol*K)", "U": "J/mol",
    "H": "J/mol", "A": "J/mol", "G": "J/mol", "Cp": "J/(mol*K)",
    "Cv": "J/(mol*K)", "alpha": "1/K", "kappaT": "1/Pa", "kappaS": "1/Pa",
    "gammaV": "Pa/K", "muJT": "K/Pa", "w": "m/s", "M": "kg/mol",
    "rho": "kg/m3", "Z": "1", "gamma": "1", "omega": "1",
}

_TOKEN = re.compile(r"\s*(?:([A-Za-z°%]+\d?|1)|\^\s*(-?\d+)|([*/()·]))")
_NUMBER = re.compile(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*(.*?)\s*$")


def _tokens(text):
    pos, out = 0, []
    text = text.strip()
    while pos < len(text):
        m = _TOKEN.match(text, pos)
        if not m:
            raise UnitError(f"cannot read unit '{text}' near '{text[pos:]}'")
        if m.group(1) is not None:
            out.append(("u", m.group(1)))
        elif m.group(2) is not None:
            out.append(("^", int(m.group(2))))
        else:
            out.append(("op", m.group(3)))
        pos = m.end()
    return out


def parse_unit(text):
    """Return ``(factor_to_SI, dimension_vector)`` for a multiplicative unit."""
    toks = _tokens(text)
    if not toks:
        return 1.0, _d()
    i = 0

    def peek():
        return toks[i] if i < len(toks) else None

    def term():
        nonlocal i
        t = peek()
        if t is None:
            raise UnitError(f"incomplete unit '{text}'")
        if t == ("op", "("):
            i += 1
            val = expr()
            if peek() != ("op", ")"):
                raise UnitError(f"unbalanced parentheses in unit '{text}'")
            i += 1
        elif t[0] == "u":
            i += 1
            if t[1] in AFFINE:
                raise UnitError(
                    f"'{t[1]}' is an offset temperature scale and cannot be "
                    f"combined with other units ('{text}')")
            if t[1] not in UNITS:
                raise UnitError(f"unknown unit '{t[1]}' in '{text}'")
            val = UNITS[t[1]]
        else:
            raise UnitError(f"cannot read unit '{text}'")
        t = peek()
        if t is not None and t[0] == "^":
            i += 1
            val = (val[0] ** t[1], tuple(e * t[1] for e in val[1]))
        return val

    def expr():
        nonlocal i
        val = term()
        while peek() is not None and peek() != ("op", ")"):
            t = peek()
            divide = t == ("op", "/")
            if t[0] == "op" and t[1] in "*/·":
                i += 1
            nxt = term()
            if divide:
                val = (val[0] / nxt[0], tuple(a - b for a, b in zip(val[1], nxt[1])))
            else:
                val = (val[0] * nxt[0], tuple(a + b for a, b in zip(val[1], nxt[1])))
        return val

    val = expr()
    if i != len(toks):
        raise UnitError(f"cannot read unit '{text}'")
    return val


def dims_str(dims):
    parts = [f"{b}^{e}" if e != 1 else b for b, e in zip(BASE, dims) if e]
    return " ".join(parts) if parts else "dimensionless"


def parse_quantity(text, var=None):
    """Convert ``'<number> <unit>'`` to an SI float.

    If ``var`` names an engine variable, the dimension is checked against it.
    """
    m = _NUMBER.match(str(text))
    if not m:
        raise UnitError(f"cannot read a number from '{text}'"
                        + (f" for {var}" if var else ""))
    value, unit = float(m.group(1)), m.group(2)
    if unit in AFFINE:
        scale, shift = AFFINE[unit]
        si, dims = scale * (value + shift), _d(K=1)
    elif unit == "":
        if var is not None and VAR_DIMS.get(var) != _d():
            raise UnitError(
                f"{var} = '{text}' has no unit; expected a unit of dimension "
                f"{dims_str(VAR_DIMS[var])} (for example {SI_UNIT[var]})")
        si, dims = value, _d()
    else:
        factor, dims = parse_unit(unit)
        si = value * factor
    if var is not None:
        want = VAR_DIMS.get(var)
        if want is None:
            raise UnitError(f"'{var}' is not a quantity known to the engine")
        if tuple(dims) != tuple(want):
            raise UnitError(
                f"{var} = '{text}' has dimension {dims_str(dims)}; expected "
                f"{dims_str(want)} (for example {SI_UNIT[var]})")
    return si


def convert(value_si, unit):
    """Convert an SI value to ``unit`` (multiplicative or offset)."""
    if unit in AFFINE:
        scale, shift = AFFINE[unit]
        return value_si / scale - shift
    factor, _ = parse_unit(unit)
    return value_si / factor
