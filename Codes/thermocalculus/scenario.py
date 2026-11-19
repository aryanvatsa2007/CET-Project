"""Scenario files: read an input file, run it, write labelled output files.

A scenario is a plain-text INI file (see Codes/README.md for the format).  Each
run writes two files:

    <name>_output.csv   quantity, value, unit, provenance, std_uncertainty, note
    <name>_log.txt      what was done, warnings, and the reason for any failure

A failing scenario does not crash: the failure type and message are written to
both files so that edge cases are part of the recorded results.
"""
import configparser
import csv
import math
import os

from . import knowledge, props
from .constants import P_REF, R, T_REF
from .eos import EoS, Fluid, IdealGasCp, get_fluid
from .errors import (InconsistentInputError, InputError, MissingDataError,
                     SingularDerivativeError, ThermoError)
from .symbolic import SYMS, derivative, parse_derivative
from .tables import PVTTable, read_table_csv
from .uncertainty import monte_carlo, propagate
from .units import SI_UNIT, convert, parse_quantity

TYPES = ("state", "process", "saturation", "derivative", "infer", "table")
_DIM_OF = {"Tc": "T", "Pc": "P", "omega": "omega", "M": "M"}


def _var(key):
    """Engine variable whose dimension a key such as 'T1' or 'Pc' must have."""
    if key in _DIM_OF:
        return _DIM_OF[key]
    return key.rstrip("0123456789")


class Result:
    def __init__(self, name, kind):
        self.name, self.kind = name, kind
        self.rows, self.log = [], []
        self.status = "ok"
        self.budget = {}          # quantity -> (value, u, {input: variance share})

    def add(self, quantity, value, unit, provenance, u=None, note=""):
        self.rows.append({"quantity": quantity, "value": value, "unit": unit,
                          "provenance": provenance, "std_uncertainty": u,
                          "mc_std_uncertainty": None, "note": note})

    def say(self, text):
        self.log.append(text)

    def get(self, quantity):
        for r in self.rows:
            if r["quantity"] == quantity:
                return r["value"]
        raise KeyError(quantity)


# ---------------------------------------------------------------------------
def _read(path):
    cp = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    cp.optionxform = str
    if not cp.read(path):
        raise InputError(f"cannot read scenario file '{path}'")
    if "scenario" not in cp:
        raise InputError("the scenario file needs a [scenario] section")
    return cp


def _build_model(cp, base, res):
    """Return (fluid, eos) from the [model] section."""
    sec = cp["model"] if "model" in cp else {}
    name = sec.get("fluid")
    if name:
        fluid = get_fluid(name)
        res.say(f"fluid {fluid.name}: critical constants from {fluid.source_critical}")
    else:
        fluid = Fluid("user fluid", "", None, None, None, None)
    for key, attr in (("Tc", "Tc"), ("Pc", "Pc"), ("omega", "omega"), ("M", "M")):
        if key in sec:
            setattr(fluid, attr, parse_quantity(sec[key], _var(key)))
    extrap = sec.get("extrapolate_cp0", "no").strip().lower() in ("yes", "true", "1")
    if "cp0_shomate" in sec:
        vals = [float(x) for x in sec["cp0_shomate"].replace(",", " ").split()]
        if len(vals) != 7:
            raise InputError("cp0_shomate needs seven numbers: A B C D E Tmin Tmax")
        fluid.cp0 = IdealGasCp("shomate", vals[:5], vals[5], vals[6])
    elif "cp0_constant" in sec:
        fluid.cp0 = IdealGasCp("constant", [parse_quantity(sec["cp0_constant"], "Cp")])
    elif "cp0_table" in sec:
        cols = read_table_csv(os.path.join(base, sec["cp0_table"]))
        if "Cp0" not in cols or "T" not in cols:
            raise InputError("a Cp0 table needs columns 'T [..]' and 'Cp0 [..]'")
        fluid.cp0 = IdealGasCp("table", table=(cols["T"], cols["Cp0"]))
        res.say(f"ideal-gas Cp from table {sec['cp0_table']} "
                f"({fluid.cp0.tmin:g} K to {fluid.cp0.tmax:g} K)")
    if fluid.cp0 is not None:
        fluid.cp0.extrapolate = extrap
        if fluid.cp0.tref != T_REF:
            res.say(f"NOTE: the Cp0 model does not cover {T_REF} K; H and S are "
                    f"referred to the ideal gas at {fluid.cp0.tref:g} K and {P_REF:g} Pa")
    kind = sec.get("eos", "PR")
    eos = EoS(kind, fluid.Tc, fluid.Pc, fluid.omega)
    res.say("equation of state " + eos.describe())
    return fluid, eos


def _state_inputs(sec, suffix=""):
    lack = [k for k in ("T", "P") if k not in sec]
    if lack:
        raise MissingDataError(
            f"section [{sec.name}] needs " + " and ".join(lack), lack)
    return (parse_quantity(sec["T"], "T"), parse_quantity(sec["P"], "P"),
            sec.get("phase", "stable").strip())


def _unc(cp):
    """Standard uncertainties {key: SI value} and Monte Carlo sample count."""
    if "uncertainty" not in cp:
        return {}, 0
    u, n = {}, 0
    for key, text in cp["uncertainty"].items():
        if key == "monte_carlo":
            n = int(text)
            continue
        var = _var(key)
        if var == "T":                       # an uncertainty is a difference
            val = parse_quantity(text.replace("degC", "K"), "T")
        else:
            val = parse_quantity(text, var)
        u[key] = abs(val)
    return u, n


def _deriv_unit(x, y):
    ux, uy = SI_UNIT[x], SI_UNIT[y]
    return f"({ux})/({uy})"


def _note(contrib):
    top = sorted(contrib.items(), key=lambda kv: -kv[1])[:3]
    return "variance share: " + ", ".join(f"{k} {100 * v:.1f}%" for k, v in top)


def _apply_uncertainty(res, func, x, u, n_mc, provenance):
    """Attach propagated uncertainties to rows already present in ``res``."""
    lack = [k for k in u if k not in x]
    if lack:
        raise InputError("uncertainty given for unknown input(s): " + ", ".join(lack))
    lin = propagate(func, x, u)
    res.budget = lin
    for row in res.rows:
        q = row["quantity"]
        if q in lin and lin[q][1] > 0:
            row["std_uncertainty"] = lin[q][1]
            row["note"] = _note(lin[q][2])
    res.say("uncertainty: first-order propagation of "
            + ", ".join(f"u({k}) = {v:g}" for k, v in u.items()) + " (SI units)")
    if n_mc:
        mc, failed = monte_carlo(func, x, u, n_mc)
        res.say(f"Monte Carlo check with {n_mc} draws ({failed} failed); "
                "ratio = Monte Carlo u / linear u:")
        for row in res.rows:
            q = row["quantity"]
            if q in mc and q in lin and lin[q][1] > 0:
                row["mc_std_uncertainty"] = mc[q][1]
                res.say(f"  {q}: linear u = {lin[q][1]:.6g}, Monte Carlo u = "
                        f"{mc[q][1]:.6g}, ratio = {mc[q][1] / lin[q][1]:.4f}")
    return lin


def _targets(cp):
    sec = cp["targets"] if "targets" in cp else {}
    ders = [d.strip() for d in sec.get("derivatives", "").split(",") if d.strip()]
    qts = [q.strip() for q in sec.get("quantities", "").split(",") if q.strip()]
    return ders, qts


def _add_state_rows(res, st, prov, prefix=""):
    for key, unit in props.STATE_UNITS:
        if key in st:
            p = "stated" if key in ("T", "P") else prov
            res.add(prefix + key, st[key], unit, p)


def _add_derivatives(res, ders, values, prov):
    for text in ders:
        x, y, z = parse_derivative(text)
        expr = derivative(x, y, z)
        try:
            val = expr.eval(values)
        except MissingDataError as exc:
            res.status = "partial"
            res.say(f"MISSING DATA: {text} = {expr} : {exc}")
            continue
        except SingularDerivativeError as exc:
            res.status = "partial"
            res.say(f"SINGULAR: {text} = {expr} : {exc}")
            continue
        label = f"(d{x}/d{y})_{z}"
        res.add(label, val, _deriv_unit(x, y), prov, note=f"= {expr}")


# ---------------------------------------------------------------------------
def _run_state(cp, base, res):
    fluid, eos = _build_model(cp, base, res)
    if "state" not in cp:
        raise InputError("a state scenario needs a [state] section")
    T, P, phase = _state_inputs(cp["state"])
    prov = f"parameterised ({eos.kind})"
    st = props.state(fluid, eos, T, P, phase)
    res.say(f"phase: {st['phase']}; volume roots at this state: "
            + ", ".join(f"{v:.6g}" for v in eos.volumes(T, P)) + " m3/mol")
    _add_state_rows(res, st, prov)
    if st["missing"]:
        res.status = "partial"
        res.say("MISSING DATA: " + ", ".join(st["missing"]) + " not available; "
                "volumetric and residual properties only")
    ders, _ = _targets(cp)
    _add_derivatives(res, ders, {s: st[s] for s in SYMS if s in st}, prov)
    if fluid.cp0 is not None:
        for w in fluid.cp0.warnings:
            res.say("WARNING: " + w)
    u, n_mc = _unc(cp)
    if u:
        x = {"T": T, "P": P, "Tc": fluid.Tc, "Pc": fluid.Pc, "omega": fluid.omega}
        x = {k: v for k, v in x.items() if v is not None}

        def func(xi):
            e = EoS(eos.kind, xi.get("Tc"), xi.get("Pc"), xi.get("omega"))
            s = props.state(fluid, e, xi["T"], xi["P"], phase)
            out = {k: v for k, v in s.items() if isinstance(v, float)}
            vals = {q: s[q] for q in SYMS if q in s}
            for text in ders:
                a, b, c = parse_derivative(text)
                try:
                    out[f"(d{a}/d{b})_{c}"] = derivative(a, b, c).eval(vals)
                except (MissingDataError, SingularDerivativeError):
                    pass
            return out
        _apply_uncertainty(res, func, x, u, n_mc, prov)


def _run_process(cp, base, res):
    fluid, eos = _build_model(cp, base, res)
    for s in ("state1", "state2"):
        if s not in cp:
            raise InputError("a process scenario needs [state1] and [state2] sections")
    T1, P1, ph1 = _state_inputs(cp["state1"])
    T2, P2, ph2 = _state_inputs(cp["state2"])
    prov = f"parameterised ({eos.kind})"
    out, s1, s2 = props.process(fluid, eos, T1, P1, T2, P2, ph1, ph2)
    res.say(f"state 1: {s1['phase']}; state 2: {s2['phase']}")
    for key, unit in (("dH", "J/mol"), ("dS", "J/(mol*K)"), ("dU", "J/mol"),
                      ("dA", "J/mol"), ("dG", "J/mol"), ("dV", "m3/mol"),
                      ("dH_ideal", "J/mol"), ("dH_residual", "J/mol"),
                      ("dS_ideal", "J/(mol*K)"), ("dS_residual", "J/(mol*K)")):
        res.add(key, out[key], unit, prov)
    _add_state_rows(res, s1, prov, "state1.")
    _add_state_rows(res, s2, prov, "state2.")
    for w in fluid.cp0.warnings:
        res.say("WARNING: " + w)
    u, n_mc = _unc(cp)
    if u:
        x = {"T1": T1, "P1": P1, "T2": T2, "P2": P2,
             "Tc": fluid.Tc, "Pc": fluid.Pc, "omega": fluid.omega}
        x = {k: v for k, v in x.items() if v is not None}

        def func(xi):
            e = EoS(eos.kind, xi.get("Tc"), xi.get("Pc"), xi.get("omega"))
            o, _, _ = props.process(fluid, e, xi["T1"], xi["P1"], xi["T2"], xi["P2"],
                                    ph1, ph2)
            return o
        _apply_uncertainty(res, func, x, u, n_mc, prov)


def _run_saturation(cp, base, res):
    fluid, eos = _build_model(cp, base, res)
    if "state" not in cp or "T" not in cp["state"]:
        raise MissingDataError("a saturation scenario needs T in [state]", ["T"])
    T = parse_quantity(cp["state"]["T"], "T")
    prov = f"parameterised ({eos.kind})"

    def calc(e, T):
        P, VL, VV = e.psat(T)
        rl, rv = e.residual(T, VL), e.residual(T, VV)
        dH = (rv["UR"] + P * VV) - (rl["UR"] + P * VL)
        return {"Psat": P, "VL": VL, "VV": VV, "dHvap": dH, "dSvap": dH / T}

    out = calc(eos, T)
    res.add("T", T, "K", "stated")
    for key, unit in (("Psat", "Pa"), ("VL", "m3/mol"), ("VV", "m3/mol"),
                      ("dHvap", "J/mol"), ("dSvap", "J/(mol*K)")):
        res.add(key, out[key], unit, prov)
    res.say("saturation from equality of fugacity in both phases; the enthalpy "
            "of vaporisation is a residual-property difference and needs no Cp0")
    u, n_mc = _unc(cp)
    if u:
        x = {"T": T, "Tc": fluid.Tc, "Pc": fluid.Pc, "omega": fluid.omega}
        x = {k: v for k, v in x.items() if v is not None}
        func = lambda xi: calc(EoS(eos.kind, xi.get("Tc"), xi.get("Pc"),
                                   xi.get("omega")), xi["T"])
        _apply_uncertainty(res, func, x, u, n_mc, prov)


def _run_derivative(cp, base, res):
    """Symbolic reduction only; numeric values if [given] supplies the symbols."""
    ders, _ = _targets(cp)
    if not ders:
        raise InputError("a derivative scenario needs 'derivatives = ...' in [targets]")
    values = {}
    if "given" in cp:
        for key, text in cp["given"].items():
            values[key] = parse_quantity(text, _var(key))
    for text in ders:
        x, y, z = parse_derivative(text)
        expr = derivative(x, y, z)
        label = f"(d{x}/d{y})_{z}"
        needs = sorted(expr.symbols())
        res.say(f"{label} = {expr}    [needs: {', '.join(needs) or 'nothing'}]")
        lack = [s for s in needs if s not in values]
        if lack:
            if values:
                res.status = "partial"
                res.say(f"MISSING DATA: {label} not evaluated; missing "
                        + ", ".join(lack))
            res.add(label, float("nan"), _deriv_unit(x, y), "symbolic only",
                    note=f"= {expr}")
        else:
            res.add(label, expr.eval(values), _deriv_unit(x, y),
                    "inferred from stated values", note=f"= {expr}")


def _run_infer(cp, base, res):
    if "given" not in cp:
        raise InputError("an infer scenario needs a [given] section")
    known = {k: parse_quantity(t, _var(k)) for k, t in cp["given"].items()}
    rtol = float(cp["scenario"].get("consistency_rtol", "0.01"))
    rep = knowledge.infer(known, rtol=rtol)
    for name, fact in rep.facts.items():
        res.add(name, fact.value, SI_UNIT.get(name, ""), fact.provenance,
                note="" if fact.how == "input" else fact.how)
    for rel, q, old, new, diff, ok in rep.checks:
        res.say(f"{'check passed' if ok else 'INCONSISTENT'}: {rel} gives {q} = "
                f"{new:.6g}, but {q} = {old:.6g} is already known "
                f"(relative difference {diff:.3g}, tolerance {rtol:g})")
    for v in rep.violations:
        res.say("INCONSISTENT: " + v)
    if not rep.consistent:
        res.status = "inconsistent"
    ders, qts = _targets(cp)
    for q in qts:
        try:
            knowledge.require(rep, q)
        except MissingDataError as exc:
            if res.status == "ok":
                res.status = "partial"
            res.say("MISSING DATA: " + str(exc))
    for text in ders:
        x, y, z = parse_derivative(text)
        try:
            expr, val = knowledge.derivative_from_facts(text, rep)
            res.add(f"(d{x}/d{y})_{z}", val, _deriv_unit(x, y), "inferred",
                    note=f"= {expr}")
        except MissingDataError as exc:
            if res.status == "ok":
                res.status = "partial"
            res.say("MISSING DATA: " + str(exc))
    u, n_mc = _unc(cp)
    if u:
        def func(xi):
            r = knowledge.infer(xi, rtol=rtol)
            out = {k: f.value for k, f in r.facts.items()}
            for text in ders:
                a, b, c = parse_derivative(text)
                try:
                    out[f"(d{a}/d{b})_{c}"] = knowledge.derivative_from_facts(text, r)[1]
                except MissingDataError:
                    pass
            return out
        _apply_uncertainty(res, func, known, u, n_mc, "inferred")


def _run_table(cp, base, res):
    if "table" not in cp or "pvt" not in cp["table"]:
        raise InputError("a table scenario needs 'pvt = <file.csv>' in [table]")
    path = os.path.join(base, cp["table"]["pvt"])
    table = PVTTable.from_csv(path)
    res.say(f"P-V-T table {cp['table']['pvt']}: {len(table.Tg)} temperatures "
            f"({table.Tg[0]:g} to {table.Tg[-1]:g} K) x {len(table.Pg)} pressures "
            f"({table.Pg[0]:g} to {table.Pg[-1]:g} Pa); bicubic spline in Z(T,P)")
    fluid = Fluid("tabulated fluid", "", None, None, None, None)
    if "model" in cp:
        sec = cp["model"]
        if sec.get("fluid"):
            fluid = get_fluid(sec["fluid"])
        if "M" in sec:
            fluid.M = parse_quantity(sec["M"], "M")
        if "cp0_table" in sec:
            cols = read_table_csv(os.path.join(base, sec["cp0_table"]))
            fluid.cp0 = IdealGasCp("table", table=(cols["T"], cols["Cp0"]))
        elif "cp0_constant" in sec:
            fluid.cp0 = IdealGasCp("constant", [parse_quantity(sec["cp0_constant"], "Cp")])
    prov = "tabulated (spline)"
    if "state1" in cp and "state2" in cp:
        T1, P1, _ = _state_inputs(cp["state1"])
        T2, P2, _ = _state_inputs(cp["state2"])
        s1 = table.state(T1, P1, fluid.cp0, fluid.M)
        s2 = table.state(T2, P2, fluid.cp0, fluid.M)
        if "H" in s1:
            for key, unit in (("H", "J/mol"), ("S", "J/(mol*K)"), ("U", "J/mol")):
                res.add("d" + key, s2[key] - s1[key], unit, prov)
        else:
            res.status = "partial"
            res.say("MISSING DATA: Cp0(T) not supplied; only residual changes reported")
        res.add("dH_residual", s2["HR"] - s1["HR"], "J/mol", prov)
        res.add("dS_residual", s2["SR"] - s1["SR"], "J/(mol*K)", prov)
        _add_state_rows(res, s1, prov, "state1.")
        _add_state_rows(res, s2, prov, "state2.")
        states = (s1, s2)
    elif "state" in cp:
        T, P, _ = _state_inputs(cp["state"])
        st = table.state(T, P, fluid.cp0, fluid.M)
        _add_state_rows(res, st, prov)
        if st["missing"]:
            res.status = "partial"
            res.say("MISSING DATA: " + ", ".join(st["missing"]) + " not supplied; "
                    "volumetric and residual properties only")
        ders, _ = _targets(cp)
        _add_derivatives(res, ders, {s: st[s] for s in SYMS if s in st}, prov)
        states = (st,)
    else:
        raise InputError("a table scenario needs [state] or [state1] and [state2]")
    for k, st in enumerate(states, 1):
        for wmsg in st.get("warnings", []):
            res.say("WARNING: " + wmsg)
        res.say(f"state {k}: low-pressure (below table) correction to HR = "
                f"{st['HR_lowP_correction']:.4g} J/mol of {st['HR']:.4g} J/mol")


_RUNNERS = {"state": _run_state, "process": _run_process,
            "saturation": _run_saturation, "derivative": _run_derivative,
            "infer": _run_infer, "table": _run_table}


def run(path):
    """Run a scenario file and return a :class:`Result` (never raises ThermoError)."""
    base = os.path.dirname(os.path.abspath(path))
    name = os.path.splitext(os.path.basename(path))[0]
    res = Result(name, "?")
    try:
        cp = _read(path)
        res.name = cp["scenario"].get("name", name)
        res.kind = cp["scenario"].get("type", "").strip()
        if cp["scenario"].get("description"):
            res.say("description: " + cp["scenario"]["description"])
        if res.kind not in _RUNNERS:
            raise InputError(f"scenario type '{res.kind}' is not one of "
                             + ", ".join(TYPES))
        _RUNNERS[res.kind](cp, base, res)
        units = dict(cp["output_units"]) if "output_units" in cp else {}
        for row in res.rows:
            key = row["quantity"].split(".")[-1]
            if key in units and isinstance(row["value"], float):
                row["value"] = convert(row["value"], units[key])
                scale = convert(1.0, units[key]) - convert(0.0, units[key])
                for col in ("std_uncertainty", "mc_std_uncertainty"):
                    if row[col] is not None:
                        row[col] *= scale
                row["unit"] = units[key]
    except ThermoError as exc:
        res.status = "error: " + type(exc).__name__
        res.say(f"ERROR ({type(exc).__name__}): {exc}")
    except (configparser.Error, ValueError) as exc:
        res.status = "error: InputError"
        res.say(f"ERROR (InputError): {exc}")
    return res


def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return "nan" if math.isnan(v) else f"{v:.10g}"
    return str(v)


def write(res, outdir, source=""):
    """Write <name>_output.csv and <name>_log.txt; return their paths."""
    os.makedirs(outdir, exist_ok=True)
    out_csv = os.path.join(outdir, f"{res.name}_output.csv")
    out_log = os.path.join(outdir, f"{res.name}_log.txt")
    head = [f"# scenario: {res.name}", f"# type: {res.kind}",
            f"# status: {res.status}"]
    if source:
        head.append(f"# input file: {source}")
    head.append("# values are in SI units unless the unit column says otherwise; "
                f"reference state: ideal gas, H = S = 0 at {T_REF} K and {P_REF:g} Pa")
    with open(out_csv, "w", newline="") as fh:
        fh.write("\n".join(head) + "\n")
        wr = csv.writer(fh)
        wr.writerow(["quantity", "value", "unit", "provenance",
                     "std_uncertainty", "mc_std_uncertainty", "note"])
        for r in res.rows:
            wr.writerow([r["quantity"], _fmt(r["value"]), r["unit"],
                         r["provenance"], _fmt(r["std_uncertainty"]),
                         _fmt(r["mc_std_uncertainty"]), r["note"]])
    with open(out_log, "w") as fh:
        fh.write("\n".join(head) + "\n\n" + "\n".join(res.log) + "\n")
    if res.budget:
        inputs = sorted({k for _, _, sh in res.budget.values() for k in sh})
        with open(os.path.join(outdir, f"{res.name}_uncertainty_budget.csv"),
                  "w", newline="") as fh:
            fh.write("\n".join(head[:3]) + "\n# fraction of the output variance "
                     "contributed by each uncertain input (first-order propagation); "
                     "SI units\n")
            wr = csv.writer(fh)
            wr.writerow(["quantity", "value", "std_uncertainty"]
                        + [f"share_{k}" for k in inputs])
            for q, (val, uq, sh) in res.budget.items():
                if uq > 0:
                    wr.writerow([q, _fmt(val), _fmt(uq)]
                                + [f"{sh.get(k, 0.0):.6f}" for k in inputs])
    return out_csv, out_log
