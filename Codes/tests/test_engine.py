"""Test suite of the thermocalculus engine.  Run from the Codes folder:

    python -m pytest -q
"""
import csv
import glob
import math
import os
import sys

import numpy as np
import pytest
from scipy.optimize import fsolve

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "Data"))

from thermocalculus import knowledge, props, scenario, symbolic, units  # noqa: E402
from thermocalculus.constants import R  # noqa: E402
from thermocalculus.eos import EOS_NAMES, EoS, IdealGasCp, get_fluid  # noqa: E402
from thermocalculus.errors import (InconsistentInputError, MissingDataError,  # noqa: E402
                                   OutOfRangeError, UnitError)
from thermocalculus.numerics import NumericalModel  # noqa: E402
from thermocalculus.tables import PVTTable  # noqa: E402
from thermocalculus.uncertainty import monte_carlo, propagate  # noqa: E402

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")
CO2 = get_fluid("CO2")
N2 = get_fluid("N2")


def eos(kind, f=CO2):
    return EoS(kind, f.Tc, f.Pc, f.omega)


# ---------------------------------------------------------------- symbolic
def test_maxwell_relations_hold_symbolically():
    assert all(ok for *_, ok in symbolic.check_maxwell())


@pytest.mark.parametrize("text, expected", [
    ("(dT/dP)_H", "V*(T*alpha - 1)/Cp"),
    ("(dT/dP)_S", "T*V*alpha/Cp"),
    ("(dP/dT)_V", "alpha/kappaT"),
    ("(dS/dP)_T", "-V*alpha"),
    ("(dH/dP)_T", "V*(1 - T*alpha)"),
    ("(dU/dV)_T", "(T*alpha - P*kappaT)/kappaT"),
    ("(dH/dT)_P", "Cp"),
    ("(dG/dP)_T", "V"),
    ("(dG/dT)_P", "-S"),
    ("(dA/dV)_T", "-P"),
    ("(dU/dS)_V", "T"),
    ("(dH/dS)_P", "T"),
])
def test_textbook_reductions(text, expected):
    assert str(symbolic.reduce_derivative(text)) == expected


def test_heat_capacity_identities():
    m = symbolic.Poly.mono
    T = symbolic.Ratio(m(T=1))
    cv = symbolic.derivative("S", "T", "V") * T          # Cv = T (dS/dT)_V
    cp = symbolic.derivative("S", "T", "P") * T
    diff = symbolic.Ratio(m(T=1, V=1, alpha=2), m(kappaT=1))
    assert cp == symbolic.Ratio(m(Cp=1))
    assert (cp - cv) == diff                              # Cp - Cv = T V alpha^2 / kappaT
    assert symbolic.derivative("U", "T", "V") == cv       # Cv = (dU/dT)_V


def test_reciprocal_and_cyclic_rules():
    d = symbolic.derivative
    one = symbolic.Ratio(symbolic.Poly.const(1))
    for x, y, z in symbolic.all_derivatives():
        assert d(x, y, z) * d(y, x, z) == one
    for x, y, z in [("P", "V", "T"), ("H", "S", "G"), ("U", "A", "T")]:
        assert d(x, y, z) * d(y, z, x) * d(z, x, y) == one * -1


def _state_functions(e, fluid):
    def f(name):
        return lambda T, P: props.state(fluid, e, T, P)[name]
    funcs = {v: f(v) for v in ("V", "S", "U", "H", "A", "G")}
    funcs["T"] = lambda T, P: T
    funcs["P"] = lambda T, P: P
    return funcs


@pytest.mark.parametrize("kind, T0, P0", [("PR", 350.0, 10e6), ("IG", 350.0, 1e6),
                                          ("vdW", 400.0, 5e6)])
def test_all_336_derivatives_against_constrained_finite_differences(kind, T0, P0):
    """Each reduced formula is compared with a direct numerical derivative taken
    along the constraint Z = const, found by solving for (T, P); this does not
    use the Jacobian identity or the differential table."""
    e = eos(kind)
    F = _state_functions(e, CO2)
    st = props.state(CO2, e, T0, P0)
    vals = {s: st[s] for s in symbolic.SYMS}

    def grad(name):
        g = F[name]
        return ((g(T0 * (1 + 1e-6), P0) - g(T0 * (1 - 1e-6), P0)) / (2e-6 * T0),
                (g(T0, P0 * (1 + 1e-6)) - g(T0, P0 * (1 - 1e-6))) / (2e-6 * P0))

    grads = {v: grad(v) for v in F}
    worst = 0.0
    for x, y, z in symbolic.all_derivatives():
        y0, z0 = F[y](T0, P0), F[z](T0, P0)
        sy = math.hypot(grads[y][0] * T0, grads[y][1] * P0)
        sz = math.hypot(grads[z][0] * T0, grads[z][1] * P0)
        jac = grads[y][0] * grads[z][1] - grads[y][1] * grads[z][0]
        cond = abs(jac) * T0 * P0 / (sy * sz)
        if cond < 1e-3:          # y and z nearly dependent here: skip ill-posed case
            continue
        formula = symbolic.derivative(x, y, z).eval(vals)
        h = 2e-5 * sy

        def solve(dy):
            def eqs(u):
                T, P = T0 * (1 + u[0]), P0 * (1 + u[1])
                return [(F[y](T, P) - y0 - dy) / sy, (F[z](T, P) - z0) / sz]
            u = fsolve(eqs, [0.0, 0.0], xtol=1e-14)
            return F[x](T0 * (1 + u[0]), P0 * (1 + u[1]))

        numeric = (solve(h) - solve(-h)) / (2 * h)
        scale = max(abs(formula), abs(numeric), 1e-300)
        err = abs(formula - numeric) / scale
        if abs(formula) < 1e-12 * math.hypot(grads[x][0] * T0, grads[x][1] * P0) / sy:
            continue             # derivative is zero to rounding (e.g. ideal-gas (dT/dP)_H)
        worst = max(worst, err)
        assert err < 2e-4, f"(d{x}/d{y})_{z}: formula {formula}, numeric {numeric}"
    assert worst < 2e-4


# ---------------------------------------------------------------- EoS
@pytest.mark.parametrize("kind", [k for k in EOS_NAMES if k != "IG"])
@pytest.mark.parametrize("T, P", [(350.0, 10e6), (320.0, 2e6), (500.0, 30e6), (280.0, 8e6)])
def test_closed_forms_agree_with_numerical_scheme(kind, T, P):
    e = eos(kind)
    cp0 = IdealGasCp("shomate", CO2.cp0.c, 0.0, 1e9)
    fl = type(CO2)("x", "", CO2.M, CO2.Tc, CO2.Pc, CO2.omega, cp0)
    a = props.state(fl, e, T, P)
    b = NumericalModel(e.pressure, cp0, CO2.M).state(T, P, a["V"])
    for k in ("V", "alpha", "kappaT", "HR", "SR", "GR", "Cv", "Cp", "w", "muJT"):
        assert a[k] == pytest.approx(b[k], rel=5e-8), k


@pytest.mark.parametrize("kind", ["vdW", "RK", "SRK", "PR"])
def test_critical_point_conditions(kind):
    """(dP/dV)_T = 0 and P = Pc at T = Tc on the critical isotherm."""
    e = eos(kind)
    Zc = {"vdW": 0.375, "RK": 1 / 3, "SRK": 1 / 3, "PR": 0.3074}[kind]
    Vc = Zc * R * e.Tc / e.Pc
    assert e.pressure(e.Tc, Vc) == pytest.approx(e.Pc, rel=2e-3)
    assert abs(e.dPdV(e.Tc, Vc)) * Vc / e.Pc < 2e-3


def test_ideal_gas_limits():
    e = eos("IG")
    st = props.state(CO2, e, 400.0, 2e5)
    assert st["Z"] == 1 and st["HR"] == 0 and st["SR"] == 0
    assert st["Cp"] - st["Cv"] == pytest.approx(R)
    assert st["alpha"] == pytest.approx(1 / 400.0)
    assert st["kappaT"] == pytest.approx(1 / 2e5)
    assert abs(st["muJT"]) < 1e-15
    real = props.state(CO2, eos("PR"), 400.0, 10.0)      # P -> 0
    assert real["Z"] == pytest.approx(1.0, abs=1e-6)
    assert real["Cp"] == pytest.approx(st["Cp"], rel=1e-6)


def test_saturation_equal_fugacity_and_reference():
    e = eos("PR")
    P, VL, VV = e.psat(280.0)
    assert e.lnphi(280.0, VL) == pytest.approx(e.lnphi(280.0, VV), abs=1e-9)
    assert e.pressure(280.0, VL) == pytest.approx(P, rel=1e-8)
    assert P == pytest.approx(4.1607e6, rel=0.01)        # NIST value at 280 K
    with pytest.raises(OutOfRangeError):
        e.psat(310.0)


def test_state_functions_are_path_independent():
    e = eos("PR")
    a, _, _ = props.process(CO2, e, 320.0, 1e6, 380.0, 12e6)
    b, _, _ = props.process(CO2, e, 320.0, 1e6, 350.0, 20e6)
    c, _, _ = props.process(CO2, e, 350.0, 20e6, 380.0, 12e6)
    for k in ("dH", "dS", "dU", "dG"):
        assert a[k] == pytest.approx(b[k] + c[k], rel=1e-10)


def test_singular_derivative_is_reported():
    """In an ideal gas H depends on T only, so T cannot vary at constant H."""
    from thermocalculus.errors import SingularDerivativeError
    st = props.state(CO2, eos("IG"), 350.0, 1e6)
    with pytest.raises(SingularDerivativeError):
        symbolic.derivative("P", "T", "H").eval({s: st[s] for s in symbolic.SYMS})


def test_missing_parameters_are_reported():
    with pytest.raises(MissingDataError) as exc:
        EoS("PR", 304.0, 7.4e6, None)
    assert exc.value.missing == ["omega"]


# ---------------------------------------------------------------- units
def test_unit_conversions():
    q = units.parse_quantity
    assert q("25 degC", "T") == pytest.approx(298.15)
    assert q("77 degF", "T") == pytest.approx(298.15)
    assert q("1 atm", "P") == 101325.0
    assert q("10 bar", "P") == 1e6
    assert q("0.19235 L/mol", "V") == pytest.approx(1.9235e-4)
    assert q("22414 cm3/mol", "V") == pytest.approx(0.022414)
    assert q("37 J/(mol*K)", "Cp") == q("37 J/mol/K", "Cp") == q("37 J mol^-1 K^-1", "Cp")
    assert q("8.84 cal/(mol*K)", "Cp") == pytest.approx(8.84 * 4.184)
    assert q("6.0435 K/MPa", "muJT") == pytest.approx(6.0435e-6)
    assert q("44.0095 g/mol", "M") == pytest.approx(0.0440095)
    assert units.convert(1e6, "bar") == pytest.approx(10.0)
    assert units.convert(298.15, "degC") == pytest.approx(25.0)


@pytest.mark.parametrize("text, var", [("10 K", "P"), ("350", "T"), ("5 furlong", "P"),
                                       ("3 degC/mol", "T"), ("abc", "T")])
def test_unit_errors(text, var):
    with pytest.raises(UnitError):
        units.parse_quantity(text, var)


# ---------------------------------------------------------------- tables
def _pr_table(Ts, Ps):
    e = eos("PR")
    T, P, V = [], [], []
    for t in Ts:
        for p in Ps:
            T.append(t), P.append(p), V.append(e.volume(t, p))
    return PVTTable(T, P, V), e


def test_table_reproduces_the_model_that_generated_it():
    Ts = np.arange(320.0, 421.0, 5.0)
    Ps = np.array([0.05, 0.1, 0.5, 1, 2, 3, 4, 5, 6, 8, 10]) * 1e6
    tab, e = _pr_table(Ts, Ps)
    ref = props.state(CO2, e, 372.0, 7e6)
    got = tab.state(372.0, 7e6, CO2.cp0, CO2.M)
    assert got["V"] == pytest.approx(ref["V"], rel=1e-5)
    assert got["alpha"] == pytest.approx(ref["alpha"], rel=2e-3)
    assert got["kappaT"] == pytest.approx(ref["kappaT"], rel=2e-3)
    assert got["HR"] == pytest.approx(ref["HR"], rel=3e-3)
    assert got["SR"] == pytest.approx(ref["SR"], rel=5e-3)
    assert got["Cp"] == pytest.approx(ref["Cp"], rel=1e-2)


def test_table_refuses_to_extrapolate_and_checks_grid():
    tab, _ = _pr_table(np.arange(320.0, 401.0, 20.0), np.array([1, 2, 4, 8]) * 1e6)
    for T, P in [(450.0, 2e6), (350.0, 0.5e6), (350.0, 9e6), (300.0, 2e6)]:
        with pytest.raises(OutOfRangeError):
            tab.state(T, P)
    from thermocalculus.errors import InputError
    with pytest.raises(InputError):
        PVTTable([300, 300, 310], [1e5, 2e5, 1e5], [0.02, 0.01, 0.02])


def test_nist_table_gives_nist_enthalpy_change():
    tab = PVTTable.from_csv(os.path.join(DATA, "inputs", "tables", "co2_pvt_table_nist.csv"))
    s1 = tab.state(400.0, 0.1e6, CO2.cp0, CO2.M)
    s2 = tab.state(350.0, 10e6, CO2.cp0, CO2.M)
    assert s2["H"] - s1["H"] == pytest.approx((20.458 - 26.286) * 1e3, rel=5e-3)
    assert s2["S"] - s1["S"] == pytest.approx(80.159 - 132.14, rel=5e-3)


# ---------------------------------------------------------------- knowledge
def test_rules_agree_with_symbolic_reductions():
    st = props.state(CO2, eos("PR"), 350.0, 10e6)
    vals = {s: st[s] for s in symbolic.SYMS}
    d = symbolic.derivative
    assert st["Cv"] == pytest.approx(350.0 * d("S", "T", "V").eval(vals))
    assert st["muJT"] == pytest.approx(d("T", "P", "H").eval(vals))
    assert st["kappaS"] == pytest.approx(-d("V", "P", "S").eval(vals) / st["V"])
    assert st["gammaV"] == pytest.approx(d("P", "T", "V").eval(vals))


def test_inference_recovers_a_full_state():
    st = props.state(CO2, eos("PR"), 350.0, 10e6)
    given = {k: st[k] for k in ("T", "P", "V", "Cp", "w", "muJT")}
    given["M"] = CO2.M
    rep = knowledge.infer(given)
    assert rep.consistent
    for k in ("alpha", "kappaT", "kappaS", "Cv", "gamma", "gammaV", "Z"):
        assert rep.facts[k].provenance == "inferred"
        assert rep.facts[k].value == pytest.approx(st[k], rel=1e-9)


def test_missing_data_alternatives():
    rep = knowledge.infer({"T": 350.0, "P": 1e7, "V": 1.9e-4, "Cp": 85.0, "kappaT": 1.5e-7})
    assert "Z" in rep.facts and "Cv" not in rep.facts
    with pytest.raises(MissingDataError) as exc:
        knowledge.require(rep, "Cv")
    assert ["alpha"] in knowledge.alternatives("Cv", rep.facts)
    assert exc.value.missing == ["alpha"]
    with pytest.raises(MissingDataError):
        knowledge.derivative_from_facts("(dT/dP)_H", rep)


def test_inconsistent_input_is_detected():
    bad = {"T": 350.0, "V": 1.9e-4, "Cp": 85.0, "Cv": 95.0, "alpha": 0.0105, "kappaT": 1.5e-7}
    rep = knowledge.infer(bad)
    assert not rep.consistent
    assert any("smaller than Cv" in v for v in rep.violations)
    assert not knowledge.infer({"T": 300.0, "V": 1e-4, "alpha": 1e-3, "kappaT": -1e-9}).consistent
    with pytest.raises(InconsistentInputError):
        raise InconsistentInputError("x")


# ---------------------------------------------------------------- uncertainty
def test_linear_propagation_matches_analytic_and_monte_carlo():
    f = lambda x: {"y": x["a"] * x["b"] ** 2}
    x, u = {"a": 2.0, "b": 3.0}, {"a": 0.02, "b": 0.03}
    val, uy, share = propagate(f, x, u)["y"]
    exact = math.sqrt((9 * 0.02) ** 2 + (2 * 2 * 3 * 0.03) ** 2)
    assert val == 18.0 and uy == pytest.approx(exact, rel=1e-8)
    assert sum(share.values()) == pytest.approx(1.0)
    assert share["b"] == pytest.approx((12 * 0.03) ** 2 / exact ** 2)
    mc, failed = monte_carlo(f, x, u, n=20000)
    assert failed == 0 and mc["y"][1] == pytest.approx(exact, rel=0.03)


# ---------------------------------------------------------------- scenarios
EXPECTED = {
    "01_state_co2_pr": "ok", "02_process_co2_pr": "ok", "03_saturation_co2_pr": "ok",
    "04_derivatives_symbolic": "ok", "05_infer_from_measurements": "ok",
    "06_infer_missing_data": "partial", "07_infer_inconsistent": "inconsistent",
    "08_table_co2_process": "ok", "09_table_out_of_range": "error: OutOfRangeError",
    "10_unit_wrong_dimension": "error: UnitError", "11_unit_missing": "error: UnitError",
    "12_missing_cp0": "partial", "13_cp0_out_of_range": "error: OutOfRangeError",
    "14_state_n2_srk": "ok", "15_state_ideal_gas_cp_table": "ok",
    "16_unknown_fluid": "error: MissingDataError",
    "17_singular_derivative": "partial",
}


@pytest.mark.parametrize("name, status", sorted(EXPECTED.items()))
def test_scenarios_run_with_expected_status(name, status, tmp_path):
    res = scenario.run(os.path.join(DATA, "inputs", name + ".txt"))
    assert res.status == status
    out_csv, out_log = scenario.write(res, str(tmp_path))
    with open(out_csv) as fh:
        rows = list(csv.reader(l for l in fh if not l.startswith("#")))
    assert rows[0][:3] == ["quantity", "value", "unit"]
    assert os.path.getsize(out_log) > 0


def test_every_input_file_has_an_expected_status():
    names = {os.path.basename(p)[:-4] for p in glob.glob(os.path.join(DATA, "inputs", "*.txt"))}
    assert names == set(EXPECTED)


def test_inference_from_nist_values_reproduces_nist_cv():
    res = scenario.run(os.path.join(DATA, "inputs", "05_infer_from_measurements.txt"))
    assert res.get("Cv") == pytest.approx(38.325, rel=2e-3)     # NIST: 38.325 J/(mol K)


def test_mixed_units_and_output_units():
    res = scenario.run(os.path.join(DATA, "inputs", "14_state_n2_srk.txt"))
    assert res.get("T") == pytest.approx(300.0)
    assert res.get("P") == pytest.approx(10.0)                  # MPa after conversion
    assert 1e3 / res.get("V") / 1e3 == pytest.approx(3.9883, rel=0.03)  # NIST mol/L


# ---- interactive prompting -------------------------------------------------
def _scripted(answers, **kw):
    from thermocalculus import interactive
    out = []
    it = iter(answers)

    def feed(prompt):
        try:
            a = next(it)
        except StopIteration:
            raise EOFError
        out.append(prompt + a)
        return a
    sess = interactive.Session(input_fn=feed, output_fn=out.append)
    if kw.get("model"):
        sess.set_model(kw["model"])
    rep = sess.run()
    return rep, "\n".join(out)


def test_interactive_session_prompts_for_missing_data_and_reproduces_nist_cv():
    rep, text = _scripted([
        "target Cv", "T = 350 K", "P = 10 MPa", "V = 0.19235 L/mol",
        "Cp = 85.729 J/(mol*K)", "muJT = 6.0435 K/MPa", "M = 44.0095 g/mol",
        "w = 248.62 m/s", "quit"])
    assert "Cannot compute Cv yet" in text and "supply any one of" in text
    assert rep.facts["Cv"].value == pytest.approx(38.325, rel=2e-3)


def test_interactive_session_rejects_wrong_unit_and_unknown_quantity():
    rep, text = _scripted(["target Cv", "T = 350 mL", "foo = 3 K", "quit"])
    assert "rejected" in text and "not a quantity the engine knows" in text
    assert "T" not in rep.facts


def test_interactive_session_uses_eos_model_when_T_and_P_are_given():
    rep, text = _scripted(["target w", "T = 350 K", "P = 10 MPa", "quit"], model="CO2 PR")
    assert rep.facts["w"].provenance.startswith("parameterised")
    assert 150 < rep.facts["w"].value < 350


def test_interactive_session_skip_gives_up_on_target():
    rep, text = _scripted(["target Cv", "T = 350 K", "skip", "quit"])
    assert "giving up on Cv" in text
