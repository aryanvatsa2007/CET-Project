"""Regenerate every Data file, every figure, the manuscript number macros and the
website data from the engine.  Nothing in Data/ (except the NIST tables saved in
Data/reference/raw_nist and the scenario files in Data/inputs) is typed by hand.

Usage (from the project folder or from Codes/):

    python Codes/run_all.py            # everything
    python Codes/run_all.py --no-figures
"""
import csv
import glob
import math
import os
import shutil
import subprocess
import sys
import warnings

import numpy as np
from scipy.optimize import fsolve

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "Data")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "tools"))

from thermocalculus import props, scenario, symbolic  # noqa: E402
from thermocalculus.constants import R  # noqa: E402
from thermocalculus.errors import ThermoError  # noqa: E402
from thermocalculus.eos import EoS, IdealGasCp, _KINDS, get_fluid  # noqa: E402
from thermocalculus.numerics import NumericalModel  # noqa: E402
from thermocalculus.tables import PVTTable  # noqa: E402

EOS_LIST = ["IG", "vdW", "RK", "SRK", "PR"]
warnings.filterwarnings("ignore", category=RuntimeWarning)


# ------------------------------------------------------------------ helpers
def wcsv(relpath, comments, header, rows):
    path = os.path.join(DATA, relpath)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        for c in comments:
            fh.write(f"# {c}\n")
        wr = csv.writer(fh)
        wr.writerow(header)
        for r in rows:
            wr.writerow([f"{x:.8g}" if isinstance(x, float) else x for x in r])
    print(f"  wrote Data/{relpath} ({len(rows)} rows)")


def rcsv(relpath):
    with open(os.path.join(DATA, relpath), newline="") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def make_eos(kind, f):
    return EoS(kind, f.Tc, f.Pc, f.omega)


# ------------------------------------------------------------------ 1 reference
def step_reference():
    print("[1] reference data and derived input tables")
    import build_reference_data
    build_reference_data.main()
    shutil.copy(os.path.join(HERE, "thermocalculus", "data", "fluids.csv"),
                os.path.join(DATA, "inputs", "fluids.csv"))
    print("  copied the fluid database to Data/inputs/fluids.csv")


# ------------------------------------------------------------------ 2 scenarios
def step_scenarios():
    print("[2] scenarios")
    out = os.path.join(DATA, "outputs")
    os.makedirs(out, exist_ok=True)
    rows = []
    for path in sorted(glob.glob(os.path.join(DATA, "inputs", "*.txt"))):
        res = scenario.run(path)
        scenario.write(res, out, "Data/inputs/" + os.path.basename(path))
        flag = [l for l in res.log if l.split(":")[0] in
                ("MISSING DATA", "INCONSISTENT", "SINGULAR") or l.startswith("ERROR")]
        rows.append([res.name, res.kind, res.status, len(res.rows),
                     flag[0] if flag else ""])
    wcsv("outputs/scenario_summary.csv",
         ["outcome of every scenario in Data/inputs (one row per input file)",
          "status: ok | partial (some targets not computable) | inconsistent | error: <type>"],
         ["scenario", "type", "status", "n_output_rows", "first_message"], rows)
    # interactive session (the engine prompts for missing data), driven by a script
    from thermocalculus import interactive
    lines = []
    interactive.run_session(os.path.join(DATA, "inputs", "sessions", "ask_session_script.txt"),
                            output_fn=lines.append)
    with open(os.path.join(out, "ask_session_transcript.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("  wrote Data/outputs/ask_session_transcript.txt (interactive session)")


# ------------------------------------------------------------------ 3 verification
def _constrained_fd(F, grads, x, y, z, T0, P0):
    """(dx/dy)_z by stepping along z = const (2x2 nonlinear solve); None if y and z
    are nearly dependent at this state."""
    y0, z0 = F[y](T0, P0), F[z](T0, P0)
    sy = math.hypot(grads[y][0] * T0, grads[y][1] * P0)
    sz = math.hypot(grads[z][0] * T0, grads[z][1] * P0)
    jac = grads[y][0] * grads[z][1] - grads[y][1] * grads[z][0]
    if abs(jac) * T0 * P0 / (sy * sz) < 1e-3:
        return None
    h = 2e-5 * sy

    def solve(dy):
        def eqs(u):
            T, P = T0 * (1 + u[0]), P0 * (1 + u[1])
            return [(F[y](T, P) - y0 - dy) / sy, (F[z](T, P) - z0) / sz]
        u = fsolve(eqs, [0.0, 0.0], xtol=1e-14)
        return F[x](T0 * (1 + u[0]), P0 * (1 + u[1]))
    return (solve(h) - solve(-h)) / (2 * h)


def step_verification():
    print("[3] verification of the mathematics")
    co2 = get_fluid("CO2")
    d = symbolic.derivative
    # 3a every reduced derivative
    rows = [[x, y, z, str(d(x, y, z)), d(x, y, z).latex(),
             " ".join(sorted(d(x, y, z).symbols()))]
            for x, y, z in symbolic.all_derivatives()]
    wcsv("verification/all_derivatives_reduced.csv",
         ["all 336 first derivatives (dX/dY)_Z of T,P,V,S,U,H,A,G reduced to "
          "T,P,V,S,Cp,alpha,kappaT by the symbolic module"],
         ["X", "Y", "Z", "expression", "latex", "symbols_needed"], rows)
    # 3b differential table and Maxwell relations
    rows = [[v, str(symbolic.Ratio(symbolic.DIFF[v][0])), str(symbolic.Ratio(symbolic.DIFF[v][1])),
             symbolic.Ratio(symbolic.DIFF[v][0]).latex(), symbolic.Ratio(symbolic.DIFF[v][1]).latex()]
            for v in symbolic.VARIABLES]
    wcsv("verification/differential_table.csv",
         ["partial derivatives of the eight variables in the (T, P) plane used by the Jacobian reduction"],
         ["X", "(dX/dT)_P", "(dX/dP)_T", "latex_dT", "latex_dP"], rows)
    mx = symbolic.check_maxwell()
    one = symbolic.Ratio(symbolic.Poly.const(1))
    trip = symbolic.all_derivatives()
    n_rec = sum(d(x, y, z) * d(y, x, z) == one for x, y, z in trip)
    n_cyc = sum(d(x, y, z) * d(y, z, x) * d(z, x, y) == one * -1 for x, y, z in trip)
    rows = [[t, o, v, "holds" if ok else "FAILS"] for t, o, v, ok in mx]
    rows.append(["reciprocity (dX/dY)_Z (dY/dX)_Z = 1", f"{len(trip)} triples", "",
                 f"{n_rec} of {len(trip)} hold"])
    rows.append(["cyclic rule (dX/dY)_Z (dY/dZ)_X (dZ/dX)_Y = -1", f"{len(trip)} triples", "",
                 f"{n_cyc} of {len(trip)} hold"])
    wcsv("verification/identity_checks.csv",
         ["exact symbolic identities checked by polynomial cross-multiplication"],
         ["identity", "origin_or_scope", "reduced_value", "result"], rows)
    # 3c 336 derivatives against constrained finite differences
    rows = []
    for label, kind, T0, P0 in (("PR supercritical", "PR", 350.0, 10e6),
                                ("PR gas", "PR", 320.0, 2e6)):
        e = make_eos(kind, co2)
        F = {v: (lambda T, P, v=v: props.state(co2, e, T, P)[v])
             for v in ("V", "S", "U", "H", "A", "G")}
        F["T"], F["P"] = (lambda T, P: T), (lambda T, P: P)
        st = props.state(co2, e, T0, P0)
        vals = {s: st[s] for s in symbolic.SYMS}
        grads = {}
        for v, g in F.items():
            grads[v] = ((g(T0 * (1 + 1e-6), P0) - g(T0 * (1 - 1e-6), P0)) / (2e-6 * T0),
                        (g(T0, P0 * (1 + 1e-6)) - g(T0, P0 * (1 - 1e-6))) / (2e-6 * P0))
        for x, y, z in trip:
            expr = d(x, y, z)
            num = _constrained_fd(F, grads, x, y, z, T0, P0)
            if num is None:
                rows.append([label, kind, T0, P0, x, y, z, str(expr), "", "", "",
                             "skipped: Y and Z nearly dependent"])
                continue
            val = expr.eval(vals)
            dev = abs(val - num) / max(abs(val), abs(num))
            rows.append([label, kind, T0, P0, x, y, z, str(expr), val, num, dev, "checked"])
    wcsv("verification/derivatives_336_check.csv",
         ["each reduced formula evaluated with CO2 state values and compared with a direct numerical",
          "derivative taken along the constraint Z = const (central difference after a 2x2 nonlinear solve);",
          "SI units; relative_deviation = |formula - numeric| / max(|formula|, |numeric|)"],
         ["state", "eos", "T_K", "P_Pa", "X", "Y", "Z", "expression", "formula_value",
          "numeric_value", "relative_deviation", "status"], rows)
    # 3d closed forms against the generic numerical scheme
    rows = []
    cp0 = IdealGasCp("shomate", co2.cp0.c, 0.0, 1e9)       # range check off for this test
    fl = type(co2)("co2", "", co2.M, co2.Tc, co2.Pc, co2.omega, cp0)
    for kind in ("vdW", "RK", "SRK", "PR"):
        e = make_eos(kind, co2)
        for T, P in ((350.0, 10e6), (320.0, 2e6), (500.0, 30e6), (280.0, 8e6)):
            a = props.state(fl, e, T, P)
            b = NumericalModel(e.pressure, cp0, co2.M).state(T, P, a["V"])
            for k in ("alpha", "kappaT", "HR", "SR", "GR", "Cv", "Cp", "w", "muJT"):
                rows.append([kind, T, P, k, a[k], b[k], abs(a[k] - b[k]) / abs(a[k])])
    wcsv("verification/closedform_vs_numerical.csv",
         ["closed-form cubic-EoS results against the model-independent scheme (Richardson finite differences",
          "and Gauss-Kronrod quadrature of the residual Helmholtz energy); CO2 parameters; SI units"],
         ["eos", "T_K", "P_Pa", "property", "closed_form", "numerical_scheme",
          "relative_difference"], rows)
    # 3e thermodynamic consistency of the state functions
    rows = []
    for kind in EOS_LIST:
        e = make_eos(kind, co2)
        for T, P in ((350.0, 10e6), (320.0, 2e6), (400.0, 20e6)):
            st = props.state(co2, e, T, P)
            g = lambda t, p: props.state(co2, e, t, p)["G"]
            h = lambda t, p: props.state(co2, e, t, p)["H"]
            gT = (g(T * (1 + 1e-5), P) - g(T * (1 - 1e-5), P)) / (2e-5 * T)
            gP = (g(T, P * (1 + 1e-5)) - g(T, P * (1 - 1e-5))) / (2e-5 * P)
            hT = (h(T * (1 + 1e-5), P) - h(T * (1 - 1e-5), P)) / (2e-5 * T)
            rows.append([kind, T, P, "(dG/dT)_P = -S", gT, -st["S"], abs(gT + st["S"]) / abs(st["S"])])
            rows.append([kind, T, P, "(dG/dP)_T = V", gP, st["V"], abs(gP - st["V"]) / st["V"]])
            rows.append([kind, T, P, "(dH/dT)_P = Cp", hT, st["Cp"], abs(hT - st["Cp"]) / st["Cp"]])
    wcsv("verification/consistency_checks.csv",
         ["numerical derivatives of the computed G and H surfaces against the independently computed S, V and Cp"],
         ["eos", "T_K", "P_Pa", "identity", "numerical_derivative", "direct_value",
          "relative_difference"], rows)
    # 3f critical point and constants
    rows = []
    for kind, (eps, sig, oa, ob) in _KINDS.items():
        e = make_eos(kind, co2)
        Zc = {"vdW": 0.375, "RK": 1 / 3, "SRK": 1 / 3, "PR": 0.3074}[kind]
        Vc = Zc * R * e.Tc / e.Pc
        rows.append([kind, eps, sig, oa, ob, getattr(e, "m", ""), e.ac, e.b,
                     e.pressure(e.Tc, Vc) / e.Pc - 1.0, e.dPdV(e.Tc, Vc) * Vc / e.Pc])
    wcsv("verification/eos_constants.csv",
         ["constants of the generic cubic form P = RT/(V-b) - a(T)/((V+eps b)(V+sig b)) and CO2 parameters;",
          "last two columns: residuals of the critical-point conditions P(Tc,Vc)=Pc and (dP/dV)_T=0"],
         ["eos", "eps", "sig", "Omega_a", "Omega_b", "m_CO2", "a_c_Pa_m6_mol-2",
          "b_m3_mol-1", "P(Tc,Vc)/Pc-1", "Vc/Pc*(dP/dV)_Tc"], rows)


# ------------------------------------------------------------------ 4 validation
PROPS = [("rho", "mol/L"), ("Cv", "J/(mol*K)"), ("Cp", "J/(mol*K)"), ("w", "m/s"),
         ("muJT", "K/MPa"), ("dH", "kJ/mol"), ("dS", "J/(mol*K)")]
NIST_COL = {"rho": "rho [mol/L]", "Cv": "Cv [J/(mol*K)]", "Cp": "Cp [J/(mol*K)]",
            "w": "w [m/s]", "muJT": "muJT [K/MPa]", "dH": "H [kJ/mol]", "dS": "S [J/(mol*K)]"}


def _model_props(st, ref):
    return {"rho": 1e-3 / st["V"], "Cv": st["Cv"], "Cp": st["Cp"],
            "w": st.get("w", float("nan")),
            "muJT": st["muJT"] * 1e6, "dH": (st["H"] - ref["H"]) * 1e-3,
            "dS": st["S"] - ref["S"]}


def _compare(rows_out, tag, T, P, phase, model, nist_row, nist_ref):
    for k, unit in PROPS:
        nv = float(nist_row[NIST_COL[k]])
        if k in ("dH", "dS"):
            nv -= float(nist_ref[NIST_COL[k]])
        dev = model[k] - nv
        pct = 100.0 * dev / nv if (k not in ("dH", "dS") and nv != 0) else float("nan")
        rows_out.append([T, P, phase, tag, k, unit, model[k], nv, dev, pct])


def _region(P):
    return "0.1-4 MPa" if P <= 4.01 else ("6-10 MPa" if P <= 10.01 else "15-20 MPa")


def step_validation():
    print("[4] validation against NIST reference data")
    co2, n2 = get_fluid("CO2"), get_fluid("N2")
    nist = rcsv("reference/nist_co2_isobars_310-400K_0.1-20MPa.csv")
    key = lambda r: (float(r["T [K]"]), float(r["P [MPa]"]))
    ref_row = next(r for r in nist if key(r) == (400.0, 0.1))
    head = ["T_K", "P_MPa", "nist_phase", "model", "property", "unit", "model_value",
            "nist_value", "deviation", "deviation_percent"]
    note = ["deviation = model - NIST; deviation_percent = 100 * deviation / NIST",
            "dH and dS are changes from the reference state named below (the engine and NIST use different",
            "zero points); for them only the absolute deviation is meaningful"]
    # 4a cubic equations of state, CO2
    rows = []
    for kind in EOS_LIST:
        e = make_eos(kind, co2)
        ref = props.state(co2, e, 400.0, 0.1e6)
        for r in nist:
            T, P = key(r)
            st = props.state(co2, e, T, P * 1e6)
            _compare(rows, kind, T, P, r["phase"], _model_props(st, ref), r, ref_row)
    wcsv("validation/co2_eos_vs_nist.csv",
         ["carbon dioxide: ideal-gas and cubic equations of state against NIST WebBook values (Span-Wagner)",
          "171 state points, 310-400 K, 0.1-20 MPa; reference state for dH, dS: 400 K, 0.1 MPa"] + note,
         head, rows)
    eos_rows = rows
    # 4a' smooth model curves for the figures (same engine calls, denser abscissa)
    rows = []
    for kind in EOS_LIST:
        e = make_eos(kind, co2)
        for T in (310.0, 350.0, 400.0):
            for P in np.arange(0.1, 20.01, 0.1):
                st = props.state(co2, e, T, P * 1e6)
                rows.append(["isotherm", kind, T, float(P), st["Z"], 1e-3 / st["V"],
                             st["Cp"], st["w"], st["muJT"] * 1e6])
        for T in np.arange(310.0, 400.01, 0.5):
            st = props.state(co2, e, float(T), 10e6)
            rows.append(["isobar", kind, float(T), 10.0, st["Z"], 1e-3 / st["V"],
                         st["Cp"], st["w"], st["muJT"] * 1e6])
    wcsv("validation/co2_model_curves.csv",
         ["model curves for Figures 3 and 5: three isotherms (0.1-20 MPa) and the 10 MPa isobar"],
         ["curve", "model", "T_K", "P_MPa", "Z", "rho_mol_per_L", "Cp_J_per_molK",
          "w_m_per_s", "muJT_K_per_MPa"], rows)
    rows = []
    for kind in ("vdW", "RK", "SRK", "PR"):
        e = make_eos(kind, co2)
        for T in np.arange(220.0, 303.01, 1.0):
            P, VL, VV = e.psat(float(T))
            rows.append([kind, float(T), P * 1e-6, 1e-3 / VL, 1e-3 / VV])
    wcsv("validation/co2_saturation_curves.csv",
         ["model saturation curves for Figure 6 (equal-fugacity solution every 1 K)"],
         ["model", "T_K", "Psat_MPa", "rhoL_mol_per_L", "rhoV_mol_per_L"], rows)
    # 4b nitrogen isotherm
    nn = rcsv("reference/nist_n2_isotherm_300K_1-31MPa.csv")
    rows = []
    for kind in EOS_LIST:
        e = make_eos(kind, n2)
        ref = props.state(n2, e, 300.0, 1e6)
        for r in nn:
            T, P = key(r)
            st = props.state(n2, e, T, P * 1e6)
            _compare(rows, kind, T, P, r["phase"], _model_props(st, ref), r, nn[0])
    wcsv("validation/n2_eos_vs_nist.csv",
         ["nitrogen at 300 K, 1-31 MPa: equations of state against NIST WebBook values (Span et al.)",
          "reference state for dH, dS: 300 K, 1 MPa"] + note, head, rows)
    n2_rows = rows
    # 4c saturation
    sat = rcsv("reference/nist_co2_saturation_220-300K.csv")
    liq = {float(r["T [K]"]): r for r in sat if r["phase"] == "liquid"}
    vap = {float(r["T [K]"]): r for r in sat if r["phase"] == "vapor"}
    rows = []
    for kind in ("vdW", "RK", "SRK", "PR"):
        e = make_eos(kind, co2)
        for T in sorted(liq):
            P, VL, VV = e.psat(T)
            rl, rv = e.residual(T, VL), e.residual(T, VV)
            dH = ((rv["UR"] + P * VV) - (rl["UR"] + P * VL)) * 1e-3
            ref = {"Psat": float(liq[T]["P [MPa]"]), "rhoL": float(liq[T]["rho [mol/L]"]),
                   "rhoV": float(vap[T]["rho [mol/L]"]),
                   "dHvap": float(vap[T]["H [kJ/mol]"]) - float(liq[T]["H [kJ/mol]"])}
            mod = {"Psat": P * 1e-6, "rhoL": 1e-3 / VL, "rhoV": 1e-3 / VV, "dHvap": dH}
            for k, unit in (("Psat", "MPa"), ("rhoL", "mol/L"), ("rhoV", "mol/L"),
                            ("dHvap", "kJ/mol")):
                rows.append([T, kind, k, unit, mod[k], ref[k], mod[k] - ref[k],
                             100 * (mod[k] / ref[k] - 1)])
    wcsv("validation/co2_saturation_vs_nist.csv",
         ["carbon dioxide vapour-liquid saturation, 220-300 K: equal-fugacity solution of each cubic",
          "equation against NIST WebBook values; dHvap needs no ideal-gas heat capacity"],
         ["T_K", "model", "property", "unit", "model_value", "nist_value", "deviation",
          "deviation_percent"], rows)
    sat_rows = rows
    # 4d tabulated route
    tab_full = PVTTable.from_csv(os.path.join(DATA, "inputs", "tables", "co2_pvt_table_nist.csv"))
    Tn = np.array([key(r)[0] for r in nist])
    Pn = np.array([key(r)[1] for r in nist]) * 1e6
    Vn = np.array([float(r["V [L/mol]"]) for r in nist]) * 1e-3
    coarse = np.isclose(Tn % 10, 0)
    tab_coarse = PVTTable(Tn[coarse], Pn[coarse], Vn[coarse])
    rows = []
    for tag, tab, pick in (("table, all 19 temperatures", tab_full, lambda T: True),
                           ("table, 10 K spacing, held-out nodes", tab_coarse,
                            lambda T: T % 10 != 0)):
        ref = tab.state(400.0, 0.1e6, co2.cp0, co2.M)
        for r in nist:
            T, P = key(r)
            if not pick(T):
                continue
            try:
                model = _model_props(tab.state(T, P * 1e6, co2.cp0, co2.M), ref)
            except ThermoError:      # the engine refuses: spline violates stability
                model = {k: float("nan") for k, _ in PROPS}
            _compare(rows, tag, T, P, r["phase"], model, r, ref_row)
    wcsv("validation/co2_table_route_vs_nist.csv",
         ["tabulated-input route: only the NIST molar volumes on a (T, P) grid are given to the engine;",
          "all other properties follow from spline derivatives/integrals and the ideal-gas Cp correlation.",
          "'all 19 temperatures': evaluated at the grid nodes (tests derivatives and integrals);",
          "model_value = nan: the engine refused the state (spline violates kappaT > 0) or could not form the",
          "property (Cv <= 0); such points are counted in summary_deviations.csv",
          "'10 K spacing, held-out nodes': table built from every second temperature and evaluated at the",
          "omitted temperatures (tests interpolation as well); reference state for dH, dS: 400 K, 0.1 MPa"]
         + note, head, rows)
    tab_rows = rows
    # 4e summary
    out = []

    def summarise(dataset, rows, region_of):
        groups = {}
        for T, P, ph, model, k, unit, mv, nv, dev, pct in rows:
            use_abs = k in ("dH", "dS")
            groups.setdefault((model, region_of(P), k, unit if use_abs else "%"), []).append(
                abs(dev) if use_abs else abs(pct))
            groups.setdefault((model, "all", k, unit if use_abs else "%"), []).append(
                abs(dev) if use_abs else abs(pct))
        for (model, region, k, unit), v in groups.items():
            ok = [x for x in v if not math.isnan(x)]     # NaN: refused or not computable
            out.append([dataset, model, region, k, len(ok), float(np.mean(ok)),
                        float(np.max(ok)), unit, len(v) - len(ok)])

    summarise("CO2 310-400 K", eos_rows, _region)
    summarise("N2 300 K", n2_rows, lambda P: "1-31 MPa")
    summarise("CO2 tabulated route", tab_rows, _region)
    groups = {}
    for T, model, k, unit, mv, nv, dev, pct in sat_rows:
        groups.setdefault((model, k), []).append(abs(pct))
    for (model, k), v in groups.items():
        out.append(["CO2 saturation 220-300 K", model, "all", k, len(v),
                    float(np.mean(v)), float(np.max(v)), "%", 0])
    out = [r for r in out if not (r[0] == "N2 300 K" and r[2] == "all")]
    wcsv("validation/summary_deviations.csv",
         ["average absolute deviation (AAD) and maximum absolute deviation from NIST, by data set, model,",
          "pressure region and property; unit '%' means percent deviation, otherwise absolute deviation in that unit"],
         ["dataset", "model", "region", "property", "n_points", "AAD", "max_abs_deviation",
          "unit", "n_refused_or_not_computable"],
         out)


def main():
    step_reference()
    step_scenarios()
    step_verification()
    step_validation()
    import products
    products.step_numbers(ROOT)
    products.step_website_data(ROOT)
    if "--no-figures" not in sys.argv:
        import figures
        figures.make_all(ROOT)
    products.step_index(ROOT)
    print("done")


if __name__ == "__main__":
    main()
