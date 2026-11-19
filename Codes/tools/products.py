"""Products derived from the Data files: LaTeX number macros and table bodies for
the manuscript, the website data file, and the Data index.  Called by run_all.py.

Every number quoted in the manuscript is a macro defined in generated_numbers.tex,
and Data/manuscript_numbers.csv records which Data file each macro comes from.
"""
import csv
import glob
import json
import math
import os
import platform
import re
import subprocess
import sys

import numpy as np


def rcsv(root, rel):
    with open(os.path.join(root, "Data", rel), newline="") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def sci(x, digits=1):
    """1.5e-05 -> 1.5 \\times 10^{-5} (math mode)."""
    m, e = f"{x:.{digits}e}".split("e")
    return f"{m} \\times 10^{{{int(e)}}}"


def tex_escape(s):
    return s.replace("_", r"\_").replace("%", r"\%").replace("&", r"\&").replace("^", r"\^{}")


# --------------------------------------------------------------------------
def step_numbers(root):
    print("[5] manuscript numbers and tables")
    N = []          # (macro, text, source, meaning)

    def add(macro, text, source, meaning):
        assert re.fullmatch(r"[A-Za-z]+", macro), macro
        N.append((macro, str(text), source, meaning))

    # ---- verification
    f = "verification/derivatives_336_check.csv"
    rows = rcsv(root, f)
    chk = [r for r in rows if r["status"] == "checked"]
    add("NDerivTotal", len(rcsv(root, "verification/all_derivatives_reduced.csv")),
        "verification/all_derivatives_reduced.csv", "number of reduced derivatives")
    allr = rcsv(root, "verification/all_derivatives_reduced.csv")
    add("NDerivNoS", sum("S" not in r["symbols_needed"].split() for r in allr),
        "verification/all_derivatives_reduced.csv", "derivatives that do not need S")
    add("NDerivAllSeven", sum(len(r["symbols_needed"].split()) == 7 for r in allr),
        "verification/all_derivatives_reduced.csv", "derivatives needing all seven symbols")
    add("NDerivChecked", len(chk), f, "derivative checks performed (two states)")
    add("NDerivSkipped", len(rows) - len(chk), f, "derivative checks skipped")
    for tag, state in (("SC", "PR supercritical"), ("Gas", "PR gas")):
        d = [float(r["relative_deviation"]) for r in chk if r["state"] == state]
        add("DerivMaxDev" + tag, sci(max(d)), f, f"max relative deviation, {state}")
        add("DerivMedDev" + tag, sci(float(np.median(d))), f, f"median relative deviation, {state}")
        add("DerivNBelow" + tag, sum(x < 1e-6 for x in d), f, f"count below 1e-6, {state}")
    idn = rcsv(root, "verification/identity_checks.csv")
    add("NMaxwell", sum(r["result"] == "holds" for r in idn), "verification/identity_checks.csv",
        "Maxwell relations confirmed")
    add("NRecip", idn[4]["result"].split()[0], "verification/identity_checks.csv",
        "reciprocity identities confirmed")
    add("NCyclic", idn[5]["result"].split()[0], "verification/identity_checks.csv",
        "cyclic-rule identities confirmed")
    f = "verification/closedform_vs_numerical.csv"
    d = [float(r["relative_difference"]) for r in rcsv(root, f)]
    add("ClosedNumMax", sci(max(d)), f, "max closed-form vs numerical difference")
    add("ClosedNumN", len(d), f, "number of closed-form vs numerical comparisons")
    f = "verification/consistency_checks.csv"
    d = [float(r["relative_difference"]) for r in rcsv(root, f)]
    add("ConsistMax", sci(max(d)), f, "max residual of G and H surface identities")
    add("ConsistN", len(d), f, "number of consistency comparisons")
    ec = {r["eos"]: r for r in rcsv(root, "verification/eos_constants.csv")}
    add("PRCritResid", sci(abs(float(ec["PR"]["P(Tc,Vc)/Pc-1"]))),
        "verification/eos_constants.csv", "PR critical-pressure residual")

    # ---- validation
    f = "validation/summary_deviations.csv"
    S = rcsv(root, f)

    def aad(ds, model, region, prop, col="AAD"):
        r = [x for x in S if x["dataset"] == ds and x["model"].startswith(model)
             and x["region"] == region and x["property"] == prop]
        assert len(r) == 1, (ds, model, region, prop, len(r))
        return float(r[0][col])

    C = "CO2 310-400 K"
    names = {"PR": "PR", "SRK": "SRK", "RK": "RK", "vdW": "vdW", "IG": "IG"}
    for m in names:
        add(f"{m}rhoAAD", f"{aad(C, m, 'all', 'rho'):.1f}", f, f"{m} density AAD %, all CO2 points")
        add(f"{m}rhoMax", f"{aad(C, m, 'all', 'rho', 'max_abs_deviation'):.0f}", f,
            f"{m} density max deviation %")
    for m in ("PR", "SRK"):
        add(f"{m}rhoAADlow", f"{aad(C, m, '0.1-4 MPa', 'rho'):.1f}", f, f"{m} density AAD %, 0.1-4 MPa")
        add(f"{m}rhoAADmid", f"{aad(C, m, '6-10 MPa', 'rho'):.1f}", f, f"{m} density AAD %, 6-10 MPa")
        add(f"{m}rhoAADhigh", f"{aad(C, m, '15-20 MPa', 'rho'):.1f}", f, f"{m} density AAD %, 15-20 MPa")
        add(f"{m}CpAAD", f"{aad(C, m, 'all', 'Cp'):.1f}", f, f"{m} Cp AAD %")
        add(f"{m}CpMax", f"{aad(C, m, 'all', 'Cp', 'max_abs_deviation'):.0f}", f, f"{m} Cp max %")
        add(f"{m}CvAAD", f"{aad(C, m, 'all', 'Cv'):.1f}", f, f"{m} Cv AAD %")
        add(f"{m}wAAD", f"{aad(C, m, 'all', 'w'):.1f}", f, f"{m} sound speed AAD %")
        add(f"{m}muAAD", f"{aad(C, m, 'all', 'muJT'):.1f}", f, f"{m} Joule-Thomson AAD %")
        add(f"{m}dHAAD", f"{1e3 * aad(C, m, 'all', 'dH'):.0f}", f, f"{m} enthalpy-change AAD, J/mol")
        add(f"{m}dSAAD", f"{aad(C, m, 'all', 'dS'):.2f}", f, f"{m} entropy-change AAD, J/(mol K)")
    add("IGdHAAD", f"{1e3 * aad(C, 'IG', 'all', 'dH'):.0f}", f, "ideal-gas enthalpy-change AAD, J/mol")
    add("vdWCvAAD", f"{aad(C, 'vdW', 'all', 'Cv'):.1f}", f, "vdW Cv AAD %")
    add("IGCvAAD", f"{aad(C, 'IG', 'all', 'Cv'):.1f}", f, "ideal-gas Cv AAD %")
    Sn = "CO2 saturation 220-300 K"
    for m in ("PR", "SRK", "RK", "vdW"):
        add(f"{m}PsatAAD", f"{aad(Sn, m, 'all', 'Psat'):.1f}", f, f"{m} saturation pressure AAD %")
        add(f"{m}rhoLAAD", f"{aad(Sn, m, 'all', 'rhoL'):.1f}", f, f"{m} saturated liquid density AAD %")
    add("PRdHvapAAD", f"{aad(Sn, 'PR', 'all', 'dHvap'):.1f}", f, "PR enthalpy of vaporisation AAD %")
    add("PRdHvapMax", f"{aad(Sn, 'PR', 'all', 'dHvap', 'max_abs_deviation'):.0f}", f,
        "PR enthalpy of vaporisation max %")
    add("PRrhoVAAD", f"{aad(Sn, 'PR', 'all', 'rhoV'):.1f}", f, "PR saturated vapour density AAD %")
    Nn = "N2 300 K"
    for m in ("PR", "SRK"):
        add(f"Ntwo{m}rhoAAD", f"{aad(Nn, m, '1-31 MPa', 'rho'):.1f}", f, f"N2 {m} density AAD %")
        add(f"Ntwo{m}CpAAD", f"{aad(Nn, m, '1-31 MPa', 'Cp'):.1f}", f, f"N2 {m} Cp AAD %")
        add(f"Ntwo{m}wAAD", f"{aad(Nn, m, '1-31 MPa', 'w'):.1f}", f, f"N2 {m} sound speed AAD %")
    Tn = "CO2 tabulated route"
    full, held = "table, all", "table, 10 K"
    add("TabdHAAD", f"{1e3 * aad(Tn, full, 'all', 'dH'):.0f}", f, "table route enthalpy AAD, J/mol")
    add("TabdSAAD", f"{aad(Tn, full, 'all', 'dS'):.2f}", f, "table route entropy AAD, J/(mol K)")
    add("TabdHAADlow", f"{1e3 * aad(Tn, full, '0.1-4 MPa', 'dH'):.0f}", f, "table route enthalpy AAD 0.1-4 MPa, J/mol")
    add("TabdHMax", f"{1e3 * aad(Tn, full, 'all', 'dH', 'max_abs_deviation'):.0f}", f, "table route enthalpy max, J/mol")
    add("TabCpAAD", f"{aad(Tn, full, 'all', 'Cp'):.0f}", f, "table route Cp AAD %")
    add("TabCvAAD", f"{aad(Tn, full, 'all', 'Cv'):.0f}", f, "table route Cv AAD %")
    add("TabCpAADlow", f"{aad(Tn, full, '0.1-4 MPa', 'Cp'):.0f}", f, "table route Cp AAD %, 0.1-4 MPa")
    add("TabHeldCpAADlow", f"{aad(Tn, held, '0.1-4 MPa', 'Cp'):.1f}", f, "held-out table Cp AAD %, 0.1-4 MPa")
    add("TabHeldrhoAAD", f"{aad(Tn, held, 'all', 'rho'):.2f}", f, "held-out table density AAD %")
    add("TabHeldrhoMax", f"{aad(Tn, held, 'all', 'rho', 'max_abs_deviation'):.1f}", f, "held-out table density max %")
    add("TabHelddHAAD", f"{1e3 * aad(Tn, held, 'all', 'dH'):.0f}", f, "held-out table enthalpy AAD, J/mol")
    tr = rcsv(root, "validation/co2_table_route_vs_nist.csv")
    refused = {(r["T_K"], r["P_MPa"]) for r in tr if r["model"].startswith(full)
               and r["property"] == "rho" and r["model_value"] == "nan"}
    add("TabNRefused", len(refused), "validation/co2_table_route_vs_nist.csv",
        "grid nodes refused by the table route (kappaT <= 0)")
    add("TabNNodes", len({(r["T_K"], r["P_MPa"]) for r in tr if r["model"].startswith(full)}),
        "validation/co2_table_route_vs_nist.csv", "grid nodes tested, full table")
    co2ref = rcsv(root, "reference/nist_co2_isobars_310-400K_0.1-20MPa.csv")
    nref = len(co2ref) + len(rcsv(root, "reference/nist_co2_saturation_220-300K.csv")) \
        + len(rcsv(root, "reference/nist_n2_isotherm_300K_1-31MPa.csv"))
    add("NNistCOtwo", len(co2ref), "reference/nist_co2_isobars_310-400K_0.1-20MPa.csv",
        "CO2 single-phase reference points")
    add("NNistTotal", nref, "reference/*.csv", "all NIST reference state points")

    # ---- scenarios
    f = "outputs/scenario_summary.csv"
    sc = rcsv(root, f)
    add("NScen", len(sc), f, "number of scenarios")
    add("NScenOK", sum(r["status"] == "ok" for r in sc), f, "scenarios completing normally")
    add("NScenPartial", sum(r["status"] == "partial" for r in sc), f, "scenarios with partial results")
    add("NScenIncons", sum(r["status"] == "inconsistent" for r in sc), f, "scenarios flagged inconsistent")
    add("NScenError", sum(r["status"].startswith("error") for r in sc), f, "scenarios refused with an error")

    def out(name):
        return {r["quantity"]: r for r in rcsv(root, f"outputs/{name}_output.csv")}

    def nist_at(T, P):
        return next(r for r in co2ref if float(r["T [K]"]) == T and float(r["P [MPa]"]) == P)

    o5, f5 = out("05_infer_from_measurements"), "outputs/05_infer_from_measurements_output.csv"
    add("InfCv", f"{float(o5['Cv']['value']):.2f}", f5, "Cv inferred from stated measurements")
    add("InfCvU", f"{float(o5['Cv']['std_uncertainty']):.2f}", f5, "its standard uncertainty")
    add("InfCvNist", nist_at(350.0, 10.0)["Cv [J/(mol*K)]"],
        "reference/nist_co2_isobars_310-400K_0.1-20MPa.csv", "NIST Cv at 350 K, 10 MPa")
    add("InfKappaT", f"{1e6 * float(o5['kappaT']['value']):.4f}", f5, "inferred kappaT, 1/MPa")
    add("InfAlpha", f"{1e3 * float(o5['alpha']['value']):.3f}", f5, "inferred alpha, 1e-3 /K")
    o8, f8 = out("08_table_co2_process"), "outputs/08_table_co2_process_output.csv"
    add("TabScendH", f"{float(o8['dH']['value']):.0f}", f8, "table scenario enthalpy change, J/mol")
    add("TabScendS", f"{float(o8['dS']['value']):.2f}", f8, "table scenario entropy change")
    a, b = nist_at(350.0, 10.0), nist_at(400.0, 0.1)
    add("TabScendHNist", f"{1e3 * (float(a['H [kJ/mol]']) - float(b['H [kJ/mol]'])):.0f}",
        "reference/nist_co2_isobars_310-400K_0.1-20MPa.csv", "NIST enthalpy change for the same states")
    add("TabScendSNist", f"{float(a['S [J/(mol*K)]']) - float(b['S [J/(mol*K)]']):.2f}",
        "reference/nist_co2_isobars_310-400K_0.1-20MPa.csv", "NIST entropy change for the same states")
    o3, f3 = out("03_saturation_co2_pr"), "outputs/03_saturation_co2_pr_output.csv"
    add("SatP", f"{1e-6 * float(o3['Psat']['value']):.3f}", f3, "PR saturation pressure at 280 K, MPa")
    add("SatPU", f"{1e-6 * float(o3['Psat']['std_uncertainty']):.3f}", f3, "its standard uncertainty, MPa")
    sat = rcsv(root, "reference/nist_co2_saturation_220-300K.csv")
    add("SatPNist", next(r["P [MPa]"] for r in sat if float(r["T [K]"]) == 280.0),
        "reference/nist_co2_saturation_220-300K.csv", "NIST saturation pressure at 280 K, MPa")
    o1, f1 = out("01_state_co2_pr"), "outputs/01_state_co2_pr_output.csv"
    add("OneRho", f"{1e-3 / float(o1['V']['value']):.3f}", f1, "PR molar density at 350 K, 10 MPa, mol/L")
    add("OneRhoNist", nist_at(350.0, 10.0)["rho [mol/L]"],
        "reference/nist_co2_isobars_310-400K_0.1-20MPa.csv", "NIST density at 350 K, 10 MPa")
    add("OneCp", f"{float(o1['Cp']['value']):.2f}", f1, "PR Cp at 350 K, 10 MPa")
    add("OneCpNist", nist_at(350.0, 10.0)["Cp [J/(mol*K)]"],
        "reference/nist_co2_isobars_310-400K_0.1-20MPa.csv", "NIST Cp at 350 K, 10 MPa")
    add("OneCpU", f"{float(o1['Cp']['std_uncertainty']):.2f}", f1, "its standard uncertainty")
    add("OneMu", f"{1e6 * float(o1['muJT']['value']):.3f}", f1, "PR Joule-Thomson coefficient, K/MPa")
    add("OneMuU", f"{1e6 * float(o1['muJT']['std_uncertainty']):.3f}", f1, "its standard uncertainty")
    ratios = []
    for nm in ("01_state_co2_pr", "02_process_co2_pr"):
        for r in rcsv(root, f"outputs/{nm}_output.csv"):
            if r["mc_std_uncertainty"] and r["provenance"] != "stated":
                ratios.append(float(r["mc_std_uncertainty"]) / float(r["std_uncertainty"]))
    add("MCRatioMin", f"{min(ratios):.2f}", "outputs/0[12]_*_output.csv", "min Monte Carlo / linear u")
    add("MCRatioMax", f"{max(ratios):.2f}", "outputs/0[12]_*_output.csv", "max Monte Carlo / linear u")
    add("MCN", len(ratios), "outputs/0[12]_*_output.csv", "outputs compared")
    o2, f2 = out("02_process_co2_pr"), "outputs/02_process_co2_pr_output.csv"
    add("TwodH", f"{float(o2['dH']['value']):.0f}", f2, "process enthalpy change, J/mol")
    add("TwodHU", f"{float(o2['dH']['std_uncertainty']):.0f}", f2, "its standard uncertainty")
    add("TwodS", f"{float(o2['dS']['value']):.2f}", f2, "process entropy change")
    add("TwodSU", f"{float(o2['dS']['std_uncertainty']):.2f}", f2, "its standard uncertainty")
    bud = {r["quantity"]: r for r in rcsv(root, "outputs/02_process_co2_pr_uncertainty_budget.csv")}
    fb = "outputs/02_process_co2_pr_uncertainty_budget.csv"
    add("TwodHShareTtwo", f"{100 * float(bud['dH']['share_T2']):.0f}", fb, "share of T2 in var(dH), %")
    add("TwodHShareOmega", f"{100 * float(bud['dH']['share_omega']):.0f}", fb, "share of omega in var(dH), %")
    add("TwodSSharePone", f"{100 * float(bud['dS']['share_P1']):.0f}", fb, "share of P1 in var(dS), %")
    add("TwodSShareTtwo", f"{100 * float(bud['dS']['share_T2']):.0f}", fb, "share of T2 in var(dS), %")
    b5 = {r["quantity"]: r for r in
          rcsv(root, "outputs/05_infer_from_measurements_uncertainty_budget.csv")}
    fb5 = "outputs/05_infer_from_measurements_uncertainty_budget.csv"
    add("InfCvShareCp", f"{100 * float(b5['Cv']['share_Cp']):.0f}", fb5, "share of Cp in var(Cv), %")
    add("InfCvShareMu", f"{100 * float(b5['Cv']['share_muJT']):.0f}", fb5, "share of muJT in var(Cv), %")

    # ---- software
    import matplotlib
    import scipy
    codes = os.path.join(root, "Codes")
    res = subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=codes,
                         capture_output=True, text=True)
    m = re.search(r"(\d+) passed", res.stdout)
    failed = re.search(r"(\d+) failed", res.stdout)
    with open(os.path.join(root, "Data", "verification", "test_suite_result.txt"), "w") as fh:
        fh.write("# result of 'python -m pytest -q' run from Codes/ by run_all.py\n")
        fh.write(res.stdout[-600:])
    add("NTests", m.group(1) if m else "0", "verification/test_suite_result.txt", "tests passed")
    add("NTestsFailed", failed.group(1) if failed else "0", "verification/test_suite_result.txt", "tests failed")
    add("PyVer", platform.python_version(), "runtime", "Python version used")
    add("NumpyVer", np.__version__, "runtime", "NumPy version used")
    add("ScipyVer", scipy.__version__, "runtime", "SciPy version used")
    add("MplVer", matplotlib.__version__, "runtime", "Matplotlib version used")
    n_lines = sum(len(open(p).read().splitlines())
                  for p in glob.glob(os.path.join(codes, "thermocalculus", "*.py")))
    add("NCodeLines", n_lines, "Codes/thermocalculus/*.py", "lines in the engine package")

    with open(os.path.join(root, "generated_numbers.tex"), "w") as fh:
        fh.write("% GENERATED by Codes/run_all.py -- do not edit.  Every number quoted in the\n"
                 "% manuscript is one of these macros; Data/manuscript_numbers.csv gives its source.\n")
        for macro, text, _, meaning in N:
            fh.write(f"\\newcommand{{\\{macro}}}{{{text}}}% {meaning}\n")
    with open(os.path.join(root, "Data", "manuscript_numbers.csv"), "w", newline="") as fh:
        fh.write("# every number quoted in manuscript.tex: LaTeX macro, value as typeset, source file "
                 "under Data/, meaning\n")
        wr = csv.writer(fh)
        wr.writerow(["macro", "value_as_typeset", "source", "meaning"])
        wr.writerows(N)
    print(f"  wrote generated_numbers.tex and Data/manuscript_numbers.csv ({len(N)} numbers)")
    _tables(root, S, sc)


def _tables(root, S, sc):
    T = {}
    # differentials
    d = rcsv(root, "verification/differential_table.csv")
    T["TabDiffBody"] = "\n".join(
        f"${r['X']}$ & ${r['latex_dT']}$ & ${r['latex_dP']}$ \\\\" for r in d)
    # selected derivatives
    sel = rcsv(root, "outputs/04_derivatives_symbolic_output.csv")
    allr = {(r["X"], r["Y"], r["Z"]): r for r in rcsv(root, "verification/all_derivatives_reduced.csv")}
    sym = {"T": "T", "P": "P", "V": "V", "S": "S", "Cp": "C_P", "alpha": r"\alpha",
           "kappaT": r"\kappa_T"}
    lines = []
    for r in sel:
        m = re.match(r"\(d(\w)/d(\w)\)_(\w)", r["quantity"])
        x, y, z = m.groups()
        a = allr[(x, y, z)]
        need = ", ".join(sym[s] for s in ("T", "P", "V", "S", "Cp", "alpha", "kappaT")
                         if s in a["symbols_needed"].split()) or "\\text{none}"
        lines.append(f"$\\left(\\partial {x}/\\partial {y}\\right)_{{{z}}}$ & "
                     f"$\\displaystyle {a['latex']}$ & ${need}$ \\\\[4pt]")
    T["TabDerivBody"] = "\n".join(lines)
    # EoS constants
    ec = rcsv(root, "verification/eos_constants.csv")
    alpha = {"vdW": "1", "RK": "T_r^{-1/2}", "SRK": r"[1+m(1-T_r^{1/2})]^2",
             "PR": r"[1+m(1-T_r^{1/2})]^2"}
    eps = {"vdW": "0", "RK": "0", "SRK": "0", "PR": r"1-\sqrt{2}"}
    sig = {"vdW": "0", "RK": "1", "SRK": "1", "PR": r"1+\sqrt{2}"}
    T["TabEosBody"] = "\n".join(
        f"{r['eos']} & ${eps[r['eos']]}$ & ${sig[r['eos']]}$ & {float(r['Omega_a']):.5f} & "
        f"{float(r['Omega_b']):.5f} & ${alpha[r['eos']]}$ \\\\" for r in ec)
    # fluids
    fl = rcsv(root, "inputs/fluids.csv")
    T["TabFluidBody"] = "\n".join(
        f"{r['name']} & {1e3 * float(r['M_kg_per_mol']):.4f} & {float(r['Tc_K']):.4f} & "
        f"{1e-6 * float(r['Pc_Pa']):.4f} & {float(r['omega']):.5f} & {float(r['cpA']):.5f} & "
        f"{float(r['cpB']):.5f} & {float(r['cpC']):.5f} & {float(r['cpD']):.5f} & "
        f"{float(r['cpE']):.6f} & {float(r['cp_Tmin_K']):.0f}--{float(r['cp_Tmax_K']):.0f} \\\\"
        for r in fl)

    # validation summary
    def cell(ds, model, region, prop, scale=1.0, fmt="{:.1f}"):
        r = [x for x in S if x["dataset"] == ds and x["model"].startswith(model)
             and x["region"] == region and x["property"] == prop]
        return fmt.format(scale * float(r[0]["AAD"])) if r else "--"

    lines = []
    C = "CO2 310-400 K"
    for region, label in (("0.1-4 MPa", "0.1--4"), ("6-10 MPa", "6--10"), ("15-20 MPa", "15--20"),
                          ("all", "all")):
        for k, m in enumerate(("IG", "vdW", "RK", "SRK", "PR")):
            first = f"\\multirow{{5}}{{*}}{{{label}}}" if k == 0 else ""
            vals = [cell(C, m, region, p) for p in ("rho", "Cv", "Cp", "w", "muJT")]
            vals += [cell(C, m, region, "dH", 1e3, "{:.0f}"), cell(C, m, region, "dS", 1.0, "{:.2f}")]
            lines.append(f"{first} & {m} & " + " & ".join(vals) + " \\\\")
        lines.append("\\midrule" if region != "all" else "")
    T["TabValBody"] = "\n".join(lines).rstrip()
    lines = []
    for m in ("IG", "vdW", "RK", "SRK", "PR"):
        vals = [cell("N2 300 K", m, "1-31 MPa", p) for p in ("rho", "Cv", "Cp", "w", "muJT")]
        vals += [cell("N2 300 K", m, "1-31 MPa", "dH", 1e3, "{:.0f}"),
                 cell("N2 300 K", m, "1-31 MPa", "dS", 1.0, "{:.2f}")]
        lines.append(f"{'1--31' if m == 'IG' else ''} & {m} & " + " & ".join(vals) + " \\\\")
    T["TabValNtwoBody"] = "\n".join(lines)
    lines = []
    for m in ("vdW", "RK", "SRK", "PR"):
        vals = [cell("CO2 saturation 220-300 K", m, "all", p) for p in ("Psat", "rhoL", "rhoV", "dHvap")]
        lines.append(f"{m} & " + " & ".join(vals) + " \\\\")
    T["TabSatBody"] = "\n".join(lines)
    lines = []
    for model, label in (("table, all", "all 19 temperatures, evaluated at the nodes"),
                         ("table, 10 K", "10 K spacing, evaluated at held-out nodes")):
        lines.append(f"\\multicolumn{{7}}{{@{{}}l}}{{\\emph{{{label}}}}} \\\\")
        for region, rl in (("0.1-4 MPa", "0.1--4"), ("6-10 MPa", "6--10"), ("15-20 MPa", "15--20")):
            vals = [cell("CO2 tabulated route", model, region, "rho", 1.0, "{:.3f}")]
            vals += [cell("CO2 tabulated route", model, region, p, 1.0, "{:.0f}") for p in ("Cv", "Cp", "w")]
            vals += [cell("CO2 tabulated route", model, region, "dH", 1e3, "{:.0f}"),
                     cell("CO2 tabulated route", model, region, "dS", 1.0, "{:.2f}")]
            lines.append(f"{rl} & " + " & ".join(vals) + " \\\\")
    T["TabRouteBody"] = "\n".join(lines)
    # scenarios
    what = {
        "01": "state, 3 derivatives, uncertainty", "02": "process, uncertainty",
        "03": "saturation at 280 K", "04": "12 derivatives, symbolic only",
        "05": "inference from stated $V$, $C_P$, $w$, $\\mu_{JT}$", "06": "$\\alpha$ not supplied",
        "07": "stated $C_V > C_P$", "08": "process from tabulated $V(T,P)$",
        "09": "450 K requested, table ends at 400 K", "10": "$P$ given in kelvin",
        "11": "$T$ given without a unit", "12": "user fluid without $C_P^{\\circ}$",
        "13": "250 K, below $C_P^{\\circ}$ range", "14": "N$_2$, SRK, $^\\circ$C and bar input",
        "15": "ideal gas, tabulated $C_P^{\\circ}$", "16": "fluid not in database",
        "17": "$(\\partial P/\\partial T)_H$ in an ideal gas"}
    resp = {"ok": "completed", "partial": "partial result; missing input named",
            "inconsistent": "flagged inconsistent"}
    lines = []
    for r in sc:
        num = r["scenario"][:2]
        st = r["status"]
        outcome = resp.get(st, "refused: " + st.replace("error: ", "").replace("Error", " error")
                           .replace("OutOfRange", "out-of-range").replace("MissingData", "missing-data")
                           .replace("Unit", "unit"))
        if num == "17":
            outcome = "partial result; singular derivative reported"
        lines.append(f"{num} & {r['type']} & {what[num]} & {outcome} \\\\")
    T["TabScenBody"] = "\n".join(lines)
    # uncertainty budget
    bud = rcsv(root, "outputs/02_process_co2_pr_uncertainty_budget.csv")
    unit = {"dH": "J mol$^{-1}$", "dS": "J mol$^{-1}$ K$^{-1}$", "dU": "J mol$^{-1}$",
            "dG": "J mol$^{-1}$"}
    lab = {"dH": "$\\Delta h$", "dS": "$\\Delta s$", "dU": "$\\Delta u$", "dG": "$\\Delta g$"}
    lines = []
    for r in bud:
        if r["quantity"] not in unit:
            continue
        v, u = float(r["value"]), float(r["std_uncertainty"])
        prec = 0 if abs(v) > 500 else 2
        shares = " & ".join(f"{100 * float(r['share_' + k]):.0f}"
                            for k in ("T1", "P1", "T2", "P2", "Tc", "Pc", "omega"))
        lines.append(f"{lab[r['quantity']]} ({unit[r['quantity']]}) & {v:.{prec}f} & {u:.{prec}f} & {shares} \\\\")
    T["TabBudgetBody"] = "\n".join(lines)
    with open(os.path.join(root, "generated_tables.tex"), "w") as fh:
        fh.write("% GENERATED by Codes/run_all.py from the CSV files in Data/ -- do not edit.\n")
        for k, v in T.items():
            fh.write(f"\\newcommand{{\\{k}}}{{%\n{v}\n}}\n")
    print(f"  wrote generated_tables.tex ({len(T)} table bodies)")


# --------------------------------------------------------------------------
def step_website_data(root):
    """Website data: the same CSV content, as one JavaScript file."""
    print("[6] website data")
    sys.path.insert(0, os.path.join(root, "Codes"))
    from thermocalculus import symbolic
    from thermocalculus.eos import get_fluid
    der = {}
    for x, y, z in symbolic.all_derivatives():
        e = symbolic.derivative(x, y, z)
        pack = lambda p: [[float(c), list(k)] for k, c in p.t.items()]
        der[f"{x}{y}{z}"] = {"s": str(e), "n": pack(e.num), "d": pack(e.den),
                             "need": sorted(e.symbols())}
    fluids = {}
    for key in ("CO2", "N2"):
        f = get_fluid(key)
        fluids[key] = {"name": f.name, "M": f.M, "Tc": f.Tc, "Pc": f.Pc, "omega": f.omega,
                       "cp": list(f.cp0.c), "cpT": [f.cp0.tmin, f.cp0.tmax]}
    nist = [[float(r["T [K]"]), float(r["P [MPa]"]), float(r["rho [mol/L]"]),
             float(r["Cv [J/(mol*K)]"]), float(r["Cp [J/(mol*K)]"]), float(r["w [m/s]"]),
             float(r["muJT [K/MPa]"]), float(r["H [kJ/mol]"]), float(r["S [J/(mol*K)]"])]
            for r in rcsv(root, "reference/nist_co2_isobars_310-400K_0.1-20MPa.csv")]
    sat = [[float(r["T [K]"]), float(r["P [MPa]"]), float(r["rho [mol/L]"]), r["phase"]]
           for r in rcsv(root, "reference/nist_co2_saturation_220-300K.csv")]
    summ = [[r["dataset"], r["model"], r["region"], r["property"], float(r["AAD"]),
             float(r["max_abs_deviation"]), r["unit"]]
            for r in rcsv(root, "validation/summary_deviations.csv")]
    scen = []
    for r in rcsv(root, "outputs/scenario_summary.csv"):
        with open(os.path.join(root, "Data", "inputs", r["scenario"] + ".txt")) as fh:
            text = fh.read()
        with open(os.path.join(root, "Data", "outputs", r["scenario"] + "_log.txt")) as fh:
            log = [l for l in fh.read().splitlines() if l and not l.startswith("#")]
        rows = [[x["quantity"], x["value"], x["unit"], x["provenance"], x["std_uncertainty"]]
                for x in rcsv(root, f"outputs/{r['scenario']}_output.csv")][:14]
        scen.append({"name": r["scenario"], "type": r["type"], "status": r["status"],
                     "input": text, "log": log[:12], "rows": rows})
    budget = {}
    for nm in ("02_process_co2_pr", "05_infer_from_measurements"):
        rows = rcsv(root, f"outputs/{nm}_uncertainty_budget.csv")
        keys = [k[6:] for k in rows[0] if k.startswith("share_")]
        budget[nm] = {"inputs": keys, "rows": [
            [r["quantity"], float(r["value"]), float(r["std_uncertainty"])]
            + [float(r["share_" + k]) for k in keys] for r in rows]}
    nums = {r["macro"]: r["value_as_typeset"] for r in rcsv(root, "manuscript_numbers.csv")}
    chk = [float(r["relative_deviation"]) for r in rcsv(root, "verification/derivatives_336_check.csv")
           if r["status"] == "checked" and r["state"] == "PR supercritical"]
    # reference values for the JavaScript port of the equations of state
    from thermocalculus import props
    from thermocalculus.eos import EoS
    js_ref = []
    for key in ("CO2", "N2"):
        f = get_fluid(key)
        for kind in ("IG", "vdW", "RK", "SRK", "PR"):
            for Tt, Pp in ((320.0, 2e6), (350.0, 10e6), (400.0, 20e6)):
                st = props.state(f, EoS(kind, f.Tc, f.Pc, f.omega), Tt, Pp)
                js_ref.append([key, kind, Tt, Pp] + [st[k] for k in
                              ("V", "Z", "alpha", "kappaT", "HR", "SR", "Cp", "Cv", "w", "muJT", "H", "S")])
    data = {"der": der, "fluids": fluids, "nist": nist, "sat": sat, "summary": summ,
            "scenarios": scen, "budget": budget, "numbers": nums,
            "derivCheck": sorted(chk), "jsRef": js_ref}
    path = os.path.join(root, "Codes", "website", "data.js")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("// GENERATED by Codes/run_all.py from the CSV files in Data/ -- do not edit.\n")
        fh.write("window.TC_DATA = " + json.dumps(data, separators=(",", ":")) + ";\n")
    print(f"  wrote Codes/website/data.js ({os.path.getsize(path) // 1024} kB)")
    engine = os.path.join(root, "Codes", "website", "engine.js")
    check = os.path.join(root, "Codes", "website", "check_engine.js")
    if os.path.exists(engine) and os.path.exists(check):
        try:
            res = subprocess.run(["node", check], capture_output=True, text=True, timeout=120)
            with open(os.path.join(root, "Data", "verification", "website_js_vs_python.csv"), "w") as fh:
                fh.write(res.stdout)
            print("  wrote Data/verification/website_js_vs_python.csv:",
                  res.stdout.strip().splitlines()[1][2:] if res.stdout else res.stderr[:200])
        except (OSError, subprocess.SubprocessError) as exc:
            print("  node not available; website cross-check skipped:", exc)


# --------------------------------------------------------------------------
DESCRIPTIONS = [
    (r"inputs/\d\d_.*\.txt", "INPUT scenario file (plain-text INI); see Codes/README.md for the format"),
    (r"inputs/sessions/ask_session_script\.txt", "INPUT scripted answers for the interactive session (python -m thermocalculus ask --script ...)"),
    (r"outputs/ask_session_transcript\.txt", "OUTPUT transcript of the interactive session: the engine asks for missing data, rejects a wrong unit, then computes the targets"),
    (r"inputs/fluids\.csv", "INPUT fluid database used by the engine (copy of Codes/thermocalculus/data/fluids.csv): M, Tc, Pc, omega, Shomate Cp0 coefficients, sources"),
    (r"inputs/tables/co2_pvt_table_nist\.csv", "INPUT table for the tabulated-input scenarios: T [K], P [MPa], V [L/mol] on a 19 x 9 grid (NIST values)"),
    (r"inputs/tables/co2_cp0_table\.csv", "INPUT table of ideal-gas Cp0(T) [J/(mol K)], 300-500 K"),
    (r"outputs/scenario_summary\.csv", "OUTPUT one line per scenario: type, status, first diagnostic message"),
    (r"outputs/.*_output\.csv", "OUTPUT of the scenario of the same name: quantity, value, unit, provenance, standard uncertainty (linear and Monte Carlo), note"),
    (r"outputs/.*_log\.txt", "OUTPUT log of the scenario of the same name: model used, warnings, missing-data / inconsistency / error messages"),
    (r"outputs/.*_uncertainty_budget\.csv", "OUTPUT variance share of every uncertain input for every output of that scenario"),
    (r"reference/raw_nist/.*", "REFERENCE NIST WebBook table as retrieved (text; columns T K, P MPa, rho mol/L, V L/mol, U kJ/mol, H kJ/mol, S J/(mol K), Cv, Cp J/(mol K), w m/s, JT K/MPa, phase)"),
    (r"reference/nist_.*\.csv", "REFERENCE NIST WebBook values with labelled columns, built from raw_nist/ after integrity checks"),
    (r"verification/all_derivatives_reduced\.csv", "OUTPUT all 336 reduced derivatives: expression, LaTeX, symbols needed"),
    (r"verification/differential_table\.csv", "OUTPUT (dX/dT)_P and (dX/dP)_T of the eight variables"),
    (r"verification/identity_checks\.csv", "OUTPUT symbolic checks: Maxwell relations, reciprocity, cyclic rule"),
    (r"verification/derivatives_336_check\.csv", "OUTPUT each reduced formula against a constrained numerical derivative, two CO2 states"),
    (r"verification/closedform_vs_numerical\.csv", "OUTPUT closed-form cubic EoS properties against the generic numerical scheme"),
    (r"verification/consistency_checks\.csv", "OUTPUT numerical derivatives of G and H against S, V, Cp"),
    (r"verification/eos_constants\.csv", "OUTPUT constants of each cubic EoS and critical-point residuals"),
    (r"verification/test_suite_result\.txt", "OUTPUT result of the automated test suite"),
    (r"verification/website_js_vs_python\.csv", "OUTPUT cross-check of the website's JavaScript equations of state against the Python engine"),
    (r"validation/co2_eos_vs_nist\.csv", "OUTPUT five models against NIST at 171 CO2 state points: rho, Cv, Cp, w, muJT, dH, dS"),
    (r"validation/n2_eos_vs_nist\.csv", "OUTPUT five models against NIST for N2 at 300 K"),
    (r"validation/co2_saturation_vs_nist\.csv", "OUTPUT saturation pressure, coexisting densities, enthalpy of vaporisation against NIST"),
    (r"validation/co2_table_route_vs_nist\.csv", "OUTPUT tabulated-input route against NIST (full table and held-out test)"),
    (r"validation/co2_model_curves\.csv", "OUTPUT dense model curves drawn in Figures 3 and 5"),
    (r"validation/co2_saturation_curves\.csv", "OUTPUT dense model saturation curves drawn in Figure 6"),
    (r"validation/summary_deviations\.csv", "OUTPUT average and maximum absolute deviations by data set, model, region, property"),
    (r"manuscript_numbers\.csv", "OUTPUT every number quoted in the manuscript with the Data file it comes from"),
]


def step_index(root):
    print("[8] data index")
    data = os.path.join(root, "Data")
    files = sorted(os.path.relpath(os.path.join(d, f), data).replace(os.sep, "/")
                   for d, _, fs in os.walk(data) for f in fs if f != "DATA_INDEX.txt")
    lines = ["DATA INDEX  (generated by Codes/run_all.py)",
             "=" * 60,
             "All files are plain text (.txt) or comma-separated (.csv).  Lines starting with '#'",
             "are comments that state content, source and units.  Values are in SI units unless a",
             "column header or the 'unit' column says otherwise.",
             "",
             "Folders:  inputs/ what the engine reads | outputs/ what the scenarios produce |",
             "          reference/ NIST data used for validation | verification/ checks of the",
             "          mathematics | validation/ comparison with reference data",
             ""]
    unknown = []
    for f in files:
        desc = next((d for pat, d in DESCRIPTIONS if re.fullmatch(pat, f)), None)
        if desc is None:
            unknown.append(f)
            desc = "(no description)"
        lines.append(f"{f}\n    {desc}")
    with open(os.path.join(data, "DATA_INDEX.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  wrote Data/DATA_INDEX.txt ({len(files)} files, {len(unknown)} without description)")
