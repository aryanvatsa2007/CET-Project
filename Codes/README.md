# thermocalculus: a computational engine for the calculus of thermodynamics

CL 202 project, Topic 2. This folder holds the complete code package.

## What it does

| Capability | Module |
|---|---|
| Reduces any first derivative (dX/dY)_Z of T, P, V, S, U, H, A, G (336 cases) exactly to T, P, V, S, Cp, alpha, kappaT | `thermocalculus/symbolic.py` |
| Ideal-gas, van der Waals, Redlich-Kwong, Soave-Redlich-Kwong and Peng-Robinson equations of state with closed-form residual properties, phase-root selection and saturation | `thermocalculus/eos.py`, `props.py` |
| Model-independent numerical differentiation and integration for any P(T, V) | `thermocalculus/numerics.py` |
| Tabulated P-V-T and Cp0(T) input (splines, no extrapolation) | `thermocalculus/tables.py` |
| Computability: infers what follows from the data given, lists the inputs that are missing, flags inconsistent inputs | `thermocalculus/knowledge.py` |
| Propagation of input uncertainty (first order, with Monte Carlo check) | `thermocalculus/uncertainty.py` |
| Interactive mode: the engine prompts the user for missing inputs until each target is computable (or the user declines) | `thermocalculus/interactive.py` |
| Units: every input is `<number> <unit>`; wrong or missing units are rejected | `thermocalculus/units.py` |
| Scenario files in, labelled CSV and log files out | `thermocalculus/scenario.py` |

## Folder structure

```
Codes/
  README.md                 this file
  run_all.py                regenerates every Data file, figure and manuscript number
  thermocalculus/           the engine (Python package)
    data/fluids.csv         fluid database (CO2, N2) with sources
  tests/test_engine.py      automated tests
  (thermocalculus/interactive.py  the prompting dialogue behind `python -m thermocalculus ask`)
  tools/
    build_reference_data.py NIST tables -> labelled CSV, with integrity checks
    figures.py              manuscript figures (drawn only from Data/ files)
    products.py             manuscript numbers and tables, website data, Data index
    audit.py                consistency audit -> ../AUDIT_REPORT.txt
  website/                  interactive web version (static; see below)
```

## Installation

Python 3.10 or newer. The engine itself needs only NumPy and SciPy; Matplotlib
is needed for the figures and pytest for the tests.

```
pip install numpy scipy matplotlib pytest
```

Tested with Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2.
No installation of the package is needed: run the commands below from `Codes/`.

## How to run, step by step

All commands are run from the `Codes` folder.

1. Reduce a derivative symbolically:

   ```
   python -m thermocalculus derive "(dT/dP)_H"
   ```
   Expected output:
   ```
   (dT/dP)_H = V*(T*alpha - 1)/Cp
   needs: Cp, T, V, alpha
   LaTeX: \frac{V\left(T \alpha - 1\right)}{C_P}
   ```

2. Verify the four Maxwell relations: `python -m thermocalculus maxwell`

3. List the fluid database: `python -m thermocalculus fluids`

4. Run one scenario (input files are in `../Data/inputs`):

   ```
   python -m thermocalculus run ../Data/inputs/01_state_co2_pr.txt -o my_outputs
   ```
   This writes `my_outputs/01_state_co2_pr_output.csv` (quantity, value, unit,
   provenance, standard uncertainty, note), `..._log.txt` and, when
   uncertainties are given, `..._uncertainty_budget.csv`. The same files are in
   `../Data/outputs` for comparison. Several files can be given at once:
   `python -m thermocalculus run ../Data/inputs/*.txt -o my_outputs`

5. Interactive session: the engine asks for the data it is missing.

   ```
   python -m thermocalculus ask                       # type at the keyboard; 'help' lists commands
   python -m thermocalculus ask --fluid CO2 --eos PR  # T and P alone are then enough for an EoS state
   python -m thermocalculus ask --script ../Data/inputs/sessions/ask_session_script.txt
   ```
   Example (typed lines follow the `>` prompt): give the targets with
   `target (dT/dP)_H; Cv; w`, then values such as `T = 350 K`. For every target
   that cannot be computed the engine prints the smallest sets of extra
   quantities that would suffice and asks for one of them
   (`Enter a value (name = value unit), or 'skip' ...`). Wrong units are
   rejected and asked again. `show` lists what is known and inferred; `quit`
   ends the session and prints the final result. The recorded scripted session
   is `../Data/outputs/ask_session_transcript.txt`.

6. Regenerate everything used in the manuscript (about one minute):

   ```
   python run_all.py
   ```
   This rebuilds `../Data/` (outputs, verification, validation, index), the
   figure PDFs, `../generated_numbers.tex`, `../generated_tables.tex` and
   `website/data.js`. Then compile the manuscript from the project folder:
   `pdflatex manuscript; bibtex manuscript; pdflatex manuscript; pdflatex manuscript`

7. Run the tests: `python -m pytest -q` (expected: all tests pass).

8. Run the consistency audit: `python tools/audit.py` (writes `../AUDIT_REPORT.txt`).

## Scenario file format

A scenario is a plain-text file with sections. Every quantity carries a unit.

```
[scenario]
name = 01_state_co2_pr
type = state            ; state | process | saturation | derivative | infer | table
description = free text

[model]
fluid = CO2             ; or give Tc, Pc, omega, M yourself
eos = PR                ; IG | vdW | RK | SRK | PR
; optional: cp0_constant = 37 J/(mol*K) | cp0_table = tables/file.csv
;           cp0_shomate = A B C D E Tmin Tmax | extrapolate_cp0 = yes

[state]                 ; process and table scenarios use [state1] and [state2]
T = 350 K
P = 10 MPa
phase = stable          ; stable | liquid | vapor

[given]                 ; infer and derivative scenarios: stated values, no EoS
V = 0.19235 L/mol

[table]                 ; table scenarios
pvt = tables/co2_pvt_table_nist.csv

[targets]
derivatives = (dT/dP)_H, (dU/dV)_T
quantities = Cv, w

[uncertainty]           ; standard uncertainties of inputs or model parameters
T = 0.1 K
monte_carlo = 2000      ; optional check

[output_units]
P = MPa
```

Units understood: K, degC, degF, degR; Pa, kPa, MPa, GPa, bar, mbar, atm, psi,
torr; m3, L, mL, cm3, dm3; J, kJ, MJ, cal, kcal; mol, kmol; kg, g; m, cm, mm; s;
combined with `*`, `/`, `^` and parentheses, for example `J/(mol*K)`.

Tables are CSV files whose header cells are written `name [unit]`, for example
`T [K],P [MPa],V [L/mol]` (a full rectangular grid) or `T [K],Cp0 [J/(mol*K)]`.

## Scenario outcomes

Each run ends with one status: `ok`; `partial` (some targets could not be
computed, the log names the missing input); `inconsistent` (stated values
contradict each other); or `error: <type>` (the request was refused: unit
error, out-of-range table or correlation, missing data, malformed input). The
17 scenarios in `../Data/inputs` cover all of these.

## Website

`website/` is a static site (HTML, CSS, JavaScript, no build step). Open
`website/index.html` in a browser, or publish the folder with GitHub Pages
(repository Settings -> Pages -> deploy from the folder that contains
`index.html`). The page loads two web fonts from Google Fonts and nothing else.

Link to the hosted copy: [ADD THE GROUP'S GITHUB PAGES LINK HERE AFTER PUBLISHING]

The equations of state are re-implemented in `website/engine.js`;
`node website/check_engine.js` compares them with the Python engine (result in
`../Data/verification/website_js_vs_python.csv`). The reduced derivatives and
all stored results shown on the site are read from `website/data.js`, which
`run_all.py` writes from the Data files.

## Limitations

One pure fluid, one phase per state; no mixtures, reactions or solids; first
derivatives only. Only CO2 and N2 are in the database. The CO2 ideal-gas Cp
correlation starts at 298 K. Cubic equations are inaccurate near the critical
point and for liquid densities. Second-derivative properties (Cp, Cv, speed of
sound) from a tabulated P-V-T grid are not reliable. Uncertainty propagation
assumes uncorrelated inputs and does not include the error of the equation of
state. See the manuscript, section "Limitations and Failure Modes".
