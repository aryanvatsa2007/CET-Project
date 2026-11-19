/* Cross-check of the website's JavaScript equations of state against reference
   values written by the Python engine (data.js -> jsRef).  Run: node check_engine.js
   Prints a CSV; run_all.py stores it as Data/verification/website_js_vs_python.csv */
const path = require("path");
global.window = {};
require(path.join(__dirname, "data.js"));
const E = require(path.join(__dirname, "engine.js"));
const D = global.window.TC_DATA;
const names = ["V", "Z", "alpha", "kappaT", "HR", "SR", "Cp", "Cv", "w", "muJT", "H", "S"];
let rows = [], worst = 0;
for (const r of D.jsRef) {
  const fl = D.fluids[r[0]], st = E.state(fl, E.makeEos(r[1], fl), r[2], r[3]);
  names.forEach((n, i) => {
    const py = r[4 + i], js = st[n];
    const diff = py === 0 ? Math.abs(js) : Math.abs(js - py) / Math.abs(py);
    worst = Math.max(worst, diff);
    rows.push([r[0], r[1], r[2], r[3], n, py, js, diff.toExponential(3)].join(","));
  });
}
// derivative records evaluated in JS at one state
const fl = D.fluids.CO2, st = E.state(fl, E.makeEos("PR", fl), 350, 10e6);
const mu = E.evalDerivative(D.der.TPH, st);
worst = Math.max(worst, Math.abs(mu - st.muJT) / Math.abs(st.muJT));
console.log("# website JavaScript engine against the Python engine (SI units)");
console.log(`# max relative difference: ${worst.toExponential(2)} over ${rows.length} comparisons`);
console.log("fluid,eos,T_K,P_Pa,property,python,javascript,relative_difference");
console.log(rows.join("\n"));
