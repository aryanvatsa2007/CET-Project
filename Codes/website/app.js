/* Thermocalculus website: interface and the WebGL scene.
   All physics comes from engine.js; all stored results come from data.js. */
(function () {
  "use strict";
  var D = window.TC_DATA, E = window.TCEngine, CO2 = D.fluids.CO2;
  var S = { eos: "PR", T: 350, P: 10, x: "T", y: "P", z: "H", vprop: "rho", veos: "PR", scen: 0 };
  var VARS = ["T", "P", "V", "S", "U", "H", "A", "G"];
  var MODELS = ["IG", "vdW", "RK", "SRK", "PR"];
  var MODEL_NAME = { IG: "Ideal gas", vdW: "van der Waals", RK: "Redlich–Kwong", SRK: "Soave–RK", PR: "Peng–Robinson" };
  var MODEL_COLOR = { PR: "#5aa2ff", SRK: "#ff8a5c", RK: "#39d6a0", vdW: "#f5c043", IG: "#f29bc0" };
  var $ = function (id) { return document.getElementById(id); };
  var eosCache = {};
  function eos(kind) { return eosCache[kind] || (eosCache[kind] = E.makeEos(kind, CO2)); }
  function el(tag, attrs, html) {
    var n = document.createElement(tag);
    for (var k in attrs || {}) n.setAttribute(k, attrs[k]);
    if (html !== undefined) n.innerHTML = html;
    return n;
  }
  function sig(x, n) {
    if (x === null || !isFinite(x)) return "–";
    var a = Math.abs(x);
    if (a !== 0 && (a < 1e-3 || a >= 1e6)) {
      var p = x.toExponential((n || 4) - 1).split("e");
      return p[0] + " × 10<sup>" + (+p[1]) + "</sup>";
    }
    return (+x.toPrecision(n || 4)).toString().replace("-", "−");
  }
  function texNum(s) { return s.replace(/ \\times 10\^\{(-?\d+)\}/, " × 10<sup>$1</sup>").replace("-", "−"); }

  /* ------------------------------------------------------------ hero facts */
  var N = D.numbers;
  [[N.NDerivTotal, "derivatives reduced exactly"],
   [texNum(N.DerivMaxDevSC), "largest gap between a formula and a direct numerical derivative"],
   [N.NNistTotal, "NIST reference states used for validation"],
   [N.NTests, "automated tests, all passing"]].forEach(function (f) {
    $("facts").appendChild(el("div", {}, "<strong>" + f[0] + "</strong><span>" + f[1] + "</span>"));
  });
  $("checked").innerHTML = "All " + N.NDerivTotal + " formulas were compared with direct numerical differentiation at two states (largest relative gap " +
    texNum(N.DerivMaxDevSC) + "). Closed-form results match an independent numerical scheme to " + texNum(N.ClosedNumMax) +
    ". Against NIST data for carbon dioxide between 310 and 400 K and up to 20 MPa, Peng–Robinson densities are off by " + N.PRrhoAAD +
    "% on average and by up to " + N.PRrhoMax + "% near the critical point.";

  /* ------------------------------------------------------------ formulas */
  function symHTML(s) {
    return s.replace(/kappaT/g, "κ<sub>T</sub>").replace(/alpha/g, "α").replace(/Cp/g, "C<sub>P</sub>")
      .replace(/\^(-?\d+)/g, "<sup>$1</sup>").replace(/\*/g, " ").replace(/ - /g, " − ").replace(/^-/, "−");
  }
  function splitTop(s) {
    var depth = 0;
    for (var i = 0; i < s.length; i++) {
      if (s[i] === "(") depth++; else if (s[i] === ")") depth--;
      else if (s[i] === "/" && depth === 0) return i;
    }
    return -1;
  }
  function unwrap(s) {
    if (s[0] !== "(") return s;
    var depth = 0;
    for (var i = 0; i < s.length; i++) {
      if (s[i] === "(") depth++; else if (s[i] === ")") { depth--; if (depth === 0) return i === s.length - 1 ? s.slice(1, -1) : s; }
    }
    return s;
  }
  function formulaHTML(s) {
    var i = splitTop(s);
    if (i < 0) return "<span class='m'>" + symHTML(s) + "</span>";
    var num = s.slice(0, i), den = s.slice(i + 1), neg = "";
    if (num[0] === "-" && splitTop(num.slice(1)) < 0 && num.indexOf(" ") < 0) { neg = "−"; num = num.slice(1); }
    return neg + "<span class='frac m'><span>" + symHTML(unwrap(num)) + "</span><span>" + symHTML(unwrap(den)) + "</span></span>";
  }
  var UNIT = { T: "K", P: "Pa", V: "m³/mol", S: "J/(mol·K)", U: "J/mol", H: "J/mol", A: "J/mol", G: "J/mol" };

  function chips(id, key) {
    var box = $(id);
    VARS.forEach(function (v) {
      var b = el("button", { type: "button", "class": "chip", "aria-pressed": "false", "aria-label": v }, v);
      b.addEventListener("click", function () { S[key] = v; fixChoice(key); renderDerivative(); });
      box.appendChild(b);
    });
  }
  function fixChoice(changed) {
    // the three choices must differ; move the others off a clash
    ["x", "y", "z"].forEach(function (k) {
      if (k === changed) return;
      var others = ["x", "y", "z"].filter(function (o) { return o !== k; }).map(function (o) { return S[o]; });
      if (others.indexOf(S[k]) >= 0) S[k] = VARS.filter(function (v) { return others.indexOf(v) < 0; })[0];
    });
  }
  function renderDerivative() {
    [["cx", "x"], ["cy", "y"], ["cz", "z"]].forEach(function (p) {
      Array.prototype.forEach.call($(p[0]).children, function (b) {
        b.setAttribute("aria-pressed", b.textContent === S[p[1]] ? "true" : "false");
      });
    });
    var rec = D.der[S.x + S.y + S.z];
    $("eq").innerHTML = "<span class='lhs'><span class='par'>(</span><span class='frac m'><span>∂" + S.x + "</span><span>∂" + S.y +
      "</span></span><span class='par'>)</span><span class='hold m'>" + S.z + "</span></span><span>=</span>" + formulaHTML(rec.s);
    $("needs").innerHTML = rec.need.length ? "<span>Needs</span>" + rec.need.map(function (n) {
      return "<span class='need'>" + symHTML(n) + "</span>"; }).join("") +
      (rec.need.indexOf("S") < 0 ? "<span>no entropy required</span>" : "") : "<span>Needs no data: it is a constant.</span>";
    var st = E.state(CO2, eos(S.eos), S.T, S.P * 1e6), val = E.evalDerivative(rec, st);
    $("dval").innerHTML = val === null
      ? "For carbon dioxide at " + S.T + " K and " + S.P + " MPa (" + MODEL_NAME[S.eos] + ") this derivative is singular: the property you vary cannot change while the third is held fixed."
      : "For carbon dioxide at " + S.T + " K and " + S.P + " MPa (" + MODEL_NAME[S.eos] + "): <b>" + sig(val, 5) + "</b> " + UNIT[S.x] + " per " + UNIT[S.y];
  }

  /* ------------------------------------------------------------ segmented controls */
  function seg(id, items, key, after) {
    var box = $(id);
    items.forEach(function (it) {
      var b = el("button", { type: "button", "aria-pressed": S[key] === it[0] ? "true" : "false" }, it[1]);
      b.addEventListener("click", function () {
        S[key] = it[0];
        Array.prototype.forEach.call(box.children, function (c) { c.setAttribute("aria-pressed", c === b ? "true" : "false"); });
        after();
      });
      box.appendChild(b);
    });
  }

  /* ------------------------------------------------------------ state explorer */
  var NIST = D.nist;    // [T, P, rho, Cv, Cp, w, muJT, H, S]
  function nistAt(T) { return NIST.filter(function (r) { return r[0] === T; }); }
  function renderState() {
    $("oT").textContent = S.T + " K";
    $("oP").textContent = S.P.toFixed(1) + " MPa";
    var st = E.state(CO2, eos(S.eos), S.T, S.P * 1e6);
    var near = nistAt(S.T).filter(function (r) { return Math.abs(r[1] - S.P) < 1e-9; })[0];
    var items = [
      ["Compressibility factor Z", sig(st.Z, 4), near ? sig(near[1] * 1e6 / (near[2] * 1e3 * E.R * S.T), 4) : null],
      ["Density, mol/L", sig(1e-3 / st.V, 4), near ? sig(near[2], 5) : null],
      ["Heat capacity C<sub>P</sub>, J/(mol·K)", sig(st.Cp, 4), near ? sig(near[4], 5) : null],
      ["Heat capacity C<sub>V</sub>, J/(mol·K)", sig(st.Cv, 4), near ? sig(near[3], 5) : null],
      ["Speed of sound, m/s", sig(st.w, 4), near ? sig(near[5], 5) : null],
      ["Joule–Thomson, K/MPa", sig(st.muJT * 1e6, 4), near ? sig(near[6], 5) : null]];
    $("read").innerHTML = items.map(function (i) {
      return "<div><span>" + i[0] + "</span><b>" + i[1] + "</b>" + (i[2] ? "<small>NIST " + i[2] + "</small>" : "<small>&nbsp;</small>") + "</div>";
    }).join("");
    $("readnote").textContent = near ? "NIST lists this exact state, so reference values are shown under each number."
      : "Move the pressure to 0.1, 1, 2, 4, 6, 8, 10, 15 or 20 MPa to compare with a NIST value at this exact state.";
    drawChart();
    renderDerivative();
    scene.setMarker(st);
  }
  function drawChart() {
    var W = 640, H = 400, m = { l: 52, r: 14, t: 10, b: 44 };
    var x = function (p) { return m.l + (W - m.l - m.r) * p / 20; };
    var y = function (z) { return m.t + (H - m.t - m.b) * (1.05 - z) / 0.95; };
    var svg = "<svg viewBox='0 0 " + W + " " + H + "' role='img' aria-label='Compressibility factor against pressure'>";
    [0.2, 0.4, 0.6, 0.8, 1.0].forEach(function (z) {
      svg += "<line x1='" + m.l + "' x2='" + (W - m.r) + "' y1='" + y(z) + "' y2='" + y(z) + "' stroke='rgba(160,170,230,.14)'/>" +
        "<text x='" + (m.l - 9) + "' y='" + (y(z) + 4.5) + "' text-anchor='end' fill='#a3a9cc' font-size='13'>" + z.toFixed(1) + "</text>";
    });
    [0, 5, 10, 15, 20].forEach(function (p) {
      svg += "<text x='" + x(p) + "' y='" + (H - m.b + 20) + "' text-anchor='middle' fill='#a3a9cc' font-size='13'>" + p + "</text>";
    });
    svg += "<text x='" + (m.l + (W - m.l - m.r) / 2) + "' y='" + (H - 4) + "' text-anchor='middle' fill='#a3a9cc' font-size='13.5'>pressure, MPa</text>";
    svg += "<text transform='translate(14 " + (m.t + (H - m.t - m.b) / 2) + ") rotate(-90)' text-anchor='middle' fill='#a3a9cc' font-size='13.5'>compressibility factor Z</text>";
    var curves = {};
    MODELS.forEach(function (k) {
      var pts = [], e = eos(k);
      for (var p = 0.1; p <= 20.0001; p += 0.1) pts.push([p, p * 1e6 * e.volume(S.T, p * 1e6) / (E.R * S.T)]);
      curves[k] = pts;
      var on = k === S.eos;
      svg += "<path d='M" + pts.map(function (q) { return x(q[0]).toFixed(1) + " " + y(q[1]).toFixed(1); }).join("L") +
        "' fill='none' stroke='" + MODEL_COLOR[k] + "' stroke-width='" + (on ? 3.2 : 1.4) + "' opacity='" + (on ? 1 : 0.5) +
        "' stroke-linejoin='round'" + (k === "IG" ? " stroke-dasharray='3 5'" : "") + "/>";
    });
    nistAt(S.T).forEach(function (r) {
      var z = r[1] * 1e6 / (r[2] * 1e3 * E.R * S.T);
      svg += "<circle cx='" + x(r[1]) + "' cy='" + y(z) + "' r='5' fill='#0b0f2a' stroke='#f4f0e6' stroke-width='2'/>";
    });
    var zc = S.P * 1e6 * eos(S.eos).volume(S.T, S.P * 1e6) / (E.R * S.T);
    svg += "<line id='cross' x1='" + x(S.P) + "' x2='" + x(S.P) + "' y1='" + m.t + "' y2='" + (H - m.b) + "' stroke='rgba(255,241,201,.45)' stroke-dasharray='4 4'/>" +
      "<circle cx='" + x(S.P) + "' cy='" + y(zc) + "' r='7' fill='#fff1c9' stroke='#0b0f2a' stroke-width='2'/>";
    svg += "<rect id='hit' x='" + m.l + "' y='" + m.t + "' width='" + (W - m.l - m.r) + "' height='" + (H - m.t - m.b) + "' fill='transparent'/></svg>";
    var box = $("chart"), tip = $("ctip");
    var old = box.querySelector("svg");
    if (old) old.remove();
    box.insertAdjacentHTML("afterbegin", svg);
    $("ctitle").textContent = "Compressibility factor of carbon dioxide at " + S.T + " K";
    $("clegend").innerHTML = MODELS.slice().reverse().map(function (k) {
      return "<span><i style='border-color:" + MODEL_COLOR[k] + (k === "IG" ? ";border-top-style:dotted" : "") + "'></i>" + MODEL_NAME[k] + "</span>";
    }).join("") + "<span><i class='dot'></i>NIST reference</span>";
    var hit = box.querySelector("#hit");
    hit.addEventListener("mousemove", function (ev) {
      var r = hit.getBoundingClientRect(), p = Math.max(0.1, Math.min(20, 20 * (ev.clientX - r.left) / r.width));
      var idx = Math.round((p - 0.1) / 0.1), rows = "<b>" + curves.PR[idx][0].toFixed(1) + " MPa</b>";
      MODELS.slice().reverse().forEach(function (k) {
        rows += "<br><span style='color:" + (k === S.eos ? "#14122a" : "#5d5a78") + "'>" + MODEL_NAME[k] + " " + curves[k][idx][1].toFixed(3) + "</span>";
      });
      tip.innerHTML = rows;
      var br = box.getBoundingClientRect();
      tip.style.left = Math.min(ev.clientX - br.left + 14, br.width - 190) + "px";
      tip.style.top = (ev.clientY - br.top + 12) + "px";
      tip.style.opacity = 1;
    });
    hit.addEventListener("mouseleave", function () { tip.style.opacity = 0; });
  }

  /* ------------------------------------------------------------ validation grid */
  var VP = { rho: ["Density", 2, function (st) { return 1e-3 / st.V; }], Cp: ["C<sub>P</sub>", 4, function (st) { return st.Cp; }],
    Cv: ["C<sub>V</sub>", 3, function (st) { return st.Cv; }], w: ["Sound speed", 5, function (st) { return st.w; }],
    muJT: ["Joule–Thomson", 6, function (st) { return st.muJT * 1e6; }] };
  var Ts = [], Ps = [];
  NIST.forEach(function (r) { if (Ts.indexOf(r[0]) < 0) Ts.push(r[0]); if (Ps.indexOf(r[1]) < 0) Ps.push(r[1]); });
  Ts.sort(function (a, b) { return a - b; }); Ps.sort(function (a, b) { return b - a; });
  function mix(a, b, t) { return a.map(function (v, i) { return Math.round(v + (b[i] - v) * t); }); }
  function divColor(v, lim) {
    var t = Math.max(-1, Math.min(1, v / lim)), mid = [38, 43, 71];
    var c = t < 0 ? mix(mid, [61, 108, 246], -t) : mix(mid, [255, 138, 61], t);
    return "rgb(" + c.join(",") + ")";
  }
  function renderGrid() {
    var prop = VP[S.vprop], e = eos(S.veos), lim = S.vprop === "rho" ? 12 : 30;
    $("ylab").innerHTML = Ps.map(function (p) { return "<span>" + p + "</span>"; }).join("");
    $("xlab").innerHTML = Ts.map(function (t, i) { return "<span>" + (i % 2 ? "" : t) + "</span>"; }).join("");
    var cells = $("cells"), tip = $("vtip");
    cells.innerHTML = "";
    Ps.forEach(function (p) {
      Ts.forEach(function (t) {
        var ref = NIST.filter(function (r) { return r[0] === t && r[1] === p; })[0];
        var model = prop[2](E.state(CO2, e, t, p * 1e6)), dev = 100 * (model / ref[prop[1]] - 1);
        var c = el("div", { "class": "cell", tabindex: "0", "aria-label": t + " K, " + p + " MPa: " + dev.toFixed(1) + " percent" });
        c.style.background = divColor(dev, lim);
        var show = function () {
          tip.innerHTML = "<b>" + t + " K, " + p + " MPa</b><br>" + MODEL_NAME[S.veos] + " " + sig(model, 5) + "<br>NIST " + sig(ref[prop[1]], 5) +
            "<br><b>" + (dev > 0 ? "+" : "−") + Math.abs(dev).toFixed(1) + "%</b>";
          var br = cells.parentNode.parentNode.getBoundingClientRect(), cr = c.getBoundingClientRect();
          tip.style.left = Math.min(cr.left - br.left + cr.width + 8, br.width - 170) + "px";
          tip.style.top = (cr.top - br.top + cr.height + 6) + "px";
          tip.style.opacity = 1;
        };
        c.addEventListener("mouseenter", show); c.addEventListener("focus", show);
        c.addEventListener("mouseleave", function () { tip.style.opacity = 0; });
        c.addEventListener("blur", function () { tip.style.opacity = 0; });
        cells.appendChild(c);
      });
    });
    var h = cells.children[0].getBoundingClientRect().height;
    Array.prototype.forEach.call($("ylab").children, function (s) { s.style.height = h + "px"; });
    $("smin").textContent = "−" + lim + "% (model too low)";
    $("smax").textContent = "+" + lim + "% (model too high)";
    $("vsum").innerHTML = MODELS.map(function (k) {
      var row = D.summary.filter(function (r) { return r[0] === "CO2 310-400 K" && r[1] === k && r[2] === "all" && r[3] === S.vprop; })[0];
      return "<div class='" + (k === S.veos ? "on" : "") + "'><span>" + MODEL_NAME[k] + "</span><b>" + row[4].toFixed(1) + "%</b></div>";
    }).join("");
  }

  /* ------------------------------------------------------------ scenarios */
  function statusChip(st) {
    var cls = st === "ok" ? "ok" : (st.indexOf("error") === 0 ? "err" : "partial");
    var txt = st === "ok" ? "Completed" : st === "partial" ? "Partial" : st === "inconsistent" ? "Inconsistent" : "Refused";
    return "<span class='st " + cls + "'>" + txt + "</span>";
  }
  function esc(s) { return s.replace(/&/g, "&amp;").replace(/</g, "&lt;"); }
  function renderScenarios() {
    var list = $("slist");
    list.innerHTML = "";
    D.scenarios.forEach(function (sc, i) {
      var title = sc.name.slice(3).replace(/_/g, " ");
      var b = el("button", { type: "button", "aria-pressed": i === S.scen ? "true" : "false" },
        "<span>" + sc.name.slice(0, 2) + "&ensp;" + title + "</span>" + statusChip(sc.status));
      b.addEventListener("click", function () { S.scen = i; renderScenarios(); });
      list.appendChild(b);
    });
    var sc = D.scenarios[S.scen];
    var flagged = sc.log.filter(function (l) { return /^(ERROR|MISSING DATA|INCONSISTENT|SINGULAR|WARNING)/.test(l); });
    var msgs = flagged.length ? flagged.map(function (l) {
      return "<p class='msg " + (/^ERROR|^INCONSISTENT/.test(l) ? "err" : "") + "'>" + esc(l) + "</p>"; }).join("")
      : "<p class='msg ok'>Completed with no warnings.</p>";
    var rows = sc.rows.length ? "<table><thead><tr><th>Quantity</th><th>Value</th><th>Unit</th><th>Where it came from</th></tr></thead><tbody>" +
      sc.rows.map(function (r) {
        var v = r[1] === "nan" ? "symbolic" : sig(+r[1], 5) + (r[4] && +r[4] > 0 ? " ± " + sig(+r[4], 2) : "");
        return "<tr><td>" + esc(r[0]) + "</td><td>" + v + "</td><td>" + esc(r[2]) + "</td><td>" + esc(r[3]) + "</td></tr>";
      }).join("") + "</tbody></table>" : "<p class='note' style='margin:0'>No values were returned.</p>";
    $("sview").innerHTML = "<h3>What came back " + statusChip(sc.status) + "</h3>" + msgs +
      "<div class='two' style='margin-top:16px'><div><h3>Input file</h3><div class='box'><pre>" + esc(sc.input) +
      "</pre></div></div><div><h3>First rows of the output file</h3><div class='box'>" + rows + "</div></div></div>";
  }

  /* ------------------------------------------------------------ uncertainty */
  var IN_COLOR = { T1: "#f5c043", P1: "#5aa2ff", T2: "#f29bc0", P2: "#ff8a5c", Tc: "#7fd66a", Pc: "#39d6a0", omega: "#a590ff",
    V: "#ff8a5c", Cp: "#5aa2ff", w: "#f5c043", muJT: "#39d6a0" };
  var IN_NAME = { T1: "initial temperature", P1: "initial pressure", T2: "final temperature", P2: "final pressure", Tc: "critical temperature",
    Pc: "critical pressure", omega: "acentric factor", V: "volume", Cp: "heat capacity", w: "sound speed", muJT: "Joule–Thomson coefficient" };
  var Q_NAME = { dH: "Δh", dS: "Δs", dU: "Δu", dG: "Δg", dV: "Δv", alpha: "α", kappaT: "κ<sub>T</sub>", Cv: "C<sub>V</sub>", gamma: "γ",
    "(dT/dP)_S": "(∂T/∂P)<sub>S</sub>" };
  function renderBudget(key, keysId, barsId, wanted) {
    var b = D.budget[key];
    $(keysId).innerHTML = b.inputs.map(function (k) {
      return "<span><i style='background:" + IN_COLOR[k] + "'></i>" + IN_NAME[k] + "</span>"; }).join("");
    $(barsId).innerHTML = wanted.map(function (q) {
      var r = b.rows.filter(function (x) { return x[0] === q; })[0];
      var segs = b.inputs.map(function (k, i) {
        var pct = 100 * r[3 + i];
        return pct < 0.4 ? "" : "<i style='flex:" + pct + " 1 0;background:" + IN_COLOR[k] + "' title='" + IN_NAME[k] + ": " + pct.toFixed(0) + "% of the variance'></i>";
      }).join("");
      var top = b.inputs.map(function (k, i) { return [k, r[3 + i]]; }).sort(function (a, c) { return c[1] - a[1]; })[0];
      return "<div class='barrow'><span class='q'>" + Q_NAME[q] + "</span><div class='stack'>" + segs + "</div><span class='u'>" +
        (100 * top[1]).toFixed(0) + "% " + IN_NAME[top[0]] + "</span></div>";
    }).join("");
  }

  /* ------------------------------------------------------------ WebGL scene */
  var scene = (function () {
    var canvas = $("scene"), gl = null, still = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    try { gl = canvas.getContext("webgl", { antialias: true, alpha: true, premultipliedAlpha: false }); } catch (err) { gl = null; }
    var api = { setMarker: function () {} };
    if (!gl) { canvas.style.display = "none"; return api; }
    var pr = E.makeEos("PR", CO2), b = pr.b;
    var LV0 = Math.log10(1.25 * b), LV1 = Math.log10(70 * b), PMAX = 24e6, T0 = 222, T1 = 438;
    function X(V) { return -1.7 + 3.4 * (Math.log10(V) - LV0) / (LV1 - LV0); }
    function Y(P) { return -0.75 + 1.9 * P / PMAX; }
    function Zc(T) { return -1.35 + 2.7 * (T - T0) / (T1 - T0); }
    var ramp = [[0.16, 0.26, 0.95], [0.56, 0.29, 1.0], [1.0, 0.31, 0.51], [1.0, 0.65, 0.23], [1.0, 0.95, 0.79]];
    function color(t) {
      t = Math.max(0, Math.min(0.9999, t)) * 4;
      var i = Math.floor(t), f = t - i;
      return [0, 1, 2].map(function (k) { return ramp[i][k] + (ramp[i + 1][k] - ramp[i][k]) * f; });
    }
    var lines = [], points = [];
    function seg(p, q, c, u1, u2, a) { lines.push(p[0], p[1], p[2], c[0], c[1], c[2], u1, a, q[0], q[1], q[2], c[0], c[1], c[2], u2, a); }
    var NV = 110, grid = [], Tlist = [];
    for (var T = 225; T <= 435; T += 5) Tlist.push(T);
    Tlist.forEach(function (T, ti) {
      var sat = T < pr.Tc ? pr.psat(T) : null, row = [], c = color((T - T0) / (T1 - T0));
      for (var i = 0; i < NV; i++) {
        var V = Math.pow(10, LV0 + (LV1 - LV0) * i / (NV - 1));
        var P = sat && V > sat.VL && V < sat.VV ? sat.P : pr.pressure(T, V);
        row.push(P > 0 && P <= PMAX ? [X(V), Y(P), Zc(T)] : null);
      }
      grid.push(row);
      for (i = 0; i < NV - 1; i++) if (row[i] && row[i + 1]) seg(row[i], row[i + 1], c, i / NV + ti * 0.37, (i + 1) / NV + ti * 0.37, 1.0);
      for (i = 0; i < NV; i++) if (row[i]) points.push(row[i][0], row[i][1], row[i][2], c[0], c[1], c[2], i / NV + ti * 0.37, 0.8);
    });
    for (var i = 0; i < NV; i += 5) for (var ti = 0; ti < Tlist.length - 1; ti++) {
      var p = grid[ti][i], q = grid[ti + 1][i];
      if (p && q) seg(p, q, color((Tlist[ti] - T0) / (T1 - T0)), 0.5, 0.5, 0.22);
    }
    var prev = null;                       // coexistence dome
    for (T = 225; T <= 303.5; T += 1.5) {
      var s = pr.psat(T), cur = [[X(s.VL), Y(s.P), Zc(T)], [X(s.VV), Y(s.P), Zc(T)]];
      if (prev) { seg(prev[0], cur[0], [1, 0.95, 0.8], 0.2, 0.2, 1); seg(prev[1], cur[1], [1, 0.95, 0.8], 0.2, 0.2, 1); }
      prev = cur;
    }
    seg(prev[0], prev[1], [1, 0.95, 0.8], 0.2, 0.2, 1);
    var marker = new Float32Array([0, 0, 0, 1, 0.97, 0.85, 0, 1]);

    function shader(type, src) { var sh = gl.createShader(type); gl.shaderSource(sh, src); gl.compileShader(sh); return sh; }
    var prog = gl.createProgram();
    gl.attachShader(prog, shader(gl.VERTEX_SHADER,
      "attribute vec3 p;attribute vec3 c;attribute vec2 u;uniform mat4 m;uniform float t;uniform float size;varying vec4 v;" +
      "void main(){vec4 q=m*vec4(p,1.0);gl_Position=q;float pulse=0.6+0.9*pow(0.5+0.5*sin(6.2832*(u.x*3.0-t)),8.0);" +
      "float fog=clamp(1.9-0.2*q.w,0.4,1.0);v=vec4(c*pulse,u.y*fog);gl_PointSize=size*(6.0/q.w);}"));
    gl.attachShader(prog, shader(gl.FRAGMENT_SHADER,
      "precision mediump float;varying vec4 v;uniform float round;void main(){float a=v.a;" +
      "if(round>0.5){float d=length(gl_PointCoord-0.5);a*=smoothstep(0.5,0.0,d);}gl_FragColor=vec4(v.rgb,a);}"));
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) { canvas.style.display = "none"; return api; }
    gl.useProgram(prog);
    function buffer(data) { var bf = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, bf); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(data), gl.DYNAMIC_DRAW); return bf; }
    var bl = buffer(lines), bp = buffer(points), bm = buffer(marker);
    var ap = gl.getAttribLocation(prog, "p"), ac = gl.getAttribLocation(prog, "c"), au = gl.getAttribLocation(prog, "u");
    var um = gl.getUniformLocation(prog, "m"), ut = gl.getUniformLocation(prog, "t"), us = gl.getUniformLocation(prog, "size"), ur = gl.getUniformLocation(prog, "round");
    function bind(bf) {
      gl.bindBuffer(gl.ARRAY_BUFFER, bf);
      gl.enableVertexAttribArray(ap); gl.vertexAttribPointer(ap, 3, gl.FLOAT, false, 32, 0);
      gl.enableVertexAttribArray(ac); gl.vertexAttribPointer(ac, 3, gl.FLOAT, false, 32, 12);
      gl.enableVertexAttribArray(au); gl.vertexAttribPointer(au, 2, gl.FLOAT, false, 32, 24);
    }
    gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE); gl.disable(gl.DEPTH_TEST);
    var mx = 0, my = 0, tx = 0, ty = 0, scroll = 0;
    window.addEventListener("mousemove", function (ev) { tx = ev.clientX / window.innerWidth - 0.5; ty = ev.clientY / window.innerHeight - 0.5; });
    window.addEventListener("scroll", function () {
      scroll = window.scrollY / Math.max(1, document.body.scrollHeight - window.innerHeight);
      document.body.classList.toggle("reading", window.scrollY > window.innerHeight * 0.55);
    }, { passive: true });
    function resize() {
      var dpr = Math.min(2, window.devicePixelRatio || 1);
      canvas.width = Math.round(window.innerWidth * dpr); canvas.height = Math.round(window.innerHeight * dpr);
      gl.viewport(0, 0, canvas.width, canvas.height);
    }
    window.addEventListener("resize", resize); resize();
    function matrix(time) {
      var asp = canvas.width / canvas.height, wide = asp > 1.15;
      var az = -0.62 + 2.2 * scroll + 0.35 * mx + (still ? 0 : 0.12 * Math.sin(time * 0.17));
      var elv = 0.42 + 0.18 * my, rad = wide ? 4.3 : 6.2;
      var ca = Math.cos(az), sa = Math.sin(az), ce = Math.cos(elv), se = Math.sin(elv);
      // view: rotate about y by az, then about x by elv, then push back
      var f = 1 / Math.tan(0.5 * 0.72), n = 0.1, fr = 40, ox = wide ? 0.36 : 0, oy = wide ? 0.02 : 0.3;
      var r = [ca, 0, -sa, sa * se, ce, ca * se, sa * ce, -se, ca * ce];    // rows of the rotation
      var v = [r[0], r[1], r[2], 0, r[3], r[4], r[5], 0, r[6], r[7], r[8], -rad];
      var m = new Float32Array(16);
      // clip = P * V, written column-major
      for (var col = 0; col < 4; col++) {
        var vx = v[col], vy = v[4 + col], vz = v[8 + col], vw = col === 3 ? 1 : 0;
        var cw = -vz;
        m[col * 4] = f / asp * vx + ox * cw;
        m[col * 4 + 1] = f * vy + oy * cw;
        m[col * 4 + 2] = (fr + n) / (n - fr) * vz + (2 * fr * n) / (n - fr) * vw;
        m[col * 4 + 3] = cw;
      }
      return m;
    }
    function frame(ms) {
      var time = ms / 1000;
      mx += (tx - mx) * 0.05; my += (ty - my) * 0.05;
      gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT);
      gl.uniformMatrix4fv(um, false, matrix(time));
      gl.uniform1f(ut, still ? 0.3 : time * 0.11);
      gl.uniform1f(ur, 0); gl.uniform1f(us, 1); bind(bl); gl.drawArrays(gl.LINES, 0, lines.length / 8);
      gl.uniform1f(ur, 1); gl.uniform1f(us, 7 * Math.min(2, window.devicePixelRatio || 1)); bind(bp); gl.drawArrays(gl.POINTS, 0, points.length / 8);
      gl.uniform1f(us, 30 * Math.min(2, window.devicePixelRatio || 1)); bind(bm); gl.drawArrays(gl.POINTS, 0, 1); gl.drawArrays(gl.POINTS, 0, 1);
      if (!still || frame.first) { frame.first = false; requestAnimationFrame(frame); }
    }
    frame.first = true;
    requestAnimationFrame(frame);
    if (still) { window.addEventListener("scroll", function () { requestAnimationFrame(frame); }, { passive: true }); }
    api.setMarker = function (st) {
      var s = E.state(CO2, pr, st.T, st.P);      // the surface is Peng–Robinson
      marker[0] = X(s.V); marker[1] = Y(st.P); marker[2] = Zc(st.T);
      gl.bindBuffer(gl.ARRAY_BUFFER, bm); gl.bufferSubData(gl.ARRAY_BUFFER, 0, marker);
      if (still) requestAnimationFrame(frame);
    };
    return api;
  })();

  /* ------------------------------------------------------------ wire up */
  chips("cx", "x"); chips("cy", "y"); chips("cz", "z");
  seg("eos", MODELS.slice().reverse().map(function (k) { return [k, MODEL_NAME[k]]; }), "eos", renderState);
  seg("vprop", Object.keys(VP).map(function (k) { return [k, VP[k][0]]; }), "vprop", renderGrid);
  seg("veos", MODELS.slice().reverse().map(function (k) { return [k, MODEL_NAME[k]]; }), "veos", renderGrid);
  $("sT").addEventListener("input", function (ev) { S.T = +ev.target.value; renderState(); });
  $("sP").addEventListener("input", function (ev) { S.P = +ev.target.value; renderState(); });
  renderState();
  renderGrid();
  renderScenarios();
  renderBudget("02_process_co2_pr", "k1", "b1", ["dH", "dS", "dU", "dG", "dV"]);
  renderBudget("05_infer_from_measurements", "k2", "b2", ["alpha", "kappaT", "Cv", "gamma", "(dT/dP)_S"]);
  window.addEventListener("resize", renderGrid);

  var links = document.querySelectorAll(".bar nav a");
  var io = new IntersectionObserver(function (entries) {
    entries.forEach(function (en) {
      if (en.isIntersecting) Array.prototype.forEach.call(links, function (a) {
        a.classList.toggle("on", a.getAttribute("href") === "#" + en.target.id);
      });
    });
  }, { rootMargin: "-45% 0px -50% 0px" });
  document.querySelectorAll("section").forEach(function (s) { io.observe(s); });
})();
