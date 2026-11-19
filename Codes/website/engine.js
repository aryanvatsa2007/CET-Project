/* JavaScript port of the equation-of-state part of the thermocalculus engine.
   It mirrors Codes/thermocalculus/eos.py and props.py so that the website shows
   the same numbers as the Python package; check_engine.js compares the two.
   Reduced derivatives are NOT re-derived here: they are read from data.js,
   which run_all.py writes from the Python symbolic module. */
(function (root) {
  "use strict";
  var R = 8.314462618, T_REF = 298.15, P_REF = 1e5;
  var CB = Math.cbrt(2);
  var KINDS = {
    vdW: [0, 0, 27 / 64, 1 / 8],
    RK: [0, 1, 1 / (9 * (CB - 1)), (CB - 1) / 3],
    SRK: [0, 1, 1 / (9 * (CB - 1)), (CB - 1) / 3],
    PR: [1 - Math.SQRT2, 1 + Math.SQRT2, 0.45724, 0.07780]
  };

  function bisect(f, lo, hi, n) {
    var flo = f(lo);
    for (var i = 0; i < (n || 200); i++) {
      var mid = 0.5 * (lo + hi), fm = f(mid);
      if ((fm > 0) === (flo > 0)) { lo = mid; flo = fm; } else { hi = mid; }
      if (Math.abs(hi - lo) <= 1e-15 * Math.abs(hi)) break;
    }
    return 0.5 * (lo + hi);
  }

  function makeEos(kind, fl) {
    var e = { kind: kind, ideal: kind === "IG", Tc: fl.Tc, Pc: fl.Pc, b: 0 };
    if (!e.ideal) {
      var k = KINDS[kind];
      e.eps = k[0]; e.sig = k[1];
      e.ac = k[2] * R * R * fl.Tc * fl.Tc / fl.Pc;
      e.b = k[3] * R * fl.Tc / fl.Pc;
      if (kind === "SRK") e.m = 0.480 + 1.574 * fl.omega - 0.176 * fl.omega * fl.omega;
      if (kind === "PR") e.m = 0.37464 + 1.54226 * fl.omega - 0.26992 * fl.omega * fl.omega;
    }
    e.a = function (T) {
      if (e.ideal) return [0, 0, 0];
      if (kind === "vdW") return [e.ac, 0, 0];
      var Tr = T / e.Tc;
      if (kind === "RK") return [e.ac / Math.sqrt(Tr), -0.5 * e.ac * Math.pow(Tr, -1.5) / e.Tc,
        0.75 * e.ac * Math.pow(Tr, -2.5) / (e.Tc * e.Tc)];
      var s = Math.sqrt(Tr), g = 1 + e.m * (1 - s);
      return [e.ac * g * g, -e.ac * e.m * g / (s * e.Tc),
        e.ac * e.m * (1 + e.m) / (2 * e.Tc * e.Tc * s * s * s)];
    };
    var den = function (V) { return (V + e.eps * e.b) * (V + e.sig * e.b); };
    e.pressure = function (T, V) {
      return e.ideal ? R * T / V : R * T / (V - e.b) - e.a(T)[0] / den(V);
    };
    e.dPdT = function (T, V) {
      return e.ideal ? R / V : R / (V - e.b) - e.a(T)[1] / den(V);
    };
    e.dPdV = function (T, V) {
      if (e.ideal) return -R * T / (V * V);
      return -R * T / ((V - e.b) * (V - e.b)) +
        e.a(T)[0] * (2 * V + (e.eps + e.sig) * e.b) / (den(V) * den(V));
    };
    e.volumes = function (T, P) {
      if (e.ideal) return [R * T / P];
      var A = e.a(T)[0] * P / (R * T * R * T), B = e.b * P / (R * T), s = e.eps, g = e.sig;
      var c2 = (s + g - 1) * B - 1, c1 = s * g * B * B - (s + g) * (B * B + B) + A,
        c0 = -(s * g * (B * B * B + B * B) + A * B);
      // real roots of Z^3 + c2 Z^2 + c1 Z + c0 (trigonometric / Cardano)
      var p = c1 - c2 * c2 / 3, q = 2 * c2 * c2 * c2 / 27 - c2 * c1 / 3 + c0, zs = [];
      var disc = q * q / 4 + p * p * p / 27;
      if (disc > 0) {
        var sq = Math.sqrt(disc);
        zs.push(Math.cbrt(-q / 2 + sq) + Math.cbrt(-q / 2 - sq) - c2 / 3);
      } else {
        var r = 2 * Math.sqrt(-p / 3), th = Math.acos(Math.max(-1, Math.min(1, 3 * q / (p * r)))) / 3;
        for (var i = 0; i < 3; i++) zs.push(r * Math.cos(th - 2 * Math.PI * i / 3) - c2 / 3);
      }
      var out = [];
      zs.forEach(function (z) {
        if (!(z > B)) return;
        var V = z * R * T / P;
        for (var it = 0; it < 60; it++) {
          var dV = (e.pressure(T, V) - P) / e.dPdV(T, V);
          if (!(V - dV > e.b)) break;
          V -= dV;
          if (Math.abs(dV) < 1e-15 * V) break;
        }
        if (V > e.b && Math.abs(e.pressure(T, V) - P) < 1e-8 * P) out.push(V);
      });
      out.sort(function (a, b) { return a - b; });
      return out.filter(function (v, i) { return i === 0 || Math.abs(v - out[i - 1]) > 1e-9 * v; });
    };
    e.q = function (V) {
      if (kind === "vdW") return 1 / V;
      return Math.log((V + e.sig * e.b) / (V + e.eps * e.b)) / (e.b * (e.sig - e.eps));
    };
    e.residual = function (T, V) {
      if (e.ideal) return { AR: 0, SR: 0, UR: 0, CvR: 0 };
      var a = e.a(T), q = e.q(V), ln = Math.log(1 - e.b / V);
      return { AR: -R * T * ln - a[0] * q, SR: R * ln + a[1] * q, UR: -(a[0] - T * a[1]) * q, CvR: T * a[2] * q };
    };
    e.lnphi = function (T, V) {
      if (e.ideal) return 0;
      var Z = e.pressure(T, V) * V / (R * T);
      return e.residual(T, V).AR / (R * T) + Z - 1 - Math.log(Z);
    };
    e.volume = function (T, P) {
      var vs = e.volumes(T, P), best = vs[0];
      vs.forEach(function (v) { if (e.lnphi(T, v) < e.lnphi(T, best)) best = v; });
      return best;
    };
    e.spinodal = function (T) {
      var prev = null, pv = 0, hits = [];
      for (var i = 0; i <= 3000; i++) {
        var V = e.b * (1 + Math.pow(10, -4 + 8 * i / 3000)), d = e.dPdV(T, V);
        if (prev !== null && (d > 0) !== (pv > 0)) hits.push([prev, V]);
        prev = V; pv = d;
      }
      if (hits.length < 2) return null;
      var f = function (V) { return e.dPdV(T, V); };
      return [bisect(f, hits[0][0], hits[0][1]), bisect(f, hits[hits.length - 1][0], hits[hits.length - 1][1])];
    };
    e.psat = function (T) {
      if (e.ideal || T >= e.Tc) return null;
      var sp = e.spinodal(T);
      if (!sp) return null;
      var phi = e.pressure(T, sp[1]), plo = Math.max(e.pressure(T, sp[0]), 1e-10 * phi);
      var roots = function (P) {
        var f = function (V) { return e.pressure(T, V) - P; };
        return [bisect(f, e.b * (1 + 1e-13), sp[0]), bisect(f, sp[1], 10 * R * T / P + sp[1])];
      };
      var gap = function (P) { var v = roots(P); return e.lnphi(T, v[0]) - e.lnphi(T, v[1]); };
      var P = bisect(gap, plo + 1e-9 * (phi - plo), phi - 1e-9 * (phi - plo), 120);
      var v = roots(P);
      return { P: P, VL: v[0], VV: v[1] };
    };
    return e;
  }

  function cp0(fl, T) {
    var c = fl.cp, t = T / 1000;
    return c[0] + c[1] * t + c[2] * t * t + c[3] * t * t * t + c[4] / (t * t);
  }
  function hInt(fl, T) {
    var c = fl.cp, t = T / 1000;
    return 1000 * (c[0] * t + c[1] * t * t / 2 + c[2] * t * t * t / 3 + c[3] * t * t * t * t / 4 - c[4] / t);
  }
  function sInt(fl, T) {
    var c = fl.cp, t = T / 1000;
    return c[0] * Math.log(t) + c[1] * t + c[2] * t * t / 2 + c[3] * t * t * t / 3 - c[4] / (2 * t * t);
  }

  function state(fl, eos, T, P) {
    var V = eos.volume(T, P), dPdT = eos.dPdT(T, V), dPdV = eos.dPdV(T, V);
    var Z = P * V / (R * T), res = eos.residual(T, V);
    var o = { T: T, P: P, V: V, Z: Z, alpha: -dPdT / (V * dPdV), kappaT: -1 / (V * dPdV) };
    o.HR = res.UR + P * V - R * T;
    o.SR = res.SR + R * Math.log(Z);
    var c0 = cp0(fl, T);
    o.Cv = c0 - R + res.CvR;
    o.Cp = o.Cv - T * dPdT * dPdT / dPdV;
    o.kappaS = o.kappaT * o.Cv / o.Cp;
    o.w = Math.sqrt(V / (fl.M * o.kappaS));
    o.muJT = V * (T * o.alpha - 1) / o.Cp;
    o.H = hInt(fl, T) - hInt(fl, T_REF) + o.HR;
    o.S = sInt(fl, T) - sInt(fl, T_REF) - R * Math.log(P / P_REF) + o.SR;
    o.U = o.H - P * V; o.A = o.U - T * o.S; o.G = o.H - T * o.S;
    o.phase = eos.ideal ? "ideal gas" : (T >= eos.Tc ? (P >= eos.Pc ? "supercritical" : "gas") : "subcritical");
    return o;
  }

  var SYMS = ["T", "P", "V", "S", "Cp", "alpha", "kappaT"];
  function evalPoly(terms, vals) {
    var tot = 0;
    terms.forEach(function (t) {
      var x = t[0];
      for (var i = 0; i < 7; i++) if (t[1][i]) x *= Math.pow(vals[SYMS[i]], t[1][i]);
      tot += x;
    });
    return tot;
  }
  /* value of a reduced derivative record from data.js; null if singular */
  function evalDerivative(rec, vals) {
    var d = evalPoly(rec.d, vals), size = 0;
    rec.d.forEach(function (t) { size += Math.abs(evalPoly([t], vals)); });
    if (Math.abs(d) <= 1e-12 * size) return null;
    return evalPoly(rec.n, vals) / d;
  }

  var api = { R: R, makeEos: makeEos, state: state, evalDerivative: evalDerivative, cp0: cp0, SYMS: SYMS };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.TCEngine = api;
})(typeof window !== "undefined" ? window : globalThis);
