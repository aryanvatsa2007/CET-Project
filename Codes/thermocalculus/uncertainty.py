"""Propagation of input uncertainty.

First-order (GUM) propagation for uncorrelated inputs,

    u_y^2 = sum_i (dy/dx_i)^2 u_i^2 ,

with sensitivities from central differences, plus a Monte Carlo check that does
not rely on linearisation.  The fraction (dy/dx_i u_i)^2 / u_y^2 is reported for
every input so the dominant contribution is visible.
"""
import math

import numpy as np


def propagate(func, x, u, rel_step=1e-5):
    """Linear propagation.

    func : callable taking a dict {name: value} and returning a dict of outputs
    x    : nominal inputs
    u    : standard uncertainties of some of the inputs (same units as x)

    Returns {output: (value, u, {input: fractional variance contribution})}.
    """
    y0 = func(dict(x))
    sens = {}
    for name, ui in u.items():
        if ui == 0:
            continue
        h = rel_step * (abs(x[name]) if x[name] != 0 else 1.0)
        hi, lo = dict(x), dict(x)
        hi[name] += h
        lo[name] -= h
        yh, yl = func(hi), func(lo)
        sens[name] = {k: (yh[k] - yl[k]) / (2 * h) for k in y0}
    out = {}
    for k, val in y0.items():
        parts = {n: (sens[n][k] * u[n]) ** 2 for n in sens}
        var = sum(parts.values())
        out[k] = (val, math.sqrt(var),
                  {n: (p / var if var > 0 else 0.0) for n, p in parts.items()})
    return out


def monte_carlo(func, x, u, n=2000, seed=12345):
    """Monte Carlo propagation with independent normal inputs.

    Returns {output: (mean, standard deviation)} and the number of failed draws.
    """
    rng = np.random.default_rng(seed)
    names = [k for k in u if u[k] != 0]
    draws = {k: rng.normal(x[k], u[k], n) for k in names}
    samples, failed = {}, 0
    for i in range(n):
        xi = dict(x)
        for k in names:
            xi[k] = float(draws[k][i])
        try:
            yi = func(xi)
        except Exception:
            failed += 1
            continue
        for k, v in yi.items():
            samples.setdefault(k, []).append(v)
    return ({k: (float(np.mean(v)), float(np.std(v, ddof=1)))
             for k, v in samples.items()}, failed)
