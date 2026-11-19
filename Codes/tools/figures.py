"""Manuscript figures.  Every figure is drawn only from CSV files in Data/, so a
figure can never disagree with the data files.  Sizes follow the ACS graphics
rules: single column 3.33 in, double column up to 7 in, lettering >= 4.5 pt,
lines >= 0.5 pt.  Called by run_all.py.
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

COL1, COL2 = 3.33, 7.0
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
# fixed colour and line style per model (never re-assigned by rank)
STYLE = {
    "PR":  dict(color="#2a78d6", ls="-"),
    "SRK": dict(color="#eb6834", ls="--"),
    "RK":  dict(color="#1baf7a", ls="-."),
    "vdW": dict(color="#eda100", ls=":"),
    "IG":  dict(color="#e87ba4", ls=(0, (1, 2.5))),
}
INPUT_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
DIVERGING = LinearSegmentedColormap.from_list(
    "div", ["#1f5fb0", "#8fb8ea", "#ecebe8", "#f3a98b", "#c8431a"])

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "axes.labelsize": 7.5,
    "axes.titlesize": 7.5, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "axes.linewidth": 0.6, "lines.linewidth": 1.1,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "legend.frameon": False, "savefig.dpi": 600, "pdf.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False,
    "mathtext.fontset": "dejavusans",
})


def rcsv(root, rel):
    with open(os.path.join(root, "Data", rel), newline="") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def save(fig, root, name):
    fig.savefig(os.path.join(root, f"{name}.pdf"), bbox_inches="tight", pad_inches=0.02, dpi=600)
    plt.close(fig)
    print(f"  wrote {name}.pdf")


def tag(ax, letter, x=-0.02, y=1.04):
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=8, fontweight="bold",
            ha="right", va="bottom")


def nist_marker(ax, x, y, label="NIST reference"):
    return ax.plot(x, y, ls="none", marker="o", ms=2.8, mfc="white", mec=INK, mew=0.7,
                   label=label, zorder=5)[0]


# --------------------------------------------------------------- Fig. 2
def fig_verification(root):
    rows = rcsv(root, "verification/derivatives_336_check.csv")
    cf = rcsv(root, "verification/closedform_vs_numerical.csv")
    fig, (a, b) = plt.subplots(2, 1, figsize=(COL1, 4.3), gridspec_kw=dict(hspace=0.55))
    for state, color, ls in (("PR supercritical", "#2a78d6", "-"), ("PR gas", "#eb6834", "--")):
        dev = sorted(max(float(r["relative_deviation"]), 1e-16) for r in rows
                     if r["state"] == state and r["status"] == "checked")
        T = next(r for r in rows if r["state"] == state)
        a.semilogy(np.arange(1, len(dev) + 1), dev, color=color, ls=ls,
                   label=f"{float(T['T_K']):.0f} K, {float(T['P_Pa']) / 1e6:.0f} MPa")
    a.set_xlabel("derivative, ranked by deviation")
    a.set_ylabel("relative deviation,\nformula vs numerical")
    a.set_xlim(0, 340)
    a.set_ylim(1e-13, 1e-3)
    a.grid(axis="y", color=GRID, lw=0.5)
    a.legend(loc="upper left", title="CO$_2$, Peng–Robinson state", title_fontsize=6.5)
    tag(a, "a", x=-0.2)
    props = ["alpha", "kappaT", "HR", "SR", "GR", "Cv", "Cp", "w", "muJT"]
    labels = [r"$\alpha$", r"$\kappa_T$", r"$H^R$", r"$S^R$", r"$G^R$", r"$C_V$",
              r"$C_P$", r"$w$", r"$\mu_{JT}$"]
    for k, kind in enumerate(("PR", "SRK", "RK", "vdW")):
        y = [max(max(float(r["relative_difference"]) for r in cf
                     if r["eos"] == kind and r["property"] == p), 1e-16) for p in props]
        b.semilogy(np.arange(len(props)) + (k - 1.5) * 0.16, y, ls="none",
                   marker="osD^"[k], ms=3.2, color=STYLE[kind]["color"], label=kind)
    b.set_xticks(range(len(props)))
    b.set_xticklabels(labels)
    b.set_ylabel("max relative difference,\nclosed form vs numerical")
    b.set_ylim(1e-15, 1e-7)
    b.grid(axis="y", color=GRID, lw=0.5)
    b.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.22), handletextpad=0.2,
             columnspacing=1.0)
    tag(b, "b", x=-0.2)
    save(fig, root, "fig2_verification")


# --------------------------------------------------------------- Fig. 3
def _eos_rows(root):
    return rcsv(root, "validation/co2_eos_vs_nist.csv")


def fig_isotherms(root):
    rows = _eos_rows(root)
    curves = rcsv(root, "validation/co2_model_curves.csv")
    fig, axes = plt.subplots(1, 3, figsize=(COL2, 2.35), sharey=True,
                             gridspec_kw=dict(wspace=0.08))
    Rg = 8.314462618
    for ax, T, letter in zip(axes, (310.0, 350.0, 400.0), "abc"):
        sel = [r for r in rows if float(r["T_K"]) == T and r["property"] == "rho"]
        P = sorted({float(r["P_MPa"]) for r in sel})
        zof = lambda rho, p: p * 1e6 / (rho * 1e3 * Rg * T)
        for kind in ("IG", "vdW", "RK", "SRK", "PR"):
            pts = sorted((float(r["P_MPa"]), float(r["Z"])) for r in curves
                         if r["model"] == kind and r["curve"] == "isotherm"
                         and float(r["T_K"]) == T)
            ax.plot(*zip(*pts), label=kind, **STYLE[kind])
        pts = sorted((float(r["P_MPa"]), float(r["nist_value"])) for r in sel
                     if r["model"] == "PR")
        nist_marker(ax, [p for p, _ in pts], [zof(v, p) for p, v in pts])
        ax.set_title(f"{T:.0f} K", pad=3)
        ax.set_xlabel("pressure, $P$ (MPa)")
        ax.set_xlim(0, 20.5)
        ax.grid(color=GRID, lw=0.5)
        tag(ax, letter, x=0.02 if letter != "a" else -0.1)
    axes[0].set_ylabel("compressibility factor, $Z$")
    axes[0].set_ylim(0.15, 1.05)
    h, l = axes[0].get_legend_handles_labels()
    order = [5, 4, 3, 2, 1, 0]
    fig.legend([h[i] for i in order], [l[i] for i in order], ncol=6, loc="upper center",
               bbox_to_anchor=(0.5, 1.08), handlelength=2.6, columnspacing=1.6)
    save(fig, root, "fig3_co2_isotherms")


# --------------------------------------------------------------- Fig. 4
def _grid(rows, model, prop, col):
    T = sorted({float(r["T_K"]) for r in rows})
    P = sorted({float(r["P_MPa"]) for r in rows})
    g = np.full((len(P), len(T)), np.nan)
    for r in rows:
        if r["model"] == model and r["property"] == prop and r[col] not in ("", "nan"):
            g[P.index(float(r["P_MPa"])), T.index(float(r["T_K"]))] = float(r[col])
    return T, P, g


def _heat(ax, T, P, g, norm, cmap=DIVERGING):
    im = ax.pcolormesh(np.arange(len(T) + 1), np.arange(len(P) + 1), g, cmap=cmap,
                       norm=norm, edgecolors="white", linewidth=0.4)
    ax.set_xticks(np.arange(len(T))[::2] + 0.5)
    ax.set_xticklabels([f"{t:.0f}" for t in T[::2]])
    ax.set_yticks(np.arange(len(P)) + 0.5)
    ax.set_yticklabels([f"{p:g}" for p in P])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlabel("temperature, $T$ (K)")
    return im


def fig_density_map(root):
    rows = _eos_rows(root)
    fig, axes = plt.subplots(1, 2, figsize=(COL2, 2.25), gridspec_kw=dict(wspace=0.12))
    norm = TwoSlopeNorm(vmin=-12, vcenter=0, vmax=12)
    for ax, kind, letter in zip(axes, ("PR", "SRK"), "ab"):
        T, P, g = _grid(rows, kind, "rho", "deviation_percent")
        im = _heat(ax, T, P, g, norm)
        worst = np.nanargmax(np.abs(g))
        i, j = np.unravel_index(worst, g.shape)
        ax.text(j + 0.5, i + 0.5, f"{g[i, j]:.0f}".replace("-", "–"), ha="center",
                va="center", fontsize=5.2, color="white")
        ax.set_title({"PR": "Peng–Robinson", "SRK": "Soave–Redlich–Kwong"}[kind], pad=3)
        tag(ax, letter, x=-0.02)
    axes[0].set_ylabel("pressure, $P$ (MPa)")
    axes[1].set_yticklabels([])
    cb = fig.colorbar(im, ax=axes, pad=0.015, aspect=22, extend="both")
    cb.set_label(r"density deviation, 100($\rho_{\mathrm{EoS}}$/$\rho_{\mathrm{NIST}}$ – 1) (%)")
    cb.outline.set_visible(False)
    save(fig, root, "fig4_co2_density_deviation")


# --------------------------------------------------------------- Fig. 5
def fig_isobar(root):
    rows = [r for r in _eos_rows(root) if float(r["P_MPa"]) == 10.0]
    curves = [r for r in rcsv(root, "validation/co2_model_curves.csv") if r["curve"] == "isobar"]
    ccol = {"Cp": "Cp_J_per_molK", "w": "w_m_per_s", "muJT": "muJT_K_per_MPa"}
    fig, axes = plt.subplots(1, 3, figsize=(COL2, 2.3), gridspec_kw=dict(wspace=0.36))
    spec = [("Cp", "isobaric heat capacity,\n$C_P$ (J mol$^{-1}$ K$^{-1}$)", (0, 360)),
            ("w", "speed of sound, $w$ (m s$^{-1}$)", (150, 420)),
            ("muJT", "Joule–Thomson coefficient,\n$\\mu_{JT}$ (K MPa$^{-1}$)", (0, 10))]
    for ax, (prop, ylabel, ylim), letter in zip(axes, spec, "abc"):
        for kind in ("vdW", "RK", "SRK", "PR"):
            pts = sorted((float(r["T_K"]), float(r[ccol[prop]])) for r in curves
                         if r["model"] == kind)
            ax.plot(*zip(*pts), label=kind, **STYLE[kind])
        pts = sorted((float(r["T_K"]), float(r["nist_value"])) for r in rows
                     if r["model"] == "PR" and r["property"] == prop)
        nist_marker(ax, *zip(*pts))
        ax.set_xlabel("temperature, $T$ (K)")
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        ax.set_xlim(308, 402)
        ax.grid(color=GRID, lw=0.5)
        tag(ax, letter, x=-0.2)
    h, l = axes[0].get_legend_handles_labels()
    order = [4, 3, 2, 1, 0]
    fig.legend([h[i] for i in order], [l[i] for i in order], ncol=5, loc="upper center",
               bbox_to_anchor=(0.5, 1.07), handlelength=2.6, columnspacing=1.6)
    save(fig, root, "fig5_co2_isobar_10MPa")


# --------------------------------------------------------------- Fig. 6
def fig_saturation(root):
    rows = rcsv(root, "validation/co2_saturation_vs_nist.csv")
    fig, (a, b) = plt.subplots(2, 1, figsize=(COL1, 4.3), gridspec_kw=dict(hspace=0.42))
    curves = rcsv(root, "validation/co2_saturation_curves.csv")
    for kind in ("vdW", "RK", "SRK", "PR"):
        c = [r for r in curves if r["model"] == kind]
        Tc = [float(r["T_K"]) for r in c]
        a.semilogy(Tc, [float(r["Psat_MPa"]) for r in c], label=kind, **STYLE[kind])
        b.plot([float(r["rhoV_mol_per_L"]) for r in c], Tc, **STYLE[kind])
        b.plot([float(r["rhoL_mol_per_L"]) for r in c], Tc, **STYLE[kind])
    pts = sorted((float(r["T_K"]), float(r["nist_value"])) for r in rows
                 if r["model"] == "PR" and r["property"] == "Psat")
    nist_marker(a, *zip(*pts))
    for prop in ("rhoL", "rhoV"):
        pts = sorted((float(r["T_K"]), float(r["nist_value"])) for r in rows
                     if r["model"] == "PR" and r["property"] == prop)
        nist_marker(b, [v for _, v in pts], [t for t, _ in pts])
    a.set_xlabel("temperature, $T$ (K)")
    a.set_ylabel("saturation pressure,\n$P^{\\mathrm{sat}}$ (MPa)")
    a.set_yticks([0.5, 1, 2, 4, 8])
    a.set_yticklabels(["0.5", "1", "2", "4", "8"])
    a.minorticks_off()
    a.set_ylim(0.45, 9)
    a.grid(color=GRID, lw=0.5)
    h, l = a.get_legend_handles_labels()
    a.legend([h[i] for i in (4, 3, 2, 1, 0)], [l[i] for i in (4, 3, 2, 1, 0)],
             loc="lower right", handlelength=2.6)
    tag(a, "a", x=-0.17)
    b.set_xlabel("saturated density, $\\rho$ (mol L$^{-1}$)")
    b.set_ylabel("temperature, $T$ (K)")
    b.set_xlim(0, 28)
    b.set_ylim(218, 306)
    b.grid(color=GRID, lw=0.5)
    b.text(3.2, 226, "vapor", color=MUTED, fontsize=6.5)
    b.text(20.5, 293, "liquid", color=MUTED, fontsize=6.5)
    tag(b, "b", x=-0.17)
    save(fig, root, "fig6_co2_saturation")


# --------------------------------------------------------------- Fig. 7
def fig_table_route(root):
    rows = [r for r in rcsv(root, "validation/co2_table_route_vs_nist.csv")
            if r["model"].startswith("table, all")]
    fig, axes = plt.subplots(2, 1, figsize=(COL1, 4.4), gridspec_kw=dict(hspace=0.5))
    spec = [("dH", "deviation", 1e3, "$\\Delta h$ deviation (J mol$^{-1}$)", 300),
            ("Cp", "deviation_percent", 1.0, "$C_P$ deviation (%)", 40)]
    for ax, (prop, col, scale, label, lim), letter in zip(axes, spec, "ab"):
        T, P, g = _grid(rows, "table, all 19 temperatures", prop, col)
        g = g * scale
        im = _heat(ax, T, P, g, TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim))
        for i, j in zip(*np.where(np.isnan(g))):
            ax.text(j + 0.5, i + 0.5, "×", ha="center", va="center", fontsize=7, color=INK)
        ax.set_ylabel("pressure, $P$ (MPa)")
        cb = fig.colorbar(im, ax=ax, pad=0.02, aspect=14, extend="both")
        cb.set_label(label)
        cb.outline.set_visible(False)
        tag(ax, letter, x=-0.12)
    save(fig, root, "fig7_table_route")


# --------------------------------------------------------------- Fig. 8
def _budget(root, rel, quantities):
    rows = rcsv(root, rel)
    inputs = [k[6:] for k in rows[0] if k.startswith("share_")]
    data = {r["quantity"]: [float(r["share_" + k]) for k in inputs] for r in rows}
    return inputs, [(q, data[q]) for q in quantities]


def _stack(ax, inputs, data, labels, input_labels):
    y = np.arange(len(data))[::-1]
    left = np.zeros(len(data))
    for k, name in enumerate(inputs):
        vals = np.array([100 * d[1][k] for d in data])
        ax.barh(y, vals, left=left, height=0.62, color=INPUT_COLORS[k],
                edgecolor="white", linewidth=0.8, label=input_labels.get(name, name))
        for yy, l0, v in zip(y, left, vals):
            if v >= 14:
                ax.text(l0 + v / 2, yy, f"{v:.0f}", ha="center", va="center",
                        fontsize=5.8, color="white", fontweight="bold")
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(0, 100)
    ax.set_xlabel("share of output variance (%)")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)


def fig_uncertainty(root):
    fig, (a, b) = plt.subplots(1, 2, figsize=(COL2, 2.2), gridspec_kw=dict(wspace=0.32))
    names = {"T1": "$T_1$", "P1": "$P_1$", "T2": "$T_2$", "P2": "$P_2$", "Tc": "$T_c$",
             "Pc": "$P_c$", "omega": "$\\omega$", "V": "$V$", "Cp": "$C_P$", "w": "$w$",
             "muJT": "$\\mu_{JT}$"}
    inputs, data = _budget(root, "outputs/02_process_co2_pr_uncertainty_budget.csv",
                           ["dH", "dS", "dU", "dG", "dV"])
    _stack(a, inputs, data, ["$\\Delta h$", "$\\Delta s$", "$\\Delta u$", "$\\Delta g$",
                             "$\\Delta v$"], names)
    a.legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0), handlelength=1.0,
             columnspacing=0.9, handletextpad=0.35)
    tag(a, "a", x=-0.07, y=1.02)
    inputs, data = _budget(root, "outputs/05_infer_from_measurements_uncertainty_budget.csv",
                           ["alpha", "kappaT", "Cv", "gamma", "(dT/dP)_S"])
    _stack(b, inputs, data, ["$\\alpha$", "$\\kappa_T$", "$C_V$", "$\\gamma$",
                             "$(\\partial T/\\partial P)_S$"], names)
    b.legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0), handlelength=1.0,
             columnspacing=0.9, handletextpad=0.35)
    tag(b, "b", x=-0.07, y=1.02)
    save(fig, root, "fig8_uncertainty_budget")


# --------------------------------------------------------------- TOC graphic
def fig_toc(root):
    rows = [r for r in _eos_rows(root) if float(r["T_K"]) == 320.0 and r["property"] == "rho"]
    fig = plt.figure(figsize=(3.33, 1.87))
    ax = fig.add_axes([0.60, 0.24, 0.37, 0.62])
    for kind in ("vdW", "RK", "SRK", "PR"):
        pts = sorted((float(r["P_MPa"]), float(r["model_value"])) for r in rows
                     if r["model"] == kind)
        ax.plot(*zip(*pts), **STYLE[kind], lw=1.0)
    pts = sorted((float(r["P_MPa"]), float(r["nist_value"])) for r in rows if r["model"] == "PR")
    nist_marker(ax, *zip(*pts))
    ax.set_xlabel("$P$ (MPa)", labelpad=1)
    ax.set_ylabel("$\\rho$ (mol L$^{-1}$)", labelpad=1)
    ax.set_title("CO$_2$, 320 K: EoS vs NIST", fontsize=6.5, pad=3)
    fig.text(0.03, 0.80, "any derivative", fontsize=7, color=MUTED)
    fig.text(0.03, 0.60, r"$\left(\frac{\partial T}{\partial P}\right)_{\!H}$", fontsize=13)
    fig.text(0.245, 0.615, r"$=\ \frac{V(T\alpha-1)}{C_P}$", fontsize=11.5)
    fig.text(0.03, 0.36, "→ what data it needs", fontsize=7, color=MUTED)
    fig.text(0.03, 0.22, "→ value ± uncertainty", fontsize=7, color=MUTED)
    fig.text(0.03, 0.08, "→ or what is missing", fontsize=7, color=MUTED)
    fig.savefig(os.path.join(root, "toc_graphic.pdf"), dpi=600)
    plt.close(fig)
    print("  wrote toc_graphic.pdf")


def make_all(root):
    print("[7] figures")
    fig_verification(root)
    fig_isotherms(root)
    fig_density_map(root)
    fig_isobar(root)
    fig_saturation(root)
    fig_table_route(root)
    fig_uncertainty(root)
    fig_toc(root)
