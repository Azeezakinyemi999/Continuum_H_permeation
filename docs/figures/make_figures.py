#!/usr/bin/env python3
"""Regenerate every figure in the documentation, from the Example Wall.

No figure here reports a real run. Each one is either computed from
:mod:`example_wall` -- a small fictitious system (see ``docs/CHECKPOINT.md`` D5) -- or
is a labelled schematic. A figure therefore keeps its meaning when the active study
changes, which a figure drawn from a real study would not.

Colours and linestyles come from the vendored production palette in
:mod:`example_wall`, which :mod:`check_equations` holds in sync with
``model_config.COLORS``. Limit-check figures reuse the notebooks' panel (D) idiom:
*this level collapses to the one below.*

Run from the repository root::

    python docs/figures/make_figures.py            # all of them
    python docs/figures/make_figures.py 4 11 13    # only these
"""

from __future__ import annotations

import math
import os
import pathlib
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import example_wall as ex  # noqa: E402

os.environ.setdefault("PERMEATION_STUDY", "Guo_etal_2025/study_config.py")

from calculations import oxide_permeation as ox  # noqa: E402
from calculations import permeability as perm  # noqa: E402
from calculations import permeation_calc as pc  # noqa: E402
from calculations import surface_kinetics as sk  # noqa: E402
from calculations.defective_metal import (  # noqa: E402
    calculate_effective_diffusivity_trapping,
    calculate_gb_enhanced_diffusivity,
)
from calculations.interface_solver import calculate_oxide_metal_system  # noqa: E402
from calculations.parallel_oxide_defect_paths import (  # noqa: E402
    calculate_parallel_path_flux,
)
from calculations.utils import arrhenius  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parent
S = ex.DOC_STYLE
C = ex.COLORS
CS = ex.CURVE_STYLES

#: Vector output. These are line drawings, so SVG stays small and sharp at any zoom.
EXT = "svg"

#: Neutral ink for schematics, annotation and reference rules. Deliberately not one of
#: the level hues -- those carry meaning and must not be spent on furniture.
INK = "#222222"
MUTED = "#6b6b6b"
RULE = "#999999"
HILITE = "#fff3cd"      # callout wash
WARN = "#b3261e"        # a trap or a wrong answer

plt.rcParams.update(
    {
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "axes.linewidth": 1.1,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "legend.framealpha": 0.95,
        "legend.edgecolor": "0.55",
        "svg.fonttype": "none",
    }
)

_captions: dict[int, tuple[str, str]] = {}


def _fig(n, slug, caption, figsize=None, ncols=1):
    """Open a figure and register its caption."""
    _captions[n] = (slug, caption)
    fig, axes = plt.subplots(
        1, ncols, figsize=figsize or (S["figsize"][0], S["figsize"][1])
    )
    return fig, axes


def _save(fig, n, slug):
    path = OUT / f"fig-{n:02d}-{slug}.{EXT}"
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  fig {n:02d}  {path.name}")
    return path


def _style_axes(ax, xlabel=None, ylabel=None, title=None):
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=S["fontsize_axis"])
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=S["fontsize_axis"])
    if title:
        ax.set_title(title, fontsize=S["fontsize_title"])
    ax.tick_params(labelsize=S["fontsize_tick"])
    ax.grid(alpha=S["grid_alpha"], linewidth=0.7)
    return ax


def _box(ax, x, y, w, h, label, facecolor="white", edgecolor=INK, fontsize=None,
         textcolor=INK, lw=1.4, ls="-"):
    ax.add_patch(
        mpatches.FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            facecolor=facecolor, edgecolor=edgecolor, linewidth=lw, linestyle=ls,
        )
    )
    ax.text(
        x + w / 2, y + h / 2, label,
        ha="center", va="center", color=textcolor,
        fontsize=fontsize or S["fontsize_annotation"], fontweight="bold", zorder=5,
    )


def _arrow(ax, xy_from, xy_to, color=INK, lw=1.6, ls="-", style="-|>"):
    ax.annotate(
        "", xy=xy_to, xytext=xy_from,
        arrowprops=dict(arrowstyle=style, color=color, lw=lw, linestyle=ls,
                        shrinkA=2, shrinkB=2),
    )


def _blank(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    return ax


# =============================================================================
# PART I -- OVERVIEW
# =============================================================================

def fig01():
    """The wall end to end, with every symbol drawn where it acts."""
    fig, ax = _fig(
        1, "the-wall",
        "The Example Wall in cross-section. Read left to right: hydrogen arrives as "
        "H2, dissociates (or does not) at a surface, crosses the oxide, crosses the "
        "metal, and leaves. Each symbol is drawn at the place it applies; the "
        "sorption law changes at the oxide/metal interface, which is the single most "
        "consequential fact about this model.",
        figsize=(8.2, 4.2),
    )
    _blank(ax)

    # Layers, drawn contiguous: the wall has no gaps in it.
    top, bot, h = 0.74, 0.40, 0.34
    _box(ax, 0.03, bot, 0.23, h, "gas\nH$_2$", facecolor="#f2f6fa", fontsize=10)
    _box(ax, 0.26, bot, 0.18, h, "oxide\n$L_{ox}$", facecolor="#fdf2ec",
         edgecolor=C["L2a"], fontsize=10)
    _box(ax, 0.44, bot, 0.36, h, "metal\n$L_m$", facecolor="#fdf0f0",
         edgecolor=C["L5"], fontsize=10)
    _box(ax, 0.80, bot, 0.17, h, "sweep\n$P_{down}$", facecolor="#f2f6fa", fontsize=10)

    # The two interfaces that matter
    ax.plot([0.26, 0.26], [bot, top], color=C["L2a"], lw=3.5, zorder=6)
    ax.plot([0.44, 0.44], [bot, top], color=C["L2b"], lw=3.5, zorder=6)

    ax.annotate("gas/oxide surface\ndissociation lives here (L6)",
                xy=(0.26, top), xytext=(0.20, 0.93),
                fontsize=8.4, color=C["L2a"], fontweight="bold", ha="center",
                arrowprops=dict(arrowstyle="-|>", color=C["L2a"], lw=1.3))
    ax.annotate("oxide/metal interface\n$P_{int}$ is solved here",
                xy=(0.44, top), xytext=(0.66, 0.93),
                fontsize=8.4, color=C["L2b"], fontweight="bold", ha="center",
                arrowprops=dict(arrowstyle="-|>", color=C["L2b"], lw=1.3))

    _arrow(ax, (0.045, 0.80), (0.21, 0.80), color=MUTED, lw=2)
    ax.text(0.128, 0.825, "$J$", ha="center", fontsize=11, color=MUTED,
            fontweight="bold")

    # Sorption laws, drawn under the layer each governs, well clear of the callout.
    ax.text(0.35, 0.355, "Henry\n$C=K_{ox}P$", ha="center", va="top",
            fontsize=8.6, color=C["L2a"], fontweight="bold")
    ax.text(0.62, 0.355, "Sieverts\n$C=K_s\\sqrt{P}$", ha="center", va="top",
            fontsize=8.6, color=C["L5"], fontweight="bold")
    ax.annotate("", xy=(0.44, 0.135), xytext=(0.44, 0.385),
                arrowprops=dict(arrowstyle="-", color=RULE, lw=1.0, linestyle=":"))
    ax.text(0.44, 0.075, "the law changes HERE", ha="center", fontsize=9,
            color=WARN, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.32", fc=HILITE, ec=WARN, lw=1.1))

    # Pressures, sitting on the boundary each belongs to
    ax.text(0.045, 0.755, "$P_{up}$", fontsize=10, color=INK, fontweight="bold")
    ax.text(0.448, 0.755, "$P_{int}$", fontsize=10, color=C["L2b"], fontweight="bold")
    ax.text(0.815, 0.755, "$P_{down}$", fontsize=10, color=INK, fontweight="bold")

    ax.set_title("Figure 1 — where each quantity lives", fontsize=S["fontsize_title"])
    return _save(fig, 1, "the-wall")


def fig02():
    """The level hierarchy, drawn as a fork rather than a ladder."""
    fig, ax = _fig(
        2, "level-fork",
        "The levels are cumulative, but they are not one ladder. Levels 1-5 share "
        "Model 1 (dissociation at the oxide/metal interface, Henry oxide). The L6 "
        "family is a different branch under Model 2 (dissociation at the gas surface, "
        "Sieverts oxide). A result does not carry across the fork: L5 and L5L6 "
        "describe the same wall with different physics and are not comparable.",
        figsize=(8.2, 4.6),
    )
    _blank(ax)

    chain = [
        ("L1\nbulk metal", C["L1"]),
        ("L2\n+ oxide", C["L2a"]),
        ("L3\n+ oxide defects", C["L3"]),
        ("L4\n+ microstructure", C["L4_both"]),
        ("L5\nfull wall", C["L5"]),
    ]
    w, h, gap = 0.165, 0.20, 0.035
    x = 0.02
    for i, (label, color) in enumerate(chain):
        _box(ax, x, 0.56, w, h, label, edgecolor=color, fontsize=9)
        if i:
            _arrow(ax, (x - gap, 0.66), (x, 0.66), color=MUTED, lw=1.4)
        x += w + gap

    # The L6 branch
    x6 = 0.02 + 2 * (w + gap)
    for i, (label, color) in enumerate(
        [("L3L6", C["L3"]), ("L4L6", C["L4_both"]), ("L5L6", C["L5"])]
    ):
        bx = x6 + i * (w + gap)
        _box(ax, bx, 0.14, w, h, label, edgecolor=color, ls="--", fontsize=9)
        if i:
            _arrow(ax, (bx - gap, 0.24), (bx, 0.24), color=MUTED, lw=1.4, ls="--")

    _arrow(ax, (x6 + w / 2, 0.55), (x6 + w / 2, 0.35), color=WARN, lw=1.8, ls="--")
    ax.text(x6 + w / 2 + 0.015, 0.45, "add surface kinetics\n= change the oxide's law",
            fontsize=8.5, color=WARN, fontweight="bold", va="center")

    ax.text(0.02, 0.80, "Model 1 — Henry oxide, dissociation at the interface",
            fontsize=9.5, color=INK, fontweight="bold")
    ax.text(0.02, 0.055, "Model 2 — Sieverts oxide, dissociation at the gas surface",
            fontsize=9.5, color=INK, fontweight="bold")

    ax.annotate("", xy=(0.99, 0.46), xytext=(0.01, 0.46),
                arrowprops=dict(arrowstyle="-", color=RULE, lw=1.2, linestyle=(0, (6, 4))))
    ax.text(0.995, 0.47, "the fork", fontsize=8.5, color=RULE, ha="right",
            fontweight="bold")

    ax.set_title("Figure 2 — two branches, not one ladder", fontsize=S["fontsize_title"])
    return _save(fig, 2, "level-fork")


def fig03():
    """The resistance network: what is in series and what is in parallel."""
    fig, ax = _fig(
        3, "resistance-network",
        "The wall as resistances. Surface, oxide and metal are in SERIES, so their "
        "resistances add. Defect paths are in PARALLEL with the intact oxide, and "
        "because flux is extensive those paths combine as an area-weighted sum of "
        "fluxes, not as the reciprocal sum used for parallel electrical resistors.",
        figsize=(8.2, 3.8),
    )
    _blank(ax)

    y = 0.62
    _box(ax, 0.04, y, 0.17, 0.17, "$R_{surf}$\n(L6 only)", edgecolor=C["L2a"],
         ls="--", fontsize=9)
    _box(ax, 0.30, y, 0.17, 0.17, "$R_{ox}$", edgecolor=C["L2a"], fontsize=9)
    _box(ax, 0.56, y, 0.17, 0.17, "$R_{metal}$", edgecolor=C["L5"], fontsize=9)
    for x0, x1 in ((0.21, 0.30), (0.47, 0.56)):
        _arrow(ax, (x0, y + 0.085), (x1, y + 0.085), color=MUTED, lw=1.4)
    ax.text(0.385, 0.86, "in series — resistances ADD",
            ha="center", fontsize=9.5, color=INK, fontweight="bold")

    # Parallel defect path
    _box(ax, 0.30, 0.22, 0.17, 0.17, "$R_{defect}$", edgecolor=C["L3"], fontsize=9)
    ax.plot([0.265, 0.265], [0.305, y + 0.085], color=MUTED, lw=1.3)
    ax.plot([0.265, 0.30], [0.305, 0.305], color=MUTED, lw=1.3)
    ax.plot([0.505, 0.505], [0.305, y + 0.085], color=MUTED, lw=1.3)
    ax.plot([0.47, 0.505], [0.305, 0.305], color=MUTED, lw=1.3)
    ax.text(0.385, 0.15, "in parallel — FLUXES add, area-weighted\n"
                         "$J=f_{i}\\,j_{i}+f_{d}\\,j_{d}$",
            ha="center", fontsize=9, color=C["L3"], fontweight="bold")

    ax.text(0.80, y + 0.085, "$J$", fontsize=12, color=INK, fontweight="bold",
            va="center")
    _arrow(ax, (0.73, y + 0.085), (0.78, y + 0.085), color=MUTED, lw=1.6)

    ax.text(0.99, 0.03,
            "$R_{surf}$ is present only in the L6 family",
            ha="right", fontsize=8.5, color=MUTED, style="italic")
    ax.set_title("Figure 3 — series and parallel, and why the sum is of fluxes",
                 fontsize=S["fontsize_title"])
    return _save(fig, 3, "resistance-network")


# =============================================================================
# PART II -- THEORY
# =============================================================================

def fig04():
    """Henry vs Sieverts: slope 1 against slope one-half."""
    fig, ax = _fig(
        4, "henry-vs-sieverts",
        "The two sorption laws on log-log axes. A Henry (molecular) layer gives flux "
        "linear in pressure, slope 1. A Sieverts (atomic) layer gives flux in the "
        "square root, slope 1/2. The slope is therefore a measurement, not a "
        "convention: it tells you which law governs the step that controls the wall.",
    )
    alpha, beta = ex.permeances()
    P = np.logspace(2, 6, 200)

    ax.loglog(P, alpha * P, color=C["L2a"], lw=S["linewidth"], ls="-",
              label=r"Henry, oxide:  $J=\alpha P$   (slope 1)")
    ax.loglog(P, beta * np.sqrt(P), color=C["L1"], lw=S["linewidth"], ls="--",
              label=r"Sieverts, metal:  $J=\beta\sqrt{P}$   (slope $\frac{1}{2}$)")

    for P0, fn, col, dy in ((1e4, lambda p: alpha * p, C["L2a"], 1.9),
                            (1e4, lambda p: beta * np.sqrt(p), C["L1"], 0.45)):
        ax.annotate("", xy=(P0 * 6, fn(P0) * dy), xytext=(P0, fn(P0) * dy),
                    arrowprops=dict(arrowstyle="-", color=col, lw=1.0, ls=":"))

    ax.text(2.5e4, alpha * 1e4 * 3.2, "slope 1", color=C["L2a"],
            fontsize=S["fontsize_annotation"], fontweight="bold")
    ax.text(2.5e4, beta * np.sqrt(1e4) * 0.3, "slope 1/2", color=C["L1"],
            fontsize=S["fontsize_annotation"], fontweight="bold")

    _style_axes(ax, "upstream pressure $P$  [Pa]", "flux $J$  [mol m$^{-2}$ s$^{-1}$]",
                "Figure 4 — the slope identifies the law")
    ax.legend(fontsize=S["fontsize_legend"], loc="upper left")
    return _save(fig, 4, "henry-vs-sieverts")


def fig05():
    """The K_ox units trap: one symbol, two unit systems, no conversion."""
    fig, ax = _fig(
        5, "kox-units-trap",
        "The same symbol K_ox carries different units in the two models, because it "
        "describes the dissolution of different species. There is no conversion "
        "factor between the two columns: they are not the same quantity in different "
        "units, they are different quantities. calculations.permeability refuses to "
        "combine them, raising UnitMismatchError, and that refusal is the physics.",
        figsize=(8.0, 4.4),
    )
    _blank(ax)

    _box(ax, 0.04, 0.30, 0.38, 0.52, "", facecolor="#fdf6f0", edgecolor=C["L2a"], lw=1.6)
    _box(ax, 0.58, 0.30, 0.38, 0.52, "", facecolor="#f0f4fd", edgecolor=C["L2b"], lw=1.6)

    ax.text(0.23, 0.76, "Model 1  (L1–L5)", ha="center", fontsize=10.5,
            color=C["L2a"], fontweight="bold")
    ax.text(0.77, 0.76, "Model 2  (L6 family)", ha="center", fontsize=10.5,
            color=C["L2b"], fontweight="bold")

    left = [
        "dissociation at the oxide/metal interface",
        "intact H$_2$ crosses the oxide",
        r"$C = K_{ox}\,P$    (Henry)",
        r"$[K_{ox}]$ = mol m$^{-3}$ Pa$^{-1}$",
        r"$[\Phi_{ox}]$ = mol m$^{-1}$ s$^{-1}$ Pa$^{-1}$",
    ]
    right = [
        "dissociation at the gas/oxide surface",
        "atomic H crosses the oxide",
        r"$C = K_{ox}\sqrt{P}$    (Sieverts)",
        r"$[K_{ox}]$ = mol m$^{-3}$ Pa$^{-1/2}$",
        r"$[\Phi_{ox}]$ = mol m$^{-1}$ s$^{-1}$ Pa$^{-1/2}$",
    ]
    for i, (a, b) in enumerate(zip(left, right)):
        yy = 0.68 - i * 0.082
        weight = "bold" if i >= 2 else "normal"
        ax.text(0.23, yy, a, ha="center", fontsize=8.6, color=INK, fontweight=weight)
        ax.text(0.77, yy, b, ha="center", fontsize=8.6, color=INK, fontweight=weight)

    # The crossed-out conversion
    ax.annotate("", xy=(0.56, 0.56), xytext=(0.44, 0.56),
                arrowprops=dict(arrowstyle="<|-|>", color=RULE, lw=1.6))
    ax.plot([0.46, 0.54], [0.52, 0.60], color=WARN, lw=3, solid_capstyle="round")
    ax.plot([0.46, 0.54], [0.60, 0.52], color=WARN, lw=3, solid_capstyle="round")
    ax.text(0.50, 0.44, "no conversion\nexists", ha="center", fontsize=8.8,
            color=WARN, fontweight="bold")

    ax.text(0.50, 0.17,
            "combine_series(metal, oxide)  →  UnitMismatchError",
            ha="center", fontsize=9.2, color=WARN, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.45", fc=HILITE, ec=WARN, lw=1.2))
    ax.text(0.50, 0.055,
            "which is why Model 1 has no single system permeability — only the four tiers",
            ha="center", fontsize=8.4, color=MUTED, style="italic")

    ax.set_title("Figure 5 — one symbol, two quantities",
                 fontsize=S["fontsize_title"])
    return _save(fig, 5, "kox-units-trap")


def fig06():
    """The interface pressure is where two flux curves cross."""
    fig, ax = _fig(
        6, "interface-root",
        "Solving Level 2b. The oxide flux falls as the interface pressure rises; the "
        "metal flux rises. They cross at exactly one interface pressure, and that "
        "crossing is the answer. The curve shapes differ because the laws differ — "
        "the oxide line is straight in P, the metal curve goes as the square root. "
        "The bracket brentq searches is the whole physical range.",
    )
    alpha, beta = ex.permeances()
    P_int = np.linspace(1.0, ex.P_UP * 0.999, 400)
    J_ox = alpha * (ex.P_UP - P_int)
    J_m = beta * np.sqrt(P_int)

    ax.plot(P_int, J_ox, color=C["L2a"], lw=S["linewidth"],
            label=r"oxide:  $\alpha\,(P_{up}-P_{int})$")
    ax.plot(P_int, J_m, color=C["L5"], lw=S["linewidth"],
            label=r"metal:  $\beta\,\sqrt{P_{int}}$")

    P_star = ex.exact_interface_pressure()
    J_star = ex.exact_level2b_flux()
    ax.plot([P_star], [J_star], "o", color=C["L2b"], ms=10, zorder=6,
            markeredgecolor="white", markeredgewidth=1.4)
    ax.annotate(
        f"the solution\n$P_{{int}}$ = {P_star:,.0f} Pa\n$J$ = {J_star:.3g}",
        xy=(P_star, J_star), xytext=(P_star * 4.2, J_star * 1.9),
        fontsize=S["fontsize_annotation"], color=C["L2b"], fontweight="bold",
        arrowprops=dict(arrowstyle="-|>", color=C["L2b"], lw=1.4),
    )
    ax.axvline(P_star, color=C["L2b"], lw=1.0, ls=":")

    ax.axvspan(0, ex.P_UP, color=RULE, alpha=0.07)
    ax.text(ex.P_UP * 0.5, max(J_ox) * 0.93,
            "bracket: $[P_{down},\\,P_{up}]$ — one root, guaranteed",
            ha="center", fontsize=8.6, color=MUTED, fontweight="bold")

    _style_axes(ax, "interface pressure $P_{int}$  [Pa]",
                "flux  [mol m$^{-2}$ s$^{-1}$]",
                "Figure 6 — the interface pressure as a root")
    ax.set_xscale("log")
    ax.set_xlim(10, ex.P_UP)
    ax.legend(fontsize=S["fontsize_legend"], loc="upper left")
    return _save(fig, 6, "interface-root")


def fig07():
    """Resistance fractions: many parameters collapse into one number."""
    fig, axes = _fig(
        7, "resistance-fractions",
        "Where the resistance sits, as the coating thickens. Six transport parameters "
        "collapse into a single fraction, and that fraction is what the regime label "
        "reports. Left: the handover as the oxide thickens. Right: the Example Wall's "
        "own split. Note the handover is gradual — the regime label puts a hard "
        "boundary on a smooth curve.",
        figsize=(8.6, 3.9), ncols=2,
    )
    ax1, ax2 = axes
    alpha0, beta = ex.permeances()

    L_ox = np.logspace(-9, -4, 160)
    fr_ox, fr_m = [], []
    for L in L_ox:
        props = ex.oxide_props()
        props["thickness"] = L
        r = calculate_oxide_metal_system(ex.P_UP, ex.P_DOWN, props, ex.metal_props())
        P_int = r["P_interface"]
        R_ox = ox.calculate_oxide_resistance(ex.D_OX, ex.K_OX, L)
        R_m = ox.calculate_metal_resistance(ex.D_M, ex.K_S_M, ex.L_M, max(P_int, 1e-12))
        fr_ox.append(R_ox / (R_ox + R_m))
        fr_m.append(R_m / (R_ox + R_m))

    ax1.semilogx(L_ox * 1e6, fr_ox, color=C["L2a"], lw=S["linewidth"], label="oxide")
    ax1.semilogx(L_ox * 1e6, fr_m, color=C["L5"], lw=S["linewidth"], label="metal")
    ax1.axhline(0.5, color=RULE, ls="--", lw=1.1)
    ax1.text(L_ox[2] * 1e6, 0.53, "the 0.5 line the regime label uses",
             fontsize=8.0, color=MUTED, fontweight="bold")
    ax1.axvline(ex.L_OX * 1e6, color=C["L2b"], ls=":", lw=1.4)
    ax1.text(ex.L_OX * 1e6 * 1.25, 0.08, "Example\nWall", fontsize=8.2,
             color=C["L2b"], fontweight="bold")
    _style_axes(ax1, "oxide thickness  [µm]", "fraction of total resistance",
                "(A) the handover is smooth")
    ax1.set_ylim(0, 1)
    ax1.legend(fontsize=S["fontsize_legend"], loc="center left")

    r = calculate_oxide_metal_system(ex.P_UP, ex.P_DOWN, ex.oxide_props(), ex.metal_props())
    P_int = r["P_interface"]
    R_ox = ox.calculate_oxide_resistance(ex.D_OX, ex.K_OX, ex.L_OX)
    R_m = ox.calculate_metal_resistance(ex.D_M, ex.K_S_M, ex.L_M, P_int)
    tot = R_ox + R_m
    ax2.barh([0], [R_ox / tot], color=C["L2a"], edgecolor="white", height=0.5, label="oxide")
    ax2.barh([0], [R_m / tot], left=[R_ox / tot], color=C["L5"], edgecolor="white",
             height=0.5, label="metal")
    ax2.text(R_ox / tot / 2, 0, f"{R_ox/tot:.0%}", ha="center", va="center",
             color="white", fontweight="bold", fontsize=11)
    ax2.text(R_ox / tot + R_m / tot / 2, 0, f"{R_m/tot:.0%}", ha="center", va="center",
             color="white", fontweight="bold", fontsize=11)
    ax2.set_yticks([])
    ax2.set_xlim(0, 1)
    _style_axes(ax2, "fraction of total resistance", None,
                "(B) the Example Wall's split")
    ax2.grid(False)
    ax2.legend(fontsize=S["fontsize_legend"], loc="lower center", ncol=2)
    fig.suptitle("Figure 7 — six parameters collapse into one fraction",
                 fontsize=S["fontsize_title"], y=1.03)
    fig.tight_layout()
    return _save(fig, 7, "resistance-fractions")


def fig08():
    """Langmuir coverage, and the (1-theta)^2 term that does the suppressing."""
    fig, axes = _fig(
        8, "langmuir-coverage",
        "Surface coverage and what it costs. Left: coverage rises towards saturation "
        "with pressure. Right: the vacancy factor (1-θ)², which multiplies the "
        "dissociation rate. Once the surface is crowded the wall cannot take hydrogen "
        "in any faster, however much pressure is applied — which is why Level 6 can "
        "make a wall surface-limited that Level 5 called metal-limited.",
        figsize=(8.6, 3.9), ncols=2,
    )
    ax1, ax2 = axes
    theta = np.linspace(1e-4, 0.995, 400)
    ax1.plot(theta, [ex.exact_g_theta(t) ** 2 for t in theta], color=C["L2a"],
             lw=S["linewidth"])
    ax1.set_yscale("log")
    _style_axes(ax1, r"coverage  $\theta$", r"$P_{virtual}=g(\theta)^2$  [Pa]",
                r"(A) coverage sets an effective pressure")
    ax1.axvline(0.5, color=RULE, ls=":", lw=1.1)
    ax1.text(0.52, ax1.get_ylim()[0] * 30, r"$\theta=\frac{1}{2}$", fontsize=8.5,
             color=MUTED, fontweight="bold")

    ax2.plot(theta, (1 - theta) ** 2, color=C["L5"], lw=S["linewidth"])
    ax2.fill_between(theta, 0, (1 - theta) ** 2, color=C["L5"], alpha=0.10)
    _style_axes(ax2, r"coverage  $\theta$", r"vacancy factor  $(1-\theta)^2$",
                "(B) a crowded surface stops accepting")
    ax2.set_ylim(0, 1)
    for t in (0.5, 0.9):
        ax2.plot([t], [(1 - t) ** 2], "o", color=C["L5"], ms=7, zorder=5,
                 markeredgecolor="white")
        ax2.annotate(f"θ={t:g} → {(1-t)**2:.2f}", xy=(t, (1 - t) ** 2),
                     xytext=(t - 0.33, (1 - t) ** 2 + 0.17), fontsize=8.4,
                     color=C["L5"], fontweight="bold",
                     arrowprops=dict(arrowstyle="-|>", color=C["L5"], lw=1.1))
    fig.suptitle("Figure 8 — coverage, and the vacancy factor it costs",
                 fontsize=S["fontsize_title"], y=1.03)
    fig.tight_layout()
    return _save(fig, 8, "langmuir-coverage")


# =============================================================================
# PART III -- THE LEVELS
# =============================================================================

def fig09():
    """Level 1: the square-root law and the Arrhenius slope."""
    fig, axes = _fig(
        9, "level1-sqrt-arrhenius",
        "Level 1's two signatures. Left: flux is straight against √P through the "
        "origin — the Sieverts signature. Right: on an Arrhenius plot the slope "
        "recovers Q = E_D + H_s exactly, which is how a measured permeability is "
        "decomposed. Both panels are limit checks: if either fails, the metal block "
        "is wrong.",
        figsize=(8.6, 3.9), ncols=2,
    )
    ax1, ax2 = axes
    P = np.linspace(0, 4e5, 120)
    _, beta = ex.permeances()
    ax1.plot(np.sqrt(P), beta * np.sqrt(P), **{k: v for k, v in CS["L1"].items()
                                               if k != "marker"},
             label=r"$J=\beta\sqrt{P}$")
    ax1.plot(np.sqrt(P[::20]), beta * np.sqrt(P[::20]), "o", color=C["L1"], ms=5)
    _style_axes(ax1, r"$\sqrt{P}$  [Pa$^{1/2}$]", "flux  [mol m$^{-2}$ s$^{-1}$]",
                "(A) straight against √P, through zero")
    ax1.legend(fontsize=S["fontsize_legend"])

    T = np.linspace(600, 1000, 120)
    # arrhenius() validates its temperature with a bare `if temperature <= 0`, so it
    # takes scalars only; mapping is deliberate, not a stylistic choice.
    D = np.array([arrhenius(ex.D_M, ex.E_D_M, t, ex.T_K) for t in T])
    Ks = np.array([arrhenius(ex.K_S_M, ex.H_S_M, t, ex.T_K) for t in T])
    J = (D * Ks / ex.L_M) * math.sqrt(ex.P_UP)
    ax2.semilogy(1000.0 / T, J, color=C["L1"], lw=S["linewidth"])
    slope = np.polyfit(1.0 / T, np.log(J), 1)[0]
    Q_fit = -slope * ex.R_GAS
    ax2.text(
        0.05, 0.08,
        f"fitted $Q$ = {Q_fit/1e3:.1f} kJ/mol\nstated $E_D+H_s$ = "
        f"{(ex.E_D_M+ex.H_S_M)/1e3:.1f} kJ/mol",
        transform=ax2.transAxes, fontsize=8.6, color=C["L1"], fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.4", fc=HILITE, ec=RULE, lw=0.9),
    )
    _style_axes(ax2, "1000/T  [K$^{-1}$]", "flux  [mol m$^{-2}$ s$^{-1}$]",
                "(B) the slope recovers $Q=E_D+H_s$")
    fig.suptitle("Figure 9 — Level 1's two signatures",
                 fontsize=S["fontsize_title"], y=1.03)
    fig.tight_layout()
    return _save(fig, 9, "level1-sqrt-arrhenius")


def fig10():
    """Level 2 limit check, in the production panel (D) idiom."""
    fig, ax = _fig(
        10, "level2-limit",
        "Limit check in the notebooks' panel (D) idiom. As the oxide is thinned the "
        "series wall must climb to the bare-metal flux, and it does. A solver that "
        "failed to recover the level below would show here as a curve that flattens "
        "at the wrong value — which is the point of drawing it rather than asserting it.",
    )
    L_ox = np.logspace(-12, -5, 110)
    J = []
    for L in L_ox:
        props = ex.oxide_props()
        props["thickness"] = L
        J.append(calculate_oxide_metal_system(ex.P_UP, ex.P_DOWN, props,
                                              ex.metal_props())["flux"])
    ax.loglog(L_ox * 1e6, J, color=C["L2b"], lw=S["linewidth"], label="Level 2b, series")
    ax.axhline(ex.exact_level1_flux(), color=C["L1"], ls="--", lw=1.8,
               label="Level 1, bare metal")
    ax.axvline(ex.L_OX * 1e6, color=RULE, ls=":", lw=1.3)
    ax.text(ex.L_OX * 1e6 * 1.3, min(J) * 3, "Example Wall", fontsize=8.4,
            color=MUTED, fontweight="bold", rotation=90)
    ax.annotate(r"$L_{ox}\rightarrow 0$  recovers L1",
                xy=(L_ox[3] * 1e6, ex.exact_level1_flux()),
                xytext=(L_ox[3] * 1e6 * 3, ex.exact_level1_flux() * 0.06),
                fontsize=S["fontsize_annotation"], color=C["L1"], fontweight="bold",
                arrowprops=dict(arrowstyle="-|>", color=C["L1"], lw=1.3))
    _style_axes(ax, "oxide thickness  [µm]", "flux  [mol m$^{-2}$ s$^{-1}$]",
                "Figure 10 — (D) limit check: L2b → L1")
    ax.legend(fontsize=S["fontsize_legend"], loc="lower left")
    return _save(fig, 10, "level2-limit")


def fig11():
    """Level 3: the defect takeover threshold, and that it is not universal."""
    fig, ax = _fig(
        11, "level3-threshold",
        "When do defects take over? The crossing — where the defect path carries as "
        "much flux as the whole intact area — moves by orders of magnitude with the "
        "quality of the coating. A threshold quoted without its oxide is meaningless. "
        "The better the oxide, the smaller the defect fraction that ruins it.",
    )
    f = np.logspace(-6, -0.5, 90)
    for mult, label, color, ls in (
        (1.0, "Example Wall oxide", C["L3"], "-"),
        (1e-3, "oxide 1000× less permeable", C["L2a"], "--"),
    ):
        props = ex.oxide_props()
        props["D_ox"] = ex.D_OX * mult
        ratio = []
        for fd in f:
            r = calculate_parallel_path_flux(
                ex.P_UP, ex.P_DOWN, props, ex.metal_props(), ex.defect_params(fd)
            )
            ratio.append(r["flux_defect_contribution"] / r["flux_intact_contribution"])
        ratio = np.asarray(ratio)
        ax.loglog(f, ratio, color=color, lw=S["linewidth"], ls=ls, label=label)
        cross = np.interp(0.0, np.log10(ratio), np.log10(f))
        ax.plot([10**cross], [1.0], "o", color=color, ms=9, zorder=6,
                markeredgecolor="white", markeredgewidth=1.3)
        ax.annotate(f"takeover at $f$ ≈ {10**cross:.1e}",
                    xy=(10**cross, 1.0), xytext=(10**cross * 0.02, 10 ** (1.1)),
                    fontsize=8.4, color=color, fontweight="bold",
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=1.2))

    ax.axhline(1.0, color=RULE, ls=":", lw=1.4)
    ax.text(f[1], 1.5, "defect path = intact path", fontsize=8.4, color=MUTED,
            fontweight="bold")
    _style_axes(ax, "defect area fraction  $f_d$",
                "defect flux ÷ intact flux",
                "Figure 11 — the takeover threshold is not a constant")
    ax.legend(fontsize=S["fontsize_legend"], loc="lower right")
    return _save(fig, 11, "level3-threshold")


def fig12():
    """Level 4: two effects, only one of which matters at structural grain sizes."""
    fig, axes = _fig(
        12, "level4-competition",
        "Level 4 adds two effects that pull opposite ways, but they are not of "
        "comparable size. Left: grain-boundary enhancement stays within a few per "
        "cent of 1 for any grain size a structural alloy actually has, and only "
        "becomes a real multiplier in the ultrafine range. Right: trapping changes "
        "the diffusivity by orders of magnitude over the same window. For an "
        "engineering wall, Level 4 is a trapping model with a grain-boundary term "
        "that rounds away — which is worth knowing before spending effort measuring "
        "grain size.",
        figsize=(8.6, 3.9), ncols=2,
    )
    ax1, ax2 = axes
    C_L = ex.K_S_M * math.sqrt(ex.P_UP)

    # (A) enhancement against grain size
    d = np.logspace(-8, -4, 120)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        enh = np.array([
            calculate_gb_enhanced_diffusivity(
                D_bulk=ex.D_M, temperature=ex.T_K, grain_size=dd,
                gb_thickness=ex.GB_THICKNESS,
            )["enhancement_ratio"]
            for dd in d
        ])
    ax1.semilogx(d * 1e9, enh, color=C["L4_gb"], lw=S["linewidth"])
    ax1.axhline(1.0, color=RULE, ls="--", lw=1.3)
    ax1.axvspan(1e3, 1e5, color=C["L4_gb"], alpha=0.08)
    ax1.text(6e3, enh.max() * 0.80, "grain sizes a structural\nalloy actually has",
             fontsize=8.0, color=C["L4_gb"], fontweight="bold", ha="center")
    ax1.axvline(ex.GRAIN_SIZE * 1e9, color=C["L4_both"], ls=":", lw=1.5)
    enh_ex = float(np.interp(ex.GRAIN_SIZE * 1e9, d * 1e9, enh))
    ax1.annotate(f"Example Wall\n{ex.GRAIN_SIZE*1e6:g} µm → {enh_ex:.4f}×",
                 xy=(ex.GRAIN_SIZE * 1e9, enh_ex),
                 xytext=(3e2, 1.0 + (enh.max() - 1) * 0.55),
                 fontsize=8.0, color=C["L4_both"], fontweight="bold",
                 arrowprops=dict(arrowstyle="-|>", color=C["L4_both"], lw=1.2))
    _style_axes(ax1, "grain size  [nm]", r"$D_{eff}/D_{lattice}$",
                "(A) grain boundaries: almost nothing")

    # (B) trapping against temperature, with the GB term drawn for scale
    T = np.linspace(500, 820, 140)     # within the GB correlation's stated data range
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        trap_ratio = np.array([
            calculate_effective_diffusivity_trapping(
                D_lattice=ex.D_M, temperature=t, trap_list=ex.trap_list(),
                lattice_concentration=C_L, lattice_density=ex.N_L,
            )["D_eff"] / ex.D_M
            for t in T
        ])
        gb_ratio = np.array([
            calculate_gb_enhanced_diffusivity(
                D_bulk=ex.D_M, temperature=t, grain_size=ex.GRAIN_SIZE,
                gb_thickness=ex.GB_THICKNESS,
            )["enhancement_ratio"]
            for t in T
        ])

    ax2.semilogy(T, trap_ratio, color=C["L4_trap"], lw=S["linewidth"],
                 label="traps alone  (slows down)")
    ax2.semilogy(T, gb_ratio, color=C["L4_gb"], lw=S["linewidth"],
                 label=f"grain boundaries at {ex.GRAIN_SIZE*1e6:g} µm")
    ax2.semilogy(T, trap_ratio * gb_ratio, color=C["L4_both"],
                 lw=S["linewidth"] + 0.6, ls="--",
                 label="both, as the code composes them")
    ax2.axhline(1.0, color=RULE, ls="--", lw=1.3)
    ax2.set_ylim(top=3.0)
    ax2.text(0.03, 0.93, "no net effect", transform=ax2.transAxes,
             fontsize=8.0, color=MUTED, fontweight="bold")
    ax2.annotate("the trapping curve and the combined curve\n"
                 "are indistinguishable at this grain size",
                 xy=(T[70], trap_ratio[70]), xytext=(T[18], 3e-3),
                 fontsize=7.8, color=C["L4_both"], fontweight="bold",
                 arrowprops=dict(arrowstyle="-|>", color=C["L4_both"], lw=1.2))
    _style_axes(ax2, "temperature  [K]", r"$D_{eff}/D_{lattice}$",
                "(B) trapping: orders of magnitude")
    ax2.legend(fontsize=7.8, loc="lower right")

    fig.suptitle("Figure 12 — two effects, one that matters",
                 fontsize=S["fontsize_title"], y=1.03)
    fig.tight_layout()
    return _save(fig, 12, "level4-competition")


def fig13():
    """Level 5: the corrections do not multiply."""
    fig, ax = _fig(
        13, "level5-not-multiplicative",
        "A trap that produces a plausible wrong answer. Having measured the oxide "
        "defect correction and the metal microstructure correction separately, it is "
        "natural to multiply them. The true Level 5 result is not that product, "
        "because the two corrections both act through the interface pressure and so "
        "are coupled. The grey bar is the wrong answer; the red bar is the right one.",
        figsize=(7.4, 4.2),
    )
    base = calculate_oxide_metal_system(
        ex.P_UP, ex.P_DOWN, ex.oxide_props(), ex.metal_props()
    )["flux"]

    f_big = 0.05
    J_defects = calculate_parallel_path_flux(
        ex.P_UP, ex.P_DOWN, ex.oxide_props(), ex.metal_props(),
        ex.defect_params(f_big),
    )["flux_total"]

    metal_slow = ex.metal_props()
    metal_slow["D_metal"] = ex.exact_D_eff_trapping()
    J_metal = calculate_oxide_metal_system(
        ex.P_UP, ex.P_DOWN, ex.oxide_props(), metal_slow
    )["flux"]

    J_naive = base * (J_defects / base) * (J_metal / base)
    J_true = calculate_parallel_path_flux(
        ex.P_UP, ex.P_DOWN, ex.oxide_props(), metal_slow, ex.defect_params(f_big)
    )["flux_total"]

    labels = ["L2b\nbaseline", "oxide defects\nonly", "slow metal\nonly",
              "naive\nproduct", "true L5"]
    values = [base, J_defects, J_metal, J_naive, J_true]
    colors = [C["L2b"], C["L3"], C["L4_both"], RULE, C["L5"]]
    bars = ax.bar(labels, values, color=colors, edgecolor="white", width=0.68)
    bars[3].set_hatch("//")

    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.03, f"{v:.3g}",
                ha="center", fontsize=8.4, fontweight="bold", color=INK)

    err = (J_naive - J_true) / J_true
    ax.annotate(
        f"the product is off by {err:+.1%}",
        xy=(3.5, max(J_naive, J_true) * 1.14), ha="center",
        fontsize=9.4, color=WARN, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.4", fc=HILITE, ec=WARN, lw=1.2),
    )
    ax.plot([3, 4], [J_naive, J_true], color=WARN, lw=1.4, ls=":", marker="o", ms=4)

    _style_axes(ax, None, "flux  [mol m$^{-2}$ s$^{-1}$]",
                "Figure 13 — corrections do not multiply")
    ax.set_ylim(0, max(values) * 1.32)
    ax.grid(axis="x", alpha=0)
    return _save(fig, 13, "level5-not-multiplicative")


def fig14():
    """Level 6: the surface suppresses the interface pressure, not the flux directly."""
    fig, axes = _fig(
        14, "level6-pressure-not-flux",
        "The single most consequential trap in this model. A slow surface does not "
        "throttle the flux by some efficiency factor; it starves the layers behind it "
        "by holding the interface pressure down. Left: the interface pressure "
        "collapses as dissociation slows. Right: the flux follows, and the wall "
        "changes which step limits it — so a conclusion drawn at Level 5 can reverse "
        "under the L6 family on the very same wall.",
        figsize=(8.6, 3.9), ncols=2,
    )
    ax1, ax2 = axes
    k = np.logspace(-13, -4, 90)
    P_int, J, limiting = [], [], []
    for kd in k:
        s = sk.solve_steady_state_flux_L1L6(
            ex.P_UP, ex.P_DOWN, ex.L_M, kd, ex.K_EQ, ex.D_M, ex.K_S_M
        )
        P_int.append(s["P_int"])
        J.append(s["J_ss"])
        limiting.append(s["rate_limiting"])

    ax1.loglog(k, P_int, color=C["L2a"], lw=S["linewidth"])
    ax1.axhline(ex.P_UP, color=C["L1"], ls="--", lw=1.6)
    ax1.text(k[1], ex.P_UP * 1.3, r"$P_{up}$ — the L1 assumption", fontsize=8.2,
             color=C["L1"], fontweight="bold")
    _style_axes(ax1, r"$k_{diss}$  [mol m$^{-2}$ s$^{-1}$ Pa$^{-1}$]",
                r"$P_{int}$  [Pa]", "(A) the surface starves the wall")

    ax2.loglog(k, J, color=C["L5"], lw=S["linewidth"], label="L1L6")
    ax2.axhline(ex.exact_level1_flux(), color=C["L1"], ls="--", lw=1.6,
                label="L1 (instant dissociation)")
    switch = None
    for i in range(1, len(limiting)):
        if limiting[i] != limiting[i - 1]:
            switch = k[i]
            break
    if switch:
        for a in (ax1, ax2):
            a.axvline(switch, color=WARN, ls=":", lw=1.6)
        ax2.annotate("surface-limited\n← | →\nmetal-limited",
                     xy=(switch, ex.exact_level1_flux() * 1e-2),
                     xytext=(switch * 60, ex.exact_level1_flux() * 2e-3),
                     fontsize=8.0, color=WARN, fontweight="bold", ha="left",
                     arrowprops=dict(arrowstyle="-|>", color=WARN, lw=1.2))
    _style_axes(ax2, r"$k_{diss}$  [mol m$^{-2}$ s$^{-1}$ Pa$^{-1}$]",
                "flux  [mol m$^{-2}$ s$^{-1}$]", "(B) and the flux follows")
    ax2.legend(fontsize=S["fontsize_legend"], loc="lower right")
    fig.suptitle("Figure 14 — the surface suppresses pressure, not flux",
                 fontsize=S["fontsize_title"], y=1.03)
    fig.tight_layout()
    return _save(fig, 14, "level6-pressure-not-flux")


# =============================================================================
# PART III-b -- PERMEABILITY
# =============================================================================

def fig15():
    """The four permeability tiers, and what each needs to be reproduced."""
    fig, ax = _fig(
        15, "permeability-tiers",
        "Why there is no single number called 'the permeability of this wall'. Four "
        "tiers, distinguished by what each needs before it can be reproduced. Tiers "
        "1-3 follow from parameters alone; tier 4 is a property of the wall AND the "
        "experiment, so it is always reported with its operating point. Under Model 1 "
        "tier 3 does not exist at all, because the two layers' units do not combine.",
        figsize=(8.4, 3.9),
    )
    _blank(ax)

    rows = [
        ("1  intrinsic", "nothing", r"$D_m K_s$,  $D_{ox}K_{ox}$", C["L1"], True),
        ("2  effective layer", "nothing", r"$D_{eff}K_s$,  $\Phi_{ox,eff}$", C["L4_both"], True),
        ("3  stack", r"$L$", r"$\Phi_{transport}$ — Model 2 only", C["L2b"], False),
        ("4  apparent", r"$L$ and $P$", r"$\Phi_{app}$ from the solved flux", C["L5"], False),
    ]
    ax.text(0.03, 0.88, "tier", fontsize=9.5, fontweight="bold", color=INK)
    ax.text(0.30, 0.88, "depends on", fontsize=9.5, fontweight="bold", color=INK)
    ax.text(0.55, 0.88, "quantity", fontsize=9.5, fontweight="bold", color=INK)
    ax.plot([0.02, 0.98], [0.85, 0.85], color=INK, lw=1.3)

    for i, (tier, dep, qty, color, reproducible) in enumerate(rows):
        y = 0.70 - i * 0.165
        ax.add_patch(mpatches.Rectangle((0.02, y - 0.045), 0.96, 0.125,
                                        facecolor=color, alpha=0.08, edgecolor="none"))
        ax.text(0.03, y, tier, fontsize=9.2, fontweight="bold", color=color)
        ax.text(0.30, y, dep, fontsize=9.0, color=INK)
        ax.text(0.55, y, qty, fontsize=9.0, color=INK)
        ax.text(0.95, y, "reproducible" if reproducible else "needs its operating point",
                fontsize=7.6, color=MUTED, ha="right", style="italic")

    ax.plot([0.02, 0.98], [0.70 - 1.5 * 0.165 - 0.045, 0.70 - 1.5 * 0.165 - 0.045],
            color=RULE, lw=1.1, ls="--")
    ax.text(0.50, 0.035,
            "above the dashed line: parameters alone.   below: wall + experiment.",
            ha="center", fontsize=8.2, color=MUTED, fontweight="bold")
    ax.set_title("Figure 15 — four tiers, not one number", fontsize=S["fontsize_title"])
    return _save(fig, 15, "permeability-tiers")


def fig16():
    """The pressure exponent as a model-free diagnostic."""
    fig, ax = _fig(
        16, "pressure-exponent",
        "The pressure exponent n = dlnJ/dlnP read off a wall, with no model "
        "internals. Pure Sieverts transport gives 1/2, a Henry layer gives 1, and a "
        "mixed wall sits between — moving as pressure shifts which step controls. "
        "This is an independent check on the regime label: if they disagree, one of "
        "them is wrong.",
    )
    P = np.logspace(2, 7, 60)
    alpha, beta = ex.permeances()

    n_series = [
        perm.pressure_exponent(
            lambda p: calculate_oxide_metal_system(
                p, ex.P_DOWN, ex.oxide_props(), ex.metal_props()
            )["flux"],
            P0,
        )
        for P0 in P
    ]
    ax.semilogx(P, n_series, color=C["L2b"], lw=S["linewidth"] + 0.4,
                label="the Example Wall (oxide + metal)")
    ax.axhline(0.5, color=C["L1"], ls="--", lw=1.6,
               label=r"pure Sieverts step,  $n=\frac{1}{2}$")
    ax.axhline(1.0, color=C["L2a"], ls="--", lw=1.6, label=r"pure Henry layer,  $n=1$")

    ax.fill_between(P, 0.5, 1.0, color=RULE, alpha=0.06)
    ax.text(P[2], 0.75, "mixed control", fontsize=8.6, color=MUTED, fontweight="bold")
    ax.annotate("oxide-controlled here\n(this wall's coating is the bottleneck)",
                xy=(1e3, float(np.interp(1e3, P, n_series))),
                xytext=(3e3, 0.62), fontsize=8.3, color=C["L2b"], fontweight="bold",
                arrowprops=dict(arrowstyle="-|>", color=C["L2b"], lw=1.2))

    _style_axes(ax, "upstream pressure  [Pa]", r"pressure exponent  $n$",
                "Figure 16 — n as a model-free diagnostic")
    ax.set_ylim(0.35, 1.12)
    ax.legend(fontsize=S["fontsize_legend"], loc="lower left")
    return _save(fig, 16, "pressure-exponent")


# =============================================================================
# PART IV -- THE SENSITIVITY WORKFLOW
# =============================================================================

def fig17():
    """The sensitivity pipeline, stage by stage."""
    fig, ax = _fig(
        17, "sa-pipeline",
        "The regime-stratified sensitivity workflow. Read top to bottom. The two "
        "shaded stages are the ones that make this workflow different from an "
        "ordinary global sensitivity analysis: the sample is partitioned by which "
        "physical process limits the flux, and each regime is then re-sampled on its "
        "own so that every regime has enough draws to support an index.",
        figsize=(7.2, 5.4),
    )
    _blank(ax)

    stages = [
        ("1 · validate the configuration", "check_against_config", False),
        ("2 · global LHS scan", "run_global_lhs_scan", False),
        ("3 · partition by regime", "partition_by_regime", True),
        ("4 · size the draws per regime", "size_draws_for_target", True),
        ("5 · targeted re-scan", "run_targeted_regime_scans", True),
        ("6 · PAWN and δ per regime", "givendata_sensitivity_by_regime", False),
        ("7 · rank and compare", "summarize_givendata", False),
    ]
    h = 0.105
    for i, (title, fn, shaded) in enumerate(stages):
        y = 0.88 - i * 0.125
        _box(ax, 0.06, y, 0.62, h, "",
             facecolor=HILITE if shaded else "white",
             edgecolor=C["L5"] if shaded else INK, lw=1.5 if shaded else 1.2)
        ax.text(0.09, y + h / 2, title, fontsize=9.2, va="center",
                fontweight="bold", color=INK)
        ax.text(0.71, y + h / 2, fn, fontsize=7.8, va="center", color=MUTED,
                family="monospace")
        if i:
            _arrow(ax, (0.37, y + h + 0.020), (0.37, y + h + 0.001), color=MUTED, lw=1.3)

    ax.add_patch(mpatches.FancyBboxPatch(
        (0.035, 0.88 - 4 * 0.125 - 0.012), 0.665, 0.38,
        boxstyle="round,pad=0.01", facecolor="none", edgecolor=C["L5"],
        lw=1.6, linestyle="--"))
    ax.text(0.70, 0.88 - 3 * 0.125 + 0.17, "what makes it\nregime-stratified",
            fontsize=8.4, color=C["L5"], fontweight="bold", va="center")

    ax.set_title("Figure 17 — the sensitivity workflow", fontsize=S["fontsize_title"])
    return _save(fig, 17, "sa-pipeline")


def fig18():
    """Why the sampling has to be aimed. Schematic."""
    fig, axes = _fig(
        18, "aimed-sampling",
        "Schematic, not a run. An unaimed Latin hypercube over the default parameter "
        "ranges lands almost entirely in one regime, leaving the others with too few "
        "draws to support any index. Aiming the sample per regime is what creates the "
        "study; without it the other regimes are not under-resolved, they are absent.",
        figsize=(8.4, 3.8), ncols=2,
    )
    ax1, ax2 = axes
    regimes = ["metal", "oxide", "surface", "defect"]
    colors = [C["L5"], C["L2a"], C["L2b"], C["L3"]]

    unaimed = [96.0, 2.4, 1.0, 0.6]
    aimed = [25.0, 25.0, 25.0, 25.0]
    floor = 300
    N = 20000

    for ax, data, title in (
        (ax1, unaimed, "(A) unaimed LHS over the default ranges"),
        (ax2, aimed, "(B) targeted, sized per regime"),
    ):
        counts = [N * d / 100.0 for d in data]
        bars = ax.bar(regimes, counts, color=colors, edgecolor="white", width=0.66)
        ax.axhline(floor, color=WARN, ls="--", lw=1.6)
        ax.set_yscale("log")
        ax.set_ylim(1, N * 2)
        for b, c in zip(bars, counts):
            ax.text(b.get_x() + b.get_width() / 2, c * 1.25, f"{c:,.0f}",
                    ha="center", fontsize=8.0, fontweight="bold", color=INK)
        _style_axes(ax, None, "samples in regime" if ax is ax1 else None, title)
        ax.grid(axis="x", alpha=0)

    ax1.text(3.45, floor * 1.5, "minimum cluster\nfor an index",
             fontsize=8.0, color=WARN, fontweight="bold", ha="right")
    for i, c in enumerate([N * d / 100.0 for d in unaimed]):
        if c < floor:
            ax1.text(i, c * 0.35, "too few", ha="center", fontsize=8.0,
                     color=WARN, fontweight="bold")
    fig.suptitle("Figure 18 — the presets are what create the study (schematic)",
                 fontsize=S["fontsize_title"], y=1.02)
    fig.tight_layout()
    return _save(fig, 18, "aimed-sampling")


def fig19():
    """The dummy-parameter significance floor. Schematic."""
    fig, ax = _fig(
        19, "dummy-floor",
        "Schematic, not a run. A dummy parameter that the model never reads is "
        "carried through the whole analysis. Whatever index it earns is noise, so it "
        "sets the floor: anything at or below it is indistinguishable from a "
        "parameter with no effect. Reading a ranking without this line drawn invites "
        "treating noise as a result.",
        figsize=(7.4, 4.2),
    )
    names = ["$D_m$", "$K_s$", "$L_m$", "$f_{pinhole}$", "$E_b$",
             "$d_{grain}$", "dummy", "$k_{diss}$", "$\\gamma_{crack}$", "$N_T$"]
    values = [0.41, 0.33, 0.27, 0.19, 0.12, 0.085, 0.052, 0.048, 0.031, 0.024]
    floor = values[names.index("dummy")]

    colors = []
    for n, v in zip(names, values):
        if n == "dummy":
            colors.append(WARN)
        elif v > floor:
            colors.append(C["L5"])
        else:
            colors.append(RULE)

    order = np.argsort(values)
    ax.barh([names[i] for i in order], [values[i] for i in order],
            color=[colors[i] for i in order], edgecolor="white", height=0.68)
    ax.axvline(floor, color=WARN, ls="--", lw=1.8)
    ax.text(floor * 1.08, 0.2, "dummy floor", color=WARN, fontsize=9,
            fontweight="bold", rotation=90, va="bottom")

    ax.axvspan(0, floor, color=WARN, alpha=0.07)
    ax.text(floor * 0.5, 8.6, "indistinguishable\nfrom no effect", ha="center",
            fontsize=8.2, color=WARN, fontweight="bold")

    _style_axes(ax, r"importance index  $\delta$", None,
                "Figure 19 — the floor a dummy parameter sets (schematic)")
    ax.grid(axis="y", alpha=0)
    return _save(fig, 19, "dummy-floor")


def fig20():
    """Preset suppression: a pinned parameter cannot show its leverage. Schematic."""
    fig, axes = _fig(
        20, "preset-suppression",
        "Schematic, not a run. A preset that pins a parameter to force a regime also "
        "removes that parameter's variance, so the analysis reports it as "
        "unimportant. The index is correct and the conclusion drawn from it is wrong: "
        "a below-floor result for a pinned parameter says nothing about its physical "
        "leverage. Always read a ranking against the preset that produced it.",
        figsize=(8.6, 3.9), ncols=2,
    )
    ax1, ax2 = axes

    ax1.axhspan(-2.2, 0.0, color=C["L3"], alpha=0.16)
    ax1.plot([0, 1], [-2.2, -2.2], color=C["L3"], lw=3)
    ax1.plot([0, 1], [0.0, 0.0], color=C["L3"], lw=3)
    ax1.annotate("", xy=(0.5, 0.0), xytext=(0.5, -2.2),
                 arrowprops=dict(arrowstyle="<|-|>", color=C["L3"], lw=2.0))
    ax1.text(0.56, -1.1, "2.2 decades\nof real leverage", fontsize=9,
             color=C["L3"], fontweight="bold", va="center")
    ax1.plot([0.5], [-1.4], "o", color=WARN, ms=12, zorder=6,
             markeredgecolor="white", markeredgewidth=1.5)
    ax1.text(0.44, -1.4, "pinned here\nby the preset", fontsize=8.4, color=WARN,
             fontweight="bold", ha="right", va="center")
    ax1.set_xlim(0, 1)
    ax1.set_ylim(-2.6, 0.5)
    ax1.set_xticks([])
    _style_axes(ax1, None, r"$\log_{10} f_{pinhole}$",
                "(A) what the parameter could do")
    ax1.grid(axis="x", alpha=0)

    names = ["$D_m$", "$K_s$", "dummy", "$f_{pinhole}$"]
    vals = [0.41, 0.33, 0.052, 0.009]
    cols = [C["L5"], C["L5"], WARN, RULE]
    ax2.barh(names, vals, color=cols, edgecolor="white", height=0.62)
    ax2.axvline(0.052, color=WARN, ls="--", lw=1.7)
    ax2.annotate("reads as unimportant —\nbut only because it was pinned",
                 xy=(0.009, 3), xytext=(0.14, 2.35), fontsize=8.3, color=WARN,
                 fontweight="bold",
                 arrowprops=dict(arrowstyle="-|>", color=WARN, lw=1.3))
    _style_axes(ax2, r"importance index  $\delta$", None,
                "(B) what the analysis reports")
    ax2.grid(axis="y", alpha=0)

    fig.suptitle("Figure 20 — a preset can blind the analysis (schematic)",
                 fontsize=S["fontsize_title"], y=1.02)
    fig.tight_layout()
    return _save(fig, 20, "preset-suppression")


# =============================================================================

FIGURES = {
    1: fig01, 2: fig02, 3: fig03, 4: fig04, 5: fig05,
    6: fig06, 7: fig07, 8: fig08, 9: fig09, 10: fig10,
    11: fig11, 12: fig12, 13: fig13, 14: fig14, 15: fig15,
    16: fig16, 17: fig17, 18: fig18, 19: fig19, 20: fig20,
}


def write_caption_index():
    """Write captions.md, the single list every chapter's figure caption comes from."""
    lines = [
        "# Figure captions",
        "",
        "Generated by `make_figures.py`. Do not edit by hand: edit the caption in the",
        "figure's own function and regenerate, so the caption and the drawing cannot",
        "drift apart.",
        "",
        "| # | file | caption |",
        "|---|---|---|",
    ]
    for n in sorted(_captions):
        slug, caption = _captions[n]
        lines.append(f"| {n} | `fig-{n:02d}-{slug}.{EXT}` | {caption} |")
    path = OUT / "captions.md"
    path.write_text("\n".join(lines) + "\n")
    print(f"\n  captions -> {path.name}")


def main(argv):
    wanted = [int(a) for a in argv[1:]] or sorted(FIGURES)
    unknown = [n for n in wanted if n not in FIGURES]
    if unknown:
        print(f"no such figure(s): {unknown}", file=sys.stderr)
        return 2

    print(f"Regenerating {len(wanted)} figure(s) from the Example Wall into {OUT}/\n")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for n in wanted:
            FIGURES[n]()

    if wanted == sorted(FIGURES):
        write_caption_index()
    print(f"\n{len(wanted)} / {len(wanted)} figures written")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
