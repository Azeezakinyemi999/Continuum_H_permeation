"""The figures of the parameter sweeps, in the house style: J and Φ_app against the swept variable.

Used by the four effect-of notebooks (one alloy, at the study's operating point or at three
temperatures or pressures) and by Comparisons/study-comparison.ipynb (one curve per alloy). Each
notebook figure becomes one figure per panel (``flux`` and ``phi``; ``<layer>-flux`` and
``<layer>-phi`` where one notebook figure held several sweeps), each a ``figures.Panel`` with a LaTeX
caption whose numbers come from the data; the notebook writes them with ``figures.save_panels``.

On each figure the marks show the limiting step at every point (colour = the step, an open mark = a
defect path carries most of the flux, a diamond = trap control; ``house_style.regime_mark``), and
each limiting step is fitted within each curve: an Arrhenius fit for temperature, a power law
otherwise, drawn as a line in the step's colour (line style per curve). The fitted exponents are on
the plot, as the second line of each limiting step's entry in the key below it. Several curves are
named at their right-hand end. Pressure figures of J carry grey reference slopes of 1 and 0.5,
labelled; the pinhole figures a grey line for the no-pinhole limit (f -> 0) of each curve.

The notebooks keep computing their own points and fit tables; nothing here changes a number.
"""
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator
from scipy import stats

from calculations import house_style as hs
from calculations.figures import Panel
from calculations.permeability import UNITS_HENRY, UNITS_SIEVERTS

MODULE = "calculations/sweep_figures.py"
R_GAS = 8.314  # J/(mol K)
LINE_STYLES = ["-", (0, (6, 3)), (0, (1.5, 1.5))]            # one per curve, in order

#: sweep key: (x label, log x, what the caption calls x, caption unit, SI -> axis factor; None = 1000/T)
SWEEP_AXES = {
    "pressure": (r"Upstream pressure $P_\mathrm{up}$ (Pa)", True, "upstream pressure", "Pa", 1.0),
    "temperature": (r"1000/$T$ (K$^{-1}$)", False, "temperature", "K", None),
    "oxide-thickness": (r"Oxide thickness $L_\mathrm{ox}$ (µm)", True, "oxide thickness", "µm", 1e6),
    "metal-thickness": (r"Metal thickness $L_\mathrm{m}$ (mm)", True, "metal thickness", "mm", 1e3),
    "total-thickness": (r"Wall thickness $L_\mathrm{total}$ (mm)", True,
                        "whole-wall thickness (oxide and metal scaled together)", "mm", 1e3),
    "pinholes": (r"Pinhole area fraction $f$", True, "pinhole area fraction", "", 1.0),
}
UNITS_TEX = {UNITS_SIEVERTS: r"mol m$^{-1}$ s$^{-1}$ Pa$^{-1/2}$", UNITS_HENRY: r"mol m$^{-1}$ s$^{-1}$ Pa$^{-1}$"}
MERGE_WORD = {"3T": "temperatures", "3P": "pressures", "studies": "alloys", "op": "conditions"}
#: The thickness sweeps of a notebook figure row, as sweep keys.
THICKNESS_KEYS = {"oxide": "oxide-thickness", "metal": "metal-thickness", "total": "total-thickness"}


def level_layers(level):
    """The layers of a level's wall: {'metal'}, {'oxide'} or both."""
    base = level[:-2] if level.endswith("L6") else level
    return {"L1": {"metal"}, "L4": {"metal"}, "L2a": {"oxide"}}.get(base, {"oxide", "metal"})


def alloy_name(metal_key):
    """The alloy's name on a figure, from a study's metal key ('316L_Guo_2025' -> '316L'); else the key."""
    return hs.ALLOY_NAMES.get(metal_key, metal_key)


def p_label(p):
    """A pressure for a curve name: '$10^{-7}$ Pa', '1 Pa', '$10^{7}$ Pa', or '316 Pa'."""
    k = int(np.floor(np.log10(p) + 1e-9))
    if np.isclose(p, 10.0 ** k):
        return "1 Pa" if k == 0 else rf"$10^{{{k}}}$ Pa"
    return f"{p:.3g} Pa"


def curve(series, x, J, phi, regime, T, P, L_ox, L_m, J0=None, phi0=None):
    """One curve of a sweep: its name, the swept values in SI (x: Pa, K, m or the area fraction), the
    flux, Φ_app and regime label at each point, and the conditions (scalars or arrays). J0 and phi0:
    the no-pinhole limit of a pinhole sweep."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    full = lambda v: np.broadcast_to(np.asarray(v, dtype=float), (n,)).copy()   # noqa: E731
    return dict(series=series, x_si=x, J=np.asarray(J, dtype=float), phi=np.asarray(phi, dtype=float),
                regime=np.asarray(regime).astype(str), T=full(T), P=full(P), L_ox=full(L_ox), L_m=full(L_m),
                J0=J0, phi0=phi0)


def _frame(curves, sweep_key):
    factor = SWEEP_AXES[sweep_key][4]
    rows = []
    for c in curves:
        d = pd.DataFrame({k: c[k] for k in ("x_si", "J", "phi", "regime", "T", "P", "L_ox", "L_m")})
        d["series"] = c["series"]
        rows.append(d)
    d = pd.concat(rows, ignore_index=True)
    d["x"] = 1000.0 / d["T"] if factor is None else d.x_si * factor
    return d


def fit(sweep_key, x_si, T, y):
    """(exponent, pre-factor, fitted y) of one group of points: Arrhenius for temperature (E_a in
    kJ/mol from ln y against 1000/T, pre-factor in the units of y), else a power law y = A x^n (SI)."""
    if sweep_key == "temperature":
        f = stats.linregress(1000.0 / T, np.log(y))
        return -f.slope * R_GAS, np.exp(f.intercept), np.exp(f.intercept + f.slope * 1000.0 / T)
    f = stats.linregress(np.log(x_si), np.log(y))
    return f.slope, np.exp(f.intercept), np.exp(f.intercept) * x_si ** f.slope


def exponent_text(sweep_key, qty, values, prefactor=None):
    """The fitted exponent(s) of one limiting step, for the key: 'E_a = 52 kJ mol-1', 'n = 0.50',
    'slope 0.00' (Φ_app of a pressure or thickness sweep), with a range over several curves."""
    if sweep_key == "temperature":
        text = rf"$E_a$ = {hs.value_range(values, '{:.0f}')} kJ mol$^{{-1}}$"
        if prefactor is not None:
            text += f", $A$ = {hs.sci(prefactor, 2)}"
        return text
    if qty == "phi" and sweep_key != "pinholes":
        return f"slope {hs.value_range(values)}"
    symbol = "n" if sweep_key in ("pressure", "pinholes") else "m"
    return f"${symbol}$ = {hs.value_range(values)}"


def _tidy_linear_y(ax):
    """A y axis set_log_y made linear (values within a factor of 1.5): at most five ticks, written
    '3.21×10⁻⁶' when the values are very small or large, with just enough digits to tell them apart."""
    if ax.get_yscale() != "linear":
        return
    lo, hi = ax.get_ylim()
    ticks = [t for t in MaxNLocator(5).tick_values(lo, hi) if lo <= t <= hi]
    if len(ticks) < 2:
        return
    ax.yaxis.set_major_locator(FixedLocator(ticks))           # only ticks inside the axis: no hidden labels
    step = np.min(np.diff(ticks))
    top = max(abs(t) for t in ticks)
    k = int(np.floor(np.log10(top))) if top > 0 else 0
    digits = max(0, int(np.ceil(-np.log10(step / 10.0 ** k) - 1e-9)))
    if -2 <= k <= 3:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.{max(0, digits - k)}f}"))
    else:
        ax.yaxis.set_major_formatter(FuncFormatter(
            lambda v, _: "0" if v == 0 else rf"${v / 10.0 ** k:.{digits}f}\times10^{{{k}}}$".replace("-", "\u2212", 1)
            if v < 0 else rf"${v / 10.0 ** k:.{digits}f}\times10^{{{k}}}$"))


def _ticks(ax, sweep_key, xlog):
    if xlog:
        lo, hi = ax.get_xlim()
        hs.log_ticks(ax.xaxis, lo, min(hi, 1.0) if sweep_key == "pinholes" else hi)   # f is at most 1
    if ax.get_yscale() == "log":
        hs.log_ticks(ax.yaxis, *ax.get_ylim())
    else:
        _tidy_linear_y(ax)


def _draw(d, sweep_key, qty, phi_units, merge_word, limits=None):
    """One figure: marks by limiting step, a fit line per step and curve, curve names at the right-hand
    end, the reference slopes (pressure, J), the no-pinhole limits; the key below the plot carries
    each step's fitted exponents. Returns (fig, ax, groups, guide_series, keyed)."""
    xlabel, xlog = SWEEP_AXES[sweep_key][:2]
    d = d.assign(y=d[qty if qty == "phi" else "J"])
    d = d[np.isfinite(d.y) & (d.y > 0)]
    names = list(dict.fromkeys(d.series))
    fig, ax = hs.new_figure()
    order = d.regime.value_counts().index                   # largest group first
    for reg in order:
        m = d.regime == reg
        colour, marker, hollow = hs.regime_mark(reg)
        ax.plot(d.x[m], d.y[m], marker, ls="none", ms=9, color=colour,
                mfc="white" if hollow else colour, mec=colour if hollow else "white",
                mew=2.0 if hollow else 1.0, zorder=3)

    # fits: one per limiting step within each curve, as the notebooks fit them
    fits = {}                                               # regime -> [(exponent, pre-factor), ...]
    for ci, name in enumerate(names):
        s = d[d.series == name]
        for reg in dict.fromkeys(s.regime):
            g = s[s.regime == reg].sort_values("x")
            if not hs.is_fitted(reg) or len(g) < 3:
                continue
            k, A, line = fit(sweep_key, g.x_si.to_numpy(), g["T"].to_numpy(), g.y.to_numpy())
            fits.setdefault(reg, []).append((k, A))
            ax.plot(g.x, line, color=hs.regime_mark(reg)[0], ls=LINE_STYLES[ci % 3], lw=hs.CURVE_PT, zorder=2)

    # the no-pinhole limits (pinhole sweeps): a grey line per curve, in its line style
    ref_y = []
    for ci, (name, y0) in enumerate((limits or {}).items()):
        if y0 is not None and np.isfinite(y0) and y0 > 0:
            ax.axhline(y0, color=hs.GREY, lw=2.4, ls=LINE_STYLES[ci % 3], zorder=1)
            ref_y.append(y0)

    hs.set_log_y(ax, np.r_[d.y.to_numpy(), ref_y], pad_decades=0.35)
    if xlog:
        ax.set_xscale("log")
        ax.set_xlim(*hs.log_limits(d.x, 0.15, frac=0.03))
    else:
        span = d.x.max() - d.x.min()
        ax.set_xlim(d.x.min() - 0.04 * span, d.x.max() + 0.04 * span)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"Flux $J$ (mol m$^{-2}$ s$^{-1}$)" if qty == "J"
                  else rf"$\Phi_\mathrm{{app}}$ ({UNITS_TEX[phi_units]})")

    ends = []
    if len(names) > 1:                                      # one curve needs no name: the caption gives it
        for name in names:                                  # curves that end together share one name
            a = d[d.series == name]
            i = a.x.idxmax()
            for e in ends:
                if np.isclose(e[0], a.x[i]) and abs(np.log10(e[1] / a.y[i])) < 0.006:
                    e[2].append(name)
                    break
            else:
                ends.append((a.x[i], a.y[i], [name]))
        ends = [(x, y, f"All {len(ns)} {merge_word}".replace("All 3 ", "All three ")
                 if len(ns) == len(names) and len(ns) > 2
                 else ("Both " + merge_word if len(ns) == len(names) == 2 else hs.tex_join(ns)))
                for x, y, ns in ends]

    _ticks(ax, sweep_key, xlog)
    top = hs.celsius_top_axis(ax) if sweep_key == "temperature" else None
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    end_texts = []
    if ends:
        hs.room_for_end_labels(ax, [n for _, _, n in ends], renderer, log=xlog)
        for _ in range(4):                                  # end labels, widening the x-axis until they fit
            _ticks(ax, sweep_key, xlog)
            fig.canvas.draw()
            end_texts = hs.end_labels(ax, ends)
            _ticks(ax, sweep_key, xlog)                     # end_labels may have extended the y-axis
            fig.canvas.draw()
            over = max(hs.bbox_frac(ax, t, renderer).x1 for t in end_texts) - 0.985
            if over <= 0:
                break
            for t in end_texts:
                t.remove()
            x0, x1 = ax.get_xlim()
            if xlog:
                u0, u1 = np.log10(x0), np.log10(x1)
                ax.set_xlim(x0, 10 ** (u1 + (over + 0.02) * (u1 - u0) * 1.15))
            else:
                ax.set_xlim(x0, x1 + (over + 0.02) * (x1 - x0) * 1.15)
    if top is not None:
        hs.set_celsius_ticks(ax, top)

    # reference slopes on the pressure figures of J: n = 1 (J ∝ P) and n = 0.5 (J ∝ √P), through the
    # middle point of the middle curve, labelled on the plot; grey, they step back from the data
    guide_series = None
    if sweep_key == "pressure" and qty == "J":
        guide_series = names[len(names) // 2]
        a = d[d.series == guide_series].sort_values("x")
        x_i, y_i = a.x.iloc[len(a) // 2], a.y.iloc[len(a) // 2]
        (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
        xs = np.logspace(np.log10(d.x.min()), np.log10(d.x.max()), 200)
        log_y = ax.get_yscale() == "log"
        fig.canvas.draw()
        points = ax.transData.transform(np.column_stack([d.x, d.y]))      # display coordinates of the marks
        taken = [t.get_window_extent(renderer) for t in end_texts]
        for n, dash in ((1.0, (0, (7, 3))), (0.5, (0, (2, 2)))):
            g = y_i * (xs / x_i) ** n
            ax.plot(xs, g, color=hs.GREY, ls=dash, lw=2.4, zorder=1)
            inside = np.flatnonzero((g > y0 * (1.6 if log_y else 1)) & (g < y1 / (1.6 if log_y else 1)))
            if not inside.size:
                continue
            best = None                                                     # (marks covered, k, dy)
            axbox = ax.get_window_extent(renderer)
            for frac in (0.8, 0.65, 0.9, 0.5, 0.35, 0.2, 0.05):           # the first spot that covers no mark
                k = inside[int(frac * (inside.size - 1))]
                for dy in ((10, 22) if n == 1 else (-18, -30)):
                    label = ax.annotate(f"slope {n:g}", (xs[k], g[k]), xytext=(6, dy), textcoords="offset points",
                                        fontsize=hs.ANNOT_PT, fontweight="bold", color="black", zorder=5,
                                        bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
                    box = label.get_window_extent(renderer).expanded(1.05, 1.15)
                    label.remove()
                    if (box.x1 > axbox.x1 or box.y0 < axbox.y0 or box.y1 > axbox.y1
                            or any(box.overlaps(t) for t in taken)):
                        continue
                    covered = box.count_contains(points)
                    if best is None or covered < best[0]:
                        best = (covered, k, dy, box)
                    if not covered:
                        break
                if best is not None and best[0] == 0:
                    break
            if best is not None:                                            # else the fewest marks covered
                _, k, dy, box = best
                ax.annotate(f"slope {n:g}", (xs[k], g[k]), xytext=(6, dy), textcoords="offset points",
                            fontsize=hs.ANNOT_PT, fontweight="bold", color="black", zorder=5,
                            bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
                taken.append(box)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)

    # the key: one entry per limiting step, its fitted exponents on the second line
    groups = {}
    for reg in order:
        short = hs.regime_short(reg).replace("\n", ", ")
        colour, marker, hollow = hs.regime_mark(reg)
        g = groups.setdefault((short, colour), {"regs": [], "marker": marker, "hollow": hollow, "fits": []})
        g["regs"].append(reg)
        g["fits"] += fits.get(reg, [])
    handles, labels, colours = [], [], []
    for (short, colour), g in groups.items():
        values = np.array([k for k, _ in g["fits"]])
        g["values"] = values
        text = short
        if len(values):
            single_A = g["fits"][0][1] if len(names) == 1 and len(g["fits"]) == 1 else None
            text += "\n" + exponent_text(sweep_key, qty, values, single_A)
        handles.append(Line2D([], [], ls="none", ms=10, marker=g["marker"], color=colour,
                              mfc="white" if g["hollow"] else colour, mec=colour, mew=2.0 if g["hollow"] else 1.0))
        labels.append(text)
        colours.append(colour)
    if limits:
        handles.append(Line2D([], [], color=hs.GREY, lw=2.4))
        labels.append(r"No pinholes ($f \rightarrow 0$)")
        colours.append("black")
    hs.key_below(fig, handles, labels, ncol=2, colours=colours)
    return fig, ax, d, groups, guide_series


# --------------------------------------------------------------------------------------------------
# Captions
# --------------------------------------------------------------------------------------------------


def _marks_text(groups, sweep_key, qty):
    """Caption words for the marks: each limiting step, its marker and its fitted exponent."""
    marks = []
    for (short, _), g in groups.items():
        words = ("no single limiting step" if short == "Mixed"
                 else hs.regime_long(g["regs"][0]).replace(" (not fitted)", "").replace(": ", ", "))
        words = words[0].lower() + words[1:]
        shape = ("open " if g["hollow"] else "") + hs.MARKER_WORDS.get(g["marker"], "marks")
        if len(g["values"]):
            fit_words = exponent_text(sweep_key, qty, g["values"]).replace("−", "$-$")
            marks.append(f"{words} ({shape}; {fit_words})")
        elif any(hs.is_fitted(r) for r in g["regs"]):
            marks.append(f"{words} ({shape}; too few points to fit)")
        else:
            marks.append(f"{words} ({shape}; not fitted)")
    return hs.tex_join(marks)


def _quantity_words(qty, phi_units):
    if qty == "J":
        return "Steady-state hydrogen flux $J$"
    drive = r"\Delta P" if phi_units == UNITS_HENRY else r"\Delta\sqrt{P}"
    return rf"Apparent permeability $\Phi_\mathrm{{app}} = J\,L_\mathrm{{total}}/{drive}$"


def _range_words(d, sweep_key):
    unit = SWEEP_AXES[sweep_key][3]
    if sweep_key == "temperature":                          # kelvin as computed: 573 to 1273 K
        return f"{d['T'].min():.0f} to {d['T'].max():.0f}~K"
    unit_tex = {"µm": r"\textmu m"}.get(unit, unit)
    return f"{hs.tex_num(d.x.min())} to {hs.tex_num(d.x.max())}{('~' + unit_tex) if unit_tex else ''}"


def _conditions(d, sweep_key, variant, layers, P_down, own):
    """Caption words for what was held fixed, read from the curves themselves; `own` names a curve
    when a held value differs between curves ('alloy': 'each alloy at its own temperature')."""
    def one(col, fmt, what):
        vals = np.unique(np.round(d[col].dropna().to_numpy(dtype=float), 12))
        return fmt(vals[0]) if len(vals) == 1 else f"each {own} at its own {what}"
    parts = []
    if variant == "3T":
        parts.append("at " + hs.tex_join([f"{t:.0f}" for t in sorted(d["T"].unique())]) + "~K (one curve each)")
    elif sweep_key != "temperature":
        parts.append(one("T", lambda v: f"at {v:.0f}~K", "temperature"))
    if variant == "3P":
        parts.append("at upstream pressures " + hs.tex_join([hs.tex_num(p) for p in sorted(d["P"].unique())])
                     + "~Pa (one curve each)")
    elif sweep_key != "pressure":
        parts.append(one("P", lambda v: f"upstream pressure {hs.tex_num(v)}~Pa", "upstream pressure"))
    micron = r"~\textmu m"
    if "oxide" in layers and sweep_key not in ("oxide-thickness", "total-thickness"):
        parts.append(one("L_ox", lambda v: f"oxide {v * 1e6:g}{micron}", "oxide thickness"))
    if "metal" in layers and sweep_key not in ("metal-thickness", "total-thickness"):
        parts.append(one("L_m", lambda v: f"metal {v * 1e3:g}~mm", "metal thickness"))
    if sweep_key == "total-thickness" and len(layers) == 2:
        r = np.unique(np.round((d.L_m / d.L_ox).to_numpy(), 6))
        parts.append(f"oxide to metal 1:{r[0]:.0f}" if len(r) == 1 else f"each {own}'s own oxide-to-metal ratio")
    pdn = np.unique(np.atleast_1d(P_down))
    parts.append(f"downstream {hs.tex_num(pdn[0])}~Pa" if len(pdn) == 1 else f"each {own}'s own downstream pressure")
    if sweep_key == "pinholes":
        parts.append("pinholes are the only oxide defects")
    return ", ".join(parts)


def _caption(level, sweep_key, qty, d, groups, phi_units, variant, alloy, P_down, guide_series, limits):
    n_pts = d.groupby("series").size()
    each = " each" if len(n_pts) > 1 else ""
    n_txt = f"{n_pts.iloc[0]} points{each}" if n_pts.nunique() == 1 else "points per curve as computed"
    own = "alloy" if variant == "studies" else "curve"
    who = (f" for {hs.tex_join([str(n).replace(' ', '~') for n in dict.fromkeys(d.series)])}" if variant == "studies"
           else (f" ({alloy.replace(' ', '~')})" if alloy else ""))
    fit_kind = "Arrhenius fits" if sweep_key == "temperature" else "power-law fits"
    multi = d.series.nunique() > 1
    text = (f"{_quantity_words(qty, phi_units)} through {hs.level_wall(level)}{who}, against "
            f"{SWEEP_AXES[sweep_key][2]} ({n_txt}, {_range_words(d, sweep_key)}); "
            f"{_conditions(d, sweep_key, variant, level_layers(level), P_down, own)}. "
            f"Marks show the limiting step at each point: {_marks_text(groups, sweep_key, qty)}. "
            f"Lines are {fit_kind} within each limiting step"
            + (", one per curve (solid, dashed and dotted in the order of the curves named at their ends)" if multi
               else "") + "; their exponents are in the key.")
    if guide_series:
        g = str(guide_series).replace(" ", "~")
        text += (rf" Grey dashed lines: slopes of 1 (long dashes, $J \propto P$) and 0.5 (short dashes, "
                 rf"$J \propto \sqrt{{P}}$) through the middle point of the {g} curve.")
    if limits:
        text += r" Grey lines: each curve's no-pinhole limit ($f \rightarrow 0$)."
    return text


# --------------------------------------------------------------------------------------------------
# The figures
# --------------------------------------------------------------------------------------------------


def sweep_panels(sweep_key, level, curves, phi_units, P_down, variant="op", alloy=None, group=""):
    """The two figures of one sweep, J and Φ_app, as ``figures.Panel``s named '<group>-flux' and
    '<group>-phi' (just 'flux' and 'phi' without a group).

    sweep_key: 'temperature', 'pressure', 'oxide-thickness', 'metal-thickness', 'total-thickness' or
    'pinholes'; curves: ``curve(...)`` records, one per curve; variant: 'op' (one curve at the
    operating point), '3T' or '3P' (three temperatures or pressures), 'studies' (one curve per alloy);
    alloy: the alloy's name, for the caption. Returns [] when no curve has a point."""
    curves = [c for c in curves if c is not None and len(c["x_si"])]
    if not curves:
        return []
    d0 = _frame(curves, sweep_key)
    merge = MERGE_WORD.get(variant, "curves")
    out = []
    for qty, name in (("J", "flux"), ("phi", "phi")):
        limits = None
        key = "J0" if qty == "J" else "phi0"
        if sweep_key == "pinholes" and any(c.get(key) is not None for c in curves):
            limits = {c["series"]: c.get(key) for c in curves}
        if not np.isfinite(d0[qty]).any():
            continue
        fig, ax, d, groups, guide = _draw(d0, sweep_key, qty, phi_units, merge, limits)
        caption = _caption(level, sweep_key, qty, d, groups, phi_units, variant, alloy, P_down, guide, limits)
        out.append(Panel(f"{group}-{name}" if group else name, fig, ax, caption, d.x.to_numpy(), d.y.to_numpy(),
                         keyed=True))
    return out


def show(panels):
    """Show each panel's figure in the notebook (no-op outside IPython), before save_panels closes it."""
    try:
        from IPython import get_ipython
        from IPython.display import display
    except ImportError:
        return panels
    if get_ipython() is None:
        return panels
    for p in panels:
        display(p.fig)
    return panels
