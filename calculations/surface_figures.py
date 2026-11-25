"""The figures of the Surface_proposal notebook: the surface-kinetics (+L6) levels, one figure per panel.

Each function takes the arrays the notebook computed (the physics stays in the notebook) and returns
a list of ``calculations.figures.Panel``; the notebook writes them with ``figures.save_panels``. Every
figure is in the house style (``calculations/house_style.py``): 180 mm wide, black text, no titles,
a key below the plot where labels do not fit beside the data.

Colours keep their meaning from the other figures of the project:
  - a rate-limiting step has its step colour (``house_style.STEP_COLOURS``): oxide surface green,
    metal surface purple, oxide diffusion vermillion, metal diffusion blue, mixed grey;
  - a model level has its level hue (``house_style.level_style``); a +L6 level keeps its base level's
    hue (L2+L6 is L2b, L2+L4+L6 is L4, L3+L4+L6 is L5);
  - defect types and other conditions are black, told apart by marker and line style;
  - fits, theory lines, limits and thresholds are grey dashed or dotted lines that step back.
Fitted values (slopes, activation energies, crossover pressures) sit on the figure, as the second line
of a key entry in the series' colour or as short check lines in a clear corner, and in the caption.

Panel names (files ``<figure>-<panel>``), the same in every figure that has the panel:
  flux        flux against upstream pressure
  coverage    surface coverage against pressure (L1+L6, L2a+L6)
  arrhenius   flux against inverse temperature at one pressure per regime (L2+L6)
  defects     flux against defect fraction, one line per defect type (L3+L6)
  diffusivity D_eff / D_m against temperature, one line per microstructure mode (L2+L4+L6)
  enhancement flux ratios against pressure (L3+L4+L6)
  resistance  resistance fractions against pressure
  limit       the limit check: the +L6 flux against the flux it must recover
"""
from itertools import groupby

import numpy as np
from matplotlib.lines import Line2D
from scipy import stats

from calculations import house_style as hs
from calculations.figures import Panel

MODULE = "calculations/surface_figures.py"
P_LABEL = r"Upstream pressure $P_\mathrm{up}$ (Pa)"
J_UNIT = r"mol m$^{-2}$ s$^{-1}$"
GREY = hs.GREY
DOTTED = (0, (1, 1.6))
DASHDOT = (0, (6, 2, 1.5, 2))
REF_DASHED = dict(color=GREY, lw=2.4, ls=(0, (5, 3)))
REF_DOTTED = dict(color=GREY, lw=2.4, ls=DOTTED)
REF_DASHDOT = dict(color=GREY, lw=2.4, ls=DASHDOT)

STEP_NAMES = {"oxide surface": "Oxide-surface limited", "metal surface": "Metal-surface limited",
              "oxide diffusion": "Oxide-diffusion limited", "metal diffusion": "Metal-diffusion limited",
              "mixed": "Mixed"}
STEP_WORDS = {k: v.lower() for k, v in STEP_NAMES.items()}
STEP_MARKERS = {"oxide surface": "P", "metal surface": "P", "oxide diffusion": "s", "metal diffusion": "o",
                "mixed": "X"}
#: resistance fractions: the step's colour, and a line style and marker each, so colour is not the only cue
FRACTION_STYLE = {"surface": ("-", "P"), "oxide": ("--", "s"), "metal": (DOTTED, "o")}
FRACTION_NAMES = {"surface": "Surface (dissociation)", "oxide": "Oxide (diffusion)", "metal": "Metal (diffusion)"}
#: conditions (defect types, defect fractions, modes) in one hue: line style and marker apart
CONDITION_MARKS = [("-", "o"), ("--", "s"), (DOTTED, "^"), (DASHDOT, "D"), ((0, (3, 1.5)), "v")]


# --------------------------------------------------------------------------------------------------
# Words and numbers
# --------------------------------------------------------------------------------------------------

def _pa(v):
    """'1e-7 Pa' style pressure for captions: $10^{-7}$~Pa, 13~Pa."""
    return f"{hs.tex_num(v, 2)}~Pa"


def _pa_key(v):
    return f"{hs.sci(v, 2)} Pa"


def _pct(v):
    """A percentage for a figure: 79, 0.49, 2.8x10^-14 (mathtext)."""
    return hs.sci(v, 2) + " %"


def _tex_pct(v):
    return hs.tex_num(v, 2) + r"~\%"


def _steps(surface):
    """regime word -> limiting step, for a figure whose surface step is `surface`."""
    return {"surface": surface, "oxide": "oxide diffusion", "metal": "metal diffusion", "mixed": "mixed"}


def regimes(frac_surface, frac_oxide=None, frac_metal=None):
    """The notebook's regime per point: the step above half of the total resistance, else 'mixed'."""
    fs = np.asarray(frac_surface, float)
    fo = np.full_like(fs, np.nan) if frac_oxide is None else np.asarray(frac_oxide, float)
    fm = np.full_like(fs, np.nan) if frac_metal is None else np.asarray(frac_metal, float)
    with np.errstate(invalid="ignore"):
        return np.where(fs > 0.5, "surface", np.where(fo > 0.5, "oxide", np.where(fm > 0.5, "metal", "mixed")))


def _ok(*arrays, positive=True):
    m = np.ones(len(arrays[0]), bool)
    for a in arrays:
        a = np.asarray(a, float)
        m &= np.isfinite(a)
        if positive:
            m &= a > 0
    return m


def segment_slopes(P, J, rl):
    """{regime: [slope of each run of three or more consecutive points]} of log J against log P."""
    P, J = np.asarray(P, float), np.asarray(J, float)
    idx = np.flatnonzero(_ok(P, J))
    out = {}
    for regime in dict.fromkeys(rl[idx]):
        sel = idx[rl[idx] == regime]
        for _, g in groupby(enumerate(sel), lambda t: t[0] - t[1]):
            grp = [k for _, k in g]
            if len(grp) >= 3:
                out.setdefault(regime, []).append(stats.linregress(np.log10(P[grp]), np.log10(J[grp])).slope)
    return out


def net_slope(P, J):
    m = _ok(P, J)
    if m.sum() < 3:
        return np.nan
    return stats.linregress(np.log10(np.asarray(P, float)[m]), np.log10(np.asarray(J, float)[m])).slope


def _slope_text(slopes, sep=", "):
    return sep.join(f"{s:.2f}" for s in slopes)


def _regime_slope_line(P, J, rl, surface):
    """'surface-limited slope 0.99; oxide-limited 0.51' for a key's second line."""
    parts = []
    for regime, slopes in segment_slopes(P, J, rl).items():
        if regime != "mixed":
            parts.append(f"{regime}-limited slope {_slope_text(slopes)}")
    return "; ".join(parts)


def coincide_words(P, curves, names, tol=0.02):
    """'The curves of B and C lie on that of A to within 0.3 %.' for curves (on one P) within tol of the
    first one, else ''."""
    base = np.asarray(curves[0], float)
    close, worst = [], 0.0
    for J, name in zip(curves[1:], names[1:]):
        J = np.asarray(J, float)
        m = _ok(base, J)
        if m.sum() < 2:
            continue
        d = float(np.max(np.abs(J[m] - base[m]) / base[m]))
        if d < tol:
            close.append(name)
            worst = max(worst, d)
    if not close:
        return ""
    return (f" The curve{'s' if len(close) > 1 else ''} of {hs.tex_join(close)} lie{'' if len(close) > 1 else 's'} "
            f"on that of {names[0]} to within {_tex_pct(100 * worst) if worst > 0 else 'rounding'}, so the "
            "markers are staggered to show each.")


def _check_lines(max_err, mean_err=None, name=None, threshold=1.0):
    """The notebook's limit check as short lines: errors, and 'validated' or 'check failed'."""
    head = f"{name}: " if name else ""
    if max_err is None or not np.isfinite(max_err):
        return [f"{head}no valid points: check failed"]
    lines = [f"{head}max error {_pct(max_err)}"]
    if mean_err is not None and np.isfinite(mean_err):
        lines.append(f"mean error {_pct(mean_err)}")
    lines.append("validated (max error under 1 %)" if max_err < threshold else "check failed: max error above 1 %")
    return lines


# --------------------------------------------------------------------------------------------------
# Drawing pieces
# --------------------------------------------------------------------------------------------------

def _place_note(ax, lines, points, avoid=()):
    """Write the check lines in the corner of the axes that covers the fewest data points."""
    if not lines:
        return None
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pts = np.vstack([hs._to_axes(ax, x, y) for x, y in points if len(x)]) if points else np.empty((0, 2))
    pts = pts[np.isfinite(pts).all(axis=1)] if len(pts) else pts
    taken = [hs.bbox_frac(ax, a, renderer) for a in avoid]
    text = "\n".join(lines)
    t = ax.text(0, 0, text, transform=ax.transAxes, fontsize=hs.ANNOT_PT, fontweight="bold", color="black",
                linespacing=1.25, zorder=7)
    best = None
    for x, y, ha, va in ((0.03, 0.96, "left", "top"), (0.97, 0.04, "right", "bottom"),
                         (0.97, 0.96, "right", "top"), (0.03, 0.04, "left", "bottom"),
                         (0.97, 0.5, "right", "center"), (0.03, 0.5, "left", "center"),
                         (0.5, 0.96, "center", "top"), (0.5, 0.04, "center", "bottom")):
        t.set_position((x, y)); t.set_ha(ha); t.set_va(va)
        b = hs.bbox_frac(ax, t, renderer)
        covered = int(((pts[:, 0] > b.x0 - 0.01) & (pts[:, 0] < b.x1 + 0.01)
                       & (pts[:, 1] > b.y0 - 0.01) & (pts[:, 1] < b.y1 + 0.01)).sum()) if len(pts) else 0
        clash = sum(not (b.x1 < o.x0 or b.x0 > o.x1 or b.y1 < o.y0 or b.y0 > o.y1) for o in taken)
        score = covered + 1000 * clash
        if best is None or score < best[0]:
            best = (score, (x, y, ha, va))
    x, y, ha, va = best[1]
    t.set_position((x, y)); t.set_ha(ha); t.set_va(va)
    t.set_bbox(dict(facecolor="white", edgecolor="none", alpha=0.9, pad=2))
    return t


def _log_axes(ax, x, y, xlog=True, ylog=True, x_pad=0.1):
    if xlog:
        lo, hi = hs.log_limits(x, x_pad, frac=0.03)
        ax.set_xscale("log")
        ax.set_xlim(lo, hi)
        hs.log_ticks(ax.xaxis, lo, hi)
    if ylog:
        hs.set_log_y(ax, y)


def _style(level=None, **over):
    s = hs.level_style(level) if level else dict(color="black", ls="-", marker="o")
    s.update(over)
    return s


def _plot(ax, x, y, style, lw=hs.CURVE_PT, ms=8, n_marks=14, zorder=3, stagger=(0, 1)):
    """A line with sparse markers, in one style dict (color, ls, marker). stagger=(k, n): the k-th of n
    curves that may lie on top of each other puts its markers between the others', so each shows."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    marker = style.get("marker")
    every = max(1, len(x) // n_marks) if marker else None
    if marker and stagger[1] > 1:
        every = max(stagger[1], len(x) // (8 if stagger[1] > 2 else 11))   # room between the curves' markers
        every = (stagger[0] * every // stagger[1], every)
    mfc = style.get("mfc", style["color"])
    return ax.plot(x, y, color=style["color"], ls=style.get("ls", "-"), lw=lw, marker=marker or "None", ms=ms,
                   markevery=every, mfc=mfc, mec=style["color"], mew=1.5, zorder=zorder)[0]


def _handle(style, lw=hs.CURVE_PT, ms=9, line=True):
    return Line2D([], [], color=style["color"], ls=style.get("ls", "-") if line else "none", lw=lw,
                  marker=style.get("marker") or "None", ms=ms, mfc=style.get("mfc", style["color"]),
                  mec=style["color"], mew=1.5)


def _ref_handle(ref):
    return Line2D([], [], **ref)


# --------------------------------------------------------------------------------------------------
# Panels shared by the figures
# --------------------------------------------------------------------------------------------------

def regime_flux_panel(P, J, frac_surface, frac_oxide=None, frac_metal=None, surface="oxide surface",
                      model="the model", what="", guides=True):
    """Flux against upstream pressure, the curve drawn in the colour of its limiting step, with the slope
    of each regime in the key and the slope-1 and slope-0.5 guides. Returns (Panel or None, words)."""
    P, J = np.asarray(P, float), np.asarray(J, float)
    ok = _ok(P, J)
    if ok.sum() < 2:
        return None, ""
    rl = regimes(frac_surface, frac_oxide, frac_metal)
    steps = _steps(surface)
    fig, ax = hs.new_figure()
    slopes = segment_slopes(P, J, rl)
    order = [r for r in ("surface", "oxide", "metal", "mixed") if np.any(ok & (rl == r))]
    handles, labels, colours = [], [], []
    idx = np.flatnonzero(ok)
    for regime in order:
        step = steps[regime]
        style = dict(color=hs.STEP_COLOURS[step], ls="-", marker=STEP_MARKERS[step])
        sel = idx[rl[idx] == regime]
        for _, g in groupby(enumerate(sel), lambda t: t[0] - t[1]):
            grp = [k for _, k in g]
            ext = grp + ([grp[-1] + 1] if grp[-1] + 1 < len(P) and ok[grp[-1] + 1] else [])
            ax.plot(P[ext], J[ext], color=style["color"], lw=hs.OUTLINE_PT, solid_capstyle="butt", zorder=3)
            marks = grp[:: max(1, len(idx) // 14)]
            ax.plot(P[marks], J[marks], ls="none", marker=style["marker"], ms=9, color=style["color"],
                    mec="white", mew=0.8, zorder=4)
        handles.append(_handle(style, lw=hs.OUTLINE_PT))
        lab = STEP_NAMES[step]
        if slopes.get(regime) and regime != "mixed":
            lab += f"\nslope {_slope_text(slopes[regime])}"
        labels.append(lab)
        colours.append(style["color"])
    _log_axes(ax, P[ok], J[ok])
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    if guides:
        Pv, Jv = P[ok], J[ok]
        ax.plot(Pv, Jv[0] * (Pv / Pv[0]) ** 1.0, zorder=2, **REF_DASHED)
        ax.plot(Pv, Jv[-1] * (Pv / Pv[-1]) ** 0.5, zorder=2, **REF_DOTTED)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        handles += [_ref_handle(REF_DASHED), _ref_handle(REF_DOTTED)]
        labels += ["Slope 1 guide", "Slope 0.5 guide"]
        colours += ["black", "black"]
    ax.set_xlabel(P_LABEL)
    ax.set_ylabel(rf"Flux $J$ ({J_UNIT})")
    hs.key_below(fig, handles, labels, ncol=3, colours=colours)
    ns = net_slope(P, J)
    _place_note(ax, [f"net slope {ns:.2f}"], [(P[ok], J[ok])])
    words = "; ".join(f"{STEP_WORDS[steps[r]]} slope {_slope_text(s)}" for r, s in slopes.items() if r != "mixed")
    caption = (f"Steady-state flux of {model} against upstream pressure{what}, drawn in the colour of the "
               f"rate-limiting step at each pressure (the step that holds more than half of the total "
               f"resistance). Fitted log--log slopes: {words or 'no regime spans three points'}; net slope "
               f"{ns:.2f}. Grey dashed and dotted lines: slope 1 (surface-limited) and slope 0.5 "
               f"(diffusion-limited) guides through the first and last points.")
    return Panel("flux", fig, ax, caption, P[ok], J[ok], keyed=True), words


def resistance_panel(P, fractions, crossovers=(), surface="oxide surface", model="the model", note=None):
    """Resistance fractions (%) against pressure, one line per step, with the 50 % threshold and the
    crossover pressures. fractions: {'surface'|'oxide'|'metal': array}; crossovers: [(name, P)]."""
    P = np.asarray(P, float)
    steps = _steps(surface)
    fig, ax = hs.new_figure()
    handles, labels, colours, xs, ys = [], [], [], [], []
    for key in ("surface", "oxide", "metal"):
        if key not in fractions or fractions[key] is None:
            continue
        f = np.asarray(fractions[key], float) * 100
        ok = np.isfinite(f) & np.isfinite(P) & (P > 0)
        if not ok.any():
            continue
        ls, marker = FRACTION_STYLE[key]
        style = dict(color=hs.STEP_COLOURS[steps[key]], ls=ls, marker=marker)
        _plot(ax, P[ok], f[ok], style)
        handles.append(_handle(style))
        labels.append(FRACTION_NAMES[key])
        colours.append(style["color"])
        xs.append(P[ok]); ys.append(f[ok])
    x_all = np.concatenate(xs) if xs else P[P > 0]
    _log_axes(ax, x_all, None, ylog=False)
    ax.set_ylim(-3, 105)
    ax.set_yticks(np.arange(0, 101, 20))
    hs.reference_line(ax, 50, **REF_DASHED)
    handles.append(_ref_handle(REF_DASHED)); labels.append("50 % threshold"); colours.append("black")
    shown = [(n, p) for n, p in crossovers if p is not None and np.isfinite(p)]
    for name, p in shown:
        hs.reference_line(ax, p, vertical=True, **REF_DOTTED)
    if shown:
        handles.append(_ref_handle(REF_DOTTED))
        labels.append("Crossover: " + ", ".join(f"{n} {_pa_key(p)}" for n, p in shown))
        colours.append("black")
    ax.set_xlabel(P_LABEL)
    ax.set_ylabel("Resistance fraction (%)")
    hs.key_below(fig, handles, labels, ncol=2, colours=colours)
    lines = list(note or [])
    if not shown:
        dominant = _dominant(fractions)
        lines.append(f"no crossover: {STEP_WORDS[steps[dominant]]} throughout" if dominant else "no crossover")
    guides = [(np.full(40, p), np.linspace(0, 100, 40)) for _, p in shown] + [(x_all, np.full(len(x_all), 50.0))]
    _place_note(ax, lines, list(zip(xs, ys)) + guides)
    cross_words = ("; ".join(f"{n} at {_pa(p)}" for n, p in shown) if shown
                   else "no crossover in the range")
    caption = (f"Share of the total resistance of {model} held by each step against upstream pressure: "
               + hs.tex_join([FRACTION_NAMES[k].lower() for k in ("surface", "oxide", "metal") if k in fractions])
               + f". The grey dashed line is the 50~\\% threshold; crossovers (grey dotted): {cross_words}.")
    return Panel("resistance", fig, ax, caption, np.concatenate(xs) if xs else None,
                 np.concatenate(ys) if ys else None, keyed=True)


def _dominant(fractions):
    best, key = -1, None
    for k, f in fractions.items():
        if f is None:
            continue
        n = int(np.nansum(np.asarray(f, float) > 0.5))
        if n > best:
            best, key = n, k
    return key if best > 0 else None


def parity_panel(series, x_label, y_label, note, words):
    """A parity plot: each series' model flux against the reference flux it must recover, with the
    perfect-agreement line. series: [(x_ref, y_model, label, style)]."""
    fig, ax = hs.new_figure()
    xs, ys, handles, labels, colours = [], [], [], [], []
    several = len(series) > 1
    for k, (x, y, label, style) in enumerate(series):
        x, y = np.asarray(x, float), np.asarray(y, float)
        ok = _ok(x, y)
        if not ok.any():
            continue
        if several:                       # series that coincide: large open markers under small filled ones
            style = {**style, "mfc": "white" if k == 0 else style["color"]}
        ms = 14 if several and k == 0 else (7 if several else 9)
        ax.plot(x[ok], y[ok], ls="none", marker=style["marker"], ms=ms, color=style["color"],
                mfc=style.get("mfc", style["color"]), mec=style["color"], mew=2.0, zorder=4 + k)
        xs.append(x[ok]); ys.append(y[ok])
        handles.append(_handle(style, line=False, ms=11 if several and k == 0 else 8)); labels.append(label)
        colours.append(style["color"])
    if xs:
        allv = np.concatenate(xs + ys)
        lo, hi = hs.log_limits(allv, 0.3, frac=0.04)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        hs.log_ticks(ax.xaxis, lo, hi); hs.log_ticks(ax.yaxis, lo, hi)
        ax.plot([lo, hi], [lo, hi], zorder=2, **REF_DASHED)
    handles.append(_ref_handle(REF_DASHED)); labels.append("Perfect agreement"); colours.append("black")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    hs.key_below(fig, handles, labels, ncol=2, colours=colours)
    _place_note(ax, note, list(zip(xs, ys)))
    caption = words + " The grey dashed line is perfect agreement."
    return Panel("limit", fig, ax, caption, np.concatenate(xs) if xs else None,
                 np.concatenate(ys) if ys else None, keyed=True)


def _limit_words(max_err, mean_err=None):
    if max_err is None or not np.isfinite(max_err):
        return "No valid points: the check failed."
    words = f"Largest relative difference {_tex_pct(max_err)}"
    if mean_err is not None and np.isfinite(mean_err):
        words += f", mean {_tex_pct(mean_err)}"
    return words + (": the check passes (under 1~\\%)." if max_err < 1.0 else
                    ": the check fails (above 1~\\%), as the notebook flags.")


def _errors(x_ref, y_model):
    x, y = np.asarray(x_ref, float), np.asarray(y_model, float)
    ok = _ok(x, y)
    if not ok.any():
        return np.nan, np.nan
    e = np.abs(y[ok] - x[ok]) / x[ok] * 100
    return float(e.max()), float(e.mean())


def _conditions(T_C, extra=""):
    return f" at {T_C:.0f}~\\textdegree C{extra}"


# --------------------------------------------------------------------------------------------------
# L1+L6 and L2a+L6: one surface and one diffusion resistance
# --------------------------------------------------------------------------------------------------

def two_resistance_panels(P, J, theta, frac_surface, frac_diffusion, K_eq, P_crossover, J_reference, J_fast,
                          layer, k_fast, T_C, material=""):
    """The four figures of L1+L6 (layer='metal': a bare metal surface over a perfect metal) or L2a+L6
    (layer='oxide': an oxide surface over a perfect oxide): flux, coverage, resistance, limit."""
    level = "L1L6" if layer == "metal" else "L2aL6"
    name = "L1+L6" if layer == "metal" else "L2a+L6"
    surface = "metal surface" if layer == "metal" else "oxide surface"
    model = f"{name} (a perfect {'metal' if layer == 'metal' else 'oxide'} with surface kinetics)"
    cond = _conditions(T_C)
    P = np.asarray(P, float)
    fr = {"surface": frac_surface, layer: frac_diffusion}
    panels = []
    flux, _ = regime_flux_panel(P, J, frac_surface, frac_diffusion if layer == "oxide" else None,
                                frac_diffusion if layer == "metal" else None, surface=surface, model=model, what=cond)
    panels.append(flux)
    panels.append(_coverage_panel(P, theta, K_eq, P_crossover, level, name, cond))
    cross = [("surface → " + layer, P_crossover)]
    panels.append(resistance_panel(P, fr, cross, surface=surface, model=f"{name}{cond}"))
    mx, mn = _errors(J_reference, J_fast)
    ref_name = "Sieverts (L1)" if layer == "metal" else "L2a"
    style = _style(level)
    panels.append(parity_panel(
        [(J_reference, J_fast, f"{name}, $k_\\mathrm{{diss}}$ = {hs.sci(k_fast, 1)}", style)],
        rf"{ref_name} flux ({J_UNIT})", rf"{name} flux ({J_UNIT})",
        [f"$k_\\mathrm{{diss}} \\to \\infty$ ({hs.sci(k_fast, 1)})"] + _check_lines(mx, mn),
        f"Limit check{cond}: with fast dissociation ($k_\\mathrm{{diss}}$ = {hs.tex_num(k_fast, 1)}) the {name} "
        f"flux must recover the {ref_name} flux without surface kinetics, "
        + ("$J = \\beta(\\sqrt{P_\\mathrm{up}} - \\sqrt{P_\\mathrm{down}})$" if layer == "metal" else
           "$J = \\alpha(\\sqrt{P_\\mathrm{up}} - \\sqrt{P_\\mathrm{down}})$")
        + f", at every pressure. {_limit_words(mx, mn)}"))
    return [p for p in panels if p is not None]


def _coverage_panel(P, theta, K_eq, P_crossover, level, name, cond):
    P, theta = np.asarray(P, float), np.asarray(theta, float)
    ok = np.isfinite(P) & np.isfinite(theta) & (P > 0)
    if not ok.any():
        return None
    fig, ax = hs.new_figure()
    style = _style(level)
    _plot(ax, P[ok], theta[ok], style)
    theta_eq = np.sqrt(K_eq * P) / (1 + np.sqrt(K_eq * P))
    _log_axes(ax, P[ok], None, ylog=False)
    ax.plot(P, theta_eq, zorder=2, **REF_DASHDOT)
    ax.set_ylim(-0.04, 1.06)
    handles = [_handle(style), _ref_handle(REF_DASHDOT)]
    labels = [r"$\theta$ at steady state", r"$\theta_\mathrm{eq}$, Langmuir isotherm"]
    colours = [style["color"], "black"]
    if P_crossover is not None and np.isfinite(P_crossover):
        hs.reference_line(ax, P_crossover, vertical=True, **REF_DOTTED)
        handles.append(_ref_handle(REF_DOTTED)); labels.append(f"Crossover, {_pa_key(P_crossover)}")
        colours.append("black")
    ax.set_xlabel(P_LABEL)
    ax.set_ylabel(r"Surface coverage $\theta$")
    hs.key_below(fig, handles, labels, ncol=3, colours=colours)
    gap = np.nanmax(theta_eq[ok] - theta[ok])
    caption = (f"Surface coverage of {name}{cond} against upstream pressure: the steady-state coverage and the "
               r"Langmuir isotherm $\theta_\mathrm{eq} = \sqrt{K_\mathrm{eq}P}/(1+\sqrt{K_\mathrm{eq}P})$ "
               f"(grey dash-dot, $K_\\mathrm{{eq}}$ = {hs.tex_num(K_eq, 3)}~Pa$^{{-1}}$). Where the two coincide the "
               "surface is in local equilibrium (diffusion-limited); below it, fast diffusion depletes the surface. "
               f"The coverage runs from {hs.tex_num(np.nanmin(theta[ok]), 2)} to {hs.tex_num(np.nanmax(theta[ok]), 2)}; "
               f"the largest shortfall below the isotherm is {hs.tex_num(gap, 2)}."
               + (f" Crossover (grey dotted) at {_pa(P_crossover)}." if P_crossover is not None and
                  np.isfinite(P_crossover) else ""))
    return Panel("coverage", fig, ax, caption, P[ok], theta[ok], keyed=True)


# --------------------------------------------------------------------------------------------------
# L2+L6: oxide surface, oxide and metal
# --------------------------------------------------------------------------------------------------

def l2l6_panels(P, J, frac_surface, frac_oxide, frac_metal, crossovers, arrhenius, J_reference, J_fast,
                k_fast, T_C):
    """L2+L6 (perfect oxide on a perfect metal, oxide-surface kinetics): flux, arrhenius, resistance, limit.
    arrhenius: {'inv_T': 1000/T, 'cases': [{'regime', 'P_fixed', 'J', 'classes', 'fit', 'E_ref'}]}, fit a
    dict with 'slope', 'intercept', 'r_sq', 'E_a' (J/mol), 'n_pure' and 'T_pure_range' (°C), or None."""
    name, cond = "L2+L6", _conditions(T_C)
    model = f"{name} (a perfect oxide on a perfect metal, with oxide-surface kinetics)"
    panels = []
    flux, _ = regime_flux_panel(P, J, frac_surface, frac_oxide, frac_metal, model=model, what=cond)
    panels.append(flux)
    panels.append(arrhenius_panel(arrhenius, name))
    panels.append(resistance_panel(P, {"surface": frac_surface, "oxide": frac_oxide, "metal": frac_metal},
                                   crossovers, model=f"{name}{cond}"))
    mx, mn = _errors(J_reference, J_fast)
    panels.append(parity_panel(
        [(J_reference, J_fast, f"{name}, $k_\\mathrm{{diss}}$ = {hs.sci(k_fast, 1)}", _style("L2bL6"))],
        rf"Series-resistance flux ({J_UNIT})", rf"{name} flux ({J_UNIT})",
        [f"$k_\\mathrm{{diss}} \\to \\infty$ ({hs.sci(k_fast, 1)})"] + _check_lines(mx, mn),
        f"Limit check{cond}: with fast dissociation ($k_\\mathrm{{diss}}$ = {hs.tex_num(k_fast, 1)}) the {name} flux "
        r"must recover the series-resistance flux $J = (\sqrt{P_\mathrm{up}} - \sqrt{P_\mathrm{down}})/(1/\alpha + "
        f"1/\\beta)$ at every pressure. {_limit_words(mx, mn)}"))
    return [p for p in panels if p is not None]


ARRHENIUS_CLASSES = {   # class -> (step regime, filled?)
    "pure_surface": ("surface", True), "mixed_surface": ("surface", False),
    "pure_oxide": ("oxide", True), "mixed_oxide": ("oxide", False),
    "pure_metal": ("metal", True), "mixed_metal": ("metal", False), "mixed": ("mixed", True),
}


def arrhenius_panel(arrhenius, name, surface="oxide surface"):
    """Flux against 1000/T at one fixed pressure per regime; each point in the colour of the step limiting
    it (filled: that step holds 99 % of the resistance, open: over half), the Arrhenius fit through the
    filled points dashed, and each fit's E_a against the expected value in the key."""
    inv_T = np.asarray(arrhenius["inv_T"], float)
    cases = [c for c in arrhenius["cases"] if c.get("J") is not None]
    if not cases:
        return None
    steps = _steps(surface)
    fig, ax = hs.new_figure()
    handles, labels, colours, xs, ys = [], [], [], [], []
    open_used = mixed_used = False
    words = []
    for c in cases:
        Jt = np.asarray(c["J"], float)
        classes = np.asarray(c["classes"])
        step = steps[c["regime"]]
        colour = hs.STEP_COLOURS[step]
        for cls, (reg, filled) in ARRHENIUS_CLASSES.items():
            m = (classes == cls) & _ok(Jt)
            if not m.any():
                continue
            st = steps[reg]
            col = hs.STEP_COLOURS[st]
            ax.plot(inv_T[m], Jt[m], ls="none", marker=STEP_MARKERS[st], ms=10 if filled else 8,
                    color=col, mfc=col if filled else "white", mec=col, mew=2.0, zorder=4)
            xs.append(inv_T[m]); ys.append(Jt[m])
            open_used |= not filled
            mixed_used |= reg == "mixed"
        fit = c.get("fit")
        lab = f"{STEP_NAMES[step]}, {_pa_key(c['P_fixed'])}"
        if fit:
            pure = np.asarray(fit.get("valid", (classes == f"pure_{c['regime']}") & _ok(Jt)))
            xf = inv_T[pure]
            ax.plot(xf, np.exp(fit["slope"] * xf + fit["intercept"]), color=colour, ls=(0, (5, 3)), lw=2.6,
                    zorder=3)
            lab += (f"\n$E_a$ = {hs.kj(fit['E_a'] / 1000, 1)} (expected {fit_ref(c)}), $R^2$ = {fit['r_sq']:.3f}")
            words.append(f"{STEP_WORDS[step]} at {_pa(c['P_fixed'])}: $E_a$ = {hs.tex_kj(fit['E_a'] / 1000, 1)} "
                         f"against {hs.tex_kj(c['E_ref'] / 1000, 1)} expected ($R^2$ = {fit['r_sq']:.4f}, "
                         f"{fit.get('n_pure', int(pure.sum()))} points)")
        else:
            lab += "\nnot extractable (under 3 pure points)"
            words.append(f"{STEP_WORDS[step]} at {_pa(c['P_fixed'])}: not extractable")
        handles.append(Line2D([], [], color=colour, ls=(0, (5, 3)) if fit else "none", lw=2.6,
                              marker=STEP_MARKERS[step], ms=9, mfc=colour, mec=colour))
        labels.append(lab)
        colours.append(colour)
    for c in arrhenius["cases"]:
        if c.get("J") is None:
            step = steps[c["regime"]]
            handles.append(Line2D([], [], ls="none", marker=STEP_MARKERS[step], ms=9, color=hs.STEP_COLOURS[step]))
            labels.append(f"{STEP_NAMES[step]}: no pressure where it dominates")
            colours.append(hs.STEP_COLOURS[step])
            words.append(f"{STEP_WORDS[step]}: no pressure where it holds 90~\\% of the resistance")
    if open_used:
        handles.append(Line2D([], [], ls="none", marker="o", ms=9, mfc="white", mec="black", mew=2.0))
        labels.append("Open: over half, not 99 %, of the resistance"); colours.append("black")
    x_all, y_all = np.concatenate(xs), np.concatenate(ys)
    pad = 0.05 * (x_all.max() - x_all.min())
    ax.set_xlim(x_all.min() - pad, x_all.max() + pad)
    hs.set_log_y(ax, y_all)
    ax.set_xlabel(r"1000/$T$ (K$^{-1}$)")
    ax.set_ylabel(rf"Flux $J$ ({J_UNIT})")
    top = hs.celsius_top_axis(ax)
    hs.set_celsius_ticks(ax, top)
    hs.key_below(fig, handles, labels, ncol=1, colours=colours)
    T_lo, T_hi = 1000 / x_all.max() - 273.15, 1000 / x_all.min() - 273.15
    caption = (f"Arrhenius plot of {name}: flux against inverse temperature ({T_lo:.0f} to {T_hi:.0f}~\\textdegree C) "
               "at one fixed pressure per regime, the median pressure of the range where that step holds over 90~\\% "
               "of the resistance. Each point is in the colour of the step limiting it there: filled where it holds "
               "99~\\% or more, open where it holds over half"
               + (", grey crosses where no step does" if mixed_used else "")
               + ". Dashed lines: Arrhenius fits through the filled points. " + "; ".join(words) + ".")
    return Panel("arrhenius", fig, ax, caption, x_all, y_all, keyed=True)


def fit_ref(case):
    return hs.kj(case["E_ref"] / 1000, 1).replace(" kJ mol$^{-1}$", "")


# --------------------------------------------------------------------------------------------------
# L3+L6: defective oxide
# --------------------------------------------------------------------------------------------------

def l3l6_panels(P, curves, slope_baseline, frac_surface, frac_oxide, frac_metal, crossovers, defects,
                J_reference, J_limit, f_limit, T_C, total_fraction):
    """L3+L6 (defective oxide on a perfect metal, surface kinetics): flux, defects, resistance, limit.
    curves: [(kind, f, J)] with kind 'baseline' (f = 0), 'fraction' or 'full' (f = 1, the dominant defect type
    in place of the whole oxide; f is then that type's name); defects: {'P', 'f' (fractions), 'J_ref',
    'types': [(type, J array)], 'f_cross': {type: f or None}}."""
    name, cond = "L3+L6", _conditions(T_C)
    panels = [_fraction_flux_panel(P, curves, slope_baseline, cond, total_fraction),
              _defect_panel(defects, cond)]
    panels.append(resistance_panel(
        P, {"surface": frac_surface, "oxide": frac_oxide, "metal": frac_metal}, crossovers,
        model=f"{name}{cond}, in the path that carries the most flux at each pressure",
        note=["dominant path"]))
    mx, mn = _errors(J_reference, J_limit)
    panels.append(parity_panel(
        [(J_reference, J_limit, f"{name}, every defect fraction {hs.sci(f_limit, 1)}", _style("L3L6"))],
        rf"L2+L6 flux ({J_UNIT})", rf"{name} flux ({J_UNIT})",
        [f"$f \\to 0$ ({hs.sci(f_limit, 1)}, all paths)"] + _check_lines(mx, mn),
        f"Limit check{cond}: with every defect fraction set to {hs.tex_num(f_limit, 1)} the {name} flux must recover "
        f"the L2+L6 flux of the intact oxide at every pressure. {_limit_words(mx, mn)}"))
    return [p for p in panels if p is not None]


def _fraction_flux_panel(P, curves, slope_baseline, cond, total_fraction):
    P = np.asarray(P, float)
    fig, ax = hs.new_figure()
    handles, labels, colours, xs, ys = [], [], [], [], []
    k = 0
    words, drawn, names = [], [], []
    n = sum(_ok(P, np.asarray(J, float)).any() for _, _, J in curves)
    for kind, f, J in curves:
        J = np.asarray(J, float)
        ok = _ok(P, J)
        if not ok.any():
            continue
        if kind == "baseline":
            style = _style("L2bL6")
            lab = "f = 0, intact oxide (L2+L6)"
            if np.isfinite(slope_baseline):
                lab += f"\nnet slope {slope_baseline:.2f}"
        elif kind == "full":
            style = dict(color="black", ls=DOTTED, marker="x")
            lab = f"f = 1, all {str(f).replace('_', ' ')}"
        else:
            ls, marker = CONDITION_MARKS[k % len(CONDITION_MARKS)]
            k += 1
            style = _style("L3", ls=ls, marker=marker)
            lab = f"f = {hs.sci(100 * f, 2)} %"
        _plot(ax, P[ok], J[ok], style, stagger=(len(drawn), n))
        handles.append(_handle(style)); labels.append(lab); colours.append(style["color"])
        xs.append(P[ok]); ys.append(J[ok])
        words.append(lab.split("\n")[0].replace("%", r"\%").replace("f = ", "$f$ = "))
        drawn.append(J); names.append(words[-1])
    _log_axes(ax, np.concatenate(xs), np.concatenate(ys))
    ax.set_xlabel(P_LABEL)
    ax.set_ylabel(rf"Flux $J$ ({J_UNIT})")
    hs.key_below(fig, handles, labels, ncol=3, colours=colours)
    caption = (f"Flux of L3+L6 (defective oxide on a perfect metal, with surface kinetics){cond} against upstream "
               f"pressure for several total defect fractions (the configured mix of defect types, "
               f"{hs.tex_num(100 * total_fraction, 2)}~\\% in total, scaled): " + hs.tex_join(words)
               + f". Net log--log slope of the intact-oxide curve: {slope_baseline:.2f}."
               + (coincide_words(P, drawn, names) if drawn else ""))
    return Panel("flux", fig, ax, caption, np.concatenate(xs), np.concatenate(ys), keyed=True)


DEFECT_NAMES = {"pinhole": "Pinholes", "crack": "Cracks", "grain_boundary": "Oxide grain boundaries"}


def _defect_panel(defects, cond):
    f = np.asarray(defects["f"], float) * 100
    fig, ax = hs.new_figure()
    handles, labels, colours, xs, ys = [], [], [], [], []
    words = []
    for k, (dtype, J) in enumerate(defects["types"]):
        J = np.asarray(J, float)
        ok = _ok(f, J)
        if not ok.any():
            continue
        ls, marker = CONDITION_MARKS[k % len(CONDITION_MARKS)]
        style = dict(color="black", ls=ls, marker=marker)
        _plot(ax, f[ok], J[ok], style, stagger=(k, len(defects["types"])))
        fc = defects.get("f_cross", {}).get(dtype)
        lab = DEFECT_NAMES.get(dtype, dtype)
        lab += f"\nabove the baseline from {hs.sci(100 * fc, 2)} %" if fc is not None else "\nnever above the baseline"
        handles.append(_handle(style)); labels.append(lab); colours.append("black")
        xs.append(f[ok]); ys.append(J[ok])
        words.append(f"{DEFECT_NAMES.get(dtype, dtype).lower()} "
                     + (f"from {hs.tex_num(100 * fc, 2)}~\\%" if fc is not None else "never"))
    x_all, y_all = np.concatenate(xs), np.concatenate(ys)
    J_ref = defects.get("J_ref")
    _log_axes(ax, x_all, np.r_[y_all, J_ref] if J_ref is not None and np.isfinite(J_ref) else y_all)
    if J_ref is not None and np.isfinite(J_ref):
        hs.reference_line(ax, J_ref, **REF_DASHED)
        handles.append(_ref_handle(REF_DASHED)); labels.append("L2+L6 baseline, f = 0"); colours.append("black")
    ax.set_xlabel("Defect fraction (%)")
    ax.set_ylabel(rf"Flux $J$ ({J_UNIT})")
    cap = 100.0
    lo, hi = ax.get_xlim()
    hs.log_ticks(ax.xaxis, lo, min(hi, cap))
    hs.key_below(fig, handles, labels, ncol=2, colours=colours)
    caption = (f"Flux of L3+L6 against the area fraction of one defect type{cond} and {_pa(defects['P'])}, each "
               "type alone (the others set to zero): " + hs.tex_join(DEFECT_NAMES.get(t, t).lower()
                                                                     for t, _ in defects["types"])
               + ". The grey dashed line is the L2+L6 flux of the intact oxide"
               + (f", {hs.tex_num(J_ref, 3)}~{J_UNIT.replace(' ', '~')}" if J_ref is not None and np.isfinite(J_ref) else "")
               + ". The flux rises above it: " + "; ".join(words) + ".")
    return Panel("defects", fig, ax, caption, x_all, y_all, keyed=True)


# --------------------------------------------------------------------------------------------------
# L2+L4+L6: perfect oxide on a defective metal
# --------------------------------------------------------------------------------------------------

MODE_STYLE = {"none": dict(hs.level_style("L2bL6")), "gb_only": hs.level_style("L4_gb"),
              "trapping_only": hs.level_style("L4_trap"), "both": hs.level_style("L4_both")}
MODE_NAMES = {"none": "mode none (L2+L6)", "gb_only": "Grain boundaries only", "trapping_only": "Trapping only",
              "both": "Both (traps and grain boundaries)"}


def l2l4l6_panels(P, mode_curves, mode, slope_net, frac_surface, frac_oxide, frac_metal, crossovers,
                  diffusivity, D_ratio, J_reference, J_none, T_C, grain_um):
    """L2+L4+L6 (perfect oxide on a defective metal, surface kinetics): flux, diffusivity, resistance, limit.
    mode_curves: [(mode, J, frac_surface, frac_oxide, frac_metal)] on P; diffusivity: {'T_C', 'ratio':
    {mode: array}, 'T_cross': {mode: °C or None}, 'T_op_C'}; D_ratio: D_eff/D_m at the operating mode."""
    name, cond = "L2+L4+L6", _conditions(T_C, f", grain size {grain_um:.0f}~\\textmu m")
    P = np.asarray(P, float)
    fig, ax = hs.new_figure()
    handles, labels, colours, xs, ys = [], [], [], [], []
    words, drawn, names = [], [], []
    n = sum(_ok(P, np.asarray(c[1], float)).any() for c in mode_curves)
    for m, J, fs, fo, fm in mode_curves:
        J = np.asarray(J, float)
        ok = _ok(P, J)
        if not ok.any():
            continue
        style = MODE_STYLE.get(m, _style("L4"))
        _plot(ax, P[ok], J[ok], style, stagger=(len(drawn), n))
        drawn.append(J); names.append(MODE_NAMES.get(m, m).lower().replace("mode none (l2+l6)", "mode none"))
        lab = MODE_NAMES.get(m, m)
        if m == mode:
            line = _regime_slope_line(P, J, regimes(fs, fo, fm), "oxide surface")
            lab += (f"\n{line}" if line else "") + (f"\nnet slope {slope_net:.2f}" if np.isfinite(slope_net) else "")
            words.append(f"for the operating mode ({mode.replace('_', ' ')}): {line}, net slope {slope_net:.2f}")
        handles.append(_handle(style)); labels.append(lab); colours.append(style["color"])
        xs.append(P[ok]); ys.append(J[ok])
    _log_axes(ax, np.concatenate(xs), np.concatenate(ys))
    ax.set_xlabel(P_LABEL)
    ax.set_ylabel(rf"Flux $J$ ({J_UNIT})")
    hs.key_below(fig, handles, labels, ncol=2, colours=colours)
    flux = Panel("flux", fig, ax,
                 f"Flux of {name} (a perfect oxide on a defective metal, with oxide-surface kinetics){cond} against "
                 "upstream pressure for each microstructure mode of the metal: none (the L2+L6 baseline), grain "
                 "boundaries only, trapping only, and both. Fitted log--log slopes " + "; ".join(words) + "."
                 + coincide_words(P, drawn, names),
                 np.concatenate(xs), np.concatenate(ys), keyed=True)
    panels = [flux, _diffusivity_panel(diffusivity, name)]
    panels.append(resistance_panel(
        P, {"surface": frac_surface, "oxide": frac_oxide, "metal": frac_metal}, crossovers,
        model=f"{name}{cond}, mode {mode.replace('_', ' ')}",
        note=[f"mode {mode.replace('_', ' ')}", f"$D_\\mathrm{{eff}}/D_m$ = {D_ratio:.4f}"]))
    mx, mn = _errors(J_reference, J_none)
    panels.append(parity_panel(
        [(J_reference, J_none, f"{name}, mode none", _style("L2bL6"))],
        rf"L2+L6 flux ({J_UNIT})", rf"{name} flux ({J_UNIT})",
        ["mode none"] + _check_lines(mx, mn),
        f"Limit check{cond}: with the microstructure switched off (mode none) the {name} flux must recover the "
        f"L2+L6 flux at every pressure. {_limit_words(mx, mn)}"))
    return [p for p in panels if p is not None]


def _diffusivity_panel(d, name):
    T = np.asarray(d["T_C"], float)
    fig, ax = hs.new_figure()
    handles, labels, colours, xs, ys = [], [], [], [], []
    words = []
    for m, ratio in d["ratio"].items():
        r = np.asarray(ratio, float)
        ok = np.isfinite(r) & (r > 0)
        if not ok.any():
            continue
        style = MODE_STYLE.get(m, _style("L4"))
        _plot(ax, T[ok], r[ok], style, stagger=(list(d["ratio"]).index(m), len(d["ratio"])))
        lab = MODE_NAMES.get(m, m).replace("mode none (L2+L6)", "mode none ($D_\\mathrm{eff} = D_m$)")
        tc = d.get("T_cross", {}).get(m)
        if tc is not None and np.isfinite(tc):
            lab += f"\nequals $D_m$ at {tc:.0f} °C"
            words.append(f"{MODE_NAMES.get(m, m).lower()} equals the lattice value at {tc:.0f}~\\textdegree C")
        handles.append(_handle(style)); labels.append(lab); colours.append(style["color"])
        xs.append(T[ok]); ys.append(r[ok])
    x_all, y_all = np.concatenate(xs), np.concatenate(ys)
    pad = 0.04 * (x_all.max() - x_all.min())
    ax.set_xlim(x_all.min() - pad, x_all.max() + pad)
    hs.set_log_y(ax, np.r_[y_all, 1.0])
    hs.reference_line(ax, 1.0, **REF_DASHED)
    handles.append(_ref_handle(REF_DASHED)); labels.append(r"$D_\mathrm{eff} = D_m$ (lattice)"); colours.append("black")
    T_op = d.get("T_op_C")
    if T_op is not None:
        hs.reference_line(ax, T_op, vertical=True, **REF_DOTTED)
        handles.append(_ref_handle(REF_DOTTED)); labels.append(f"Operating temperature, {T_op:.0f} °C")
        colours.append("black")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel(r"$D_\mathrm{eff}/D_m$")
    hs.key_below(fig, handles, labels, ncol=2, colours=colours)
    lo = {m: np.nanmin(np.asarray(r, float)) for m, r in d["ratio"].items()}
    caption = (f"Effective diffusivity of the defective metal relative to the lattice, $D_\\mathrm{{eff}}/D_m$, against "
               f"temperature ({x_all.min():.0f} to {x_all.max():.0f}~\\textdegree C) for each microstructure mode of "
               f"{name}. Grey dashed: the lattice value; grey dotted: the operating temperature. The lowest values are "
               + hs.tex_join(f"{hs.tex_num(v, 2)} ({MODE_NAMES.get(m, m).lower()})" for m, v in lo.items() if np.isfinite(v))
               + (". " + "; ".join(words) + "." if words else "."))
    return Panel("diffusivity", fig, ax, caption, x_all, y_all, keyed=True)


# --------------------------------------------------------------------------------------------------
# L3+L4+L6: the full model
# --------------------------------------------------------------------------------------------------

HIERARCHY = {   # case -> (key name, style): +L6 levels in their base level's hue, each its own line and marker
    "L2L6": ("L2+L6: f = 0, mode none", _style("L2bL6", ls="--", marker="s")),
    "L3L6": ("L3+L6: configured f, mode none", _style("L3L6", ls=DASHDOT, marker="D")),
    "L2L4L6": ("L2+L4+L6: f = 0, operating mode", _style("L4L6", ls=DOTTED, marker="^")),
    "L3L4L6": ("L3+L4+L6: both active", _style("L5L6", ls="-", marker="o")),
}


def l3l4l6_panels(P, hierarchy, fractions_full, slope_net, frac_surface, frac_oxide, frac_metal, crossovers,
                  n_unresolved, J_limit1, J_ref1, J_limit2, J_ref2, mode, T_C, total_fraction):
    """L3+L4+L6 (defective oxide on a defective metal, surface kinetics): flux, enhancement, resistance, limit.
    hierarchy: {case: J} on P for 'L2L6', 'L3L6', 'L2L4L6', 'L3L4L6'; fractions_full: (fs, fo, fm) of the full
    model on P, for its regime slopes; limit 1: f -> 0 against L2+L4+L6; limit 2: mode none against L3+L6."""
    name = "L3+L4+L6"
    cond = _conditions(T_C, f", mode {mode.replace('_', ' ')}, defect fraction {hs.tex_num(100 * total_fraction, 2)}~\\%")
    P = np.asarray(P, float)
    fig, ax = hs.new_figure()
    handles, labels, colours, xs, ys = [], [], [], [], []
    slope_words = ""
    drawn, names = [], []
    n = sum(_ok(P, np.asarray(hierarchy.get(c, np.full(len(P), np.nan)), float)).any() for c in HIERARCHY)
    for case, (lab, style) in HIERARCHY.items():
        J = np.asarray(hierarchy.get(case, np.full(len(P), np.nan)), float)
        ok = _ok(P, J)
        if not ok.any():
            continue
        _plot(ax, P[ok], J[ok], style, stagger=(len(drawn), n))
        drawn.append(J); names.append(lab.split(":")[0])
        if case == "L3L4L6":
            slope_words = _regime_slope_line(P, J, regimes(*fractions_full), "oxide surface")
            lab += (f"\n{slope_words}" if slope_words else "") + \
                   (f"\nnet slope {slope_net:.2f}" if np.isfinite(slope_net) else "")
        handles.append(_handle(style)); labels.append(lab); colours.append(style["color"])
        xs.append(P[ok]); ys.append(J[ok])
    _log_axes(ax, np.concatenate(xs), np.concatenate(ys))
    ax.set_xlabel(P_LABEL)
    ax.set_ylabel(rf"Flux $J$ ({J_UNIT})")
    hs.key_below(fig, handles, labels, ncol=2, colours=colours)
    flux = Panel("flux", fig, ax,
                 f"Flux against upstream pressure{cond} through the hierarchy of surface-kinetics levels: L2+L6 "
                 "(intact oxide, perfect metal), L3+L6 (defective oxide), L2+L4+L6 (defective metal) and L3+L4+L6 "
                 f"(both). Fitted log--log slopes of the full model: {slope_words or 'none'}; net slope "
                 f"{slope_net:.2f}." + coincide_words(P, drawn, names), np.concatenate(xs), np.concatenate(ys),
                 keyed=True)
    panels = [flux, _enhancement_panel(P, hierarchy, cond)]
    note = [f"dominant path, mode {mode.replace('_', ' ')}"]
    if n_unresolved:
        note.append(f"{n_unresolved} points with no dominant path")
    panels.append(resistance_panel(
        P, {"surface": frac_surface, "oxide": frac_oxide, "metal": frac_metal}, crossovers,
        model=f"{name}{cond}, flux-weighted over the paths", note=note))
    mx1, mn1 = _errors(J_ref1, J_limit1)
    mx2, mn2 = _errors(J_ref2, J_limit2)
    panels.append(parity_panel(
        [(J_ref1, J_limit1, r"$f \to 0$, against L2+L4+L6", _style("L4", marker="o")),
         (J_ref2, J_limit2, "mode none, against L3+L6", _style("L3", marker="s"))],
        rf"Reference flux ({J_UNIT})", rf"{name} flux ({J_UNIT})",
        _check_lines(mx1, None, r"$f \to 0$")[:1] + _check_lines(mx2, None, "mode none")[:1]
        + (["both checks pass (under 1 %)"] if np.isfinite(mx1) and np.isfinite(mx2) and max(mx1, mx2) < 1 else
           ["check failed: a max error above 1 %"]),
        f"Two limit checks{cond}: with every defect fraction near zero the {name} flux must recover L2+L4+L6 "
        f"(largest difference {_tex_pct(mx1)}), and with the microstructure switched off (mode none) it must "
        f"recover L3+L6 (largest difference {_tex_pct(mx2)})."
        + (" Both pass (under 1~\\%)." if np.isfinite(mx1) and np.isfinite(mx2) and max(mx1, mx2) < 1 else
           " At least one fails (above 1~\\%).")))
    return [p for p in panels if p is not None]


def _enhancement_panel(P, hierarchy, cond):
    P = np.asarray(P, float)
    base = np.asarray(hierarchy.get("L2L6"), float)
    ok_b = _ok(base)
    E = {}
    for key, case in (("oxide", "L3L6"), ("metal", "L2L4L6"), ("both", "L3L4L6")):
        J = np.asarray(hierarchy.get(case, np.full(len(P), np.nan)), float)
        e = np.full(len(P), np.nan)
        e[ok_b] = J[ok_b] / base[ok_b]
        E[key] = e
    E["product"] = E["oxide"] * E["metal"]
    styles = {"oxide": (_style("L3", ls="--", marker="D"), r"$E_\mathrm{oxide}$ = L3+L6 / L2+L6"),
              "metal": (_style("L4", ls=DASHDOT, marker="^"), r"$E_\mathrm{metal}$ = L2+L4+L6 / L2+L6"),
              "both": (_style("L5", ls="-", marker="o"), r"$E_\mathrm{both}$ = L3+L4+L6 / L2+L6"),
              "product": (dict(color="black", ls=DOTTED, marker="x"), r"$E_\mathrm{oxide} \times E_\mathrm{metal}$")}
    fig, ax = hs.new_figure()
    handles, labels, colours, xs, ys = [], [], [], [], []
    for key, (style, lab) in styles.items():
        e = E[key]
        ok = np.isfinite(e) & (e > 0) & (P > 0)
        if not ok.any():
            continue
        _plot(ax, P[ok], e[ok], style)
        handles.append(_handle(style)); labels.append(lab); colours.append(style["color"])
        xs.append(P[ok]); ys.append(e[ok])
    x_all, y_all = np.concatenate(xs), np.concatenate(ys)
    _log_axes(ax, x_all, None, ylog=False)
    hs.set_linear_y(ax, y_all, include=[1.0], frac=0.1)
    hs.reference_line(ax, 1.0, **REF_DASHED)
    handles.append(_ref_handle(REF_DASHED)); labels.append("E = 1, no enhancement"); colours.append("black")
    ax.set_xlabel(P_LABEL)
    ax.set_ylabel(r"Flux ratio $E$")
    hs.key_below(fig, handles, labels, ncol=2, colours=colours)
    both, prod = E["both"], E["product"]
    m = np.isfinite(both) & np.isfinite(prod) & (prod > 0)
    behaviour, words = None, ""
    if m.any():
        ratio = np.nanmedian(both[m] / prod[m])
        if ratio > 0.8:
            behaviour = r"multiplicative: $E_\mathrm{both} \approx E_\mathrm{oxide} \times E_\mathrm{metal}$"
            words = "multiplicative"
        elif np.nanmedian(both[m]) > np.nanmedian(np.maximum(E["oxide"][m], E["metal"][m])):
            behaviour, words = "super-multiplicative", "super-multiplicative"
        else:
            behaviour, words = "sub-multiplicative: one effect dominates", "sub-multiplicative"
        _place_note(ax, [behaviour], list(zip(xs, ys)))
    rng = lambda e: f"{hs.tex_num(np.nanmin(e), 3)} to {hs.tex_num(np.nanmax(e), 3)}"   # noqa: E731
    caption = (f"Flux enhancement against upstream pressure{cond}: each defective level divided by the L2+L6 flux "
               f"of the intact system, $E_\\mathrm{{oxide}}$ ({rng(E['oxide'])}), $E_\\mathrm{{metal}}$ "
               f"({rng(E['metal'])}) and $E_\\mathrm{{both}}$ ({rng(E['both'])}), with the product "
               r"$E_\mathrm{oxide} \times E_\mathrm{metal}$ (black dotted) that a multiplicative combination would give"
               + (f"; the combination is {words}" if words else "") + ". Grey dashed: $E$ = 1.")
    return Panel("enhancement", fig, ax, caption, x_all, y_all, keyed=True)
