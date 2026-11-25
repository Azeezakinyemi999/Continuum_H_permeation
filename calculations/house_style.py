"""The house plotting style of every figure the package and its notebooks draw.

The research-plots rules (180 x 128 mm at final size, TeX Gyre Heros or its URW original Nimbus
Sans, 17 pt ticks and 20 pt axis titles, a 2.2 pt black frame with ticks inside on all four sides,
curves of 3 pt or more, black text, no grid, no titles), colours that keep one meaning in every
figure of the project, and every figure written as PDF, 600 dpi PNG and a LaTeX caption whose
numbers come from the results (``calculations.figures.save_panels``). Each panel of a multi-panel
figure is its own figure.

Colours (confirmed by the user on 2026-11-24):
  limiting step  oxide diffusion vermillion, metal diffusion blue, metal surface purple, oxide surface
                 green, mixed and undefined grey (marks only; their labels are black); marker shape and
                 open/filled come from ``wall_labels.describe``
  model level    L1 black, L2a orange, L2b blue, L3 green, L4 (all variants) sky blue told apart by line
                 style and marker, L5 vermillion; a +L6 level keeps its level's hue, dashed with squares
  The two maps reuse hues (vermillion is oxide diffusion and L5), so one figure never uses both.
Fits, theory lines, limits and thresholds are grey dashed lines that step back; conditions
(temperatures, defect types) are black, told apart by marker and line style.
"""
import math
import os
import re

import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, MaxNLocator, NullFormatter

import matplotlib as mpl

from calculations.research_style import (ANNOT_PT, COLOURS, CURVE_PT, MM, OUTLINE_PT,  # noqa: F401
                                         check_figure, use_research_style)
from calculations.wall_labels import describe, fitted, parse

# --------------------------------------------------------------------------------------------------
# Names
# --------------------------------------------------------------------------------------------------

#: The comparison studies, in the order of study-comparison.ipynb.
STUDIES = ["Guo_etal_2025", "Fuerst_et_al_2024", "incoloy802_cr2o3"]

#: Metal key in a study_config.py -> the alloy's name on a figure.
ALLOY_NAMES = {
    "316L_Guo_2025": "316L",
    "Hastelloy_N_fuerst_2024": "Hastelloy N",
    "metal_X40_NiCrAlTi_31_19_Incoloy802_Schmidt1985": "Incoloy 802",
}

#: Model level -> what the wall is, for captions.
LEVEL_WALLS = {
    "L1": "a perfect metal (lattice diffusion)",
    "L2a": "a perfect oxide alone (molecular H$_2$)",
    "L2b": "a perfect oxide on a perfect metal",
    "L3": "a defective oxide on a perfect metal",
    "L4": "a defective metal (traps and grain boundaries)",
    "L5": "a defective oxide on a defective metal",
}


def level_wall(level):
    """Caption words for a level key such as 'L3' or 'L3L6'."""
    base, surface = (level[:-2], True) if level.endswith("L6") else (level, False)
    return LEVEL_WALLS[base] + (", with surface kinetics" if surface else "")


# --------------------------------------------------------------------------------------------------
# Colours: one meaning each, in every paper figure
# --------------------------------------------------------------------------------------------------

#: Limiting step -> colour. Used where colour means the limiting step (regime figures).
STEP_COLOURS = {
    "oxide diffusion": COLOURS["vermilion"],
    "metal diffusion": COLOURS["blue"],
    "metal surface": COLOURS["purple"],
    "oxide surface": COLOURS["green"],
    "mixed": COLOURS["grey_points"],          # not fitted: the marks step back
    "undefined": COLOURS["grey_points"],
}

#: Model level -> line and marker style. Used where colour means the level (level figures).
#: One hue per level family; the three L4 variants share sky blue and differ by line and marker;
#: a +L6 (surface kinetics) variant keeps its level's hue and is dashed with square markers.
LEVEL_STYLE = {
    "L1": dict(color=COLOURS["black"], ls="-", marker="o"),
    "L2a": dict(color=COLOURS["orange"], ls="-", marker="o"),
    "L2b": dict(color=COLOURS["blue"], ls="-", marker="o"),
    "L3": dict(color=COLOURS["green"], ls="-", marker="o"),
    "L4": dict(color=COLOURS["sky"], ls="-", marker="D"),
    "L4_gb": dict(color=COLOURS["sky"], ls="-", marker="^"),
    "L4_trap": dict(color=COLOURS["sky"], ls="--", marker="v"),
    "L4_both": dict(color=COLOURS["sky"], ls=":", marker="D"),
    "L5": dict(color=COLOURS["vermilion"], ls="-", marker="o"),
}


#: Sensitivity-analysis regime (the cluster a sample falls in) -> (colour, marker, name). The colours
#: follow STEP_COLOURS: oxide-limited is vermilion, metal-limited blue, surface-limited green;
#: defect-dominated (L5 only) takes orange, which no limiting step uses.
SA_REGIME_STYLE = {
    "oxide": (COLOURS["vermilion"], "s", "Oxide-limited"),
    "metal": (COLOURS["blue"], "o", "Metal-limited"),
    "surface": (COLOURS["green"], "P", "Surface-limited"),
    "defect": (COLOURS["orange"], "D", "Defect-dominated"),
}

def level_style(level):
    """Line and marker style for a level key; '+L6' variants are dashed with squares."""
    if level.endswith("L6") and level[:-2] in LEVEL_STYLE:
        return {**LEVEL_STYLE[level[:-2]], "ls": "--", "marker": "s"}
    return dict(LEVEL_STYLE[level])


# --------------------------------------------------------------------------------------------------
# Regimes (limiting steps)
# --------------------------------------------------------------------------------------------------

_PATH_WORDS = {"intact": "Intact oxide", "pinholes": "Pinholes", "cracks": "Cracks", "oxide GBs": "Oxide GBs"}
_STEP_WORDS = {
    "oxide diffusion": "oxide-diffusion limited", "metal diffusion": "metal-diffusion limited",
    "oxide surface": "oxide-surface limited", "metal surface": "metal-surface limited",
    "mixed": "mixed", "undefined": "solve failed",
}
MARKER_WORDS = {"o": "circles", "s": "squares", "X": "crosses", "x": "crosses", "P": "plus signs",
                "D": "diamonds", "^": "triangles", "v": "triangles"}


def regime_mark(label):
    """(colour, marker, hollow) for a regime label from wall_labels."""
    _, colour_old, marker, fill = describe(label)
    step = parse(label)[1]
    return STEP_COLOURS.get(step, COLOURS["grey_points"]), marker, fill == "hollow"


def regime_short(label):
    """A short label for the figure: 'Pinholes\\nmetal-surface limited', 'Oxide-diffusion limited'."""
    path, step, detail = parse(label)
    if step == "mixed":            # no single limiting step, whichever path dominates: one 'Mixed' group
        return "Mixed"
    words = _STEP_WORDS.get(step, step)
    if detail == "traps":
        words += ", trap-controlled"
    if path:
        return f"{_PATH_WORDS.get(path, path)}\n{words}"
    return words[0].upper() + words[1:]


def regime_long(label):
    """The full regime name, for captions (wall_labels.describe's words)."""
    return describe(label)[0]


def regime_fit_name(label):
    """The regime's name as study-comparison.ipynb writes it in its fit tables ('H$_2$' -> 'H₂')."""
    return regime_long(label).replace("$_2$", "₂")


is_fitted = fitted

# --------------------------------------------------------------------------------------------------
# Axes
# --------------------------------------------------------------------------------------------------


def log_ticks(axis, lo, hi):
    """Major ticks on a log axis. Plain numbers (0.1, 1, 10) when every label lies between 0.001 and
    10,000 and the axis spans six decades or fewer; otherwise powers of ten every 1, 2 or 3 decades,
    so the labels never crowd."""
    k0, k1 = int(np.ceil(np.log10(lo) - 1e-9)), int(np.floor(np.log10(hi) + 1e-9))
    span = k1 - k0
    if span < 2:                                 # under two decades: 1, 2 and 5 times powers of ten
        for subs in ((1, 2, 5), (1, 2, 3, 5, 7), (1, 2, 3, 4, 5, 6, 7, 8, 9)):
            ticks = [m * 10.0 ** k for k in range(k0 - 1, k1 + 1) for m in subs if lo <= m * 10.0 ** k <= hi]
            if len(ticks) >= 3:                  # a short axis: at least three labelled ticks
                break
        axis.set_major_locator(FixedLocator(ticks))
        axis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=200))
        axis.set_minor_formatter(NullFormatter())

        plain = all(-3 <= int(np.floor(np.log10(t) + 1e-9)) <= 3 for t in ticks)   # one format per axis

        def mantissa(v, _):
            k = int(np.floor(np.log10(v) + 1e-9))
            m = round(v / 10 ** k)
            if plain:
                return f"{v:,.{max(0, -k)}f}"
            return rf"$10^{{{k}}}$" if m == 1 else rf"${m}\times10^{{{k}}}$"
        axis.set_major_formatter(FuncFormatter(mantissa))
        return
    step = 1 if span <= 6 else (2 if span <= 12 else 3)
    ks = [k for k in range(k0, k1 + 1) if k % step == 0] or [k0]
    axis.set_major_locator(FixedLocator([10.0 ** k for k in ks]))
    axis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=200) if step > 1
                           else LogLocator(base=10, subs=np.arange(2, 10), numticks=200))
    axis.set_minor_formatter(NullFormatter())
    plain = span <= 6 and min(ks) >= -3 and max(ks) <= 4
    if plain:
        axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.{max(0, -int(np.floor(np.log10(v))))}f}"))
    else:
        axis.set_major_formatter(FuncFormatter(lambda v, _: rf"$10^{{{int(round(np.log10(v)))}}}$"))


def celsius_top_axis(ax, label="Temperature (°C)", min_gap=0.11):
    """A top axis in °C over a 1000/T axis, with round ticks spaced at least min_gap (axes fraction) apart."""
    ax.tick_params(axis="x", top=False)
    top = ax.secondary_xaxis("top", functions=(lambda v: 1000 / np.maximum(v, 1e-9) - 273.15,
                                               lambda c: 1000 / (c + 273.15)))
    top.set_xlabel(label)
    return top


def set_celsius_ticks(ax, top, min_gap=0.11):
    """Choose the °C ticks after the x-limits are final: every 100 or 200 °C, none closer than min_gap."""
    x0, x1 = ax.get_xlim()
    lo_c, hi_c = sorted(1000 / np.array([x0, x1]) - 273.15)
    step = 200 if hi_c - lo_c > 500 else 100
    chosen, last = [], None
    for c in sorted((c for c in range(0, 2001, step) if lo_c <= c <= hi_c), reverse=True):   # left to right
        f = (1000 / (c + 273.15) - x0) / (x1 - x0)
        if last is None or f - last >= min_gap:
            chosen.append(c)
            last = f
    top.set_xticks(chosen)


def log_limits(values, pad_decades=0.3, frac=0.05):
    """(lo, hi) around positive values, padded in log space by pad_decades or by frac of the span,
    whichever is larger, so end markers never touch the frame on an axis of many decades."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v) & (v > 0)]
    lo, hi = np.log10(v.min()), np.log10(v.max())
    pad = max(pad_decades, frac * (hi - lo))
    return 10 ** (lo - pad), 10 ** (hi + pad)


def assert_inside(ax, x, y):
    """Fail if any finite data point falls outside the axis limits."""
    (x0, x1), (y0, y1) = sorted(ax.get_xlim()), sorted(ax.get_ylim())
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    inside = (x[ok] >= x0) & (x[ok] <= x1) & (y[ok] >= y0) & (y[ok] <= y1)
    assert inside.all(), f"{(~inside).sum()} data points fall outside the axes"


def _to_axes(ax, x, y):
    """Data coordinates -> axes fractions."""
    pts = ax.transData.transform(np.column_stack([x, y]))
    return ax.transAxes.inverted().transform(pts)


# --------------------------------------------------------------------------------------------------
# Labels next to the data
# --------------------------------------------------------------------------------------------------


def room_for_end_labels(ax, names, renderer, offset_pt=10, log=True, side="right"):
    """Widen the x-axis on one side so labels written beyond each curve's end on that side fit."""
    fig = ax.figure
    probe = ax.text(0, 0, max(names, key=len), fontsize=ANNOT_PT, fontweight="bold", transform=ax.transAxes)
    width_px = probe.get_window_extent(renderer).width + offset_pt * fig.dpi / 72 + 6
    probe.remove()
    frac = min(0.45, width_px / ax.get_window_extent(renderer).width)
    x0, x1 = ax.get_xlim()
    grow = frac / (1 - frac)
    if log:
        u0, u1 = np.log10(x0), np.log10(x1)
        ax.set_xlim(*((x0, 10 ** (u1 + grow * (u1 - u0))) if side == "right" else (10 ** (u0 - grow * (u1 - u0)), x1)))
    else:
        ax.set_xlim(*((x0, x1 + grow * (x1 - x0)) if side == "right" else (x0 - grow * (x1 - x0), x1)))


def end_labels(ax, ends, gap=0.075, offset_pt=10, floor=0.045, ceiling=0.95, side="right"):
    """Name each curve beyond its end on ``side``. ``ends`` is [(x, y, name), ...] in data units.

    The top label stays level with its own curve; each label below is pushed down where two curves
    nearly coincide, so no label sits level with a different curve. When the stack would run out of
    the bottom of the axes, the y-axis is extended downwards (never the labels moved up into the
    data). Returns the text artists.
    """
    def stack():
        frac = sorted(((_to_axes(ax, [x], [y])[0], name) for x, y, name in ends), key=lambda t: -t[0][1])
        placed = []
        for (_, yf), _ in frac:
            placed.append(yf if not placed else min(yf, placed[-1] - gap))
        return frac, placed

    frac, placed = stack()
    for _ in range(8):                                 # room at the top: the top label stays level with its curve
        if not placed or max(placed) <= ceiling:
            break
        y0, y1 = ax.get_ylim()
        if ax.get_yscale() == "log":
            u0, u1 = np.log10(y0), np.log10(y1)
            ax.set_ylim(y0, 10 ** (u1 + (max(placed) - ceiling + 0.01) * (u1 - u0) * 1.2))
        else:
            ax.set_ylim(y0, y1 + (max(placed) - ceiling + 0.01) * (y1 - y0) * 1.2)
        frac, placed = stack()
    for _ in range(8):
        if not placed or min(placed) >= floor:
            break
        y0, y1 = ax.get_ylim()
        if ax.get_yscale() == "log":
            u0, u1 = np.log10(y0), np.log10(y1)
            ax.set_ylim(10 ** (u0 - (floor - min(placed) + 0.01) * (u1 - u0) * 1.2), y1)
        else:
            ax.set_ylim(y0 - (floor - min(placed) + 0.01) * (y1 - y0) * 1.2, y1)
        frac, placed = stack()
    texts = []
    for ((xf, _), name), y_label in zip(frac, placed):
        texts.append(ax.annotate(name, (xf, y_label), xycoords="axes fraction",
                                 xytext=(offset_pt if side == "right" else -offset_pt, 0),
                                 textcoords="offset points", va="center", ha="left" if side == "right" else "right",
                                 fontsize=ANNOT_PT, fontweight="bold", color=COLOURS["black"], zorder=6))
    return texts


def bbox_frac(ax, artist, renderer):
    """An artist's extent in axes fractions."""
    return artist.get_window_extent(renderer).transformed(ax.transAxes.inverted())


_bbox_frac = bbox_frac


def place_labels(ax, groups, renderer, avoid=(), margin=0.03, pad=0.018, glyph_w=0.035):
    """Write each group's label in its colour, in free space next to its points.

    ``groups`` is [(text, colour, x, y, glyph), ...]: x and y are arrays in data units, and glyph is
    None or a marker style (``marker``, ``mfc``, ``mec``, ``mew``) drawn just left of the label, for
    groups that share a colour and differ only by marker. ``avoid`` holds artists (such as the end
    labels) that no label may cover. Each label is tried on rings around its group's centre, below
    and to the right first, and goes to the nearest spot that covers no data point and no other
    label.

    Returns ``(artists, ok)``: every text and glyph drawn, and whether all labels found a free spot.
    When ``ok`` is False the caller should remove the artists and use a key instead.
    """
    all_pts = np.vstack([_to_axes(ax, g[2], g[3]) for g in groups]) if groups else np.empty((0, 2))
    all_pts = all_pts[np.isfinite(all_pts).all(axis=1)]
    taken = [_bbox_frac(ax, a, renderer) for a in avoid]
    angles = np.deg2rad([-90, -60, -120, -30, -150, 0, 180, 30, 150, -75, -105, 60, 120, 90])
    artists, ok = [], True
    for text, colour, x, y, glyph in groups:
        pts = _to_axes(ax, x, y)
        pts = pts[np.isfinite(pts).all(axis=1)]
        centre = np.median(pts, axis=0) if len(pts) else np.array([0.5, 0.5])
        t = ax.text(0, 0, text, transform=ax.transAxes, color=colour, fontsize=ANNOT_PT,
                    fontweight="bold", zorder=6, linespacing=1.15)
        artists.append(t)
        gw = glyph_w if glyph else 0.0
        t.set_position((0.5, 0.5)); t.set_ha("left"); t.set_va("bottom")
        size = _bbox_frac(ax, t, renderer)
        w, h = size.width, size.height
        found = None
        for r in np.arange(0.05, 1.2, 0.025):
            for a in angles:
                dx, dy = np.cos(a), np.sin(a)
                cx, cy = centre[0] + r * dx, centre[1] + r * dy
                left = cx + gw if dx > 0.3 else (cx - w if dx < -0.3 else cx - (w + gw) / 2 + gw)
                bottom = cy if dy > 0.3 else (cy - h if dy < -0.3 else cy - h / 2)
                x0, x1, y0, y1 = left - gw - pad, left + w + pad, bottom - pad, bottom + h + pad
                if x0 + pad < margin or x1 - pad > 1 - margin or y0 + pad < margin or y1 - pad > 1 - margin:
                    continue
                if len(all_pts) and ((all_pts[:, 0] > x0) & (all_pts[:, 0] < x1)
                                     & (all_pts[:, 1] > y0) & (all_pts[:, 1] < y1)).any():
                    continue
                if any(not (x1 < o.x0 or x0 > o.x1 or y1 < o.y0 or y0 > o.y1) for o in taken):
                    continue
                t.set_position((left, bottom))
                t.set_multialignment("right" if (dx < -0.3 and not glyph) else "left")
                found = _bbox_frac(ax, t, renderer)
                break
            if found is not None:
                break
        if found is None:
            ok = False
            continue
        taken.append(found)
        if glyph:
            line_h = found.height / (text.count("\n") + 1)
            artists += ax.plot([found.x0 - gw / 2], [found.y1 - line_h / 2], ls="none", ms=9,
                               transform=ax.transAxes, clip_on=False, zorder=6, **glyph)
    return artists, ok


# --------------------------------------------------------------------------------------------------
# Captions
# --------------------------------------------------------------------------------------------------


def tex_num(v, sig=2):
    """A number for a LaTeX caption: 0.66, 18, 490, or $10^{-7}$ and $2.5\\times10^{4}$ outside 0.01-1000;
    true minus signs."""
    v = float(v)
    if not np.isfinite(v):
        return "n/a"
    if v == 0:
        return "0"
    if 0.01 <= abs(v) < 1000:
        return f"{float(f'{v:.{sig}g}'):g}".replace("-", "$-$")
    k = int(np.floor(np.log10(abs(v))))
    a = float(f"{v / 10 ** k:.{sig}g}")
    if abs(a) >= 10:                                   # 9.996e-14 rounds to 10: write 1e-13
        a, k = a / 10, k + 1
    body = f"10^{{{k}}}" if abs(a) == 1 else f"{abs(a):g}\\times10^{{{k}}}"
    return f"${'-' if a < 0 else ''}{body}$"


def value_range(values, fmt="{:.2f}"):
    """'1.00', '0.51–0.52', or '−0.99 to −0.62' (true minus signs; 'to' when a value is negative)."""
    vals = sorted({float(fmt.format(v + 0.0)) for v in values})
    vals = [0.0 if v == 0 else v for v in vals]                    # no '−0.00'
    text = [fmt.format(v).replace("-", "\u2212") for v in vals]
    if len(text) == 1:
        return text[0]
    return f"{text[0]} to {text[-1]}" if vals[0] < 0 else f"{text[0]}–{text[-1]}"


def tex_join(items):
    """'a', 'a and b', 'a, b and c'."""
    items = list(items)
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


# --------------------------------------------------------------------------------------------------
# Figures, keys, limits and captions
# --------------------------------------------------------------------------------------------------

FIGSIZE = (180 * MM, 128 * MM)
FALLBACK_FACE = "Nimbus Sans"          # URW's Helvetica design; TeX Gyre Heros is built from it
_FACE = None
GREY = COLOURS["grey_points"]
#: Fits, theory lines, limits and thresholds: grey dashed, stepping back behind the data.
REF_LINE = dict(color=GREY, lw=2.4, ls=(0, (5, 3)))


def use_style():
    """The research style, once per process. Where neither TeX Gyre Heros, Helvetica nor Arial is
    installed (a cluster without TeX), the same design from URW, Nimbus Sans, instead of DejaVu Sans."""
    global _FACE
    if _FACE is not None:
        use_research_style()
        if _FACE == FALLBACK_FACE:
            _fallback_rc(_FACE)
        return _FACE
    face = use_research_style()
    if face == "DejaVu Sans":
        for root in ("/usr/share/fonts", "/usr/local/share/fonts"):
            for dirpath, _dirs, files in os.walk(root):
                for f in files:
                    if f.startswith("NimbusSans-") and f.endswith(".otf"):
                        font_manager.fontManager.addfont(os.path.join(dirpath, f))
        if FALLBACK_FACE in {f.name for f in font_manager.fontManager.ttflist}:
            face = FALLBACK_FACE
            _fallback_rc(face)
    _FACE = face
    return face


def _fallback_rc(face):
    mpl.rcParams.update({"font.sans-serif": [face, "DejaVu Sans"], "mathtext.rm": f"{face}:regular",
                         "mathtext.it": f"{face}:italic", "mathtext.bf": f"{face}:bold"})


def new_figure():
    """A 180 x 128 mm figure and its axes, in the house style."""
    import matplotlib.pyplot as plt
    use_style()
    return plt.subplots(figsize=FIGSIZE)


def line_kw(level, lw=CURVE_PT, ms=9, markers=True):
    """Keyword arguments for ax.plot of one model level."""
    s = level_style(level)
    kw = dict(color=s["color"], ls=s["ls"], lw=lw, zorder=3 if s["ls"] == "-" else 4)
    if markers:
        kw.update(marker=s["marker"], ms=ms, mfc=s["color"], mec=s["color"], mew=1.2)
    return kw


def set_log_y(ax, values, **kw):
    """A log y axis around the values; values all within a factor of 1.5 get a tight linear axis with
    plain numbers instead, so small differences stay visible."""
    v = np.asarray(values, dtype=float).ravel()
    v = v[np.isfinite(v) & (v > 0)]
    if len(v) and v.max() / v.min() < 1.5:
        ax.set_yscale("linear")
        set_linear_y(ax, v, frac=0.12)
        lo, hi = ax.get_ylim()
        ax.set_yticks([t for t in MaxNLocator(6).tick_values(lo, hi) if lo <= t <= hi])
        ax.set_ylim(lo, hi)
        if v.max() < 0.01 or v.max() >= 1e4:               # tiny or huge values: ×10ⁿ once, at the axis top
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True, useOffset=False)
        else:
            ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        return
    lo, hi = log_limits(values, **kw)
    ax.set_yscale("log")
    ax.set_ylim(lo, hi)
    log_ticks(ax.yaxis, lo, hi)


def set_log_x(ax, values, **kw):
    lo, hi = log_limits(values, **kw)
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    log_ticks(ax.xaxis, lo, hi)


def set_linear_y(ax, values, frac=0.07, include=(), floor=None, ceiling=None):
    """y limits around the values (and `include`), padded; floor/ceiling clip the padding only."""
    v = np.concatenate([np.asarray(values, float).ravel(), np.asarray(include, float).ravel()])
    v = v[np.isfinite(v)]
    if not len(v):
        return
    lo, hi = v.min(), v.max()
    pad = frac * ((hi - lo) or abs(hi) or 1.0)
    lo, hi = lo - pad, hi + pad
    if floor is not None and v.min() >= floor:
        lo = max(lo, floor)
    if ceiling is not None and v.max() <= ceiling:
        hi = min(hi, ceiling)
    ax.set_ylim(lo, hi)


def reference_line(ax, y, vertical=False, **kw):
    """A grey dashed reference line (a fit, limit or threshold), named in the key by the caller."""
    style = {**REF_LINE, **kw}
    return ax.axvline(y, zorder=2, **style) if vertical else ax.axhline(y, zorder=2, **style)


def text_label(ax, x, y, text, color="black", ha="left", va="center", coords="data", backing=True, **kw):
    """A bold direct label on a white backing (grey becomes black: text is never grey)."""
    return ax.annotate(text, (x, y), xycoords="axes fraction" if coords == "axes" else "data",
                       ha=ha, va=va, fontsize=ANNOT_PT, fontweight="bold", color="black" if color == GREY else color,
                       zorder=7, linespacing=1.15,
                       bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none") if backing else None, **kw)


def points_outside(ax, x, y):
    """The number of finite data points outside the axis limits (0 when all are inside)."""
    if x is None or y is None:
        return 0
    (x0, x1), (y0, y1) = sorted(ax.get_xlim()), sorted(ax.get_ylim())
    x, y = np.asarray(x, float).ravel(), np.asarray(y, float).ravel()
    if x.size != y.size:
        return 0
    ok = np.isfinite(x) & np.isfinite(y)
    if ax.get_yscale() == "log":
        ok &= y > 0
    if ax.get_xscale() == "log":
        ok &= x > 0
    return int((~((x[ok] >= x0) & (x[ok] <= x1) & (y[ok] >= y0) & (y[ok] <= y1))).sum())


def key_below(fig, handles, labels, ncol=4, colours=None):
    """A frameless key below the plot; the figure grows by its rows so the plot keeps its size.
    colours: one text colour per entry (the series colour; grey becomes black, as text is never grey).
    An entry's second line carries the series' fitted values. ncol is the most columns: the key is
    measured as drawn and loses a column until it fits the figure's width."""
    if not labels:
        return None
    ncol = max(1, min(ncol, len(labels)))
    lines = max(1, max(str(lab).count("\n") + 1 for lab in labels))
    height = fig.get_figheight()
    while True:
        rows = math.ceil(len(labels) / ncol)
        fig.set_figheight(height + (8.5 * rows * lines + 2.5 * rows * (lines > 1) + 6) * MM)
        leg = fig.legend(handles, labels, loc="outside lower center", ncol=ncol, frameon=False,
                         handlelength=2.2, handletextpad=0.45, columnspacing=1.1, borderaxespad=0.2,
                         labelspacing=0.55 if lines > 1 else 0.5)
        if ncol == 1:
            break
        try:
            wide = leg.get_window_extent(fig.canvas.get_renderer()).width > 0.98 * fig.bbox.width
        except Exception:                                  # a backend without a renderer: keep the estimate
            break
        if not wide:
            break
        leg.remove()
        ncol -= 1
    if colours:
        for text, c in zip(leg.get_texts(), colours):
            text.set_color("black" if c in (GREY, None, "white") else c)
    return leg


def level_key(fig, levels, names=None, values=None, extra=(), markers=True, ncol=4):
    """The key of the levels drawn, plus extra (handle, label) rows. names: {level: label};
    values: {level: text} printed as the entry's second line, in the level's colour."""
    handles = [Line2D([], [], **line_kw(lv, markers=markers)) for lv in levels] + [h for h, _ in extra]
    labels = [(names or {}).get(lv, lv) + (f"\n{values[lv]}" if values and values.get(lv) else "")
              for lv in levels] + [lab for _, lab in extra]
    colours = [level_style(lv)["color"] for lv in levels] + ["black"] * len(extra)
    return key_below(fig, handles, labels, ncol=ncol, colours=colours)


def sci(v, sig=3):
    """A number for a key in mathtext: '0.426', '$5.98\\times10^{-7}$', with a true minus sign."""
    v = float(v)
    if v == 0 or not np.isfinite(v):
        return f"{v:g}"
    if 0.01 <= abs(v) < 1000:
        return f"{float(f'{v:.{sig}g}'):g}".replace("-", "−")
    k = int(np.floor(np.log10(abs(v))))
    a = float(f"{v / 10 ** k:.{sig}g}")
    if abs(a) >= 10:                                   # 9.996e-14 rounds to 10: write 1e-13
        a, k = a / 10, k + 1
    sign = "−" if a < 0 else ""
    return f"{sign}$10^{{{k}}}$" if abs(a) == 1 else f"{sign}${abs(a):g}\\times10^{{{k}}}$"


def kj(v, digits=None, unit=True):
    """'52 kJ mol$^{-1}$' for keys, from kJ/mol, with a true minus sign; one decimal below 10."""
    if v is None or not np.isfinite(v):
        return "n/a"
    digits = (1 if abs(v) < 10 else 0) if digits is None else digits
    return f"{v:.{digits}f}".replace("-", "−") + (r" kJ mol$^{-1}$" if unit else "")


def tex_kj(v, digits=0):
    return "n/a" if v is None or not np.isfinite(v) else f"{v:.{digits}f}".replace("-", "$-$") + r"~kJ~mol$^{-1}$"


def caption_tex(name, module, text):
    return rf"""% Written by {module}. Do not edit the numbers by hand.
\begin{{figure}}
  \centering
  \includegraphics[width=180mm]{{{name}.pdf}}
  \caption{{{text}}}
  \label{{fig:{name}}}
\end{{figure}}
"""


def finish(fig, ax, outstem, module, caption, x=None, y=None, keyed=False):
    """Frame on top, the checks, then <outstem>.pdf, <outstem>.png (600 dpi) and <outstem>_caption.tex.
    Problems are returned, never raised: a figure must not stop a notebook. Returns (files, problems)."""
    import matplotlib.pyplot as plt
    for a in fig.axes:
        for spine in a.spines.values():
            spine.set_zorder(10)
    name = os.path.basename(outstem)
    problems = [p for p in check_figure(fig) if not (keyed and "legend" in p)]
    if _FACE == FALLBACK_FACE:          # the check knows only TeX Gyre Heros, Helvetica and Arial
        problems = [p for p in problems if not (p.startswith("text not in") and "DejaVu" not in p)]
    n_out = points_outside(ax, x, y) if ax is not None else 0
    if n_out:
        problems.append(f"{n_out} data points fall outside the axes")
    os.makedirs(os.path.dirname(outstem) or ".", exist_ok=True)
    files = []
    for ext in ("pdf", "png"):                 # not pathlib's with_suffix: '1e-7Pa' is not an extension
        fig.savefig(f"{outstem}.{ext}", dpi=600, facecolor="white")
        files.append(f"{outstem}.{ext}")
    with open(outstem + "_caption.tex", "w") as fh:
        fh.write(caption_tex(name, module, caption))
    files.append(outstem + "_caption.tex")
    plt.close(fig)
    return files, problems


_FORMULA = re.compile(r"(?<=[A-Z])(\d+)|(?<=[A-Z][a-z])(\d+)")


def sub(s):
    """'H2O and Cr2O3' -> 'H$_{2}$O and Cr$_{2}$O$_{3}$' for matplotlib text; text with $ is left alone."""
    s = str(s)
    return s if "$" in s else _FORMULA.sub(lambda m: f"$_{{{m.group(0)}}}$", s)
