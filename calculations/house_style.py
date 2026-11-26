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


def _rgba(c):
    try:
        return tuple(np.round(mpl.colors.to_rgba(c), 3))
    except (ValueError, TypeError):
        return None


def _drawn_line(a):
    """Whether a Line2D is drawn with a visible line (not only markers)."""
    return a.get_linestyle() not in ("None", "none", "", " ") and a.get_linewidth() > 0


def _drawn_marker(a):
    return a.get_marker() not in (None, "None", "none", "", " ")


def _draws(handle, a):
    """Whether the artist a on the axes is (part of) what the key handle stands for: the same colour,
    and the handle's marker (with its fill) or its line style. A scatter matches by colour when the
    handle has a marker."""
    from matplotlib.collections import PathCollection
    if not isinstance(handle, Line2D) or not a.get_visible():
        return False
    if isinstance(a, PathCollection):
        colours = {_rgba(c) for c in list(a.get_facecolor()) + list(a.get_edgecolor())}
        rgb = _rgba(handle.get_color())
        return _drawn_marker(handle) and any(c is not None and c[:3] == rgb[:3] for c in colours)
    if not isinstance(a, Line2D) or _rgba(handle.get_color()) != _rgba(a.get_color()):
        return False
    marker = _drawn_marker(handle) and _drawn_marker(a) and a.get_marker() == handle.get_marker() \
        and _rgba(a.get_markerfacecolor()) == _rgba(handle.get_markerfacecolor())
    line = _drawn_line(handle) and _drawn_line(a) and a.get_linestyle() == handle.get_linestyle()
    return marker or line


def _axes_points(ax, a, step=0.004):
    """An artist's points in axes fractions; along a drawn line, every ``step`` of the axes, so a label
    can be kept off the line between the points as well as off the points."""
    if isinstance(a, Line2D):
        xy, trans = np.asarray(a.get_xydata(), float), a.get_transform()
    else:                                                   # a scatter collection
        xy, trans = np.asarray(a.get_offsets(), float), a.get_offset_transform()
    if xy.ndim != 2 or not len(xy):
        return np.empty((0, 2))
    with np.errstate(all="ignore"):
        pts = ax.transAxes.inverted().transform(trans.transform(xy))
    pts = pts[np.isfinite(pts).all(axis=1)]
    if isinstance(a, Line2D) and _drawn_line(a) and len(pts) > 1:
        dense = [pts[:1]]
        for p, q in zip(pts[:-1], pts[1:]):
            n = int(np.ceil(np.hypot(*(q - p)) / step))
            dense.append(p + (q - p) * np.linspace(0, 1, max(n, 1) + 1)[1:, None])
        pts = np.vstack(dense)
    return pts[(pts[:, 0] > -0.02) & (pts[:, 0] < 1.02) & (pts[:, 1] > -0.02) & (pts[:, 1] < 1.02)]


def _obstacles(ax, renderer):
    """What a label inside the plot must not cover: the points of every line, mark and scatter (in axes
    fractions), and the boxes of every text and filled patch. None when an image or mesh fills the
    plot (a heatmap): there is no free space to write in."""
    from matplotlib.collections import PathCollection, QuadMesh
    if ax.images or any(isinstance(c, QuadMesh) for c in ax.collections):
        return None
    pts = [_axes_points(ax, a) for a in ax.lines]
    pts += [_axes_points(ax, c) for c in ax.collections if isinstance(c, PathCollection)]
    boxes = [bbox_frac(ax, t, renderer) for t in ax.texts if t.get_visible() and t.get_text().strip()]
    for patch in ax.patches:
        b = bbox_frac(ax, patch, renderer)
        if patch.get_fill():
            boxes.append(b)
        else:                                               # an outline: its edges only
            u = np.linspace(0, 1, 60)
            pts.append(np.vstack([np.c_[b.x0 + u * b.width, np.full_like(u, b.y0)],
                                  np.c_[b.x0 + u * b.width, np.full_like(u, b.y1)],
                                  np.c_[np.full_like(u, b.x0), b.y0 + u * b.height],
                                  np.c_[np.full_like(u, b.x1), b.y0 + u * b.height]]))
    pts = [p for p in pts if len(p)]
    return (np.vstack(pts) if pts else np.empty((0, 2))), boxes


def _break_points(line):
    """Where a line may break for a label inside the plot: after ', ', '; ' or ': ' and before ' (',
    never inside $...$ math. Returns [(piece, separator before it)], the first separator ''."""
    pieces, cur, sep, math, i = [], "", "", False, 0
    while i < len(line):
        ch = line[i]
        if ch == "$":
            math = not math
        hit = next((t for t in (", ", "; ", ": ", " (") if not math and line.startswith(t, i)), None)
        if hit:
            pieces.append((cur, sep))
            cur, sep = ("(" if hit == " (" else ""), hit
            i += len(hit)
            continue
        cur += ch
        i += 1
    return pieces + [(cur, sep)]


def _inside_text(label, width=26):
    """A key label rewritten for the plot: each line longer than about ``width`` characters breaks at
    its commas, semicolons, colons and brackets ('Intact oxide, metal-diffusion limited' -> 'Intact
    oxide' over 'metal-diffusion limited'), so the label is narrower and finds room more easily. The
    separator at a break is dropped; pieces that stay on one line keep theirs."""
    visible = lambda t: len(re.sub(r"\$[^$]*\$", "xx", t))   # noqa: E731
    out = []
    for line in str(label).split("\n"):
        if visible(line) <= width:
            out.append(line)
            continue
        broken = []
        for part, sep in _break_points(line):
            joined = (broken[-1] + (" " if sep == " (" else sep) + part) if broken else None
            if joined is not None and visible(joined) <= width:
                broken[-1] = joined
            else:
                broken.append(part)
        out += broken
    return "\n".join(out)


def _glyph(ax, handle, x_right, y, gw):
    """The handle's mark (and line) drawn just left of a label, in axes fractions."""
    kw = dict(color=handle.get_color(), transform=ax.transAxes, clip_on=False, zorder=6)
    out = []
    if _drawn_line(handle):
        out += ax.plot([x_right - gw, x_right - 0.008], [y, y], ls=handle.get_linestyle(),
                       lw=handle.get_linewidth(), solid_capstyle="butt", **kw)
    if _drawn_marker(handle):
        out += ax.plot([x_right - (gw + 0.008) / 2 if _drawn_line(handle) else x_right - gw / 2], [y], ls="none",
                       marker=handle.get_marker(), ms=min(handle.get_markersize(), 10),
                       mfc=handle.get_markerfacecolor(), mec=handle.get_markeredgecolor(),
                       mew=handle.get_markeredgewidth(), **kw)
    return out


def _occupancy(pts, boxes, n):
    """A summed-area table of an n x n grid over the axes: a cell counts when a point or box is in it."""
    occ = np.zeros((n, n))
    if len(pts):
        ij = np.clip((pts * n).astype(int), 0, n - 1)
        occ[ij[:, 1], ij[:, 0]] = 1
    for b in boxes:
        i0, i1 = int(np.clip(np.floor(b.x0 * n), 0, n)), int(np.clip(np.ceil(b.x1 * n), 0, n))
        j0, j1 = int(np.clip(np.floor(b.y0 * n), 0, n)), int(np.clip(np.ceil(b.y1 * n), 0, n))
        occ[j0:j1, i0:i1] = 1
    return np.pad(occ.cumsum(0).cumsum(1), ((1, 0), (1, 0)))


def place_labels(ax, entries, renderer, margin=0.025, pad=0.014, n=200, step=0.01, keep=400, max_tries=40000):
    """Write each label in its colour, in bold, in free space next to its own data.

    ``entries`` is [(text, colour, home, glyph), ...]: home holds the entry's points in axes fractions
    (None for an entry nothing on the plot draws: it goes to the clearest corner), glyph is None or the
    key handle whose mark is drawn left of the label (for entries that share a colour, whose text colour
    is not the series colour, or that draw nothing). Every position on a grid over the plot that
    covers no line, mark or text is a free spot; its cost is its distance to the nearest point of the
    label's own data. The labels take the free spots with the least total cost that do not overlap each
    other. Returns (artists, ok); when ok is False nothing is left on the axes and the caller uses a key
    below instead."""
    found_obstacles = _obstacles(ax, renderer)
    if found_obstacles is None:
        return [], False
    pts, boxes = found_obstacles
    sat = _occupancy(pts, boxes, n)
    texts, spots = [], []
    for text, colour, home, glyph in entries:
        t = ax.text(0.5, 0.5, text, transform=ax.transAxes, color=colour, fontsize=ANNOT_PT, fontweight="bold",
                    ha="left", va="bottom", zorder=6, linespacing=1.15)
        texts.append(t)
        size = bbox_frac(ax, t, renderer)
        w, h = size.width, size.height
        gw = 0.0 if glyph is None else (0.075 if _drawn_line(glyph) else 0.035)
        lefts = np.arange(margin + gw, 1 - margin - w + 1e-9, step)
        bottoms = np.arange(margin, 1 - margin - h + 1e-9, step)
        if not len(lefts) or not len(bottoms):
            spots.append([])
            continue
        L, B = (a.ravel() for a in np.meshgrid(lefts, bottoms))
        x0, x1, y0, y1 = L - gw - pad, L + w + pad, B - pad, B + h + pad
        i0, i1 = np.clip(np.floor(x0 * n), 0, n).astype(int), np.clip(np.ceil(x1 * n), 0, n).astype(int)
        j0, j1 = np.clip(np.floor(y0 * n), 0, n).astype(int), np.clip(np.ceil(y1 * n), 0, n).astype(int)
        free = (sat[j1, i1] - sat[j0, i1] - sat[j1, i0] + sat[j0, i0]) == 0
        if not free.any():
            spots.append([])
            continue
        L, B, x0, x1, y0, y1 = (a[free] for a in (L, B, x0, x1, y0, y1))
        hp = (np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]]) if home is None
              else home[:: max(1, len(home) // 400)] if len(home) else np.array([[0.5, 0.5]]))
        dx = np.maximum(np.maximum(x0[:, None] - hp[None, :, 0], hp[None, :, 0] - x1[:, None]), 0)
        dy = np.maximum(np.maximum(y0[:, None] - hp[None, :, 1], hp[None, :, 1] - y1[:, None]), 0)
        cost = np.hypot(dx, dy).min(axis=1)
        mine = []
        for k in np.argsort(cost, kind="stable"):
            if any(abs(L[k] - m[0]) < 0.025 and abs(B[k] - m[1]) < 0.025 for m in mine):
                continue                                    # nearly the same spot as a better one
            mine.append((L[k], B[k], x0[k], x1[k], y0[k], y1[k], cost[k]))
            if len(mine) >= keep:
                break
        spots.append(np.array(mine) if mine else np.empty((0, 7)))

    if not all(len(sp) for sp in spots):
        for t in texts:
            t.remove()
        return [], False

    # The least total cost with no two labels overlapping: a depth-first search that places the label
    # with the fewest spots left first, drops every spot that would overlap it from the other labels'
    # lists (backing off at once when a list empties), and abandons a branch that cannot beat the best.
    best, tries = [np.inf, None], [0]

    def search(remaining, chosen, total):
        if not remaining:
            best[:] = [total, dict(chosen)]
            return
        floor = {j: spots[j][idx[0], 6] for j, idx in remaining.items()}
        if total + sum(floor.values()) >= best[0] or tries[0] > max_tries:
            return
        i = min(remaining, key=lambda j: len(remaining[j]))
        others = sum(v for j, v in floor.items() if j != i)
        for m in remaining[i]:
            tries[0] += 1
            c = total + spots[i][m, 6]
            if c + others >= best[0] or tries[0] > max_tries:
                break                                       # spots are sorted by cost: none better follows
            bx0, bx1, by0, by1 = spots[i][m, 2:6]
            left, ok = {}, True
            for j, idx in remaining.items():
                if j == i:
                    continue
                o = spots[j][idx]
                keep_j = idx[(o[:, 3] < bx0) | (o[:, 2] > bx1) | (o[:, 5] < by0) | (o[:, 4] > by1)]
                if not len(keep_j):
                    ok = False
                    break
                left[j] = keep_j
            if ok:
                chosen[i] = m
                search(left, chosen, c)
                del chosen[i]

    search({i: np.arange(len(sp)) for i, sp in enumerate(spots)}, {}, 0.0)
    if best[1] is None:
        for t in texts:
            t.remove()
        return [], False
    artists = list(texts)
    for i, (text, colour, home, glyph) in enumerate(entries):
        left, bottom = spots[i][best[1][i], :2]
        texts[i].set_position((left, bottom))
        if glyph is not None:
            box = bbox_frac(ax, texts[i], renderer)
            gw = 0.075 if _drawn_line(glyph) else 0.035
            line_h = box.height / (text.count("\n") + 1)
            artists += _glyph(ax, glyph, box.x0 - 0.006, box.y1 - line_h / 2, gw - 0.006)
    return artists, True


def labels_or_key(fig, ax, handles, labels, ncol=4, colours=None, fallback=True):
    """Each key entry's label inside the plot next to its own data (research-plots: labels next to the
    data, not a key), or, when one does not find free space, the whole key below the plot instead.

    The data of an entry are the artists on ``ax`` its handle draws (same colour, marker and line
    style); an entry nothing draws goes to the clearest corner. Call after the axes limits, ticks and any
    notes are final. Returns ("inside", artists) or ("key", legend); with fallback=False a figure whose
    labels do not fit gets no key and ("none", None) comes back, so the caller can make room and retry."""
    if not labels:
        return "inside", []
    if ax is None:
        ax = fig.axes[0]
    texts = ["black" if c in (GREY, None, "white") else c
             for c in (colours or [getattr(h, "get_color", lambda: "black")() for h in handles])]
    texts = [c if isinstance(c, str) else mpl.colors.to_hex(c) for c in texts]
    fig.canvas.draw()                                       # the constrained layout has placed the axes
    renderer = fig.canvas.get_renderer()
    entries = []
    for h, lab, ink in zip(handles, labels, texts):
        home = [_axes_points(ax, a) for a in list(ax.lines) + list(ax.collections) if _draws(h, a)]
        home = [p for p in home if len(p)]
        series = mpl.colors.to_hex(h.get_color())
        shared = sum(t.lower() == ink.lower() for t in texts) > 1 or ink.lower() != series.lower()
        # an entry nothing on the plot draws (an explanation such as 'Open: over half of the resistance')
        # carries its mark and goes to a clear corner: its home is the four corners of the plot
        home = np.vstack(home) if home else None
        entries.append((_inside_text(lab), ink, home, h if (shared or home is None) else None))
    artists, ok = place_labels(ax, entries, renderer)
    if ok:
        return "inside", artists
    if not fallback:
        return "none", None
    return "key", key_below(fig, handles, labels, ncol=ncol, colours=colours)


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


def grow_y(ax, frac):
    """Extend the y-axis by ``frac`` of its height above and below (in decades on a log axis), with
    the ticks redone for the new range; room for labels inside the plot."""
    lo, hi = ax.get_ylim()
    if ax.get_yscale() == "log":
        u0, u1 = np.log10(lo), np.log10(hi)
        lo, hi = 10 ** (u0 - frac * (u1 - u0)), 10 ** (u1 + frac * (u1 - u0))
        ax.set_ylim(lo, hi)
        log_ticks(ax.yaxis, lo, hi)
    else:
        span = hi - lo
        lo, hi = lo - frac * span, hi + frac * span
        ax.set_ylim(lo, hi)
        if isinstance(ax.yaxis.get_major_locator(), FixedLocator):
            ax.set_yticks([t for t in MaxNLocator(6).tick_values(lo, hi) if lo <= t <= hi])
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
    return labels_or_key(fig, fig.axes[0], handles, labels, ncol=ncol, colours=colours)


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
