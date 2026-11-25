"""House-style figures of the Proposal notebook: the validation of each model level, one figure per panel.

The notebook computes every sweep (the physics stays there) and passes the arrays to one function per
level; each returns a list of ``calculations.figures.Panel`` for ``save_panels``. Each panel is a
180 x 128 mm figure in the house style (``calculations.house_style``): a level in its level hue
(``LEVEL_STYLE``; the defective-metal modes share L4's sky blue and differ by line style and marker),
conditions (temperatures, defect fractions) told apart by marker and line style, fits, theory lines,
limits and thresholds grey dashed, a key below the plot, no titles. The check the notebook makes on a
panel (fitted slope or energy against the expected value, R², the error of a limit check) is written on
the panel in a clear corner, and stated as failed when it fails. Every number in a caption comes from
the arrays passed in.

    from calculations import level_figures as LF
    from calculations.figures import save_panels
    save_panels(LF.single_layer_panels(...), "Level1_Perfect_Metal", FIGS, module=LF.MODULE)
"""
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from calculations import house_style as hs
from calculations.figures import Panel
from calculations.research_style import ANNOT_PT, COLOURS, CURVE_PT

MODULE = "calculations/level_figures.py"
BLACK = COLOURS["black"]
GREY = hs.GREY
#: conditions of one figure (temperatures, defect fractions): marker and line style, in order
CONDITION_MARKS = [("o", "-"), ("s", (0, (4, 2))), ("^", ":"), ("D", "-."), ("v", (0, (1, 3)))]

FLUX = r"Flux $J$ (mol m$^{-2}$ s$^{-1}$)"
P_UP = r"Upstream pressure $P_\mathrm{up}$ (Pa)"
UNITS_SIEVERTS = r"(mol m$^{-1}$ s$^{-1}$ Pa$^{-1/2}$)"
UNITS_HENRY = r"(mol m$^{-1}$ s$^{-1}$ Pa$^{-1}$)"


# --------------------------------------------------------------------------------------------------
# Small pieces
# --------------------------------------------------------------------------------------------------

def _tex(s):
    """A name for a LaTeX caption: underscores and per cents escaped."""
    return str(s).replace("\\", "").replace("_", r"\_").replace("%", r"\%")


def _alloy(name):
    return _tex(hs.ALLOY_NAMES.get(name, name))


def _length(m):
    """A length for a caption in nm, µm or mm: '500~nm', '100~$\\mu$m', '1~mm'."""
    for scale, unit in ((1e9, "nm"), (1e6, r"$\mu$m"), (1e3, "mm")):
        if m * scale < 1000 or unit == "mm":
            return f"{hs.tex_num(m * scale, 3)}~{unit}"


def _degC(T):
    return f"{T - 273.15:.0f}~$^\\circ$C"


def _minus(text):
    return text.replace("-", "−")


def _style(level, **over):
    return {**hs.level_style(level), **over}


def _every(x):
    return max(1, len(np.atleast_1d(x)) // 10)


def _draw(ax, x, y, st, ms=9, lw=CURVE_PT, markers=True, z=3):
    """One series: a line in the style's colour and line style, markers every tenth of the points.
    Returns the keyword arguments of its key handle."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    kw = dict(color=st["color"], ls=st.get("ls", "-"), lw=lw)
    marker = st.get("marker")
    if markers and marker not in (None, "None", ""):
        kw.update(marker=marker, ms=ms, mfc=st.get("mfc", st["color"]), mec=st["color"], mew=1.6)
    ax.plot(x, y, zorder=z, markevery=_every(x) if "marker" in kw else None, **kw)
    return kw


def _points(ax, x, y, st, ms=10, z=4):
    """Points without a line (a parity check)."""
    kw = dict(ls="none", marker=st.get("marker", "o"), ms=ms, color=st["color"],
              mfc=st.get("mfc", st["color"]), mec=st["color"], mew=1.6)
    ax.plot(np.asarray(x, float), np.asarray(y, float), zorder=z, **kw)
    return kw


def _handle(kw):
    return Line2D([], [], **kw)


def _ref(ax, x, y, z=2, **over):
    """A grey dashed fit, theory line or limit; z=6 draws a fit over the data it fits, so its dashes show."""
    kw = {**hs.REF_LINE, **over}
    ax.plot(np.asarray(x, float), np.asarray(y, float), zorder=z, **kw)
    return kw


def _set_y(ax, values, **kw):
    """house_style.set_log_y, with powers of ten on a narrow linear axis of very small or large values
    (a flux of 2.9157e-6 to 2.9158e-6 reads 2.9157 ×10⁻⁶, not 0.0000029157)."""
    hs.set_log_y(ax, values, **kw)
    if ax.get_yscale() == "linear":
        lo, hi = ax.get_ylim()
        if max(abs(lo), abs(hi)) < 1e-2 or max(abs(lo), abs(hi)) >= 1e4:
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
        ax.set_yticks([t for t in ax.get_yticks() if lo <= t <= hi])     # no tick label beyond the limits
        ax.set_ylim(lo, hi)


def _ratio(r):
    """A ratio for a key: '1.0000', '0.7381', or 2.8×10⁻⁴ (mathtext)."""
    return f"{r:.4f}" if 0.01 <= abs(r) < 100 else hs.sci(r, 3)


def _local_slope(x, y, end, n=4):
    """d log y / d log x over the first or last n points."""
    x, y = np.log10(np.asarray(x, float)), np.log10(np.asarray(y, float))
    sl = slice(0, n) if end == "low" else slice(-n, None)
    return float(np.polyfit(x[sl], y[sl], 1)[0])


def _arrhenius_x(ax, T_K, label=r"1000/$T$ (K$^{-1}$)"):
    """x = 1000/T with a top axis in °C, padded so the end markers stay inside."""
    x = 1000.0 / np.asarray(T_K, float)
    span = x.max() - x.min()
    ax.set_xlim(x.min() - 0.05 * span, x.max() + 0.05 * span)
    ax.set_xlabel(label)
    top = hs.celsius_top_axis(ax)
    hs.set_celsius_ticks(ax, top)
    top.tick_params(axis="x", direction="in", which="both")
    return x


def _log_axes(ax, x, y, xpad=0.15, ypad=0.3):
    hs.set_log_x(ax, x, pad_decades=xpad)
    _set_y(ax, y, pad_decades=ypad)


def _note(ax, lines):
    """The panel's check result, in the corner that covers the fewest data points, on a white backing.
    Call after the axes limits and the key are final."""
    lines = [ln for ln in lines if ln]
    if not lines:
        return None
    fig = ax.figure
    fig.canvas.draw()                                  # the constrained layout has placed the axes
    renderer = fig.canvas.get_renderer()
    pts = []
    for line in ax.get_lines():
        xd, yd = line.get_xdata(orig=False), line.get_ydata(orig=False)
        try:
            xd, yd = np.asarray(xd, float), np.asarray(yd, float)
        except (TypeError, ValueError):
            continue
        if xd.size and xd.size == yd.size:
            ok = np.isfinite(xd) & np.isfinite(yd)
            if ax.get_xscale() == "log":
                ok &= xd > 0
            if ax.get_yscale() == "log":
                ok &= yd > 0
            if ok.any():
                pts.append(hs._to_axes(ax, xd[ok], yd[ok]))
    pts = np.vstack(pts) if pts else np.empty((0, 2))
    t = ax.text(0, 0, "\n".join(lines), transform=ax.transAxes, fontsize=ANNOT_PT, fontweight="bold",
                color=BLACK, linespacing=1.25, zorder=7)
    best = None
    for x, y, ha, va in ((0.03, 0.96, "left", "top"), (0.97, 0.04, "right", "bottom"),
                         (0.97, 0.96, "right", "top"), (0.03, 0.04, "left", "bottom")):
        t.set_position((x, y))
        t.set_ha(ha)
        t.set_va(va)
        b = hs.bbox_frac(ax, t, renderer)
        covered = int(((pts[:, 0] > b.x0 - 0.01) & (pts[:, 0] < b.x1 + 0.01)
                       & (pts[:, 1] > b.y0 - 0.01) & (pts[:, 1] < b.y1 + 0.01)).sum()) if len(pts) else 0
        if best is None or covered < best[0]:
            best = (covered, (x, y, ha, va))
    x, y, ha, va = best[1]
    t.set_position((x, y))
    t.set_ha(ha)
    t.set_va(va)
    t.set_bbox(dict(facecolor="white", edgecolor="none", alpha=0.9, pad=2))
    return t


def _key(fig, entries, ncol=3):
    """A key below the plot from [(handle kwargs or artist, label)], text in the series colour."""
    handles, labels, colours = [], [], []
    for h, lab in entries:
        if isinstance(h, dict):
            colours.append(h.get("color", BLACK))
            h = _handle(h)
        else:
            colours.append(getattr(h, "get_facecolor", lambda: BLACK)() if isinstance(h, Patch) else BLACK)
        handles.append(h)
        labels.append(lab)
    colours = [c if isinstance(c, str) else BLACK for c in colours]
    return hs.key_below(fig, handles, labels, ncol=ncol, colours=colours)


def _rel_err(model, ref):
    model, ref = np.asarray(model, float), np.asarray(ref, float)
    return float(np.nanmax(np.abs(model - ref) / np.abs(ref) * 100))


def _err_text(v):
    """A percentage for the figure: 0, 0.49 or 2.8×10⁻¹⁴ (mathtext)."""
    return hs.sci(v, 2)


def _parity(name, series, xlabel, ylabel, caption, check_lines):
    """Model against a reference (an analytical solution or a lower level): points per series and the
    grey dashed line of perfect agreement, equal log axes. series: [(label, style, x_ref, y_model)]."""
    fig, ax = hs.new_figure()
    entries, xs, ys = [], [], []
    for label, st, xr, ym in series:
        kw = _points(ax, xr, ym, st)
        entries.append((kw, label))
        xs.append(np.asarray(xr, float))
        ys.append(np.asarray(ym, float))
    allv = np.concatenate(xs + ys)
    lo, hi = hs.log_limits(allv, pad_decades=0.3)
    ref = _ref(ax, [lo, hi], [lo, hi], z=6)
    entries.append((ref, "Perfect agreement"))
    for setter, axis in ((ax.set_xscale, ax.xaxis), (ax.set_yscale, ax.yaxis)):
        setter("log")
        hs.log_ticks(axis, lo, hi)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    _key(fig, entries, ncol=3)
    _note(ax, check_lines)
    return Panel(name, fig, ax, caption, np.concatenate(xs), np.concatenate(ys), keyed=True)


# --------------------------------------------------------------------------------------------------
# Level 1 and Level 2a: one layer (a perfect metal, or a perfect oxide alone)
# --------------------------------------------------------------------------------------------------

def single_layer_panels(level, pressure, temperature, thickness, parity, conditions):
    """Level 1 (perfect metal, Sieverts + Fick) or Level 2a (perfect oxide alone, Henry + Fick).

    pressure     dict(P, J, slope, intercept, r2, expected): the pressure sweep and its log-log fit
    temperature  dict(T, phi, E, E_expected, r2, slope, intercept): the permeability (J/mol energies;
                 ln Φ = intercept + slope · 1000/T)
    thickness    dict(L, J, J_theory, slope, r2, expected): the thickness sweep (L in m) and J ∝ 1/L
    parity       [(T, J_analytical, J_model)]: the model against the analytical flux at each T
    conditions   dict(name, T, P_ref, L, P_down): the layer's name, the operating temperature (K), the
                 pressure of the thickness sweep (Pa), the layer thickness (m), the downstream pressure

    Panels: flux_vs_pressure, permeability_vs_T, flux_vs_thickness, analytical_check."""
    metal = level == "L1"
    st = hs.level_style(level)
    c = conditions
    law = "Sieverts" if metal else "Henry"
    what = (f"the perfect metal {_alloy(c['name'])} (Level~1, Sieverts' law and Fick's law)" if metal else
            f"the perfect oxide {_tex(c['name'])} alone (Level~2a, Henry's law and Fick's law)")
    thick = (f"{c['L'] * 1e3:.3g}~mm thick" if metal else f"{c['L'] * 1e6:.3g}~$\\mu$m thick")
    out = []

    # flux against pressure
    p = pressure
    fig, ax = hs.new_figure()
    kw = _draw(ax, p["P"], p["J"], st)
    fit = 10 ** (p["slope"] * np.log10(p["P"]) + p["intercept"])
    ref = _ref(ax, p["P"], fit, z=6)
    _log_axes(ax, p["P"], np.r_[p["J"], fit])
    ax.set_xlabel(P_UP)
    ax.set_ylabel(FLUX)
    _key(fig, [(kw, "Model"), (ref, "Log-log fit")])
    ok = abs(p["slope"] - p["expected"]) < 0.01
    _note(ax, [_minus(f"slope {p['slope']:.3f}, expected {p['expected']:.3f}"), f"$R^2$ = {p['r2']:.4f}",
               None if ok else "check failed"])
    out.append(Panel("flux_vs_pressure", fig, ax,
                     f"Steady-state flux through {what}, {thick}, against the upstream pressure at {_degC(c['T'])} "
                     f"({hs.tex_num(p['P'][0])} to {hs.tex_num(p['P'][-1])}~Pa, {len(p['P'])} points). The grey dashed "
                     f"line is the log-log fit: slope {p['slope']:.4f} against the {p['expected']:.1f} of {law}'s law "
                     f"($R^2$ = {p['r2']:.6f}).",
                     p["P"], np.r_[p["J"]], keyed=True))

    # permeability against 1000/T
    t = temperature
    fig, ax = hs.new_figure()
    x = 1000.0 / np.asarray(t["T"], float)
    kw = _draw(ax, x, t["phi"], st)
    fit = np.exp(t["slope"] * x + t["intercept"])
    ref = _ref(ax, x, fit, z=6)
    _arrhenius_x(ax, t["T"])
    _set_y(ax, np.r_[t["phi"], fit], pad_decades=0.3)
    sym = r"$\Phi$" if metal else r"$\Phi_\mathrm{ox}$"
    ax.set_ylabel(f"Permeability {sym} {UNITS_SIEVERTS if metal else UNITS_HENRY}")
    _key(fig, [(kw, "Model" if metal else r"$\Phi_\mathrm{ox} = D_\mathrm{ox} K_\mathrm{ox}$"), (ref, "Arrhenius fit")])
    E, Ex = t["E"] / 1e3, t["E_expected"] / 1e3
    ok = abs(E - Ex) < 0.01 * max(1.0, abs(Ex))
    esym = r"$E_\Phi$" if metal else r"$E_{\Phi,\mathrm{ox}}$"
    _note(ax, [f"{esym} = {hs.kj(E, 2)}", _minus(f"expected {Ex:.2f}, ") + f"$R^2$ = {t['r2']:.4f}",
               None if ok else "check failed"])
    out.append(Panel("permeability_vs_T", fig, ax,
                     f"Permeability of {what} against inverse temperature ({t['T'][0]:.0f} to {t['T'][-1]:.0f}~K, "
                     f"{len(t['T'])} points). The grey dashed line is the Arrhenius fit: "
                     f"{esym} = {hs.tex_kj(E, 2)} against {hs.tex_kj(Ex, 2)} expected from "
                     + ("$E_D + \\Delta H_s$" if metal else "$E_{D,\\mathrm{ox}} + \\Delta H_\\mathrm{sol,ox}$")
                     + f" ($R^2$ = {t['r2']:.6f}).",
                     x, np.asarray(t["phi"], float), keyed=True))

    # flux against thickness
    k = thickness
    scale, unit = (1e3, "mm") if metal else (1e6, r"$\mu$m")
    L = np.asarray(k["L"], float) * scale
    fig, ax = hs.new_figure()
    kw = _draw(ax, L, k["J"], st)
    ref = _ref(ax, L, k["J_theory"], z=6)
    _log_axes(ax, L, np.r_[k["J"], k["J_theory"]])
    ax.set_xlabel((r"Metal thickness $L$" if metal else r"Oxide thickness $\delta$") + f" ({unit})")
    ax.set_ylabel(FLUX)
    _key(fig, [(kw, "Model"), (ref, r"Theory, $J \propto 1/L$" if metal else r"Theory, $J \propto 1/\delta$")])
    ok = abs(k["slope"] - k["expected"]) < 0.01
    _note(ax, [_minus(f"slope {k['slope']:.3f}, expected {k['expected']:.3f}"), f"$R^2$ = {k['r2']:.4f}",
               None if ok else "check failed"])
    out.append(Panel("flux_vs_thickness", fig, ax,
                     f"Steady-state flux through {what} against its thickness at {_degC(c['T'])} and "
                     f"{hs.tex_num(c['P_ref'])}~Pa upstream ({L[0]:.3g} to {L[-1]:.3g}~{unit}). The grey dashed line is "
                     f"$J \\propto 1/{'L' if metal else chr(92) + 'delta'}$ through the middle point; the log-log slope is "
                     f"{k['slope']:.4f} against {k['expected']:.0f} ($R^2$ = {k['r2']:.6f}).",
                     L, np.asarray(k["J"], float), keyed=True))

    # model against the analytical solution
    series, errs = [], []
    for i, (T, Ja, Jm) in enumerate(parity):
        marker, _ = CONDITION_MARKS[i % len(CONDITION_MARKS)]
        series.append((f"{T - 273.15:.0f} °C", dict(color=BLACK, marker=marker, mfc="white" if i % 2 else BLACK), Ja, Jm))
        errs.append(_rel_err(Jm, Ja))
    emax = max(errs)
    formula = (r"$J = \Phi \sqrt{P}/L$" if metal else r"$J = D_\mathrm{ox} K_\mathrm{ox} \Delta P / \delta$")
    out.append(_parity("analytical_check", series, r"Analytical flux (mol m$^{-2}$ s$^{-1}$)",
                       r"Model flux (mol m$^{-2}$ s$^{-1}$)",
                       f"The model's flux through {what} against the analytical solution {formula}, at "
                       f"{hs.tex_join([_degC(T) for T, _, _ in parity])} over the pressure range, {thick}. The grey dashed "
                       f"line is perfect agreement; the largest relative difference is {hs.tex_num(emax)}~\\%.",
                       ["exact match" if emax < 1e-6 else f"max error {_err_text(emax)} %",
                        None if emax < 0.01 else "check failed"]))
    return out


# --------------------------------------------------------------------------------------------------
# Level 2b: perfect oxide + perfect metal in series
# --------------------------------------------------------------------------------------------------

def level2b_panels(pressure, temperature, thickness, limit, conditions):
    """Level 2b (perfect oxide on a perfect metal, coupled at the interface pressure).

    pressure     dict(P, J2b, J1, J2a): Level 2b against Level 1 and Level 2a over the pressure range
    temperature  dict(T, J2b, J1, P, E_app, r2): against 1000/T at pressure P (E_app in J/mol)
    thickness    dict(L, J2b, J1, J_oxide_limit, delta, P, L_transition): against the metal thickness (m),
                 for an oxide delta (m) thin enough to show both regimes; L_transition (m) or None
    limit        dict(delta, P_norm, delta_99): the normalised interface pressure as the oxide thins (m);
                 delta_99 (m) where it passes 0.99, or None
    conditions   dict(metal, oxide, T, L_oxide, L_metal)

    Panels: flux_vs_pressure, flux_vs_T, flux_vs_metal_thickness, limit_check."""
    c = conditions
    s1, s2a, s2b = hs.level_style("L1"), hs.level_style("L2a"), hs.level_style("L2b")
    system = (f"a {c['L_oxide'] * 1e6:.3g}~$\\mu$m {_tex(c['oxide'])} oxide on {c['L_metal'] * 1e3:.3g}~mm of "
              f"{_alloy(c['metal'])}")
    out = []

    p = pressure
    fig, ax = hs.new_figure()
    k2b = _draw(ax, p["P"], p["J2b"], s2b, z=5)
    k1 = _draw(ax, p["P"], p["J1"], {**s1, "ls": (0, (4, 2))})
    k2a = _draw(ax, p["P"], p["J2a"], {**s2a, "ls": ":"})
    _log_axes(ax, p["P"], np.r_[p["J2b"], p["J1"], p["J2a"]])
    ax.set_xlabel(P_UP)
    ax.set_ylabel(FLUX)
    _key(fig, [(k2b, "L2b, oxide + metal"), (k1, "L1, metal only"), (k2a, "L2a, oxide only")])
    lo, hi = _local_slope(p["P"], p["J2b"], "low"), _local_slope(p["P"], p["J2b"], "high")
    _note(ax, [_minus(f"L2b slope {lo:.2f} at low pressure"), _minus(f"{hi:.2f} at high pressure")])
    out.append(Panel("flux_vs_pressure", fig, ax,
                     f"Steady-state flux through {system} (Level~2b) against the upstream pressure at {_degC(c['T'])}, "
                     f"with the metal alone (Level~1, slope 0.5) and the oxide alone (Level~2a, slope 1). The Level~2b "
                     f"log-log slope is {lo:.2f} over the lowest pressures and {hi:.2f} over the highest: the oxide "
                     f"limits the flux at low pressure and the metal at high pressure.",
                     np.r_[p["P"], p["P"], p["P"]], np.r_[p["J2b"], p["J1"], p["J2a"]], keyed=True))

    t = temperature
    fig, ax = hs.new_figure()
    x = 1000.0 / np.asarray(t["T"], float)
    k2b = _draw(ax, x, t["J2b"], s2b, z=5)
    k1 = _draw(ax, x, t["J1"], {**s1, "ls": (0, (4, 2))})
    _arrhenius_x(ax, t["T"])
    _set_y(ax, np.r_[t["J2b"], t["J1"]], pad_decades=0.3)
    ax.set_ylabel(FLUX)
    E = t["E_app"] / 1e3
    _key(fig, [(k2b, f"L2b, oxide + metal\n$E_\\mathrm{{app}}$ = {hs.kj(E, 1)}"), (k1, "L1, metal only")], ncol=2)
    _note(ax, [f"$E_\\mathrm{{app}}$ = {hs.kj(E, 2)}", f"$R^2$ = {t['r2']:.4f}"])
    out.append(Panel("flux_vs_T", fig, ax,
                     f"Steady-state flux through {system} (Level~2b) and through the metal alone (Level~1) against "
                     f"inverse temperature at {hs.tex_num(t['P'])}~Pa upstream ({t['T'][0]:.0f} to {t['T'][-1]:.0f}~K). "
                     f"A single Arrhenius fit to Level~2b gives an apparent activation energy of {hs.tex_kj(E, 2)} "
                     f"($R^2$ = {t['r2']:.4f}); the curvature shows the limiting layer changing with temperature.",
                     np.r_[x, x], np.r_[t["J2b"], t["J1"]], keyed=True))

    k = thickness
    L = np.asarray(k["L"], float) * 1e3
    fig, ax = hs.new_figure()
    k2b = _draw(ax, L, k["J2b"], s2b, z=5)
    k1 = _draw(ax, L, k["J1"], {**s1, "ls": (0, (4, 2))})
    lim = _ref(ax, [L[0], L[-1]], [k["J_oxide_limit"]] * 2, color=s2a["color"], ls=":", lw=CURVE_PT)
    entries = [(k2b, "L2b, oxide + metal"), (k1, "L1, metal only"), (lim, "Oxide limit")]
    if k.get("L_transition"):
        Lt = k["L_transition"] * 1e3
        tr = dict(color=GREY, ls="-.", lw=2.4)
        ax.axvline(Lt, zorder=2, **tr)
        entries.append((tr, f"Transition, {Lt:.3g} mm"))
    _log_axes(ax, L, np.r_[k["J2b"], k["J1"], k["J_oxide_limit"]])
    ax.set_xlabel(r"Metal thickness $L$ (mm)")
    ax.set_ylabel(FLUX)
    _key(fig, entries, ncol=2)
    out.append(Panel("flux_vs_metal_thickness", fig, ax,
                     f"Steady-state flux against the metal thickness for a thin {k['delta'] * 1e9:.3g}~nm {_tex(c['oxide'])} "
                     f"oxide on {_alloy(c['metal'])} at {_degC(c['T'])} and {hs.tex_num(k['P'])}~Pa, chosen to show both "
                     f"regimes: Level~2b, the metal alone (Level~1) and the oxide-limited flux "
                     f"{hs.tex_num(k['J_oxide_limit'])}~mol~m$^{{-2}}$~s$^{{-1}}$ (dotted)."
                     + (f" Level~2b falls 10\\,\\% below the oxide limit at $L \\approx$ {k['L_transition'] * 1e3:.3g}~mm "
                        "(grey dash-dot line)." if k.get("L_transition") else
                        " Level~2b stays within 10\\,\\% of the oxide limit over the whole range."),
                     np.r_[L, L], np.r_[k["J2b"], k["J1"]], keyed=True))

    m = limit
    d = np.asarray(m["delta"], float) * 1e6
    fig, ax = hs.new_figure()
    kn = _draw(ax, d, m["P_norm"], s2b, z=5)
    one = _ref(ax, [d.min(), d.max()], [1.0, 1.0], color=s1["color"])
    entries = [(kn, r"L2b, $(P_\mathrm{int} - P_\mathrm{down})/(P_\mathrm{up} - P_\mathrm{down})$"),
               (one, r"L1 limit, $P_\mathrm{int} = P_\mathrm{up}$")]
    if m.get("delta_99"):
        d99 = m["delta_99"] * 1e6
        th = dict(color=GREY, ls=":", lw=2.4)
        ax.axvline(d99, zorder=2, **th)
        entries.append((th, f"Ratio above 0.99 below {hs.sci(d99, 2)} µm"))
    hs.set_log_x(ax, d, pad_decades=0.15)
    ax.set_ylim(0, 1.1)
    ax.invert_xaxis()                                       # thin oxide on the right: delta -> 0
    ax.set_xlabel(r"Oxide thickness $\delta$ ($\mu$m)")
    ax.set_ylabel(r"$P_\mathrm{int}$, normalised")
    _key(fig, entries, ncol=1)
    final = float(np.asarray(m["P_norm"])[-1])
    ok = final > 0.99
    _note(ax, [r"$\delta \to 0$: " + f"ratio {final:.4f}", "metal dominates, as Level 1" if ok else "check failed: expected 1"])
    out.append(Panel("limit_check", fig, ax,
                     f"Limit check of Level~2b: the interface pressure, normalised as "
                     f"$(P_\\mathrm{{int}} - P_\\mathrm{{down}})/(P_\\mathrm{{up}} - P_\\mathrm{{down}})$, as the "
                     f"{_tex(c['oxide'])} oxide thins from {d[0]:.3g} to {hs.tex_num(d[-1])}~$\\mu$m on "
                     f"{c['L_metal'] * 1e3:.3g}~mm of {_alloy(c['metal'])} at {_degC(c['T'])}. As $\\delta \\to 0$ it "
                     f"reaches {final:.4f}" + (" (the dashed line is 1): the metal becomes the only resistance, as in "
                                              "Level~1." if ok else ", short of the 1 expected: the check fails.")
                     + (f" It passes 0.99 below $\\delta$ = {hs.tex_num(m['delta_99'] * 1e6)}~$\\mu$m." if m.get("delta_99") else ""),
                     d, np.asarray(m["P_norm"], float), keyed=True))
    return out


# --------------------------------------------------------------------------------------------------
# Level 3: defective oxide + perfect metal (parallel paths)
# --------------------------------------------------------------------------------------------------

def level3_panels(pressure, fraction, prf, limit, conditions):
    """Level 3 (a defective oxide on a perfect metal: intact and pinhole paths in parallel).

    pressure    dict(P, curves=[(f, J)], J1): the flux for each defect fraction f (f = 0 is Level 2b) and
                Level 1, against pressure
    fraction    dict(f, J_total, J_defect, J_intact, J_L2b, J_L1, f_crossover, P): against the defect
                fraction at pressure P; f_crossover (fraction) where the defect path takes over, or None
    prf         dict(f, PRF, PRF_perfect, f_half): the permeation reduction factor J(L1)/J
    limit       dict(T, J_L2b, J_L3, f_check): Level 3 at a vanishing defect fraction against Level 2b
    conditions  dict(metal, oxide, T, L_oxide, L_metal, defect_type)

    Panels: flux_vs_pressure, flux_vs_defect_fraction, prf, limit_check."""
    c = conditions
    s1, s2b, s3 = hs.level_style("L1"), hs.level_style("L2b"), hs.level_style("L3")
    system = (f"a {c['L_oxide'] * 1e6:.3g}~$\\mu$m {_tex(c['oxide'])} oxide with {c.get('defect_type', 'pinhole')}s on "
              f"{c['L_metal'] * 1e3:.3g}~mm of {_alloy(c['metal'])}")
    out = []

    p = pressure
    fig, ax = hs.new_figure()
    entries, ys = [], []
    k = 0
    for f, J in p["curves"]:
        if f == 0:
            kw = _draw(ax, p["P"], J, s2b, z=5)
            entries.append((kw, "f = 0 (L2b)"))
        else:
            marker, ls = CONDITION_MARKS[(k + 1) % len(CONDITION_MARKS)]
            kw = _draw(ax, p["P"], J, {**s3, "marker": marker, "ls": ls}, z=4)
            entries.append((kw, f"f = {100 * f:g} %"))
            k += 1
        ys.append(np.asarray(J, float))
    k1 = _draw(ax, p["P"], p["J1"], {**s1, "ls": (0, (4, 2))})
    entries.append((k1, "L1, metal only"))
    ys.append(np.asarray(p["J1"], float))
    _log_axes(ax, p["P"], np.concatenate(ys))
    ax.set_xlabel(P_UP)
    ax.set_ylabel(FLUX)
    _key(fig, entries, ncol=3)
    slopes = [(f, _local_slope(p["P"], J, "low")) for f, J in p["curves"]]
    out.append(Panel("flux_vs_pressure", fig, ax,
                     f"Steady-state flux through {system} (Level~3) against the upstream pressure at {_degC(c['T'])}, for "
                     f"defect area fractions {hs.tex_join([f'{100 * f:g}' for f, _ in p['curves']])}~\\% ($f$ = 0 is "
                     f"Level~2b), with the metal alone (Level~1). Over the lowest pressures the log-log slope is "
                     + hs.tex_join([f"{s:.2f} ($f$ = {100 * f:g}~\\%)" for f, s in slopes])
                     + ": defects move the flux from the oxide's slope 1 toward the metal's 0.5.",
                     np.tile(np.asarray(p["P"], float), len(ys)), np.concatenate(ys), keyed=True))

    q = fraction
    f100 = np.asarray(q["f"], float) * 100
    fig, ax = hs.new_figure()
    kt = _draw(ax, f100, q["J_total"], s3, z=5)
    kd = _draw(ax, f100, q["J_defect"], {**s3, "ls": (0, (4, 2)), "marker": "s"})
    ki = _draw(ax, f100, q["J_intact"], {**s2b, "ls": ":", "marker": "^"})
    r2 = _ref(ax, [f100[0], f100[-1]], [q["J_L2b"]] * 2, color=s2b["color"], ls="-.", lw=2.4)
    r1 = _ref(ax, [f100[0], f100[-1]], [q["J_L1"]] * 2, color=s1["color"], ls="-.", lw=2.4)
    entries = [(kt, "Total flux"), (kd, "Defect path"), (ki, "Intact path"), (r2, "L2b (f = 0)"), (r1, "L1 (f = 1)")]
    if q.get("f_crossover"):
        fc = q["f_crossover"] * 100
        cr = dict(color=GREY, ls=(0, (5, 3)), lw=2.4)
        ax.axvline(fc, zorder=2, **cr)
        entries.append((cr, f"Crossover, f = {fc:.3g} %"))
    yall = np.r_[q["J_total"], q["J_defect"], q["J_intact"], q["J_L2b"], q["J_L1"]]
    _log_axes(ax, f100, yall)
    ax.set_xlabel("Defect area fraction (%)")
    ax.set_ylabel(FLUX)
    _key(fig, entries, ncol=3)
    out.append(Panel("flux_vs_defect_fraction", fig, ax,
                     f"Steady-state flux through {system} (Level~3) against the defect area fraction at {_degC(c['T'])} "
                     f"and {hs.tex_num(q['P'])}~Pa ({hs.tex_num(f100[0])} to {f100[-1]:.3g}~\\%): the total, the "
                     f"defect-path and the intact-path contributions, with Level~2b ($f$ = 0, "
                     f"{hs.tex_num(q['J_L2b'])}) and Level~1 ($f$ = 1, {hs.tex_num(q['J_L1'])}~mol~m$^{{-2}}$~s$^{{-1}}$) "
                     "as dash-dot lines."
                     + (f" The defect path carries more than the intact path above $f \\approx$ {q['f_crossover'] * 100:.3g}~\\%."
                        if q.get("f_crossover") else ""),
                     np.tile(f100, 3), np.r_[q["J_total"], q["J_defect"], q["J_intact"]], keyed=True))

    r = prf
    f100 = np.asarray(r["f"], float) * 100
    fig, ax = hs.new_figure()
    kp = _draw(ax, f100, r["PRF"], s3, z=5)
    rp = _ref(ax, [f100[0], f100[-1]], [r["PRF_perfect"]] * 2, color=s2b["color"], lw=2.4)
    r1 = _ref(ax, [f100[0], f100[-1]], [1.0, 1.0])
    hs.set_log_x(ax, f100, pad_decades=0.15)
    hs.set_linear_y(ax, np.r_[r["PRF"], r["PRF_perfect"], 1.0], include=[0.0], floor=0.0)
    ax.set_xlabel("Defect area fraction (%)")
    ax.set_ylabel(r"PRF = $J_\mathrm{L1}/J$")
    _key(fig, [(kp, "With defects"), (rp, f"Perfect oxide, PRF = {hs.sci(r['PRF_perfect'], 3)}"), (r1, "No barrier, PRF = 1")],
         ncol=1)
    halves = bool(np.nanmin(r["PRF"]) <= r["PRF_perfect"] / 2)
    _note(ax, [f"PRF max {hs.sci(r['PRF_perfect'], 3)}",
               f"halved at $f \\approx$ {hs.sci(r['f_half'] * 100, 2)} %" if halves else "never halved in this range"])
    out.append(Panel("prf", fig, ax,
                     f"Permeation reduction factor of the oxide, PRF = $J_\\mathrm{{L1}}/J$ (the bare metal's flux over "
                     f"the coated one), against the defect area fraction for {system} at {_degC(c['T'])}. The perfect "
                     f"oxide gives PRF = {hs.tex_num(r['PRF_perfect'], 4)}"
                     + (f"; it halves at $f \\approx$ {hs.tex_num(r['f_half'] * 100)}~\\%." if halves else
                        f"; it falls to {hs.tex_num(np.nanmin(r['PRF']), 4)} and does not halve over this range.")
                     + " PRF = 1 is no barrier (grey dashed).",
                     f100, np.asarray(r["PRF"], float), keyed=True))

    m = limit
    err = _rel_err(m["J_L3"], m["J_L2b"])
    out.append(_parity("limit_check", [(f"L3, f = {hs.sci(m['f_check'], 1)}", s3, m["J_L2b"], m["J_L3"])],
                       r"Level 2b flux (mol m$^{-2}$ s$^{-1}$)", r"Level 3 flux (mol m$^{-2}$ s$^{-1}$)",
                       f"Limit check of Level~3: with a vanishing defect fraction ($f$ = {hs.tex_num(m['f_check'])}) its "
                       f"flux against Level~2b at each temperature from {m['T'][0]:.0f} to {m['T'][-1]:.0f}~K, for "
                       f"{system}. The grey dashed line is perfect agreement; the largest relative difference is "
                       f"{hs.tex_num(err)}~\\%" + (": Level~3 recovers Level~2b." if err < 0.01 else ", so the check fails."),
                       [r"$f \to 0$: " + f"max error {_err_text(err)} %", None if err < 0.01 else "check failed"]))
    return out


# --------------------------------------------------------------------------------------------------
# Level 4: defective metal (grain boundaries and traps)
# --------------------------------------------------------------------------------------------------

MODE_STYLE = {"none": "L1", "gb_only": "L4_gb", "trapping_only": "L4_trap", "both": "L4_both"}
MODE_WORDS = {"none": "Perfect lattice (L1)", "gb_only": "Grain boundaries only", "trapping_only": "Trapping only",
              "both": "Both (L4)"}


def level4_panels(pressure, temperature, grain, limit, conditions):
    """Level 4 (a defective metal: grain boundaries speed diffusion, traps slow it).

    pressure     dict(P, J1, Jgb, Jtrap, Jboth, T): each mode against pressure at T
    temperature  dict(T, mod_gb, mod_trap, mod_both, T_crossover): D_eff/D_lattice against temperature
    grain        dict(d, J1, Jgb, enhancement, T, P): against the grain size (m), grain boundaries only
    limit        dict(T, J1, Jnone, P): mode 'none' against Level 1
    conditions   dict(metal, L_metal, grain_size, gb_thickness, trap_density, binding_eV)

    Panels: flux_vs_pressure, modification_vs_T, flux_vs_grain_size, limit_check."""
    c = conditions
    st = {m: hs.level_style(l) for m, l in MODE_STYLE.items()}
    metal = f"{c['L_metal'] * 1e3:.3g}~mm of {_alloy(c['metal'])}"
    micro = (f"grain size {_length(c['grain_size'])}, grain-boundary width {_length(c['gb_thickness'])}, "
             f"traps of {c['binding_eV']:.2f}~eV at {hs.tex_num(c['trap_density'])}~m$^{{-3}}$")
    out = []

    p = pressure
    fig, ax = hs.new_figure()
    ratios = {"gb_only": np.mean(p["Jgb"] / p["J1"]), "trapping_only": np.mean(p["Jtrap"] / p["J1"]),
              "both": np.mean(p["Jboth"] / p["J1"])}
    entries = [(_draw(ax, p["P"], p["J1"], st["none"], z=5), "L1, perfect lattice")]
    for mode, key, z in (("gb_only", "Jgb", 4), ("trapping_only", "Jtrap", 4), ("both", "Jboth", 4)):
        kw = _draw(ax, p["P"], p[key], st[mode], z=z)
        entries.append((kw, f"L4, {MODE_WORDS[mode].replace(' (L4)', '').lower()}\n× {_ratio(ratios[mode])} of L1"))
    _log_axes(ax, p["P"], np.r_[p["J1"], p["Jgb"], p["Jtrap"], p["Jboth"]])
    ax.set_xlabel(P_UP)
    ax.set_ylabel(FLUX)
    _key(fig, entries, ncol=2)
    out.append(Panel("flux_vs_pressure", fig, ax,
                     f"Steady-state flux through {metal} against the upstream pressure at {_degC(p['T'])}, for each "
                     f"defective-metal mode (Level~4; {micro}) and the perfect lattice (Level~1). On average over the "
                     f"pressures, grain boundaries alone give {hs.tex_num(ratios['gb_only'], 4)} times the Level~1 flux, "
                     f"trapping alone {hs.tex_num(ratios['trapping_only'], 3)} times and both "
                     f"{hs.tex_num(ratios['both'], 3)} times (the key gives the same factors).",
                     np.tile(np.asarray(p["P"], float), 4), np.r_[p["J1"], p["Jgb"], p["Jtrap"], p["Jboth"]], keyed=True))

    t = temperature
    TC = np.asarray(t["T"], float) - 273.15
    fig, ax = hs.new_figure()
    entries = []
    for mode, key in (("gb_only", "mod_gb"), ("trapping_only", "mod_trap"), ("both", "mod_both")):
        entries.append((_draw(ax, TC, t[key], st[mode]), MODE_WORDS[mode]))
    one = _ref(ax, [TC[0], TC[-1]], [1.0, 1.0])
    entries.append((one, "No modification"))
    span = TC[-1] - TC[0]
    ax.set_xlim(TC[0] - 0.04 * span, TC[-1] + 0.04 * span)
    _set_y(ax, np.r_[t["mod_gb"], t["mod_trap"], t["mod_both"], 1.0], pad_decades=0.2)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel(r"$D_\mathrm{eff}/D_\mathrm{lattice}$")
    _key(fig, entries, ncol=2)
    both = np.asarray(t["mod_both"], float)
    crosses = bool(np.nanmin(both) < 1 < np.nanmax(both))
    _note(ax, [f"combined = 1 at {t['T_crossover'] - 273.15:.0f} °C" if crosses else
               f"combined {hs.sci(np.nanmin(both), 3)} to {hs.sci(np.nanmax(both), 3)}"])
    out.append(Panel("modification_vs_T", fig, ax,
                     f"Effective diffusivity of {_alloy(c['metal'])} relative to its lattice diffusivity, "
                     f"$D_\\mathrm{{eff}}/D_\\mathrm{{lattice}}$, against temperature ({TC[0]:.0f} to {TC[-1]:.0f}~$^\\circ$C; "
                     f"{micro}): grain boundaries alone, trapping alone and both. The combined factor runs from "
                     f"{hs.tex_num(np.nanmin(both), 3)} to {hs.tex_num(np.nanmax(both), 3)}"
                     + (f" and crosses 1 near {t['T_crossover'] - 273.15:.0f}~$^\\circ$C: trapping dominates below and "
                        "grain boundaries above." if crosses else
                        f"; it does not cross 1 in this range (closest to 1 at {t['T_crossover'] - 273.15:.0f}~$^\\circ$C)."),
                     np.tile(TC, 3), np.r_[t["mod_gb"], t["mod_trap"], t["mod_both"]], keyed=True))

    g = grain
    d = np.asarray(g["d"], float) * 1e6
    fig, ax = hs.new_figure()
    k1 = _draw(ax, d, np.full_like(d, float(g["J1"])), {**st["none"], "ls": (0, (4, 2))}, markers=False)
    kg = _draw(ax, d, g["Jgb"], st["gb_only"], z=5)
    hs.set_log_x(ax, d, pad_decades=0.1)
    _set_y(ax, np.r_[g["Jgb"], g["J1"]], pad_decades=0.2)
    ax.set_xlabel(r"Grain size ($\mu$m)")
    ax.set_ylabel(FLUX)
    enh = np.asarray(g["enhancement"], float)
    _key(fig, [(k1, "L1, perfect lattice"),
               (kg, f"L4, grain boundaries only\n× {enh[0]:.5f} to {enh[-1]:.5f} of L1")],
         ncol=2)
    out.append(Panel("flux_vs_grain_size", fig, ax,
                     f"Steady-state flux through {metal} with grain boundaries only (no traps) against the grain size "
                     f"({hs.tex_num(d[0])} to {hs.tex_num(d[-1])}~$\\mu$m) at {_degC(g['T'])} and {hs.tex_num(g['P'])}~Pa, "
                     f"with the perfect lattice (Level~1). The grain-boundary volume fraction is $3\\delta_\\mathrm{{gb}}/d$, "
                     f"so smaller grains give more fast paths; the flux is {hs.tex_num(enh[0], 4)} times Level~1 at the "
                     f"smallest grains and {hs.tex_num(enh[-1], 4)} times at the largest.",
                     np.r_[d, d], np.r_[g["Jgb"], np.full_like(d, float(g["J1"]))], keyed=True))

    m = limit
    err = _rel_err(m["Jnone"], m["J1"])
    out.append(_parity("limit_check", [("L4, mode 'none'", st["none"], m["J1"], m["Jnone"])],
                       r"Level 1 flux (mol m$^{-2}$ s$^{-1}$)", r"L4 flux, no defects (mol m$^{-2}$ s$^{-1}$)",
                       f"Limit check of Level~4: with no microstructure (mode `none') its flux against Level~1 at each "
                       f"temperature from {m['T'][0]:.0f} to {m['T'][-1]:.0f}~K, {metal}, {hs.tex_num(m['P'])}~Pa. The "
                       f"grey dashed line is perfect agreement; the largest relative difference is {hs.tex_num(err)}~\\%"
                       + (": Level~4 recovers Level~1." if err < 0.01 else ", so the check fails."),
                       [f"max error {_err_text(err)} %", None if err < 0.01 else "check failed"]))
    return out


# --------------------------------------------------------------------------------------------------
# Level 5: defective oxide + defective metal
# --------------------------------------------------------------------------------------------------

def level5_panels(pressure, temperature, modes, limit, conditions):
    """Level 5 (the full system: a defective oxide on a defective metal).

    pressure     dict(P, J1, J2b, J3, J4, J5, T): every level against pressure at T
    temperature  dict(T, J1, J2b, J5, E1, E5, P): against 1000/T; E1, E5 the Arrhenius energies (kJ/mol)
    modes        dict(f, none, gb, trap, both, T, P): Level 5 for each oxide defect fraction and metal mode
    limit        dict(J2b, J5_f0, J1, J5_f1, f0, f1, T): Level 5 at f -> 0 and f -> 1 (mode 'none')
                 against Level 2b and Level 1
    conditions   dict(metal, oxide, L_oxide, L_metal, defect_fraction, grain_size, trap_density)

    Panels: flux_vs_pressure, flux_vs_T, defects_and_microstructure, limit_check."""
    c = conditions
    system = (f"a {c['L_oxide'] * 1e6:.3g}~$\\mu$m {_tex(c['oxide'])} oxide with {100 * c['defect_fraction']:.3g}~\\% "
              f"pinholes on {c['L_metal'] * 1e3:.3g}~mm of {_alloy(c['metal'])} (grain size {_length(c['grain_size'])}, "
              f"traps at {hs.tex_num(c['trap_density'])}~m$^{{-3}}$)")
    out = []

    p = pressure
    fig, ax = hs.new_figure()
    levels = [("L1", "J1", "L1, metal only"), ("L2b", "J2b", "L2b, oxide + metal"), ("L3", "J3", "L3, defective oxide"),
              ("L4_both", "J4", "L4, defective metal"), ("L5", "J5", "L5, full system")]
    entries = []
    for lv, key, label in levels:
        entries.append((_draw(ax, p["P"], p[key], hs.level_style(lv), z=5 if lv == "L5" else 3), label))
    yall = np.r_[tuple(np.asarray(p[k], float) for _, k, _ in levels)]
    _log_axes(ax, p["P"], yall)
    ax.set_xlabel(P_UP)
    ax.set_ylabel(FLUX)
    _key(fig, entries, ncol=3)
    i = int(np.argmin(np.abs(np.log10(np.asarray(p["P"], float)) - 5)))
    out.append(Panel("flux_vs_pressure", fig, ax,
                     f"Steady-state flux of every level against the upstream pressure at {_degC(p['T'])}, for {system}. "
                     f"At {hs.tex_num(p['P'][i])}~Pa the levels give "
                     + hs.tex_join([f"{hs.tex_num(p[k][i])} ({lab.split(',')[0]})" for _, k, lab in levels])
                     + "~mol~m$^{-2}$~s$^{-1}$.",
                     np.tile(np.asarray(p["P"], float), len(levels)), yall, keyed=True))

    t = temperature
    fig, ax = hs.new_figure()
    x = 1000.0 / np.asarray(t["T"], float)
    k1 = _draw(ax, x, t["J1"], hs.level_style("L1"))
    k2 = _draw(ax, x, t["J2b"], hs.level_style("L2b"))
    k5 = _draw(ax, x, t["J5"], hs.level_style("L5"), z=5)
    _arrhenius_x(ax, t["T"])
    _set_y(ax, np.r_[t["J1"], t["J2b"], t["J5"]], pad_decades=0.3)
    ax.set_ylabel(FLUX)
    _key(fig, [(k1, f"L1, metal only\n$E_a$ = {hs.kj(t['E1'], 1)}"), (k2, "L2b, oxide + metal"),
               (k5, f"L5, full system\n$E_a$ = {hs.kj(t['E5'], 1)}")], ncol=3)
    out.append(Panel("flux_vs_T", fig, ax,
                     f"Steady-state flux of Levels~1, 2b and 5 against inverse temperature at {hs.tex_num(t['P'])}~Pa "
                     f"({t['T'][0]:.0f} to {t['T'][-1]:.0f}~K), for {system}. Single Arrhenius fits give "
                     f"$E_a$ = {hs.tex_kj(t['E1'], 2)} for Level~1 and {hs.tex_kj(t['E5'], 2)} for Level~5.",
                     np.r_[x, x, x], np.r_[t["J1"], t["J2b"], t["J5"]], keyed=True))

    m = modes
    f = np.asarray(m["f"], float)
    fig, ax = hs.new_figure()
    xs = np.arange(len(f))
    w = 0.2
    bars = [("none", "No microstructure", hs.level_style("L1")["color"], None),
            ("gb", "Grain boundaries only", hs.level_style("L4_gb")["color"], None),
            ("trap", "Trapping only", hs.level_style("L4_trap")["color"], "//"),
            ("both", "Both (L5)", hs.level_style("L5")["color"], None)]
    entries, vals = [], []
    for j, (key, label, colour, hatch) in enumerate(bars):
        v = np.asarray(m[key], float)
        ax.bar(xs + (j - 1.5) * w, v, w, color=colour, hatch=hatch, edgecolor="black" if hatch is None else "white",
               linewidth=1.2 if hatch is None else 0, zorder=3)
        if hatch:
            ax.bar(xs + (j - 1.5) * w, v, w, fill=False, edgecolor="black", linewidth=1.2, zorder=3)
        entries.append((Patch(facecolor=colour, hatch=hatch, edgecolor="black", linewidth=1.2), label))
        vals.append(v)
    allv = np.concatenate(vals)
    if np.nanmax(allv) / np.nanmin(allv[allv > 0]) < 10:      # bars on a short log axis would mislead: from zero
        ax.set_ylim(0, np.nanmax(allv) * 1.12)
        if np.nanmax(allv) < 1e-2:
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    else:
        lo, hi = hs.log_limits(allv, pad_decades=0.3)
        ax.set_yscale("log")
        ax.set_ylim(lo, hi)
        hs.log_ticks(ax.yaxis, lo, hi)
    ax.set_xlim(-0.6, len(f) - 0.4)
    ax.set_xticks(xs, [f"{100 * v:g}" for v in f])
    ax.tick_params(axis="x", which="both", bottom=False, top=False)
    ax.set_xlabel("Oxide defect area fraction (%)")
    ax.set_ylabel(FLUX)
    hs.key_below(fig, [h for h, _ in entries], [lab for _, lab in entries], ncol=2)
    out.append(Panel("defects_and_microstructure", fig, ax,
                     f"Level~5 flux at {_degC(m['T'])} and {hs.tex_num(m['P'])}~Pa for each oxide defect area fraction "
                     f"({hs.tex_join([f'{100 * v:g}' for v in f])}~\\%) and each treatment of the metal: no "
                     f"microstructure, grain boundaries only, trapping only (hatched) and both, for a "
                     f"{c['L_oxide'] * 1e6:.3g}~$\\mu$m {_tex(c['oxide'])} oxide on {c['L_metal'] * 1e3:.3g}~mm of "
                     f"{_alloy(c['metal'])}. The flux runs from {hs.tex_num(allv.min())} to {hs.tex_num(allv.max())}"
                     "~mol~m$^{-2}$~s$^{-1}$.",
                     keyed=True))

    q = limit
    e0, e1 = _rel_err(q["J5_f0"], q["J2b"]), _rel_err(q["J5_f1"], q["J1"])
    out.append(_parity("limit_check",
                       [(r"L5, $f \to 0$, against L2b", hs.level_style("L2b"), q["J2b"], q["J5_f0"]),
                        (r"L5, $f \to 1$, against L1", {**hs.level_style("L1"), "marker": "s"}, q["J1"], q["J5_f1"])],
                       r"Lower-level flux (mol m$^{-2}$ s$^{-1}$)", r"Level 5 flux (mol m$^{-2}$ s$^{-1}$)",
                       f"Limit checks of Level~5 with no metal microstructure (mode `none'), over the pressure range at "
                       f"{_degC(q['T'])}: with a vanishing defect fraction ($f$ = {hs.tex_num(q['f0'])}) against Level~2b, "
                       f"and with the oxide almost all defect ($f$ = {q['f1']:g}) against Level~1. The grey dashed line is "
                       f"perfect agreement; the largest relative differences are {hs.tex_num(e0)}~\\% and "
                       f"{hs.tex_num(e1)}~\\%.",
                       [r"$f \to 0$: " + f"error {_err_text(e0)} %", r"$f \to 1$: " + f"error {_err_text(e1)} %"]))
    return out
