r"""Fixed text of the study report (LaTeX snippets, not escaped).

Each section has a 'title', a 'goal', 'docs' (where to read more: (path from the
repository root, heading or None, what it explains)), 'columns' (what each table
column means) and 'figures': one entry per kind of figure ('key', used by
:mod:`calculations.report`) with a 'caption' and 'panels': (panel, what it shows, how
to read it, [(what you see, what it means), ...]).

The wording follows the pages it points to and the notebooks' own explanations
(the first cell of each effect-of notebook). Every pointer is checked when the report
is built and by tests/test_report.py, so a renamed page or heading is caught rather
than left as a dead link.
"""

# the levels, in the order of the hierarchy; Model 1 is L1-L5, Model 2 the L6 family
LEVELS = [
    ('L1', 'Level 1 --- perfect metal'),
    ('L2a', 'Level 2a --- perfect oxide'),
    ('L2b', 'Level 2b --- oxide and metal in series'),
    ('L3', 'Level 3 --- defective oxide'),
    ('L4', 'Level 4 --- defective metal (grain boundaries and traps)'),
    ('L5', 'Level 5 --- full system'),
    ('L1L6', 'Level 1 + surface kinetics'),
    ('L2aL6', 'Level 2a + surface kinetics'),
    ('L2L6', 'Level 2 + surface kinetics'),
    ('L2bL6', 'Level 2b + surface kinetics'),
    ('L2L4L6', 'Levels 2 + 4 + surface kinetics'),
    ('L3L6', 'Level 3 + surface kinetics'),
    ('L3L4L6', 'Levels 3 + 4 + surface kinetics'),
    ('L4L6', 'Level 4 + surface kinetics'),
    ('L5L6', 'Level 5 + surface kinetics (the full system)'),
]

CONVENTIONS = r"""
\begin{itemize}[leftmargin=*]
\item \textbf{Units.} Flux $J$ in mol\,m$^{-2}$\,s$^{-1}$; permeability $\Phi$ and apparent permeability
  $\Phi_\mathrm{app}$ in mol\,m$^{-1}$\,s$^{-1}$\,Pa$^{-1/2}$; pressures in Pa; diffusivities in m$^2$/s;
  activation energies of the sweep fits in kJ/mol.
\item \textbf{Two models.} Levels 1--5 are Model 1 (a Henry-law oxide, no surface step); the L6 family is
  Model 2 (surface dissociation and recombination, a Sieverts oxide). Compare across the two with care: a
  difference mixes the surface kinetics with the change of oxide law.
\item \textbf{Regimes} are named by the largest resistance fraction (surface, oxide, metal, or defect paths).
  The resistance fractions of a level sum to 1.
\item \textbf{Apparent permeability} $\Phi_\mathrm{app} = J\,L/\Delta\sqrt{P}$ ($\Delta P$ for L2a) is a material
  property only where thickness-proportional Sieverts transport controls; elsewhere it depends on the operating point.
\item \textbf{PRF}, the permeation reduction factor: at Level 5, bare-metal flux over the coated wall's flux; in the
  pinhole sweep, $J_p/J_i$, pinhole path over intact oxide.
\item \textbf{Fits.} The figures draw one fit per regime, over the window where that regime holds, named by its
  dominant path and limiting step. Two whole-sweep values are in the tables only, never drawn: the \emph{single fit}
  of the temperature sweep (one Arrhenius fit over every temperature) and the \emph{net} $n$ of the pressure sweep
  ($\Delta\ln J/\Delta\ln P$ between the two ends of the range, the average local exponent, not a fitted line).
\item \textbf{``--''} in a table is a value the study does not have.
\end{itemize}
"""

SUMMARY = {
    'title': 'Run summary',
    'docs': [
        ('docs/index.md', 'Where to start', 'the documentation, by question'),
        ('docs/getting-started.md', 'The model hierarchy', 'the levels L1--L6 the sections below follow'),
        ('docs/theory/two-models.md', 'The fork', 'why L1--L5 and the L6 family are two models, not six levels'),
        ('docs/config-file.md', 'One study, one file', r'the \texttt{study\_config.py} this report is about'),
        ('docs/config-file.md', '`CONDITIONS`', 'the operating point and the sweep ranges'),
        ('docs/how-to/switch-study.md', 'Then validate', r'\texttt{check\_against\_config}, whose findings are listed below'),
        ('docs/getting-started.md', 'Read these caveats before trusting a number', 'the caveats to read first'),
        ('docs/caveats.md', None, 'every known caveat, by level'),
        ('OUTPUTS.md', None, 'every returned quantity, function by function'),
        ('PARAMETERS.md', None, 'every input, with units and provenance'),
    ],
}

_LIMIT = ('(D)', 'A limit check: this level reduced to a lower one, or to a closed form.', '',
          [('points off the line', 'the limit is not recovered; the panel states the error')])

SECTIONS = {
    'levels': {
        'title': 'Levels 1--5: the transport hierarchy (Model 1)',
        'goal': r'Build the wall up one physical effect at a time --- bulk metal, an oxide layer, defects in the '
                r'oxide, grain boundaries and traps in the metal, then all together --- and see which effect '
                r'controls the flux at the operating point (\path{Proposal.ipynb}).',
        'docs': [
            ('docs/theory/level1-metal.md', 'What it models', "Sieverts' law and Fick diffusion"),
            ('docs/theory/level2-oxide.md', 'Level 2a — the perfect oxide', "a Henry-law oxide"),
            ('docs/theory/level2-oxide.md', 'Level 2b — oxide and metal in series', 'the interface-pressure solve'),
            ('docs/theory/level2-oxide.md', 'Reading the regime', 'the resistance ratio'),
            ('docs/theory/level3-defective-oxide.md', 'The parallel-path model', 'pinholes, cracks, grain boundaries'),
            ('docs/theory/level4-microstructure.md', 'Two effects, applied in sequence', 'grain-boundary enhancement, then trapping'),
            ('docs/theory/level5-full-system.md', 'The assembly', 'how Levels 3 and 4 combine'),
            ('docs/theory/level5-full-system.md', 'Reading the resistance fractions', 'what the fractions and regime mean'),
            ('docs/theory/level5-full-system.md', 'Values at the default operating point', 'the documented values, to compare with'),
            ('docs/caveats.md', 'Level 5 is metal-limited, so the oxide barely matters — but Level 5L6 is not',
             'why the two full systems disagree'),
        ],
        'columns': [
            ('level, quantity', r'the model level and the returned quantity (\texttt{OUTPUTS.md} names them all).'),
            ('value, units', 'at the operating point of this study.'),
        ],
        'figures': [
            {'key': 'L1', 'caption': 'Level 1: perfect metal (Sieverts + Fick).',
             'panels': [('(A)--(C)', 'Flux against pressure, permeability against temperature, and flux against '
                         'thickness, each with its fit and the expected value.', '',
                         [('slope 0.5 in $P$ and $-1$ in $L$', 'Sieverts dissolution and Fick diffusion, as the theory says')]),
                        ('(D)', 'The code against the closed form $J = \\Phi\\sqrt{P}/L$.', '',
                         [('points off the line', 'the implementation disagrees with the closed form')])]},
            {'key': 'L2a', 'caption': 'Level 2a: perfect oxide.',
             'panels': [('(A)--(C)', 'The oxide alone against pressure, temperature and thickness.', '',
                         [('slope 1 in $P$', "Henry's law: the oxide dissolves H$_2$ as molecules")]), _LIMIT]},
            {'key': 'L2b', 'caption': 'Level 2b: oxide and metal in series.',
             'panels': [('(A)--(C)', 'The coupled wall against pressure, temperature and oxide thickness, beside Level 1.', '',
                         [('a slope moving between 1 and 0.5', 'control passing between the oxide and the metal')]), _LIMIT]},
            {'key': 'L3', 'caption': 'Level 3: defective oxide.',
             'panels': [('(A)--(C)', 'The wall with oxide defects, against pressure, temperature and defect fraction.', '',
                         [('the curve meeting Level 1', 'the defects bypass the oxide')]), _LIMIT]},
            {'key': 'L4', 'caption': 'Level 4: defective metal (grain boundaries and traps).',
             'panels': [('(A)--(C)', 'The metal with grain-boundary enhancement and trapping, against its parameters.', '',
                         [('flux below Level 1', 'traps slow transport more than grain boundaries speed it')]), _LIMIT]},
            {'key': 'L5', 'caption': 'Level 5: the full system (defective oxide + defective metal).',
             'panels': [('(A)', 'Every level against pressure at the operating temperature.', '',
                         [('Level 5 on top of one lower level', 'the effect of that level controls the system there')]),
                        ('(B)', 'Arrhenius plot of Levels 1, 2b and 5, with activation energies.', '',
                         [('a bend', 'control passing from one step to another with temperature')]),
                        ('(C)', 'Oxide against metal effects: flux with no microstructure, grain boundaries only, traps '
                         'only, and both, at each oxide defect fraction.', '',
                         [('traps-only and combined bars equal', 'trapping dominates the microstructure effect')]),
                        ('(D)', r'Limit checks: Level 5 with $f \to 0$ against Level 2b, and $f \to 1$ against Level 1.', '',
                         [('points off the line', 'a limit is not recovered')])]},
        ],
    },
    'surface': {
        'title': 'Level 6: surface kinetics (Model 2)',
        'goal': r'Add the surface step --- H$_2$ dissociation and recombination with Langmuir coverage --- to each '
                r'level, and find where the surface, rather than the bulk, controls the flux '
                r'(\path{Surface_proposal.ipynb}).',
        'docs': [
            ('docs/theory/level6-surface-kinetics.md', 'Three fluxes, one unknown', 'the steady-state balance'),
            ('docs/theory/level6-surface-kinetics.md', 'Flux matching is the correctness check', 'what the limit checks test'),
            ('docs/theory/level6-surface-kinetics.md', 'The result that matters: the surface suppresses pressure, not flux directly',
             'how to read the interface pressure'),
            ('docs/theory/level6-surface-kinetics.md', 'Reading coverage', r'$\theta$ against the isotherm'),
            ('docs/theory/level6-surface-kinetics.md', 'The full system: L5L6', 'the full Model 2 system'),
            ('docs/theory/equilibrium-models.md', 'Why Langmuir and Sieverts do not conflict', 'the two equilibrium laws'),
            ('docs/caveats.md', '`rate_limiting` and `fraction_surface` do not measure what the surface cost',
             'why a small surface fraction does not mean a small surface effect'),
            ('docs/caveats.md', '`PRF` is `nan` from the Level 5L6 wrapper', 'the missing PRF of L5L6'),
            ('docs/caveats.md', 'L5-vs-L5L6 differences mix surface kinetics with a change of oxide law',
             'comparing L5 with L5L6'),
        ],
        'columns': [
            ('level, quantity', 'as in the previous section; $\\theta$ is the surface coverage, $\\eta_\\mathrm{surf}$ the '
             'surface efficiency of the L5L6 wrapper.'),
        ],
        'figures': [
            {'key': 'core_validation', 'caption': 'Core validation of one surface-kinetics level at the operating temperature.',
             'panels': [('(A) Flux against pressure', 'With the slopes of the surface-limited and diffusion-limited parts marked.', '',
                         [('slope 1 at low $P$, 0.5 at high $P$', 'surface-limited below the crossover, diffusion-limited above')]),
                        ('(B) What the level adds', 'As titled in the panel: the surface coverage beside the Langmuir isotherm '
                         'for the simplest levels; for the others, the effect of the level\\textquotesingle s own ingredient '
                         '--- the oxide against temperature, the defect fraction, or the metal microstructure.', '',
                         [(r'$\theta$ below $\theta_\mathrm{eq}$', 'fast diffusion empties the surface faster than it fills'),
                          (r'$D_\mathrm{eff}/D_m$ below 1', 'trapping slows transport more than grain boundaries speed it'),
                          ('a curve that stays flat across the defect fraction', 'the defects do not open a faster path at this point')]),
                        ('(C) Rate-limiting step', 'The resistance fraction of each step against pressure.', '',
                         [('the curves crossing at 50\\,\\%', 'the crossover pressure between surface and bulk control')]),
                        ('(D) Limit check', 'The level reduced to a simpler one: fast dissociation against no surface step, '
                         'or the microstructure switched off against the level without it.', '',
                         [('points off the line, or a warning box', 'the limit is not reached: read the box before trusting the panel')])]},
        ],
    },
    'temperature': {
        'title': 'Effect of temperature',
        'goal': r'Sweep the temperature over \texttt{T\_range} at the operating pressure, fit $J = A\,e^{-E_a/RT}$ '
                r'over each regime\textquotesingle s window, and see how the controlling step changes with temperature. '
                r'The figures show the regime fits; the single fit over every temperature is reported in the tables '
                r'only (\path{effect-of-temperature.ipynb}).',
        'docs': [
            ('docs/how-to/parameter-sweeps.md', 'The temperature sweep: regime fits and the single fit',
             'what the figures draw, and what only the table gives'),
            ('docs/how-to/parameter-sweeps.md', 'Reading a table', 'every column of the fits tables'),
            ('docs/theory/level1-metal.md', 'Temperature dependence', 'the activation energy of the perfect metal'),
            ('docs/theory/permeability.md', 'Starting point: which model are you in?', 'what an apparent permeability means in each model'),
            ('docs/theory/two-models.md', 'What may be compared with what', 'which levels may be compared'),
            ('docs/caveats.md', '`permeability` is apparent and needs its operating point', r'why $\Phi_\mathrm{app}$ is not a material property'),
            ('docs/caveats.md', 'Both materials are extrapolated at the default operating point', 'temperatures outside the data'),
        ],
        'figures': [
            {'key': 'main', 'caption': 'Temperature sweep at the operating pressure: regime-resolved Arrhenius fits.',
             'panels': [('Flux (top) and $\\Phi_\\mathrm{app}$ (bottom)', r'Against $1000/T$, coloured by regime, with one '
                         'Arrhenius fit per regime window.',
                         'Marker and colour name the dominant path and the limiting step; each fit, and the single fit '
                         'over every temperature, is in the table.',
                         [('a steeper segment', 'a step with a higher activation energy controls there'),
                          ('a change of marker', 'the dominant path changes, for instance from intact oxide to pinholes'),
                          ('a single-fit $E_a$ far from the regime values', 'it averages over regimes; quote the regime fits')])]},
            {'key': '3P', 'caption': 'The temperature sweep at the low, middle and high end of \\texttt{P\\_range}.',
             'panels': [('Each panel', 'One curve per pressure; colour and marker are the regime, line style the pressure.', '',
                         [('parallel curves', 'the same step controls at every pressure'),
                          ('regime changes that move with pressure', 'the crossover temperature depends on pressure')])]},
        ],
    },
    'pressure': {
        'title': 'Effect of pressure',
        'goal': r'Sweep the upstream pressure over \texttt{P\_range} at the operating temperature and fit $J = A\,P^{\,n}$ '
                r'over each regime\textquotesingle s window: the slope $n$ names the controlling step, $n \approx 1$ for a '
                r'Henry oxide (Model 1) or the surface rate (Model 2), $n = 0.5$ for Sieverts transport. The tables also '
                r'give the \emph{net} $n$ of the whole trend, $\Delta\ln J/\Delta\ln P$ between the two ends of the range '
                r'--- the average of the local exponent, not a fitted line, so it is not drawn '
                r'(\path{effect-of-pressure.ipynb}).',
        'docs': [
            ('docs/how-to/parameter-sweeps.md', 'The pressure exponent: regime n and net n',
             'the two exponents of the tables, defined'),
            ('docs/how-to/parameter-sweeps.md', 'Why not one least-squares line over the whole range',
             'why the net $n$ replaces a single fit'),
            ('docs/how-to/parameter-sweeps.md', 'Reading the two together', 'what a net $n$ between the regime values means'),
            ('docs/theory/level1-metal.md', 'The √P signature', '$n = 0.5$ for Sieverts transport'),
            ('docs/theory/level2-oxide.md', "Henry's law, not Sieverts", '$n = 1$ for a Henry oxide'),
            ('docs/theory/level6-surface-kinetics.md', 'The result that matters: the surface suppresses pressure, not flux directly',
             'the surface in a pressure sweep'),
            ('docs/theory/permeability.md', 'Tier 4 — the wall, and it genuinely drifts with pressure',
             r'why $\Phi_\mathrm{app}$ of the wall changes with pressure'),
            ('docs/caveats.md', '`permeability` is apparent and needs its operating point', r'quoting $\Phi_\mathrm{app}$'),
        ],
        'figures': [
            {'key': 'main', 'caption': 'Pressure sweep at the operating temperature: regime-resolved power-law fits.',
             'panels': [('Flux (top)', '$J$ against $P_\\mathrm{up}$ on log--log axes, coloured by regime, with one fit per '
                         'regime and the labelled reference slopes $n = 1$ (long dashes) and $n = 0.5$ (short dashes).',
                         'Compare each regime\\textquotesingle s line with the references; the net $n$ of the whole range is '
                         'in the table.',
                         [('$n \\approx 1$', 'a Henry oxide (Model 1) or the surface rate (Model 2) controls'),
                          ('$n \\approx 0.5$', 'Sieverts transport controls'),
                          ('$n$ in between', 'mixed control, or a regime change inside the window')]),
                        ('Apparent permeability (bottom)', r'$\Phi_\mathrm{app}$ from the same flux; its slope is $n - 0.5$ ($n - 1$ for L2a).', '',
                         [('a flat plateau', 'Sieverts transport controls: $\\Phi_\\mathrm{app}$ is a true permeability there'),
                          ('a rise with $\\sqrt{P}$', 'it is not a material property in that window')])]},
            {'key': '3T', 'caption': 'The pressure sweep at the low, middle and high end of \\texttt{T\\_range}.',
             'panels': [('Each panel', 'One curve per temperature; colour and marker are the regime, line style the '
                         'temperature; the reference slopes are drawn through the middle temperature.', '',
                         [('the slope changing with temperature', 'the controlling step depends on temperature')])]},
        ],
    },
    'thickness': {
        'title': 'Effect of thickness',
        'goal': r'Sweep the oxide thickness, the metal thickness, and the whole wall scaled at a fixed ratio, at the '
                r'operating point, and fit $J \propto L^{\,m}$: $m = -1$ when the swept thickness controls the flux, '
                r'$m \approx 0$ when something else does. The scaled wall separates bulk from surface: surface '
                r'kinetics do not depend on thickness, so they pull $m$ toward 0 (\path{effect-of-thickness.ipynb}).',
        'docs': [
            ('docs/how-to/parameter-sweeps.md', 'The thickness sweep: one figure per swept layer',
             'the three sweeps, and which levels get which'),
            ('docs/theory/level1-metal.md', 'The resistance analogy, and where it breaks', 'resistances in series'),
            ('docs/theory/level2-oxide.md', 'The limiting cases', 'oxide-limited and metal-limited walls'),
            ('docs/theory/permeability.md', 'Why the sorption law decides whether permeability is intrinsic',
             r'when $\Phi_\mathrm{app}$ is a material property'),
            ('docs/theory/permeability.md', 'The surface is a driving force, not a resistance', 'why surface kinetics pull $m$ toward 0'),
        ],
        'figures': [
            {'key': 'main', 'caption': 'A thickness sweep at the operating point: one figure per swept layer.',
             'panels': [('Left and right', '$J$ and $\\Phi_\\mathrm{app}$ against the swept thickness: the oxide (metal '
                         'fixed), the metal (oxide fixed), or the whole wall scaled at a fixed ratio. A level has only the '
                         'sweeps of the layers it has; a one-layer level (L1, L2a, L4, L1L6, L2aL6) has just its own '
                         'layer, since scaling its whole wall would repeat that sweep.', 'The fits table gives $m$ for each regime.',
                         [('$m = -1$', 'the swept thickness controls the flux'),
                          ('$m \\approx 0$', 'something else controls: the other layer, or the surface'),
                          ('the whole wall with $m$ above $-1$ and $\\Phi_\\mathrm{app}$ growing with $L$', 'surface kinetics take part')])]},
            {'key': '3P', 'caption': 'A thickness sweep at the low, middle and high end of \\texttt{P\\_range}.',
             'panels': [('Left and right', 'As above, one curve per pressure.', '', [])]},
            {'key': '3T', 'caption': 'A thickness sweep at the low, middle and high end of \\texttt{T\\_range}.',
             'panels': [('Left and right', 'As above, one curve per temperature.', '', [])]},
        ],
    },
    'pinholes': {
        'title': 'Effect of pinholes',
        'goal': r'Sweep the pinhole area fraction $f$ with cracks and oxide grain boundaries at zero. Intact oxide and '
                r'bare metal are two paths side by side, $J = (1-f)\,J_i + f\,J_p$, and one number sets the curve: '
                r'$\mathrm{PRF} = J_p/J_i$. Pinholes carry half the flux at $f_{1/2} = 1/(1+\mathrm{PRF})$, so the better '
                r'the oxide, the fewer pinholes it takes to undo it (\path{effect-of-pinholes.ipynb}).',
        'docs': [
            ('docs/how-to/parameter-sweeps.md', 'The pinhole sweep: PRF and the half-flux fraction',
             'PRF, $f_{1/2}$ and the paths table'),
            ('docs/theory/level3-defective-oxide.md', 'When do defects actually matter?', 'the conditions for a defect to matter'),
            ('docs/theory/level3-defective-oxide.md', 'Permeation Reduction Factor', 'PRF'),
            ('docs/theory/permeability.md', 'The coupled wall: pinholes are conditional', 'pinholes in the coupled wall'),
            ('docs/theory/permeability.md', 'Contrast, not area fraction, is the lever', r'why PRF, not $f$, sets the curve'),
            ('docs/theory/permeability.md', "The pinhole path carries the metal's surface kinetics", 'the pinhole path in Model 2'),
            ('docs/theory/equilibrium-models.md', 'The pinhole special case: a third surface to worry about', 'the bare-metal surface'),
        ],
        'figures': [
            {'key': 'main', 'caption': 'Pinhole sweep at the operating point.',
             'panels': [('Flux and apparent permeability', 'Against the pinhole area fraction $f$; grey: no pinholes.',
                         'Colour and marker are the regime; lines are power-law fits $J \\propto f^{\\,n}$.',
                         [('flat until a knee', r'the intact oxide carries the flux until $f$ nears $f_{1/2}$'),
                          ('a rise after the knee', 'pinholes carry most of the flux: the oxide is bypassed'),
                          ('a fall', r'the pinhole path passes less than intact oxide ($\mathrm{PRF} < 1$), for instance '
                           'when dissociation on the bare metal is slow')])]},
            {'key': '3P', 'caption': 'The pinhole sweep at the low, middle and high end of \\texttt{P\\_range}.',
             'panels': [('Each panel', 'One curve per pressure; colour and marker are the regime, line style the condition.', '', [])]},
            {'key': '3T', 'caption': 'The pinhole sweep at the low, middle and high end of \\texttt{T\\_range}.',
             'panels': [('Each panel', 'One curve per temperature.', '', [])]},
        ],
    },
}

_SA_DOCS = [
    ('docs/how-to/sensitivity-analysis.md', 'What this answers, and why it is stratified', 'why the analysis is split by regime'),
    ('docs/how-to/sensitivity-analysis.md', 'The two estimators', r'Borgonovo $\delta$ and PAWN'),
    ('docs/how-to/sensitivity-analysis.md', 'The dummy parameter sets the significance threshold', r'what $\delta$/floor $> 1$ means'),
    ('docs/how-to/sensitivity-analysis.md', 'Below the floor does not mean "no physical effect"', 'reading a parameter below the floor'),
    ('docs/how-to/sensitivity-analysis.md', 'Sizing the draws', 'the cluster sizes'),
    ('docs/config-file.md', '`REGIME_PRESETS` and `REGIME_PRESETS_L5`', "the presets that aim each regime's draws"),
    ('docs/caveats.md', 'Presets suppress the very parameters they pin', 'a pinned parameter cannot rank'),
    ('docs/caveats.md', 'The default `n_dummy=3` is too small and inflates every margin', 'the floor itself is uncertain'),
    ('docs/caveats.md', "The sweep range extrapolates past the oxide's validated window", 'the ranges of the scan'),
]
_SA_COLUMNS = [
    ('regime', 'the cluster: draws whose largest resistance is in that part of the wall.'),
    (r'$\delta$/floor at $T$', r'Borgonovo $\delta$ of the parameter for the flux, over the dummy parameter\textquotesingle s '
     r'$\delta$; above 1 the parameter is distinguishable from noise at that temperature.'),
    ('best', 'the largest of those ratios; the table lists parameters with best $\\ge 1$, highest first.'),
    ('draws', 'how many scan draws fell in the regime at that temperature.'),
]
_SA_FIGURES = [
    {'key': 'regime_exploration', 'caption': 'The outputs of each regime cluster, at one temperature.',
     'panels': [('Panels', 'The distribution of each output (flux, permeability, coverage) within each regime.', '',
                 [('well separated distributions', 'the regimes are distinct states, so splitting by regime is meaningful')])]},
    {'key': 'givendata', 'caption': 'Given-data sensitivity of the flux within one regime.',
     'panels': [('Left', r'Borgonovo $\delta$ per parameter, with its spread.', '',
                 [('a bar no longer than the dummy', 'indistinguishable from noise at this sample size, which is not the same as no effect')]),
                ('Right', 'PAWN median per parameter.', '',
                 [('the same order as the left panel', 'the ranking does not depend on the estimator')])]},
    {'key': 'comparison_heatmap', 'caption': r'Cross-regime sensitivity ($\delta$) of the flux, at one temperature.',
     'panels': [('Cells', r'Colour and the large number: $\delta$ as a fraction of its column maximum; corner: the raw $\delta$.', '',
                 [('a parameter dark in one column only', 'it matters only in that regime')])]},
    {'key': 'regime_geometry', 'caption': "Where the clusters lie against one regime's two strongest drivers.",
     'panels': [('Axes', r'Each cluster\textquotesingle s draws against the top two $\delta$ drivers; boxes outline the clusters.', '',
                 [('overlapping clusters', 'those two parameters alone do not decide the regime')])]},
]
SECTIONS['sensitivity_L5'] = {
    'title': 'Sensitivity by regime: Level 5',
    'goal': r'Find which inputs control the Level 5 flux within each regime (oxide, metal, defect): a Latin-hypercube '
            r'scan aimed at each regime, split into clusters, and given-data sensitivity within each cluster '
            r'(\path{sensitivity_regime_L5.ipynb}).',
    'docs': _SA_DOCS, 'columns': _SA_COLUMNS, 'figures': _SA_FIGURES}
SECTIONS['sensitivity_L5L6'] = {
    'title': 'Sensitivity by regime: Level 5 + surface kinetics',
    'goal': r'As for Level 5, for the full Model 2 system, whose regimes are surface, oxide and metal '
            r'(\path{sensitivity_regime_L5L6.ipynb}).',
    'docs': _SA_DOCS, 'columns': _SA_COLUMNS, 'figures': _SA_FIGURES}
SECTIONS['regimes'] = {
    'title': 'Regime clusters in parallel coordinates',
    'goal': r'Interactive parallel-coordinates views of the regime clusters of the sensitivity scans '
            r'(\path{regime_parallel_coords.ipynb}, \path{regime_parallel_coords_L5.ipynb}). They are HTML pages, '
            r'to open in a browser; a PDF cannot hold them, so this section lists them.',
    'docs': [
        ('docs/how-to/sensitivity-analysis.md', 'Notebooks and outputs', 'which notebook writes what'),
        ('docs/caveats.md', 'Run the SA notebooks before the parallel-coordinates ones', 'the order the notebooks must run in'),
    ],
}
