r"""One PDF per study, built from the study folder's notebook outputs.

The report reads what the notebooks wrote beside them --- ``<Notebook>-results/``
(figures, and the fit tables of the effect-of notebooks) and ``sa_results*/`` (the
sensitivity tables) --- plus ``report/data/``, which :mod:`calculations.report_data`
writes first: the model at the operating point, level by level, and the facts of
the study. Every notebook with outputs gets a section: its goal, where in ``docs/``
to read more, what was computed, key numbers, its tables, its figures, and how to
read them. The fixed text is in :mod:`calculations.report_text`; the LaTeX blocks are
in :mod:`calculations.report_latex`.

    python -m calculations.report Guo_etal_2025                 # -> Guo_etal_2025/report/Guo_etal_2025_report.pdf
    python -m calculations.report Guo_etal_2025 --tex-only      # no LaTeX here: keep the .tex
    python -m calculations.report Guo_etal_2025 --no-data       # reuse report/data/ as it is
    python -m calculations.report Guo_etal_2025 --run-notebooks # run every notebook first, in order
    python -m calculations.report Guo_etal_2025 --run-notebooks Proposal,effect-of-pinholes

Importing this module binds no study; the data step runs in its own process with
``PERMEATION_STUDY`` set, because the model resolves one study per process.
"""
import argparse
import csv
import datetime
import glob
import hashlib
import json
import os
import re
import subprocess
import sys
from collections import Counter

from calculations import report_latex as L
from calculations import report_text as RT

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SECTION_ORDER = ['levels', 'surface', 'temperature', 'pressure', 'thickness', 'pinholes',
                 'sensitivity_L5', 'sensitivity_L5L6', 'regimes']

# the folders of each section; a section is written when one of them has files
SECTION_FOLDERS = {
    'levels': ['Proposal-results'],
    'surface': ['Surface_proposal-results'],
    'temperature': ['effect-of-temperature-results'],
    'pressure': ['effect-of-pressure-results'],
    'thickness': ['effect-of-thickness-results'],
    'pinholes': ['effect-of-pinholes-results'],
    'sensitivity_L5': ['sensitivity_regime_L5-results', 'sa_results_L5'],
    'sensitivity_L5L6': ['sensitivity_regime_L5L6-results', 'sa_results'],
    'regimes': ['regime_parallel_coords-results', 'regime_parallel_coords_L5-results'],
}
# the effect-of sweeps: section -> (results folder, the word in the file names, figure variants)
SWEEPS = {
    'temperature': ('effect-of-temperature-results', 'arrhenius', ('3P',)),
    'pressure': ('effect-of-pressure-results', 'pressure', ('3T',)),
    'thickness': ('effect-of-thickness-results', 'thickness', ('3P', '3T')),
    'pinholes': ('effect-of-pinholes-results', 'pinholes', ('3P', '3T')),
}
LEVEL_NAMES = dict(RT.LEVELS)
LEVEL_RANK = {lv: i for i, (lv, _) in enumerate(RT.LEVELS)}
PROPOSAL_STEMS = {'L1': 'Level1_Perfect_Metal', 'L2a': 'Level2a_Perfect_Oxide', 'L2b': 'Level2b_Oxide_Metal',
                  'L3': 'Level3_Defective_Oxide', 'L4': 'Level4_Defective_Metal', 'L5': 'Level5_Full_System'}
NEAR = 0.05                         # |slope - ideal| within this reads as the ideal value


# ── values ───────────────────────────────────────────────────────────────────
def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if v == v else None


def level_of(stem):
    """'L5L6_arrhenius_3P' -> 'L5L6'; 'Level2b_Oxide_Metal' -> 'L2b'; else None."""
    m = re.match(r'Level(\d[ab]?)_', stem)
    if m:
        return 'L' + m.group(1)
    m = re.match(r'(L\d[ab]?(?:L\d[ab]?)*)_', stem)
    return m.group(1) if m else None


def level_sort(levels):
    return sorted(set(levels), key=lambda lv: (LEVEL_RANK.get(lv, 99), lv))


def level_name(lv):
    return LEVEL_NAMES.get(lv, lv)


def model_of(lv):
    return 'Model 2' if lv.endswith('L6') else 'Model 1'


def git(*args):
    try:
        return subprocess.run(['git', '-C', REPO_ROOT, *args], capture_output=True, text=True,
                              timeout=20).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ''


def slope_word(v, ideals):
    """The ideal a fitted slope sits on, from [(value, words)], or 'between'."""
    if v is None:
        return '--'
    for ideal, words in ideals:
        if abs(v - ideal) <= NEAR:
            return words
    return 'between the ideal values'


# ── the study folder ─────────────────────────────────────────────────────────
class StudyReport(L.LatexReport):
    """The report of one study folder."""

    def __init__(self, study_dir, log=print, web_base='auto'):
        super().__init__(os.path.abspath(study_dir), docs_root=REPO_ROOT,
                         web_base=L.git_web_base(REPO_ROOT) if web_base == 'auto' else web_base, log=log)
        self.name = os.path.basename(self.out.rstrip('/'))
        self.config = os.path.join(self.out, 'study_config.py')
        try:
            with open(os.path.join(self.dir, 'data', 'study_facts.json')) as fh:
                self.facts = json.load(fh)
        except (OSError, ValueError):
            self.facts = {}
        self.op = self.rows('report/data/operating_point.csv', mark=False)
        self.built = []

    def rows(self, rel, mark=True):
        path = os.path.join(self.out, rel)
        if not os.path.isfile(path):
            return []
        if mark:
            self.used.add(rel)
        with open(path, newline='', encoding='utf-8', errors='replace') as fh:
            return list(csv.DictReader(fh))

    def files(self, folder, pattern='*'):
        return sorted(os.path.relpath(p, self.out) for p in glob.glob(os.path.join(self.out, folder, pattern))
                      if os.path.isfile(p))

    def figures(self, folder):
        """{stem: file} of a results folder, the PDF when there is one (vector), else the PNG."""
        out = {}
        for f in self.files(folder):
            stem, ext = os.path.splitext(os.path.basename(f))
            if ext.lower() == '.pdf' or (ext.lower() == '.png' and stem not in out):
                out[stem] = f
        return out

    def outputs(self, name):
        return [f for folder in SECTION_FOLDERS[name] for f in self.files(folder)]

    def op_value(self, level, quantity):
        for r in self.op:
            if r['level'] == level and r['quantity'] == quantity:
                return r['value']
        return None


def op_num(rep, level, quantity):
    return num(rep.op_value(level, quantity))


def op_text(rep, level, quantity):
    v = rep.op_value(level, quantity)
    return '--' if v in (None, '') else L.tex_num(v)


# ── title page and run summary ───────────────────────────────────────────────
def title_page(rep):
    f = rep.facts
    commit = git('rev-parse', '--short', 'HEAD')
    dirty = bool(git('status', '--porcelain', '--', 'calculations'))
    sha = f.get('study_sha256')
    if not sha and os.path.isfile(rep.config):
        with open(rep.config, 'rb') as fh:
            sha = hashlib.sha256(fh.read()).hexdigest()
    cond = f.get('conditions') or {}
    point = ', '.join(s for s in (
        f"$T$ = {cond['T_operating']:g} K" if 'T_operating' in cond else '',
        f"$P_\\mathrm{{up}}$ = {L.tex_num(cond['P_upstream'])} Pa" if 'P_upstream' in cond else '',
        f"$L_\\mathrm{{metal}}$ = {L.tex_num(cond['L_metal'])} m" if 'L_metal' in cond else '',
        f"$L_\\mathrm{{oxide}}$ = {L.tex_num(cond['L_oxide'])} m" if 'L_oxide' in cond else '') if s)
    facts = [('study', L.tex_escape(f.get('study_name') or rep.name)), ('folder', L.tex_path(rep.out)),
             ('input file', L.tex_file('study_config.py') + (f', sha256 \\texttt{{{sha[:16]}}}\\ldots' if sha else ''))]
    if f.get('metal_name'):
        facts.append(('metal, oxide', L.tex_escape(f['metal_name']) + ', ' + L.tex_escape(f.get('oxide_name') or '--')))
    if point:
        facts.append(('operating point', point))
    facts += [('operating point computed', L.tex_escape((f.get('computed_at') or 'not computed').replace('T', ' '))),
              ('report built', datetime.datetime.now().strftime('%Y-%m-%d %H:%M')),
              ('code', ('commit ' + L.tex_escape(commit) + (r', \textbf{uncommitted changes in calculations/}' if dirty else ''))
               if commit else 'not a git checkout'),
              ('documentation', 'linked on GitHub at that commit' if rep.web_base else r'\texttt{docs/} in the repository')]
    title = L.tex_escape(f.get('study_name') or rep.name) + ': hydrogen permeation report'
    footer = (r'Built from the notebook outputs of the study folder and the operating-point table of '
              r'\texttt{calculations.report\_data}; no notebook was run to make it. Rebuild with '
              r'\texttt{python -m calculations.report ' + L.tex_escape(rep.name) + r'}.')
    return title, L.title_page(title, r'study folder \texttt{' + L.tex_escape(rep.name) + '}', facts, footer)


def results_age(rep):
    """Per results folder: when it was last written, and whether that was before the last change
    to study_config.py or the last commit to calculations/."""
    cfg_t = os.path.getmtime(rep.config) if os.path.isfile(rep.config) else None
    code_t = num(git('log', '-1', '--format=%ct', '--', 'calculations'))
    rows = []
    for folder in sorted({x for v in SECTION_FOLDERS.values() for x in v}):
        fs = rep.files(folder)
        if not fs:
            continue
        t = max(os.path.getmtime(os.path.join(rep.out, x)) for x in fs)
        rows.append({'folder': folder + '/', 'n': len(fs),
                     'made': datetime.datetime.fromtimestamp(t).strftime('%Y-%m-%d %H:%M'),
                     'cfg': None if cfg_t is None else t < cfg_t, 'code': None if code_t is None else t < code_t})
    return rows


def run_summary(rep):
    f = rep.facts
    out = [r'\section{Run summary}', r'\label{sec:summary}',
           r'\paragraph{Where to read more.}', rep.pointers_block('Run summary', RT.SUMMARY['docs'])]

    # the operating point
    cond = f.get('conditions') or {}
    if cond:
        keys = ['T_operating', 'P_upstream', 'P_downstream', 'L_metal', 'L_oxide', 'T_range', 'P_range',
                'L_metal_range', 'L_oxide_range', 'defect_fraction_range']
        rows = [{'k': k, 'v': ', '.join(L.tex_num(x) for x in cond[k]) if isinstance(cond[k], list) else L.tex_num(cond[k])}
                for k in keys if k in cond]
        out += [r'\subsection{The operating point}',
                rep.table(rows, [('k', 'entry of \\texttt{CONDITIONS}', 'l', L.tex_file), ('v', 'value', 'l', lambda s: s)],
                          r'The operating point and sweep ranges of \texttt{study\_config.py} (SI units: K, Pa, m).')]
        props = f.get('properties_at_T') or {}
        if props:
            out.append(r'At $T_\mathrm{operating}$: ' + ', '.join(
                f'{L.tex_escape(k)} = {L.tex_num(v)}' for k, v in props.items() if v is not None)
                + f' (metal {L.tex_escape(f.get("metal_name"))}, oxide {L.tex_escape(f.get("oxide_name"))}).')

    # the model at the operating point
    if rep.op:
        rep.used.add('report/data/operating_point.csv')
        out += [r'\subsection{The model at the operating point}',
                'Each level evaluated at the operating point with the functions the smoke test and the sensitivity '
                'wrappers use (the source column names them). The sections below quote these values.',
                rep.table(rep.op, [('level', 'level', 'l', L.tex_escape), ('quantity', 'quantity', 'l', L.tex_file),
                                   ('value', 'value', 'r', None), ('units', 'units', 'l', L.tex_escape),
                                   ('source', 'source', 'l', L.tex_file)],
                          'The model at the operating point.', 'report/data/operating_point.csv', group_by='level')]
    else:
        out += [r'\subsection{The model at the operating point}',
                r'\textbf{Not computed}: \texttt{report/data/operating\_point.csv} is missing. Rebuild without '
                r'\texttt{--no-data}, or see the build log for why \texttt{calculations.report\_data} failed.']

    # caveats
    items = [r'\texttt{check\_against\_config}: ' + L.tex_md(p) for p in f.get('config_problems') or []]
    items += [r'The model while computing the operating point: ' + L.tex_md(w) for w in f.get('model_warnings') or []]
    items += [r'\texttt{report\_data} failed for ' + L.tex_md(e) for e in f.get('errors') or []]
    out += [r'\subsection{Read these before quoting any number}',
            (L.enumerated(items) if items else 'The configuration check found nothing, and the model raised no warnings.')
            + r'\par See also \texttt{docs/caveats.md}, linked above.']

    # results and their age
    rows = results_age(rep)
    out.append(r'\subsection{Notebook outputs, and the inputs they were made from}')
    if rows:
        dirty = bool(git('status', '--porcelain', '--', 'calculations'))
        out.append('A results folder is written when its notebook runs. One older than \\texttt{study\\_config.py}, '
                   'or than the last commit to \\texttt{calculations/}, may not reflect them: re-run its notebook.'
                   + (r' \textbf{\texttt{calculations/} has uncommitted changes}, which no saved output can be known to reflect.'
                      if dirty else ''))
        out.append(rep.table(rows, [('folder', 'folder', 'l', L.tex_file), ('n', 'files', 'r', None),
                                    ('made', 'last written', 'l', L.tex_escape),
                                    ('cfg', r'older than\\study\_config.py', 'l', L.tex_flag),
                                    ('code', r'older than the last\\commit to calculations/', 'l', L.tex_flag)],
                             'Results folders of this study.'))
    else:
        out.append('No notebook of this study has written results yet.')

    rows = [{'name': RT.SECTIONS[s]['title'], 'n': len(rep.outputs(s)),
             'ref': (r'\ref{sec:' + s + '}') if s in rep.built else '--'} for s in SECTION_ORDER]
    out += [r'\subsection{Sections}',
            rep.table(rows, [('name', 'section', 'l', lambda s: s), ('n', r'files\\found', 'r', None),
                             ('ref', 'section', 'l', lambda s: s)],
                      'The sections, one per notebook, in the order of the hierarchy. A section is shown for each '
                      'notebook with outputs.')]
    out += [r'\subsection{Conventions used in this report}', RT.CONVENTIONS]
    return '\n\n'.join(out)


# ── Levels 1-5 and Level 6 ───────────────────────────────────────────────────
def op_rows(rep, levels):
    return [r for r in rep.op if r['level'] in levels]


def op_table(rep, levels, caption):
    rows = op_rows(rep, levels)
    if not rows:
        return ''
    return rep.table(rows, [('level', 'level', 'l', L.tex_escape), ('quantity', 'quantity', 'l', L.tex_file),
                            ('value', 'value', 'r', None), ('units', 'units', 'l', L.tex_escape)],
                     caption, 'report/data/operating_point.csv', group_by='level')


def section_levels(rep):
    figs = rep.figures('Proposal-results')
    keys = []
    l1, l2b, l5 = op_num(rep, 'L1', 'flux'), op_num(rep, 'L2b', 'flux'), op_num(rep, 'L5', 'flux')
    if l1 is not None:
        keys.append(f'Level 1, the bare metal: $J$ = {L.tex_num(l1)} mol\\,m$^{{-2}}$\\,s$^{{-1}}$.')
    if l2b is not None:
        rr = op_num(rep, 'L2b', 'resistance_ratio')
        keys.append(f'Level 2b, oxide on metal: $J$ = {L.tex_num(l2b)}'
                    + (f', {L.tex_num(l1 / l2b)} times below Level 1' if l1 and l2b else '')
                    + f'; interface pressure {op_text(rep, "L2b", "P_interface")} Pa; oxide-to-metal resistance ratio '
                    + (L.tex_num(rr) + (' (the metal controls)' if rr is not None and rr < 1 else ' (the oxide controls)')
                       if rr is not None else '--') + '.')
    if l5 is not None:
        keys.append(f'Level 5, the full system: $J$ = {L.tex_num(l5)}, PRF = {op_text(rep, "L5", "PRF")}; regime '
                    f'{L.tex_file(rep.op_value("L5", "regime") or "--")} '
                    f'({L.tex_file(rep.op_value("L5", "regime_hierarchy") or "--")}), dominant path '
                    f'{L.tex_file(rep.op_value("L5", "dominant_path") or "--")}; resistance split oxide '
                    f'{op_text(rep, "L5", "frac_oxide")}, metal {op_text(rep, "L5", "frac_metal")}, defects '
                    f'{op_text(rep, "L5", "frac_defect")}.')
        dm = op_num(rep, 'L5', 'D_modification')
        if dm is not None:
            keys.append(f'The microstructure (Level 4, inside Level 5) changes the diffusivity by a factor '
                        f'{L.tex_num(dm)}: $D_\\mathrm{{eff}}$ = {op_text(rep, "L5", "D_eff")} m$^2$/s.')
    parts = [('What was computed', 'The levels one at a time, each against pressure, temperature and its own parameter, '
              'with a limit check against a lower level; and, for this report, each level at the operating point.'),
             ('Key numbers', L.bullets(keys) or 'The operating-point table was not computed.'),
             ('At the operating point', op_table(rep, ('L1', 'L2b', 'L5'), 'Levels 1, 2b and 5 at the operating point.'))]
    for lv, stem in PROPOSAL_STEMS.items():
        if stem in figs:
            parts.append(('figure', lv, figs[stem]))
    for stem in sorted(set(figs) - set(PROPOSAL_STEMS.values())):
        parts.append(('figure', 'L5', figs[stem], f'({L.tex_escape(stem)}.)'))
    return rep.section('levels', RT.SECTIONS['levels'], parts)


def section_surface(rep):
    figs = rep.figures('Surface_proposal-results')
    keys = []
    l1, l1l6 = op_num(rep, 'L1', 'flux'), op_num(rep, 'L1L6', 'flux')
    if l1l6 is not None:
        keys.append(f'Level 1 + surface: $J$ = {L.tex_num(l1l6)}' + (f', {L.tex_num(l1 / l1l6)} times below the bare metal'
                                                                    if l1 and l1l6 else '')
                    + f'; coverage $\\theta$ = {op_text(rep, "L1L6", "theta")}; interface pressure '
                    f'{op_text(rep, "L1L6", "P_interface")} Pa; rate-limiting '
                    f'{L.tex_file(rep.op_value("L1L6", "rate_limiting") or "--")}.')
    l5l6, l5 = op_num(rep, 'L5L6', 'flux'), op_num(rep, 'L5', 'flux')
    if l5l6 is not None:
        keys.append(f'Level 5 + surface, the full Model 2 system: $J$ = {L.tex_num(l5l6)}'
                    + (f' ({L.tex_num(l5 / l5l6)} times below Level 5)' if l5 and l5l6 else '')
                    + f'; regime {L.tex_file(rep.op_value("L5L6", "regime") or "--")}, system rate-limiting '
                    f'{L.tex_file(rep.op_value("L5L6", "system_rate_limiting") or "--")}; resistance split surface '
                    f'{op_text(rep, "L5L6", "frac_surface")}, oxide {op_text(rep, "L5L6", "frac_oxide")}, metal '
                    f'{op_text(rep, "L5L6", "frac_metal")}; $\\theta$ = {op_text(rep, "L5L6", "theta")}.')
        if rep.op_value('L5', 'regime') and rep.op_value('L5', 'regime') != rep.op_value('L5L6', 'regime'):
            keys.append('Level 5 and Level 5 + surface name different regimes. The difference mixes the surface step '
                        'with the change of oxide law between the two models (the caveat linked above).')
        if op_num(rep, 'L5L6', 'PRF') is None:
            keys.append(r'PRF of Level 5 + surface is not defined by its wrapper (\texttt{nan}; the caveat linked above).')
    parts = [('What was computed', 'One core-validation figure per level with the surface step: flux, coverage (or the '
              'microstructure effect) and the rate-limiting step against pressure, and a limit check against a simpler '
              'level; and Level 1 + surface and Level 5 + surface at the operating point.'),
             ('Key numbers', L.bullets(keys) or 'The operating-point table was not computed.'),
             ('At the operating point', op_table(rep, ('L1L6', 'L5L6'), 'The surface-kinetics levels at the operating point.'))]
    by_level = {}
    for stem, f in figs.items():
        by_level.setdefault(level_of(stem) or 'other', []).append((stem, f))
    for lv in level_sort(by_level):
        parts.append(('+' + level_name(lv), f'{level_name(lv)} ({model_of(lv)}).'))
        for stem, f in sorted(by_level[lv], key=lambda x: ('core_validation' not in x[0], x[0])):
            # one figure per level is left (the notebook's FIGURE 1, core validation), whatever its file name
            parts.append(('figure', 'core_validation', f, f'({L.tex_escape(stem)}.)', 0.55))
    return rep.section('surface', RT.SECTIONS['surface'], parts)


# ── the effect-of sweeps ─────────────────────────────────────────────────────
def value_column(header):
    """The fitted slope or activation-energy column: the one after 'points'."""
    if 'points' in header and header.index('points') + 1 < len(header):
        return header[header.index('points') + 1]
    return None


def fits_table(rep, rel, caption):
    rows = rep.rows(rel)
    if not rows:
        return ''
    cols = []
    for h in rows[0].keys():
        if h == 'fit':
            cols.append((h, 'fit', 'P4.4', L.tex_escape))
        elif h in ('quantity', 'condition', 'sweep', 'A units') or 'window' in h:
            cols.append((h, L.tex_escape(h).replace(' (', r'\\('), 'l', L.tex_escape))
        elif h.startswith('R'):
            cols.append((h, L.tex_escape(h), 'r', lambda v: L.tex_num(v, 4)))
        else:
            cols.append((h, L.tex_escape(h).replace(' (', r'\\('), 'r', None))
    return rep.table(rows, cols, caption, rel)


def paths_segments(rows):
    """The regimes along f of a pinhole sweep: one row per run of equal labels, per condition."""
    out, groups = [], {}
    for r in rows:
        groups.setdefault((num(r.get('T (K)')), num(r.get('P_up (Pa)'))), []).append(r)
    for (T, P), grp in sorted(groups.items(), key=lambda kv: (kv[0][0] or 0, kv[0][1] or 0)):
        grp = sorted(grp, key=lambda r: num(r.get('f')) or 0)
        first = grp[0]
        start = 0
        for i in range(1, len(grp) + 1):
            if i == len(grp) or grp[i].get('label') != grp[start].get('label'):
                out.append({'T': T, 'P': P, 'label': grp[start].get('label'), 'f0': grp[start].get('f'),
                            'f1': grp[i - 1].get('f'), 'n': i - start, 'f_half': first.get('f_half'),
                            'PRF': first.get('PRF')})
                start = i
    return out


# the rows of a fits table that cover the whole sweep rather than one regime: the Arrhenius single fit of the
# temperature sweep, and the net exponent of the pressure sweep (the average local slope, not a fitted line)
WHOLE = ('single fit', 'net (whole range)')


def sweep_glance(rep, folder, word, levels):
    """The main result of each level: its whole-sweep row of J (single fit, or net n), or the per-regime fits
    when there is none."""
    rows = []
    for lv in levels:
        fits = rep.rows(f'{folder}/{lv}_{word}_fits.csv')
        if not fits:
            continue
        vcol = value_column(list(fits[0].keys()))
        j = [r for r in fits if r.get('quantity') == 'J']
        whole = [r for r in j if r.get('fit') in WHOLE]
        for r in (whole or j):
            same = [x for x in j if x.get('sweep', '') == r.get('sweep', '')]       # thickness: per swept layer
            rows.append({'level': lv, 'sweep': r.get('sweep', ''), 'condition': r.get('condition', ''),
                         'fit': r.get('fit', ''), 'value': num(r.get(vcol)) if vcol else None,
                         'n_regimes': len({x.get('fit') for x in same} - set(WHOLE))})
    return rows


def level_intro(name, lv, glance):
    """One line under a level's heading: its model, its main flux fit, how many regimes the sweep crosses."""
    out = f'{model_of(lv)}.'
    whole = [g for g in glance if g['fit'] in WHOLE and g['value'] is not None]
    what = {'temperature': '$E_a$', 'pressure': '$n$', 'thickness': '$m$', 'pinholes': '$n$'}[name]
    unit = ' kJ/mol' if name == 'temperature' else ''
    if whole:
        out += (' Net exponent of the flux over the whole range: ' if whole[0]['fit'] == 'net (whole range)'
                else ' Single fit of the flux over every point: ') + '; '.join(
            (f"{L.tex_escape(g['sweep'])} " if g['sweep'] else '') + f"{what} = {L.tex_num(g['value'])}{unit}"
            for g in whole) + '.'
    n = max((g['n_regimes'] for g in glance), default=0)
    if n > 1:
        out += f' The sweep crosses {n} regimes; their fits are in the tables below.'
    elif n == 1 and not whole:
        out += ' One regime over the whole sweep.'
    return out


SLOPES = {
    'pressure': [(1.0, 'about 1: a Henry oxide or the surface rate'), (0.5, 'about 0.5: Sieverts transport')],
    'thickness': [(-1.0, r'about $-1$: this thickness controls'), (0.0, 'about 0: something else controls')],
    'pinholes': [(0.0, 'about 0: pinholes do not change the flux here'), (1.0, 'about 1: pinholes carry the flux')],
}
VALUE_HEAD = {'temperature': r'$E_a$\\(kJ/mol)', 'pressure': r'$n$', 'thickness': r'$m$', 'pinholes': r'$n$'}
SWEEP_NAMES = {'oxide': 'oxide thickness swept', 'metal': 'metal thickness swept', 'total': 'whole wall scaled'}


def section_sweep(rep, name):
    folder, word, variants = SWEEPS[name]
    figs = rep.figures(folder)
    levels = level_sort(lv for lv in (level_of(s) for s in figs) if lv)
    levels += [lv for lv in level_sort(level_of(os.path.basename(f)) for f in rep.files(folder, '*.csv'))
               if lv and lv not in levels]
    levels = level_sort(levels)
    glance = sweep_glance(rep, folder, word, levels)
    keys = [f'{len(levels)} levels: ' + ', '.join(levels) + '.'] if levels else []
    vals = [g for g in glance if g['value'] is not None]
    if name == 'temperature' and vals:
        lo, hi = min(vals, key=lambda g: g['value']), max(vals, key=lambda g: g['value'])
        keys.append(f"Single-fit $E_a$ of the flux from {L.tex_num(lo['value'])} kJ/mol ({lo['level']}) to "
                    f"{L.tex_num(hi['value'])} kJ/mol ({hi['level']}).")
    elif vals:
        for g in vals:
            if g['fit'] in WHOLE or name in ('pinholes', 'thickness'):
                g['word'] = slope_word(g['value'], SLOPES[name])
        named = [g for g in vals if g.get('word')]
        what = 'Net $n$ of the whole range' if name == 'pressure' else 'Slope'
        for word_ in dict.fromkeys(g['word'] for g in named):
            lv = [g['level'] + (f" ({L.tex_escape(g['sweep'])})" if g['sweep'] else '') for g in named if g['word'] == word_]
            keys.append(f'{what} {word_}: ' + ', '.join(dict.fromkeys(lv)) + '.')
    multi = sorted({g['level'] for g in glance if g['n_regimes'] > 1}, key=lambda lv: LEVEL_RANK.get(lv, 99))
    if multi:
        keys.append('The regime changes inside the sweep for ' + ', '.join(multi)
                    + '; a whole-range value averages over those regimes, so quote the regime fits there too.')
    if name == 'pinholes':
        for lv in levels:
            seg = paths_segments(rep.rows(f'{folder}/{lv}_{word}_paths.csv', mark=False))
            if seg:
                s = seg[0]
                prf = num(s['PRF'])
                keys.append(f"{lv}: PRF = {L.tex_num(prf)}, so pinholes carry half the flux at $f_{{1/2}}$ = "
                            f"{L.tex_num(num(s['f_half']))}"
                            + (' --- the pinhole path passes \\emph{less} than intact oxide.' if prf is not None and prf < 1 else '.'))
    glance_table = rep.table(glance, [
        ('level', 'level', 'l', L.tex_escape), *([('sweep', 'sweep', 'l', L.tex_escape)] if name == 'thickness' else []),
        ('condition', 'condition', 'l', L.tex_escape), ('fit', 'fit', 'P5.6', L.tex_escape),
        ('value', VALUE_HEAD[name], 'r', None)],
        'The flux of each level at the operating point: ' + {
            'temperature': 'the single Arrhenius fit over every temperature.',
            'pressure': r'the net $n$ of the whole pressure range, $\Delta\ln J/\Delta\ln P$ between its ends.',
            'thickness': 'the fit of each regime, per swept layer.', 'pinholes': 'the fit of each regime.'}[name],
        group_by='level') if glance else ''
    parts = [('What was computed', f'For each level, the sweep at the operating point and at the low, middle and high end '
              'of the other ranges; a power-law or Arrhenius fit over the whole sweep and over each regime window. '
              'Each level is a subsection below, with its figures and fits.'),
             ('Key numbers', L.bullets(keys)), ('At a glance', glance_table)]
    for lv in levels:
        parts.append(('+' + level_name(lv), level_intro(name, lv, [g for g in glance if g['level'] == lv])))
        base = f'{lv}_{word}'
        # thickness: one figure per swept layer (<base>_oxide, _metal, _total); the others one per sweep
        layers = ['', '_oxide', '_metal', '_total'] if name == 'thickness' else ['']
        for v in ('', *(f'_{x}' for x in variants)):
            for s in layers:
                stem = f'{base}{s}{v}'
                if stem in figs:
                    what = f'{lv}, {SWEEP_NAMES[s[1:]]}' if s else lv
                    parts.append(('figure', v[1:] or 'main', figs[stem], f'({what}.)', 0.4 if s else 0.5))
        parts.append((None, fits_table(rep, f'{folder}/{base}_fits.csv', f'Fits of {lv} at the operating point.')))
        for v in variants:
            parts.append((None, fits_table(rep, f'{folder}/{base}_{v}_fits.csv', f'Fits of {lv}, three conditions ({v}).')))
        if name == 'pinholes':
            for suffix, what in (('', 'at the operating point'), *((f'_{v}', f'three conditions ({v})') for v in variants)):
                rel = f'{folder}/{base}{suffix}_paths.csv'
                seg = paths_segments(rep.rows(rel))
                if seg:
                    parts.append((None, rep.table(seg, [
                        ('T', '$T$ (K)', 'r', None), ('P', r'$P_\mathrm{up}$\\(Pa)', 'r', None),
                        ('label', 'regime', 'P5.0', L.tex_escape), ('f0', '$f$ from', 'r', None), ('f1', '$f$ to', 'r', None),
                        ('n', 'points', 'r', None), ('f_half', '$f_{1/2}$', 'r', None), ('PRF', 'PRF', 'r', None)],
                        f'Regimes of {lv} along $f$, {what}.', rel, group_by='T')))
    text = {**RT.SECTIONS[name], 'columns': [
        ('quantity', r'$J$, the flux, or $\Phi_\mathrm{app}$, the apparent permeability from the same flux.'),
        ('condition, window, points', 'where the sweep was made, and the range and number of points the fit used.'),
        ('fit', {'temperature': '\\emph{single fit}: one Arrhenius fit over every temperature, in the table only; '
                                'otherwise one regime, named by its dominant path and limiting step.',
                 'pressure': r'\emph{net (whole range)}: the net exponent of the whole trend, '
                             r'$\Delta\ln J/\Delta\ln P$ between the two ends of the range --- the average of the local '
                             r'exponent, not a fitted line, so it has no prefactor and is not drawn; otherwise one regime, '
                             'named by its dominant path and limiting step.'}.get(
            name, 'one regime, named by its dominant path and limiting step.')),
        (VALUE_HEAD[name].replace(r'\\', ' '), {'temperature': 'the activation energy of the Arrhenius fit.',
                                                'pressure': '$n$ of $J = A\\,P^{\\,n}$.',
                                                'thickness': '$m$ of $J \\propto L^{\\,m}$, for the swept layer.',
                                                'pinholes': '$n$ of $J \\propto f^{\\,n}$.'}[name]),
        ('A', 'the prefactor of the fit, in SI units.'),
    ] + ([('regime, $f$ from, $f$ to', 'one row per run of the same regime along the sweep.'),
          ('$f_{1/2}$, PRF', 'the area fraction at which pinholes carry half the flux, $1/(1+\\mathrm{PRF})$, and '
           '$\\mathrm{PRF} = J_p/J_i$.')] if name == 'pinholes' else [])}
    return rep.section(name, text, parts)


# ── sensitivity ──────────────────────────────────────────────────────────────
def section_sensitivity(rep, name):
    figs_dir, data_dir = SECTION_FOLDERS[name]
    figs = rep.figures(figs_dir)
    clusters = rep.rows(f'{data_dir}/master_clusters.csv')
    drivers = rep.rows(f'{data_dir}/delta_over_floor_by_T.csv')
    keys, parts = [], []
    if clusters:
        count = Counter((r.get('regime'), num(r.get('temperature_K'))) for r in clusters)
        regimes = sorted({k[0] for k in count})
        Ts = sorted({k[1] for k in count if k[1] is not None})
        keys.append(f'{len(clusters)} scan draws in {len(regimes)} regimes (' + ', '.join(
            f'{L.tex_escape(g)} {sum(v for k, v in count.items() if k[0] == g)}' for g in regimes)
            + (f') at {", ".join(f"{t:g}" for t in Ts)} K.' if Ts else ').'))
        grid = [{'regime': g, 'T': t, 'n': count.get((g, t), 0)} for g in regimes for t in Ts]
        parts.append(('Cluster sizes', rep.pivot(grid, 'regime', 'T', 'n', 'Draws per regime and temperature.', 'regime',
                                                 f'{data_dir}/master_clusters.csv', col_head=lambda t: f'{float(t):g} K')))
    if drivers:
        tcols = [c for c in drivers[0].keys() if re.fullmatch(r'\d+(\.\d+)?K', c)]
        above = [r for r in drivers if (num(r.get('best')) or 0) >= 1]
        for g in dict.fromkeys(r['regime'] for r in drivers):
            top = sorted((r for r in above if r['regime'] == g), key=lambda r: -num(r['best']))
            n_all = sum(r['regime'] == g for r in drivers)
            keys.append(f'{L.tex_escape(g)}: {len(top)} of {n_all} parameters above the dummy floor'
                        + (': ' + ', '.join(f"{L.tex_file(r['parameter'])} ({L.tex_num(num(r['best']))})" for r in top[:4])
                           if top else '') + '.')
        above = sorted(above, key=lambda r: (r['regime'], -num(r['best'])))
        parts.append(('Drivers above the dummy floor', rep.table(above, [
            ('regime', 'regime', 'l', L.tex_escape), ('parameter', 'parameter', 'l', L.tex_file),
            *[(c, r'$\delta$/floor\\' + c.replace('K', ' K'), 'r', None) for c in tcols], ('best', 'best', 'r', None)],
            r'Parameters whose $\delta$ for the flux exceeds the dummy floor at some temperature, by regime.',
            f'{data_dir}/delta_over_floor_by_T.csv', group_by='regime')))
    head = [('What was computed', 'A Latin-hypercube scan of the wrapper aimed at each regime with its preset, the '
             'draws split by regime, then given-data Borgonovo $\\delta$ and PAWN per parameter within each regime, '
             'against a dummy parameter\\textquotesingle s floor.'),
            ('Key numbers', L.bullets(keys) or r'The tables of \texttt{' + L.tex_escape(data_dir) + '/} are missing.')]
    order = [('regime_exploration', r'regime_exploration_.*'), ('givendata', r'givendata_.*'),
             ('comparison_heatmap', r'comparison_heatmap_.*'), ('regime_geometry', r'regime_geometry_.*')]

    def t_of(stem):
        m = re.search(r'_(\d+)K$', stem)
        return (int(m.group(1)) if m else 0, stem)
    for key, pat in order:
        for stem in sorted((s for s in figs if re.fullmatch(pat, s)), key=t_of):
            parts.append(('figure', key, figs[stem], f'({L.tex_escape(stem)}.)', 0.5))
    return rep.section(name, RT.SECTIONS[name], head + parts)


def section_regimes(rep):
    files = rep.outputs('regimes')
    rows = [{'file': f, 'kb': os.path.getsize(os.path.join(rep.out, f)) / 1024} for f in files]
    for f in files:
        rep.used.add(f)
    parts = [('Interactive figures', rep.table(rows, [('file', 'file', 'P10', L.tex_escape), ('kb', 'size (kB)', 'r', None)],
                                               'Parallel-coordinates pages: open each in a browser.', full=True))]
    return rep.section('regimes', RT.SECTIONS['regimes'], parts)


def build_section(rep, name):
    if name == 'levels':
        return section_levels(rep)
    if name == 'surface':
        return section_surface(rep)
    if name in SWEEPS:
        return section_sweep(rep, name)
    if name.startswith('sensitivity'):
        return section_sensitivity(rep, name)
    return section_regimes(rep)


# ── appendices ───────────────────────────────────────────────────────────────
def appendix_files(rep):
    rows = []
    for name in SECTION_ORDER:
        for folder in SECTION_FOLDERS[name]:
            fs = rep.files(folder)
            if not fs:
                continue
            kinds = Counter(os.path.splitext(f)[1].lower().lstrip('.') for f in fs)
            figs = rep.figures(folder)
            rows.append({'folder': folder + '/', 'kinds': ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
                         'shown': f'{sum(f in rep.used for f in fs)} of {len(fs)}'
                                  + (f' ({sum(f in rep.used for f in figs.values())} of {len(figs)} figures)' if figs else ''),
                         'section': (r'\ref{sec:' + name + '}') if name in rep.built else '--'})
    return '\n\n'.join([r'\section{Notebook outputs}', r'\label{app:files}',
                        'Every results folder of the study, ' + L.tex_path(rep.out) + '. A figure is shown as its PDF '
                        'when the notebook wrote one, else as its PNG; the subfolder \\texttt{sa\\_results/scans/} and '
                        'other nested folders are not listed.',
                        rep.table(rows, [('folder', 'folder', 'P5.4', L.tex_escape), ('kinds', 'files', 'P3.6', L.tex_escape),
                                         ('shown', 'shown', 'P4.2', L.tex_escape), ('section', 'section', 'l', lambda s: s)],
                                  'Results folders.', full=True)])


def appendix_config(rep):
    if not os.path.isfile(rep.config):
        return ''
    os.makedirs(rep.dir, exist_ok=True)
    with open(rep.config, errors='replace') as fh:
        text = L.ascii_text(fh.read())
    with open(os.path.join(rep.dir, 'study_config_listing.py'), 'w') as fh:
        fh.write(text)
    return '\n\n'.join([r'\section{The input file of this study}', r'\label{app:config}',
                        r'\texttt{study\_config.py} as it is now (non-ASCII characters transliterated).',
                        r'\lstinputlisting{study_config_listing.py}'])


def appendix_docs(rep):
    pages = {}
    for title, path, heading in rep.pointers:
        p = pages.setdefault(path, {'headings': [], 'sections': []})
        if heading and heading not in p['headings']:
            p['headings'].append(heading)
        if title not in p['sections']:
            p['sections'].append(title)
    rows = [{'page': rep.pointer_link(path), 'headings': '; '.join(L.tex_heading(h) for h in p['headings']) or '(whole page)',
             'sections': '; '.join(p['sections'])} for path, p in sorted(pages.items())]
    out = [r'\section{Documentation map}', r'\label{app:docs}',
           'Where each part of this report is explained. Paths are from the repository root'
           + ('; each links to the page on GitHub at the commit the report was built from.' if rep.web_base else '.'),
           rep.table(rows, [('page', 'page', 'P5.2', lambda s: s), ('headings', 'sections of the page used', 'P6.4', lambda s: s),
                            ('sections', 'pointed to from', 'P4.4', lambda s: s)], 'Documentation pages.', full=True)]
    if rep.broken:
        out.append(r'\textbf{Pointers that did not resolve} (a page or heading was renamed; update report\_text.py):'
                   + L.bullets([L.tex_escape(b) for b in sorted(set(rep.broken))]))
    return '\n\n'.join(out)


# ── running the notebooks ────────────────────────────────────────────────────
# The order matters: the parallel-coordinates notebooks read the scans the sensitivity
# notebooks write (docs/caveats.md, "Run the SA notebooks before the parallel-coordinates ones").
NOTEBOOK_ORDER = ['Proposal', 'Surface_proposal', 'effect-of-temperature', 'effect-of-pressure',
                  'effect-of-thickness', 'effect-of-pinholes', 'sensitivity_regime_L5', 'sensitivity_regime_L5L6',
                  'regime_parallel_coords_L5', 'regime_parallel_coords']


def kernel_for(python):
    """The Jupyter kernel whose interpreter is `python`, or None (the notebooks' own 'python3'
    may be another environment without the package; tests/test_notebook_execution.py)."""
    try:
        out = subprocess.run([python, '-m', 'jupyter', 'kernelspec', 'list', '--json'],
                             capture_output=True, text=True, timeout=120)
        specs = json.loads(out.stdout).get('kernelspecs', {}) if out.returncode == 0 else {}
    except (OSError, subprocess.SubprocessError, ValueError):
        return None
    want = os.path.realpath(python)
    for name, spec in specs.items():
        argv = spec.get('spec', {}).get('argv') or []
        if argv and os.path.realpath(argv[0]) == want:
            return name
    return None


def run_notebooks(study_dir, only=None, log=print, python=None):
    """Execute the study's notebooks in order, headless, in the study folder. The executed copies
    go to report/executed/, so the committed notebooks are not rewritten; their figures and tables
    go to the usual results folders. Returns the notebooks that failed."""
    python = python or sys.executable
    kernel = kernel_for(python)
    if kernel is None:
        raise L.ReportError(f'no Jupyter kernel runs {python}; register one with '
                            f'`{python} -m ipykernel install --user --name <env>`')
    names = [n for n in NOTEBOOK_ORDER if os.path.isfile(os.path.join(study_dir, n + '.ipynb'))]
    names += sorted(os.path.splitext(os.path.basename(p))[0] for p in glob.glob(os.path.join(study_dir, '*.ipynb'))
                    if os.path.splitext(os.path.basename(p))[0] not in names)
    if only:
        names = [n for n in names if n in only]
    out_dir = os.path.join(study_dir, 'report', 'executed')
    os.makedirs(out_dir, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if k != 'PERMEATION_STUDY'}     # the notebooks insist on their folder
    env['MPLBACKEND'] = 'Agg'
    failed = []
    for n in names:
        log(f'report: running {n}.ipynb ...')
        proc = subprocess.run([python, '-m', 'jupyter', 'nbconvert', '--execute', '--to', 'notebook',
                               f'--output-dir={out_dir}', f'--ExecutePreprocessor.kernel_name={kernel}',
                               '--ExecutePreprocessor.timeout=-1', n + '.ipynb'],
                              cwd=study_dir, env=env, capture_output=True, text=True, errors='replace')
        if proc.returncode != 0:
            failed.append(n)
            log(f'report: {n}.ipynb failed:\n' + proc.stderr[-2000:])
    return failed


# ── build ────────────────────────────────────────────────────────────────────
def run_data(study_dir, log=print, python=None):
    """calculations.report_data in its own process, with this study bound. True if it worked."""
    cfg = os.path.join(study_dir, 'study_config.py')
    rel = os.path.relpath(cfg, REPO_ROOT)
    env = dict(os.environ, PERMEATION_STUDY=cfg if rel.startswith('..') else rel)
    log(f'report: computing the operating point of {os.path.basename(study_dir)} ...')
    try:
        proc = subprocess.run([python or sys.executable, '-m', 'calculations.report_data', study_dir], cwd=REPO_ROOT,
                              env=env, capture_output=True, text=True, errors='replace', timeout=900)
    except subprocess.TimeoutExpired:
        log('report: the operating-point step took longer than 900 s; building without it')
        return False
    if proc.returncode != 0:
        log('report: the operating-point step failed; building without it:\n' + (proc.stderr or proc.stdout)[-3000:])
        return False
    return True


def write_tex(rep):
    title, title_tex = title_page(rep)
    rep.built = [s for s in SECTION_ORDER if rep.outputs(s)]
    sections = [build_section(rep, s) for s in rep.built]
    summary = run_summary(rep)
    appendices = [a for a in (appendix_files(rep), appendix_config(rep)) if a]
    appendices.append(appendix_docs(rep))
    parts = [L.preamble(title), r'\begin{document}', title_tex, r'\tableofcontents', r'\clearpage', summary,
             *[r'\clearpage' + '\n' + s for s in sections],
             r'\clearpage', r'\appendix', *[r'\clearpage' + '\n' + a for a in appendices], r'\end{document}']
    os.makedirs(rep.dir, exist_ok=True)
    path = os.path.join(rep.dir, f'{rep.name}_report.tex')
    with open(path, 'w') as fh:
        fh.write('\n\n'.join(parts) + '\n')
    return path


def build(study_dir, latex=True, data=True, latexmk=None, log=print, web_base='auto', python=None):
    """Write <study>/report/<study>_report.tex and build the PDF. data: first compute the
    operating point (report_data, in its own process). Returns the PDF path, or the .tex path
    when the PDF is not built."""
    study_dir = os.path.abspath(study_dir)
    if not os.path.isfile(os.path.join(study_dir, 'study_config.py')):
        raise L.ReportError(f'{study_dir} is not a study folder: it has no study_config.py')
    if data:
        run_data(study_dir, log, python)
    rep = StudyReport(study_dir, log=log, web_base=web_base)
    tex = write_tex(rep)
    exe = L.latexmk_path(latexmk) if latex else None
    if exe is None:
        log(f'report: wrote {tex}' + ('' if not latex else
                                       f'; LaTeX (latexmk) not found. Build it with: cd {rep.dir} && latexmk -pdf '
                                       f'{os.path.basename(tex)}'))
        return tex
    pdf = L.build_pdf(tex, exe, log)
    log(f'report: {pdf}')
    return pdf


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('studies', nargs='+', help='study folders, e.g. Guo_etal_2025')
    ap.add_argument('--tex-only', action='store_true', help='write the .tex but do not run LaTeX')
    ap.add_argument('--no-data', action='store_true', help='do not recompute report/data/ (the operating point)')
    ap.add_argument('--latexmk', default=None, help='path of latexmk (default: on PATH or /Library/TeX/texbin)')
    ap.add_argument('--run-notebooks', nargs='?', const='all', default=None, metavar='NAMES',
                    help='first execute the notebooks (all, in dependency order, or a comma-separated list of '
                         'names without .ipynb); the sensitivity scans take long')
    a = ap.parse_args(argv)
    failed = []
    for study in a.studies:
        try:
            if a.run_notebooks:
                only = None if a.run_notebooks == 'all' else set(a.run_notebooks.split(','))
                bad = run_notebooks(os.path.abspath(study), only)
                if bad:
                    print(f'report: {study}: notebooks that failed: {", ".join(bad)} (their previous outputs, '
                          'if any, are used)')
            build(study, latex=not a.tex_only, data=not a.no_data, latexmk=a.latexmk)
        except L.ReportError as exc:
            print(f'report: {study}: {exc}')
            failed.append(study)
    if failed:
        print(f'report: failed for {", ".join(failed)}')
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
