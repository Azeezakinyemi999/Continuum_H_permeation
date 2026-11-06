"""The study report (calculations.report): its documentation pointers resolve, it builds
from a study folder made here, a section appears only when its notebook has outputs, and
the operating-point step computes every level for a real study.

The PDF test runs only where latexmk is installed; everything else needs no LaTeX.
"""
import csv
import json
import os
import struct
import subprocess
import sys
import zlib

import pytest

from conftest import ROOT

from calculations import report
from calculations import report_latex as L
from calculations import report_text as RT


# ── a made-up study folder ──────────────────────────────────────────────────
def tiny_png(path):
    """A valid 1 x 1 white PNG, without PIL."""
    def chunk(kind, data):
        return struct.pack('>I', len(data)) + kind + data + struct.pack('>I', zlib.crc32(kind + data) & 0xffffffff)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'wb') as fh:
        fh.write(b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', 1, 1, 8, 2, 0, 0, 0))
                 + chunk(b'IDAT', zlib.compress(b'\x00\xff\xff\xff')) + chunk(b'IEND', b''))


def write_csv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


FIT_HEAD = ['quantity', 'condition', 'fit', 'T window (K)', 'points', 'E_a (kJ/mol)', 'A', 'A units']
P_HEAD = ['quantity', 'condition', 'fit', 'P window (Pa)', 'points', 'n (J) / slope (Φ_app)', 'A (SI)']
L_HEAD = ['quantity', 'sweep', 'condition', 'fit', 'L window (m)', 'points', 'm (J) / slope (Φ_app)', 'A (SI)']


def study(d, sensitivity=True):
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, 'study_config.py'), 'w') as fh:
        fh.write("STUDY_NAME = 'Toy_study'\nCONDITIONS = {'T_operating': 873}\n")
    tiny_png(os.path.join(d, 'Proposal-results', 'Level1_Perfect_Metal.png'))
    tiny_png(os.path.join(d, 'Surface_proposal-results', 'L1L6_core_validation.png'))
    tiny_png(os.path.join(d, 'effect-of-temperature-results', 'L1_arrhenius.png'))
    tiny_png(os.path.join(d, 'effect-of-temperature-results', 'L5L6_arrhenius.png'))
    write_csv(os.path.join(d, 'effect-of-temperature-results', 'L1_arrhenius_fits.csv'), FIT_HEAD,
              [['J', '1e5 Pa', 'single fit', '573–1273', 40, 61.75, 0.1, 'mol/m^2/s']])
    write_csv(os.path.join(d, 'effect-of-temperature-results', 'L5L6_arrhenius_fits.csv'), FIT_HEAD,
              [['J', '1e5 Pa', 'single fit', '573–1273', 40, 136.2, 57.3, 'mol/m^2/s'],
               ['J', '1e5 Pa', 'Pinholes dominant: metal-diffusion limited, trap-controlled', '573–717', 9, 88.8, 0.0058, 'mol/m^2/s'],
               ['J', '1e5 Pa', 'Intact oxide dominant: oxide-diffusion limited', '730–900', 14, 210.0, 1e3, 'mol/m^2/s']])
    tiny_png(os.path.join(d, 'effect-of-pressure-results', 'L5L6_pressure.png'))
    write_csv(os.path.join(d, 'effect-of-pressure-results', 'L5L6_pressure_fits.csv'), P_HEAD,
              [['J', '873 K', 'net (whole range)', '1.0e-07–1.0e+07', 40, 0.7588, ''],
               ['J', '873 K', 'Intact oxide dominant: oxide-surface limited', '1.0e-07–6.6e-01', 20, 0.991, 2.1e-9],
               ['J', '873 K', 'Intact oxide dominant: oxide-diffusion limited', '3.5e+00–1.0e+07', 19, 0.513, 2.5e-9]])
    # Level 1 has one layer: its metal sweep only. Level 2b has both: oxide, metal and the whole wall.
    tiny_png(os.path.join(d, 'effect-of-thickness-results', 'L1_thickness_metal.png'))
    write_csv(os.path.join(d, 'effect-of-thickness-results', 'L1_thickness_fits.csv'), L_HEAD,
              [['J', 'metal', '873 K, 1e5 Pa', 'Metal-diffusion limited', '1.00e-04–5.00e-03', 40, -1.0, 4.3e-9]])
    for sweep in ('oxide', 'metal', 'total'):
        tiny_png(os.path.join(d, 'effect-of-thickness-results', f'L2b_thickness_{sweep}.png'))
    write_csv(os.path.join(d, 'effect-of-thickness-results', 'L2b_thickness_fits.csv'), L_HEAD,
              [['J', sweep, '873 K, 1e5 Pa', 'Oxide-diffusion limited', '1.00e-09–1.00e-05', 40, -0.9, 4e-10]
               for sweep in ('oxide', 'metal', 'total')])
    write_csv(os.path.join(d, 'effect-of-pinholes-results', 'L3L6_pinholes_paths.csv'),
              ['level', 'T (K)', 'P_up (Pa)', 'f', 'label', 'dominant_path', 'J_total', 'PRF', 'f_half'],
              [['L3L6', 873, 1e5, 1e-5, 'intact/oxide', 'intact', 1e-6, 0.5, 0.667],
               ['L3L6', 873, 1e5, 1e-2, 'intact/oxide', 'intact', 1e-6, 0.5, 0.667],
               ['L3L6', 873, 1e5, 0.5, 'pinholes/metal', 'pinholes', 8e-7, 0.5, 0.667]])
    write_csv(os.path.join(d, 'effect-of-pinholes-results', 'L3L6_pinholes_fits.csv'),
              ['quantity', 'condition', 'fit', 'f window', 'points', 'n (J, Φ_app ∝ f^n)', 'A (SI)'],
              [['J', '873 K, 1e5 Pa', 'Intact oxide dominant: oxide-diffusion limited', '1e-5–0.5', 40, -0.0004, 9.6e-7]])
    if sensitivity:
        write_csv(os.path.join(d, 'sa_results', 'delta_over_floor_by_T.csv'),
                  ['regime', 'parameter', '773K', '1073K', 'best'],
                  [['metal', 'D_ref', 1.4, 1.3, 1.4], ['metal', 'K_s_ref', 0.5, 0.6, 0.6], ['surface', 'k_diss_metal_ref', 2.2, 2.0, 2.2]])
        write_csv(os.path.join(d, 'sa_results', 'master_clusters.csv'), ['D_ref', 'regime', 'temperature_K'],
                  [[1, 'metal', 773], [2, 'metal', 1073], [3, 'surface', 773]])
        tiny_png(os.path.join(d, 'sensitivity_regime_L5L6-results', 'comparison_heatmap_773K.png'))
    data = os.path.join(d, 'report', 'data')
    write_csv(os.path.join(data, 'operating_point.csv'), ['level', 'quantity', 'value', 'units', 'source'],
              [['L1', 'flux', 4.3e-6, 'mol m-2 s-1', 'x'], ['L2b', 'flux', 4.2e-6, 'mol m-2 s-1', 'x'],
               ['L2b', 'resistance_ratio', 0.005, '-', 'x'], ['L5', 'flux', 3.2e-6, 'mol m-2 s-1', 'x'],
               ['L5', 'regime', 'metal', '', 'x'], ['L5L6', 'flux', 9.1e-7, 'mol m-2 s-1', 'x'],
               ['L5L6', 'regime', 'oxide', '', 'x'], ['L5L6', 'PRF', '', '-', 'x']])
    with open(os.path.join(data, 'study_facts.json'), 'w') as fh:
        json.dump({'study_name': 'Toy_study', 'conditions': {'T_operating': 873, 'P_upstream': 1e5,
                                                             'T_range': [573, 1273]},
                   'config_problems': ['L5: default D_ref=2.8e-11 lies OUTSIDE its sweep range'],
                   'model_warnings': ['Warning: Temperature 873K outside validated range [468, 823]K'],
                   'computed_at': '2026-11-05T12:00:00'}, fh)
    return d


def tex_of(path):
    with open(path) as fh:
        return fh.read()


def build_tex(d):
    return tex_of(report.build(d, latex=False, data=False, web_base=None, log=lambda *_: None))


# ── documentation pointers ──────────────────────────────────────────────────
def all_pointers():
    for block in (RT.SUMMARY, *RT.SECTIONS.values()):
        for path, heading, _ in block['docs']:
            yield block['title'], path, heading


@pytest.mark.parametrize('title,path,heading', list(all_pointers()))
def test_every_documentation_pointer_resolves(title, path, heading):
    assert L.check_pointer(str(ROOT), path, heading) is None, f'{title}: {path} / {heading}'


def test_every_section_points_into_the_docs():
    for name in report.SECTION_ORDER:
        assert any(p.startswith('docs/') for p, _, _ in RT.SECTIONS[name]['docs']), name


@pytest.mark.parametrize('heading,anchor', [
    ('`PRF` is `nan` from the Level 5L6 wrapper', 'prf-is-nan-from-the-level-5l6-wrapper'),
    ("Henry's law, not Sieverts", 'henrys-law-not-sieverts'),
    ('Below the floor does not mean "no physical effect"', 'below-the-floor-does-not-mean-no-physical-effect'),
    ('Level 2a — the perfect oxide', 'level-2a--the-perfect-oxide'),
])
def test_github_slug(heading, anchor):
    assert L.github_slug(heading) == anchor


# ── helpers ─────────────────────────────────────────────────────────────────
def test_level_of():
    assert report.level_of('L5L6_arrhenius_3P') == 'L5L6'
    assert report.level_of('L2aL6_core_validation') == 'L2aL6'
    assert report.level_of('Level2b_Oxide_Metal') == 'L2b'
    assert report.level_of('regime_geometry_oxide') is None
    assert report.level_sort(['L5L6', 'L1', 'L2b', 'L1L6']) == ['L1', 'L2b', 'L1L6', 'L5L6']


def test_paths_segments_split_the_sweep_by_regime():
    rows = [{'T (K)': '873', 'P_up (Pa)': '1e5', 'f': f, 'label': lab, 'PRF': '0.5', 'f_half': '0.667'}
            for f, lab in (('1e-5', 'a'), ('1e-3', 'a'), ('1e-1', 'b'), ('0.5', 'b'))]
    seg = report.paths_segments(rows)
    assert [(s['label'], s['f0'], s['f1'], s['n']) for s in seg] == [('a', '1e-5', '1e-3', 2), ('b', '1e-1', '0.5', 2)]


# ── building ────────────────────────────────────────────────────────────────
def test_study_folder_gives_its_sections(tmp_path):
    body = build_tex(study(str(tmp_path / 'Toy_study')))
    for name in ('levels', 'surface', 'temperature', 'pressure', 'thickness', 'pinholes', 'sensitivity_L5L6'):
        assert r'\label{sec:' + name + '}' in body, name
    for name in ('sensitivity_L5', 'regimes'):
        assert r'\label{sec:' + name + '}' not in body, name
    assert body.count(r'\paragraph{Where to read more.}') == 7 + 1        # + the run summary
    # no R² anywhere; the pressure sweep's whole-range value is the net n; thickness has one figure per layer
    assert '$R^2$' not in body and 'R²' not in body
    assert r'Net exponent of the flux over the whole range: $n$ = $0.759$' in body
    assert '(L1, metal thickness swept.)' in body
    assert '(L1, oxide thickness swept.)' not in body and '(L1, whole wall scaled.)' not in body
    for sweep in ('oxide thickness swept', 'metal thickness swept', 'whole wall scaled'):
        assert f'(L2b, {sweep}.)' in body, sweep
    assert 'Pointers that did not resolve' not in body
    # key numbers from the tables, and the caveats from the facts
    assert r'Single-fit $E_a$ of the flux from $61.8$ kJ/mol (L1) to $136$ kJ/mol (L5L6)' in body
    assert 'The regime changes inside the sweep for L5L6' in body
    assert r'pinholes carry half the flux at $f_{1/2}$ = $0.667$' in body and r'passes \emph{less} than intact oxide' in body
    assert r'metal: 1 of 2 parameters above the dummy floor' in body
    assert 'lies OUTSIDE its sweep range' in body and 'outside validated range' in body
    assert 'Level 5 and Level 5 + surface name different regimes' in body


def test_without_the_operating_point_the_report_says_so(tmp_path):
    d = study(str(tmp_path / 'Toy_study'), sensitivity=False)
    os.remove(os.path.join(d, 'report', 'data', 'operating_point.csv'))
    body = build_tex(d)
    assert r'\textbf{Not computed}' in body


def test_a_folder_without_study_config_is_refused(tmp_path):
    with pytest.raises(L.ReportError):
        report.build(str(tmp_path), latex=False, data=False, web_base=None, log=lambda *_: None)


def test_operating_point_of_a_real_study(tmp_path):
    """report_data in its own process, the way report.build runs it, for the study the suite gates on."""
    env = dict(os.environ, PERMEATION_STUDY='Guo_etal_2025/study_config.py')
    code = ('import sys; from calculations import report_data as d; '
            'd.write(sys.argv[1], sys.argv[2])')
    proc = subprocess.run([sys.executable, '-c', code, str(ROOT / 'Guo_etal_2025'), str(tmp_path)],
                          cwd=ROOT, env=env, capture_output=True, text=True, timeout=900)
    assert proc.returncode == 0, proc.stderr[-3000:]
    with open(tmp_path / 'operating_point.csv', newline='') as fh:
        rows = list(csv.DictReader(fh))
    flux = {r['level']: float(r['value']) for r in rows if r['quantity'] == 'flux'}
    assert set(flux) == {'L1', 'L2b', 'L1L6', 'L5', 'L5L6'}
    assert all(v > 0 for v in flux.values())
    assert flux['L2b'] <= flux['L1']                  # an oxide cannot raise the flux of the bare metal
    facts = json.load(open(tmp_path / 'study_facts.json'))
    assert facts['study_name'] and facts['conditions']['T_operating'] > 0 and not facts['errors']


@pytest.mark.skipif(L.latexmk_path() is None, reason='LaTeX (latexmk) not installed')
def test_pdf_builds(tmp_path):
    d = study(str(tmp_path / 'Toy_study'))
    pdf = report.build(d, latex=True, data=False, web_base=None, log=lambda *_: None)
    assert pdf.endswith('.pdf') and os.path.getsize(pdf) > 10000
