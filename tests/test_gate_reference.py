"""The container gate's `--emit` mode, which regenerates its own reference values.

`container/model_smoke_test.py` asserts 13 pinned values per study. Whenever the
operating point changes those values have to be re-measured, and there is no other tool
that does it — `container/revalidate_reference.py` is a cross-interpreter diff script
with its own hardcoded inputs.

`--emit` avoids a second copy of the model calls by turning `REF` into a recorder: each
lookup returns a sentinel and remembers its key, and `close()`/`equal()` store the value
they were handed. That pairing is the subtle part — it relies on the lookup happening
immediately before the call it belongs to — so the test below checks the emitted values
against the committed ones, not merely that something was printed.
"""
import ast
import os
import re
import subprocess
import sys

import pytest

from conftest import KNOWN_STUDIES, ROOT

GATE = ROOT / "container" / "model_smoke_test.py"


def run_gate(study, *args):
    env = dict(os.environ)
    env["MHI_STUDY"] = f"{study}/study_config.py"
    env["MHI_STUDY_QUIET"] = "1"
    env["PYTHONPATH"] = os.pathsep.join([str(ROOT), env.get("PYTHONPATH", "")])
    return subprocess.run([sys.executable, str(GATE), *args, "-w", "."],
                          cwd=str(ROOT), env=env, capture_output=True, text=True)


def emitted_block(stdout):
    """The `'study': {...},` block `--emit` prints, as a dict."""
    match = re.search(r"^    '([^']+)': \{$(.*?)^    \},$", stdout, re.S | re.M)
    assert match, f"no REFERENCE block in output:\n{stdout[-600:]}"
    return match.group(1), ast.literal_eval("{" + match.group(2) + "}")


def committed_reference(study_name):
    source = GATE.read_text()
    block = source[source.index("REFERENCE = {"):source.index("failures:")]
    return ast.literal_eval(block[block.index("{"):block.rindex("}") + 1])[study_name]


# Guo_etal_2025 stands in for all of them: the gate runs the same calls whichever study
# is active, and running it four times would quadruple the suite's runtime for no signal.
STUDY = "Guo_etal_2025"


def test_emit_prints_a_parseable_reference_block():
    proc = run_gate(STUDY, "--emit")
    assert proc.returncode == 0, proc.stdout[-800:] + proc.stderr[-400:]
    name, emitted = emitted_block(proc.stdout)
    assert name == KNOWN_STUDIES[STUDY][0]
    assert emitted, "emitted an empty block"


def test_emit_covers_exactly_the_committed_keys():
    _, emitted = emitted_block(run_gate(STUDY, "--emit").stdout)
    assert set(emitted) == set(committed_reference(KNOWN_STUDIES[STUDY][0]))


def test_emitted_values_match_the_committed_ones():
    """Guards the sentinel pairing: a value recorded against the wrong key shows up here.

    It also means the committed `REFERENCE` is, by definition, what the model produces
    today — so a re-baseline is incomplete until this passes again.
    """
    _, emitted = emitted_block(run_gate(STUDY, "--emit").stdout)
    committed = committed_reference(KNOWN_STUDIES[STUDY][0])
    for key, got in emitted.items():
        want = committed[key]
        if isinstance(want, float):
            assert got == pytest.approx(want, rel=1e-9), key
        else:
            assert got == want, key


def test_emit_does_not_assert_and_the_gate_still_does():
    """The two modes must stay distinct: emit measures, the bare gate checks."""
    assert "nothing is being asserted" in run_gate(STUDY, "--emit").stdout
    plain = run_gate(STUDY)
    assert plain.returncode == 0, plain.stdout[-800:]
    assert "SCIENTIFIC SMOKE TEST PASSED" in plain.stdout
    assert "Paste into REFERENCE" not in plain.stdout
