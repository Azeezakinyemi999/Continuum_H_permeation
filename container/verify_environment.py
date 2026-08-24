#!/usr/bin/env python3
"""Gate 1 — is the software environment the one we validated?

Checks the *whole* locked environment, not a hand-maintained subset, by reading
the same requirements.lock.txt the image was built from. That way this gate can
never drift out of step with the build.

Note versions are read via importlib.metadata rather than module __version__
attributes: SALib has no __version__ at all, so `SALib.__version__` raises
AttributeError. Reading distribution metadata also reports what pip actually
installed rather than what a package claims about itself.

Exit 0 = pass, 1 = fail. Prints every check so a failure is diagnosable from a
copy-pasted log.
"""
import sys
from importlib.metadata import PackageNotFoundError, version as dist_version
from pathlib import Path

LOCK = Path(__file__).with_name("requirements.lock.txt")
if not LOCK.exists():                      # in-image location
    LOCK = Path("/opt/requirements.lock.txt")

EXPECTED_PYTHON = (3, 12)

# Results depend on these exact versions; the reference values in
# PACKAGING_1.0.0.md are only meaningful against them. Listed separately so a
# mismatch here is reported as scientifically significant, not just drift.
NUMERIC_CRITICAL = {"numpy", "scipy", "pandas", "SALib"}

failures: list[str] = []


def check(ok: bool, label: str, detail: str = "") -> bool:
    print(f"{'OK  ' if ok else 'FAIL'}  {label}{('  ' + detail) if detail else ''}")
    if not ok:
        failures.append(label)
    return ok


print("=" * 78)
print("HYDROGEN MODEL — ENVIRONMENT VERIFICATION (gate 1 of 2)")
print("=" * 78)

print(f"\npython   {sys.version.split()[0]}  ({sys.executable})")
check(sys.version_info[:2] == EXPECTED_PYTHON,
      "python version",
      f"expected {'.'.join(map(str, EXPECTED_PYTHON))}.x, got "
      f"{sys.version_info.major}.{sys.version_info.minor}")

# ---------------------------------------------------------------- locked deps
if not LOCK.exists():
    check(False, "requirements.lock.txt present", f"looked in {LOCK}")
else:
    pins = {}
    for line in LOCK.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "==" not in line:
            check(False, "lock line parseable", repr(line))
            continue
        name, _, ver = line.partition("==")
        pins[name] = ver

    print(f"\nlocked packages: {len(pins)}  (from {LOCK})")
    mismatched, missing = [], []
    for name in sorted(pins, key=str.lower):
        want = pins[name]
        try:
            got = dist_version(name)
        except PackageNotFoundError:
            missing.append(name)
            continue
        if got != want:
            mismatched.append((name, want, got))

    for name in missing:
        tag = " [NUMERIC]" if name in NUMERIC_CRITICAL else ""
        check(False, f"{name} installed{tag}", "NOT INSTALLED")
    for name, want, got in mismatched:
        tag = " [NUMERIC]" if name in NUMERIC_CRITICAL else ""
        check(False, f"{name} version{tag}", f"expected {want}, got {got}")

    if not missing and not mismatched:
        check(True, f"all {len(pins)} locked versions match")

    # Always echo the ones that matter, pass or fail, so logs are comparable.
    print("\nnumerically significant:")
    for name in sorted(NUMERIC_CRITICAL, key=str.lower):
        try:
            print(f"        {name:<12} {dist_version(name)}")
        except PackageNotFoundError:
            print(f"        {name:<12} NOT INSTALLED")

# ------------------------------------------------------------- functional checks
print("\nfunctional:")

try:
    from scipy.optimize import brentq
    root = brentq(lambda x: x * x - 4.0, 0.0, 5.0)
    check(abs(root - 2.0) < 1e-12, "scipy.optimize.brentq", f"root={root!r}")
except Exception as exc:
    check(False, "scipy.optimize.brentq", f"{type(exc).__name__}: {exc}")

try:
    from SALib.analyze import delta as _d, pawn as _p   # noqa: F401
    from SALib.sample import latin as _l                # noqa: F401
    check(True, "SALib sample/analyze imports")
except Exception as exc:
    check(False, "SALib sample/analyze imports", f"{type(exc).__name__}: {exc}")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import tempfile, os
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.loglog([1e2, 1e5], [1e-12, 1e-7])
    ax.add_patch(Rectangle((1e3, 1e-11), 1e3, 1e-10, alpha=0.2))
    out = os.path.join(tempfile.gettempdir(), "_verify.png")
    fig.savefig(out, dpi=100, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    size = os.path.getsize(out)
    os.unlink(out)
    check(size > 1000, "matplotlib headless savefig",
          f"backend={matplotlib.get_backend()}, {size} bytes")
except Exception as exc:
    check(False, "matplotlib headless savefig", f"{type(exc).__name__}: {exc}")

try:
    import plotly.graph_objects as go
    html = go.Figure(go.Parcoords(dimensions=[dict(label="a", values=[1, 2, 3])])
                     ).to_html(include_plotlyjs="cdn")
    check(len(html) > 500, "plotly Parcoords to_html", f"{len(html)} chars")
except Exception as exc:
    check(False, "plotly Parcoords to_html", f"{type(exc).__name__}: {exc}")

try:
    from IPython.display import display, HTML   # noqa: F401
    check(True, "IPython.display import")
except Exception as exc:
    check(False, "IPython.display import", f"{type(exc).__name__}: {exc}")

try:
    import io, pandas as pd
    df = pd.DataFrame({"flux": [1.888003770552079e-07, 3.1e-12]})
    buf = io.StringIO(); df.to_csv(buf, index=False); buf.seek(0)
    check(bool((pd.read_csv(buf)["flux"] == df["flux"]).all()),
          "pandas CSV round-trip is exact")
except Exception as exc:
    check(False, "pandas CSV round-trip is exact", f"{type(exc).__name__}: {exc}")

# GUI modules must be absent — an accidental `import turtle` would otherwise
# only fail once the image is in the company's hands. Mirrors ruff's TID251.
for mod in ("turtle", "tkinter"):
    try:
        __import__(mod)
        check(False, f"{mod} absent (headless image)", "unexpectedly importable")
    except ImportError:
        check(True, f"{mod} absent (headless image)")

print("=" * 78)
if failures:
    print(f"ENVIRONMENT VERIFICATION FAILED — {len(failures)} check(s):")
    for f in failures:
        print(f"  - {f}")
    sys.exit(1)
print("ENVIRONMENT VERIFICATION PASSED")
