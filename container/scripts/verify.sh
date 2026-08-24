#!/bin/bash
# Both acceptance gates, with networking disabled. Run this first, and after any
# transfer of the image archive.
set -uo pipefail

IMAGE="${IMAGE:-hydrogen-model:1.0.1}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$HERE/.." && pwd)"

# Check the package layout before mounting. Docker silently CREATES a missing
# bind-mount source as an empty directory, which would surface as a confusing
# gate-2 import failure rather than "your package is incomplete".
if [ ! -d "$PKG/workspace/calculations" ]; then
  echo "Package looks incomplete: no workspace/calculations/ next to scripts/." >&2
  echo "  expected: $PKG/workspace/calculations" >&2
  echo "Run this script from inside the delivered package directory." >&2
  exit 1
fi

fail=0

echo "### Gate 1/2 — environment (offline)"
if ! docker run --rm --network none "$IMAGE" python /opt/verify_environment.py; then
  echo ">>> GATE 1 FAILED: the image does not contain the validated environment." >&2
  fail=1
fi

echo
echo "### Gate 2/2 — scientific smoke test (offline)"
if ! docker run --rm --network none \
      -v "$PKG/workspace:/workspace" -w /workspace \
      --user "$(id -u):$(id -g)" \
      "$IMAGE" python /opt/model_smoke_test.py; then
  echo ">>> GATE 2 FAILED: the model does not reproduce its reference values." >&2
  fail=1
fi

echo
if [ "$fail" -ne 0 ]; then
  echo "VERIFICATION FAILED — do not use this image for analysis." >&2
  exit 1
fi
echo "VERIFICATION PASSED — $IMAGE is the validated environment."
