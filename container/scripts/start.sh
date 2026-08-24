#!/bin/bash
# Start JupyterLab against the workspace. Ctrl-C to stop.
set -euo pipefail

IMAGE="${IMAGE:-hydrogen-model:1.0.0}"
PORT="${PORT:-8888}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$HERE/.." && pwd)"

if [ ! -d "$PKG/workspace" ]; then
  echo "No workspace/ directory next to scripts/ — is the package complete?" >&2
  exit 1
fi

cat <<BANNER

  JupyterLab is starting.

  Copy the FULL URL printed below, including the ?token=... part, into your
  browser. Without the token the page will refuse to load.

      http://127.0.0.1:$PORT/lab?token=<shown below>

  Your files are in workspace/. Results are written next to the notebook that
  produces them. Press Ctrl-C here to stop the server.

BANNER

# --user "$(id -u):0" makes files written into workspace/ belong to you rather
# than to root. Group 0 is required: the image's HOME is group-writable so any
# host UID gets a writable home, which JupyterLab needs.
exec docker run --rm -it \
  --name hydrogen-jupyter \
  -p "127.0.0.1:$PORT:8888" \
  -v "$PKG/workspace:/workspace" \
  --user "$(id -u):0" \
  "$IMAGE"
