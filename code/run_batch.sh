#!/bin/zsh
# Detached batch launcher. caffeinate dies with this script and holds off idle
# sleep for the whole run; python replaces this shell (same PID) so caffeinate -w
# keeps waiting on it until the batch exits.
set -e
REPO="$(cd "$(dirname "$0")/.." && pwd)"          # repo root = parent of this script's dir
ENV="${MORPH_ENV:-/Users/cafierom/miniforge3/envs/openmm-md}"
cd "$REPO"
caffeinate -w $$ &
exec "$ENV/bin/python" code/run_batch.py
