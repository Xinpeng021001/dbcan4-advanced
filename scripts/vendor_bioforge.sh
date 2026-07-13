#!/usr/bin/env bash
# =============================================================================
# sync_bioforge.sh — refresh the src/bioforge web layer from a working copy
# =============================================================================
# The `bioforge` web UI + ingest package lives in this repo (src/bioforge) so the
# whole product runs from a single `git clone`. If you maintain the web layer in a
# separate working copy, run this to re-sync the package + migrations here.
#
# Usage:
#     scripts/vendor_bioforge.sh [BIOFORGE_WORKDIR]
#   BIOFORGE_WORKDIR defaults to $BIOFORGE_SRC or /array1/xinpeng/biodb.
# =============================================================================
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="${1:-${BIOFORGE_SRC:-/array1/xinpeng/biodb}}"

[ -d "$SRC/src/bioforge" ] || { echo "ERROR: no bioforge package at $SRC/src/bioforge" >&2; exit 2; }

echo "Syncing bioforge from: $SRC"
echo "                   to: $REPO"

# 1. package
rm -rf "$REPO/src/bioforge"
rsync -a --exclude='__pycache__' --exclude='*.pyc' --exclude='*.egg-info' \
      "$SRC/src/bioforge/" "$REPO/src/bioforge/"

# 2. alembic migrations + config (they sit at the repo root, not inside the package)
rm -rf "$REPO/db"
rsync -a --exclude='__pycache__' --exclude='*.pyc' "$SRC/db/" "$REPO/db/"
cp "$SRC/alembic.ini" "$REPO/alembic.ini"

echo "Done. Reinstall so console scripts + package-data update:  pip install -e ."
