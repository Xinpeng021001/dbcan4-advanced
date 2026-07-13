#!/usr/bin/env bash
# =============================================================================
# vendor_bioforge.sh — refresh the vendored BioForge web layer from upstream biodb
# =============================================================================
# The `bioforge` web UI + ingest package is VENDORED (copied) into this repo so
# the whole product runs from a single `git clone` (see src/bioforge/VENDORED.md).
# Run this to re-sync it from a local biodb checkout after upstream changes.
#
# Usage:
#     scripts/vendor_bioforge.sh [BIODB_CHECKOUT]
#   BIODB_CHECKOUT defaults to $BIODB_UPSTREAM or ../../biodb.
# =============================================================================
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="${1:-${BIODB_UPSTREAM:-/array1/xinpeng/biodb}}"

[ -d "$SRC/src/bioforge" ] || { echo "ERROR: no bioforge package at $SRC/src/bioforge" >&2; exit 2; }

echo "Vendoring bioforge from: $SRC"
echo "                    to: $REPO"

# 1. package
rm -rf "$REPO/src/bioforge"
rsync -a --exclude='__pycache__' --exclude='*.pyc' --exclude='*.egg-info' \
      "$SRC/src/bioforge/" "$REPO/src/bioforge/"

# 2. alembic migrations + config (live OUTSIDE the package upstream)
rm -rf "$REPO/db"
rsync -a --exclude='__pycache__' --exclude='*.pyc' "$SRC/db/" "$REPO/db/"
cp "$SRC/alembic.ini" "$REPO/alembic.ini"

# 3. provenance stamp
COMMIT="$(git -C "$SRC" rev-parse --short HEAD 2>/dev/null || echo unknown)"
SUBJECT="$(git -C "$SRC" log -1 --format='%s' 2>/dev/null || echo '')"
cat > "$REPO/src/bioforge/VENDORED.md" <<EOF
# Vendored BioForge web layer

This \`bioforge\` package is **vendored** (copied) into dbCAN4-advanced so the whole
product — pipeline **and** the web UI/ingest layer — runs from a single
\`git clone\` with no second repository to fetch.

- Upstream: https://github.com/Xinpeng021001/biodb  (branch \`feature/advanced-cazyme-integration\`)
- Snapshot: commit \`$COMMIT\` — "$SUBJECT"

Refresh with: \`scripts/vendor_bioforge.sh [BIODB_CHECKOUT]\`
EOF

echo "Done. Vendored commit $COMMIT."
echo "Reinstall so console scripts + package-data update:  pip install -e ."
