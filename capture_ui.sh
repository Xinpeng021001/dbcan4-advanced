#!/usr/bin/env bash
# capture_ui.sh — true-browser screenshots of the live BioForge web UI.
# Uses headless Google Chrome (real CSS + JS + the 3Dmol WebGL viewer), NOT
# WeasyPrint, so the capture matches what a user actually sees in a browser.
#
# The 3Dmol structure viewer needs WebGL; headless Chrome has no GPU, so we run
# ANGLE's SwiftShader software-WebGL backend and give the async fetch+render a
# generous virtual-time budget. A tall window keeps the structure card in frame.
#
# Usage: capture_ui.sh <base_url> <out_dir> [gene_ids...]
set -euo pipefail
BASE="${1:?base url e.g. http://127.0.0.1:8808}"
OUT="${2:?output dir}"
shift 2 || true
GENES=("$@"); [ ${#GENES[@]} -eq 0 ] && GENES=(1 2 3)
mkdir -p "$OUT"
CHROME=$(command -v google-chrome-stable || command -v google-chrome || command -v chromium)

# Flat pages (dashboard/browse) — no WebGL needed.
shot() { # url file
  timeout 90 "$CHROME" --headless=new --no-sandbox --hide-scrollbars \
    --window-size=1440,2000 --force-device-scale-factor=2 \
    --virtual-time-budget=9000 --screenshot="$2" "$1" 2>/dev/null || true
  [ -s "$2" ] && echo "  captured $(basename "$2") ($(wc -c <"$2") B)" || echo "  FAILED $1"
}

# Gene pages — enable software WebGL so the 3Dmol structure viewer renders.
# Warm the PDB URL first so 3Dmol's fetch() hits a warm cache within the budget.
shot_gene() { # url file
  local pdb
  pdb=$(curl -s "$1" | grep -oE 'data-pdb="[^"]+"' | head -1 | sed 's/data-pdb="//;s/"//' || true)
  [ -n "$pdb" ] && curl -s "$BASE$pdb" >/dev/null || true
  timeout 120 "$CHROME" --headless=new --no-sandbox --in-process-gpu \
    --use-gl=angle --use-angle=swiftshader-webgl --enable-unsafe-swiftshader \
    --ignore-gpu-blocklist --hide-scrollbars \
    --window-size=1400,4600 --force-device-scale-factor=1 \
    --run-all-compositor-stages-before-draw --virtual-time-budget=30000 \
    --screenshot="$2" "$1" 2>/dev/null || true
  [ -s "$2" ] && echo "  captured $(basename "$2") ($(wc -c <"$2") B)" || echo "  FAILED $1"
}

echo "### capturing $BASE -> $OUT"
shot      "$BASE/"       "$OUT/01_dashboard.png"
shot      "$BASE/browse" "$OUT/02_browse.png"
for g in "${GENES[@]}"; do shot_gene "$BASE/genes/$g" "$OUT/gene_$g.png"; done
echo "### done"
