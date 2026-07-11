#!/usr/bin/env bash
# Run CLEAN maxsep EC inference on a protein FASTA, emitting a CLEAN maxsep CSV.
# Usage: run_clean.sh <input.faa> <out.csv> [clean_dir] [clean_venv]
# Encodes the three real gotchas discovered on met:
#   (1) CLEAN_infer_fasta.py -d ARG prepends 'inputs/' internally -> pass bare name.
#   (2) ESM-1b alphabet throws KeyError '*' on trailing stop codons -> strip '*'.
#   (3) sequences >1022 aa exceed ESM-1b -> --truncation_seq_length 1022.
set -euo pipefail
IN="$1"; OUT="$2"
CLEAN_DIR="${3:-/array1/xinpeng/scratch/CLEAN}"
CLEAN_VENV="${4:-/array1/xinpeng/scratch/venv_clean}"
APP="$CLEAN_DIR/app"
NAME="clean_$$"                      # unique per-run basename

# stage cleaned FASTA (strip '*' and whitespace from sequence lines)
mkdir -p "$APP/data/inputs" "$APP/results/inputs" "$APP/data/esm_data"
awk '/^>/{print; next}{gsub(/[*[:space:]]/,""); print}' "$IN" > "$APP/data/inputs/$NAME.fasta"

pushd "$APP" >/dev/null
# STEP 1: extract ESM-1b mean embeddings (truncate long seqs)
"$CLEAN_VENV/bin/python" esm/scripts/extract.py esm1b_t33_650M_UR50S \
    "data/inputs/$NAME.fasta" data/esm_data --include mean --truncation_seq_length 1022
# STEP 2: maxsep inference (script prepends 'inputs/')
"$CLEAN_VENV/bin/python" CLEAN_infer_fasta.py -d "$NAME"
popd >/dev/null

cp "$APP/results/inputs/${NAME}_maxsep.csv" "$OUT"
echo "[run_clean] wrote $OUT"
