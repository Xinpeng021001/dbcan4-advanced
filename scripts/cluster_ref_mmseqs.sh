#!/usr/bin/env bash
# Redundancy-reduce a labeled reference FASTA with MMseqs2 linclust.
# Representative sequences keep their original ">ID|FAM[,FAM...]" headers, so the
# output *_rep_seq.fasta is directly consumable by embed_esmc.py.
#
# Usage: cluster_ref_mmseqs.sh <in.faa> <out_prefix> <tmp_dir> <min_seq_id> <cov> <threads>
# Both references are clustered with IDENTICAL parameters so the fungi-vs-all-kingdom
# comparison reflects taxonomic scope, not differing redundancy.
set -eo pipefail
IN=$1; OUT=$2; TMP=$3; SID=${4:-0.5}; COV=${5:-0.8}; THREADS=${6:-64}
mkdir -p "$TMP"
mmseqs easy-linclust "$IN" "$OUT" "$TMP" \
    --min-seq-id "$SID" -c "$COV" --cov-mode 0 \
    --threads "$THREADS" --split-memory-limit 200G -v 2
echo "reps: $(grep -c '^>' ${OUT}_rep_seq.fasta)"
