#!/usr/bin/env bash
# Launch N one-process-per-GPU ESM-C embedding shards over a FASTA, wait for all,
# and write a DONE (or FAILED) marker. Designed to run inside a screen session so
# the user can `screen -r esmc_embed` to watch, and a separate watcher can poll the marker.
#
# Usage: embed_all_gpus.sh <fasta> <out_prefix> <ngpu> <logdir> [model]
set -uo pipefail
FASTA=$1; OUT=$2; NGPU=${3:-8}; LOGDIR=${4:-.}; MODEL=${5:-esmc_600m}
PY=/array1/xinpeng/dbcan4-advanced/venv/bin/python
export HF_HOME=/array1/xinpeng/dbcan4-advanced/hf_cache
export TOKENIZERS_PARALLELISM=false
mkdir -p "$LOGDIR"
rm -f "${OUT}.DONE" "${OUT}.FAILED"
echo "[wrapper] start $(date) fasta=$FASTA out=$OUT ngpu=$NGPU" > "$LOGDIR/wrapper.log"
pids=()
for ((s=0; s<NGPU; s++)); do
  CUDA_VISIBLE_DEVICES=$s $PY embed_esmc.py --fasta "$FASTA" --out-prefix "$OUT" \
       --shard $s --nshards $NGPU --model "$MODEL" > "$LOGDIR/embed.shard$s.log" 2>&1 &
  pids+=($!)
  echo "[wrapper] launched shard $s pid ${pids[-1]} on GPU $s" >> "$LOGDIR/wrapper.log"
done
fail=0
for i in "${!pids[@]}"; do
  if ! wait "${pids[$i]}"; then
    echo "[wrapper] shard $i FAILED (pid ${pids[$i]})" >> "$LOGDIR/wrapper.log"
    fail=1
  else
    echo "[wrapper] shard $i done" >> "$LOGDIR/wrapper.log"
  fi
done
if [ $fail -eq 0 ]; then
  echo "[wrapper] ALL DONE $(date)" >> "$LOGDIR/wrapper.log"
  touch "${OUT}.DONE"
else
  echo "[wrapper] FAILED $(date)" >> "$LOGDIR/wrapper.log"
  touch "${OUT}.FAILED"
fi
