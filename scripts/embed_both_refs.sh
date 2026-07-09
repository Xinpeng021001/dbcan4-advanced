#!/usr/bin/env bash
# Driver: embed BOTH clustered references (all-kingdom then fungal) across 8 GPUs,
# writing shard npz to the stable emb_refscope dir. Meant to run inside a detached
# screen session (screen -r esmc_embed to watch). Touches ALL.DONE / ALL.FAILED.
set -uo pipefail
REF=/array1/xinpeng/dbcan4-advanced/data_refscope
EMB=/array1/xinpeng/dbcan4-advanced/emb_refscope
SCR=/array1/xinpeng/dbcan4-advanced/scripts_refscope
mkdir -p "$EMB/logs"
cd "$SCR"
rm -f "$EMB/ALL.DONE" "$EMB/ALL.FAILED"
echo "[driver] start $(date)"

bash embed_all_gpus.sh "$REF/allking_2024_c50_rep_seq.fasta" "$EMB/allking_c50" 8 "$EMB/logs/ak" esmc_600m
bash embed_all_gpus.sh "$REF/fungi_2024_c50_rep_seq.fasta"   "$EMB/fungi_c50"   8 "$EMB/logs/fu" esmc_600m

if [ -f "$EMB/allking_c50.DONE" ] && [ -f "$EMB/fungi_c50.DONE" ]; then
  echo "[driver] BOTH DONE $(date)"; touch "$EMB/ALL.DONE"
else
  echo "[driver] one or both FAILED $(date)"; touch "$EMB/ALL.FAILED"
fi
