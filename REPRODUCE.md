# Reproducibility — exact commands & parameters

All jobs run **directly on `met`** (no SLURM). Long jobs run inside a named GNU `screen`
session (`screen -dmS dbcan4_<step> …`; attach with `screen -r dbcan4_<step>`).

- Project root on met: `/array1/xinpeng/dbcan4-advanced` (`$BASE`)
- Python venv: `source $BASE/venv/bin/activate` (py3.11, torch 2.12.1+cu130, esm 3.2.1 ESM-C, faiss-cpu, scikit-learn)
- Data: `$BASE/data/` — `reference_2024.faa` (337,759 fungal CAZymes), `eval_2025.faa` (4,726),
  `eval_2025_labels.tsv`, `reference_labels_2024.tsv`
- dbCAN dev: 5.0.7.dev50+g6250e4e79 (swapped into venv site-packages over released 5.2.9)
- dbCAN DB: `/array1/xinpeng/dbcan_db` (downloaded via `run_dbcan database`, current release)

---

## Step 4 — build reference + evaluation splits
```bash
python scripts/build_reference.py \
  --cazy-2024 /array1/xinpeng/fungi-cazyme-project/CAZyDB.07142024.fungi.faa \
  --cazy-2025 /array1/xinpeng/fungi-cazyme-project/CAZyDB.07242025.fungi.faa \
  --outdir $BASE/data --eval-per-bucket 4000 --min-len 30 --max-len 1500
# novelty = exact sequence-MD5 identity vs 2024. Buckets: carried_over / novel_seq / novel_family.
```

## Step 5 — baselines

### (a) dbCAN dev, CURRENT DB — canonical baseline (HMMER + dbCAN-sub + DIAMOND in one call)
This is the flagship baseline. `run_dbcan CAZyme_annotation` runs all three tiers and writes
`overview.tsv` (columns: Gene ID, EC#, **dbCAN_hmm**, dbCAN_sub, DIAMOND, #ofTools, Recommend Results).
Do NOT hand-roll hmmscan/diamond — this single command is the reference.
```bash
run_dbcan CAZyme_annotation \
  --mode protein --input_raw_data $BASE/data/eval_2025.faa \
  --db_dir /array1/xinpeng/dbcan_db --output_dir $BASE/results/dbcan_eval2025 \
  --methods diamond,hmm,dbCANsub --threads 32
```
**Exact cutoffs used (dbCAN dev 5.0.7.dev50 defaults, read from `dbcan/parameter.py`):**

| tier | parameter | default |
|---|---|---|
| DIAMOND (CAZy) | E-value | **1e-102** |
| dbCAN HMM | E-value / coverage | **1e-15** / **0.35** |
| dbCAN-sub HMM | E-value / coverage | **1e-15** / **0.35** |
| (TF DIAMOND) | E-value / query-cover | 1e-4 / 35% |
| (TC/TF HMM) | E-value / coverage | 1e-4 / 0.35 |

Scored with `norm_fams()` (strips `(coords)` and `_eNNN` subfamily suffixes, `+`-joins) →
`dbcan_eval2025_scored.json`. dbCAN found ≥1 CAZyme for 4,680 / 4,726 eval proteins.
> **DB vintage caveat:** this DB is the *current* CAZy release, which already contains all 20
> new-2025 families (all 20 confirmed as NAME entries in dbCAN.hmm). So the high dbCAN recall on
> the "novel_family" bucket is expected and is **not** a temporal generalization result — it
> measures "what today's dbCAN finds," not "what dbCAN trained on 2024 would find."

### (b) DIAMOND temporal baseline — 2024 fungal reference (identical set the pLM methods retrieve over)
```bash
diamond makedb --in $BASE/data/reference_2024.faa -d $BASE/db/ref2024_fungi -p 32
diamond blastp -q $BASE/data/eval_2025.faa -d $BASE/db/ref2024_fungi \
  -o $BASE/results/diamond_fungiref_hits.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --max-target-seqs 5 --evalue 1e-3 -p 32 --quiet
python scripts/diamond_baseline.py \
  --hits $BASE/results/diamond_fungiref_hits.tsv --labels $BASE/data/eval_2025_labels.tsv \
  --out-pred $BASE/results/diamond_fungiref_pred.tsv \
  --out-summary $BASE/results/diamond_fungiref_summary.json --evalues 1e-15,1e-10,1e-5,1e-3
# best-hit family = top bitscore among hits with evalue <= cut.
```
> Note: an earlier DIAMOND run used the all-kingdom `CAZyDB.07142024.fa` (~1.8 GB); results were
> within 0.3% of the fungal-reference run, but the fungal-reference version above is the fair,
> apples-to-apples comparison against ESM-C (identical reference set).
> **Cutoff note:** this *standalone* DIAMOND baseline uses `--evalue 1e-3` (permissive, to measure
> raw retrieval ceiling by best-bitscore hit). dbCAN's production DIAMOND tier uses a far stricter
> `1e-102` (see table above). The two answer different questions: standalone = "is there any 2024
> fungal homolog?"; dbCAN = "is there a high-confidence CAZy assignment?".

### (c) Temporal HMMER (2024-only) — SUPPLEMENTARY fair-temporal check (not canonical dbCAN)
Because dbCAN's DB is the current release, it cannot answer "what would an HMM tier trained on 2024
do on 2025?". This one-off builds per-base-family profile HMMs from the 2024 reference only:
```bash
python scripts/hmm_baseline_2024.py \
  --ref data/reference_2024.faa --eval data/eval_2025.faa \
  --workdir hmm2024 --out-pred results/hmm2024_pred.tsv \
  --min-seqs 5 --max-seqs 300 --threads 48 --seed 0
# per family >=5 seqs (223 of 266 base families): mafft --auto -> hmmbuild --amino;
# cat -> hmmpress -> hmmscan -E 1e-3; best full-seq e-value per query -> predicted base family.
```
Result (parent-family recall): novel_seq 0.876, novel_subfamily 0.949, **truly-novel base 0.000** —
same blind spot as DIAMOND and the pLM methods on genuinely new families.

## Step 6 — ESM-C embeddings + retrieval
```bash
# embed (8-GPU sharded, ESM-C 600M, mean-pooled 1152-dim)
for s in $(seq 0 7); do
  CUDA_VISIBLE_DEVICES=$s python scripts/embed_esmc.py \
    --fasta data/reference_2024.faa --out-prefix emb/ref2024 --shard $s --nshards 8 --max-len 1500 &
done; wait
# (same for eval_2025 -> emb/eval2025)

# retrieval: kNN (k=15 majority vote) + nearest-centroid (403 family prototypes)
python scripts/retrieval_esmc.py \
  --ref-prefix emb/ref2024 --eval-prefix emb/eval2025 --labels data/eval_2025_labels.tsv \
  --out-summary results/esmc_retrieval_summary.json --out-pred results/esmc_retrieval_pred.tsv \
  --k 15 --op-threshold 0.5
```

## Step 7 — trained heads (on FROZEN ESM-C embeddings)
```bash
CUDA_VISIBLE_DEVICES=0 python scripts/train_heads.py \
  --ref-prefix emb/ref2024 --eval-prefix emb/eval2025 --labels data/eval_2025_labels.tsv \
  --outdir results/heads --epochs 30 --proj-dim 256 --hidden 1024 \
  --batch 4096 --lr 1e-3 --temp 0.1 --k 15 --min-count 2 --seed 0
# SupCon projection head + softmax classifier; val classifier acc ~0.966 in ~30 s on one A5500.
```

## Analysis
```bash
# dual-granularity scoring (exact subfamily vs parent family) + truly-novel-base split, run locally
# from downloaded pred TSVs; see scripts/*.py and results/unified_scores.json.
```

### AUROC note
Novelty-detection AUROC is computed **tie-aware** (average ranks). Vote-purity is discretized to
k+1 values, so a naive argsort-rank AUROC over-/under-counts ties — always use the average-rank form
(fixed in `retrieval_esmc.py` / `train_heads.py`).