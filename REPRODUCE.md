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

## Step 8 — reference scope: fungi-specific vs all-kingdom (MMseqs-clustered)

Tests whether the ESM-C reference should be fungi-only or include all kingdoms, and
whether MMseqs redundancy reduction is safe. Design **separates two effects** by
holding the eval_2025 set + truth fixed and swapping only the reference:
redundancy effect = unclustered→clustered *fungal*; scope effect = clustered
*fungal*→clustered *all-kingdom* (redundancy matched). All on `met`, in a screen
session. Outputs live in `$BASE/data_refscope/`, `$BASE/emb_refscope/`,
`$BASE/results_refscope/`.

### (a) build both labeled references (2024 cutoff + eval-MD5 leak guard)
`build_labeled_ref.py` keeps records with ≥1 CAZy family, drops exact-dup
sequences, drops any sequence whose MD5 matches an eval_2025 sequence (temporal /
self-match guard), and length-filters 30–1500.
```bash
# all-kingdom (2024)
python scripts/build_labeled_ref.py \
  --in-faa /array1/xinpeng/fungi-cazyme-project/CAZyDB.07142024.fa \
  --exclude-faa $BASE/data/eval_2025.faa \
  --out-prefix $BASE/data_refscope/allking_2024 --min-len 30 --max-len 1500
# fungal (2024)
python scripts/build_labeled_ref.py \
  --in-faa /array1/xinpeng/fungi-cazyme-project/CAZyDB.07142024.fungi.faa \
  --exclude-faa $BASE/data/eval_2025.faa \
  --out-prefix $BASE/data_refscope/fungi_2024 --min-len 30 --max-len 1500
# scale: all-kingdom 3,613,700 labeled -> 2,150,909 kept (46 excluded by eval-MD5 guard), 820 families
#        fungal       457,694  labeled ->   398,271 kept (0 eval-MD5 hits),               421 families
```

### (b) MMseqs redundancy reduction (50% identity, 80% coverage)
`cluster_ref_mmseqs.sh <in.faa> <out_prefix> <tmp_dir> <min_seq_id> <cov> <threads>`
wraps `mmseqs easy-linclust`; `*_rep_seq.fasta` keeps `ID|FAM` headers.
```bash
scripts/cluster_ref_mmseqs.sh $BASE/data_refscope/allking_2024.faa \
  $BASE/data_refscope/allking_2024_c50 /array1/xinpeng/scratch/mm_ak 0.5 0.8 64
scripts/cluster_ref_mmseqs.sh $BASE/data_refscope/fungi_2024.faa \
  $BASE/data_refscope/fungi_2024_c50   /array1/xinpeng/scratch/mm_fu 0.5 0.8 64
# reps: all-kingdom 2,150,909 -> 465,117 (4.6x); fungal 398,271 -> 110,299 (3.6x)
```
> **MMseqs gotcha:** `easy-linclust` writes `*_rep_seq.fasta` headers with a
> **trailing space** (`>ID|GT4 `). `embed_esmc.py` `.strip()`s the header so the
> family label is not corrupted to `'GT4 '` — an un-stripped label silently fails
> every exact match and produces near-zero scores.

### (c) embed both clustered references (8-GPU sharded, from a screen session)
`embed_all_gpus.sh <fasta> <out_prefix> <ngpu> <logdir> [model]` launches one
process per GPU and writes a `<out_prefix>.DONE` marker; `embed_both_refs.sh` drives
both references and touches `emb_refscope/ALL.DONE` when finished.
```bash
mkdir -p $BASE/scripts_refscope && cp scripts/embed_esmc.py scripts/embed_all_gpus.sh \
  scripts/embed_both_refs.sh $BASE/scripts_refscope/
screen -dmS esmc_embed bash $BASE/scripts_refscope/embed_both_refs.sh   # attach: screen -r esmc_embed
# ESM-C 600M, mean-pooled 1152-dim, ~27 seq/s/GPU, ~49 min wall for 575,416 reps on 8x A5500
# integrity: allking_c50 = 465,117 rows / 8 shards; fungi_c50 = 110,299 / 8; dim 1152
```

### (d) retrieval + trained heads vs the FIXED eval_2025 embeddings
Reuses `emb/eval2025.shard*.npz` from Step 6 (eval set never changes).
```bash
for R in allking fungi; do
  python scripts/retrieval_esmc.py \
    --ref-prefix $BASE/emb_refscope/${R}_c50 --eval-prefix $BASE/emb/eval2025 \
    --labels $BASE/data/eval_2025_labels.tsv --k 15 \
    --out-summary $BASE/results_refscope/retr_${R}_c50_summary.json \
    --out-pred    $BASE/results_refscope/retr_${R}_c50_pred.tsv
  CUDA_VISIBLE_DEVICES=0 python scripts/train_heads.py \
    --ref-prefix $BASE/emb_refscope/${R}_c50 --eval-prefix $BASE/emb/eval2025 \
    --labels $BASE/data/eval_2025_labels.tsv \
    --outdir $BASE/results_refscope/heads_${R} --epochs 30
done
```
> `train_heads.py` batches the reference projection and the eval kNN (query blocks of
> 512), so peak memory is independent of reference size — required for the 465K
> all-kingdom reference; numerically identical to the unbatched form.

### (e) three-level scoring + effect decomposition (run locally from pred TSVs)
Score class / family / subfamily (multidomain via set operations) for each of
5 methods × {unclustered-fungal, clustered-fungal, clustered-all-kingdom}. The
unclustered-fungal predictions are the Step 6/7 `esmc_retrieval_pred.tsv` /
`head_eval_pred.tsv` (redundancy baseline). Effect sizes:
`redundancy = clustered_fungal − unclustered_fungal`, `scope = clustered_allking −
clustered_fungal`. Outputs `benchmarks/refscope_threelevel_all.tsv` (270 rows) and
`benchmarks/refscope_effects.tsv`; figure `docs/figures/refscope_effect.png`.

**Headline (known families, novel-seq subfamily overlap):** contrastive kNN
0.973→0.970→0.970 (redundancy −0.003, scope 0.000). Redundancy reduction is nearly
free for trained retrieval; all-kingdom does not help fungal prediction; the raw
ESM-C centroid is the one method that degrades under clustering (0.497→0.338). See
report §4.7.

## Analysis
```bash
# dual-granularity scoring (exact subfamily vs parent family) + truly-novel-base split, run locally
# from downloaded pred TSVs; see scripts/*.py and results/unified_scores.json.
```

### AUROC note
Novelty-detection AUROC is computed **tie-aware** (average ranks). Vote-purity is discretized to
k+1 values, so a naive argsort-rank AUROC over-/under-counts ties — always use the average-rank form
(fixed in `retrieval_esmc.py` / `train_heads.py`).