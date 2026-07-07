# Daily log — dbCAN4-advanced

One-week prototype sprint. One dated entry per working day: progress, blockers, next.

---

## Day 1 — 2026-07-07

**Done**
- Method survey completed and grounded in retrieved full text / records: CAZyLingua
  (ProtT5 + two-stage RF+FFN), CLEAN (contrastive, *Science* 2023), Foldseek (3Di,
  *Nat Biotech* 2023), CAZyme3D (870,740 AF structures; ID50 = 188,574),
  SaProt (structure-aware vocab, ICLR 2024), GraphEC (ProtT5 + ESMFold geometric GNN),
  ProstT5 (sequence→3Di, enabling folding-free structure search), EZSpecificity (*Nature* 2025).
- Wrote **design document** `docs/design_dbcan4_advanced.md`: motivation, full method
  survey, two-tier architecture (ESM-C sequence tier + Foldseek/SaProt structure tier +
  fusion), temporal-holdout evaluation design, one-week map, BioForge integration, risks.
- Rendered publication-style **architecture diagram** `docs/architecture.png`.
- **Retargeted** compute from leu → **met** (128 CPU, 8× A5500, 386 GB RAM).
- Created **GitHub repo** `github.com/Xinpeng021001/dbcan4-advanced` (public) and scaffolded it.

**Decisions**
- **ESM-C** (not ESM-2) as the primary sequence pLM — 600M default on 24 GB GPUs.
- Retrieval is not limited to kNN: benchmark kNN vs nearest-centroid vs contrastive vs MLP.
- **SaProt** added as a structure-aware-embedding ablation arm alongside Foldseek.
- **Primary evaluation = 2024→2025 temporal holdout** (novel-sequence + novel-family subsets).
- EZSpecificity-style substrate prediction scoped as **Phase 2 roadmap**, not week 1.

**Blockers**
- DEFT could not be pinned to a citable record; fusion described from verified analogs
  (CLEAN-Contact, Phold, CAZyLingua's Foldseek validation) and implemented as our own rule.

**Next (Day 2)**
- Port envs + CAZy DB from leu to met; smoke-test run_dbcan / Foldseek / ESMFold.
- Build 2024/2025 reference+test splits; run dbCAN baseline; isolate the missed set.

---

## Day 2 — 2026-07-07

**Done**
- **Ported compute to met** and built the project GPU env at `/array1/xinpeng/dbcan4-advanced/venv`:
  torch 2.12.1 (CUDA, 8 GPUs), EvolutionaryScale `esm` 3.2.1 (ESM-C `esmc_600m` → 1152-dim),
  faiss/sklearn/biopython. Verified ESM-C loads on GPU and embeds. See `docs/env_setup.md`.
- Confirmed submission model: met has **no scheduler** — jobs run directly (no SLURM).
- **Built labeled reference + temporal-holdout splits** (`scripts/build_reference.py`):
  - Reference (2024): 337,759 labeled fungal CAZymes, 411 families.
  - 2025 novelty by **exact-sequence identity** vs 2024: 329,724 carried-over,
    **51,704 novel-sequence**, **728 novel-family** (20 families new in 2025).
  - Bounded eval set `eval_2025.faa`: 4,726 seqs (4,000 novel-seq + 726 novel-family).
- Corrected design-doc reference provenance (separated retrieved-full-text vs
  title/record-only citations); retrieved CLEAN-Contact full text; confirmed EZSpecificity DOI.

**Decisions**
- Novelty defined by **exact-sequence md5**, not accession — JGI numeric IDs collide across
  genomes (407,338 unique seqs vs 337,759 unique pids in 2024), so accession matching is unreliable.

**Next (Day 3)**
- Run dbCAN baseline (HMMER + DIAMOND) on eval_2025 → the "missed by sequence-similarity" subset.
- Compute ESM-C embeddings for reference + eval; kNN / centroid retrieval baselines.

---

## Day 2 (cont.) — baselines (Step 5)

**Done**
- Installed the **dev dbCAN** (5.0.7.dev50+g6250e4e79, from leu's run_dbcan_new env) into the
  project venv on met, over released 5.2.9. Downloaded the dbCAN database (7.4 GB) into
  `/array1/xinpeng/dbcan_db` (dbCAN.hmm, CAZy.dmnd, dbCAN-sub.hmm, +TF/TCDB/STP; STP.hmm needed a
  direct-URL fetch after a 403).
- **Two baselines on the 4,726-seq eval set:**
  1. **Temporally-clean sequence similarity** — DIAMOND vs the **2024** CAZy DB (best-hit family):
     - novel-sequence (known family): **96.5% exact family**
     - novel-family (new in 2025): **0.1% exact family** — the sequence-similarity blind spot.
  2. **Production dbCAN, current DB** (dev, hmm+dbCAN_sub+DIAMOND): novel-family 97% exact — but
     only because the **current DB already contains all 20 new-2025 families** (all 20 confirmed
     present as NAME entries in dbCAN.hmm). This is a "what today's dbCAN finds" reference, **not** a temporal test.

**Takeaway**
- The fair 2024->2025 comparison is baseline #1. Its collapse on novel families (0.1%) is exactly
  the gap ESM-C retrieval + structure similarity must close. Figure: `docs/figures/baseline_novelty_recall.png`.

**Next**
- ESM-C embeddings for reference (2024) + eval; kNN and nearest-centroid retrieval; measure
  novel-family recall vs the 0.1% sequence-similarity floor.

---

## Day 2 (cont.) — ESM-C retrieval baselines (Step 6)

**Done**
- Embedded reference_2024 (337,759 seqs) + eval_2025 (4,726) with **ESM-C 600M**, mean-pooled
  1152-dim, sharded across all 8 GPUs (scripts/embed_esmc.py). Embeddings cached on met at
  `/array1/xinpeng/dbcan4-advanced/emb/`.
- Two retrieval schemes (scripts/retrieval_esmc.py): **kNN** (k=15 majority vote) and
  **nearest-centroid** (CLEAN-style prototype, one L2-normalized mean per family over 403 families).

**Results (honest, and design-shaping)**
- Known families (novel-seq, n=4,000): ESM-C **kNN 85.8% exact**, centroid 45.6% — both *below*
  DIAMOND-2024 (96.5%). Off-the-shelf ESM-C retrieval does **not** beat sequence similarity here.
- Novel families (n=726): **0% exact** for both — the family isn't in the 2024 reference, so
  retrieval cannot assign it (same ceiling as DIAMOND).
- Novelty detection (flag novel_family vs novel_seq): raw top-1 cosine AUROC ~0.52 (~chance);
  best single signal centroid-margin AUROC 0.66. Mean-pooled ESM-C cosine is globally saturated
  (~0.99 both buckets), so absolute similarity is a poor novelty cue.

**Takeaway → Step 7**
- Off-the-shelf ESM-C embeddings are not tuned to separate CAZy families. This is exactly what a
  **CAZy-supervised contrastive head** should fix: pull same-family together, push different-family
  apart, so (a) known-family recall clears DIAMOND and (b) novel families land in low-density regions
  that a margin/energy score can flag. Figure: docs/figures/esmc_retrieval_findings.png.

---

## Day 2 (cont.) — trained heads (Step 7) + dbCAN HMMER baseline featured

**Trained heads on frozen ESM-C** (scripts/train_heads.py): SupCon projection head + softmax
classifier, val classifier acc 0.966 in ~30 s (1 GPU). Known-family (novel-seq) exact subfamily:
contrastive kNN **89.8%** (up from off-the-shelf 85.8%), classifier 89.1% — still below DIAMOND 96.2%.

**dbCAN HMMER featured** (user priority). Canonical baseline = `run_dbcan CAZyme_annotation` →
overview.tsv (dbCAN_hmm / dbCAN_sub / DIAMOND). Exact dbCAN dev cutoffs (from parameter.py):
DIAMOND E-value 1e-102; dbCAN HMM 1e-15 / cov 0.35; dbCAN-sub 1e-15 / cov 0.35. All cutoffs now in
REPRODUCE.md. dbCAN HMMER (current DB) subfamily recall: overall 0.822, novel-seq 0.794,
novel-family 0.975 — the 0.975 is DB-vintage (current DB already has all 20 new-2025 families).

**Two-granularity + truly-novel split** (the honest evaluation):
- 14/20 "novel families" are new *subfamilies* of parents already in 2024; only **6 are truly-novel
  base families** (CBM104, CBM3, CBM8, GT109, GT119, PL29; n=21 eval seqs).
- Parent-family recall on novel subfamilies: DIAMOND 98.9%, contrastive kNN 97.9%, HMMER-2024 94.9%.
- **Truly-novel base families: 0.000 for ALL methods** (DIAMOND, HMMER-2024, off-the-shelf & trained
  ESM-C) — every method makes a confident wrong call, none abstain. This is THE blind spot and the
  motivation for the structure tier (Step 8) + a working novelty detector.

**Supplementary temporal HMMER** (scripts/hmm_baseline_2024.py): 223 base-family HMMs from 2024 only
(mafft --auto + hmmbuild, hmmscan -E 1e-3). Confirms HMMER has the same truly-novel blind spot.

**Corrections**: AUROC now tie-aware (avg-rank) in both retrieval + train scripts; vote-purity AUROC
0.506 (was mis-reported 0.60); head_metrics.json v2 carries corrected values. Fair DIAMOND rerun uses
the identical fungal 2024 reference (96.2% novel-seq, within 0.3% of all-kingdom run).

Figures: docs/figures/benchmark_dbcan_vs_plm.png (featured), trained_heads_comparison.png.
