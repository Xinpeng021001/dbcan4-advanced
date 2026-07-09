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
- Parent-family recall on the 705 novel-subfamily-of-known-parent seqs (identical denominator):
  DIAMOND 98.9%, contrastive kNN 97.9%, HMMER-2024 97.7%, classifier 96.9%, off-the-shelf ESM-C kNN 89.9%.
  (For the FULL 726-seq novel_family bucket incl. the 21 truly-novel, these become 0.960 / 0.950 /
  0.949 in unified_scores.json — diluted by the 21 impossible cases; see novel_family_decomposition.json.)
- **Truly-novel base families: 0.000 for ALL methods** (DIAMOND, HMMER-2024, off-the-shelf & trained
  ESM-C) — every method makes a confident wrong call, none abstain. This is THE blind spot and the
  motivation for the structure tier (Step 8) + a working novelty detector.

**Supplementary temporal HMMER** (scripts/hmm_baseline_2024.py): 223 base-family HMMs from 2024 only
(mafft --auto + hmmbuild, hmmscan -E 1e-3). Confirms HMMER has the same truly-novel blind spot.

**Corrections**: AUROC now tie-aware (avg-rank) in both retrieval + train scripts; vote-purity AUROC
0.506 (was mis-reported 0.60); head_metrics.json v2 carries corrected values. Fair DIAMOND rerun uses
the identical fungal 2024 reference (96.2% novel-seq, within 0.3% of all-kingdom run).

Figures: docs/figures/benchmark_dbcan_vs_plm.png (featured), trained_heads_comparison.png.

---

## Day 2 (cont.) — IMPORTANT novelty-definition correction

User flagged that "0.0000 on truly-novel families might be errors — 2025 added new fungi data /
newly-annotated CAZymes that are NOT necessarily novel families, just proteins not annotated before."

**They were right.** My "truly-novel base family" was defined against the FUNGAL 2024 subset only.
Re-checked the 6 supposed truly-novel families against the FULL all-kingdom 2024 CAZy:
- CBM104: 0 headers (genuinely absent from 2024 CAZy)
- CBM3: 2,884 | CBM8: 136 | GT109: 75 | GT119: 65,201 | PL29: 288  -> all EXISTED in 2024 CAZy
  (in bacteria/plants/other kingdoms), just absent or sparse in fungi.

Re-scored all 726 novel_family seqs against all-kingdom 2024 CAZy (DIAMOND):
- **95.9% (696/726)** have parent family recoverable from all-kingdom 2024 -> NEW-TO-FUNGI, not new families
- only **6 seqs** have no parent hit even in all-kingdom 2024 -> genuinely-novel-to-CAZy candidates
  (3x CBM104->AA9, 2x GH2_10->GH114, 1x GT109->GT0; all confident WRONG calls)

**Implication:** this fungal 2024->2025 holdout is mostly a cross-kingdom-transfer test, NOT a
strong novel-family-discovery test (only ~6 genuinely novel seqs). The honest framing: sequence
methods recover new-to-fungi families well via cross-kingdom homology; the structure tier's real
target is the handful of genuinely-novel + the confident-wrong-call cases.
See benchmarks/novelty_stratification_corrected.json.

---

## Day 2 (cont.) — Step 8: structure tier (ESMFold + Foldseek)

Self-contained temporal design (user choice): fold our own 2024 reference + eval set on met's
8 GPUs with ESMFold (facebook/esmfold_v1 via transformers), build a Foldseek DB, TM-align search.

**Pipeline built:**
- scripts/fold_esmfold.py    — shardable ESMFold folder (8-GPU), sanitizes non-standard residues (*->X)
- scripts/select_fold_sets.py — picks eval + per-family reference reps to fold
- scripts/foldseek_search.py  — foldseek easy-search (--alignment-type 1 TM-align), best hit -> family

**Folded:** 699 eval structures + 917 reference structures (817 initial <=600aa + 100 supplementary
for large families GH2/GH3/GH31/GH36 etc. that the 600aa cap had removed). ESMFold is O(L^2);
long fungal CAZymes (>700aa) dominate cost, so we length-capped and sampled.

**Results (Foldseek vs 917-structure 2024 reference):**
- Known families (novel_seq, n=347): **71.8% exact subfamily, 75.2% parent** — structure similarity
  DOES recover CAZy families from *predicted* structures. Below DIAMOND's 96% because our reference
  is a length-capped SAMPLE, not the full CAZyDB (a complete structural ref = CAZyme3D would close this).
- Novel families (n=350): 0% — same ceiling as all methods (family absent from reference).

**Complementarity (the fusion justification), on 347 known-fam eval proteins with structures:**
- DIAMOND alone 93.9% | Foldseek alone 75.2% | **Union (fusion ceiling) 95.4%**
- Structure rescues **5 proteins DIAMOND misses** — different errors -> motivates Step 9 fusion.

**Honest limitations documented:** (1) sampled/length-capped reference underperforms full CAZyDB;
(2) large multi-domain families under-represented; (3) predicted (not experimental) structures.
The value shown is COMPLEMENTARITY, not head-to-head beating DIAMOND on a partial reference.
See benchmarks/foldseek_summary.json, struct_seq_complementarity.json; docs/figures/structure_tier_findings.png.

---

## Day 2 (cont.) — Step 9: DEFT-style fusion + consensus scoring

scripts/fusion_consensus.py: confidence-weighted consensus over all retrieval axes —
sequence (DIAMOND, HMMER-2024), pLM (contrastive kNN + softmax classifier), structure (Foldseek).
Each method emits (family, confidence); votes weighted by per-method reliability prior; winning
family's normalized score is the fusion confidence; below tau=0.35 the call ABSTAINS (flags novel/uncertain).

**Fusion beats every single method:**
- Parent-family (known, n=4000): fusion 98.4% vs ContrastiveKNN 97.6%, DIAMOND 97.3%
- Exact-subfamily (known, n=4000): fusion 98.1% vs ContrastiveKNN 97.3%, DIAMOND 97.0%

**Abstention correctly flags the genuinely-hard cases (the actionable result):**
- known family:            abstain 1.3%
- new-to-fungi:            abstain 0.3%  (correctly confident - recoverable by cross-kingdom homology)
- genuinely-novel-to-CAZy: abstain 66.7% (mean conf 0.51 vs 0.84) - fusion flags CBM104, GT109 for review

This is the DEFT-style payoff: consensus stays confident on placeable proteins (incl. new-to-fungi)
but flags genuinely-unplaceable CAZymes. See benchmarks/fusion_summary.json,
fusion_abstention_analysis.json; docs/figures/fusion_findings.png.

---

## Day 2 (cont.) — Step 11: master benchmark + final report

- scripts assembled a unified master benchmark (results/master_benchmark_v2.{tsv,json}) scoring
  all 9 methods at exact-subfamily + parent-family across 3 novelty tiers.
- docs/benchmark_report.md: full report (motivation, temporal-holdout design, the novelty-definition
  correction, per-tier results, structure complementarity, fusion abstention, dbCAN4 recommendations,
  honest limitations, reproducibility pointer).
- docs/figures/master_benchmark_figure.png: two-panel headline figure (known-family recall + novelty cliff).

**Headline findings:**
- Known families: fusion 98.1% subfamily (beats DIAMOND 98.1% marginally, every single method);
  off-the-shelf ESM-C trails DIAMOND, trained heads close the gap.
- Genuinely-novel families (n=6): 0.33 parent recall for ALL retrieval methods = irreducible cliff.
- Fusion abstention flags 66.7% of genuinely-novel vs 0.3% of new-to-fungi = the actionable signal.

**Recommendations for dbCAN4:** keep sequence tier-1; add trained pLM head tier-2; add structure
tier-3 gated on a complete reference (CAZyme3D + ProstT5); ship fusion confidence + abstention as
the headline feature; substrate-specificity (EZSpecificity-style) as Phase 2.

---

## Day 2 (cont.) — correction: fusion-vs-DIAMOND parent tie + overlap-recall labels

Review caught two presentation errors in the Step-11 deliverables (fixed):
- master_benchmark_v2.tsv columns were named "exact_subfamily" but hold OVERLAP recall
  (prediction shares >=1 family with truth). Renamed to *_overlap and clarified in the report table.
- benchmark_report.md 4.1 claimed fusion "edges past every single method (98.1% sub, 98.4% parent)".
  True at subfamily (98.10 vs DIAMOND 98.05) but NOT parent — DIAMOND 98.45 > fusion 98.38.
  Corrected: methods are effectively TIED on known families; fusion's differentiator is novelty
  abstention (4.4), not a higher known-family number.

---

## Day 4 — 2026-07-08: fair-database rerun (fixing the DIAMOND vintage confound)

**Motivation (from user).** `run_dbcan download` pulls the *current* database,
which includes HMMs and CAZy sequences added in 2025 — the same year as the eval
set. For a 2024→2025 temporal holdout that leaks the answer into the annotator's
reference, "especially DIAMOND." User built a temporally-clean 2024 DB at
`/array1/xinpeng/dbcan_db_2024` (old dbCAN HMM, old DIAMOND, old dbCAN-sub HMM)
and asked to rerun run_dbcan against it.

**Done**
- Reran `run_dbcan CAZyme_annotation --methods diamond,hmm,dbCANsub` on all 4,726
  eval_2025 proteins vs the 2024 DB (job on met, exit 0, 4,575 annotated).
- Scored the 2024-DB run and re-scored the 2025-DB run with **one identical
  parser** (norm_fams: strip coords + subfamily `_eNNN`, split multi-domain on
  `+`/`|`, family regex). Extended both to parent-level metrics.
- Re-scored **every** method (ESM-C kNN, contrastive kNN, classifier, Foldseek,
  fusion, custom DIAMOND) with the same scorer so the master table is internally
  consistent. Prior numbers all reproduced.
- New artifacts: benchmarks/master_benchmark_v3.{tsv,json},
  dbcan_db_2024_vs_2025_comparison.tsv, dbcan2024_eval2025_scored.json,
  dbcan2025_eval2025_rescored.json, allmethods_rescored.json;
  docs/figures/fair_benchmark_2024db.png; benchmark_report.md §4.5.

**Headline findings (subfamily-exact recall, 2024 → 2025 DB):**
- **HMMER known-family byte-identical: 0.794 → 0.794.** Old family HMMs are not
  rebuilt between releases, so a query matching a pre-2024 family gets the same
  answer either way. Confirms the user's intuition — only novel families move.
- **DIAMOND contaminated in BOTH tiers**: known/new-seq 0.896 → 0.980 *and*
  novel-family 0.001 → 0.992. The eval sequences are themselves in the 2025 CAZy
  release the current .dmnd is built from → near-self hits. User's concern
  confirmed; DIAMOND is the most affected method.
- **Novel-family is the smoking gun**: ≈0.00 (fair) → 0.97–0.99 (current) for all
  tiers at subfamily level. But parent-level even the fair DB recovers ~95% →
  cross-kingdom transfer, not genuine novelty (restates the §2.1 correction).
- pLM + structure + fusion tiers were already fair (train/retrieve on 2024 only);
  no rerun needed. Foldseek likewise clean.

**Verification / provenance (auditor-prompted).**
- 2024-DB vintage verified by HMMER `hmmbuild` DATE stamps inside the .hmm files
  (dbCAN.hmm newest Aug 2024; dbCAN-sub.hmm all May 2022) + absence of the 12
  genuinely-new-in-2025 family profiles — NOT by download-path date. The
  dbCAN-sub.hmm was served from a 2025-dated URL (`db_v5-1_3-11-2025`) but its
  content is a 2022 build; profile count alone is not proof of vintage.

**Blocker hit + fixed.**
- dbCAN-sub silently produced empty output on the 2024 DB. Root cause: dev
  run_dbcan writes the subfamily DB as `dbCAN_sub.hmm` (underscore) but its
  annotator config reads `dbCAN-sub.hmm` (hyphen); the wrapper masked the real
  `FileNotFoundError` as "empty output generated for downstream compatibility."
  Fix: `ln -s dbCAN_sub.hmm dbCAN-sub.hmm` in the DB dir. (Worth reporting
  upstream as a naming-consistency bug in the dev build.)

**Next**
- Optional: hold out entire families across all kingdoms for a stronger
  novel-family-discovery benchmark (current holdout is mostly cross-kingdom).
- Optional: wire the 2024-DB fair numbers into the Nextflow benchmark contract.

---

## Day 4 (cont.) — three-level scoring + fair pLM-vs-run_dbcan comparison

User asked to (1) verify ESM / other methods compare *fairly* with run_dbcan, and
(2) score at three CAZy levels — class (AA/GT/…), family (GH2), subfamily (GH13_1) —
plus multidomain.

**Fairness verification.** Re-scored all 14 methods through one scorer on the
identical truth. Confirmed ID alignment: all prediction files map cleanly onto the
4,726 truth keys (142 GenBank + 4,584 JGI), 0 not-in-truth. No method is scored on a
different protein set → the comparison is fair in denominator and temporal setting.
The pLM heads were already temporally clean (train/retrieve on 2024 only).

**Three-level results (overlap recall, overall, temporally clean):**
- Class level is near-saturated (0.90–1.00) for all sequence + pLM methods.
- Trained Contrastive-kNN (0.978 family / 0.830 subfamily) beats every run_dbcan tier
  (best 0.920 / 0.773) and trails only custom DIAMOND — the fair win we needed to
  justify the pLM tier for dbCAN4.
- Monotonicity class ≥ family ≥ subfamily holds with 0 violations across 15 methods.

**dbCAN-sub subfamily fix.** Its native output is ECAMI cluster IDs (GH13_e122), a
different namespace than CAZy subfamilies (GH13_1); the raw parser scored it a
spurious 0.000 at subfamily. Mapped each cluster → dominant CAZy subfamily via the
composition column: subfamily-exact 0.000 → 0.427 (on n=1,741 has-subfamily truth).
Still trails others because not every ECAMI cluster = one CAZy subfamily.

**Multidomain (n=332), family exact set-match.** DIAMOND/HMMER recover all domains
(custom DIAMOND 0.885, run_dbcan DIAMOND 0.876). Single-label pLM/fusion score ≈0.006
EXACT (one call per protein by construction) but overlap ~0.98–1.00, Jaccard ~0.49 —
architectural, not accuracy: pLM tier needs per-domain segmentation for dbCAN4.

**Fixed** a 0.982→0.980 rounding error in the §4.5 DIAMOND known/new-seq cell
(scored value is 0.9802) in both report and this log.

**Artifacts:** benchmarks/master_benchmark_v4_threelevel.tsv (270 rows: 15 methods ×
6 buckets × 3 levels), threelevel_all_methods.json, foldseek_structure_subset_scored.json;
docs/figures/threelevel_benchmark.png; benchmark_report.md §4.6.
