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
