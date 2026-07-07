# dbCAN4 Advanced-Methods Module — Design & Architecture

**Project:** Bringing protein-language-model and structure-similarity CAZyme annotation into dbCAN.
**Status:** Design document (Step 1 of the one-week prototype plan). Living document — updated as benchmark results arrive.
**Primary host:** `met.unl.edu` (128 CPU, 8× RTX A5500 24 GB, 386 GB RAM).
**Author/maintainer:** dbCAN development team.

---

## 1. Motivation — why sequence similarity is not enough

The current dbCAN stack (dbCAN3 / `run_dbcan` v4–v5) assigns CAZy families by three sequence-similarity engines run in parallel and reconciled into a consensus `overview.tsv`:

1. **HMMER** against the dbCAN HMM database (family-level profile HMMs built from CAZy);
2. **dbCAN_sub** (HMMER against subfamily-level profiles with EC/substrate mapping);
3. **DIAMOND** BLASTp against the CAZyDB sequence database.

All three detect homology in **sequence space**. This is fast, precise, and interpretable, but it has a structural blind spot: CAZymes that are **remote homologs** of known families — sharing fold, catalytic mechanism, and active-site geometry but having drifted below the sequence-identity detection threshold — are missed or called with low confidence. This is not a hypothetical failure mode:

- In the CAZy ontology, enzymes **within a class share a conserved fold, mechanism, and catalytic residues**, even where sequence has diverged substantially — so fold is conserved well past the point where sequence alignment fails.
- CAZyLingua (Thurimella et al., *BMC Bioinformatics* 2025), applying pLM embeddings to a mother/infant metagenomic gene catalog, **identified over 27,000 putative CAZymes missed by other tools**, including horizontally-transferred enzymes, and functionally validated a CE17 enzyme predicted to be overabundant in Crohn's disease.
- Structure-based annotation frameworks in adjacent domains show the same gain: Phold (phage annotation) reports that a ProstT5+Foldseek structural approach **outperforms sequence-based homology in functional-annotation sensitivity** while remaining fast, precisely because structure resolves homology **below the ~25% sequence-identity "twilight zone."**

The field is therefore at a transition from **sequence-similarity** to **representation-similarity** (pLM embeddings) and **structure-similarity** (fold search) methods. CAZyLingua (pLM-only) and DEFT (pLM + structure) are early movers. As dbCAN developers, our goal is not to wrap those tools but to build our **own** advanced-methods module, benchmarked against the dbCAN baseline, and to fold the winning components into a future **dbCAN4**.

### The dbCAN4 recall target, stated precisely

> Recover CAZymes that HMMER/DIAMOND miss (or call with low confidence), assign them the correct CAZy family/subfamily where the family exists, and abstain gracefully (backing off to class level) where the query belongs to a fold/family not represented in the reference — all while not degrading the precision users rely on from dbCAN3.

---

## 2. Survey of advanced methods

Two axes organize the design space: **what representation** carries the homology signal (sequence embedding vs. explicit structure), and **how** a family label is assigned from it (nonparametric retrieval vs. trained head vs. fold search). All entries below are grounded in retrieved literature.

### 2.1 Protein-language-model (pLM) similarity

The premise: a pLM trained by masked-language-modeling over hundreds of millions of sequences maps each protein to a vector where **functional/structural relatedness shows up as embedding proximity**, even when sequences are too diverged for alignment. Assignment schemes on top of embeddings:

| Scheme | How family is assigned | Strengths | Weaknesses |
|---|---|---|---|
| **kNN retrieval** | Nearest reference embeddings vote | Zero training; works with 2–3 refs; interpretable ("nearest = GH5 @ d") | Sensitive to distance threshold; a single noisy neighbor can mislead |
| **Nearest-centroid / prototype** | Distance to per-family mean embedding | Robust to single noisy neighbors; this is CLEAN's mechanism | Assumes families are roughly unimodal in embedding space |
| **Supervised contrastive head** | Fine-tuned projection pulls same-family together, then centroid/kNN in learned space | State of the art for enzyme EC (CLEAN); improves separation of confusable families | Needs labeled training; calibration effort |
| **MLP / linear probe** | Softmax classifier on frozen embeddings | Fast; strong on well-populated families | Weak on long-tail families; fixed label space |
| **Hierarchical (class→family→subfamily)** | Predict coarse level first, back off when fine level uncertain | Matches CAZy ontology; graceful abstention | More moving parts |
| **Retriever↔predictor refinement** | Retriever and trained predictor refine each other (ProtIR) | Combines strengths of both | Heaviest; likely beyond 1-week scope |

**Reference points:**
- **CAZyLingua** — the first tool to use pLMs for CAZyme family/subfamily classification. Architecture: **ProtT5** embeddings → a two-stage classifier (a random-forest CAZyme/non-CAZyme gate, then a feed-forward neural-network multiclass head over the CAZy family/subfamily ontology). Precision and recall comparable to HMM-based methods while outperforming purely sequence-based approaches.
- **CLEAN** (Yu et al., *Science* 2023) — **contrastive learning-enabled enzyme annotation**; assigns EC numbers with "better accuracy, reliability, and sensitivity compared with the state-of-the-art tool BLASTp," and can annotate understudied enzymes, correct mislabels, and flag promiscuous enzymes. Mechanism: supervised contrastive training + distance-to-cluster-centre assignment. This is the model for our contrastive head.
- **EnzymeHunter**, **ProtIR** — hierarchically-aware contrastive and retriever↔predictor refinement variants, cited in the design as the next tier.

**Our pLM choice — ESM-C (ESM Cambrian), not ESM-2.** ESM-C is the newer successor to ESM-2 with better scaling at matched parameter counts, available in 300M / 600M / 6B sizes. On the 24 GB A5500s the **600M** model is the comfortable default for batched embedding; 6B is reachable with care. ESM-2 is retained only as an optional ablation. CAZyLingua used ProtT5 and GraphEC used ProtT5-XL-U50 — using ESM-C is a deliberate, defensible upgrade rather than a re-implementation of either.

### 2.2 Structure similarity

The premise stated in §1: fold is conserved far past sequence. Two ways to bring structure in:

**(a) Explicit fold search — Foldseek.** Foldseek (van Kempen et al., *Nature Biotechnology* 2023) describes tertiary residue interactions as a **3Di structural alphabet**, turning structure search into fast sequence-style search: it is **four to five orders of magnitude faster** than Dali/TM-align/CE while retaining 86–133% of their sensitivities. This is the enabling technology for structure-based CAZyme annotation at scale.

**Reference database — CAZyme3D (our own group's resource).** CAZyme3D (Zheng, Yin et al., bioRxiv 2024, `https://pro.unl.edu/CAZyme3D/`) is a dedicated 3D-structure database for CAZymes: **870,740 AlphaFold-predicted structures** (Whole dataset) and a **188,574-sequence non-redundant ID50 subset** organized by a hierarchical classification that extends the CAZy class/clan/family/subfamily levels with new structure-defined levels (subclasses, structural-cluster groups, SCs). Using CAZyme3D as the Foldseek target database means a query protein's fold is matched against **CAZyme structures specifically**, with family labels attached — exactly what we need.

**(b) Structure-aware embeddings — SaProt.** SaProt (Su et al., *ICLR* 2024) builds a **structure-aware vocabulary**: each residue token is fused with a **Foldseek 3Di structure token**, so the pLM sees explicit local geometry rather than inferring it. Trained on ~40M structures, it surpasses sequence-only baselines across 10 downstream tasks. SaProt gives us a **second, orthogonal structure signal** distinct from Foldseek's geometric search: a structure-aware *embedding* we can run retrieval/classification on exactly as we do for ESM-C.

**The folding dependency, and a way around it.** Both structure routes need a structure (or its 3Di tokens). Options: (i) **ESMFold** on the query (GPU, ~seconds–minutes per protein) → real 3Di; (ii) **ProstT5** predicts 3Di tokens *directly from sequence* — the Bilingual/ProstT5 work reports predicting 3Di for a 1,787-protein proteome in **~44 s on GPU vs. ~48 h for ColabFold structure prediction**, an over-three-orders-of-magnitude speedup, enabling Foldseek-style structure search *without folding*. Phold uses exactly ProstT5+Foldseek in production. **Design consequence:** the structure tier can run in a fast "predicted-3Di" mode (ProstT5) for whole proteomes and a slow "true-fold" mode (ESMFold) for the hardest cases. For the one-week prototype we fold a bounded subset with ESMFold and note ProstT5 as the scale path.

### 2.3 Geometric / active-site models (roadmap tier)

- **GraphEC** (Song et al., *Nature Communications* 2024) — ProtT5-XL-U50 embeddings on **ESMFold-predicted structures**, geometric graph learning to predict EC number **and active sites**; explicitly argues prior EC predictors under-use active-site and structural characteristics.
- **CLEAN-Contact** (*Communications Biology* 2024) — CLEAN's contrastive framework augmented with structural inference for improved functional annotation.

These predict at residue/active-site resolution and are heavier to train; scoped as **Phase-2**, not week-1.

### 2.4 Substrate specificity (roadmap tier)

- **EZSpecificity** (*Nature* 2025, "Enzyme specificity prediction using cross-attention graph neural networks") — a cross-attention, SE(3)-equivariant GNN over enzyme + substrate graphs predicting **which substrates an enzyme acts on**, for the millions of enzymes lacking specificity data.

This answers dbCAN's *next* question after family assignment. Today dbCAN maps family→substrate via `fam-substrate-mapping.tsv` and CGC rules (a lookup). An EZSpecificity-style per-enzyme model is a compelling **Phase-2 substrate module** for dbCAN4, cited here as roadmap, not built in week 1.

---

## 3. Proposed dbCAN4 architecture

A **two-tier** design: a fast sequence-only tier that runs on every protein, and a structure tier reserved for hard/ambiguous cases — with a fusion layer producing one calibrated call.

![dbCAN4 advanced-methods architecture]({{artifact:art_4947927e-ba6d-477a-a081-b183ca5ab94b}})

```
                              query proteins (FASTA)
                                       │
                    ┌──────────────────┴───────────────────┐
                    │                                       │
         ┌──────────▼───────────┐              ┌────────────▼─────────────┐
         │  dbCAN3 BASELINE     │              │  TIER 1 — SEQUENCE pLM   │
         │  (unchanged)         │              │  ESM-C 600M embeddings   │
         │  HMMER · dbCAN_sub · │              │   ├─ kNN retrieval       │
         │  DIAMOND → overview  │              │   ├─ nearest-centroid    │
         └──────────┬───────────┘              │   └─ contrastive head    │
                    │                          │      (+ MLP compare)     │
                    │                          │   hierarchical backoff   │
                    │                          └────────────┬─────────────┘
                    │                                       │
                    │           gate: low-confidence / baseline-missed / disagreement
                    │                                       │
                    │                          ┌────────────▼─────────────┐
                    │                          │  TIER 2 — STRUCTURE      │
                    │                          │  ESMFold (subset) or     │
                    │                          │  ProstT5 3Di (scale)     │
                    │                          │   ├─ Foldseek vs         │
                    │                          │   │   CAZyme3D           │
                    │                          │   └─ SaProt embedding    │
                    │                          │       retrieval          │
                    │                          └────────────┬─────────────┘
                    │                                       │
                    └───────────────┬───────────────────────┘
                                    │
                       ┌────────────▼─────────────┐
                       │  FUSION / CONSENSUS       │
                       │  weighted evidence →      │
                       │  family call + confidence │
                       │  + abstain option         │
                       └────────────┬─────────────┘
                                    │
                       ┌────────────▼─────────────┐
                       │  BioForge DB + web UI     │
                       │  advanced vs baseline,    │
                       │  "missed-CAZyme" flag,    │
                       │  versioned release        │
                       └───────────────────────────┘
```

### 3.1 Tier 1 — sequence pLM (runs on everything)
- **Embed** each protein with ESM-C 600M (mean-pooled per-sequence vector; per-residue retained for future domain-level work).
- **Assign** family by three interchangeable schemes benchmarked head-to-head: kNN, nearest-centroid, and a supervised-contrastive projection (+ an MLP head as a parametric comparison).
- **Hierarchical backoff:** if no family passes threshold, report the CAZy class (GH/GT/PL/CE/AA/CBM) if *that* is confident, else abstain.

### 3.2 Tier 2 — structure (hard cases only)
- **Gate in** proteins that are baseline-missed, low-confidence in Tier 1, or where schemes disagree — keeps the expensive folding bounded.
- **Fold** with ESMFold (subset) or predict 3Di with ProstT5 (scale path).
- **Two orthogonal signals:** Foldseek 3Di search vs. CAZyme3D (fold-similarity, family label from best structural hit) **and** SaProt structure-aware embedding retrieval.

### 3.3 Fusion
- Combine Tier-1 and Tier-2 evidence into one call with a transparent weighted score and a confidence value; abstain below threshold. Design detail in Step 7; the principle is that **agreement across orthogonal signals (embedding + fold) is the strongest evidence for a genuine remote-homolog CAZyme**, which is precisely the DEFT premise implemented as our own rule.

---

## 4. Evaluation design

**Primary: temporal holdout (2024 → 2025).** Train/reference on the **2024** CAZy release; test on **2025**. Diff the releases to isolate 2025-only entries and split them into:
- **novel-sequence** — sequences newly added to families that already existed in 2024 → *can we assign the right family to a protein the model never saw?*
- **novel-family** — families/subfamilies appearing only in 2025 → the label is out-of-training, so the correct behavior is **class-level backoff or abstention**; this tests calibration, not raw accuracy.

This mirrors real prospective deployment far better than a random split: the 2025 additions are genuinely new content, and recovering them — especially entries the 2024-era HMMER/DIAMOND profiles miss — is direct evidence the advanced methods find *novel* CAZymes.

**Secondary: identity-controlled split** (≤30% identity via MMseqs2/DIAMOND clustering) — the classic remote-homolog stress test, complementary to the temporal split.

**Fungal proteome set** from `/array1/xinpeng/fungi-cazyme-project` provides a realistic application-domain test and the "HMMER/DIAMOND-missed" subset that motivates the whole project.

**Metrics:** family-level precision / recall / F1 (macro and weighted); recovery rate of baseline-missed and 2025-novel CAZymes; long-tail-family performance; calibration (does confidence track correctness?); abstention correctness on novel-family.

---

## 5. One-week execution map

| Day | Steps | Deliverable |
|---|---|---|
| 1 | Survey + this doc; GitHub repo; port envs to met | design doc, repo, `env_setup.md` |
| 2 | Reference/eval splits (2024/2025 diff); dbCAN baseline | split manifests, `baseline_overview.tsv`, missed-set |
| 3 | ESM-C embeddings; kNN + centroid retrieval | embeddings checkpoint, `retrieval_predictions.tsv` |
| 4 | Contrastive + MLP heads | head weights, `clf_predictions.tsv`, curves |
| 5 | ESMFold subset; Foldseek vs CAZyme3D; SaProt ablation | structures, `foldseek_predictions.tsv`, `saprot_predictions.tsv` |
| 6 | Fusion; BioForge integration | `fusion_predictions.tsv`, adapter, web demo |
| 7 | Benchmark, figures, final report | `benchmark_report.md`, figures, roadmap |

Daily: commit + push, dated `daily-log.md` entry.

---

## 6. Integration with BioForge (biodb)

The existing **BioForge** stack (`github.com/Xinpeng021001/biodb`) already ingests dbCAN `overview.tsv` into a versioned SQL DB with a browsable web UI, genome browser, provenance, and REST API. Its `cazyme_annotations` table keys each call by a **`tool`** field (`HMMER`/`dbCAN_sub`/`DIAMOND`) with an `n_tools_support` count — exactly the seam for advanced predictions. Integration (Step 10): extend the `tool` vocabulary (`ESM-C-kNN`, `ESM-C-centroid`, `ESM-C-contrastive`, `Foldseek-CAZyme3D`, `SaProt`, `fusion`), add a **confidence-score** column, write an ingester adapter that loads our prediction TSVs as a **new versioned release**, and surface a web view that flags **advanced-only ("missed by HMMER/DIAMOND")** CAZymes side by side with the baseline. This delivers the database + visualization output directly, reusing BioForge's versioning and provenance.

---

## 7. Risks & honest limitations

- **One week ⇒ proof-of-concept scale.** Trained heads run on *frozen* embeddings; ESMFold structures and the SaProt ablation cover a *bounded subset*, not all references. Numbers are indicative, not final.
- **Long-tail families** have few references; kNN/centroid help but rare families remain hard, and the temporal novel-family subset is deliberately a hard/abstention case.
- **Structure tier throughput** is folding-bound; ProstT5 is the noted scale path but week-1 uses ESMFold on a subset.
- **Precision guard-rail:** every recall gain is reported alongside its precision cost, since not degrading dbCAN3 precision is a hard requirement.
- **DEFT** could not be pinned to a citable record during the survey; the pLM+structure fusion is therefore described from verified analogs (CLEAN-Contact, Phold, CAZyLingua's own Foldseek validation step) and implemented as our own rule, not as a re-implementation of DEFT.

---

## References (verified against retrieved full text / records)

1. Thurimella K. *et al.* **Protein language models uncover carbohydrate-active enzyme function in metagenomics** (CAZyLingua). *BMC Bioinformatics* (2025). doi:10.1186/s12859-025-06286-y.
2. Yu T. *et al.* **Enzyme function prediction using contrastive learning** (CLEAN). *Science* 379:1358–1363 (2023). doi:10.1126/science.adf2465.
3. van Kempen M. *et al.* **Fast and accurate protein structure search with Foldseek.** *Nature Biotechnology* (2023). doi:10.1038/s41587-023-01773-0.
4. Zheng J., Yin Y. *et al.* **CAZyme3D: a database of 3D structures for carbohydrate-active enzymes.** bioRxiv (2024). doi:10.1101/2024.12.27.630555.
5. Su J., Han C., Zhou Y., Shan J., Zhou X., Yuan F. **SaProt: Protein Language Modeling with Structure-aware Vocabulary.** *ICLR* 2024.
6. Song Y. *et al.* **Accurately predicting enzyme functions through geometric graph learning on ESMFold-predicted structures** (GraphEC). *Nature Communications* (2024). doi:10.1038/s41467-024-52533-w.
7. **CLEAN-Contact: improved enzyme functional annotation prediction using contrastive learning with structural inference.** *Communications Biology* (2024). doi:10.1038/s42003-024-07359-z.
8. Heinzinger M. *et al.* **Bilingual language model for protein sequence and structure** (ProstT5). *NAR Genomics and Bioinformatics* 6(4):lqae150 (2024).
9. **Enzyme specificity prediction using cross-attention graph neural networks** (EZSpecificity). *Nature* (2025). doi:10.1038/s41586-025-09697-2.
10. ESM Cambrian (ESM-C) — EvolutionaryScale, 2024 (300M/600M/6B protein language models).

*Bibliographic details above are taken from retrieved full-text PDFs (refs 1, 3, 4, 6, 8), verified abstracts/records (refs 2, 9), and OpenReview/DBLP records (ref 5). Volume/page numbers not confirmed from a retrieved record are omitted rather than guessed.*
