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

**(b) Structure-aware embeddings — SaProt.** SaProt (SaProt: Protein Language Modeling with Structure-aware Vocabulary, *ICLR* 2024; Westlake University) builds a **structure-aware vocabulary**: each residue token is fused with a **Foldseek 3Di structure token**, so the pLM sees explicit local geometry rather than inferring it, and it reports gains over sequence-only baselines on downstream function tasks. SaProt gives us a **second, orthogonal structure signal** distinct from Foldseek's geometric search: a structure-aware *embedding* we can run retrieval/classification on exactly as we do for ESM-C. *(Training-corpus size and per-task benchmark numbers to be confirmed against the paper before they appear in any writeup.)*

**The folding dependency, and a way around it.** Both structure routes need a structure (or its 3Di tokens). Options: (i) **ESMFold** on the query (GPU, ~seconds–minutes per protein) → real 3Di; (ii) **ProstT5** predicts 3Di tokens *directly from sequence* — the Bilingual/ProstT5 work reports a large speedup over full structure prediction for generating Foldseek 3Di tokens, enabling Foldseek-style structure search *without folding* (the specific benchmark figures are from a search snippet and should be confirmed against the paper before quoting). Phold uses exactly ProstT5+Foldseek in production. **Design consequence:** the structure tier can run in a fast "predicted-3Di" mode (ProstT5) for whole proteomes and a slow "true-fold" mode (ESMFold) for the hardest cases. For the one-week prototype we fold a bounded subset with ESMFold and note ProstT5 as the scale path.

### 2.3 Geometric / active-site models (roadmap tier)

- **GraphEC** (Song et al., *Nature Communications* 2024) — ProtT5-XL-U50 embeddings on **ESMFold-predicted structures**, geometric graph learning to predict EC number **and active sites**; explicitly argues prior EC predictors under-use active-site and structural characteristics.
- **CLEAN-Contact** (*Communications Biology* 2024; details to be verified against the primary source) — reported as CLEAN's contrastive framework augmented with structural inference for improved functional annotation.

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

## 8. Gray-zone adjudication and the structure tier — resolved (2026-07-09)

Two questions raised earlier in the project's benchmark review are resolved here: (a) how to
build a defensible non-CAZyme negative set given that dbCAN.hmm/DIAMOND/dbCAN_sub each find
more candidate CAZymes than CAZy itself has labeled, and (b) how to get the structure tier
(§2.2, §3.2) actually running, since CAZyme3D was previously "not available on met" and the
structure tier was a self-contained POC.

### 8.1 The gray zone is a genuine mixture, not a labeling artifact

Streaming all 2,226 Mycocosm genomes in `/array1/xinpeng/all_genome` (dbCAN-annotated:
`overview.tsv`, `diamond.out`, `dbCAN_hmm_results.tsv`, `dbCANsub_hmm_results.tsv`,
`non_CAZyme.faa`) and recomputing the tool-hit count directly from the raw per-tool tables
(each already filtered to that tool's own significance cutoff) rather than trusting the
pre-filtered `overview.tsv`/`non_CAZyme.faa` split yields three tiers over 28,192,456 total
protein rows:

| Tier | Rule | Count | % |
|---|---|---:|---:|
| `high_confidence_cazyme` | ≥2 tools agree | 946,270 | 3.36% |
| `gray_zone` | exactly 1 tool hits a CAZy family | 2,844,297 | 10.09% |
| `high_confidence_non_cazyme` | 0 tool hits | 24,401,889 | 86.55% |

A single-tool hit is real homology evidence that simply misses dbCAN's own ≥2-tool consensus
bar — treating it as automatically non-CAZyme (as the current `non_CAZyme.faa` output does) or
automatically CAZyme would both be wrong. Structure-similarity scoring against CAZyme3D_id50
(§8.2) on a 2,483-protein stratified sample shows `structure_evidence_score` cleanly separates
known CAZymes (mean 0.653) from known non-CAZymes (mean 0.463), with the gray zone falling
exactly in between (mean 0.566) — confirming it is a real mixture population. Adjudicating the
sampled gray zone with this signal: **33.4%** structure-supports-CAZyme (candidate recall gain
— real CAZymes dbCAN's consensus rule currently misses), **38.3%** structure-supports-non-CAZyme
(likely sequence-level false positives), **28.3%** remain ambiguous even with an orthogonal
signal (a legitimate abstention population, not a modeling failure). Full detail, caveats, and
the adjudicated TSV are in `synthesis_report.md` / `gray_zone_adjudicated_structure_validated.tsv`.
**This is the negative/decoy set called for in the benchmark review's gap 1** (no
precision/negative set) and **gap 2** (no CAZyme/non-CAZyme gate) — the
`high_confidence_non_cazyme` tier (24.4M proteins) is the negative population for training a
gate; the gray zone is the population to hold out for calibrating abstention (§5.6/gap 1).

### 8.2 Structure tier: unblocked

CAZyme3D_id50 (178,356 AlphaFold structures, the group's own resource, §2.2) is now
downloaded and extracted on met at `/array1/xinpeng/cazyme3d/extracted/cazyme_id50/`. It maps
to Mycocosm proteins by **homology, not accession** (accessions are disjoint UniProt/RefSeq
IDs vs. Mycocosm's JGI gene IDs) — an mmseqs2 search at ≥30% identity/50% coverage finds a
CAZyme3D_id50 homolog for 93.3% of a Mycocosm CAZyme sample, confirming it is usable as the
Foldseek/embedding reference database as designed in §2.2/§3.2.

**ProstT5** (Rostlab/ProstT5) is installed and validated on met: sequence→3Di prediction with
no folding step, run across all 8 GPUs on the full 2,483-protein validation sample
(~15–35s/protein). This is the "scale path" §2.2 anticipated in place of folding every protein.

**SaProt** (westlake-repl/SaProt_650M_AF2) is installed and validated: AA+3Di structure-aware
embeddings via the vendored `foldseek_util.get_struc_seq` helper, producing the orthogonal
structural signal alongside Foldseek/mmseqs2 3Di-string search.

**Bulk predicted-structure sources evaluated and deprioritized:** ESM Metagenomic Atlas 2
(85.9% sequence-hit coverage on a CAZyme sample but median only 51.2% identity — mostly
distant homology, consistent with its bacterial/archaeal/environmental MGnify origin) and
AF3db (only reachable via UniProt accessions Mycocosm proteins don't have; 75.9% coverage of
the indirectly-queryable subset). Neither is worth further integration investment over the
CAZyme3D_id50 + ProstT5 + SaProt + local-ESMFold-for-gaps stack, which is now the recommended
and running structure-tier architecture for §3.2.

### 8.3 Updated architecture note

§3's two-tier diagram is unchanged in shape; the structure tier (Tier 2) is no longer a
"self-contained POC" (per the benchmark review's gap 6) — CAZyme3D_id50, ProstT5, and SaProt
are installed, validated, and running on met, gated on the gray-zone population identified in
§8.1 as the natural "hard cases" input to Tier 2.

---

## References

**Group A — retrieved and read in this session (full text or abstract confirmed):**

1. Thurimella K. *et al.* **Protein language models uncover carbohydrate-active enzyme function in metagenomics** (CAZyLingua). *BMC Bioinformatics* (2025). doi:10.1186/s12859-025-06286-y. *(full-text PDF read)*
2. Yu T. *et al.* **Enzyme function prediction using contrastive learning** (CLEAN). *Science* (2023). doi:10.1126/science.adf2465. *(abstract confirmed; full text paywalled)*
3. van Kempen M. *et al.* **Fast and accurate protein structure search with Foldseek.** *Nature Biotechnology* (2023). doi:10.1038/s41587-023-01773-0. *(full-text PDF read)*
4. Zheng J., Yin Y. *et al.* **CAZyme3D: a database of 3D structures for carbohydrate-active enzymes.** bioRxiv (2024). doi:10.1101/2024.12.27.630555. *(full text read; the 870,740 / 188,574 figures are quoted from it)*
5. Song Y. *et al.* **Accurately predicting enzyme functions through geometric graph learning on ESMFold-predicted structures** (GraphEC). *Nature Communications* (2024). doi:10.1038/s41467-024-52533-w. *(full-text PDF read)*
6. **Improved enzyme functional annotation prediction using contrastive learning with structural inference** (CLEAN-Contact). *Communications Biology* (2024). doi:10.1038/s42003-024-07359-z. *(full-text PDF retrieved; title/venue/DOI confirmed)*
7. **Enzyme specificity prediction using cross-attention graph neural networks** (EZSpecificity). *Nature* (2025). doi:10.1038/s41586-025-09697-2. *(DOI + title confirmed via resolver; full text paywalled, author list/pages not verified here)*

**Group B — identified by title/venue from search results but NOT independently retrieved this session; details to be verified before any publication or external writeup:**

8. **SaProt: Protein Language Modeling with Structure-aware Vocabulary.** *ICLR* 2024 (Westlake University; `github.com/westlake-repl/SaProt`, OpenReview `6MRm3G4NiU`). *(venue/repo from search hits; author list and benchmark specifics not confirmed here)*
9. Heinzinger M. *et al.* **Bilingual language model for protein sequence and structure** (ProstT5). *NAR Genomics and Bioinformatics* (2024), article lqae150. *(from search-result title/record; full text not fetched this session)*
10. ESM Cambrian (ESM-C) — EvolutionaryScale (300M/600M/6B protein language models). *(model family known from the loaded `esmfold2` skill; cite the official release before publication)*

> **Provenance note.** Group A (refs 1–7) references were retrieved this session — full text read for refs 1, 3, 4, 5, 6; abstract confirmed for ref 2; DOI+title confirmed via resolver for ref 7 (full text paywalled). Quoted figures come directly from the retrieved text. Group B (refs 8–10) were located by title/venue in web-search results (whose body text was not preserved in full) or known from the loaded skill; their author lists, volume/page, and benchmark specifics have **not** been independently verified here and must be checked against the primary source before this document is used externally. This note replaces an earlier, overstated "all verified" claim.
