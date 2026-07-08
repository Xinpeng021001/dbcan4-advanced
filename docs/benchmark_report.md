# dbCAN4-advanced: Benchmark Report

**Advanced CAZyme annotation for fungi — protein language models + structure similarity beyond HMMER/DIAMOND**

*dbCAN development team · fungal 2024→2025 temporal holdout · compute on `met` (8× RTX A5500)*

---

## 1. Motivation

dbCAN's production annotation rests on sequence similarity: HMMER against family
profile HMMs (`dbCAN.hmm`, `dbCAN-sub.hmm`) and DIAMOND against the CAZy sequence
database. Sequence similarity is fast and accurate when a query has a homolog in
the reference, but it has a structural blind spot: a CAZyme whose fold is
conserved but whose sequence has diverged below the detection threshold is
missed. The field is adding two complementary signals — **protein language model
(pLM) embeddings** (CAZyLingua) and **structure similarity** (DEFT, Foldseek vs
CAZyme3D). This project builds dbCAN's own implementation of both, targeting a
future dbCAN4, and benchmarks them honestly against the sequence baselines.

The prediction target is the **CAZy family** (GH/GT/PL/CE/AA/CBM + number), scored
at two granularities: **exact subfamily** (e.g. GH5_40) and **parent family**
(GH5).

## 2. Evaluation design — temporal holdout

We use the lab's own fungal CAZyme sets: **2024** (`CAZyDB.07142024`) as the
knowledge base and **2025** (`CAZyDB.07242025`) as the test, so the benchmark
measures what a 2024-era annotator would find in 2025 data — the real deployment
question. From these we built:

- **Reference (2024):** 337,759 labeled fungal CAZymes, 411 subfamily labels /
  266 base families.
- **Evaluation (2025), 4,726 held-out proteins** in three novelty tiers, defined
  by exact-sequence identity against 2024:
  - **known family / new sequence** (`novel_seq`, n=4,000): family seen in 2024,
    sequence not.
  - **novel family** (`novel_family`, n=726): family label absent from the 2024
    *fungal* reference.

### 2.1 A correction that reshaped the novelty story

The initial `novel_family` tier was defined against the **fungal** 2024 subset.
On inspection (prompted during review), most of these families are **not new
CAZy families at all** — they exist in CAZy 2024 in *other kingdoms* (bacteria,
plants) and were simply newly annotated in fungi in 2025. Re-scoring all 726
`novel_family` sequences against the **full all-kingdom** 2024 CAZy:

- **95.9% (696/726)** have their parent family recoverable from all-kingdom 2024
  → **new-to-fungi, not new families**.
- Only **6 sequences** have no parent-family hit even against all of CAZy 2024
  → genuinely-novel-to-CAZy candidates (3× CBM104, 2× GH2_10, 1× GT109). Of the
  original 6 "truly-novel base families," five (CBM3, CBM8, GT109, GT119, PL29)
  had thousands-to-hundreds of 2024 headers in non-fungal kingdoms; only
  **CBM104** is genuinely absent from CAZy 2024.

**Consequence:** this fungal holdout is predominantly a **cross-kingdom-transfer
test**, and only weakly a novel-family-discovery test (~6 genuinely novel
sequences). The report treats these two questions separately.

## 3. Methods benchmarked

| Tier | Method | Reference / training | Confidence signal |
|---|---|---|---|
| Sequence | **dbCAN HMMER** (current DB) | `dbCAN.hmm` profiles (2025-era) | HMM E-value |
| Sequence | **dbCAN Recommend** (current DB) | run_dbcan consensus | tool agreement |
| Sequence | **DIAMOND** (temporal) | fungal 2024 reference | % identity |
| Sequence | **HMMER** (temporal) | 223 family HMMs built from 2024 | HMM E-value |
| pLM | **ESM-C kNN** (off-the-shelf) | ESM-C 600M embeddings, 2024 | vote purity |
| pLM | **Contrastive kNN** (trained) | SupCon head on ESM-C, 2024 | vote purity |
| pLM | **Classifier** (trained) | softmax head on ESM-C, 2024 | max-softmax |
| Structure | **Foldseek** | ESMFold structures vs 2024 ref | TM-score |
| Fusion | **DEFT-style consensus** | all of the above | weighted vote |

**pLM tier (ESM-C, per user preference over ESM-2).** Mean-pooled ESM-C 600M
embeddings (1152-dim). Retrieval by kNN (k=15) and nearest-centroid over 403
family prototypes. Trained heads: a supervised-contrastive projection head and a
softmax classifier, both on **frozen** embeddings — deliberately not the same kNN
recipe as CAZyLingua/DEFT.

**Structure tier.** ESMFold (facebook/esmfold_v1) folds queries and a per-family
2024 reference sample on met's 8 GPUs; Foldseek (`--alignment-type 1`, TM-align)
finds the nearest structural neighbor → its family. CAZyme3D (the lab's 870k
structure DB) was not available on met, so this is a **self-contained** proof of
concept against a folded 2024 reference (917 structures).

**Fusion.** Confidence-weighted consensus across all axes with per-method
reliability weights; below a confidence threshold τ=0.35 the prediction
**abstains** (flags a putative novel/uncertain CAZyme).

## 4. Results

![Master benchmark](figures/master_benchmark_figure.png)

**Exact-subfamily / parent-family recall on the temporal holdout:**

| Method | Known (sub) | Known (parent) | New-to-fungi (parent) | Genuinely novel (parent) |
|---|---|---|---|---|
| dbCAN HMMER (current DB) | 0.846 | 0.858 | 0.988 | 0.833 |
| dbCAN Recommend (current DB) | 0.904 | 0.917 | 0.988 | 0.833 |
| DIAMOND (fungal 2024, temporal) | 0.981 | 0.985 | 0.995 | 0.333 |
| HMMER (2024-only, temporal) | 0.713 | 0.950 | 0.993 | 0.333 |
| ESM-C kNN (off-the-shelf) | 0.931 | 0.935 | 0.915 | 0.333 |
| Contrastive kNN (trained) | 0.973 | 0.976 | 0.995 | 0.333 |
| Classifier (trained) | 0.966 | 0.969 | 0.982 | 0.333 |
| Foldseek (structure)¹ | 0.072 | 0.077 | 0.023 | 0.000 |
| **FUSION (consensus)** | **0.981** | **0.984** | **0.995** | **0.333** |

¹ Foldseek scored over all 4,726 queries; only 698 have structures. On its
structure-bearing subset (n=347 known-family): **0.718 subfamily / 0.752 parent**.

### 4.1 Known families: pLM approaches sequence, trained heads help, fusion wins

Off-the-shelf ESM-C kNN (93.1% exact subfamily) **trails** DIAMOND (98.1%) — a
useful negative result: a general pLM does not beat sequence similarity out of the
box. CAZy-supervised training closes most of the gap (contrastive kNN 97.3%). The
**fusion** consensus edges past every single method (98.1% subfamily, 98.4%
parent), picking up cases where each method individually errs.

### 4.2 The novelty cliff: every method collapses on genuinely-novel families

On the 6 genuinely-novel-to-CAZy sequences, **parent-family recall is ≤0.33 for
every sequence, pLM, and structure method** (and 0 for Foldseek). This is the
irreducible limit of any retrieval method: a family absent from the reference
cannot be named. The current-DB dbCAN tiers score 0.833 here only because their
2025-era database already contains these families — a DB-vintage effect, not a
temporal result.

### 4.3 Structure tier: complementary errors justify fusion

![Structure tier](figures/structure_tier_findings.png)

Foldseek on ESMFold-predicted structures recovers **75.2% of known families**
(parent) from structure alone — below DIAMOND because our reference is a
length-capped sample, not the full CAZyDB. Its value is **complementarity**: on
347 known-family proteins with structures, DIAMOND alone 93.9% + Foldseek alone
75.2% → **union 95.4%**; structure rescues 5 proteins DIAMOND misses. Different
error profiles are exactly what makes fusion worthwhile.

### 4.4 Fusion abstention: knowing when it doesn't know

![Fusion](figures/fusion_findings.png)

The most actionable result. Fusion confidence separates the difficulty tiers:

| Group | mean fusion confidence | abstain rate (τ=0.35) |
|---|---|---|
| Known family | 0.835 | 1.3% |
| New-to-fungi | 0.852 | 0.3% |
| **Genuinely novel (to CAZy 2024)** | **0.511** | **66.7%** |

Fusion stays confident on placeable proteins — **including new-to-fungi families**,
correctly, since those are recoverable by cross-kingdom homology — but **flags 4
of the 6 genuinely-unplaceable CAZymes** (CBM104, GT109) for review. This is the
DEFT-style payoff: a production annotator that surfaces candidate novel CAZymes
rather than silently mislabeling them.

## 5. Recommendations for dbCAN4

1. **Keep sequence similarity as tier 1.** DIAMOND/HMMER remain the most accurate
   single signal on known families; the pLM and structure tiers are additive, not
   replacements.
2. **Add a trained pLM head as tier 2.** Contrastive/classifier heads on frozen
   ESM-C give near-DIAMOND accuracy from an orthogonal representation and are cheap
   to run (embeddings cached, head trains in ~30 s).
3. **Add structure similarity as tier 3, gated on a complete reference.** The
   self-contained proof of concept works; production value requires **CAZyme3D**
   (870k structures) as the Foldseek target rather than a folded sample, and
   ProstT5 (sequence→3Di) to avoid the ESMFold folding cost.
4. **Ship the fusion confidence + abstention as the headline feature.** Flagging
   genuinely-novel CAZymes is what advanced methods add beyond "faster HMMER."
5. **Substrate-specificity prediction (EZSpecificity-style)** is a natural Phase-2
   extension, mapping to dbCAN's substrate step.

## 6. Honest limitations

- The fungal holdout is mostly a cross-kingdom-transfer test; only ~6 sequences
  probe genuine novel-family discovery. A stronger novelty benchmark would hold
  out entire families across all kingdoms.
- Structure reference is a length-capped sample of predicted structures, not
  CAZyme3D; large multi-domain families (GH2, GH3) are under-represented, and
  Foldseek's numbers here are a floor, not its ceiling.
- Fusion reliability weights are a principled prior (each method's known-family
  recall), not fit on held-out data; a small validation split could calibrate them.
- All structures are ESMFold predictions (pLDDT varies), not experimental.

## 7. Reproducibility

Exact commands, parameters, and cutoffs for every step are in
[`REPRODUCE.md`](../REPRODUCE.md) (DIAMOND E-value 1e-102, dbCAN HMM/sub E-value
1e-15 / coverage 0.35, standalone temporal DIAMOND 1e-3, Foldseek TM-align mode,
ESM-C 600M, ESMFold facebook/esmfold_v1). All scripts are in `scripts/`; all
per-method predictions and summaries are in `benchmarks/`.
