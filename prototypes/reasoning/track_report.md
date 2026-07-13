# dbCAN4-advanced Reasoning-Assistance Track — Prototype Report

This is a standalone reasoning-assistance prototype built on top of the dbCAN4-advanced fusion
pipeline (sequence-similarity baselines HMMER/DIAMOND + an ESM-C protein-language-model tier of 3
retrieval heads — kNN, centroid, contrastive/classifier — combined by a confidence-weighted fusion
layer). It adds four capabilities the base pipeline does not have on its own: (1) grounded
plain-language reasoning reports per protein, (2) an explicit disagreement/review-triage rule,
(3) conformal prediction sets with a coverage guarantee, and (4) an evaluated OOD/novelty score. All
four are demonstrated on 3 real held-out 2025 fungal CAZymes and evaluated where applicable against
the full n=4,726-protein eval slice.

## 1. Grounded per-protein reasoning

For each of the 3 example proteins (267317, 602276, 169208), a plain-language rationale was drafted
via `host.llm` under a hard grounding constraint: every quantitative claim in the narrative must be
traceable to a value in the staged evidence bundle (per-head predictions/confidences, Pfam domains,
DeepTMHMM topology, CLEAN EC confidence, physicochemistry, ESMFold pLDDT), and no family/EC
assignment may be invented or overridden. The three worked cases:

- **602276 (GH11):** all 4 signals (kNN, centroid, contrastive-classifier, contrastive-kNN) and the
  Pfam domain (PF00457 Glyco_hydro_11, full coverage) agree unanimously; EC (3.2.1.8, conf 0.99),
  secretion signal, and structure (pLDDT 82.3) all corroborate. Clean consensus, no review needed.
- **267317 (GH78, true multidomain GH28+GH78):** the ESM-C kNN and contrastive heads call GH78
  (conf 0.995/0.950), the centroid head dissents to GH92 (conf 0.977). Pfam finds *two* domains —
  PF00295 (GH28-like, N-terminal, high coverage) and PF17389 (GH78-diagnostic, C-terminal) —
  confirming the protein genuinely has both domain signatures, which the fusion call (GH78, but at a
  reduced confidence of 0.667) only partially resolves. EC confidence is very low (0.105).
- **169208 (GH183, the designed failure/disagreement case):** kNN confidently (0.986) calls the
  WRONG family (GH43_6, purity only 0.267); only the centroid head recovers the true GH183
  (conf 0.984); the contrastive classifier calls a third, low-confidence family (PL42, 0.297). Fusion
  followed the confident-but-wrong kNN vote (GH43_6, final confidence only 0.308). The only Pfam
  domain (DUF4185) is uninformative, EC confidence is near-zero (0.0014), and localization carries an
  internal inconsistency (called Extracellular despite `sp_prediction=NO_SP`) that the reasoning
  narrative correctly flagged as a data-quality issue rather than silently resolving.

Full reports: `report_267317.md`, `report_602276.md`, `report_169208.md`.

## 2. Disagreement detection + review triage

A fixed, interpretable score — `0.40·(1−agreement_frac) + 0.35·(1−fusion_conf) + 0.25·min(conf_collapse,1)`
— composes signals the fusion layer already emits (vote agreement, fusion confidence, and the gap
between the best single-head confidence and the fused confidence) into three tiers: ACCEPT / WATCH /
FLAG. On the 3 example proteins, the rule reproduces exactly the triage a human reviewer would want:

| Protein | Score | Tier |
|---|---|---|
| 602276 | 0.010 | **ACCEPT** — 4/4 agreement, high confidence |
| 267317 | 0.299 | **WATCH** — 3/4 agreement, real multidomain complexity, correct call |
| 169208 | 0.612 | **FLAG** — 2/4 agreement, fusion followed a confident-but-wrong head |

Applied to the full eval slice (using the classifier's own confidence as the closest available
fusion-like proxy, since a genuine fusion-confidence field only exists for the 3 demo proteins),
FLAG-tier accuracy (33.9%) is far below
ACCEPT-tier accuracy (83.8%), and the
score achieves AUROC 0.611 for ranking wrong vs. correct calls — a real, if modest, discriminative
signal, honestly reported as a ranking aid for routing the worst ~3% of calls to a human, not a strong
standalone classifier. Full detail: `triage_report.md`, `triage_rule.py`, `triage_eval_slice.csv`.

## 3. Conformal prediction sets

A split-conformal procedure builds a per-protein candidate menu (union of predictions across
classifier / ESM-C-kNN / ESM-C-centroid / ESM-C-contrastive-kNN / DIAMOND, each with its own
confidence) and calibrates a nonconformity threshold on a held-out half of the eval slice
(n=2,363 calibration / 2,363 test). Below the achievable ceiling, calibration is well-behaved: a 70%
coverage target achieves 79.9%
empirical coverage with a mean set size of 0.95
(95.2% singleton sets); an 80% target
achieves 81.1%. **Above the
ceiling, coverage saturates**: the true family is entirely absent from every head's candidate menu for
15.6% of test proteins (driven almost entirely by `novel_family` cases), so 90%/95%
targets cannot be met by any threshold — reported honestly as a real limitation rather than hidden.
On the 3 demo proteins, the procedure surfaces exactly the same ambiguity the triage rule flagged: at
a 90% (saturated) target, 267317's set widens to {GH78, GH92} and 169208's set widens to
{GH43_6, GH183, PL42} — both correctly containing the true family, but only alongside the competing
wrong call(s), while 602276 stays a clean singleton {GH11} at every target. Detail: `conformal_report.md`,
`conformal.py`, `conformal_calibration_results.csv`, `conformal_demo_predictions.json`.

## 4. OOD / novelty scoring — evaluated against baseline

Reproducing the current baseline exactly (off-the-shelf ESM-C `centroid_margin`,
AUROC **0.6548** vs. the reported 0.6549), two alternative novelty scores were
built and evaluated for detecting `novel_family` vs. `novel_seq` on the same eval slice:

- **Unsupervised energy score** (min-max-combining 5 confidence/margin signals, no label fitting):
  AUROC **0.7386**.
- **Supervised logistic-regression score** (10 signals, evaluated by 5-fold cross-validation to avoid
  overfitting inflation): AUROC **0.7856 ± 0.0249**, vs.
  **0.6540** for the baseline signal on the identical folds.

Both scores exceed the 0.6548 baseline, and the supervised improvement holds
across all 5 CV folds — a genuine, reproducible gain. This is reported as an honest **positive**
result on the specific task evaluated (novel_family vs. novel_seq, within the space of
already-computed per-head scalar signals), while flagging the real scope limits: no true ESM-C
embedding-space Mahalanobis/energy score was computed (raw embeddings are not in this staged bundle,
only derived scalars), and none of this addresses the harder, unaddressed question of detecting
proteins outside CAZy family space entirely. On the 3 demo proteins, the energy score ranks 169208
(FLAG case) at the 100th percentile of the eval-slice novelty distribution and 267317 (WATCH case) at
the 97.8th, consistent with the disagreement-triage findings, while 602276 (ACCEPT case) sits at the
34th percentile. Detail: `ood_report.md`, `ood_novelty.py`, `ood_novelty_results.csv`,
`ood_eval_slice_scored.csv`, `ood_demo_scores.json`.

## Cross-cutting takeaway

All four components independently converge on the same read of the 3 example proteins: 602276 is a
clean, trustworthy call at every level (agreement, coverage, novelty score); 169208 is the pipeline's
genuine failure mode (confident-wrong fusion call, flagged by triage, widened by conformal sets,
ranked highest for novelty); and 267317 sits in between — correct but architecturally complex, caught
by WATCH-tier triage and a conformal set that (at high coverage targets) correctly surfaces the
dissenting centroid call rather than hiding it. This consistency across four independently-derived
signals is itself a useful validation that the underlying per-head confidence/agreement information
in the fusion layer, even without any new model training, already carries enough signal to build a
practical reasoning-assistance layer around it.

## Known limitations (repeated from each section for visibility)

1. **Triage rule** is a fixed-weight heuristic (not fit against a review-cost objective) with modest
   standalone AUROC (0.611) — useful for ranking, not a strong binary classifier.
2. **Conformal sets** have a real coverage ceiling (84.4% on this eval slice) driven by
   `novel_family` cases whose true label is never proposed by any head; 90%/95% coverage targets
   cannot be achieved on this data no matter how the threshold is tuned.
3. **OOD/novelty scoring** improvements are real but scoped to the novel_family-vs-novel_seq task
   using only existing scalar per-head signals — genuine embedding-space geometric OOD detection
   (Mahalanobis, energy-based on raw ESM-C vectors) was not possible without the raw embeddings, which
   are not part of this staged data bundle, and the harder open question of detecting non-CAZyme /
   entirely-out-of-taxonomy proteins remains unaddressed.
4. **Demo-protein feature parity**: the 3 example proteins carry a slightly reduced feature set (no
   DIAMOND run, no `contr_cent_margin`/`menu_size` fields) relative to the eval slice, so cross-scale
   comparisons (conformal menus, energy scores) use honestly-reduced versions of each procedure,
   documented in each report.

## Artifact inventory

| File | Contents |
|---|---|
| `report_267317.md`, `report_602276.md`, `report_169208.md` | Grounded per-protein reasoning reports |
| `triage_rule.py`, `triage_report.md`, `triage_eval_slice.csv` | Disagreement-detection rule, demo + eval-slice validation |
| `conformal.py`, `conformal_report.md`, `conformal_calibration_results.csv`, `conformal_demo_predictions.json` | Split-conformal implementation, coverage/set-size results, demo predictions |
| `ood_novelty.py`, `ood_report.md`, `ood_novelty_results.csv`, `ood_eval_slice_scored.csv`, `ood_demo_scores.json` | OOD/novelty scoring implementation, evaluation vs. baseline, demo scores |
| `track_report.md` | This document |
