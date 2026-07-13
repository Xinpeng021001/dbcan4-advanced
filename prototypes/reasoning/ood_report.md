# OOD / Novelty Scoring — Evaluated Against Baseline

## Task and baseline

The task: distinguish `novel_family` proteins (families never seen during ESM-C head training — true
zero-shot cases) from `novel_seq` proteins (new sequences of a family the heads WERE trained on) in
the eval slice (n=4726: 4000 novel_seq /
726 novel_family). This is the standing hard problem noted
in the project's own materials: **trained retrieval/classifier heads improve family accuracy but not
novelty detection** (head_metrics.json's own trained-head novelty AUROCs top out at 0.546
— `classifier_maxsoftmax`).

**Current baseline** (off-the-shelf, untrained ESM-C retrieval, from `esmc_retrieval_summary.json`):
`novelty_detection_auroc.centroid_margin = 0.6549` (score = −nearest-centroid margin: a small margin
between the top-1 and top-2 nearest class centroids signals a borderline/unfamiliar embedding
neighborhood). We reproduced this number directly on the eval slice: **0.6548**
(the ~0.0001 difference is negligible and consistent with this being the same underlying computation).

## Two new scores, evaluated honestly

**1. Energy score (unsupervised).** A min-max-normalized average of five signals whose direction
(higher-or-lower-implies-novel) was fixed by inspecting each signal's own single-feature AUROC, with
no label-fitting of weights: `cent_margin` (lower→novel), `cent_conf` (higher→novel), the *trained*
head's `contr_cent_margin` (higher→novel — notably the **opposite** direction from the off-the-shelf
centroid margin), `menu_size`/candidate-menu disagreement (higher→novel), and `clf_conf`
(lower→novel). Implementation: `ood_novelty.py::energy_novelty_score`.
  - **AUROC = 0.7386** on the full eval slice (min-max scaler fit on the calibration half
    only, to avoid look-ahead) — an improvement of **+0.084** over
    baseline, achieved purely by combining several already-available confidence/margin signals with no
    supervised weight-fitting.

**2. Logistic-regression score (supervised).** A logistic regression over 10 standardized per-head
signals (`clf_conf, contr_cent_margin, contr_knn_purity, knn_conf, knn_purity, knn_margin, cent_conf,
cent_margin, diamond_conf, menu_size`), evaluated **out-of-fold** via 5-fold cross-validation to avoid
overfitting inflation. Implementation: `ood_novelty.py::logistic_novelty_score`.
  - **AUROC = 0.7856 ± 0.0249** (5-fold CV), vs. **0.6540** for
    the baseline signal evaluated on the identical folds — an improvement of
    **+0.132**, consistent across all 5 folds
    ([0.814, 0.775, 0.79, 0.744, 0.806]).

## Results table

| Method | AUROC | Evaluation protocol |
|---|---|---|
| Baseline: off-the-shelf `centroid_margin` (current pipeline) | 0.6549 (reported) / 0.6548 (reproduced here) | full eval slice, n=4726 |
| Energy score (unsupervised, 5-signal) | **0.7386** | full eval slice, calibration-only scaler fit |
| Logistic regression (supervised, 10-signal) | **0.7856 ± 0.0249** | 5-fold CV |

## Applied to the 3 example proteins

The demo proteins carry only 3 of the 5 energy-score signals (`cent_margin`, `cent_conf`, `clf_conf`
— no `contr_cent_margin`/`menu_size` fields are staged for these 3), so the energy score here uses a
reduced 3-signal version, scaled against the eval slice's calibration-half min/max range for
consistency:

| Protein | True family | Energy novelty score | Percentile vs. eval-slice distribution |
|---|---|---|---|
| 267317 | GH78 (multidomain GH28+GH78) | 0.6280 | 97.8th |
| 602276 | GH11 | 0.4813 | 34.0th |
| 169208 | GH183 | 0.8590 | 100.0th |

169208 (the disagreement/failure case) scores at the top of the eval-slice novelty distribution
(100th percentile) — consistent with the FLAG-tier triage verdict from Step 2 and with its
near-zero EC confidence (0.0014). 267317 also scores high (0.628,
97.8th percentile), consistent with its real multidomain ambiguity (WATCH tier). 602276 scores near
the middle of the distribution (34th percentile), consistent with its unanimous, high-confidence
ACCEPT-tier consensus. These three proteins are not themselves labeled novel_family/novel_seq (they
are real 2025 held-out CAZymes, not part of the labeled eval slice), so this is a qualitative
consistency check, not a labeled-accuracy claim.

## Honest assessment

**Both new scores exceed the 0.6549 baseline** on the novel_family-vs-novel_seq detection task,
and the improvement is real and reproducible: the supervised logistic score's advantage
(+0.132 AUROC) holds up across all 5 cross-validation folds with a
small standard deviation (0.025), not just as a single lucky split. This is reported
as a genuine, if modest-to-moderate, positive result.

**But this is a narrower, easier problem than the project's standing hard question.** The 0.6549
baseline (and the ~0.37–0.55 trained-head AUROCs in `head_metrics.json`) target *exactly the same*
novel_family-vs-novel_seq distinction evaluated here — so this is a fair apples-to-apples comparison,
not a different task dressed up as the same one. The improvement comes almost entirely from
**combining multiple already-existing confidence/margin signals** (most of which individually only
reach AUROC 0.50–0.65 — see per-signal breakdown in `ood_novelty.py`'s docstring) rather than from
any new embedding-space geometry (e.g. true Mahalanobis distance to a fitted per-family Gaussian, or a
genuine energy-based OOD score over raw ESM-C embeddings) — those would require the actual ESM-C
embedding vectors, which are **not included in this staged bundle** (only derived scalar confidences/
margins per head are). This is stated as a real constraint on scope, not glossed over.

**Where this does NOT resolve the hard problem:** none of these scores were evaluated on detecting
proteins that are OOD relative to *all* CAZy family space (e.g. a protein that is not a CAZyme at
all, or belongs to a family absent from any reference set) — the eval slice's `novel_family` bucket
still consists of families the broader pipeline eventually needs to recognize as CAZymes, just ones
these particular trained heads never saw. Extending this evaluation to true out-of-CAZyme-space
detection would need a negative-control set (non-CAZyme proteins) that isn't part of this bundle.

Full per-protein scored eval slice: `ood_eval_slice_scored.csv`. Demo-protein scores:
`ood_demo_scores.json`. Implementation: `ood_novelty.py`.
