# dbCAN4 Precision Gate: CAZyme vs non-CAZyme (ESM-C embedding classifier)

## Motivation
The temporal-holdout benchmark measured recall/accuracy on proteins already known
to be CAZymes. It never measured **precision against a realistic non-CAZyme
background** — the review's #1 gap. A real proteome is ~97-98% non-CAZyme, so a
tool that looks excellent on a balanced set can still drown true calls in false
positives. This track builds a decoy-backed negative set and measures the gate
under realistic imbalance.

## Data
- **Positives (train):** 337,759 reference-2024 fungal CAZymes (ESM-C 600M, 1152-dim, mean-pooled).
- **Negatives (train):** 112,818 = 97,818 natural non-CAZyme fungal proteins + 15,000 shuffled-domain decoys.
- **Held-out realistic slice:** 32,320 proteins from 5 whole genomes across 5 taxonomic
  classes (Chytridiomycetes, Cryptomycota, Entomophthoromycotina, Kickxellomycotina,
  Microsporidia), **zero genome overlap** with training. Truth tiers from Track A
  multi-evidence tiering: 669 CAZyme (2.07%), 1,815 gray, 29,836 non-CAZyme.
- Zero protein-ID overlap across train positives / train negatives / realistic slice.

## Model
L2-normalized ESM-C embeddings → logistic regression (C=1.0). Saved: `cazyme_gate_model.joblib`.

## Results

### Balanced held-out (positives vs negatives, ~equal)
- AUROC = **0.9872**, AUPRC = **0.9857**

### Realistic imbalance (2.2% CAZyme, gray excluded, n=30,505) — the honest number
- AUROC = **0.9845** but AUPRC = **0.6616**
- The AUROC/AUPRC gap is the finding: high ranking quality, but precision collapses
  at the realistic base rate.

Operating points (CAZyme vs non-CAZyme):

| threshold | n_flagged | precision | recall | F1 |
|-----------|-----------|-----------|--------|-----|
| 0.50 | 2673 | 0.238 | 0.952 | 0.381 |
| 0.70 | 1643 | 0.373 | 0.915 | 0.529 |
| 0.80 | 1237 | 0.464 | 0.858 | 0.602 |
| 0.90 |  810 | 0.572 | 0.692 | 0.626 |
| 0.95 |  493 | 0.692 | 0.510 | 0.587 |
| 0.99 |  110 | 0.836 | 0.138 | 0.236 |

### Gate score cleanly ranks the three evidence tiers
- CAZyme:     mean 0.893, median 0.951
- gray:       mean 0.234, median 0.135
- non-CAZyme: mean 0.137, median 0.059

### Gray-zone adjudication (the gate's practical contribution)
Among 1,815 gray-zone proteins (sequence evidence ambiguous; **no ground-truth label**,
so these are triage suggestions, not validated calls):
- **41.4% (752) score <0.1** → gate suggests likely non-CAZyme
- 2.0% (36) score >0.9 → gate suggests likely CAZyme
- 56.6% (1,027) remain uncertain (0.1-0.9)

## Calibrated open-set triage
Thresholds t_lo=0.866, t_hi=0.994 (t_hi set for ≥90% precision on positive calls):
- **CONFIDENT_CAZYME:** precision 0.909 (recall 7.5%)
- **CONFIDENT_NON:** purity 0.995 (retires 96.8% of the clean slice; 146/669 true CAZymes lost)
- **ABSTAIN → downstream sequence/structure evidence:** 3.0%

## Design conclusion for dbCAN4
The ESM-C embedding gate is a strong **ranking / triage** signal, **not** a standalone
precision classifier at realistic base rates. Its correct role in the dbCAN4 stack:
1. **Retire confident non-CAZymes cheaply.** At a single low threshold t=0.10 on the
   clean cazyme-vs-non slice, proteins scoring <0.10 are 99.99% truly non-CAZyme
   (negative purity), which retires 63.7% of all true non-CAZymes while losing only
   1/669 true CAZymes. Applying the *same* t=0.10 to the gray zone (which has no ground
   truth) flags 41.4% of gray proteins as likely non-CAZyme — a triage suggestion, not a
   validated call. (Coverage rises with the threshold: t=0.20 retires 79% of true
   non-CAZymes at 99.99% purity, losing 3/669 CAZymes.) This cuts downstream
   structure/adjudication workload substantially at near-zero false-negative cost.
2. Provide a calibrated abstention layer feeding the fusion module, **not** replace
   HMMER/DIAMOND sequence evidence. Precision on positive calls must come from combining
   the gate with sequence homology, which is the fusion design already in place.

## Artifacts
- `cazyme_gate_model.joblib` — trained gate (sklearn LogisticRegression + L2 norm)
- `gate_eval_realistic_slice.tsv` — per-protein scores + triage on 32,320 held-out proteins
- `gate_operating_points.tsv` — precision/recall/F1 across thresholds
- `gate_calibration.json` — calibrated abstention thresholds + summary metrics
- `gate_grayzone_adjudication.json` — gray-zone triage counts
- `gate_pr_curve.npz` — full PR curve + raw scores
