# dbCAN4-advanced — Prototypes

Two standalone, review-stage prototypes built on top of the dbCAN4-advanced ESM-C tier. They run
entirely from the small **`example_data/`** bundle checked in here — no GPU, no reference index, no
live met compute, no network access. Each is a candidate for later integration into the main
pipeline / BioForge UI (see `SYNTHESIS_REPORT.md` for the integration recommendation).

```
prototypes/
├── viz/                 Training / fine-tuning + embedding visualization
│   ├── build_visualizations.py     regenerates all viz deliverables
│   ├── requirements.txt
│   └── METHODS.md                  per-deliverable data sources + methods
├── reasoning/           Grounded reasoning + conformal / OOD abstention
│   ├── run_example.py              regenerates all reasoning deliverables
│   ├── triage_rule.py              disagreement / review-triage rule
│   ├── conformal.py                split-conformal prediction sets
│   ├── ood_novelty.py              OOD / novelty scoring
│   ├── requirements.txt
│   ├── report_267317.md            grounded per-protein reasoning reports
│   ├── report_602276.md
│   ├── report_169208.md            (the disagreement / failure case)
│   ├── triage_report.md
│   ├── conformal_report.md
│   └── ood_report.md
├── example_data/        Pre-staged inputs (2.1 MB) — everything both scripts need
├── outputs/             Committed reference outputs
│   ├── viz/*.png                   static figure previews (the HTML is git-ignored, see below)
│   └── reasoning/*.csv,*.json       reproduced numeric results
└── SYNTHESIS_REPORT.md  cross-track review + what to integrate into the product
```

## Quick start

```bash
# Visualization track  (Plotly + matplotlib)
cd prototypes/viz
pip install -r requirements.txt
python build_visualizations.py --assets ../example_data --outdir ../outputs/viz
#   -> embedding_explorer.html, training_dashboard.html, calibration_confusion_explorer.html,
#      training_summary.png, reliability_diagram.png

# Reasoning track  (pandas + scikit-learn)
cd ../reasoning
pip install -r requirements.txt
PYTHONPATH="$PWD" python run_example.py --assets ../example_data --outdir ../outputs/reasoning
#   -> triage_eval_slice.csv, conformal_calibration_results.csv, conformal_demo_predictions.json,
#      ood_novelty_results.csv, ood_eval_slice_scored.csv, ood_demo_scores.json
```

Both scripts are deterministic and reproduce the numbers in the reports exactly.

## What the reasoning example prints (verified)

```
--- Disagreement-detection triage ---
AUROC (flag wrong classifier calls, n=4726): 0.6108
  ACCEPT 4437  acc 0.838      WATCH 162  acc 0.784      FLAG 127  acc 0.339
  demo:  602276 -> ACCEPT   267317 -> WATCH   169208 -> FLAG

--- Conformal prediction sets (split-conformal, 50/50 calib/test) ---
  target 70% -> empirical 0.799   80% -> 0.811   90%/95% -> 0.844 (menu-coverage ceiling)

--- OOD / novelty (novel_family vs novel_seq, n=4726) ---
  baseline off-the-shelf centroid_margin  AUROC 0.6548   (matches the pipeline's 0.6549)
  energy score (unsupervised)             AUROC 0.7386
  logistic regression (5-fold CV)         AUROC 0.7856 +/- 0.025
```

The supervised OOD score is a **+0.13 AUROC** gain over the current novelty baseline using only
signals the pipeline already computes — no new embeddings, no retraining. See
`reasoning/ood_report.md` and `SYNTHESIS_REPORT.md`.

## Interactive HTML explorers

The three HTML explorers are ~5–6 MB each (Plotly.js embedded inline so they open offline with no
CDN). They are **regenerable in one command** (above) and are therefore *not* committed to the repo
to keep git history lean — run `build_visualizations.py` to produce them locally. The static PNG
previews in `outputs/viz/` are committed for quick reference.

## Provenance and honesty notes

- `example_data/` is a pre-staged, read-only snapshot of already-computed pipeline outputs
  (reference-embedding UMAP subsample, the n=4,726 2024→2025 temporal-holdout eval predictions, the
  real 7-checkpoint head-training log, and the 3 worked-example evidence bundles). Nothing here was
  re-trained or re-embedded.
- The grounded per-protein reasoning *reports* (`report_*.md`) were drafted with an LLM under a
  strict grounding constraint (every quantitative claim traces to a value in `example_data/`).
  `run_example.py` reproduces only the **deterministic** pieces (triage, conformal, OOD); it does not
  re-run the LLM narrative, which is non-deterministic by nature.
- `novel_family` exact accuracy is 0.0 **by construction** (those families are absent from head
  training) — this is a property of the temporal-holdout design, not a model failure.
