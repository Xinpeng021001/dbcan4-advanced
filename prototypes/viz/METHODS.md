# dbCAN4-advanced — Viz-track methods note

Standalone visualization prototype for the ESM-C protein-language-model tier of dbCAN4-advanced,
built from a pre-staged, self-contained data bundle (`viz_track_assets.tar.gz`). No git repo access,
no live compute — every number below is read directly from the staged files; nothing is invented or
interpolated beyond what is stated.

All three interactive deliverables are self-contained HTML (Plotly.js embedded inline, no CDN /
external network calls at view time) and open directly in a browser.

## 1. Interactive embedding explorer — `embedding_explorer.html`

**Source:** `embedding_umap_coords.json`.

**What it shows:** a stratified reference subsample of 9,815 proteins spanning 814 fungal CAZyme
families (capped at 22 proteins/family), with 2D UMAP coordinates precomputed in two spaces:

- **Raw ESM-C** — the 1152-dim mean-pooled ESM-C embedding, cosine-metric UMAP (`raw_umap_x/y`).
- **Trained-head projection** — the contrastive head's 256-dim output (1152→1024→256 MLP),
  L2-normalized, cosine-metric UMAP (`trained_umap_x/y`).

A button toggle (top-left) switches the entire scatter between the two spaces. Points are colored by
family; the 24 most numerous families in the subsample are given distinct colors and shown in the
legend, all others are grouped as `other` (light grey, drawn behind) to keep the legend legible against
814 total families. Hover shows protein id and full family string (families can be multi-domain, e.g.
`CBM1,GH5_7`).

The 3 named query proteins (`169208`/GH183, `602276`/GH11, `267317`/GH78,GH28) are overlaid as black
star markers with direct `Q1`/`Q2`/`Q3` text labels, using `query_raw_umap_x/y` and
`query_trained_umap_x/y`.

**Regenerate:** load `embedding_umap_coords.json`, split `ref_ids`/`ref_fams` into a DataFrame with
`raw_x/y` and `trained_x/y` columns, take `value_counts()` on `fam` for the top-24 grouping, build one
`Scattergl` trace per family per view (50 traces total: 25 raw + 25 trained, each holding all family
groups plus the query overlay), and wire a `updatemenus` button pair that toggles trace `visible`
between the two 25-trace blocks.

## 2. Training dashboard — `training_dashboard.html` (+ static `training_summary.png`)

**Sources:** `train_heads.log` (training curve), `head_eval_pred.tsv` (per-family accuracy, margin/purity
distributions), `head_metrics.json` and `esmc_retrieval_summary.json` (static-figure comparison bars).

Four panels:

1. **Loss & val-accuracy vs epoch** — the training log contains exactly **7 logged epochs**
   (0, 5, 10, 15, 20, 25, 29) out of 30 total; the curve is a dual-axis line plot through exactly these
   7 points and is not interpolated to a finer grid. Loss falls from 9.9417 to 3.8565; val classifier
   accuracy rises from 0.1642 to 0.9655.
2. **Per-family classifier accuracy** — from `head_eval_pred.tsv`, grouped by the first token of
   `true_families` (multi-label rows like `CBM91,GH43_14` are grouped under their first family to avoid
   combinatorial explosion), restricted to families with **≥5 eval proteins** (177 of 269 families with
   any eval support), bars colored by novelty bucket (`novel_seq` / `novel_family` / `mixed`, the last
   meaning the family's eval proteins span both buckets).
3. **Contrastive-centroid margin distribution**, split by whether `contr_cent_pred` matches any token in
   `true_families` (overlap-match).
4. **Contrastive-kNN purity distribution**, split the same way using `contr_knn_pred`.

The static companion PNG (`figure-style`-compliant) puts the training curve beside a bar comparison of
overall exact/overlap accuracy for the untrained ESM-C kNN-retrieval baseline (0.726/0.794, from
`esmc_retrieval_summary.json`, `schemes.knn.overall["thr=0.0"]`) against the three trained-head schemes
(contrastive kNN 0.760/0.830, contrastive centroid 0.755/0.822, classifier 0.754/0.823, all from
`head_metrics.json`).

**Regenerate:** parse `train_heads.log` with the regex `epoch\s+(\d+)\s+loss\s+([\d.]+)\s+val_clf_acc\s+([\d.]+)`;
load `head_eval_pred.tsv`, derive `true_primary = true_families.split(',')[0]` and an overlap-correctness
boolean per scheme (`pred in set(true_families.split(','))`); group/filter/plot as above.

## 3. Calibration + confusion explorer — `calibration_confusion_explorer.html` (+ static `reliability_diagram.png`)

**Source:** `head_eval_pred.tsv` (all 4,726 rows).

**Reliability diagram:** `clf_conf` binned into 10 fixed-width bins [0,1]; empirical accuracy per bin is
the overlap-match rate (`clf_pred` in the true-family token set) of the proteins in that bin. Overall
Expected Calibration Error (weighted mean |confidence − accuracy| across bins) = **0.146**. The dominant
[0.9, 1.0] confidence bin holds 4,253 of 4,726 proteins (90%) and is itself overconfident (mean confidence
0.992 vs. empirical accuracy 0.838) — this single bin drives most of the aggregate ECE.

**Confusion matrix:** family × family heatmap (log10(n+1) color scale) restricted to true families with
**≥15 eval proteins** (97 of 269 families with any eval support, covering 3,866 of 4,726 proteins);
predictions to a family below that threshold are bucketed as `other` on the column axis. Clicking any
cell lists the protein ids in that (true, predicted) bucket (capped at 50 ids/cell) in a side panel, each
annotated with its novelty bucket (`novel_seq`/`novel_family`) and full true-family string.

**Regenerate:** build the confusion matrix with `pandas.crosstab(true_primary, pred_bucket)` after the
same support-threshold filter; export a `{true_fam}|||{pred_fam}: [protein_ids]` JSON lookup for every
non-empty cell and wire a `plotly_click` event listener in vanilla JS (no framework) to render the
looked-up list — see the embedded `<script>` block in the HTML file for the exact implementation.

## Known caveats / things NOT fabricated

- `train_heads.log` has only 7 logged epochs — no finer per-step or per-batch curve exists in the bundle.
- Novel-family exact accuracy is 0.0 for all trained-head schemes by construction (evaluating on families
  absent from training); the per-family accuracy panel labels these rows `novel_family` rather than
  presenting them as an undifferentiated failure.
- Multi-label `true_families` strings (e.g. `GH78,GH28`) are reduced to their first token for grouping in
  the per-family accuracy chart and confusion matrix; `overlap`-style correctness elsewhere in the bundle
  already accounts for the full label set.
- The commonly-cited "novelty AUROC ~0.63" project figure refers to `esmc_retrieval_summary.json`'s
  `novelty_detection_auroc.centroid_margin = 0.6549` (untrained ESM-C retrieval), which is higher than
  every novelty-detection signal computed on the trained heads in `head_metrics.json` — trained heads
  improve family-assignment accuracy but not OOD/novelty detection. This dashboard does not re-derive or
  re-plot that number since no novelty-score-vs-label file was included in this bundle for the trained
  heads' ROC curve itself — only `head_metrics.json`'s final scalar AUROC values are available, which are
  quoted directly.
- `eval2025_overview.tsv` (dbCAN3 baseline overview: HMMER/dbCAN_sub/DIAMOND calls) was inspected but not
  used in any deliverable in this round — no chart in the plan required joining it against
  `head_eval_pred.tsv`, and doing so wasn't necessary to satisfy the requested visualizations.

## Files in this deliverable set

| File | Type | Description |
|---|---|---|
| `embedding_explorer.html` | interactive | Deliverable 1 |
| `training_dashboard.html` | interactive | Deliverable 2 |
| `training_summary.png` | static | Figure-style companion to deliverable 2 |
| `calibration_confusion_explorer.html` | interactive | Deliverable 3 |
| `reliability_diagram.png` | static | Figure-style companion to deliverable 3 |
| `METHODS.md` | markdown | This note |
