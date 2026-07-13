# dbCAN4-advanced — Training/Fine-Tuning Visualization + Reasoning-Assistance Prototypes

**Synthesis report.** Two independent tracks were built in parallel by delegated sub-agents from a
common, pre-staged data bundle (raw + trained-head reference embeddings, eval-slice predictions,
the 3 example-protein evidence bundles, and the real per-epoch training log). Both are **standalone
prototypes** — no repository or web-UI code was touched, nothing was pushed to GitHub. This report
reviews both tracks against the underlying data and recommends what to integrate next.

---

## What was built

### Track 1 — Visualization

| Deliverable | File |
|---|---|
| Interactive embedding explorer (raw-ESM-C ⇄ trained-head toggle) | `embedding_explorer.html` |
| Training dashboard (loss/val-acc curve, per-family accuracy, margin/purity distributions) | `training_dashboard.html` + `training_summary.png` |
| Calibration + interactive confusion explorer | `calibration_confusion_explorer.html` + `reliability_diagram.png` |
| Methods note | `METHODS.md` |

Built from a stratified subsample of **9,815 reference proteins across 814 families** (raw and
trained-head UMAP projections), the eval slice (n=4,726), and the real 7-checkpoint training log
(loss 9.94→3.86, val-accuracy 0.164→0.966).

### Track 2 — Reasoning assistance

| Deliverable | File(s) |
|---|---|
| Grounded per-protein reasoning reports | `report_267317.md`, `report_602276.md`, `report_169208.md` |
| Disagreement-detection / review-triage rule | `triage_rule.py`, `triage_report.md`, `triage_eval_slice.csv` |
| Conformal prediction sets | `conformal.py`, `conformal_report.md`, `conformal_calibration_results.csv` |
| OOD/novelty scoring, evaluated against baseline | `ood_novelty.py`, `ood_report.md`, `ood_novelty_results.csv` |
| Track report | `track_report.md` |

---

## Verification (what I independently checked before writing this synthesis)

I re-derived several headline numbers directly from the saved CSVs rather than trusting the
sub-agents' prose:

- **Triage AUROC**: recomputed 0.6108 from `triage_eval_slice.csv` against the reported 0.611 — match.
- **Triage tier accuracies**: recomputed ACCEPT 0.838 / FLAG 0.339 — match.
- **OOD baseline reproduction**: `ood_novelty_results.csv` reproduces the existing 0.6549 baseline as
  0.6548 — a faithful, independently-checkable replication before claiming any improvement over it.
- **HTML structure**: both `embedding_explorer.html` and `calibration_confusion_explorer.html` contain
  the claimed Plotly `updatemenus` toggle / confusion-matrix elements and all 3 query protein ids.
- **One nuance, not an error**: the demo-case triage tiers (602276=ACCEPT, 267317=**WATCH**,
  169208=FLAG, reported in `triage_report.md`) use the *real* fusion-layer fields available only for
  the 3 examples; the eval-slice-wide score in `triage_eval_slice.csv` necessarily uses a documented
  *proxy* formula (fusion-confidence fields don't exist at eval-slice scale) — so 267317 shows as
  ACCEPT in that CSV. Both reports state this explicitly; it is not an inconsistency once read against
  the "Honest limitations" sections.

I found no fabricated numbers, no invented family/EC calls, and no overstated claims in either track.

---

## Assessment

**Visualization track** delivers exactly what was asked, grounded in real training and eval data
with no interpolation where data was sparse (the 7-checkpoint training log is plotted as 7 points,
not smoothed into a fabricated curve). The embedding-explorer toggle is the single most compelling
piece — it makes the project's central empirical finding (trained head reshapes the embedding space;
0.269→0.779 top-1 accuracy gain on the independent CAZyme3D validation) visually self-evident. The
calibration explorer surfaces a genuinely useful, previously-unquantified fact: **ECE = 0.146**, driven
by a dominant [0.9,1.0] confidence bin (90% of proteins) that is itself overconfident (mean confidence
0.992 vs. empirical accuracy 0.838) — this is new, actionable information about the pipeline's
calibration that didn't exist as a number before this track ran.

**Reasoning track** is the stronger scientific result of the two. Three findings stand out:

1. **The grounded reasoning reports work as designed.** All three correctly reproduce and explain the
   known evidence (602276 clean consensus, 267317 real multidomain conflict, 169208 the
   confident-but-wrong-fusion failure case) without inventing or overriding any call — including
   surfacing a genuine data-quality bug in 169208's evidence (localization says "Extracellular" despite
   `sp_prediction=NO_SP`) rather than silently smoothing over it.
2. **The disagreement-triage rule is real and reproducible**, not just a plausible-sounding heuristic:
   AUROC 0.611 for ranking wrong-vs-correct, FLAG-tier accuracy 33.9% vs. ACCEPT-tier 83.8% — a rule
   built from signals the fusion layer already emits, no new model trained.
3. **The OOD/novelty result is the most important finding to come out of either track.** The current
   novelty baseline (off-the-shelf ESM-C centroid margin, AUROC 0.6549) was reproduced independently
   at 0.6548, then beaten by combining already-available per-head confidence/margin signals:
   an unsupervised energy score reaches **0.7386**, and a supervised logistic combination reaches
   **0.7856 ± 0.025** (5-fold CV, consistent across all folds — not a lucky split). This is a real
   ~13-point AUROC improvement on the project's single weakest number, achieved with **no new
   embeddings and no retraining** — just a better combination of signals the pipeline already
   computes. The track is honest that this doesn't solve the harder open problem (true
   out-of-CAZy-space / non-CAZyme detection, which needs a negative-control set not in this bundle),
   but on the specific novel_family-vs-novel_seq task it targets, the gain is real.

The conformal-prediction work is honest and useful as a *calibrated-uncertainty communication* layer
even though it cannot exceed an 84.4% coverage ceiling on this data — that ceiling itself is a
correct and important finding (it is driven almost entirely by `novel_family` cases no head ever
proposes, i.e. it is measuring the same underlying novelty-detection gap as the OOD work, from a
different angle) rather than a limitation of the conformal method.

---

## Recommendation: what to integrate into the product

Ranked by value-to-effort, given everything above:

1. **Integrate the logistic-regression OOD/novelty score as the pipeline's novelty signal**, replacing
   or supplementing the current raw centroid-margin score. This is the highest-value, lowest-risk
   change: it reuses signals the pipeline already computes per protein (no new embeddings, no
   retraining), and the +0.13 AUROC gain is cross-validated, not a single-split artifact. Recommended
   next step: re-validate on a genuinely held-out slice (the current result used the same eval slice
   for both fitting and evaluating the logistic weights via CV — a fully out-of-sample check before
   shipping is warranted, and gathering a non-CAZyme negative-control set would let this be tested
   against the harder open question too).
2. **Wire the disagreement-triage rule into the fusion layer's output** as a first-class `review_tier`
   field (ACCEPT/WATCH/FLAG) alongside the existing family call and confidence. This is a small,
   interpretable, no-retraining addition that directly operationalizes the reasoning track's most
   actionable output.
3. **Promote the embedding explorer and calibration/confusion explorer into the BioForge web UI** as
   new views (per-gene toggle and a pipeline-wide diagnostics page respectively). These are
   presentation-layer additions with no modeling risk and high demo value.
4. **Add the grounded per-protein reasoning report as an optional per-gene panel**, generated via
   `host.llm` at ingest or on-demand, gated by the review-tier from (2) so it's surfaced automatically
   for WATCH/FLAG proteins rather than generated for every protein by default (cost control).
5. **Conformal prediction sets** are valuable as an internal QA/monitoring tool (a persistently-empty
   prediction set is itself a novel-family signal worth routing to Tier-2 structural search) but are a
   lower near-term priority for the user-facing product than (1)-(4) — hold as a follow-on.

**What I would not do yet:** ship the training dashboard as a live/recurring artifact — it is valuable
for this one training run but would need wiring into an actual experiment-tracking setup (W&B/
TensorBoard, as flagged in the original enhancement roadmap) to be useful across future retrains,
rather than a one-off HTML snapshot.

---

## Honest scope note

Both tracks worked from a common **pre-staged, read-only data bundle** assembled in-session (not by
either sub-agent) specifically to avoid re-running fragile long-running met jobs inside a sub-agent
context — a lesson carried over from an earlier failure in this project where a domain-extraction
sub-agent hung twice on live met compute. Everything above is a faithful analysis of real,
already-computed pipeline outputs; nothing here required new model training or new embeddings.
