# Disagreement Detection + Review Triage

## Rule

`triage_score = 0.40*(1 - agreement_frac) + 0.35*(1 - fusion_conf) + 0.25*min(conf_collapse, 1.0)`

where:
- **agreement_frac** = fraction of votes (out of the fusion layer's own vote tally) that match the
  fusion's final family call (e.g. 2/4, 3/4, 4/4 as reported in `fusion_raw.tsv`'s `agreement` column)
- **fusion_conf** = the fusion layer's own reported confidence for its final call
- **conf_collapse** = max(0, best single-head confidence − fusion confidence), capped at 1.0 — flags
  cases where at least one head was near-certain but the fused output ended up far less confident,
  a sign the fusion step is averaging away real disagreement rather than resolving it

Tiers: **ACCEPT** (score ≤ 0.15, no action) / **WATCH** (0.15 < score ≤ 0.35, accept but log for
periodic spot-check) / **FLAG** (score > 0.35, route to manual curation before accepting).

This is a fixed, interpretable composition of signals the fusion layer already emits — no new model
was fit. Thresholds were set by inspection of the 3 worked examples and checked for face validity
against the eval slice (below).

## Demonstration on the 3 example proteins

|   protein_id | true_family                  | fusion_call   |   fusion_conf |   agreement |   max_head_conf |   conf_collapse |   triage_score | triage_tier                                                     |
|-------------:|:-----------------------------|:--------------|--------------:|------------:|----------------:|----------------:|---------------:|:----------------------------------------------------------------|
|       267317 | GH78 (multidomain GH28+GH78) | GH78          |        0.6671 |           3 |          0.9953 |          0.3282 |         0.2986 | WATCH: intermediate — accept fusion call but log for spot-check |
|       602276 | GH11                         | GH11          |        0.9823 |           4 |          0.9986 |          0.0163 |         0.0103 | ACCEPT: fusion call trusted                                     |
|       169208 | GH183                        | GH43_6        |        0.3082 |           2 |          0.9861 |          0.6779 |         0.6116 | FLAG: mandatory manual review                                   |

**169208 (true GH183) — FLAGGED, score 0.612.** Agreement is only 2/4 (fusion's own vote tally), and
the fusion call (GH43_6, conf 0.308) sits far below the best individual head (ESM-C kNN, conf 0.986),
a confidence collapse of 0.678 — the single largest of the three proteins. The rule catches this case
even though nothing about the fusion output *looks* superficially wrong (it isn't an abstention, it's
a confident-looking single family). The manual-review payoff is real: the centroid head alone had
already recovered the true family (GH183, conf 0.984), and only the fusion layer's vote-following of
the confident-but-wrong kNN call obscured that.

**602276 (true GH11) — ACCEPTED, score 0.010.** Full 4/4 agreement, fusion confidence 0.982, and a
tiny confidence collapse (0.016, since the top head, the contrastive classifier, was itself the
fusion's own basis). No manual review warranted — this is the case the pipeline is designed to
resolve outright.

**267317 (true GH78, real multidomain complexity) — WATCH, score 0.299.** This sits between the other
two: agreement is 3/4 (only the centroid head dissents, calling GH92), fusion confidence is moderate
(0.667), and there's a real confidence collapse (0.328, since two heads independently reported
>0.94 confidence for GH78 while fusion landed at 0.667). The rule places this in the intermediate WATCH
tier rather than FLAG — the fusion call is directionally right and 2/3 heads directly support it, but
the reduced fusion confidence combined with the dissenting centroid head, together with the underlying
Pfam evidence for two distinct domains (PF00295 GH28-like N-terminal, PF17389 GH78-diagnostic
C-terminal — see `report_267317.md`), means this protein is exactly the kind of case worth periodic
audit even though it doesn't require blocking manual review before use. This demonstrates the rule
correctly differentiating "confidently wrong" (169208, FLAG) from "correct but architecturally
complex" (267317, WATCH) from "clean consensus" (602276, ACCEPT).

## Evaluation on the eval slice (n=4726)

Applying the same score formula (using the trained classifier's own confidence as the closest
available proxy for a "fused" decision, since a fusion-layer confidence field is not present in the
eval-slice bundle) to the full held-out evaluation set of 4726 proteins:

| triage_tier_v2   |    n |   accuracy |
|:-----------------|-----:|-----------:|
| ACCEPT           | 4437 |   0.838179 |
| FLAG             |  127 |   0.338583 |
| WATCH            |  162 |   0.783951 |

FLAG-tier accuracy (33.9%) is far
below ACCEPT-tier accuracy (83.8%),
confirming the score is picking out genuinely harder cases rather than flagging at random. As a
continuous ranking signal, `triage_score` achieves an AUROC of **0.611** for distinguishing wrong
vs. correct classifier calls on this slice — a modest but real discriminative signal (0.5 = random),
consistent with disagreement/confidence-collapse being informative but not a strong standalone
predictor of correctness at this scale. The full per-protein scored table is saved as
`triage_eval_slice.csv`.

## Honest limitations

- The eval-slice bundle lacks a genuine fusion-layer confidence/vote-agreement field (that
  infrastructure exists only for the 3 worked examples in `fusion_raw.tsv`), so the eval-slice version
  of the rule substitutes the trained classifier's own confidence and a 3-head vote among
  clf/centroid/contrastive-kNN. This is a reasonable proxy but not an exact replay of the demo-case
  formula — the two should be read as the same *design*, evaluated on the best data available at
  each scale, not as numerically identical computations.
- AUROC 0.611 means the score is a useful *ranking* signal for review triage (send the worst-scoring
  ~3% to a human, as the demo shows), not a strong binary correct/incorrect classifier on its own.
- Thresholds (0.15 / 0.35) were chosen by inspection of the 3 worked examples, not fit by optimizing
  precision/recall on the eval slice; a production system would want to tune them against a labeled
  review-cost/benefit curve.
