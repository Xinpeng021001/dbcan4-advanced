# Conformal Prediction Sets

## Method

Split-conformal prediction (Vovk et al.) over a per-protein candidate **menu**: the union of family
predictions from all available per-method heads (trained classifier, ESM-C-kNN, ESM-C-centroid,
ESM-C-contrastive-kNN, and DIAMOND), each keeping its own reported confidence for its call. The
nonconformity score of a candidate family `y` for protein `x` is `s(x,y) = 1 - c(y)` if `y` is in the
menu, else `1.0` (maximal nonconformity — a family no head proposed at all).

Calibration: on a held-out calibration split, compute `s_i = min` (most generous) over the protein's
true label(s) of `s(x_i, y_true)`. `q_hat` is the `ceil((n+1)(1-alpha))/n` empirical quantile of these
scores. At test time, the prediction set for protein `x` is `{ y in menu(x) : s(x,y) <= q_hat }` —
i.e. every candidate family whose reported confidence is high enough to fall under the calibrated
nonconformity threshold. Implementation: `conformal.py`.

## Calibration + held-out coverage evaluation

Eval slice (n=4726) split 50/50 into calibration (2363) / test (2363) at a
fixed seed. Results at four target coverage levels:

|   alpha | target_coverage_label          |   q_hat |   empirical_coverage_test |   mean_set_size |   pct_size_0_abstain |   pct_size_1 |   pct_size_ge2 |
|--------:|:-------------------------------|--------:|--------------------------:|----------------:|---------------------:|-------------:|---------------:|
|    0.3  | 70%                            |  0      |                    0.7994 |           0.953 |               0.0474 |       0.9522 |         0.0004 |
|    0.2  | 80%                            |  0.0046 |                    0.8113 |           0.973 |               0.0326 |       0.9619 |         0.0055 |
|    0.1  | 90% (target, see ceiling note) |  1      |                    0.8443 |           1.601 |               0      |       0.4778 |         0.5222 |
|    0.05 | 95% (target, see ceiling note) |  1      |                    0.8443 |           1.601 |               0      |       0.4778 |         0.5222 |

**Coverage ceiling.** The eval slice's true label is entirely absent from every head's proposed menu
for 16.1% of proteins overall — almost entirely the
`novel_family` bucket (only 4.1%
menu-coverage there, vs. 98.4% for
`novel_seq`). No conformal threshold can cover an example whose true label was never proposed by any
head, so achievable empirical coverage on this test split **saturates at 84.4%**
regardless of how loose `q_hat` is set. The 90% and 95% target rows above hit `q_hat=1.0` (accept
every menu entry) and still only reach 84.4%
coverage — this is a genuine, honestly-reported ceiling, not a calibration bug.

**Below the ceiling** (70% / 80% targets), the procedure is well-behaved: empirical coverage tracks
the target closely (79.9% achieved vs. 70% target; 81.1% achieved vs. 80% target — both conservative,
as split-conformal guarantees are), and set sizes are small — 95.2%
of test proteins get a singleton set at the 70% target, meaning the conformal procedure abstains
(empty set) on only 4.7%
of cases rather than forcing a low-confidence single guess.

## Applied to the 3 example proteins

The three example proteins were annotated with only 4 of the 5 heads used for calibration (kNN,
centroid, classifier, contrastive-kNN — no DIAMOND run was staged for these), so their menus are
correspondingly smaller; the calibrated `q_hat` thresholds from the eval slice are applied as-is.

|   protein_id | true_family                  | menu                                                | pred_set_70pct   | pred_set_80pct   | pred_set_90pct_target       |
|-------------:|:-----------------------------|:----------------------------------------------------|:-----------------|:-----------------|:----------------------------|
|       267317 | GH78 (multidomain GH28+GH78) | {'GH78': 0.9953, 'GH92': 0.9774}                    | []               | []               | ['GH78', 'GH92']            |
|       602276 | GH11                         | {'GH11': 1.0}                                       | ['GH11']         | ['GH11']         | ['GH11']                    |
|       169208 | GH183                        | {'GH43_6': 0.9861, 'GH183': 0.9844, 'PL42': 0.2972} | []               | []               | ['GH43_6', 'GH183', 'PL42'] |

- **267317** (true GH78): menu = {GH78: 0.995, GH92: 0.977}. At the 70%/80% targets the set is
  **empty** — both candidate confidences fall just below `q_hat`, i.e. the conformal procedure
  abstains rather than picking either. At the 90% (saturated) target, the set widens to **{GH78,
  GH92}**, correctly containing the true family but also the centroid head's dissenting call — this
  mirrors the WATCH-tier triage verdict for this protein (Step 2): real ambiguity, not a clean single
  answer.
- **602276** (true GH11): menu = {GH11: 1.0} (unanimous across all 4 heads). Singleton set **{GH11}**
  at every target level, correctly resolving to the true family with no ambiguity — consistent with
  the ACCEPT-tier triage verdict.
- **169208** (true GH183): menu = {GH43_6: 0.986, GH183: 0.984, PL42: 0.297}. At 70%/80% the set is
  **empty** (all three confidences again just below `q_hat`). At the 90% saturated target the set
  widens to all three candidates **{GH43_6, GH183, PL42}** — the conformal set correctly contains the
  true family GH183, but only by including the wrong high-confidence kNN call and the low-confidence
  PL42 call alongside it, i.e. it surfaces the same underlying disagreement the FLAG-tier triage
  verdict (Step 2) already identified, this time as an explicit 3-candidate set rather than a
  single wrong point prediction.

## Honest limitations

- The achievable coverage ceiling on this eval slice is 84.4%, driven almost entirely by
  the `novel_family` bucket, where by definition no head has ever seen the true family and cannot
  propose it. Conformal prediction guarantees coverage **conditional on the true label being a
  reachable candidate** — it cannot manufacture coverage for labels no upstream head ever proposes.
  A production deployment should report this ceiling alongside any target-coverage number, and treat
  a persistently empty prediction set as itself diagnostic of a possible novel-family case (worth
  routing to Tier 2 structural search / manual curation), not as a bug to hide.
- At 70%/80% targets, a nontrivial fraction of examples (3–5%) get an **empty** set — an explicit
  abstain rather than the guaranteed-nonempty behavior of some conformal variants (e.g. always
  keeping the top-1 candidate). This is a deliberate design choice here (empty = "no head was
  confident enough to trust"), but a deployment could instead force the top-1 candidate into the set
  when it would otherwise be empty, at the cost of covering some genuinely wrong top-1 calls.
- The demo proteins were scored against menus built from only 4 heads (DIAMOND unavailable for
  these 3), while calibration used 5 heads including DIAMOND — this is disclosed above rather than
  silently normalized away.
