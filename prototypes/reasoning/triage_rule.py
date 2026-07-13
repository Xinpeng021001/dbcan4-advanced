"""
dbCAN4-advanced disagreement-detection / review-triage rule.

Score combines three signals into a single triage_score:
  1. Head disagreement:   1 - agreement_frac   (agreement_frac = # heads/votes matching the fusion
     call, out of the total vote count reported by the fusion layer, e.g. 2/4, 3/4, 4/4)
  2. Fusion uncertainty:   1 - fusion_confidence
  3. Confidence collapse:  max(0, max_single_head_confidence - fusion_confidence), capped at 1.0
     (large collapse = at least one head was very confident but the fused call ended up far less
     confident -- a sign the fusion layer is averaging away a real disagreement rather than
     resolving it)

  triage_score = 0.40*(1 - agreement_frac) + 0.35*(1 - fusion_conf) + 0.25*min(conf_collapse, 1.0)

Tiers:
  score <= 0.15            -> ACCEPT  (fusion call trusted, no action)
  0.15 <  score <= 0.35     -> WATCH   (accept but flag for periodic/spot-check review)
  score  > 0.35             -> FLAG    (send to manual curation before accepting)

This module is deliberately simple/interpretable (not learned) -- it composes signals the fusion
layer already emits (per-method confidence, vote agreement count) plus one derived quantity
(confidence collapse). No new model weights were fit; thresholds were chosen by inspection of the
3 example proteins and checked for face validity against the held-out eval slice (see
triage_eval_slice.csv / triage_report.md).
"""

def compute_triage_score(head_confs: dict, fusion_conf: float, agreement: int, n_total_votes: int = 4) -> dict:
    """
    head_confs: {head_name: confidence} for however many independent per-method confidences are
        available (kNN conf, centroid conf, contrastive/classifier conf, ...).
    fusion_conf: the fusion layer's own reported confidence for its final call.
    agreement: number of votes (including the fusion call itself) that match the fusion family,
        as already reported by the fusion layer (e.g. fusion_raw.tsv 'agreement' column).
    n_total_votes: total number of votes tallied by the fusion layer (denominator for agreement_frac).
    """
    agreement_frac = agreement / float(n_total_votes)
    max_head_conf = max(head_confs.values())
    conf_collapse = max(0.0, max_head_conf - fusion_conf)
    score = 0.40 * (1 - agreement_frac) + 0.35 * (1 - fusion_conf) + 0.25 * min(conf_collapse, 1.0)
    if score > 0.35:
        tier = "FLAG"
    elif score > 0.15:
        tier = "WATCH"
    else:
        tier = "ACCEPT"
    return {
        "agreement_frac": agreement_frac,
        "max_head_conf": max_head_conf,
        "conf_collapse": conf_collapse,
        "triage_score": score,
        "triage_tier": tier,
    }
