#!/usr/bin/env python3
"""Fusion layer — combine per-method normalized predictions into one call.

Reads several standard §2.1 TSVs (ESM-C-kNN, ESM-C-centroid, ESM-C-contrastive,
Foldseek-CAZyme3D, SaProt) and, per protein, produces a single fused CAZy-family
call with a calibrated confidence. This is real, GPU-free logic — it runs
identically in `-stub-run` and in production.

Rule (design_dbcan4_advanced.md §3.3): a transparent weighted vote. Each method
votes its predicted family with weight = method_weight × its own confidence.
The winning family's score is normalized by the total possible weight to give a
[0,1] fused confidence. Agreement across *orthogonal* signals (sequence pLM +
structure) is rewarded: we record how many distinct signal types backed the call
and which methods voted, in `extra`, so the DB/UI can show the consensus.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict

# Which evidence signal each method contributes (for orthogonality bookkeeping).
SIGNAL = {
    "ESM-C-kNN": "sequence", "ESM-C-centroid": "sequence",
    "ESM-C-contrastive": "sequence",
    "Foldseek-CAZyme3D": "structure", "SaProt": "structure",
}
# Default per-method weights (structure signals weighted a touch higher because
# they are orthogonal to the sequence baseline that already ran).
DEFAULT_WEIGHTS = {
    "ESM-C-kNN": 1.0, "ESM-C-centroid": 1.0, "ESM-C-contrastive": 1.2,
    "Foldseek-CAZyme3D": 1.5, "SaProt": 1.3,
}


def _tool_from_name(path: str) -> str:
    base = path.split("/")[-1]
    return base[:-4] if base.endswith(".tsv") else base


def load(path):
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            yield row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True,
                    help="normalized §2.1 TSVs (filename stem = registry tool key)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--weights", default="", help="JSON dict tool->weight (optional)")
    ap.add_argument("--min-confidence", type=float, default=0.0)
    args = ap.parse_args()

    weights = dict(DEFAULT_WEIGHTS)
    if args.weights:
        try:
            weights.update(json.loads(args.weights))
        except (json.JSONDecodeError, TypeError):
            pass

    # per protein: family -> summed weighted score; and bookkeeping
    scores: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    votes: dict[str, dict[str, str]] = defaultdict(dict)
    signals: dict[str, set] = defaultdict(set)
    maxposs: dict[str, float] = defaultdict(float)

    for path in args.inputs:
        tool = _tool_from_name(path)
        w = weights.get(tool, 1.0)
        for row in load(path):
            pid = row["protein_id"]
            fam = row.get("family", "-")
            conf_s = row.get("confidence", "-")
            try:
                conf = float(conf_s) if conf_s not in ("-", "", None) else 0.0
            except ValueError:
                conf = 0.0
            # every present method contributes to the max possible weight
            maxposs[pid] += w
            if fam in ("-", "", None):
                continue
            scores[pid][fam] += w * conf
            votes[pid][tool] = fam
            if tool in SIGNAL:
                signals[pid].add(SIGNAL[tool])

    with open(args.out, "w", newline="") as fout:
        wri = csv.writer(fout, delimiter="\t")
        wri.writerow(["protein_id", "family", "confidence", "ec",
                      "all_families", "votes", "agreement", "signals"])
        for pid in sorted(maxposs):
            fam_scores = scores.get(pid, {})
            if not fam_scores:
                wri.writerow([pid, "-", "-", "-", "-",
                              json.dumps(votes.get(pid, {}), separators=(",", ":")),
                              0, json.dumps(sorted(signals.get(pid, [])))])
                continue
            best_fam = max(fam_scores, key=fam_scores.get)
            fused_conf = fam_scores[best_fam] / maxposs[pid] if maxposs[pid] else 0.0
            agreement = sum(1 for t, f in votes[pid].items() if f == best_fam)
            all_fams = ",".join(sorted(fam_scores, key=fam_scores.get, reverse=True))
            if fused_conf < args.min_confidence:
                best_fam = "-"
            wri.writerow([
                pid, best_fam, f"{fused_conf:.4f}", "-", all_fams,
                json.dumps(votes[pid], separators=(",", ":")),
                agreement, json.dumps(sorted(signals.get(pid, []))),
            ])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
