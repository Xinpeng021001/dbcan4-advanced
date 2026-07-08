#!/usr/bin/env python3
"""Assemble the cazyme_advanced/manifest.json (the OUTPUT_CONTRACT entry point).

Called by the Nextflow COLLATE process after per-method normalized TSVs and
feature TSVs have been published. It scans the published tree for a sample and
writes a manifest describing exactly what exists, mapping each file to its
registry tool. The BioForge ingester reads this manifest and nothing else.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os

# Registry tool keys that emit §2.1 CAZyme predictions.
PRED_TOOLS = ["ESM-C-kNN", "ESM-C-centroid", "ESM-C-contrastive",
              "Foldseek-CAZyme3D", "SaProt", "fusion"]
FEATURES = [
    ("signal_peptide", "SignalP6", "signalp6.tsv"),
    ("tm_topology", "DeepTMHMM", "deeptmhmm.tsv"),
    ("structure", "ESMFold", "structures.tsv"),
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="pipeline outdir (contains cazyme_advanced/)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--release-label", default=None)
    ap.add_argument("--release-notes", default="ESM-C + structure advanced CAZyme calls")
    ap.add_argument("--pipeline-version", default="0.1.0")
    ap.add_argument("--tool-versions", default="{}", help="JSON dict of tool->version")
    args = ap.parse_args()

    base = os.path.join(args.outdir, "cazyme_advanced")
    pred_dir = os.path.join(base, "predictions", args.sample)
    feat_dir = os.path.join(base, "features", args.sample)

    preds = []
    for tool in PRED_TOOLS:
        rel = os.path.join("predictions", args.sample, f"{tool}.tsv")
        if os.path.exists(os.path.join(base, rel)):
            preds.append({"tool": tool, "path": rel})

    feats = []
    for ftype, tool, fname in FEATURES:
        rel = os.path.join("features", args.sample, fname)
        if os.path.exists(os.path.join(base, rel)):
            feats.append({"feature_type": ftype, "tool": tool, "path": rel})

    label = args.release_label or f"advanced-{dt.date.today().isoformat()}"
    manifest = {
        "contract_version": "1.0",
        "pipeline": "dbcan4-advanced",
        "pipeline_version": args.pipeline_version,
        "created": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "release_label": label,
        "release_notes": args.release_notes,
        "tool_versions": json.loads(args.tool_versions),
        "samples": [
            {"sample_key": args.sample,
             "cazyme_predictions": preds,
             "protein_features": feats}
        ],
    }
    out = os.path.join(base, "manifest.json")
    os.makedirs(base, exist_ok=True)
    with open(out, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"[manifest] {out}: {len(preds)} prediction file(s), {len(feats)} feature file(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
