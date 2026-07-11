#!/usr/bin/env python3
"""Assemble the cazyme_advanced/manifest.json (the OUTPUT_CONTRACT v1.1 entry point).

Called by the Nextflow COLLATE process. Every per-method prediction TSV and every
per-protein feature TSV is STAGED into this task's working directory (Nextflow
`path` inputs), so we enumerate the files present in CWD (or an explicit
--stage-dir) and map each to its registry tool/feature-type. This is robust to
publishDir path resolution: we never scan the publish tree, only staged inputs.

The BioForge ingester reads this manifest and nothing else.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os

# §2.1 CAZyme-prediction files: filename stem -> registry tool key.
PRED_TOOLS = ["ESM-C-kNN", "ESM-C-centroid", "ESM-C-contrastive",
              "Foldseek-CAZyme3D", "SaProt", "fusion"]

# v1.1 per-protein features: filename -> (feature_type, registry tool).
FEATURE_FILES = {
    # signal_peptide + localization tools are resolved at runtime from each TSV's
    # own `method`/`extra` column (SignalP-6.0 + DeepLoc are license-gated; when
    # not installed the pipeline honestly sources both from DeepTMHMM). The label
    # here is the DEFAULT registry tool; the per-row provenance is authoritative.
    "signalp6.tsv":        ("signal_peptide", "SignalP6/DeepTMHMM"),
    "deeptmhmm.tsv":       ("tm_topology",    "DeepTMHMM"),
    "structures.tsv":      ("structure",      "ESMFold"),
    "domains.tsv":         ("domains",        "Pfam/hmmscan"),
    "structure_hits.tsv":  ("structure_hits", "Foldseek-CAZyme3D"),
    "localization.tsv":    ("localization",   "DeepLoc/derived"),
    "physicochem.tsv":     ("physicochem",    "Biopython"),
    "ec_prediction.tsv":   ("ec_prediction",  "CLEAN"),
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="pipeline outdir (for manifest relative paths)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--stage-dir", default=".", help="dir holding staged prediction/feature TSVs (default CWD)")
    ap.add_argument("--release-label", default=None)
    ap.add_argument("--release-notes", default="ESM-C + structure advanced CAZyme calls")
    ap.add_argument("--pipeline-version", default="0.1.0")
    ap.add_argument("--tool-versions", default="{}", help="JSON dict of tool->version")
    args = ap.parse_args()

    stage = args.stage_dir
    present = set(os.listdir(stage))

    # predictions: contract-relative path is predictions/<sample>/<tool>.tsv
    preds = []
    for tool in PRED_TOOLS:
        if f"{tool}.tsv" in present:
            preds.append({"tool": tool,
                          "path": os.path.join("predictions", args.sample, f"{tool}.tsv")})

    # features: contract-relative path is features/<sample>/<fname>
    feats = []
    for fname, (ftype, tool) in FEATURE_FILES.items():
        if fname in present:
            feats.append({"feature_type": ftype, "tool": tool,
                          "path": os.path.join("features", args.sample, fname)})

    label = args.release_label or f"advanced-{dt.date.today().isoformat()}"
    manifest = {
        "contract_version": "1.1",
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
    base = os.path.join(args.outdir, "cazyme_advanced")
    os.makedirs(base, exist_ok=True)
    out = os.path.join(base, "manifest.json")
    with open(out, "w") as fh:
        json.dump(manifest, fh, indent=2)
    # also write a copy in CWD so the process can publish it directly
    with open("manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"[manifest] {out}: {len(preds)} prediction file(s), {len(feats)} feature file(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
