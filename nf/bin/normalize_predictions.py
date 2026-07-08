#!/usr/bin/env python3
"""Normalize a method's raw prediction output to the standard §2.1 schema.

The dbCAN4-advanced prototype scripts (retrieval_esmc.py, train_heads.py,
diamond_baseline.py) each emit their own wide per-query TSV. The Nextflow
pipeline runs the real tool, then calls this to project that raw TSV onto the
ONE normalized schema the OUTPUT_CONTRACT defines:

    protein_id  family  confidence  ec  all_families  extra(json)

so every downstream consumer (BioForge ingester, web app) reads the same shape
regardless of which method produced the row. Column mappings per method are
declared here in METHOD_MAP — adding a method is a dict entry, not new code.

Usage:
    normalize_predictions.py --tool ESM-C-kNN --in raw.tsv --out ESM-C-kNN.tsv
    normalize_predictions.py --tool fusion    --in fus.tsv --out fusion.tsv

Unknown/abstain values ("-", "", "NA") become "-" in family/confidence.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys

NULL = {"", "-", "na", "none", "nan"}


def _clean(v: str | None) -> str:
    return "-" if v is None or v.strip().lower() in NULL else v.strip()


def _fnum(v: str | None) -> str:
    v = _clean(v)
    if v == "-":
        return "-"
    try:
        return f"{float(v):.4f}"
    except ValueError:
        return "-"


def _coerce(v):
    """Turn a raw TSV cell into a native type for the `extra` json object:
    parse JSON-looking cells (dicts/lists the upstream tool already emitted),
    ints, floats; else keep the string. Prevents double-encoded JSON in `extra`.
    """
    if not isinstance(v, str):
        return v
    s = v.strip()
    if s == "" or s == "-":
        return s
    if s[0] in "{[":
        try:
            return json.loads(s)
        except (json.JSONDecodeError, ValueError):
            pass
    if s.isdigit() or (s[0] == "-" and s[1:].isdigit()):
        try:
            return int(s)
        except ValueError:
            pass
    try:
        return float(s)
    except ValueError:
        return s


def _pick(row, keys):
    """Build an `extra` dict from selected columns, coercing JSON/number cells."""
    return {k: _coerce(row[k]) for k in keys if k in row and row[k] not in ("", "-")}


# For each registry tool: which raw columns carry (id, family, confidence),
# plus a builder for the method-specific `extra` json. Multiple raw layouts
# (e.g. the esmc retrieval TSV feeds BOTH kNN and centroid) are handled by
# pointing several tools at the same raw file with different column picks.
METHOD_MAP = {
    "ESM-C-kNN": dict(
        id="query_id", family="knn_pred", conf="knn_conf",
        extra=lambda r: _pick(r, ("knn_purity", "knn_margin")),
    ),
    "ESM-C-centroid": dict(
        id="query_id", family="cent_pred", conf="cent_conf",
        extra=lambda r: _pick(r, ("cent_margin",)),
    ),
    "ESM-C-contrastive": dict(
        id="query_id", family="clf_pred", conf="clf_conf",
        extra=lambda r: _pick(r, ("contr_cent_pred", "contr_cent_margin",
                                  "contr_knn_pred", "contr_knn_purity")),
    ),
    "Foldseek-CAZyme3D": dict(
        id="protein_id", family="family", conf="prob",
        extra=lambda r: _pick(r, ("target", "tmscore", "lddt")),
    ),
    "SaProt": dict(
        id="protein_id", family="family", conf="cosine",
        extra=lambda r: _pick(r, ("nn_id",)),
    ),
    "fusion": dict(
        id="protein_id", family="family", conf="confidence",
        extra=lambda r: _pick(r, ("votes", "agreement", "signals")),
    ),
}


def normalize(tool: str, in_path: str, out_path: str,
              id_col=None, family_col=None, conf_col=None) -> int:
    spec = METHOD_MAP.get(tool, {})
    id_col = id_col or spec.get("id", "protein_id")
    family_col = family_col or spec.get("family", "family")
    conf_col = conf_col or spec.get("conf", "confidence")
    extra_fn = spec.get("extra", lambda r: {})

    n = 0
    with open(in_path, newline="") as fin, open(out_path, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["protein_id", "family", "confidence", "ec",
                         "all_families", "extra"])
        for row in reader:
            pid = _clean(row.get(id_col))
            if pid == "-":
                continue
            fam = _clean(row.get(family_col))
            # A comma/;-joined multi-family cell: first is primary, rest -> all_families
            all_fams = "-"
            if "," in fam or ";" in fam:
                parts = [p.strip() for p in fam.replace(";", ",").split(",") if p.strip()]
                fam = parts[0] if parts else "-"
                all_fams = ",".join(parts) if parts else "-"
            try:
                extra = extra_fn(row)
            except Exception:
                extra = {}
            writer.writerow([
                pid, fam, _fnum(row.get(conf_col)),
                _clean(row.get("ec")), all_fams,
                json.dumps(extra, separators=(",", ":")) if extra else "{}",
            ])
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True)
    ap.add_argument("--in", dest="in_path", required=True)
    ap.add_argument("--out", dest="out_path", required=True)
    ap.add_argument("--id-col", default=None)
    ap.add_argument("--family-col", default=None)
    ap.add_argument("--conf-col", default=None)
    args = ap.parse_args()
    n = normalize(args.tool, args.in_path, args.out_path,
                  args.id_col, args.family_col, args.conf_col)
    print(f"[normalize] {args.tool}: wrote {n} rows -> {args.out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
