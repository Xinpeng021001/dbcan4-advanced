"""Parse dbCAN substrate prediction `substrate/<sample>_substrate_prediction.tsv`.

Substrate is predicted per CGC (not per gene), so this returns
{cluster_key: (substrate, score)}. Header-keyed: locate the cluster-id column and
the substrate column by name, tolerating the column-order/label drift between
run_dbcan versions.
"""
from __future__ import annotations

import csv
from pathlib import Path


def _find(idx: dict[str, int], *needles: str) -> int | None:
    for name, i in idx.items():
        if all(n in name for n in needles):
            return i
    return None


def _first(*vals: int | None) -> int | None:
    """First non-None index (index 0 is valid, so don't use `a or b`)."""
    for v in vals:
        if v is not None:
            return v
    return None


def parse_substrate(path: str | Path) -> dict[str, tuple[str | None, float | None]]:
    path = Path(path)
    out: dict[str, tuple[str | None, float | None]] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return out
        idx = {name.strip().lower().lstrip("#").strip(): i for i, name in enumerate(header)}
        cid_i = _first(_find(idx, "cluster"), _find(idx, "cgc"))
        sub_i = _find(idx, "substrate")
        score_i = _find(idx, "score")
        if cid_i is None or sub_i is None:
            raise ValueError(
                f"substrate file {path} missing cluster/substrate columns; got {header!r}"
            )
        for row in reader:
            if not row or len(row) <= max(cid_i, sub_i):
                continue
            key = row[cid_i].strip()
            if not key:
                continue
            sub = row[sub_i].strip() or None
            if sub in ("-", "N", "NA"):
                sub = None
            score = None
            if score_i is not None and score_i < len(row):
                try:
                    score = float(row[score_i].strip())
                except ValueError:
                    score = None
            out[key] = (sub, score)
    return out
