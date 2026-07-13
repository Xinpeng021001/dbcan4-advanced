"""Parse dbCAN CGC output `cgc/<sample>_cgc_standard_out.tsv`.

run_dbcan v4 emits one row per gene belonging to a CAZyme Gene Cluster (CGC):

    CGC#  Gene Type  Contig ID  Protein ID  Gene Start  Gene Stop  Direction  Protein Family

`Gene Type` is CAZyme / TC / TF / STP / null; `Protein ID` is the locus tag that
joins to our genes. Header-keyed (fails loud on an unrecognised layout). The
loader aggregates these per-gene rows into one CazymeCluster per CGC#.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


@dataclass
class CgcMember:
    cluster_key: str
    gene_key: str
    gene_type: str | None   # CAZyme | TC | TF | STP | null
    contig: str | None
    start: int | None
    end: int | None
    family: str | None


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


def _int(row: list[str], i: int | None) -> int | None:
    if i is None or i >= len(row):
        return None
    try:
        return int(row[i].strip())
    except ValueError:
        return None


def parse_cgc_standard(path: str | Path) -> Iterator[CgcMember]:
    path = Path(path)
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return
        idx = {name.strip().lower(): i for i, name in enumerate(header)}
        cgc_i = _find(idx, "cgc")
        pid_i = _first(_find(idx, "protein", "id"), _find(idx, "gene", "id"))
        if cgc_i is None or pid_i is None:
            raise ValueError(
                f"CGC standard-out {path} missing CGC#/Protein ID columns; got {header!r}"
            )
        type_i = _find(idx, "type")
        contig_i = _find(idx, "contig")
        start_i = _find(idx, "start")
        end_i = _first(_find(idx, "stop"), _find(idx, "end"))
        fam_i = _first(_find(idx, "family"), _find(idx, "annotation"))

        for row in reader:
            if not row or len(row) <= max(cgc_i, pid_i):
                continue
            cluster_key = row[cgc_i].strip()
            gene_key = row[pid_i].strip()
            if not cluster_key or not gene_key:
                continue

            def cell(i: int | None) -> str | None:
                if i is None or i >= len(row):
                    return None
                v = row[i].strip()
                return v or None

            yield CgcMember(
                cluster_key=cluster_key,
                gene_key=gene_key,
                gene_type=cell(type_i),
                contig=cell(contig_i),
                start=_int(row, start_i),
                end=_int(row, end_i),
                family=cell(fam_i),
            )
