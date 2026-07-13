"""Parser for `protein_annotation/interproscan/<sample>_interproscan_faa.tsv`.

InterProScan's TSV format is HEADERLESS and positional (see the InterProScan
docs). Columns 1-11 are always present; 12-15 are optional:

    1  Protein accession      (== our gene_key / locus_tag)
    2  MD5
    3  Sequence length
    4  Analysis               (Pfam, PANTHER, ...)
    5  Signature accession
    6  Signature description
    7  Start
    8  Stop
    9  Score / e-value
    10 Status
    11 Date
    12 InterPro accession     (optional)
    13 InterPro description   (optional)
    14 GO annotations         (optional, '|'-separated)
    15 Pathways               (optional)

We validate on column count (>= 11) and index positionally, since there is no
header to key on. Rows with fewer than 11 columns are skipped.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator


@dataclass
class InterproHit:
    gene_key: str
    analysis: str | None
    signature_acc: str | None
    signature_desc: str | None
    interpro_acc: str | None
    interpro_desc: str | None
    start: int | None
    end: int | None
    evalue: float | None
    go_terms: list[str] = field(default_factory=list)


def _opt(row: list[str], i: int) -> str | None:
    if i < len(row):
        v = row[i].strip()
        if v and v not in ("-", "None"):
            return v
    return None


def _int(row: list[str], i: int) -> int | None:
    v = _opt(row, i)
    try:
        return int(v) if v is not None else None
    except ValueError:
        return None


def _float(row: list[str], i: int) -> float | None:
    v = _opt(row, i)
    try:
        return float(v) if v is not None else None
    except ValueError:
        return None


def parse_interproscan(path: str | Path) -> Iterator[InterproHit]:
    path = Path(path)
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 11:
                continue  # not a valid IPS data row
            gene_key = row[0].strip()
            if not gene_key:
                continue
            go_raw = _opt(row, 13)
            go_terms = [g for g in go_raw.split("|")] if go_raw else []
            yield InterproHit(
                gene_key=gene_key,
                analysis=_opt(row, 3),
                signature_acc=_opt(row, 4),
                signature_desc=_opt(row, 5),
                interpro_acc=_opt(row, 11),
                interpro_desc=_opt(row, 12),
                start=_int(row, 6),
                end=_int(row, 7),
                evalue=_float(row, 8),
                go_terms=go_terms,
            )
