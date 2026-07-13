"""GFF3 parser for Prokka output.

Prokka appends the genome FASTA after a `##FASTA` line; we stop there. We parse
the 9 standard GFF columns and the key=value attribute column into a dict. Only
feature rows are yielded; comment/directive lines are skipped. This is a focused
stdlib parser (equivalent to gffutils/BCBio.GFF for our needs) so the ingester
has no heavy dependency and the parsing is transparent and testable.
"""
from __future__ import annotations

import urllib.parse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator


@dataclass
class GffFeature:
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    strand: str | None
    attributes: dict[str, str] = field(default_factory=dict)

    @property
    def gene_key(self) -> str:
        """Stable per-feature id: prefer locus_tag, else ID, else synthesized."""
        return (
            self.attributes.get("locus_tag")
            or self.attributes.get("ID")
            or f"{self.seqid}:{self.start}-{self.end}"
        )

    @property
    def product(self) -> str | None:
        return self.attributes.get("product")


def _parse_attributes(col9: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in col9.strip().split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, _, val = part.partition("=")
        out[key.strip()] = urllib.parse.unquote(val.strip())
    return out


def parse_gff3(path: str | Path) -> Iterator[GffFeature]:
    """Yield features from a GFF3 file, stopping at the ``##FASTA`` section."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                # Malformed row — skip rather than crash the whole ingest.
                continue
            seqid, source, ftype, start, end, _score, strand, _phase, attrs = cols[:9]
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            yield GffFeature(
                seqid=seqid,
                source=source,
                feature_type=ftype,
                start=start_i,
                end=end_i,
                strand=strand if strand in ("+", "-") else None,
                attributes=_parse_attributes(attrs),
            )
