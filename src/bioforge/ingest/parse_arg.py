"""Parse funcscan's hAMRonization combined ARG report.

hAMRonization harmonises the output of many AMR tools (AMRFinderPlus, RGI/CARD,
ResFinder, …) into one TSV with a shared column vocabulary. This is a run-level
file covering all samples; each row is one resistance-gene hit. Header-keyed and
tolerant: we map the harmonised column names we care about and skip the rest,
grouping hits by sample.

Returns {sample_key: [ArgHit, …]}.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ArgHit:
    sample_key: str
    gene_key: str | None       # protein/locus id to join to a Gene, if present
    gene_symbol: str | None
    gene_name: str | None
    drug_class: str | None
    resistance_mechanism: str | None
    identity: float | None
    coverage: float | None
    tool: str | None
    reference_db: str | None
    accession: str | None


def _find(idx: dict[str, int], *cands: str) -> int | None:
    """First column whose name equals or contains one of the candidate strings."""
    for c in cands:
        if c in idx:
            return idx[c]
    for name, i in idx.items():
        if any(c in name for c in cands):
            return i
    return None


def _num(row: list[str], i: int | None) -> float | None:
    if i is None or i >= len(row):
        return None
    v = row[i].strip().rstrip("%")
    try:
        return float(v)
    except ValueError:
        return None


def parse_arg_report(path: str | Path) -> dict[str, list[ArgHit]]:
    path = Path(path)
    out: dict[str, list[ArgHit]] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return out
        idx = {name.strip().lower(): i for i, name in enumerate(header)}
        sample_i = _find(idx, "input_file_name", "sample", "input_file")
        if sample_i is None:
            raise ValueError(
                f"hAMRonization report {path} has no sample column; got {header!r}"
            )
        pid_i = _find(idx, "input_protein_id", "input_gene_id", "protein_id", "input_sequence_id")
        sym_i = _find(idx, "gene_symbol")
        name_i = _find(idx, "gene_name")
        drug_i = _find(idx, "drug_class")
        mech_i = _find(idx, "resistance_mechanism")
        id_i = _find(idx, "sequence_identity", "identity")
        cov_i = _find(idx, "coverage_percentage", "coverage")
        tool_i = _find(idx, "analysis_software_name", "tool")
        db_i = _find(idx, "reference_database_name", "reference_database", "database")
        acc_i = _find(idx, "reference_accession", "accession")

        def cell(row: list[str], i: int | None) -> str | None:
            if i is None or i >= len(row):
                return None
            v = row[i].strip()
            return v or None

        for row in reader:
            if not row or len(row) <= sample_i:
                continue
            sample_key = row[sample_i].strip()
            if not sample_key:
                continue
            # Sample column is often a filename like "sampleA.tsv"; take the stem.
            sample_key = Path(sample_key).name.split(".")[0]
            hit = ArgHit(
                sample_key=sample_key,
                gene_key=cell(row, pid_i),
                gene_symbol=cell(row, sym_i),
                gene_name=cell(row, name_i),
                drug_class=cell(row, drug_i),
                resistance_mechanism=cell(row, mech_i),
                identity=_num(row, id_i),
                coverage=_num(row, cov_i),
                tool=cell(row, tool_i),
                reference_db=cell(row, db_i),
                accession=cell(row, acc_i),
            )
            out.setdefault(sample_key, []).append(hit)
    return out
