"""Provenance helpers: file checksums, tool versions, and the idempotency check.

Idempotency strategy (see DESIGN.md): before loading a source file for a sample,
we compute its sha256 and look for an existing Provenance row with the same
(sample_key, file_checksum). If found, that file is unchanged since a prior
ingest and is skipped — so re-running on identical output adds no rows.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import yaml
from sqlalchemy import select
from sqlalchemy.orm import Session

from ..models import Provenance


def sha256_file(path: str | Path, _bufsize: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while chunk := fh.read(_bufsize):
            h.update(chunk)
    return h.hexdigest()


def file_already_ingested(session: Session, sample_key: str, checksum: str) -> bool:
    stmt = select(Provenance.id).where(
        Provenance.sample_key == sample_key,
        Provenance.file_checksum == checksum,
    )
    return session.execute(stmt).first() is not None


def parse_software_versions(results_root: Path) -> dict[str, str]:
    """Flatten pipeline_info/software_versions.yml into {tool_lower: version}."""
    path = results_root / "pipeline_info" / "software_versions.yml"
    versions: dict[str, str] = {}
    if not path.exists():
        return versions
    try:
        data = yaml.safe_load(path.read_text()) or {}
    except yaml.YAMLError:
        return versions
    for _process, tools in data.items():
        if isinstance(tools, dict):
            for tool, ver in tools.items():
                versions[str(tool).lower()] = str(ver)
    return versions


def read_samplesheet_metadata(results_root: Path) -> dict[str, dict]:
    """Return {sample_key: {row columns}} from pipeline_info/samplesheet.valid.csv."""
    import csv

    path = results_root / "pipeline_info" / "samplesheet.valid.csv"
    out: dict[str, dict] = {}
    if not path.exists():
        return out
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            key = row.get("sample") or row.get("sample_id")
            if key:
                out[key] = dict(row)
    return out
