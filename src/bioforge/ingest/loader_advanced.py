"""Load the dbCAN4-advanced STANDARD OUTPUT CONTRACT as a new BioForge release.

`bioforge-ingest-advanced <manifest.json>` reads the pipeline's manifest (see
`dbcan4-advanced/nf/OUTPUT_CONTRACT.md`) and loads its advanced CAZyme calls +
per-protein features into the SAME database the baseline lives in, as a *new
versioned Release* — never overwriting the baseline.

Linkage model (the key design point): advanced calls attach to the genes the
baseline ingest already created. We do NOT create new Sample/Gene rows — we look
up the latest Sample for each manifest sample_key and join predictions onto its
genes by protein_id == gene_key. Because BioForge's `latest_sample_ids` picks the
newest Sample per sample_key (and we add none), the baseline sample stays the
"current" one, its genes stay visible, and the advanced CazymeAnnotation rows —
tagged method_family='advanced' + release_id=<advanced release> — appear
alongside the baseline calls on every gene. Advanced-vs-baseline is then a query
across `method_family`, never a mutation. Idempotent via the same
(sample_key, file_checksum) provenance short-circuit as the baseline loader.
"""
from __future__ import annotations

import hashlib
import shutil
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy import func, select
from sqlalchemy.orm import Session

from ..models import CazymeAnnotation, Gene, ProteinFeature, Release, Sample
from .parse_advanced import (
    AdvancedManifest,
    parse_manifest,
    read_features,
    read_predictions,
)
from .provenance import file_already_ingested, sha256_file


@dataclass
class AdvancedIngestReport:
    release_id: int | None
    release_label: str | None = None
    samples_matched: int = 0
    samples_missing: list[str] = field(default_factory=list)
    advanced_calls_added: int = 0
    features_added: int = 0
    structures_copied: int = 0
    unmatched_calls: int = 0
    unmatched_features: int = 0
    per_tool: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        d = self.__dict__.copy()
        return d


def _genes_by_key(session: Session, sample_id: int) -> dict[str, Gene]:
    genes = session.scalars(select(Gene).where(Gene.sample_id == sample_id)).all()
    return {g.gene_key: g for g in genes}


def _latest_sample(session: Session, sample_key: str) -> Sample | None:
    """The newest Sample row for a sample_key (max release_id) — the baseline
    sample the advanced calls should attach to."""
    stmt = (
        select(Sample)
        .where(Sample.sample_key == sample_key)
        .order_by(Sample.release_id.desc())
        .limit(1)
    )
    return session.scalars(stmt).first()


def _manifest_signature(man: AdvancedManifest) -> tuple[str, dict[str, str]]:
    """Combined sha256 over every file the manifest references (+ the manifest
    itself), so re-ingesting an unchanged manifest is a no-op."""
    per_file: dict[str, str] = {"manifest": sha256_file(man.manifest_path)}
    for s in man.samples:
        for p in s.predictions:
            if p.path.exists():
                per_file[f"pred:{s.sample_key}:{p.tool}"] = sha256_file(p.path)
        for fe in s.features:
            if fe.path.exists():
                per_file[f"feat:{s.sample_key}:{fe.feature_type}"] = sha256_file(fe.path)
    joined = "|".join(f"{k}:{v}" for k, v in sorted(per_file.items()))
    return hashlib.sha256(joined.encode()).hexdigest(), per_file


def ingest_advanced_manifest(
    session: Session,
    manifest_path: str | Path,
    label: str | None = None,
    notes: str | None = None,
    structures_dir: str | Path | None = None,
    make_current: bool = True,
) -> AdvancedIngestReport:
    man = parse_manifest(manifest_path)
    combined, per_file = _manifest_signature(man)

    # Idempotency: manifest-level signature keyed under a synthetic sample_key.
    idem_key = f"__advanced__:{man.release_label}"
    if file_already_ingested(session, idem_key, combined):
        return AdvancedIngestReport(release_id=None, release_label=man.release_label)

    release = Release(
        label=label or man.release_label,
        pipeline_version=man.pipeline_version,
        funcscan_version=None,
        source_dir=str(man.base_dir.resolve()),
        notes=notes or man.release_notes,
        is_current=make_current,
    )
    session.add(release)
    session.flush()

    report = AdvancedIngestReport(release_id=release.id, release_label=release.label)

    # Where to copy served structure PDBs (so the web layer can serve them).
    struct_out = Path(structures_dir) if structures_dir else None

    for sm in man.samples:
        sample = _latest_sample(session, sm.sample_key)
        if sample is None:
            report.samples_missing.append(sm.sample_key)
            continue
        report.samples_matched += 1
        gbk = _genes_by_key(session, sample.id)

        # --- advanced CAZyme predictions ---
        for pf in sm.predictions:
            if not pf.path.exists():
                continue
            n_tool = 0
            for call in read_predictions(pf):
                gene = gbk.get(call.protein_id)
                if gene is None:
                    report.unmatched_calls += 1
                    continue
                session.add(CazymeAnnotation(
                    gene_id=gene.id,
                    cazy_family=call.cazy_family,
                    ec_number=call.ec,
                    tool=call.tool,
                    n_tools_support=None,
                    substrate=None,
                    confidence=call.confidence,
                    method_family="advanced",
                    method_kind=call.method_kind,
                    release_id=release.id,
                    raw={"all_families": call.all_families, **call.extra},
                ))
                report.advanced_calls_added += 1
                n_tool += 1
            report.per_tool[pf.tool] = report.per_tool.get(pf.tool, 0) + n_tool
            _prov(session, release, sm.sample_key, f"advanced:{pf.tool}",
                  sample.id, str(pf.path), per_file.get(f"pred:{sm.sample_key}:{pf.tool}", ""),
                  pf.tool, man.tool_versions.get(pf.tool.split("-")[0].lower()))

        # --- per-protein features (SignalP6 / DeepTMHMM / structure) ---
        for fe in sm.features:
            if not fe.path.exists():
                continue
            for feat in read_features(fe):
                gene = gbk.get(feat.protein_id)
                if gene is None:
                    report.unmatched_features += 1
                    continue
                served_path = None
                if feat.feature_type == "structure" and feat.structure_rel_path:
                    src = (fe.path.parent / feat.structure_rel_path)
                    if src.exists() and struct_out is not None:
                        dst_rel = f"{sm.sample_key}/{feat.protein_id}.pdb"
                        dst = struct_out / dst_rel
                        dst.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copyfile(src, dst)
                        served_path = f"structures/{dst_rel}"
                        report.structures_copied += 1
                    elif src.exists():
                        served_path = str(src)
                session.add(ProteinFeature(
                    gene_id=gene.id,
                    feature_type=feat.feature_type,
                    tool=feat.tool,
                    label=feat.label,
                    score=feat.score,
                    start=feat.start,
                    end=feat.end,
                    structure_path=served_path,
                    attributes=feat.attributes,
                    release_id=release.id,
                ))
                report.features_added += 1
            _prov(session, release, sm.sample_key, f"feature:{fe.feature_type}",
                  sample.id, str(fe.path),
                  per_file.get(f"feat:{sm.sample_key}:{fe.feature_type}", ""),
                  fe.tool, None)

    if report.samples_matched == 0:
        # Nothing attached (e.g. baseline not ingested yet) — don't leave an
        # empty release. Surface the missing keys to the caller.
        session.delete(release)
        session.commit()
        report.release_id = None
        return report

    # Record the manifest-level idempotency provenance row.
    _prov(session, release, idem_key, "advanced_manifest", None,
          str(man.manifest_path), combined, man.pipeline, man.pipeline_version)

    if make_current:
        for other in session.query(Release).filter(Release.id != release.id).all():
            other.is_current = False
    session.commit()
    return report


def _prov(session, release, sample_key, entity_type, entity_id,
          source_file, checksum, tool_name, tool_version) -> None:
    from ..models import Provenance
    session.add(Provenance(
        release_id=release.id,
        sample_key=sample_key,
        entity_type=entity_type,
        entity_id=entity_id,
        source_file=source_file,
        file_checksum=checksum,
        tool_name=tool_name,
        tool_version=tool_version,
    ))
