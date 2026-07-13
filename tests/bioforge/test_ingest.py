"""Ingestion CLI tests: correct row counts + idempotency (the required test)."""
from __future__ import annotations

from sqlalchemy import func, select

from bioforge.ingest.loader import ingest_directory
from bioforge.models import (
    CazymeAnnotation,
    Gene,
    InterproDomain,
    Provenance,
    Release,
    Sample,
)

# Expected counts, hand-derived from sample_data/ (see fixture files):
#   genes:   sampleA 5 features + sampleB 2 = 7
#   cazymes: per (gene, tool, family) row  = 15
#   domains: IPS hits                       = 6
EXPECTED_GENES = 7
EXPECTED_CAZYMES = 15
EXPECTED_DOMAINS = 6


def _counts(session):
    return (
        session.scalar(select(func.count(Gene.id))),
        session.scalar(select(func.count(CazymeAnnotation.id))),
        session.scalar(select(func.count(InterproDomain.id))),
    )


def test_ingest_row_counts(Session, sample_data):
    with Session() as s:
        report = ingest_directory(s, sample_data, label="t1")
    with Session() as s:
        genes, cazymes, domains = _counts(s)
    assert report.samples_added == 2
    assert genes == EXPECTED_GENES
    assert cazymes == EXPECTED_CAZYMES
    assert domains == EXPECTED_DOMAINS
    # Report and DB agree.
    assert report.genes_added == EXPECTED_GENES
    assert report.cazymes_added == EXPECTED_CAZYMES
    assert report.domains_added == EXPECTED_DOMAINS


def test_ingest_is_idempotent(Session, sample_data):
    # First ingest.
    with Session() as s:
        ingest_directory(s, sample_data, label="t1")
    with Session() as s:
        before = _counts(s)
        releases_before = s.scalar(select(func.count(Release.id)))

    # Second ingest of identical output — must add nothing.
    with Session() as s:
        report2 = ingest_directory(s, sample_data, label="t2")

    with Session() as s:
        after = _counts(s)
        releases_after = s.scalar(select(func.count(Release.id)))

    assert report2.samples_added == 0
    assert report2.samples_skipped == 2
    assert report2.release_id is None  # empty release rolled back
    assert before == after  # no duplicated rows
    assert releases_before == releases_after == 1  # no phantom release created


def test_changed_file_creates_new_version(Session, sample_data, tmp_path):
    import shutil

    # Copy fixture so we can mutate one sample's GFF.
    work = tmp_path / "run"
    shutil.copytree(sample_data, work)

    with Session() as s:
        ingest_directory(s, work, label="v1")
    with Session() as s:
        v1_samples = s.scalar(select(func.count(Sample.id)))

    # Change sampleB's GFF (add a gene) -> should create a new sample version.
    gff = work / "results/annotation/prokka/all/sampleB/sampleB.gff"
    text = gff.read_text().replace(
        "##FASTA",
        "contig1\tProdigal:002006\tCDS\t1000\t1200\t.\t+\t0\tID=PROKKA_00003;locus_tag=PROKKA_00003;product=New enzyme\n##FASTA",
    )
    gff.write_text(text)

    with Session() as s:
        report = ingest_directory(s, work, label="v2")

    with Session() as s:
        v2_samples = s.scalar(select(func.count(Sample.id)))
        releases = s.scalar(select(func.count(Release.id)))

    assert report.samples_added == 1  # only sampleB re-ingested
    assert report.samples_skipped == 1  # sampleA unchanged
    assert v2_samples == v1_samples + 1  # old version retained (append, not overwrite)
    assert releases == 2


def test_provenance_recorded(Session, sample_data):
    with Session() as s:
        ingest_directory(s, sample_data, label="t1")
    with Session() as s:
        prov = s.scalars(select(Provenance)).all()
        # Every provenance row has a checksum + release; sample-level rows carry
        # the combined signature used for idempotency.
        assert all(p.file_checksum for p in prov)
        assert {"sample", "cazyme", "interpro"} <= {p.entity_type for p in prov}
        # Tool versions captured from software_versions.yml.
        sample_rows = [p for p in prov if p.entity_type == "sample"]
        assert all(p.tool_version == "1.14.6" for p in sample_rows)  # prokka


def test_baseline_cazymes_stamped_with_method_metadata(Session, sample_data):
    """The v3 schema extension: baseline dbCAN calls must be stamped
    method_family='baseline', a registry-derived method_kind, and the release_id
    of the ingest that produced them — the seam the advanced ingester compares
    against. No baseline call carries a confidence (only advanced calls do)."""
    from bioforge.methods import kind_of

    with Session() as s:
        report = ingest_directory(s, sample_data, label="t1")
    with Session() as s:
        calls = s.scalars(select(CazymeAnnotation)).all()
        assert calls, "expected baseline CAZyme calls"
        # All baseline family + release stamped, no confidence, kind from registry.
        assert all(c.method_family == "baseline" for c in calls)
        assert all(c.release_id == report.release_id for c in calls)
        assert all(c.confidence is None for c in calls)
        assert all(c.method_kind == kind_of(c.tool) for c in calls)
        # Only the known baseline vocabulary is present at this stage.
        assert {c.tool for c in calls} <= {"HMMER", "dbCAN_sub", "DIAMOND"}


def test_protein_features_table_available(Session, sample_data):
    """protein_features is created by the migration/metadata and is queryable
    (empty until the advanced ingester loads SignalP6/DeepTMHMM/structure rows)."""
    from bioforge.models import ProteinFeature

    with Session() as s:
        ingest_directory(s, sample_data, label="t1")
    with Session() as s:
        assert s.scalar(select(func.count(ProteinFeature.id))) == 0
