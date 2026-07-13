"""ORM models — the 6 core tables from DESIGN.md.

Design principle: append + version, never destructive reload. Each ingest run
creates one `Release`; samples/genes/annotations are inserted linked to it, and
`Provenance` records per-record lineage + source-file checksum. Re-ingesting the
same file is a no-op via the checksum short-circuit in ingest/provenance.py.

Portability: JSON columns use SQLAlchemy's generic JSON so the same schema runs
on both SQLite (dev) and Postgres (prod). Lists (e.g. GO terms) are stored as
JSON arrays for the same reason (avoids Postgres-only ARRAY type).
"""
from __future__ import annotations

from datetime import datetime, timezone

from sqlalchemy import (
    JSON,
    Boolean,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


def _utcnow() -> datetime:
    return datetime.now(timezone.utc)


class Base(DeclarativeBase):
    pass


class Release(Base):
    """One row per ingest invocation — the unit of versioning."""

    __tablename__ = "releases"

    id: Mapped[int] = mapped_column(primary_key=True)
    label: Mapped[str] = mapped_column(String(128), index=True)
    pipeline_version: Mapped[str | None] = mapped_column(String(64))
    funcscan_version: Mapped[str | None] = mapped_column(String(64))
    source_dir: Mapped[str | None] = mapped_column(Text)
    notes: Mapped[str | None] = mapped_column(Text)
    is_current: Mapped[bool] = mapped_column(Boolean, default=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=_utcnow)

    samples: Mapped[list["Sample"]] = relationship(back_populates="release")


class Sample(Base):
    """One row per sample per ingested version."""

    __tablename__ = "samples"
    __table_args__ = (
        UniqueConstraint("sample_key", "release_id", name="uq_sample_release"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    sample_key: Mapped[str] = mapped_column(String(255), index=True)
    release_id: Mapped[int] = mapped_column(ForeignKey("releases.id", ondelete="CASCADE"))
    annotation_tool: Mapped[str | None] = mapped_column(String(32))  # prokka | bakta
    n_contigs: Mapped[int] = mapped_column(Integer, default=0)
    n_genes: Mapped[int] = mapped_column(Integer, default=0)
    sample_metadata: Mapped[dict] = mapped_column("metadata", JSON, default=dict)
    track_path: Mapped[str | None] = mapped_column(Text)  # relative JBrowse track dir
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=_utcnow)

    release: Mapped["Release"] = relationship(back_populates="samples")
    genes: Mapped[list["Gene"]] = relationship(
        back_populates="sample", cascade="all, delete-orphan"
    )


class Gene(Base):
    """One row per feature (CDS/gene/RNA) parsed from the Prokka GFF3."""

    __tablename__ = "genes"
    __table_args__ = (
        Index("ix_genes_locus", "sample_id", "contig", "start", "end"),
        Index("ix_genes_key", "gene_key"),
        Index("ix_genes_product", "product"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id", ondelete="CASCADE"))
    gene_key: Mapped[str] = mapped_column(String(255))  # GFF ID / locus_tag
    contig: Mapped[str] = mapped_column(String(255))    # GFF seqid
    start: Mapped[int] = mapped_column(Integer)
    end: Mapped[int] = mapped_column(Integer)
    strand: Mapped[str | None] = mapped_column(String(1))
    feature_type: Mapped[str] = mapped_column(String(32))
    product: Mapped[str | None] = mapped_column(Text)
    attributes: Mapped[dict] = mapped_column(JSON, default=dict)  # raw GFF col 9
    protein_seq: Mapped[str | None] = mapped_column(Text)  # from <sample>_cleaned.faa
    # CGC membership (dbCAN CAZyme gene cluster); SET NULL so deleting a cluster
    # doesn't cascade into its genes.
    cluster_id: Mapped[int | None] = mapped_column(
        ForeignKey("cazyme_clusters.id", ondelete="SET NULL")
    )

    sample: Mapped["Sample"] = relationship(back_populates="genes")
    cazymes: Mapped[list["CazymeAnnotation"]] = relationship(
        back_populates="gene", cascade="all, delete-orphan"
    )
    domains: Mapped[list["InterproDomain"]] = relationship(
        back_populates="gene", cascade="all, delete-orphan"
    )
    args: Mapped[list["ArgAnnotation"]] = relationship(
        back_populates="gene", cascade="all, delete-orphan"
    )
    go_annotations: Mapped[list["GoAnnotation"]] = relationship(
        back_populates="gene", cascade="all, delete-orphan"
    )
    features: Mapped[list["ProteinFeature"]] = relationship(
        back_populates="gene", cascade="all, delete-orphan"
    )
    cluster: Mapped["CazymeCluster | None"] = relationship(back_populates="genes")

    @hybrid_property
    def length(self) -> int:
        """Feature span in bp (inclusive, GFF-style 1-based coordinates).

        Exposed as a hybrid so templates read ``gene.length`` and queries can
        ``order_by(Gene.length)`` without materialising every row in Python.
        """
        return self.end - self.start + 1

    @length.expression  # type: ignore[no-redef]
    def length(cls):  # noqa: N805
        return cls.end - cls.start + 1

    @property
    def location(self) -> str:
        """Human-readable ``contig:start-end(strand)`` label used across the UI."""
        return f"{self.contig}:{self.start:,}-{self.end:,} ({self.strand or '.'})"


class CazymeAnnotation(Base):
    """One row per gene per predicting tool call.

    Baseline calls come from dbCAN `_overview.tsv` (tool ∈ HMMER/dbCAN_sub/DIAMOND,
    ``method_family='baseline'``). The dbCAN4-advanced module adds protein-language-
    model and structure-based calls (tool ∈ ESM-C-kNN / ESM-C-centroid /
    ESM-C-contrastive / Foldseek-CAZyme3D / SaProt / fusion, ``method_family=
    'advanced'``) each carrying a calibrated ``confidence``. The full tool vocabulary,
    its display names and colours live in ``bioforge.methods`` (one source of truth).

    ``release_id`` records which ingest Release produced the call, so advanced
    predictions load as a *new versioned release* alongside the baseline without
    overwriting it — advanced-vs-baseline is then a query, not a mutation.
    """

    __tablename__ = "cazyme_annotations"
    __table_args__ = (
        Index("ix_cazyme_family", "cazy_family"),
        Index("ix_cazyme_method_family", "method_family"),
        Index("ix_cazyme_release", "release_id"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    gene_id: Mapped[int] = mapped_column(ForeignKey("genes.id", ondelete="CASCADE"))
    cazy_family: Mapped[str] = mapped_column(String(64))  # e.g. GH13, GT2, CBM48
    ec_number: Mapped[str | None] = mapped_column(String(64))
    # Canonical tool key; see bioforge.methods.REGISTRY for the full vocabulary.
    tool: Mapped[str | None] = mapped_column(String(32))
    n_tools_support: Mapped[int | None] = mapped_column(Integer)
    substrate: Mapped[str | None] = mapped_column(Text)
    # --- advanced-method extension (all nullable → baseline rows unaffected) ---
    confidence: Mapped[float | None] = mapped_column(Float)          # 0–1 calibrated
    method_family: Mapped[str] = mapped_column(String(16), default="baseline")
    method_kind: Mapped[str | None] = mapped_column(String(24))      # sequence-plm/…
    release_id: Mapped[int | None] = mapped_column(
        ForeignKey("releases.id", ondelete="CASCADE")
    )
    raw: Mapped[dict] = mapped_column(JSON, default=dict)

    gene: Mapped["Gene"] = relationship(back_populates="cazymes")

    @property
    def is_advanced(self) -> bool:
        return self.method_family == "advanced"


class CazymeCluster(Base):
    """A dbCAN CAZyme Gene Cluster (CGC) — a run of neighbouring CAZyme/TC/TF/STP
    genes, from `cgc/<sample>_cgc_standard_out.tsv`, with a predicted substrate
    (from `substrate/<sample>_substrate_prediction.tsv`). Member genes point back
    via ``Gene.cluster_id``. Mirrors dbCAN-PUL / CAZy's gene-cluster concept."""

    __tablename__ = "cazyme_clusters"
    __table_args__ = (
        Index("ix_cgc_locus", "sample_id", "contig", "start", "end"),
        Index("ix_cgc_key", "cluster_key"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id", ondelete="CASCADE"))
    cluster_key: Mapped[str] = mapped_column(String(64))  # e.g. CGC1
    contig: Mapped[str] = mapped_column(String(255))
    start: Mapped[int] = mapped_column(Integer)
    end: Mapped[int] = mapped_column(Integer)
    composition: Mapped[str | None] = mapped_column(String(128))  # e.g. "CAZyme×3, TF×1"
    n_genes: Mapped[int] = mapped_column(Integer, default=0)
    substrate: Mapped[str | None] = mapped_column(Text)
    substrate_score: Mapped[float | None] = mapped_column(Float)
    raw: Mapped[dict] = mapped_column(JSON, default=dict)

    sample: Mapped["Sample"] = relationship()
    genes: Mapped[list["Gene"]] = relationship(back_populates="cluster")


class InterproDomain(Base):
    """One row per InterProScan domain hit, from `_interproscan_faa.tsv`."""

    __tablename__ = "interpro_domains"
    __table_args__ = (
        Index("ix_interpro_acc", "interpro_acc"),
        Index("ix_interpro_sig", "signature_acc"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    gene_id: Mapped[int] = mapped_column(ForeignKey("genes.id", ondelete="CASCADE"))
    analysis: Mapped[str | None] = mapped_column(String(64))  # Pfam, PANTHER, ...
    signature_acc: Mapped[str | None] = mapped_column(String(128))
    signature_desc: Mapped[str | None] = mapped_column(Text)
    interpro_acc: Mapped[str | None] = mapped_column(String(64))
    interpro_desc: Mapped[str | None] = mapped_column(Text)
    start: Mapped[int | None] = mapped_column(Integer)
    end: Mapped[int | None] = mapped_column(Integer)
    evalue: Mapped[float | None] = mapped_column(Float)
    go_terms: Mapped[list] = mapped_column(JSON, default=list)
    raw: Mapped[dict] = mapped_column(JSON, default=dict)

    gene: Mapped["Gene"] = relationship(back_populates="domains")


class ArgAnnotation(Base):
    """An antimicrobial-resistance gene hit, from funcscan's hAMRonization report.

    Mirrors CARD / AMRFinderPlus / ResFinder content: a resistance gene call with
    its drug class, mechanism, %identity/coverage and source tool/DB. ``gene_id``
    is nullable so a hit that doesn't join a called gene still attaches to the
    sample (counted as unmatched at ingest)."""

    __tablename__ = "arg_annotations"
    __table_args__ = (
        Index("ix_arg_drug_class", "drug_class"),
        Index("ix_arg_symbol", "gene_symbol"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    gene_id: Mapped[int | None] = mapped_column(ForeignKey("genes.id", ondelete="CASCADE"))
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id", ondelete="CASCADE"))
    gene_symbol: Mapped[str | None] = mapped_column(String(128))
    gene_name: Mapped[str | None] = mapped_column(Text)
    drug_class: Mapped[str | None] = mapped_column(String(255))
    resistance_mechanism: Mapped[str | None] = mapped_column(String(255))
    identity: Mapped[float | None] = mapped_column(Float)
    coverage: Mapped[float | None] = mapped_column(Float)
    tool: Mapped[str | None] = mapped_column(String(64))
    reference_db: Mapped[str | None] = mapped_column(String(64))
    accession: Mapped[str | None] = mapped_column(String(64))
    raw: Mapped[dict] = mapped_column(JSON, default=dict)

    gene: Mapped["Gene | None"] = relationship(back_populates="args")
    sample: Mapped["Sample"] = relationship()


class GoAnnotation(Base):
    """A Gene Ontology term on a gene, normalised out of the InterProScan GO column
    so GO can be browsed/faceted (AmiGO/QuickGO style). One row per (gene, GO)."""

    __tablename__ = "go_annotations"
    __table_args__ = (Index("ix_go_id", "go_id"),)

    id: Mapped[int] = mapped_column(primary_key=True)
    gene_id: Mapped[int] = mapped_column(ForeignKey("genes.id", ondelete="CASCADE"))
    go_id: Mapped[str] = mapped_column(String(20))  # GO:0005975
    source: Mapped[str | None] = mapped_column(String(64))  # analysis that assigned it

    gene: Mapped["Gene"] = relationship(back_populates="go_annotations")


class ProteinFeature(Base):
    """A per-protein functional/structural feature beyond CAZy-family calls —
    the "deeper analysis" layer for a (candidate) fungal CAZyme.

    One row per (gene, feature). Deliberately generic so the same table holds
    heterogeneous signals produced by the advanced annotation pipeline:

      feature_type   tool          meaning of the columns
      ------------   ----          -----------------------
      signal_peptide SignalP6      label='SP'/'NO_SP'/'LIPO'/…, score=P(SP),
                                   start/end = cleavage site; attributes: per-class
                                   probabilities, cleavage motif.
      tm_topology    DeepTMHMM     label='TM'/'SP+TM'/'Globular', score=n_helices,
                                   attributes: topology string, helix spans.
      structure      ESMFold/AF    label=source, score=mean pLDDT, structure_path =
                                   served PDB/mmCIF; attributes: model, len, per-…
      disorder/…     (future)      free to extend without a migration.

    ``release_id`` ties the feature to the ingest that produced it (same
    versioning contract as CazymeAnnotation)."""

    __tablename__ = "protein_features"
    __table_args__ = (
        Index("ix_pf_gene", "gene_id"),
        Index("ix_pf_type", "feature_type"),
        Index("ix_pf_release", "release_id"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    gene_id: Mapped[int] = mapped_column(ForeignKey("genes.id", ondelete="CASCADE"))
    feature_type: Mapped[str] = mapped_column(String(32))  # signal_peptide|tm_topology|structure
    tool: Mapped[str | None] = mapped_column(String(48))   # SignalP6 | DeepTMHMM | ESMFold …
    label: Mapped[str | None] = mapped_column(String(128))
    score: Mapped[float | None] = mapped_column(Float)
    start: Mapped[int | None] = mapped_column(Integer)
    end: Mapped[int | None] = mapped_column(Integer)
    structure_path: Mapped[str | None] = mapped_column(Text)  # served PDB/mmCIF (structure rows)
    attributes: Mapped[dict] = mapped_column(JSON, default=dict)
    release_id: Mapped[int | None] = mapped_column(
        ForeignKey("releases.id", ondelete="CASCADE")
    )

    gene: Mapped["Gene"] = relationship(back_populates="features")


class Provenance(Base):
    """Per-record lineage: which source file (+ checksum, tool, version) produced it."""

    __tablename__ = "provenance"
    __table_args__ = (
        Index("ix_prov_entity", "entity_type", "entity_id"),
        Index("ix_prov_checksum", "sample_key", "file_checksum"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    release_id: Mapped[int] = mapped_column(ForeignKey("releases.id", ondelete="CASCADE"))
    sample_key: Mapped[str | None] = mapped_column(String(255))
    entity_type: Mapped[str] = mapped_column(String(32))  # sample|gene|cazyme|interpro
    entity_id: Mapped[int | None] = mapped_column(Integer)
    source_file: Mapped[str] = mapped_column(Text)
    file_checksum: Mapped[str] = mapped_column(String(64))
    tool_name: Mapped[str | None] = mapped_column(String(64))
    tool_version: Mapped[str | None] = mapped_column(String(64))
    ingested_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=_utcnow)
