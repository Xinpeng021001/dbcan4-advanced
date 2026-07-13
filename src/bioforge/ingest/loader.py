"""Core ETL: load discovered sample inputs into the database with provenance.

Idempotency model (robust + simple): for each sample we compute a combined
signature = sha256 over its source-file checksums (GFF + dbCAN + IPS). If a
Provenance 'sample' row already exists for (sample_key, combined_signature), the
sample is unchanged since a prior ingest and is skipped entirely — so a re-run on
identical output adds zero rows. If any source file changed, a NEW Sample version
is created under the current Release (append + version, never overwrite history).

All of a sample's genes/annotations therefore live in one Release, keeping each
version internally consistent.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy.orm import Session

from ..methods import kind_of as method_kind_of
from ..models import (
    ArgAnnotation,
    CazymeAnnotation,
    CazymeCluster,
    Gene,
    GoAnnotation,
    InterproDomain,
    Provenance,
    Release,
    Sample,
)
from .discover import SampleInputs, _results_root, discover_samples
from .jbrowse import prepare_sample_track
from .parse_arg import parse_arg_report
from .parse_cgc import parse_cgc_standard
from .parse_dbcan import parse_dbcan_overview
from .parse_fasta import parse_protein_fasta
from .parse_gff import parse_gff3
from .parse_interpro import parse_interproscan
from .parse_substrate import parse_substrate
from .provenance import (
    file_already_ingested,
    parse_software_versions,
    read_samplesheet_metadata,
    sha256_file,
)


@dataclass
class IngestReport:
    release_id: int | None
    samples_added: int = 0
    samples_skipped: int = 0
    genes_added: int = 0
    cazymes_added: int = 0
    domains_added: int = 0
    sequences_added: int = 0
    clusters_added: int = 0
    args_added: int = 0
    go_added: int = 0
    unmatched_cazyme: int = 0
    unmatched_domains: int = 0
    unmatched_arg: int = 0

    def as_dict(self) -> dict:
        return self.__dict__.copy()


def _combined_signature(inputs: SampleInputs) -> tuple[str, dict[str, str]]:
    """Return (combined_sha, {role: file_checksum}) for the sample's sources.

    Every source file a sample draws from is hashed in, so a change to any of
    them (a new CGC call, an added sequence, …) yields a new combined signature
    and therefore a new versioned Sample on re-ingest.
    """
    per_file: dict[str, str] = {}
    for role, path in (
        ("gff", inputs.gff_path),          # OPTIONAL: only when coordinates were provided
        ("dbcan", inputs.dbcan_path),
        ("interpro", inputs.interpro_path),
        ("faa", inputs.faa_path),
        ("cgc", inputs.cgc_path),
        ("cgc_std", inputs.cgc_std_path),
        ("substrate", inputs.substrate_path),
        ("arg", inputs.arg_path),
    ):
        if path:
            per_file[role] = sha256_file(path)
    joined = "|".join(f"{k}:{v}" for k, v in sorted(per_file.items()))
    return hashlib.sha256(joined.encode()).hexdigest(), per_file


def ingest_directory(
    session: Session,
    pipeline_output: str | Path,
    label: str | None = None,
    notes: str | None = None,
    tracks_dir: str | Path | None = None,
) -> IngestReport:
    root = _results_root(pipeline_output)
    samples = discover_samples(pipeline_output)
    versions = parse_software_versions(root)
    metadata = read_samplesheet_metadata(root)

    label = label or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    release = Release(
        label=label,
        pipeline_version=versions.get("nf-core/funcscan") or versions.get("nextflow"),
        funcscan_version=versions.get("nf-core/funcscan"),
        source_dir=str(Path(pipeline_output).resolve()),
        notes=notes,
        is_current=True,
    )
    session.add(release)
    session.flush()

    report = IngestReport(release_id=release.id)

    for inp in samples:
        combined, per_file = _combined_signature(inp)
        if file_already_ingested(session, inp.sample_key, combined):
            report.samples_skipped += 1
            continue
        _load_sample(session, release, inp, combined, per_file, versions,
                     metadata, report, tracks_dir)

    if report.samples_added == 0:
        # Nothing new — don't leave an empty release lying around.
        session.delete(release)
        session.commit()
        report.release_id = None
        return report

    # Mark this the current release; demote older ones (display convenience).
    for other in session.query(Release).filter(Release.id != release.id).all():
        other.is_current = False
    session.commit()
    return report


def _load_sample(
    session: Session,
    release: Release,
    inp: SampleInputs,
    combined: str,
    per_file: dict[str, str],
    versions: dict[str, str],
    metadata: dict[str, dict],
    report: IngestReport,
    tracks_dir: str | Path | None = None,
) -> None:
    # dbCAN4 is protein-input by default (no genome, no Prokka). Only when the
    # user supplied a GFF do we run in coordinate mode and stamp a genome tool.
    has_gff = inp.gff_path is not None
    sample = Sample(
        sample_key=inp.sample_key,
        release_id=release.id,
        annotation_tool=("prokka" if has_gff else "dbcan4-protein"),
        sample_metadata=metadata.get(inp.sample_key, {}),
    )
    session.add(sample)
    session.flush()

    gene_by_key: dict[str, Gene] = {}
    contigs: set[str] = set()

    if has_gff:
        # --- coordinate mode: genes from a user-provided GFF3 ---
        for feat in parse_gff3(inp.gff_path):
            gene = Gene(
                sample_id=sample.id,
                gene_key=feat.gene_key,
                contig=feat.seqid,
                start=feat.start,
                end=feat.end,
                strand=feat.strand,
                feature_type=feat.feature_type,
                product=feat.product,
                attributes=feat.attributes,
            )
            session.add(gene)
            contigs.add(feat.seqid)
            # Last feature with a given key wins for annotation joins (CDS expected).
            gene_by_key[feat.gene_key] = gene
        session.flush()

    # --- protein sequences: build genes here in protein mode; else join by key ---
    if inp.faa_path:
        seqs = parse_protein_fasta(inp.faa_path)
        matched = 0
        for key, seq in seqs.items():
            gene = gene_by_key.get(key)
            if gene is None and not has_gff:
                # Protein mode: each protein IS a gene. No genomic coordinates —
                # use its own residue span (1..L, aa) and mark it coordinate-free.
                gene = Gene(
                    sample_id=sample.id,
                    gene_key=key,
                    contig="protein",          # sentinel: not a genomic contig
                    start=1,
                    end=len(seq),               # residue span (aa), not bp
                    strand=None,
                    feature_type="protein",
                    product=None,
                    attributes={"coordinate_system": "protein_residues"},
                )
                session.add(gene)
                contigs.add("protein")
                gene_by_key[key] = gene
            if gene is not None:
                gene.protein_seq = seq
                matched += 1
        report.sequences_added += matched
        session.flush()
        _provenance(session, release, inp.sample_key, "sequence", sample.id,
                    str(inp.faa_path), per_file.get("faa", ""),
                    "interproscan", versions.get("interproscan"))

    n_genes = session.query(Gene).filter(Gene.sample_id == sample.id).count()
    report.genes_added += n_genes
    sample.n_genes = n_genes
    sample.n_contigs = len(contigs)

    # --- JBrowse 2 track files: only meaningful with real genomic coordinates ---
    if tracks_dir is not None and has_gff:
        try:
            track = prepare_sample_track(
                inp.gff_path, inp.sample_key, tracks_dir, sorted(contigs)
            )
            sample.track_path = track["track_dir"]
        except Exception as exc:  # track prep must not abort a DB ingest
            report_note = f"track prep failed for {inp.sample_key}: {exc}"
            print(f"[warn] {report_note}")

    _provenance(session, release, inp.sample_key, "sample", sample.id,
                str(inp.gff_path) if has_gff else str(inp.faa_path or inp.dbcan_path),
                combined, ("prokka" if has_gff else "dbcan4-protein"),
                versions.get("prokka") if has_gff else None)

    # --- CAZyme annotations from dbCAN overview ---
    if inp.dbcan_path:
        for call in parse_dbcan_overview(inp.dbcan_path):
            gene = gene_by_key.get(call.gene_key)
            if gene is None:
                report.unmatched_cazyme += 1
                continue
            ann = CazymeAnnotation(
                gene_id=gene.id,
                cazy_family=call.cazy_family,
                ec_number=call.ec_number,
                tool=call.tool,
                n_tools_support=call.n_tools_support,
                substrate=None,  # per-CGC substrate join deferred (see DESIGN.md)
                # Baseline provenance: stamp the release + method metadata so the
                # advanced ingester (a separate release) can be compared against it.
                # method_kind comes from the registry (one source of truth).
                method_family="baseline",
                method_kind=method_kind_of(call.tool),
                release_id=release.id,
                raw={"tool": call.tool},
            )
            session.add(ann)
            report.cazymes_added += 1
        session.flush()
        _provenance(session, release, inp.sample_key, "cazyme", sample.id,
                    str(inp.dbcan_path), per_file.get("dbcan", ""),
                    "run_dbcan", versions.get("run_dbcan"))

    # --- CAZyme gene clusters (CGC) + substrate ---
    if inp.cgc_std_path:
        _load_clusters(session, release, inp, sample, gene_by_key, per_file,
                       versions, report)

    # --- antimicrobial-resistance genes (hAMRonization) ---
    if inp.arg_path:
        _load_args(session, release, inp, sample, gene_by_key, per_file, report)

    # --- InterPro domains (+ normalised GO terms) ---
    if inp.interpro_path:
        go_seen: set[tuple[int, str]] = set()  # de-dup (gene_id, go_id)
        for hit in parse_interproscan(inp.interpro_path):
            gene = gene_by_key.get(hit.gene_key)
            if gene is None:
                report.unmatched_domains += 1
                continue
            dom = InterproDomain(
                gene_id=gene.id,
                analysis=hit.analysis,
                signature_acc=hit.signature_acc,
                signature_desc=hit.signature_desc,
                interpro_acc=hit.interpro_acc,
                interpro_desc=hit.interpro_desc,
                start=hit.start,
                end=hit.end,
                evalue=hit.evalue,
                go_terms=hit.go_terms,
            )
            session.add(dom)
            report.domains_added += 1
            for go_id in hit.go_terms or []:
                key = (gene.id, go_id)
                if key in go_seen:
                    continue
                go_seen.add(key)
                session.add(GoAnnotation(
                    gene_id=gene.id, go_id=go_id, source=hit.analysis))
                report.go_added += 1
        session.flush()
        _provenance(session, release, inp.sample_key, "interpro", sample.id,
                    str(inp.interpro_path), per_file.get("interpro", ""),
                    "interproscan", versions.get("interproscan"))

    report.samples_added += 1


def _load_clusters(session, release, inp, sample, gene_by_key, per_file,
                   versions, report) -> None:
    """Aggregate CGC per-gene rows into CazymeCluster rows, wire member genes via
    Gene.cluster_id, and propagate the per-CGC substrate onto member CAZyme rows."""
    from collections import Counter

    members = list(parse_cgc_standard(inp.cgc_std_path))
    if not members:
        return
    substrates = parse_substrate(inp.substrate_path) if inp.substrate_path else {}

    by_key: dict[str, list] = {}
    for m in members:
        by_key.setdefault(m.cluster_key, []).append(m)

    for key, ms in by_key.items():
        member_genes = [gene_by_key.get(m.gene_key) for m in ms]
        starts = [m.start for m in ms if m.start is not None] + \
                 [g.start for g in member_genes if g is not None]
        ends = [m.end for m in ms if m.end is not None] + \
               [g.end for g in member_genes if g is not None]
        contig = next((m.contig for m in ms if m.contig), None) or next(
            (g.contig for g in member_genes if g is not None), "?")
        types = Counter((m.gene_type or "null") for m in ms)
        composition = ", ".join(f"{t}×{c}" for t, c in types.items())
        sub, score = substrates.get(key, (None, None))

        cluster = CazymeCluster(
            sample_id=sample.id,
            cluster_key=key,
            contig=contig,
            start=min(starts) if starts else 0,
            end=max(ends) if ends else 0,
            composition=composition,
            n_genes=len(ms),
            substrate=sub,
            substrate_score=score,
            raw={"members": [m.gene_key for m in ms]},
        )
        session.add(cluster)
        session.flush()
        for g in member_genes:
            if g is not None:
                g.cluster_id = cluster.id
                if sub:  # propagate per-CGC substrate onto the gene's CAZyme rows
                    session.query(CazymeAnnotation).filter(
                        CazymeAnnotation.gene_id == g.id
                    ).update({CazymeAnnotation.substrate: sub})
        report.clusters_added += 1
    session.flush()
    _provenance(session, release, inp.sample_key, "cluster", sample.id,
                str(inp.cgc_std_path), per_file.get("cgc_std", ""),
                "run_dbcan", versions.get("run_dbcan"))


def _load_args(session, release, inp, sample, gene_by_key, per_file, report) -> None:
    """Load this sample's ARG hits from the run-level hAMRonization report."""
    hits = parse_arg_report(inp.arg_path).get(inp.sample_key, [])
    if not hits:
        return
    for h in hits:
        gene = gene_by_key.get(h.gene_key) if h.gene_key else None
        if gene is None:
            report.unmatched_arg += 1
        session.add(ArgAnnotation(
            gene_id=gene.id if gene is not None else None,
            sample_id=sample.id,
            gene_symbol=h.gene_symbol,
            gene_name=h.gene_name,
            drug_class=h.drug_class,
            resistance_mechanism=h.resistance_mechanism,
            identity=h.identity,
            coverage=h.coverage,
            tool=h.tool,
            reference_db=h.reference_db,
            accession=h.accession,
            raw={"tool": h.tool},
        ))
        report.args_added += 1
    session.flush()
    _provenance(session, release, inp.sample_key, "arg", sample.id,
                str(inp.arg_path), per_file.get("arg", ""),
                "hamronization", None)


def _provenance(session, release, sample_key, entity_type, entity_id,
                source_file, checksum, tool_name, tool_version) -> None:
    session.add(
        Provenance(
            release_id=release.id,
            sample_key=sample_key,
            entity_type=entity_type,
            entity_id=entity_id,
            source_file=source_file,
            file_checksum=checksum,
            tool_name=tool_name,
            tool_version=tool_version,
        )
    )
