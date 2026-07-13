"""Read queries backing the API and pages.

"Current" data = the latest sample version per sample_key (highest release_id).
We compute it at query time rather than trusting a flag, so history in older
releases stays queryable while the UI shows current state by default.

Query helpers here fall into three groups: row fetchers (`search_genes`,
`sample_genes`, `get_gene`), aggregate facets/stats that back the dashboard and
filter dropdowns (`dataset_stats`, `*_facets`), and small serialisers
(`gene_annotation_summary`) shared by the JSON API.
"""
from __future__ import annotations

from dataclasses import dataclass

from sqlalchemy import func, select
from sqlalchemy.orm import Session, selectinload

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

# Sort keys the browse page understands, mapped to ORDER BY columns. Keeping this
# an allow-list (rather than trusting a raw column name from the query string)
# means an arbitrary ?sort= value can never reach SQL.
GENE_SORTS = {
    "location": (Sample.sample_key, Gene.contig, Gene.start),
    "name": (Gene.gene_key,),
    "product": (Gene.product,),
    "length": (Gene.length.desc(),),
}


def latest_sample_ids(session: Session) -> list[int]:
    """IDs of the newest Sample row per sample_key."""
    sub = (
        select(Sample.sample_key, func.max(Sample.release_id).label("mr"))
        .group_by(Sample.sample_key)
        .subquery()
    )
    stmt = select(Sample.id).join(
        sub, (Sample.sample_key == sub.c.sample_key) & (Sample.release_id == sub.c.mr)
    )
    return list(session.scalars(stmt))


def list_samples(session: Session) -> list[Sample]:
    ids = latest_sample_ids(session)
    if not ids:
        return []
    stmt = select(Sample).where(Sample.id.in_(ids)).order_by(Sample.sample_key)
    return list(session.scalars(stmt))


def get_sample(session: Session, sample_id: int) -> Sample | None:
    return session.get(Sample, sample_id)


def sample_genes(session: Session, sample_id: int) -> list[Gene]:
    stmt = (
        select(Gene)
        .where(Gene.sample_id == sample_id)
        .options(selectinload(Gene.cazymes), selectinload(Gene.domains))
        .order_by(Gene.contig, Gene.start)
    )
    return list(session.scalars(stmt))


def get_gene(session: Session, gene_id: int) -> Gene | None:
    """A single gene with its sample + all annotations eagerly loaded."""
    stmt = (
        select(Gene)
        .where(Gene.id == gene_id)
        .options(
            selectinload(Gene.sample),
            selectinload(Gene.cazymes),
            selectinload(Gene.domains),
            selectinload(Gene.args),
            selectinload(Gene.cluster),
            selectinload(Gene.features),
            selectinload(Gene.go_annotations),
        )
    )
    return session.scalars(stmt).first()


def gene_cazyme_comparison(gene: Gene) -> dict:
    """Split a gene's CAZyme calls into baseline vs advanced and compute the
    advanced-only families — the data behind the side-by-side view.

    Returns a dict with:
      baseline   : {family: [calls]}          calls from HMMER/dbCAN_sub/DIAMOND
      advanced   : {family: [calls]}          calls from ESM-C/Foldseek/SaProt/fusion
      families   : sorted union of all families (family-level, subfamily collapsed)
      advanced_only : set of families ONLY the advanced methods called
      rows       : per-family rows for the comparison table (family, baseline
                   tools, advanced tools+confidence, advanced_only flag)
    """
    import re
    from ..methods import display_of, colour_of, is_advanced as tool_is_advanced

    def fam_level(f: str) -> str:
        m = re.match(r"^([A-Za-z]+\d+)", str(f))
        return m.group(1) if m else str(f)

    def is_placeholder(f: str) -> bool:
        # GH0 / GT0 / CBM0 = 'unclassified' placeholder, not a real family call.
        return bool(re.match(r"^[A-Za-z]+0$", str(f)))

    baseline: dict[str, list] = {}
    advanced: dict[str, list] = {}
    for c in gene.cazymes:
        bucket = advanced if (c.method_family == "advanced" or tool_is_advanced(c.tool)) else baseline
        bucket.setdefault(fam_level(c.cazy_family), []).append(c)

    base_fams = set(baseline)
    adv_fams = set(advanced)
    all_fams = sorted(base_fams | adv_fams)

    def advanced_only_flag(fam: str, a_calls: list) -> bool:
        """A family is 'advanced-only' only if (a) no baseline tool called it,
        (b) it is not an unclassified placeholder, and (c) it has CONSENSUS
        support — >=2 distinct advanced methods, or a fusion call. This avoids
        flagging a lone low-agreement centroid guess (the weakest method) as a
        discovery the baseline 'missed'."""
        if fam in base_fams or is_placeholder(fam):
            return False
        tools = {c.tool for c in a_calls}
        has_fusion = "fusion" in tools
        n_methods = len(tools)
        return has_fusion or n_methods >= 2

    rows = []
    advanced_only = set()
    for fam in all_fams:
        b_calls = baseline.get(fam, [])
        a_calls = advanced.get(fam, [])
        confs = [c.confidence for c in a_calls if c.confidence is not None]
        adv_only = advanced_only_flag(fam, a_calls)
        if adv_only:
            advanced_only.add(fam)
        rows.append({
            "family": fam,
            "baseline_tools": sorted({c.tool for c in b_calls if c.tool}),
            "baseline_subfamilies": sorted({c.cazy_family for c in b_calls}),
            "advanced_calls": [
                {"tool": c.tool, "display": display_of(c.tool),
                 "colour": colour_of(c.tool), "confidence": c.confidence,
                 "subfamily": c.cazy_family} for c in a_calls
            ],
            "advanced_best_conf": max(confs) if confs else None,
            "in_baseline": fam in base_fams,
            "in_advanced": fam in adv_fams,
            "advanced_only": adv_only,
        })
    # advanced-only families first, then by best confidence
    rows.sort(key=lambda r: (not r["advanced_only"],
                             -(r["advanced_best_conf"] or 0), r["family"]))
    return {
        "baseline": baseline, "advanced": advanced,
        "families": all_fams, "advanced_only": advanced_only,
        "rows": rows,
        "n_baseline": len(base_fams), "n_advanced": len(adv_fams),
        "n_advanced_only": len(advanced_only),
    }


def gene_row_advanced_only(gene: Gene) -> list[str]:
    """The advanced-only families for one gene (consensus-filtered), reusing the
    same logic as the detail view. Cheap: operates on already-loaded cazymes."""
    return sorted(gene_cazyme_comparison(gene)["advanced_only"])


def genes_advanced_only_map(genes: list[Gene]) -> dict[int, list[str]]:
    """{gene_id: [advanced-only families]} for a list of genes (browse rows).
    Requires gene.cazymes to be loaded (search_genes eager-loads them)."""
    return {g.id: gene_row_advanced_only(g) for g in genes}


def gene_features_by_type(gene: Gene) -> dict:
    """Group a gene's ProteinFeature rows by feature_type for the deep-dive view.

    Covers the OUTPUT_CONTRACT v1.1 feature vocabulary: the original
    signal_peptide/tm_topology/structure plus domain (Pfam), structure_hit
    (Foldseek), localization, physicochem, ec_prediction. Domain and
    structure_hit rows are sorted (domains N->C by start; hits by score desc)
    so the template can render them directly."""
    out: dict[str, list] = {"signal_peptide": [], "tm_topology": [], "structure": [],
                            "domain": [], "structure_hit": [], "localization": [],
                            "physicochem": [], "ec_prediction": [], "other": []}
    # Normalize feature_type spelling across ingesters: some emit plural
    # ("domains", "structure_hits") while the contract key is singular. Map
    # known variants onto the canonical key so the template sees them.
    _alias = {"domains": "domain", "structure_hits": "structure_hit",
              "signal_peptides": "signal_peptide", "ec": "ec_prediction"}
    for f in gene.features:
        key = _alias.get(f.feature_type, f.feature_type)
        out.get(key, out["other"]).append(f)
    # stable, meaningful ordering for the multi-row types
    out["domain"].sort(key=lambda f: (f.start if f.start is not None else 0))
    out["structure_hit"].sort(key=lambda f: (f.score if f.score is not None else 0),
                              reverse=True)
    return out


# Palette for the structure viewer's colour-by-domain mode (domains drawn in order).
_DOMAIN_PALETTE = ["#3b82f6", "#10b981", "#ef4444", "#8b5cf6",
                   "#0ea5e9", "#f97316", "#14b8a6", "#e11d48"]


def _coord(*vals):
    """First value coercible to int, else None (coords may be str or nested)."""
    for v in vals:
        if v is None:
            continue
        try:
            return int(float(v))
        except (TypeError, ValueError):
            continue
    return None


def structure_domain_spec(features: dict) -> list[dict]:
    """Build the colour-by-domain spec the gene-page 3D viewer consumes.

    Emits [{name, family, start, end, color}, ...] — signal peptide first (amber),
    then each Pfam/CAZyme domain in N->C order in a distinct palette colour. Coordinates
    and labels may live in top-level columns OR nested under `attributes` depending on the
    ingester, so both are tried; any span without resolvable start/end is skipped. Consumed
    by the template as `structure_domains | tojson`, so JSON escaping is handled downstream.
    """
    spec: list[dict] = []
    sp = (features.get("signal_peptide") or [None])[0]
    if sp is not None:
        a = getattr(sp, "attributes", None) or {}
        span = a.get("sp_span") if isinstance(a, dict) else None
        s = _coord(getattr(sp, "start", None), span[0] if span else None)
        e = _coord(getattr(sp, "end", None), span[1] if span else None)
        if s is not None and e is not None:
            spec.append({"name": "Signal peptide", "family": "SP (SignalP/DeepTMHMM)",
                         "start": s, "end": e, "color": "#f59e0b"})
    for i, d in enumerate(features.get("domain") or []):
        a = getattr(d, "attributes", None) or {}
        a = a if isinstance(a, dict) else {}
        s = _coord(getattr(d, "start", None), a.get("start"))
        e = _coord(getattr(d, "end", None), a.get("end"))
        if s is None or e is None:
            continue
        family = a.get("acc") or a.get("name") or getattr(d, "label", None) or "domain"
        name = a.get("name") or getattr(d, "label", None) or family
        spec.append({"name": str(name), "family": str(family),
                     "start": s, "end": e,
                     "color": _DOMAIN_PALETTE[len(spec) % len(_DOMAIN_PALETTE)]})
    return spec


def _search_stmt(
    ids: list[int],
    q: str | None,
    cazy_family: str | None,
    sample_key: str | None,
    ec: str | None = None,
    drug_class: str | None = None,
    go: str | None = None,
):
    """Shared WHERE-building for search_genes / count_genes (one source of truth).

    New filters are added here once and both search_genes and count_genes (and
    therefore the browse page + CSV/FASTA exports) pick them up automatically.
    """
    stmt = select(Gene).join(Sample, Gene.sample_id == Sample.id).where(
        Gene.sample_id.in_(ids)
    )
    if q:
        like = f"%{q}%"
        stmt = stmt.where(Gene.product.ilike(like) | Gene.gene_key.ilike(like))
    if sample_key:
        stmt = stmt.where(Sample.sample_key == sample_key)
    if cazy_family:
        stmt = stmt.where(
            Gene.id.in_(
                select(CazymeAnnotation.gene_id).where(
                    CazymeAnnotation.cazy_family.ilike(cazy_family)
                )
            )
        )
    if ec:
        stmt = stmt.where(
            Gene.id.in_(
                select(CazymeAnnotation.gene_id).where(
                    CazymeAnnotation.ec_number == ec
                )
            )
        )
    if drug_class:
        stmt = stmt.where(
            Gene.id.in_(
                select(ArgAnnotation.gene_id).where(
                    ArgAnnotation.drug_class.ilike(drug_class),
                    ArgAnnotation.gene_id.is_not(None),
                )
            )
        )
    if go:
        stmt = stmt.where(
            Gene.id.in_(select(GoAnnotation.gene_id).where(GoAnnotation.go_id == go))
        )
    return stmt


def search_genes(
    session: Session,
    q: str | None = None,
    cazy_family: str | None = None,
    sample_key: str | None = None,
    ec: str | None = None,
    drug_class: str | None = None,
    go: str | None = None,
    limit: int = 200,
    offset: int = 0,
    sort: str = "location",
) -> list[Gene]:
    """Search across current samples by gene name/product, CAZy family, sample."""
    ids = latest_sample_ids(session)
    if not ids:
        return []
    stmt = _search_stmt(ids, q, cazy_family, sample_key, ec, drug_class, go).options(
        selectinload(Gene.cazymes),
        selectinload(Gene.domains),
        selectinload(Gene.sample),
    )
    stmt = stmt.order_by(*GENE_SORTS.get(sort, GENE_SORTS["location"]))
    stmt = stmt.limit(limit).offset(offset)
    return list(session.scalars(stmt))


def count_genes(
    session: Session,
    q: str | None = None,
    cazy_family: str | None = None,
    sample_key: str | None = None,
    ec: str | None = None,
    drug_class: str | None = None,
    go: str | None = None,
) -> int:
    """Total matches for a search (for pagination), ignoring limit/offset."""
    ids = latest_sample_ids(session)
    if not ids:
        return 0
    inner = _search_stmt(ids, q, cazy_family, sample_key, ec, drug_class, go).subquery()
    return session.scalar(select(func.count()).select_from(inner)) or 0


def go_facets(session: Session, limit: int = 40) -> list[tuple[str, int]]:
    """Distinct GO terms across current samples, with gene counts (browse-by-GO)."""
    ids = latest_sample_ids(session)
    if not ids:
        return []
    stmt = (
        select(GoAnnotation.go_id, func.count(func.distinct(GoAnnotation.gene_id)))
        .join(Gene, GoAnnotation.gene_id == Gene.id)
        .where(Gene.sample_id.in_(ids))
        .group_by(GoAnnotation.go_id)
        .order_by(func.count(func.distinct(GoAnnotation.gene_id)).desc(), GoAnnotation.go_id)
        .limit(limit)
    )
    return [(go, n) for go, n in session.execute(stmt)]


def ec_facets(session: Session, limit: int = 40) -> list[tuple[str, int]]:
    """Distinct EC numbers across current samples, with gene counts (browse-by-EC)."""
    ids = latest_sample_ids(session)
    if not ids:
        return []
    stmt = (
        select(
            CazymeAnnotation.ec_number,
            func.count(func.distinct(CazymeAnnotation.gene_id)),
        )
        .join(Gene, CazymeAnnotation.gene_id == Gene.id)
        .where(Gene.sample_id.in_(ids), CazymeAnnotation.ec_number.is_not(None))
        .group_by(CazymeAnnotation.ec_number)
        .order_by(func.count(func.distinct(CazymeAnnotation.gene_id)).desc())
        .limit(limit)
    )
    return [(ec, n) for ec, n in session.execute(stmt)]


def arg_drug_class_facets(session: Session) -> list[tuple[str, int]]:
    """Distinct ARG drug classes across current samples, with hit counts."""
    ids = latest_sample_ids(session)
    if not ids:
        return []
    stmt = (
        select(ArgAnnotation.drug_class, func.count())
        .where(ArgAnnotation.sample_id.in_(ids),
               ArgAnnotation.drug_class.is_not(None))
        .group_by(ArgAnnotation.drug_class)
        .order_by(func.count().desc())
    )
    return [(dc, n) for dc, n in session.execute(stmt)]


def sample_args(session: Session, sample_id: int) -> list[ArgAnnotation]:
    stmt = (
        select(ArgAnnotation)
        .where(ArgAnnotation.sample_id == sample_id)
        .options(selectinload(ArgAnnotation.gene))
        .order_by(ArgAnnotation.drug_class, ArgAnnotation.gene_symbol)
    )
    return list(session.scalars(stmt))


def cazyme_family_facets(session: Session) -> list[tuple[str, int]]:
    """Distinct CAZy families across current samples, with gene counts."""
    ids = latest_sample_ids(session)
    if not ids:
        return []
    stmt = (
        select(
            CazymeAnnotation.cazy_family,
            func.count(func.distinct(CazymeAnnotation.gene_id)),
        )
        .join(Gene, CazymeAnnotation.gene_id == Gene.id)
        .where(Gene.sample_id.in_(ids))
        .group_by(CazymeAnnotation.cazy_family)
        .order_by(CazymeAnnotation.cazy_family)
    )
    return [(fam, n) for fam, n in session.execute(stmt)]


def feature_type_facets(session: Session) -> list[tuple[str, int]]:
    """Gene feature-type breakdown (CDS/tRNA/…) across current samples."""
    ids = latest_sample_ids(session)
    if not ids:
        return []
    stmt = (
        select(Gene.feature_type, func.count())
        .where(Gene.sample_id.in_(ids))
        .group_by(Gene.feature_type)
        .order_by(func.count().desc())
    )
    return [(t, n) for t, n in session.execute(stmt)]


def top_cazyme_families(session: Session, limit: int = 10) -> list[tuple[str, int]]:
    """The most common CAZy families by gene count — powers the dashboard chart."""
    fams = cazyme_family_facets(session)
    return sorted(fams, key=lambda kv: kv[1], reverse=True)[:limit]


@dataclass
class DatasetStats:
    n_samples: int
    n_genes: int
    n_cazymes: int
    n_domains: int
    n_families: int
    n_releases: int
    n_clusters: int = 0
    n_args: int = 0


def dataset_stats(session: Session) -> DatasetStats:
    """Headline counts over the *current* sample set (plus total release count)."""
    ids = latest_sample_ids(session)
    if not ids:
        return DatasetStats(0, 0, 0, 0, 0, session.scalar(
            select(func.count()).select_from(Release)) or 0)
    n_genes = session.scalar(
        select(func.count()).select_from(Gene).where(Gene.sample_id.in_(ids))
    ) or 0
    n_cazymes = session.scalar(
        select(func.count()).select_from(CazymeAnnotation)
        .join(Gene, CazymeAnnotation.gene_id == Gene.id)
        .where(Gene.sample_id.in_(ids))
    ) or 0
    n_domains = session.scalar(
        select(func.count()).select_from(InterproDomain)
        .join(Gene, InterproDomain.gene_id == Gene.id)
        .where(Gene.sample_id.in_(ids))
    ) or 0
    n_families = session.scalar(
        select(func.count(func.distinct(CazymeAnnotation.cazy_family)))
        .join(Gene, CazymeAnnotation.gene_id == Gene.id)
        .where(Gene.sample_id.in_(ids))
    ) or 0
    n_releases = session.scalar(select(func.count()).select_from(Release)) or 0
    n_clusters = session.scalar(
        select(func.count()).select_from(CazymeCluster)
        .where(CazymeCluster.sample_id.in_(ids))
    ) or 0
    n_args = session.scalar(
        select(func.count()).select_from(ArgAnnotation)
        .where(ArgAnnotation.sample_id.in_(ids))
    ) or 0
    return DatasetStats(
        n_samples=len(ids),
        n_genes=n_genes,
        n_cazymes=n_cazymes,
        n_domains=n_domains,
        n_families=n_families,
        n_releases=n_releases,
        n_clusters=n_clusters,
        n_args=n_args,
    )


def sample_stats(session: Session, sample_id: int) -> dict:
    """Per-sample rollups for the sample-detail summary: annotated-gene counts,
    feature-type breakdown, and the sample's own top CAZy families."""
    ft = [
        (t, n)
        for t, n in session.execute(
            select(Gene.feature_type, func.count())
            .where(Gene.sample_id == sample_id)
            .group_by(Gene.feature_type)
            .order_by(func.count().desc())
        )
    ]
    n_cazyme_genes = session.scalar(
        select(func.count(func.distinct(CazymeAnnotation.gene_id)))
        .join(Gene, CazymeAnnotation.gene_id == Gene.id)
        .where(Gene.sample_id == sample_id)
    ) or 0
    n_domain_genes = session.scalar(
        select(func.count(func.distinct(InterproDomain.gene_id)))
        .join(Gene, InterproDomain.gene_id == Gene.id)
        .where(Gene.sample_id == sample_id)
    ) or 0
    families = [
        (fam, n)
        for fam, n in session.execute(
            select(
                CazymeAnnotation.cazy_family,
                func.count(func.distinct(CazymeAnnotation.gene_id)),
            )
            .join(Gene, CazymeAnnotation.gene_id == Gene.id)
            .where(Gene.sample_id == sample_id)
            .group_by(CazymeAnnotation.cazy_family)
            .order_by(func.count(func.distinct(CazymeAnnotation.gene_id)).desc())
            .limit(12)
        )
    ]
    n_clusters = session.scalar(
        select(func.count()).select_from(CazymeCluster)
        .where(CazymeCluster.sample_id == sample_id)
    ) or 0
    n_args = session.scalar(
        select(func.count()).select_from(ArgAnnotation)
        .where(ArgAnnotation.sample_id == sample_id)
    ) or 0
    return {
        "feature_types": ft,
        "n_cazyme_genes": n_cazyme_genes,
        "n_domain_genes": n_domain_genes,
        "n_clusters": n_clusters,
        "n_args": n_args,
        "top_families": families,
    }


def sample_clusters(session: Session, sample_id: int) -> list[CazymeCluster]:
    """CAZyme gene clusters for a sample, contig-ordered, with member genes loaded."""
    stmt = (
        select(CazymeCluster)
        .where(CazymeCluster.sample_id == sample_id)
        .options(selectinload(CazymeCluster.genes))
        .order_by(CazymeCluster.contig, CazymeCluster.start)
    )
    return list(session.scalars(stmt))


def get_cluster(session: Session, cluster_id: int) -> CazymeCluster | None:
    stmt = (
        select(CazymeCluster)
        .where(CazymeCluster.id == cluster_id)
        .options(
            selectinload(CazymeCluster.genes).selectinload(Gene.cazymes),
            selectinload(CazymeCluster.sample),
        )
    )
    return session.scalars(stmt).first()


def gene_provenance(session: Session, gene: Gene) -> list[Provenance]:
    """Provenance rows covering a gene's sample (source files, tools, checksums)."""
    stmt = (
        select(Provenance)
        .where(Provenance.entity_id == gene.sample_id)
        .order_by(Provenance.entity_type)
    )
    return list(session.scalars(stmt))


def gene_neighbors(
    session: Session, gene: Gene, flank_bp: int = 8000, max_genes: int = 60
) -> tuple[list[Gene], int, int]:
    """Genes on the same contig within ``flank_bp`` of ``gene`` (the genomic
    neighbourhood), ordered by start. Returns (genes, region_start, region_end)
    where the region is clamped to the returned features' extent."""
    lo = gene.start - flank_bp
    hi = gene.end + flank_bp
    stmt = (
        select(Gene)
        .where(
            Gene.sample_id == gene.sample_id,
            Gene.contig == gene.contig,
            Gene.end >= lo,
            Gene.start <= hi,
        )
        .options(selectinload(Gene.cazymes))
        .order_by(Gene.start)
        .limit(max_genes)
    )
    genes = list(session.scalars(stmt))
    if not genes:
        return [], max(lo, 1), hi
    region_start = max(min(g.start for g in genes), 1)
    region_end = max(g.end for g in genes)
    return genes, region_start, region_end


def genes_by_contig(genes: list[Gene]) -> list[tuple[str, list[Gene], int, int]]:
    """Group a sample's genes by contig for the per-contig genome track.

    Returns [(contig, genes, region_start, region_end), …] contig-sorted. The
    region spans from 1 to the largest feature end on that contig (we don't store
    contig lengths, so the rightmost feature is the best available extent)."""
    by_contig: dict[str, list[Gene]] = {}
    for g in genes:
        by_contig.setdefault(g.contig, []).append(g)
    out = []
    for contig in sorted(by_contig):
        grp = sorted(by_contig[contig], key=lambda g: g.start)
        out.append((contig, grp, 1, max(g.end for g in grp)))
    return out


def comparative_matrix(session: Session, by: str = "family"):
    """Presence/count matrix of a facet across current samples (pangenome/IMG style).

    Returns (samples, rows) where each row is (label, [count per sample], total),
    sorted by total desc. ``by`` is 'family' (CAZy families × samples, gene counts)
    or 'drug_class' (ARG drug classes × samples, hit counts)."""
    ids = latest_sample_ids(session)
    samples = list_samples(session)
    sample_ids = [s.id for s in samples]
    if not sample_ids:
        return [], []

    if by == "drug_class":
        stmt = (
            select(ArgAnnotation.drug_class, ArgAnnotation.sample_id, func.count())
            .where(ArgAnnotation.sample_id.in_(ids),
                   ArgAnnotation.drug_class.is_not(None))
            .group_by(ArgAnnotation.drug_class, ArgAnnotation.sample_id)
        )
    else:
        by = "family"
        stmt = (
            select(CazymeAnnotation.cazy_family, Gene.sample_id,
                   func.count(func.distinct(CazymeAnnotation.gene_id)))
            .join(Gene, CazymeAnnotation.gene_id == Gene.id)
            .where(Gene.sample_id.in_(ids))
            .group_by(CazymeAnnotation.cazy_family, Gene.sample_id)
        )

    agg: dict[str, dict[int, int]] = {}
    for label, sid, n in session.execute(stmt):
        agg.setdefault(label, {})[sid] = n
    rows = []
    for label, per in agg.items():
        counts = [per.get(sid, 0) for sid in sample_ids]
        rows.append((label, counts, sum(counts)))
    rows.sort(key=lambda r: (-r[2], r[0]))
    return samples, rows


def search_all(session: Session, q: str, per_type: int = 10) -> list[dict]:
    """Unified full-text-ish search across entity types (UniProt/EBI-style).

    Portable ILIKE over the indexed text columns; returns typed hits with a URL.
    (For large corpora, back this with SQLite FTS5 / Postgres tsvector.)"""
    if not q or not q.strip():
        return []
    ids = latest_sample_ids(session)
    if not ids:
        return []
    like = f"%{q.strip()}%"
    out: list[dict] = []

    genes = session.scalars(
        select(Gene).join(Sample, Gene.sample_id == Sample.id)
        .where(Gene.sample_id.in_(ids),
               Gene.product.ilike(like) | Gene.gene_key.ilike(like))
        .options(selectinload(Gene.sample))
        .order_by(Gene.gene_key).limit(per_type)
    ).all()
    for g in genes:
        out.append({"type": "gene", "label": g.gene_key,
                    "sub": f"{g.product or ''} · {g.sample.sample_key}",
                    "url": f"/genes/{g.id}"})

    fams = session.execute(
        select(CazymeAnnotation.cazy_family, func.count(func.distinct(CazymeAnnotation.gene_id)))
        .join(Gene, CazymeAnnotation.gene_id == Gene.id)
        .where(Gene.sample_id.in_(ids), CazymeAnnotation.cazy_family.ilike(like))
        .group_by(CazymeAnnotation.cazy_family).limit(per_type)
    ).all()
    for fam, n in fams:
        out.append({"type": "CAZy family", "label": fam, "sub": f"{n} gene(s)",
                    "url": f"/browse?family={fam}"})

    doms = session.execute(
        select(InterproDomain.interpro_acc, InterproDomain.interpro_desc)
        .join(Gene, InterproDomain.gene_id == Gene.id)
        .where(Gene.sample_id.in_(ids), InterproDomain.interpro_acc.is_not(None),
               InterproDomain.interpro_acc.ilike(like)
               | InterproDomain.interpro_desc.ilike(like))
        .distinct().limit(per_type)
    ).all()
    for acc, desc in doms:
        out.append({"type": "InterPro", "label": acc, "sub": desc or "",
                    "url": f"/browse?q={acc}"})

    args = session.execute(
        select(ArgAnnotation.gene_symbol, ArgAnnotation.drug_class)
        .where(ArgAnnotation.sample_id.in_(ids),
               ArgAnnotation.gene_symbol.ilike(like)
               | ArgAnnotation.drug_class.ilike(like))
        .distinct().limit(per_type)
    ).all()
    for sym, dc in args:
        out.append({"type": "resistance", "label": sym or dc,
                    "sub": dc or "", "url": f"/browse?drug_class={dc}" if dc else "/browse"})

    return out


def list_releases(session: Session) -> list[Release]:
    return list(
        session.scalars(select(Release).order_by(Release.created_at.desc()))
    )


def gene_annotation_summary(gene: Gene) -> dict:
    """Compact dict for API responses."""
    return {
        "id": gene.id,
        "gene_key": gene.gene_key,
        "contig": gene.contig,
        "start": gene.start,
        "end": gene.end,
        "length": gene.length,
        "strand": gene.strand,
        "feature_type": gene.feature_type,
        "product": gene.product,
        "cazy_families": sorted({c.cazy_family for c in gene.cazymes}),
        "interpro": sorted(
            {d.interpro_acc for d in gene.domains if d.interpro_acc}
        ),
    }


# ---------------------------------------------------------------------------
# AI report (on-demand, from the ingested DB) — reuses the pipeline's grounded
# report generator (bioforge.api.ai_report.build_report) via a duck-typed
# assets adapter, so the served report is identical in structure to the one the
# one-command pipeline writes to results/cazyme_advanced/ai_reports/.
# ---------------------------------------------------------------------------
import json as _json
from collections import OrderedDict as _OD
from . import ai_report as _air


def _load_attrs(f):
    """A ProteinFeature.attributes may be a dict or a JSON string; return dict."""
    a = getattr(f, "attributes", None)
    if isinstance(a, str):
        try:
            a = _json.loads(a)
        except (ValueError, TypeError):
            a = {}
    return a if isinstance(a, dict) else {}


class DBAssets:
    """Duck-typed stand-in for build_ai_report.Assets, backed by DB rows for ONE gene.

    The generator only calls .one(name,key,pid), .tsv(name,key), .comprehensive(pid),
    .manifest() and .has(name); we reconstruct exactly the row dicts each evidence
    builder reads, from the gene's ProteinFeature + CazymeAnnotation rows.
    """

    def __init__(self, gene, pid, comprehensive=None, manifest=None):
        self.gene = gene
        self.pid = str(pid)
        self._comprehensive = comprehensive
        self._manifest = manifest or {}
        self._features = gene_features_by_type(gene)
        self._caz = list(getattr(gene, "cazymes", []) or [])
        self._tables = self._build_tables()

    # -- helpers -----------------------------------------------------------
    def _caz_raw(self, tool):
        for c in self._caz:
            if c.tool == tool:
                raw = c.raw
                if isinstance(raw, str):
                    try:
                        raw = _json.loads(raw)
                    except (ValueError, TypeError):
                        raw = {}
                return c, (raw if isinstance(raw, dict) else {})
        return None, {}

    def _build_tables(self):
        pid = self.pid
        t = {}

        # ESM-C heads + fusion from cazyme_annotations
        c, raw = self._caz_raw("ESM-C-kNN")
        if c:
            t["raw_knn.tsv"] = [{"query_id": pid, "knn_pred": c.cazy_family,
                                 "knn_conf": c.confidence,
                                 "knn_purity": raw.get("knn_purity"),
                                 "knn_margin": raw.get("knn_margin")}]
        c, raw = self._caz_raw("ESM-C-centroid")
        if c:
            t["raw_centroid.tsv"] = [{"query_id": pid, "cent_pred": c.cazy_family,
                                      "cent_conf": c.confidence,
                                      "cent_margin": raw.get("cent_margin")}]
        c, raw = self._caz_raw("ESM-C-contrastive")
        if c:
            t["raw_contrastive.tsv"] = [{"query_id": pid, "clf_pred": c.cazy_family,
                                         "clf_conf": c.confidence,
                                         "contr_knn_pred": raw.get("contr_knn_pred"),
                                         "contr_knn_purity": raw.get("contr_knn_purity")}]
        c, raw = self._caz_raw("fusion")
        if c:
            fams = raw.get("all_families")
            votes = dict(raw.get("votes") or {})
            # The pipeline's fusion_raw.tsv counts agreement over the four voters
            # INCLUDING the fusion call's own vote (e.g. 3/4), and the review-flag
            # thresholds are tuned for that convention. The ingester stores votes
            # over the three ESM-C heads only and an agreement count to match. To
            # keep the served report consistent with the pipeline report, add the
            # fusion self-vote and count agreement as voters == final_family / 4.
            votes.setdefault("fusion", c.cazy_family)
            agreement = sum(1 for v in votes.values() if v == c.cazy_family)
            t["fusion_raw.tsv"] = [{"protein_id": pid, "family": c.cazy_family,
                                    "confidence": c.confidence,
                                    "agreement": agreement,
                                    "all_families": (",".join(fams) if isinstance(fams, list) else fams),
                                    "votes": _json.dumps(votes),
                                    "signals": _json.dumps(raw.get("signals") or [])}]

        # Pfam domains
        drows = []
        for f in self._features.get("domain", []):
            a = _load_attrs(f)
            extra = a.get("extra")
            if isinstance(extra, dict):
                extra = _json.dumps(extra)
            drows.append({"protein_id": pid,
                          "acc": a.get("acc"), "name": a.get("name") or getattr(f, "label", None),
                          "start": a.get("start") if a.get("start") is not None else getattr(f, "start", None),
                          "end": a.get("end") if a.get("end") is not None else getattr(f, "end", None),
                          "evalue": a.get("evalue"), "score": a.get("score"),
                          "extra": extra})
        if drows:
            t["domains.tsv"] = drows

        # EC
        ecs = self._features.get("ec_prediction", [])
        if ecs:
            f = ecs[0]; a = _load_attrs(f)
            t["ec_prediction.tsv"] = [{"protein_id": pid,
                                       "ec_number": getattr(f, "label", None),
                                       "confidence": getattr(f, "score", None),
                                       "rank": a.get("rank"), "tool": getattr(f, "tool", None),
                                       "extra": _json.dumps({"confidence_type": a.get("confidence_type")})}]

        # Structure
        sts = self._features.get("structure", [])
        if sts:
            f = sts[0]; a = _load_attrs(f)
            t["structures.tsv"] = [{"protein_id": pid, "plddt": getattr(f, "score", None),
                                    "length": a.get("length"), "source": getattr(f, "label", None) or a.get("source"),
                                    "extra": _json.dumps({"model": a.get("model")})}]

        # DeepTMHMM topology
        tms = self._features.get("tm_topology", [])
        if tms:
            f = tms[0]; a = _load_attrs(f)
            t["deeptmhmm.tsv"] = [{"protein_id": pid, "prediction": getattr(f, "label", None),
                                   "n_tm": a.get("n_tm") or 0, "topology": a.get("topology"),
                                   "extra": _json.dumps({"tool": a.get("tool", "DeepTMHMM"),
                                                         "has_signal_peptide": a.get("has_signal_peptide")})}]

        # Localization
        locs = self._features.get("localization", [])
        if locs:
            f = locs[0]; a = _load_attrs(f)
            t["localization.tsv"] = [{"protein_id": pid, "localization": getattr(f, "label", None),
                                      "confidence": getattr(f, "score", None), "method": a.get("method"),
                                      "extra": _json.dumps({"basis": a.get("basis"),
                                                           "sp_prediction": a.get("sp_prediction")})}]

        # Physicochem
        pcs = self._features.get("physicochem", [])
        if pcs:
            f = pcs[0]; a = _load_attrs(f)
            t["physicochem.tsv"] = [{"protein_id": pid, "mw_da": getattr(f, "score", None),
                                     "pi": a.get("pi"), "instability": a.get("instability"),
                                     "gravy": a.get("gravy"), "aromaticity": a.get("aromaticity"),
                                     "extra": _json.dumps({"n_glyc_sequons": a.get("n_glyc_sequons")})}]
        return t

    # -- Assets interface --------------------------------------------------
    def has(self, name):
        return name in self._tables or (name == "manifest.json" and bool(self._manifest))

    def tsv(self, name, key):
        rows = self._tables.get(name, [])
        out = {}
        for r in rows:
            out.setdefault(str(r.get(key)), []).append(r)
        return out

    def one(self, name, key, pid):
        rows = self.tsv(name, key).get(str(pid))
        return rows[0] if rows else None

    def comprehensive(self, pid):
        return self._comprehensive

    def manifest(self):
        return self._manifest


def build_ai_report_for_gene(gene, comprehensive=None, manifest=None):
    """Assemble the grounded AI report dict for a gene, straight from the DB."""
    pid = gene.gene_key
    A = DBAssets(gene, pid, comprehensive=comprehensive, manifest=manifest)
    return _air.build_report(A, str(pid))
