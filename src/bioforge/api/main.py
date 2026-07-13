"""FastAPI app: JSON REST API + server-rendered Jinja2 pages + static tracks.

Routes
  Pages:  /                dashboard (dataset stats + charts + samples)
          /browse          gene browse + search (paginated, sortable, CSV export)
          /samples/{id}    sample detail (embeds JBrowse, annotation table)
          /genes/{id}      single-gene detail (annotations + provenance)
          /releases        release / changelog list
          /jbrowse/{id}    standalone JBrowse 2 embed (iframe target)
  API:    /api/stats, /api/samples, /api/genes/search, /api/samples/{id},
          /api/releases
  Export: /browse.csv      current gene search as CSV
  Static: /tracks/...      per-sample bgzip/tabix GFF + config.json (for JBrowse)
"""
from __future__ import annotations

import csv
import io
import math
from pathlib import Path

from fastapi import Depends, FastAPI, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from sqlalchemy.orm import Session

from ..config import tracks_dir
from ..db import make_session_factory
from ..ingest.parse_fasta import to_fasta
from . import queries as Q
from . import xrefs
from .track import render_track

BASE = Path(__file__).resolve().parent
templates = Jinja2Templates(directory=str(BASE / "templates"))
# Expose the cross-reference resolver to every template/macro so accessions can
# render an outbound link without the template importing Python.
templates.env.globals["xref"] = xrefs.xref

PER_PAGE = 50  # gene rows per browse page

SessionLocal = make_session_factory()


def get_session() -> Session:
    with SessionLocal() as s:
        yield s


def create_app() -> FastAPI:
    app = FastAPI(title="BioForge", version="0.2.0")

    tdir = tracks_dir()
    app.mount("/tracks", StaticFiles(directory=str(tdir)), name="tracks")
    static_dir = BASE / "static"
    static_dir.mkdir(exist_ok=True)
    app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")
    # Served protein structures (PDBs copied here by the advanced ingester) — the
    # 3D viewer on the gene deep-dive fetches from /structures/<sample>/<pid>.pdb.
    struct_dir = tdir.parent / "structures"
    struct_dir.mkdir(parents=True, exist_ok=True)
    app.mount("/structures", StaticFiles(directory=str(struct_dir)), name="structures")

    # ---------------- Pages ----------------
    @app.get("/", response_class=HTMLResponse)
    def dashboard(request: Request, session: Session = Depends(get_session)):
        return templates.TemplateResponse(
            request,
            "dashboard.html",
            {
                "stats": Q.dataset_stats(session),
                "samples": Q.list_samples(session),
                "top_families": Q.top_cazyme_families(session, limit=12),
                "feature_types": Q.feature_type_facets(session),
                "releases": Q.list_releases(session)[:5],
            },
        )

    @app.get("/browse", response_class=HTMLResponse)
    def browse(
        request: Request,
        q: str | None = None,
        family: str | None = None,
        sample: str | None = None,
        ec: str | None = None,
        drug_class: str | None = None,
        go: str | None = None,
        advanced_only: int = 0,
        sort: str = "location",
        page: int = 1,
        session: Session = Depends(get_session),
    ):
        searched = bool(q or family or sample or ec or drug_class or go or advanced_only)
        page = max(1, page)
        if advanced_only:
            # Advanced-only filter needs per-gene consensus computation, so fetch
            # the full current matching set, filter, then paginate in Python.
            all_matches = Q.search_genes(
                session, q=q, cazy_family=family, sample_key=sample, ec=ec,
                drug_class=drug_class, go=go, limit=100000, offset=0, sort=sort)
            adv_map = Q.genes_advanced_only_map(all_matches)
            matches = [g for g in all_matches if adv_map.get(g.id)]
            total = len(matches)
            pages = max(1, math.ceil(total / PER_PAGE))
            page = min(page, pages)
            genes = matches[(page - 1) * PER_PAGE: page * PER_PAGE]
            advonly_map = {g.id: adv_map[g.id] for g in genes}
        else:
            total = Q.count_genes(session, q=q, cazy_family=family, sample_key=sample,
                                  ec=ec, drug_class=drug_class, go=go)
            pages = max(1, math.ceil(total / PER_PAGE))
            page = min(page, pages)
            genes = Q.search_genes(
                session, q=q, cazy_family=family, sample_key=sample, ec=ec,
                drug_class=drug_class, go=go,
                limit=PER_PAGE, offset=(page - 1) * PER_PAGE, sort=sort,
            )
            advonly_map = Q.genes_advanced_only_map(genes)
        return templates.TemplateResponse(
            request,
            "browse.html",
            {
                "genes": genes,
                "advonly_map": advonly_map,
                "advanced_only": 1 if advanced_only else 0,
                "q": q or "", "family": family or "", "sample": sample or "",
                "ec": ec or "", "drug_class": drug_class or "", "go": go or "",
                "sort": sort,
                "families": Q.cazyme_family_facets(session),
                "ecs": Q.ec_facets(session),
                "drug_classes": Q.arg_drug_class_facets(session),
                "gos": Q.go_facets(session),
                "samples": Q.list_samples(session),
                "searched": searched,
                "total": total, "page": page, "pages": pages,
                "per_page": PER_PAGE,
            },
        )

    @app.get("/browse.csv")
    def browse_csv(
        q: str | None = None,
        family: str | None = None,
        sample: str | None = None,
        ec: str | None = None,
        drug_class: str | None = None,
        go: str | None = None,
        sort: str = "location",
        session: Session = Depends(get_session),
    ):
        genes = Q.search_genes(
            session, q=q, cazy_family=family, sample_key=sample, ec=ec,
            drug_class=drug_class, go=go, limit=100_000, sort=sort,
        )
        buf = io.StringIO()
        w = csv.writer(buf)
        w.writerow([
            "gene_key", "sample", "contig", "start", "end", "length",
            "strand", "feature_type", "product", "cazy_families", "interpro",
        ])
        for g in genes:
            w.writerow([
                g.gene_key, g.sample.sample_key, g.contig, g.start, g.end,
                g.length, g.strand or "", g.feature_type, g.product or "",
                ";".join(sorted({c.cazy_family for c in g.cazymes})),
                ";".join(sorted({d.interpro_acc for d in g.domains if d.interpro_acc})),
            ])
        buf.seek(0)
        return StreamingResponse(
            iter([buf.getvalue()]),
            media_type="text/csv",
            headers={"Content-Disposition": "attachment; filename=bioforge_genes.csv"},
        )

    def _fasta_response(genes, filename: str) -> StreamingResponse:
        records = [
            (f"{g.gene_key} {g.sample.sample_key} {g.product or ''}".strip(),
             g.protein_seq)
            for g in genes if g.protein_seq
        ]
        return StreamingResponse(
            iter([to_fasta(records)]),
            media_type="text/x-fasta",
            headers={"Content-Disposition": f"attachment; filename={filename}"},
        )

    @app.get("/browse.faa")
    def browse_faa(
        q: str | None = None,
        family: str | None = None,
        sample: str | None = None,
        ec: str | None = None,
        drug_class: str | None = None,
        go: str | None = None,
        sort: str = "location",
        session: Session = Depends(get_session),
    ):
        genes = Q.search_genes(
            session, q=q, cazy_family=family, sample_key=sample, ec=ec,
            drug_class=drug_class, go=go, limit=100_000, sort=sort,
        )
        return _fasta_response(genes, "bioforge_proteins.faa")

    @app.get("/genes/{gene_id}/protein.faa")
    def gene_faa(gene_id: int, session: Session = Depends(get_session)):
        gene = Q.get_gene(session, gene_id)
        if gene is None or not gene.protein_seq:
            return HTMLResponse("No protein sequence for this gene.", status_code=404)
        return _fasta_response([gene], f"{gene.gene_key}.faa")

    @app.get("/genes/{gene_id}/ai_report.json")
    def gene_ai_report(gene_id: int, session: Session = Depends(get_session)):
        """Downloadable, grounded 'prompt pack' JSON for one protein.

        Built on demand from the gene's ingested evidence, reusing the pipeline's
        report generator so the served report matches the one dbcan4_workup.sh writes.
        Paste/upload into an LLM to describe the protein and answer questions about it.
        """
        gene = Q.get_gene(session, gene_id)
        if gene is None:
            return JSONResponse({"error": "gene not found"}, status_code=404)
        try:
            report = Q.build_ai_report_for_gene(gene)
        except Exception as exc:  # never 500 the download; surface a diagnostic
            return JSONResponse({"error": "could not build AI report",
                                 "detail": str(exc)}, status_code=500)
        fname = f"{gene.gene_key}_ai_report.json"
        return JSONResponse(
            report,
            headers={"Content-Disposition": f'attachment; filename="{fname}"'},
        )

    @app.get("/samples/{sample_id}/proteins.faa")
    def sample_faa(sample_id: int, session: Session = Depends(get_session)):
        sample = Q.get_sample(session, sample_id)
        if sample is None:
            return HTMLResponse("Sample not found.", status_code=404)
        genes = Q.sample_genes(session, sample_id)
        return _fasta_response(genes, f"{sample.sample_key}_proteins.faa")

    @app.get("/samples/{sample_id}", response_class=HTMLResponse)
    def sample_detail(
        request: Request, sample_id: int, session: Session = Depends(get_session)
    ):
        sample = Q.get_sample(session, sample_id)
        if sample is None:
            return _not_found(request, "Sample not found")
        genes = Q.sample_genes(session, sample_id)
        clusters = Q.sample_clusters(session, sample_id)
        by_contig_cl: dict[str, list] = {}
        for cl in clusters:
            by_contig_cl.setdefault(cl.contig, []).append(cl)
        contig_tracks = [
            (contig, render_track(grp, contig, rs, re,
                                  clusters=by_contig_cl.get(contig)))
            for contig, grp, rs, re in Q.genes_by_contig(genes)
        ]
        return templates.TemplateResponse(
            request,
            "sample_detail.html",
            {"sample": sample, "genes": genes,
             "stats": Q.sample_stats(session, sample_id),
             "clusters": clusters,
             "args": Q.sample_args(session, sample_id),
             "contig_tracks": contig_tracks},
        )

    @app.get("/genes/{gene_id}", response_class=HTMLResponse)
    def gene_detail(
        request: Request, gene_id: int, session: Session = Depends(get_session)
    ):
        gene = Q.get_gene(session, gene_id)
        if gene is None:
            return _not_found(request, "Gene not found")
        # Protein-input mode (dbcan4, no genome): no genomic neighbourhood.
        if gene.contig == "protein":
            neighbors, context_track = [], None
        else:
            neighbors, r_start, r_end = Q.gene_neighbors(session, gene)
            context_track = render_track(
                neighbors, gene.contig, r_start, r_end, focal_id=gene.id
            ) if neighbors else None
        features = Q.gene_features_by_type(gene)
        return templates.TemplateResponse(
            request,
            "gene_detail.html",
            {"gene": gene, "provenance": Q.gene_provenance(session, gene),
             "context_track": context_track,
             "comparison": Q.gene_cazyme_comparison(gene),
             "features": features,
             "structure_domains": Q.structure_domain_spec(features),
             "n_neighbors": len(neighbors) - 1 if neighbors else 0},
        )

    @app.get("/compare", response_class=HTMLResponse)
    def compare(
        request: Request, by: str = "family",
        session: Session = Depends(get_session),
    ):
        by = "drug_class" if by == "drug_class" else "family"
        samples, rows = Q.comparative_matrix(session, by=by)
        return templates.TemplateResponse(
            request, "compare.html",
            {"samples": samples, "rows": rows, "by": by},
        )

    @app.get("/search", response_class=HTMLResponse)
    def search_page(
        request: Request, q: str | None = None,
        session: Session = Depends(get_session),
    ):
        hits = Q.search_all(session, q) if q else []
        return templates.TemplateResponse(
            request, "search.html", {"q": q or "", "hits": hits},
        )

    @app.get("/releases", response_class=HTMLResponse)
    def releases_page(request: Request, session: Session = Depends(get_session)):
        return templates.TemplateResponse(
            request, "releases.html", {"releases": Q.list_releases(session)}
        )

    @app.get("/jbrowse/{sample_id}", response_class=HTMLResponse)
    def jbrowse_embed(
        request: Request, sample_id: int, session: Session = Depends(get_session)
    ):
        sample = Q.get_sample(session, sample_id)
        if sample is None or not sample.track_path:
            return HTMLResponse(
                "<p style='font-family:sans-serif;padding:1rem'>"
                "No genome track available for this sample.</p>"
            )
        return templates.TemplateResponse(
            request,
            "jbrowse_embed.html",
            {"sample": sample,
             "config_url": f"/tracks/{sample.track_path}/config.json"},
        )

    def _not_found(request: Request, message: str) -> HTMLResponse:
        return templates.TemplateResponse(
            request, "not_found.html", {"message": message}, status_code=404
        )

    # ---------------- JSON API ----------------
    @app.get("/api/stats")
    def api_stats(session: Session = Depends(get_session)):
        s = Q.dataset_stats(session)
        return {
            **s.__dict__,
            "top_families": [
                {"family": f, "genes": n}
                for f, n in Q.top_cazyme_families(session, limit=12)
            ],
            "feature_types": [
                {"type": t, "count": n} for t, n in Q.feature_type_facets(session)
            ],
        }

    @app.get("/api/samples")
    def api_samples(session: Session = Depends(get_session)):
        return [
            {
                "id": s.id,
                "sample_key": s.sample_key,
                "annotation_tool": s.annotation_tool,
                "n_genes": s.n_genes,
                "n_contigs": s.n_contigs,
                "release_id": s.release_id,
                "metadata": s.sample_metadata,
            }
            for s in Q.list_samples(session)
        ]

    @app.get("/api/samples/{sample_id}")
    def api_sample(sample_id: int, session: Session = Depends(get_session)):
        s = Q.get_sample(session, sample_id)
        if s is None:
            return JSONResponse({"error": "not found"}, status_code=404)
        genes = Q.sample_genes(session, sample_id)
        return {
            "id": s.id,
            "sample_key": s.sample_key,
            "annotation_tool": s.annotation_tool,
            "n_genes": s.n_genes,
            "n_contigs": s.n_contigs,
            "release_id": s.release_id,
            "track_path": s.track_path,
            "genes": [Q.gene_annotation_summary(g) for g in genes],
        }

    @app.get("/api/genes/search")
    def api_search(
        q: str | None = None,
        family: str | None = Query(None),
        sample: str | None = None,
        ec: str | None = None,
        drug_class: str | None = None,
        go: str | None = None,
        limit: int = 200,
        offset: int = 0,
        sort: str = "location",
        session: Session = Depends(get_session),
    ):
        genes = Q.search_genes(
            session, q=q, cazy_family=family, sample_key=sample, ec=ec,
            drug_class=drug_class, go=go, limit=limit, offset=offset, sort=sort,
        )
        return [Q.gene_annotation_summary(g) for g in genes]

    @app.get("/api/args/search")
    def api_args(
        drug_class: str | None = None,
        session: Session = Depends(get_session),
    ):
        genes = Q.search_genes(session, drug_class=drug_class, limit=100_000)
        out = []
        for g in genes:
            for a in g.args:
                if drug_class and (a.drug_class or "").lower() != drug_class.lower():
                    continue
                out.append({
                    "sample": g.sample.sample_key, "gene_key": g.gene_key,
                    "gene_symbol": a.gene_symbol, "drug_class": a.drug_class,
                    "resistance_mechanism": a.resistance_mechanism,
                    "identity": a.identity, "coverage": a.coverage,
                    "tool": a.tool, "reference_db": a.reference_db,
                    "accession": a.accession,
                })
        return out

    @app.get("/api/matrix")
    def api_matrix(by: str = "family", session: Session = Depends(get_session)):
        samples, rows = Q.comparative_matrix(session, by=by)
        return {
            "by": "drug_class" if by == "drug_class" else "family",
            "samples": [s.sample_key for s in samples],
            "rows": [
                {"label": label, "counts": counts, "total": total}
                for label, counts, total in rows
            ],
        }

    @app.get("/api/search")
    def api_search_all(q: str = "", session: Session = Depends(get_session)):
        return Q.search_all(session, q)

    @app.get("/api/releases")
    def api_releases(session: Session = Depends(get_session)):
        return [
            {
                "id": r.id,
                "label": r.label,
                "funcscan_version": r.funcscan_version,
                "pipeline_version": r.pipeline_version,
                "is_current": r.is_current,
                "created_at": r.created_at.isoformat() if r.created_at else None,
                "notes": r.notes,
            }
            for r in Q.list_releases(session)
        ]

    return app


app = create_app()
