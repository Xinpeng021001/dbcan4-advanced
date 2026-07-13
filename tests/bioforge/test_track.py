"""Tests for the native genomic-context track (queries + SVG renderer)."""
from __future__ import annotations

from bioforge.api import queries as Q
from bioforge.api.track import render_track
from bioforge.ingest.loader import ingest_directory


def _ingest(Session, sample_data, tmp_path):
    with Session() as s:
        ingest_directory(s, sample_data, label="t", tracks_dir=tmp_path / "tracks")


def test_neighbors_share_contig(Session, sample_data, tmp_path):
    _ingest(Session, sample_data, tmp_path)
    with Session() as s:
        # PROKKA_00001..00004 are on contig1; PROKKA_00005 on contig2.
        gene = next(g for g in Q.sample_genes(s, 1) if g.gene_key == "PROKKA_00001")
        neighbors, r_start, r_end = Q.gene_neighbors(s, gene)
        contigs = {g.contig for g in neighbors}
        assert contigs == {"contig1"}                 # never crosses a contig
        assert gene.id in {g.id for g in neighbors}    # focal gene is included
        assert r_start >= 1 and r_end >= r_start


def test_render_track_svg(Session, sample_data, tmp_path):
    _ingest(Session, sample_data, tmp_path)
    with Session() as s:
        gene = next(g for g in Q.sample_genes(s, 1) if g.gene_key == "PROKKA_00001")
        neighbors, rs, re = Q.gene_neighbors(s, gene)
        svg = render_track(neighbors, gene.contig, rs, re, focal_id=gene.id)
    assert svg.startswith("<svg") and svg.endswith("</svg>")
    assert svg.count("<polygon") == len(neighbors)     # one arrow per feature
    assert "focal" in svg                              # focal gene styled
    assert 'href="/genes/' in svg                      # arrows link to gene pages


def test_genes_by_contig_grouping(Session, sample_data, tmp_path):
    _ingest(Session, sample_data, tmp_path)
    with Session() as s:
        groups = Q.genes_by_contig(Q.sample_genes(s, 1))
    contigs = [c for c, _, _, _ in groups]
    assert contigs == sorted(contigs)                  # contig-sorted
    for _contig, grp, rstart, rend in groups:
        assert rstart == 1
        assert rend == max(g.end for g in grp)
