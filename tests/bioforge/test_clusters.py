"""CGC cluster + substrate ingest, queries, and track band rendering (M3)."""
from __future__ import annotations

from bioforge.api import queries as Q
from bioforge.api.track import render_track
from bioforge.ingest.loader import ingest_directory
from bioforge.ingest.parse_cgc import parse_cgc_standard
from bioforge.ingest.parse_substrate import parse_substrate


def test_parsers_handle_zero_index_columns(tmp_path):
    # Regression: column at index 0 must be found (not treated as falsy).
    cgc = tmp_path / "s_cgc_standard_out.tsv"
    cgc.write_text(
        "CGC#\tGene Type\tContig ID\tProtein ID\tGene Start\tGene Stop\tDirection\tProtein Family\n"
        "CGC1\tCAZyme\tc1\tG1\t10\t99\t+\tGH1\n"
    )
    members = list(parse_cgc_standard(cgc))
    assert members[0].cluster_key == "CGC1" and members[0].gene_key == "G1"
    assert members[0].start == 10 and members[0].family == "GH1"

    sub = tmp_path / "s_substrate_prediction.tsv"
    sub.write_text("Cluster ID\tBest hit PUL\tSubstrate\tScore\nCGC1\tPUL1\txylan\t9.0\n")
    assert parse_substrate(sub) == {"CGC1": ("xylan", 9.0)}


def test_cluster_ingest(Session, sample_data, tmp_path):
    with Session() as s:
        report = ingest_directory(s, sample_data, label="t", tracks_dir=tmp_path / "t")
        assert report.clusters_added == 2
        clusters = Q.sample_clusters(s, 1)
        assert len(clusters) == 1
        cgc = clusters[0]
        assert cgc.cluster_key == "CGC1"
        assert cgc.substrate == "starch"
        assert {g.gene_key for g in cgc.genes} == {
            "PROKKA_00001", "PROKKA_00002", "PROKKA_00003"}
        # substrate propagated onto member CAZyme rows
        g1 = next(g for g in Q.sample_genes(s, 1) if g.gene_key == "PROKKA_00001")
        assert all(c.substrate == "starch" for c in g1.cazymes)
        assert g1.cluster_id == cgc.id


def test_cluster_band_in_track(Session, sample_data, tmp_path):
    with Session() as s:
        ingest_directory(s, sample_data, label="t", tracks_dir=tmp_path / "t")
        genes = Q.sample_genes(s, 1)
        clusters = Q.sample_clusters(s, 1)
        svg = render_track(genes, "contig1", 1, 2200, clusters=clusters)
    assert "cluster-band" in svg
    assert "CGC1" in svg
