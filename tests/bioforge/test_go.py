"""GO-term normalisation + browse-by-GO / browse-by-EC (M5, M1)."""
from __future__ import annotations

from bioforge.api import queries as Q
from bioforge.ingest.loader import ingest_directory


def test_go_normalised_and_browsable(Session, sample_data, tmp_path):
    with Session() as s:
        report = ingest_directory(s, sample_data, label="t", tracks_dir=tmp_path / "t")
        assert report.go_added == 7          # de-duped (gene, go) pairs across IPS rows
        facets = dict(Q.go_facets(s))
        assert facets["GO:0004553"] == 4     # shared by 4 genes across both samples
        assert facets["GO:0005975"] == 1
        hits = Q.search_genes(s, go="GO:0005975")
        assert [g.gene_key for g in hits] == ["PROKKA_00001"]


def test_ec_browse(Session, sample_data, tmp_path):
    with Session() as s:
        ingest_directory(s, sample_data, label="t", tracks_dir=tmp_path / "t")
        facets = dict(Q.ec_facets(s))
        assert "3.2.1.1" in facets
        hits = Q.search_genes(s, ec="3.2.1.1")
        assert [g.gene_key for g in hits] == ["PROKKA_00001"]
