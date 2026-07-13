"""ARG (antimicrobial resistance) ingest, filter, and API (M4)."""
from __future__ import annotations

from bioforge.api import main as api_main
from bioforge.api import queries as Q
from bioforge.ingest.loader import ingest_directory
from bioforge.ingest.parse_arg import parse_arg_report


def test_parse_arg_report(tmp_path):
    p = tmp_path / "hamronization_combined_report.tsv"
    p.write_text(
        "input_file_name\tinput_protein_id\tgene_symbol\tdrug_class\t"
        "sequence_identity\tanalysis_software_name\treference_database_name\n"
        "sampleX.tsv\tG7\tblaTEM\tbeta-lactam\t99.9%\trgi\tCARD\n"
    )
    hits = parse_arg_report(p)
    assert set(hits) == {"sampleX"}          # filename stem used as sample key
    h = hits["sampleX"][0]
    assert h.gene_key == "G7" and h.gene_symbol == "blaTEM"
    assert h.identity == 99.9 and h.reference_db == "CARD"


def test_arg_ingest_and_facets(Session, sample_data, tmp_path):
    with Session() as s:
        report = ingest_directory(s, sample_data, label="t", tracks_dir=tmp_path / "t")
        assert report.args_added == 3 and report.unmatched_arg == 0
        facets = dict(Q.arg_drug_class_facets(s))
        assert facets.get("macrolide") == 1 and facets.get("beta-lactam") == 1
        # ARG joins to its gene
        args = Q.sample_args(s, 1)
        assert any(a.gene is not None and a.gene.gene_key == "PROKKA_00003" for a in args)


def test_browse_drug_class_filter(Session, sample_data, monkeypatch, tmp_path):
    with Session() as s:
        ingest_directory(s, sample_data, label="t", tracks_dir=tmp_path / "t")
    monkeypatch.setattr(api_main, "SessionLocal", Session)
    from starlette.testclient import TestClient
    c = TestClient(api_main.create_app())

    r = c.get("/api/genes/search", params={"drug_class": "macrolide"})
    assert [g["gene_key"] for g in r.json()] == ["PROKKA_00003"]
    r = c.get("/api/args/search", params={"drug_class": "beta-lactam"})
    assert len(r.json()) == 1 and r.json()[0]["gene_symbol"] == "blaOXA-48"
