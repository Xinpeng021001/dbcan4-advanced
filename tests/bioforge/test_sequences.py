"""Protein sequence ingest + FASTA endpoints (M2)."""
from __future__ import annotations

from bioforge.api import main as api_main
from bioforge.api import queries as Q
from bioforge.ingest.loader import ingest_directory
from bioforge.ingest.parse_fasta import parse_protein_fasta, to_fasta


def test_parse_fasta_and_serialise(tmp_path):
    faa = tmp_path / "x.faa"
    faa.write_text(">g1 sampleA amylase\nMKAA\nQQ\n>g2\nMMMM\n")
    parsed = parse_protein_fasta(faa)
    assert parsed == {"g1": "MKAAQQ", "g2": "MMMM"}
    fasta = to_fasta([("g1 sampleA amylase", "MKAA")], width=2)
    assert fasta.splitlines() == [">g1 sampleA amylase", "MK", "AA"]


def test_sequences_ingested(Session, sample_data, tmp_path):
    with Session() as s:
        report = ingest_directory(s, sample_data, label="t", tracks_dir=tmp_path / "t")
        assert report.sequences_added == 6          # 4 in sampleA CDS + 2 in sampleB
        gene = next(g for g in Q.sample_genes(s, 1) if g.gene_key == "PROKKA_00001")
        assert gene.protein_seq and gene.protein_seq.startswith("MKFL")
        trna = next(g for g in Q.sample_genes(s, 1) if g.feature_type == "tRNA")
        assert trna.protein_seq is None            # tRNA has no protein


def test_fasta_endpoints(Session, sample_data, monkeypatch, tmp_path):
    with Session() as s:
        ingest_directory(s, sample_data, label="t", tracks_dir=tmp_path / "t")
    monkeypatch.setattr(api_main, "SessionLocal", Session)
    from starlette.testclient import TestClient
    c = TestClient(api_main.create_app())

    r = c.get("/genes/1/protein.faa")
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/x-fasta")
    assert r.text.startswith(">PROKKA_00001")

    r = c.get("/samples/1/proteins.faa")
    assert r.text.count(">") == 4                  # sampleA has 4 proteins
    r = c.get("/browse.faa", params={"family": "GH13"})
    assert r.text.count(">") == 1
