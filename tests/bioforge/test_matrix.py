"""Comparative matrix + unified search (M6)."""
from __future__ import annotations

from bioforge.api import main as api_main
from bioforge.api import queries as Q
from bioforge.ingest.loader import ingest_directory


def _client(Session, sample_data, monkeypatch, tmp_path):
    with Session() as s:
        ingest_directory(s, sample_data, label="t", tracks_dir=tmp_path / "t")
    monkeypatch.setattr(api_main, "SessionLocal", Session)
    from starlette.testclient import TestClient
    return TestClient(api_main.create_app())


def test_family_matrix(Session, sample_data, tmp_path):
    with Session() as s:
        ingest_directory(s, sample_data, label="t", tracks_dir=tmp_path / "t")
        samples, rows = Q.comparative_matrix(s, by="family")
        keys = [s_.sample_key for s_ in samples]
        assert keys == ["sampleA", "sampleB"]
        row = {label: counts for label, counts, _ in rows}
        assert row["GH18"] == [0, 1]     # GH18 only in sampleB
        assert row["GH13"] == [1, 0]     # GH13 only in sampleA


def test_drug_class_matrix(Session, sample_data, tmp_path):
    with Session() as s:
        ingest_directory(s, sample_data, label="t", tracks_dir=tmp_path / "t")
        _samples, rows = Q.comparative_matrix(s, by="drug_class")
        labels = {label for label, _c, _t in rows}
        assert {"macrolide", "tetracycline", "beta-lactam"} <= labels


def test_unified_search(Session, sample_data, monkeypatch, tmp_path):
    c = _client(Session, sample_data, monkeypatch, tmp_path)
    types = {h["type"] for h in c.get("/api/search", params={"q": "amylase"}).json()}
    assert "gene" in types
    fam = c.get("/api/search", params={"q": "GH13"}).json()
    assert any(h["type"] == "CAZy family" and h["label"] == "GH13" for h in fam)
    assert c.get("/api/matrix", params={"by": "drug_class"}).json()["by"] == "drug_class"
