"""API smoke tests over an ingested fixture DB (in-memory-ish per tmp_path)."""
from __future__ import annotations

import pytest
from starlette.testclient import TestClient

from bioforge.api import main as api_main
from bioforge.ingest.loader import ingest_directory


@pytest.fixture()
def client(Session, sample_data, monkeypatch, tmp_path):
    with Session() as s:
        ingest_directory(s, sample_data, label="t1",
                         tracks_dir=tmp_path / "tracks")
    # Point the app's session factory at the test DB.
    monkeypatch.setattr(api_main, "SessionLocal", Session)
    app = api_main.create_app()
    return TestClient(app)


def test_list_samples(client):
    r = client.get("/api/samples")
    assert r.status_code == 200
    keys = {s["sample_key"] for s in r.json()}
    assert keys == {"sampleA", "sampleB"}


def test_search_by_family(client):
    r = client.get("/api/genes/search", params={"family": "GH13"})
    assert r.status_code == 200
    data = r.json()
    assert len(data) == 1
    assert data[0]["gene_key"] == "PROKKA_00001"
    assert "GH13" in data[0]["cazy_families"]


def test_search_by_product(client):
    r = client.get("/api/genes/search", params={"q": "cellulase"})
    assert [g["gene_key"] for g in r.json()] == ["PROKKA_00005"]


def test_sample_detail_endpoint(client):
    r = client.get("/api/samples/1")
    assert r.status_code == 200
    assert r.json()["n_genes"] == 5


def test_pages_render(client):
    for path in ["/", "/samples/1", "/releases", "/jbrowse/1"]:
        assert client.get(path).status_code == 200
