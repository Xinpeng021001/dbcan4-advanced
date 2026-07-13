"""Tests for the advanced-vs-baseline comparison + deep-dive feature helpers
that back the side-by-side UI."""
from __future__ import annotations

import json
from pathlib import Path

from bioforge.api import queries as Q
from bioforge.ingest.loader import ingest_directory
from bioforge.ingest.loader_advanced import ingest_advanced_manifest
from bioforge.models import Gene
from sqlalchemy import select


def _advanced_fixture(tmp: Path) -> Path:
    base = tmp / "cazyme_advanced"
    (base / "predictions/sampleA").mkdir(parents=True)
    (base / "features/sampleA/structures").mkdir(parents=True)

    def tsv(p, hdr, rows):
        with open(p, "w") as fh:
            fh.write("\t".join(hdr) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")

    h = ["protein_id", "family", "confidence", "ec", "all_families", "extra"]
    # PROKKA_00001: baseline also calls GH13 -> agreement (not advanced-only)
    # PROKKA_00005: AA9 advanced-only WITH consensus (kNN + fusion)
    # PROKKA_00002: GT99 called ONLY by centroid (single method) -> must NOT flag
    tsv(base / "predictions/sampleA/ESM-C-kNN.tsv", h, [
        ["PROKKA_00001", "GH13", "0.98", "-", "-", "{}"],
        ["PROKKA_00005", "AA9", "0.95", "-", "-", "{}"]])
    tsv(base / "predictions/sampleA/ESM-C-centroid.tsv", h, [
        ["PROKKA_00002", "GT99", "0.55", "-", "-", "{}"]])  # lone weak call
    tsv(base / "predictions/sampleA/fusion.tsv", h, [
        ["PROKKA_00001", "GH13", "0.97", "-", "-", '{"agreement":2}'],
        ["PROKKA_00005", "AA9", "0.90", "-", "-", '{"agreement":2}']])
    tsv(base / "features/sampleA/signalp6.tsv",
        ["protein_id", "prediction", "sp_prob", "cs_position", "extra"],
        [["PROKKA_00001", "SP", "0.97", "20", "{}"]])
    tsv(base / "features/sampleA/deeptmhmm.tsv",
        ["protein_id", "prediction", "n_tm", "topology", "extra"],
        [["PROKKA_00001", "Globular", "0", "-", "{}"]])
    (base / "features/sampleA/structures/PROKKA_00005.pdb").write_text(
        "ATOM      1  CA  ALA A   1       0.0     0.0     0.0  1.00 90.00           C\nEND\n")
    tsv(base / "features/sampleA/structures.tsv",
        ["protein_id", "source", "plddt", "length", "path", "extra"],
        [["PROKKA_00005", "ESMFold", "88.5", "129", "structures/PROKKA_00005.pdb", "{}"]])

    (base / "manifest.json").write_text(json.dumps({
        "contract_version": "1.0", "pipeline": "dbcan4-advanced",
        "release_label": "adv-ui-test",
        "samples": [{"sample_key": "sampleA",
                     "cazyme_predictions": [
                         {"tool": "ESM-C-kNN", "path": "predictions/sampleA/ESM-C-kNN.tsv"},
                         {"tool": "ESM-C-centroid", "path": "predictions/sampleA/ESM-C-centroid.tsv"},
                         {"tool": "fusion", "path": "predictions/sampleA/fusion.tsv"}],
                     "protein_features": [
                         {"feature_type": "signal_peptide", "tool": "SignalP6",
                          "path": "features/sampleA/signalp6.tsv"},
                         {"feature_type": "tm_topology", "tool": "DeepTMHMM",
                          "path": "features/sampleA/deeptmhmm.tsv"},
                         {"feature_type": "structure", "tool": "ESMFold",
                          "path": "features/sampleA/structures.tsv"}]}]}))
    return base / "manifest.json"


def _load(Session, sample_data, tmp_path):
    with Session() as s:
        ingest_directory(s, sample_data, label="baseline")
    man = _advanced_fixture(tmp_path)
    with Session() as s:
        ingest_advanced_manifest(s, man, structures_dir=tmp_path / "served")


def _gene(session, key):
    g = session.scalars(select(Gene).where(Gene.gene_key == key)).first()
    return Q.get_gene(session, g.id)


def test_comparison_flags_advanced_only_with_consensus(Session, sample_data, tmp_path):
    _load(Session, sample_data, tmp_path)
    with Session() as s:
        cmp = Q.gene_cazyme_comparison(_gene(s, "PROKKA_00005"))
    # AA9 recovered by kNN + fusion, absent from baseline -> advanced-only.
    assert "AA9" in cmp["advanced_only"]
    aa9 = [r for r in cmp["rows"] if r["family"] == "AA9"][0]
    assert aa9["advanced_only"] is True
    assert not aa9["in_baseline"] and aa9["in_advanced"]
    assert aa9["advanced_best_conf"] is not None
    # advanced-only rows sort first.
    assert cmp["rows"][0]["advanced_only"] is True


def test_comparison_agreement_not_flagged(Session, sample_data, tmp_path):
    _load(Session, sample_data, tmp_path)
    with Session() as s:
        cmp = Q.gene_cazyme_comparison(_gene(s, "PROKKA_00001"))
    # GH13 called by BOTH baseline and advanced -> not advanced-only.
    assert "GH13" not in cmp["advanced_only"]
    gh13 = [r for r in cmp["rows"] if r["family"] == "GH13"][0]
    assert gh13["in_baseline"] and gh13["in_advanced"]
    assert gh13["advanced_only"] is False


def test_comparison_lone_weak_call_not_flagged(Session, sample_data, tmp_path):
    """A family called by a SINGLE advanced method (no fusion, no 2nd method)
    must NOT be flagged advanced-only — avoids surfacing centroid noise."""
    _load(Session, sample_data, tmp_path)
    with Session() as s:
        cmp = Q.gene_cazyme_comparison(_gene(s, "PROKKA_00002"))
    # GT99 came only from centroid -> present as advanced call, but NOT flagged.
    assert "GT99" not in cmp["advanced_only"]
    gt99 = [r for r in cmp["rows"] if r["family"] == "GT99"][0]
    assert gt99["in_advanced"] and gt99["advanced_only"] is False


def test_features_by_type(Session, sample_data, tmp_path):
    _load(Session, sample_data, tmp_path)
    with Session() as s:
        feats = Q.gene_features_by_type(_gene(s, "PROKKA_00001"))
        assert len(feats["signal_peptide"]) == 1
        assert feats["signal_peptide"][0].label == "SP"
        assert len(feats["tm_topology"]) == 1
        feats5 = Q.gene_features_by_type(_gene(s, "PROKKA_00005"))
        assert len(feats5["structure"]) == 1
        assert feats5["structure"][0].structure_path.endswith(".pdb")


def test_browse_advanced_only_map(Session, sample_data, tmp_path):
    _load(Session, sample_data, tmp_path)
    with Session() as s:
        genes = list(s.scalars(select(Gene)))
        # reload with cazymes
        genes = [Q.get_gene(s, g.id) for g in genes]
        amap = Q.genes_advanced_only_map(genes)
    # exactly the AA9 gene (PROKKA_00005) has an advanced-only family
    flagged = {gid for gid, fams in amap.items() if fams}
    with Session() as s:
        g5 = s.scalars(select(Gene).where(Gene.gene_key == "PROKKA_00005")).first()
    assert flagged == {g5.id}
