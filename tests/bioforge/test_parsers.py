"""Unit tests for the format parsers (header-keyed dbCAN, positional IPS, GFF3)."""
from __future__ import annotations

from bioforge.ingest.parse_dbcan import _families_in_cell, parse_dbcan_overview
from bioforge.ingest.parse_gff import parse_gff3
from bioforge.ingest.parse_interpro import parse_interproscan

GFF = "sample_data/results/annotation/prokka/all/sampleA/sampleA.gff"
DBCAN = "sample_data/results/cazyme/dbcan/cazyme_annotation/sampleA_overview.tsv"
IPS = "sample_data/results/protein_annotation/interproscan/sampleA_interproscan_faa.tsv"


def test_gff_stops_at_fasta_and_parses_attrs():
    feats = list(parse_gff3(GFF))
    assert len(feats) == 5  # FASTA section excluded
    f = feats[0]
    assert f.seqid == "contig1" and f.start == 12 and f.strand == "+"
    assert f.gene_key == "PROKKA_00001"
    assert f.product == "Alpha-amylase"


def test_dbcan_family_normalisation():
    assert _families_in_cell("GH13_31(12-390)") == ["GH13"]
    assert _families_in_cell("CBM48(20-70)+GH9(80-820)") == ["CBM48", "GH9"]
    assert _families_in_cell("-") == []


def test_dbcan_emits_per_tool_family_rows():
    calls = list(parse_dbcan_overview(DBCAN))
    # PROKKA_00005: HMMER(CBM48,GH9)+dbCAN_sub(GH9)+DIAMOND(GH9) = 4 rows
    p5 = [c for c in calls if c.gene_key == "PROKKA_00005"]
    fams = sorted((c.tool, c.cazy_family) for c in p5)
    assert fams == [
        ("DIAMOND", "GH9"),
        ("HMMER", "CBM48"),
        ("HMMER", "GH9"),
        ("dbCAN_sub", "GH9"),
    ]


def test_interpro_positional_parse():
    hits = list(parse_interproscan(IPS))
    assert len(hits) == 5
    h = hits[0]
    assert h.gene_key == "PROKKA_00001"
    assert h.analysis == "Pfam"
    assert h.interpro_acc == "IPR006047"
    assert "GO:0004553" in h.go_terms
    # Row without InterPro acc / GO still parses.
    panther = [x for x in hits if x.analysis == "PANTHER"][0]
    assert panther.interpro_acc is None
