"""Cross-reference URL resolver (api/xrefs.py)."""
from __future__ import annotations

from bioforge.api import xrefs


def test_cazy_family():
    assert xrefs.cazy("GH13") == "http://www.cazy.org/GH13.html"
    assert xrefs.cazy("CBM48") == "http://www.cazy.org/CBM48.html"
    assert xrefs.cazy("notafamily") is None
    assert xrefs.cazy(None) is None


def test_interpro_and_signature():
    assert xrefs.interpro("IPR006047").endswith("/entry/InterPro/IPR006047/")
    assert xrefs.interpro("PF00128") is None  # not an InterPro accession
    assert "pfam/PF00128" in xrefs.signature("Pfam", "PF00128")
    assert "panther/PTHR10357" in xrefs.signature("PANTHER", "PTHR10357")
    assert xrefs.signature("MysteryDB", "X1") is None


def test_ec_go_arg():
    assert xrefs.ec("3.2.1.1") == "https://enzyme.expasy.org/EC/3.2.1.1"
    assert xrefs.ec("3.2.1.-") is None            # partial EC not resolvable
    assert "QuickGO/term/GO:0005975" in xrefs.go("GO:0005975")
    assert xrefs.arg("ARO:3000001").endswith("/ontology/3000001")
    assert xrefs.arg("not-an-aro") is None


def test_dispatch():
    assert xrefs.xref("cazy", "GT2") == xrefs.cazy("GT2")
    assert xrefs.xref("ec", "2.4.1.1") == xrefs.ec("2.4.1.1")
    assert xrefs.xref("unknown", "x") is None
