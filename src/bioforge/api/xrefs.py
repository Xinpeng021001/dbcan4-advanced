"""Resolve accessions to their canonical public database pages.

Every serious bio database cross-links out to the source authorities (UniProt,
InterPro, CAZy, …). This module turns the accessions we already store into
resolvable URLs so the UI can render an outbound "↗" next to each tag. Pure
functions, no dependencies, no network calls — just URL templates.

Registered as a Jinja global (`xref`) in api/main.py so templates/macros can call
``xref('cazy', 'GH13')`` without importing Python.
"""
from __future__ import annotations

import re

# InterPro member-database sub-paths, keyed by the analysis name we store in
# InterproDomain.analysis (lower-cased). Anything not listed falls back to a
# plain InterPro entry lookup, which resolves most signatures anyway.
_MEMBER_DB = {
    "pfam": "pfam",
    "panther": "panther",
    "smart": "smart",
    "prosite": "prosite",
    "prosite_profiles": "profile",
    "prosite_patterns": "prosite",
    "cdd": "cdd",
    "tigrfam": "ncbifam",
    "ncbifam": "ncbifam",
    "superfamily": "ssf",
    "prints": "prints",
    "hamap": "hamap",
    "pirsf": "pirsf",
    "gene3d": "cathgene3d",
    "cath-gene3d": "cathgene3d",
}

_CAZY_FAMILY = re.compile(r"^(GH|GT|PL|CE|AA|CBM|SLH|cohesin|dockerin)", re.IGNORECASE)


def cazy(family: str | None) -> str | None:
    """CAZy family page, e.g. GH13 → cazy.org/GH13.html."""
    if not family:
        return None
    fam = family.strip()
    if not _CAZY_FAMILY.match(fam):
        return None
    return f"http://www.cazy.org/{fam}.html"


def interpro(acc: str | None) -> str | None:
    if not acc or not acc.upper().startswith("IPR"):
        return None
    return f"https://www.ebi.ac.uk/interpro/entry/InterPro/{acc}/"


def signature(analysis: str | None, acc: str | None) -> str | None:
    """A member-database signature (Pfam/PANTHER/…) on the InterPro site."""
    if not acc:
        return None
    sub = _MEMBER_DB.get((analysis or "").strip().lower())
    if sub is None:
        return None
    return f"https://www.ebi.ac.uk/interpro/entry/{sub}/{acc}/"


def ec(number: str | None) -> str | None:
    """ExPASy ENZYME entry, e.g. 3.2.1.1. Skip partial ECs (with '-')."""
    if not number:
        return None
    n = number.strip()
    if not re.match(r"^\d+\.\d+\.\d+\.\d+$", n):
        return None
    return f"https://enzyme.expasy.org/EC/{n}"


def go(go_id: str | None) -> str | None:
    if not go_id:
        return None
    g = go_id.strip()
    if not g.upper().startswith("GO:"):
        return None
    return f"https://www.ebi.ac.uk/QuickGO/term/{g}"


def arg(accession: str | None) -> str | None:
    """CARD ontology (ARO) page, e.g. ARO:3000001 → card.mcmaster.ca/ontology/3000001."""
    if not accession:
        return None
    m = re.match(r"^ARO:?(\d+)$", accession.strip(), re.IGNORECASE)
    if not m:
        return None
    return f"https://card.mcmaster.ca/ontology/{m.group(1)}"


def xref(kind: str, acc: str | None, analysis: str | None = None) -> str | None:
    """Dispatch used by templates: xref('cazy'|'interpro'|'signature'|'ec'|'go'|'arg', acc)."""
    if kind == "cazy":
        return cazy(acc)
    if kind == "interpro":
        return interpro(acc)
    if kind == "signature":
        return signature(analysis, acc)
    if kind == "ec":
        return ec(acc)
    if kind == "go":
        return go(acc)
    if kind == "arg":
        return arg(acc)
    return None
