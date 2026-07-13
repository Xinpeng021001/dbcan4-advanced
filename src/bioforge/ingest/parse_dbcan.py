"""Parser for dbCAN `cazyme_annotation/<sample>_overview.tsv` (run_dbcan v4).

Real columns (header-keyed, tab-separated):
    Gene ID | EC# | HMMER | dbCAN_sub | DIAMOND | #ofTools

Each of the three tool columns (HMMER, dbCAN_sub, DIAMOND) may contain one or
more CAZy family calls, optionally with coordinates and subfamily suffixes, and
'+'-joined when several families hit one gene, e.g. "CBM48(20-70)+GH9(80-820)".
We normalise each token to its family (strip coords + subfamily suffix) and emit
one record per (gene, tool, family) — matching the schema's "row per predicting
tool call" design.

Header-keyed: we look up columns by name and raise on an unrecognised layout so
a run_dbcan version change fails loudly instead of loading garbage positionally.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

# Column names we depend on (case-insensitive match).
_GENE_COL = "gene id"
_EC_COL = "ec#"
_NTOOLS_COL = "#oftools"
_TOOL_COLS = ["hmmer", "dbcan_sub", "diamond"]

# A family token: letters (GH/GT/PL/CE/CBM/AA/SLH ...) + optional number,
# optionally with a subfamily suffix like _31 or _e123 that we drop.
_FAMILY_RE = re.compile(r"^([A-Za-z]+\d*)")


@dataclass
class CazymeCall:
    gene_key: str
    cazy_family: str
    ec_number: str | None
    tool: str
    n_tools_support: int | None


def _clean_tool_name(raw: str) -> str:
    mapping = {"hmmer": "HMMER", "dbcan_sub": "dbCAN_sub", "diamond": "DIAMOND"}
    return mapping.get(raw.lower(), raw)


def _families_in_cell(cell: str) -> list[str]:
    """Extract normalised family names from one tool cell."""
    cell = (cell or "").strip()
    if cell in ("", "-", "N", "NA"):
        return []
    families: list[str] = []
    for token in cell.split("+"):
        token = token.strip()
        # Drop coordinate annotations like "(12-390)".
        token = re.sub(r"\(.*?\)", "", token).strip()
        if not token or token == "-":
            continue
        m = _FAMILY_RE.match(token)
        if m:
            families.append(m.group(1).upper())
    # De-dup while preserving order.
    seen: set[str] = set()
    out = []
    for f in families:
        if f not in seen:
            seen.add(f)
            out.append(f)
    return out


def parse_dbcan_overview(path: str | Path) -> Iterator[CazymeCall]:
    path = Path(path)
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return
        idx = {name.strip().lower(): i for i, name in enumerate(header)}
        if _GENE_COL not in idx:
            raise ValueError(
                f"dbCAN overview {path} missing '{_GENE_COL}' column; got {header!r}"
            )
        present_tools = [t for t in _TOOL_COLS if t in idx]
        if not present_tools:
            raise ValueError(
                f"dbCAN overview {path} has no recognised tool columns "
                f"({_TOOL_COLS}); got {header!r}"
            )

        for row in reader:
            if not row or len(row) <= idx[_GENE_COL]:
                continue
            gene_key = row[idx[_GENE_COL]].strip()
            if not gene_key:
                continue
            ec = row[idx[_EC_COL]].strip() if _EC_COL in idx else None
            if ec in ("-", ""):
                ec = None
            n_tools = None
            if _NTOOLS_COL in idx and idx[_NTOOLS_COL] < len(row):
                try:
                    n_tools = int(row[idx[_NTOOLS_COL]].strip())
                except ValueError:
                    n_tools = None
            for tool in present_tools:
                col_i = idx[tool]
                cell = row[col_i] if col_i < len(row) else ""
                for fam in _families_in_cell(cell):
                    yield CazymeCall(
                        gene_key=gene_key,
                        cazy_family=fam,
                        ec_number=ec,
                        tool=_clean_tool_name(tool),
                        n_tools_support=n_tools,
                    )
