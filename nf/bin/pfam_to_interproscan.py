#!/usr/bin/env python3
"""Derive an InterProScan-format TSV (with GO terms) from the Pfam domains the
pipeline already computes — a self-contained, offline fallback for the real
InterProScan tool.

Why this exists
---------------
BioForge's gene page renders a Gene Ontology card and an InterPro-domains table,
but both are data-gated: they only appear when the funcscan tree contains
`protein_annotation/interproscan/<sample>_interproscan_faa.tsv`
(parsed by `bioforge.ingest.parse_interpro`). Running the real InterProScan is a
tens-of-GB install, which fights the "git clone and run" goal.

Instead we reuse the Pfam/hmmscan domains (`domains.tsv`, produced by
PFAM_DOMAINS / step 4 of dbcan4_workup.sh) and join them against the small
`pfam2go` mapping (InterPro's official Pfam->GO table). That yields real GO terms
for every Pfam-annotated protein with no extra tool and no network at run time.

Output columns are the standard headerless, positional InterProScan v5 TSV:

    1 protein_acc  2 md5  3 length  4 analysis  5 signature_acc  6 signature_desc
    7 start  8 stop  9 score/evalue  10 status  11 date
    12 interpro_acc  13 interpro_desc  14 GO (|-separated)  15 pathways

`interpro_acc`/`interpro_desc` are filled only if an optional Pfam->InterPro map
is supplied (`--pfam2interpro`); otherwise they are '-' (rendered as "—" in the
web UI). The real InterProScan step, when installed, supersedes this file and
fills those columns natively.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

PFAM_RE = re.compile(r"PF\d{5}")


def load_pfam2go(path: Path) -> dict[str, list[str]]:
    """Parse the GO `pfam2go` external2go file into {PFxxxxx: [GO:id, ...]}.

    Lines look like:  `Pfam:PF00001 7tm_1 > GO:G protein-coupled ... ; GO:0004930`
    """
    mapping: dict[str, list[str]] = {}
    if not path or not path.exists():
        return mapping
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("Pfam:"):
                continue
            m = PFAM_RE.search(line)
            go = re.search(r"(GO:\d{7})\s*$", line.strip())
            if not m or not go:
                continue
            mapping.setdefault(m.group(0), [])
            gid = go.group(1)
            if gid not in mapping[m.group(0)]:
                mapping[m.group(0)].append(gid)
    return mapping


def load_pfam2interpro(path: Path | None) -> dict[str, tuple[str, str]]:
    """Optional {PFxxxxx: (IPRxxxxxx, description)} map.

    Accepts a 2- or 3-column TSV: `PFxxxxx <TAB> IPRxxxxxx [<TAB> description]`.
    """
    mapping: dict[str, tuple[str, str]] = {}
    if not path or not Path(path).exists():
        return mapping
    with Path(path).open(encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[0].startswith("PF"):
                continue
            pf = parts[0].split(".")[0]
            ipr = parts[1].strip()
            desc = parts[2].strip() if len(parts) > 2 else "-"
            if ipr.startswith("IPR"):
                mapping[pf] = (ipr, desc or "-")
    return mapping


def read_fasta(path: Path | None) -> dict[str, str]:
    seqs: dict[str, str] = {}
    if not path or not Path(path).exists():
        return seqs
    pid, buf = None, []
    with Path(path).open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if pid:
                    seqs[pid] = "".join(buf)
                pid = line[1:].strip().split()[0].split("|")[0]
                buf = []
            else:
                buf.append(line.strip())
    if pid:
        seqs[pid] = "".join(buf)
    return seqs


def domains_from_contract(path: Path):
    """Yield (protein_id, pf_acc, name, start, end, evalue) from a domains.tsv
    (the pipeline's §2.5 contract file)."""
    with path.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("acc") or "").split(".")[0]
            if not acc.startswith("PF"):
                continue
            yield (row.get("protein_id"), acc, row.get("name") or acc,
                   row.get("start") or "", row.get("end") or "",
                   row.get("evalue") or "-")


def domains_from_domtbl(path: Path):
    """Yield the same tuples straight from an `hmmscan --domtblout` file."""
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 19:
                continue
            tname, tacc = f[0], f[1]
            acc = (tacc.split(".")[0] if tacc and tacc != "-" else tname)
            if not acc.startswith("PF"):
                continue
            yield (f[3], acc, tname, f[17], f[18], f[12])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--domains-tsv", help="pipeline §2.5 domains.tsv (from PFAM_DOMAINS)")
    src.add_argument("--domtbl", help="raw hmmscan --domtblout (alternative to --domains-tsv)")
    ap.add_argument("--pfam2go", required=True, help="GO external2go pfam2go mapping file")
    ap.add_argument("--pfam2interpro", default=None,
                    help="optional PF->IPR TSV to also fill the InterPro accession column")
    ap.add_argument("--faa", default=None, help="protein FASTA (for md5 + length columns)")
    ap.add_argument("--out", required=True, help="output <sample>_interproscan_faa.tsv")
    ap.add_argument("--date", default="01-01-2026",
                    help="value for the (unused-by-ingest) date column")
    a = ap.parse_args()

    p2go = load_pfam2go(Path(a.pfam2go))
    p2ipr = load_pfam2interpro(a.pfam2interpro)
    seqs = read_fasta(Path(a.faa) if a.faa else None)
    if not p2go:
        print(f"[pfam_to_interproscan] WARNING: empty pfam2go map at {a.pfam2go} — "
              "GO column will be blank", file=sys.stderr)

    rows = (domains_from_contract(Path(a.domains_tsv)) if a.domains_tsv
            else domains_from_domtbl(Path(a.domtbl)))

    n_rows = n_go = 0
    with open(a.out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        for pid, pf, name, start, end, evalue in rows:
            if not pid:
                continue
            seq = seqs.get(pid, "")
            md5 = hashlib.md5(seq.encode()).hexdigest() if seq else "-"
            length = str(len(seq)) if seq else "-"
            gos = p2go.get(pf, [])
            go_field = "|".join(gos) if gos else "-"
            ipr, ipr_desc = p2ipr.get(pf, ("-", "-"))
            w.writerow([pid, md5, length, "Pfam", pf, name,
                        start, end, evalue, "T", a.date,
                        ipr, ipr_desc, go_field, "-"])
            n_rows += 1
            n_go += len(gos)
    print(f"[pfam_to_interproscan] {a.out}: {n_rows} Pfam rows, {n_go} GO assignments "
          f"({'with' if p2ipr else 'no'} InterPro-accession map)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
