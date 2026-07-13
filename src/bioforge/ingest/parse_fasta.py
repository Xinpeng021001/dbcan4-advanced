"""Parse a protein FASTA (``<sample>_cleaned.faa``) into {gene_key: sequence}.

InterProScan's cleaned protein FASTA uses the Prokka/Bakta locus tag as the
record id (the same key our GFF genes carry), so the sequence joins straight back
to a Gene by ``gene_key``. Only the first whitespace-delimited token of the
header is used as the id, matching how the GFF/dbCAN/IPS parsers key genes.
"""
from __future__ import annotations

from pathlib import Path


def parse_protein_fasta(path: str | Path) -> dict[str, str]:
    """Return {gene_key: protein_sequence} (uppercase, no newlines)."""
    seqs: dict[str, list[str]] = {}
    current: str | None = None
    with Path(path).open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0].strip() or None
                if current is not None:
                    seqs.setdefault(current, [])
            elif current is not None:
                seqs[current].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items() if v}


def to_fasta(records: list[tuple[str, str]], width: int = 60) -> str:
    """Serialise [(header, sequence), …] to a wrapped FASTA string."""
    out: list[str] = []
    for header, seq in records:
        out.append(f">{header}")
        for i in range(0, len(seq), width):
            out.append(seq[i : i + width])
    return "\n".join(out) + ("\n" if out else "")
