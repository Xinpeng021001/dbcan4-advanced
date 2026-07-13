"""Discover per-sample input files inside a dbCAN4 / funcscan-style output dir.

dbCAN4 is a **fungal, protein-input** tool: the pipeline takes a protein FASTA
and there is no genome assembly, gene prediction, or Prokka step (Prokka is a
bacterial/archaeal annotator and has no place here). So a sample is anchored on
its **dbCAN overview + protein FASTA**, and there is **no GFF** in the default
path.

A GFF is treated as an **optional** input: if the user supplies genomic
coordinates (``annotation/prokka/…/<sample>.gff`` or a plain ``*.gff``), the
sample is discovered with those coordinates; otherwise genes are built directly
from the protein FASTA records (each protein is its own coordinate-free unit,
residue span 1..L).

The entry point accepts either the pipeline output root (containing ``results/``)
or the ``results/`` dir itself.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class SampleInputs:
    sample_key: str
    gff_path: Path | None = None          # OPTIONAL: only when the user provides coordinates
    dbcan_path: Path | None = None
    interpro_path: Path | None = None
    faa_path: Path | None = None          # protein FASTA (interproscan/<s>_cleaned.faa)
    cgc_path: Path | None = None          # CAZyme gene clusters (cgc/<s>_cgc.gff)
    cgc_std_path: Path | None = None      # CGC metadata (cgc/<s>_cgc_standard_out.tsv)
    substrate_path: Path | None = None    # substrate/<s>_substrate_prediction.tsv
    arg_path: Path | None = None          # run-level hAMRonization combined report


def _results_root(path: str | Path) -> Path:
    p = Path(path)
    if (p / "results").is_dir():
        return p / "results"
    return p


def discover_samples(pipeline_output: str | Path) -> list[SampleInputs]:
    """Discover samples, preferring the protein-FASTA anchor (fungal default).

    Two discovery modes, unioned by sample name (protein-mode wins on overlap):

    1. **Protein mode (default, no GFF):** every sample with a dbCAN overview
       ``cazyme/dbcan/cazyme_annotation/<sample>_overview.tsv``. Genes are built
       from the cleaned protein FASTA. ``gff_path`` stays ``None``.
    2. **Coordinate mode (optional GFF):** if the user supplied a Prokka-style
       GFF under ``annotation/prokka/…``, that sample carries ``gff_path`` so
       genes get genomic coordinates. Purely opt-in.
    """
    root = _results_root(pipeline_output)
    by_name: dict[str, SampleInputs] = {}

    # --- mode 1: protein-FASTA / dbCAN anchor (the fungal default) ---
    dbcan_dir = root / "cazyme" / "dbcan" / "cazyme_annotation"
    if dbcan_dir.is_dir():
        for ov in sorted(dbcan_dir.glob("*_overview.tsv")):
            name = ov.name[: -len("_overview.tsv")]
            by_name[name] = SampleInputs(
                sample_key=name,
                gff_path=None,
                dbcan_path=ov,
                interpro_path=_interpro_for(root, name),
                faa_path=_faa_for(root, name),
                cgc_path=_cgc_for(root, name),
                cgc_std_path=_cgc_std_for(root, name),
                substrate_path=_substrate_for(root, name),
                arg_path=_arg_for(root),
            )

    # --- mode 2: OPTIONAL user-provided GFF coordinates ---
    prokka_root = root / "annotation" / "prokka"
    if prokka_root.is_dir():
        category_dirs = [d for d in ("all", "long") if (prokka_root / d).is_dir()] or ["."]
        seen: set[str] = set()
        for cat in category_dirs:
            cat_dir = prokka_root if cat == "." else prokka_root / cat
            for sample_dir in sorted(p for p in cat_dir.iterdir() if p.is_dir()):
                name = sample_dir.name
                if name in seen:
                    continue
                gff = sample_dir / f"{name}.gff"
                if not gff.exists():
                    hits = list(sample_dir.glob("*.gff"))
                    if not hits:
                        continue
                    gff = hits[0]
                seen.add(name)
                if name in by_name:
                    # already discovered in protein mode — just attach the coords
                    by_name[name].gff_path = gff
                else:
                    by_name[name] = SampleInputs(
                        sample_key=name,
                        gff_path=gff,
                        dbcan_path=_dbcan_for(root, name),
                        interpro_path=_interpro_for(root, name),
                        faa_path=_faa_for(root, name),
                        cgc_path=_cgc_for(root, name),
                        cgc_std_path=_cgc_std_for(root, name),
                        substrate_path=_substrate_for(root, name),
                        arg_path=_arg_for(root),
                    )

    return [by_name[k] for k in sorted(by_name)]


def _arg_for(root: Path) -> Path | None:
    """The run-level hAMRonization combined ARG report (shared by all samples)."""
    return _exists(root / "arg" / "hamronization" / "hamronization_combined_report.tsv")


def _exists(p: Path) -> Path | None:
    return p if p.exists() else None


def _dbcan_for(root: Path, sample: str) -> Path | None:
    return _exists(root / "cazyme" / "dbcan" / "cazyme_annotation" / f"{sample}_overview.tsv")


def _interpro_for(root: Path, sample: str) -> Path | None:
    return _exists(
        root / "protein_annotation" / "interproscan" / f"{sample}_interproscan_faa.tsv"
    )


def _faa_for(root: Path, sample: str) -> Path | None:
    return _exists(
        root / "protein_annotation" / "interproscan" / f"{sample}_cleaned.faa"
    )


def _cgc_for(root: Path, sample: str) -> Path | None:
    return _exists(root / "cazyme" / "dbcan" / "cgc" / f"{sample}_cgc.gff")


def _cgc_std_for(root: Path, sample: str) -> Path | None:
    return _exists(root / "cazyme" / "dbcan" / "cgc" / f"{sample}_cgc_standard_out.tsv")


def _substrate_for(root: Path, sample: str) -> Path | None:
    return _exists(
        root / "cazyme" / "dbcan" / "substrate" / f"{sample}_substrate_prediction.tsv"
    )
