"""`bioforge-ingest` — the entire update workflow in one command.

Usage:
    bioforge-ingest <path-to-funcscan-output> [--tag LABEL] [--notes TEXT]
                                              [--no-tracks] [--database-url URL]

Idempotent: safe to re-run on the same output (unchanged samples are skipped).
Re-running after a new pipeline run ingests only new/changed samples as a new
versioned release, and refreshes their JBrowse tracks.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ..config import database_url as default_db_url
from ..config import tracks_dir as default_tracks_dir
from ..db import make_engine, make_session_factory
from .loader import ingest_directory


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="bioforge-ingest",
        description="Ingest an nf-core/funcscan (4.0.0+) output directory into BioForge.",
    )
    p.add_argument("pipeline_output", help="Path to funcscan output (root or results/).")
    p.add_argument("--tag", dest="label", default=None, help="Release label.")
    p.add_argument("--notes", default=None, help="Free-text release notes.")
    p.add_argument("--no-tracks", action="store_true", help="Skip JBrowse track prep.")
    p.add_argument(
        "--database-url",
        default=None,
        help="Override DATABASE_URL (else env var, else SQLite dev fallback).",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    out_path = Path(args.pipeline_output)
    if not out_path.exists():
        print(f"error: path not found: {out_path}", file=sys.stderr)
        return 2

    url = args.database_url or default_db_url()
    engine = make_engine(url)
    Session = make_session_factory(engine)

    tracks = None if args.no_tracks else default_tracks_dir()

    with Session() as session:
        report = ingest_directory(
            session,
            out_path,
            label=args.label,
            notes=args.notes,
            tracks_dir=tracks,
        )

    r = report
    if r.release_id is None:
        print("No changes: all samples already ingested (idempotent no-op).")
    else:
        print(f"Release #{r.release_id} ingested.")
    print(
        f"  samples: +{r.samples_added} added, {r.samples_skipped} skipped\n"
        f"  genes: +{r.genes_added}  cazymes: +{r.cazymes_added}  "
        f"domains: +{r.domains_added}\n"
        f"  sequences: +{r.sequences_added}  clusters: +{r.clusters_added}  "
        f"args: +{r.args_added}  go: +{r.go_added}"
    )
    if r.unmatched_cazyme or r.unmatched_domains or r.unmatched_arg:
        print(
            f"  [warn] unmatched to a gene — cazyme: {r.unmatched_cazyme}, "
            f"domains: {r.unmatched_domains}, arg: {r.unmatched_arg}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
