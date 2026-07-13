"""`bioforge-ingest-advanced` — load dbCAN4-advanced predictions as a new release.

Usage:
    bioforge-ingest-advanced <manifest.json | cazyme_advanced_dir>
        [--tag LABEL] [--notes TEXT] [--database-url URL]
        [--structures-dir DIR] [--not-current]

Reads the STANDARD OUTPUT CONTRACT the Nextflow pipeline publishes and inserts
its advanced CAZyme calls + per-protein features against the genes the baseline
ingest already loaded, as a new versioned Release. Idempotent: re-running on an
unchanged manifest adds nothing. Requires the baseline release to be ingested
first (advanced calls attach to existing genes).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ..config import database_url as default_db_url
from ..config import tracks_dir as default_tracks_dir
from ..db import make_engine, make_session_factory
from .loader_advanced import ingest_advanced_manifest


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="bioforge-ingest-advanced",
        description="Ingest dbCAN4-advanced standard outputs into BioForge as a new release.",
    )
    p.add_argument("manifest", help="Path to cazyme_advanced/manifest.json (or its dir).")
    p.add_argument("--tag", dest="label", default=None, help="Release label override.")
    p.add_argument("--notes", default=None, help="Free-text release notes.")
    p.add_argument("--database-url", default=None, help="Override DATABASE_URL.")
    p.add_argument("--structures-dir", default=None,
                   help="Where to copy served structure PDBs "
                        "(default: <tracks_dir>/../structures).")
    p.add_argument("--not-current", action="store_true",
                   help="Do not mark this the current release.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    man_path = Path(args.manifest)
    if not man_path.exists():
        print(f"error: manifest not found: {man_path}", file=sys.stderr)
        return 2

    url = args.database_url or default_db_url()
    engine = make_engine(url)
    Session = make_session_factory(engine)

    # Default structures dir sits beside the JBrowse tracks under web_static/.
    if args.structures_dir:
        struct_dir = Path(args.structures_dir)
    else:
        struct_dir = default_tracks_dir().parent / "structures"
    struct_dir.mkdir(parents=True, exist_ok=True)

    with Session() as session:
        report = ingest_advanced_manifest(
            session, man_path,
            label=args.label, notes=args.notes,
            structures_dir=struct_dir,
            make_current=not args.not_current,
        )

    r = report
    if r.release_id is None and r.samples_matched == 0 and not r.samples_missing:
        print("No changes: manifest already ingested (idempotent no-op).")
        return 0
    if r.release_id is None and r.samples_missing:
        print(f"error: no baseline samples matched {r.samples_missing}. "
              f"Ingest the baseline (bioforge-ingest) first.", file=sys.stderr)
        return 1
    print(f"Advanced release #{r.release_id} ('{r.release_label}') ingested.")
    print(f"  samples matched: {r.samples_matched}"
          + (f"  (missing: {r.samples_missing})" if r.samples_missing else ""))
    print(f"  advanced CAZyme calls: +{r.advanced_calls_added}")
    for tool, n in sorted(r.per_tool.items()):
        print(f"      {tool:<20} +{n}")
    print(f"  protein features: +{r.features_added}  "
          f"(structures copied: {r.structures_copied})")
    if r.unmatched_calls or r.unmatched_features:
        print(f"  [warn] unmatched to a gene — calls: {r.unmatched_calls}, "
              f"features: {r.unmatched_features}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
