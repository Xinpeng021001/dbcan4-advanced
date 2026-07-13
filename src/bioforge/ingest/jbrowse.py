"""Prepare JBrowse 2 track files per sample, wired into the ingester.

Steps (equivalent to `jbrowse sort-gff | bgzip | tabix`):
  1. sort GFF features by (seqid, start)  — required for tabix
  2. bgzip-compress                        — via pysam (no system htslib needed)
  3. tabix-index (.tbi)                     — via pysam
  4. write a per-sample JBrowse 2 config.json referencing the track

Using pysam for bgzip/tabix means the one-command setup has no hard dependency on
the Node `@jbrowse/cli` or a system htslib build. The optional `jbrowse text-index`
(name search) step is skipped gracefully if the CLI isn't installed.

JBrowse reads the .gz/.tbi over HTTP range requests, so these files just need to
be served statically — no JBrowse backend.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

from .parse_gff import parse_gff3

# pysam is imported lazily inside prepare_sample_track: it carries a compiled
# htslib and is only needed when actually building tracks. Keeping it out of
# module import means `bioforge-ingest --no-tracks` and the API run with no
# pysam installed at all.


def _sorted_gff_lines(gff_path: Path) -> list[str]:
    """Return GFF3 data lines sorted by (seqid, start), with a version pragma."""
    rows = []
    for feat in parse_gff3(gff_path):
        attrs = ";".join(f"{k}={v}" for k, v in feat.attributes.items())
        rows.append(
            (
                feat.seqid,
                feat.start,
                "\t".join(
                    [
                        feat.seqid,
                        feat.source,
                        feat.feature_type,
                        str(feat.start),
                        str(feat.end),
                        ".",
                        feat.strand or ".",
                        ".",
                        attrs,
                    ]
                ),
            )
        )
    rows.sort(key=lambda r: (r[0], r[1]))
    return ["##gff-version 3"] + [r[2] for r in rows]


def prepare_sample_track(
    gff_path: str | Path,
    sample_key: str,
    tracks_root: str | Path,
    contigs: list[str] | None = None,
) -> dict:
    """Build sorted/bgzipped/tabixed GFF + config.json for one sample.

    Returns a dict describing the track (relative paths, assembly name). The
    caller stores the relative track dir on the Sample row.
    """
    import pysam  # lazy: only needed when tracks are actually built

    gff_path = Path(gff_path)
    tracks_root = Path(tracks_root)
    out_dir = tracks_root / sample_key
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1 + 2: sorted, bgzip-compressed GFF
    sorted_plain = out_dir / f"{sample_key}.sorted.gff"
    sorted_plain.write_text("\n".join(_sorted_gff_lines(gff_path)) + "\n")
    gz_path = out_dir / f"{sample_key}.sorted.gff.gz"
    pysam.tabix_compress(str(sorted_plain), str(gz_path), force=True)
    sorted_plain.unlink(missing_ok=True)

    # 3: tabix index
    pysam.tabix_index(str(gz_path), preset="gff", force=True)

    # 4 (optional): text index for name search, only if jbrowse CLI present
    text_indexed = _maybe_text_index(out_dir, sample_key)

    # config.json — one assembly + one feature track for this sample
    config = _build_config(sample_key, gz_path.name, contigs or [])
    (out_dir / "config.json").write_text(json.dumps(config, indent=2))

    return {
        "track_dir": sample_key,  # relative to tracks_root
        "gff_gz": gz_path.name,
        "tbi": gz_path.name + ".tbi",
        "config": "config.json",
        "text_indexed": text_indexed,
    }


def _maybe_text_index(out_dir: Path, sample_key: str) -> bool:
    if shutil.which("jbrowse") is None:
        return False
    try:
        subprocess.run(
            ["jbrowse", "text-index", "--out", str(out_dir), "--force"],
            cwd=out_dir,
            check=True,
            capture_output=True,
        )
        return True
    except (subprocess.CalledProcessError, OSError):
        return False


def _build_config(sample_key: str, gff_gz_name: str, contigs: list[str]) -> dict:
    """Minimal JBrowse 2 config: one assembly, one GFF feature track.

    We use an inline sequence adapter placeholder assembly. For real genomes the
    ingester would also bgzip/faidx the FASTA; for the browsable MVP the feature
    track over the contigs is what matters. Region list drives the default view.
    """
    assembly_name = sample_key
    return {
        "assemblies": [
            {
                "name": assembly_name,
                "sequence": {
                    "type": "ReferenceSequenceTrack",
                    "trackId": f"{assembly_name}-ref",
                    "adapter": {
                        "type": "FromConfigSequenceAdapter",
                        "features": [
                            {
                                "refName": c,
                                "uniqueId": f"{assembly_name}-{c}",
                                "start": 0,
                                "end": 100000,
                            }
                            for c in (contigs or [assembly_name])
                        ],
                    },
                },
            }
        ],
        "tracks": [
            {
                "type": "FeatureTrack",
                "trackId": f"{sample_key}-annotation",
                "name": f"{sample_key} annotation",
                "assemblyNames": [assembly_name],
                "adapter": {
                    "type": "Gff3TabixAdapter",
                    "gffGzLocation": {"uri": gff_gz_name, "locationType": "UriLocation"},
                    "index": {
                        "location": {
                            "uri": gff_gz_name + ".tbi",
                            "locationType": "UriLocation",
                        }
                    },
                },
            }
        ],
        "defaultSession": {
            "name": f"{sample_key} view",
            "view": {
                "id": "linearGenomeView",
                "type": "LinearGenomeView",
                "tracks": [
                    {
                        "type": "FeatureTrack",
                        "configuration": f"{sample_key}-annotation",
                        "displays": [
                            {
                                "type": "LinearBasicDisplay",
                                "configuration": f"{sample_key}-annotation-LinearBasicDisplay",
                            }
                        ],
                    }
                ],
            },
        },
    }
