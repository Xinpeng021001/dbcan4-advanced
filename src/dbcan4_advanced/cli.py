#!/usr/bin/env python3
"""dbcan4 — command-line interface for the dbCAN4-advanced annotation engine.

Subcommands
-----------
  dbcan4 embed     FASTA -> ESM-C embeddings (.npz)                 [GPU]
  dbcan4 infer     embeddings -> label-free family calls (TSVs)     [CPU/GPU]
  dbcan4 annotate  FASTA -> family calls in one step (embed+infer)  [GPU]
  dbcan4 run       FASTA -> full Nextflow pipeline (baseline + advanced + features)
                   -> standard v1.1 output contract  [--serve to ingest + launch web UI]
  dbcan4 info      show resolved asset paths + versions

The heavy scientific steps shell out to the validated stage scripts
(embed_esmc.py, infer_esmc.py) shipped alongside the pipeline, so the engine
code that already runs on GPU is reused verbatim rather than reimplemented.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__
from .config import Assets, resolve_assets


def _run(cmd: list[str], **kw) -> int:
    print(f"[dbcan4] $ {' '.join(str(c) for c in cmd)}", file=sys.stderr, flush=True)
    return subprocess.call([str(c) for c in cmd], **kw)


def cmd_embed(args, A: Assets) -> int:
    out_prefix = args.out_prefix or (Path(args.fasta).with_suffix("").name + ".esmc")
    rc = _run([A.engine_python, A.script("embed_esmc.py"),
               "--fasta", args.fasta, "--out-prefix", out_prefix,
               "--model", args.model, "--nshards", "1", "--shard", "0"])
    if rc == 0:
        shard = Path(f"{out_prefix}.shard0.npz")
        canon = Path(f"{out_prefix}.npz")
        if shard.exists():
            shard.replace(canon)
        print(f"[dbcan4] embeddings -> {canon}")
    return rc


def cmd_infer(args, A: Assets) -> int:
    cmd = [A.engine_python, A.script("infer_esmc.py"),
           "--emb", args.emb,
           "--ref-prefix", args.ref_prefix or A.ref_emb_prefix,
           "--k", str(args.k),
           "--out-knn", args.out_knn,
           "--out-centroid", args.out_centroid]
    heads = args.heads or A.heads_pt
    proj_ref = args.proj_ref or A.proj_ref_npz
    if args.out_contrastive and heads and Path(heads).exists():
        cmd += ["--heads", heads, "--proj-ref", proj_ref, "--out-contrastive", args.out_contrastive]
    return _run(cmd)


def cmd_annotate(args, A: Assets) -> int:
    """FASTA -> family calls in one step (embed then label-free infer)."""
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    prefix = str(outdir / "query.esmc")
    # embed
    e = argparse.Namespace(fasta=args.fasta, out_prefix=prefix, model=args.model)
    rc = cmd_embed(e, A)
    if rc != 0:
        return rc
    # infer (label-free)
    i = argparse.Namespace(
        emb=f"{prefix}.npz", ref_prefix=args.ref_prefix, k=args.k,
        heads=args.heads, proj_ref=args.proj_ref,
        out_knn=str(outdir / "ESM-C-kNN.raw.tsv"),
        out_centroid=str(outdir / "ESM-C-centroid.raw.tsv"),
        out_contrastive=str(outdir / "ESM-C-contrastive.raw.tsv"))
    rc = cmd_infer(i, A)
    if rc == 0:
        print(f"[dbcan4] annotate complete -> {outdir}/ESM-C-*.raw.tsv")
    return rc


def cmd_run(args, A: Assets) -> int:
    """Run the full Nextflow pipeline; optionally ingest + serve."""
    nf_main = A.nf_dir / "main.nf"
    if not nf_main.exists():
        print(f"[dbcan4] ERROR: pipeline not found at {nf_main}", file=sys.stderr)
        return 2
    # build a samplesheet from a single FASTA if --input not given
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    samplesheet = args.input
    if not samplesheet:
        samplesheet = str(outdir / "samplesheet.csv")
        sample = args.sample or Path(args.fasta).with_suffix("").name
        with open(samplesheet, "w") as fh:
            fh.write("sample,faa\n")
            fh.write(f"{sample},{os.path.abspath(args.fasta)}\n")
        print(f"[dbcan4] wrote samplesheet -> {samplesheet}")
    nxf = shutil.which("nextflow") or "nextflow"
    cmd = [nxf, "run", str(nf_main), "-profile", args.profile,
           "--input", samplesheet, "--outdir", str(outdir)]
    if args.stub:
        cmd.append("-stub-run")
    if args.resume:
        cmd.append("-resume")
    rc = _run(cmd)
    if rc != 0:
        return rc
    manifest = outdir / "cazyme_advanced" / "manifest.json"
    print(f"[dbcan4] pipeline complete -> {manifest}")
    if args.serve:
        return _serve(args, A, outdir, manifest)
    return 0


def _serve(args, A: Assets, outdir: Path, manifest: Path) -> int:
    """Ingest the published contract into BioForge SQLite + launch the web UI."""
    db = os.path.abspath(args.db or (outdir / "dbcan4.db"))
    env = dict(os.environ, DATABASE_URL=f"sqlite:///{db}")
    funcscan = outdir / "funcscan"
    print(f"[dbcan4] ingesting into {db}")
    if _run(["alembic", "upgrade", "head"], env=env, cwd=args.biodb) != 0:
        print("[dbcan4] WARN: alembic upgrade failed", file=sys.stderr)
    if funcscan.exists():
        _run(["bioforge-ingest", str(funcscan)], env=env)
    _run(["bioforge-ingest-advanced", str(manifest)], env=env)
    print(f"[dbcan4] serving web UI on http://{args.host}:{args.port}  (Ctrl-C to stop)")
    return _run(["uvicorn", "bioforge.api.main:app", "--host", args.host,
                 "--port", str(args.port)], env=env)


def cmd_info(args, A: Assets) -> int:
    print(f"dbcan4-advanced {__version__}")
    print(f"  package dir : {Path(__file__).resolve().parent}")
    print(f"  pipeline    : {A.nf_dir}")
    print(f"  scripts     : {A.scripts_dir}")
    print(f"  ref index   : {A.ref_emb_prefix}")
    print(f"  heads.pt    : {A.heads_pt}  ({'exists' if A.heads_pt and Path(A.heads_pt).exists() else 'MISSING'})")
    print(f"  proj_ref    : {A.proj_ref_npz}")
    print(f"  nextflow    : {shutil.which('nextflow') or 'NOT ON PATH'}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="dbcan4", description="dbCAN4-advanced annotation engine")
    p.add_argument("--version", action="version", version=f"dbcan4-advanced {__version__}")
    p.add_argument("--assets", default=None, help="override asset root (dir with nf/ + results/)")
    sub = p.add_subparsers(dest="cmd", required=True)

    pe = sub.add_parser("embed", help="FASTA -> ESM-C embeddings")
    pe.add_argument("fasta"); pe.add_argument("--out-prefix", default=None)
    pe.add_argument("--model", default="esmc_600m"); pe.set_defaults(fn=cmd_embed)

    pi = sub.add_parser("infer", help="embeddings -> label-free family calls")
    pi.add_argument("--emb", required=True)
    pi.add_argument("--ref-prefix", default=None); pi.add_argument("--k", type=int, default=15)
    pi.add_argument("--heads", default=None); pi.add_argument("--proj-ref", default=None)
    pi.add_argument("--out-knn", default="ESM-C-kNN.raw.tsv")
    pi.add_argument("--out-centroid", default="ESM-C-centroid.raw.tsv")
    pi.add_argument("--out-contrastive", default="ESM-C-contrastive.raw.tsv")
    pi.set_defaults(fn=cmd_infer)

    pa = sub.add_parser("annotate", help="FASTA -> family calls (embed + infer)")
    pa.add_argument("fasta"); pa.add_argument("--outdir", default="dbcan4_annotate")
    pa.add_argument("--model", default="esmc_600m"); pa.add_argument("--k", type=int, default=15)
    pa.add_argument("--ref-prefix", default=None); pa.add_argument("--heads", default=None)
    pa.add_argument("--proj-ref", default=None); pa.set_defaults(fn=cmd_annotate)

    pr = sub.add_parser("run", help="full Nextflow pipeline (baseline + advanced + features)")
    pr.add_argument("--fasta", default=None, help="single protein FASTA (auto-builds samplesheet)")
    pr.add_argument("--input", default=None, help="samplesheet.csv (sample,faa); overrides --fasta")
    pr.add_argument("--sample", default=None)
    pr.add_argument("--outdir", default="dbcan4_results")
    pr.add_argument("--profile", default="met")
    pr.add_argument("--stub", action="store_true", help="add -stub-run (prove DAG, no tools/GPU)")
    pr.add_argument("--resume", action="store_true")
    pr.add_argument("--serve", action="store_true", help="ingest + launch web UI after the run")
    pr.add_argument("--db", default=None); pr.add_argument("--biodb", default=None, help="biodb repo dir (for alembic)")
    pr.add_argument("--host", default="127.0.0.1"); pr.add_argument("--port", type=int, default=8000)
    pr.set_defaults(fn=cmd_run)

    ps = sub.add_parser("info", help="show resolved asset paths + versions")
    ps.set_defaults(fn=cmd_info)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    A = resolve_assets(args.assets)
    return args.fn(args, A)


if __name__ == "__main__":
    raise SystemExit(main())
