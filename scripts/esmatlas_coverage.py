#!/usr/bin/env python
"""
Evaluate sequence-level coverage of a sample of Mycocosm fungal CAZyme
proteins against the ESM Metagenomic Atlas "highquality_clust30" sequence
set, using mmseqs2 easy-search.

Usage:
    python esmatlas_coverage.py --query_faa sample.faa \
        --esmatlas_fasta /array1/xinpeng/esmatlas/highquality_clust30.fasta \
        --out_dir /array1/xinpeng/dbcan4-advanced/structure/esmatlas_coverage \
        --threads 64 --min_id 0.3 --min_cov 0.5
"""
import argparse
import json
import os
import subprocess
import sys


def sh(cmd):
    print(f"[esmatlas] $ {cmd}", file=sys.stderr)
    subprocess.run(cmd, shell=True, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query_faa", required=True)
    ap.add_argument("--esmatlas_fasta", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--threads", type=int, default=64)
    ap.add_argument("--min_id", type=float, default=0.3, help="min fraction seq identity")
    ap.add_argument("--min_cov", type=float, default=0.5, help="min coverage")
    ap.add_argument("--sensitivity", type=float, default=5.7)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    tmp_dir = os.path.join(args.out_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    result_m8 = os.path.join(args.out_dir, "esmatlas_hits.m8")

    cmd = (
        f"mmseqs easy-search {args.query_faa} {args.esmatlas_fasta} {result_m8} {tmp_dir} "
        f"--threads {args.threads} --min-seq-id {args.min_id} -c {args.min_cov} "
        f"-s {args.sensitivity} --format-output query,target,pident,alnlen,evalue,bits,qcov,tcov"
    )
    sh(cmd)

    # count query hits
    n_query_total = sum(1 for line in open(args.query_faa) if line.startswith(">"))
    hit_queries = set()
    best_hit = {}
    with open(result_m8) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            q, t, pident = parts[0], parts[1], float(parts[2])
            hit_queries.add(q)
            if q not in best_hit or pident > best_hit[q][1]:
                best_hit[q] = (t, pident)

    n_hit = len(hit_queries)
    coverage_fraction = n_hit / n_query_total if n_query_total else 0.0

    summary = {
        "n_query_total": n_query_total,
        "n_query_with_esmatlas_hit": n_hit,
        "coverage_fraction": coverage_fraction,
        "min_seq_id": args.min_id,
        "min_cov": args.min_cov,
    }
    with open(os.path.join(args.out_dir, "coverage_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))

    with open(os.path.join(args.out_dir, "best_hits.tsv"), "w") as fh:
        fh.write("query_id\tbest_esmatlas_hit\tpident\n")
        for q, (t, pident) in best_hit.items():
            fh.write(f"{q}\t{t}\t{pident}\n")


if __name__ == "__main__":
    main()
