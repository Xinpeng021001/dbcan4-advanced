#!/usr/bin/env python3
"""Reshape a run_dbcan V5 protein run into the funcscan layout BioForge ingests.

BioForge's baseline discovery (bioforge.ingest.discover) keys on a Prokka-style
per-sample GFF plus cazyme/dbcan/cazyme_annotation/<sample>_overview.tsv. run_dbcan
in --mode protein does not emit a genome GFF (no genomic coordinates for bare
proteins), so we synthesize a minimal Prokka-compatible GFF that lays each protein
on one demo contig — every gene ID matches the protein_id in the overview, which is
what the ingester joins on. The dbCAN_hmm column is renamed to the literal 'HMMER'
that parse_dbcan.py expects (run_dbcan V5 writes 'dbCAN_hmm').

This is the in-pipeline replacement for a separate funcscan run: FASTA -> baseline.

BioForge funcscan layout (see bioforge.ingest.discover):
  <root>/annotation/prokka/all/<sample>/<sample>.gff        (sample anchor)
  <root>/cazyme/dbcan/cazyme_annotation/<sample>_overview.tsv
  <root>/protein_annotation/interproscan/<sample>_cleaned.faa
"""
from __future__ import annotations
import argparse, csv, os


def read_fasta_ids_lengths(faa):
    ids, order = {}, []
    cur, seq = None, []
    for line in open(faa):
        if line.startswith(">"):
            if cur is not None:
                ids[cur] = len("".join(seq))
            cur = line[1:].strip().split()[0].split("|")[0]
            order.append(cur); seq = []
        else:
            seq.append(line.strip())
    if cur is not None:
        ids[cur] = len("".join(seq))
    return order, ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overview", required=True, help="run_dbcan overview.tsv")
    ap.add_argument("--faa", required=True, help="input protein FASTA")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True, help="funcscan root to populate")
    ap.add_argument("--contig", default=None, help="synthetic contig name (default <sample>_contig_1)")
    args = ap.parse_args()

    contig = args.contig or f"{args.sample}_contig_1"
    order, lengths = read_fasta_ids_lengths(args.faa)
    root = args.outdir

    # --- 1. rewrite overview.tsv with dbCAN_hmm -> HMMER, strip any |families id suffix ---
    dbcan_dir = os.path.join(root, "cazyme", "dbcan", "cazyme_annotation")
    os.makedirs(dbcan_dir, exist_ok=True)
    ov_out = os.path.join(dbcan_dir, f"{args.sample}_overview.tsv")
    with open(args.overview) as fi:
        rdr = csv.reader(fi, delimiter="\t")
        header = next(rdr)
        header = ["HMMER" if h == "dbCAN_hmm" else h for h in header]
        gid_col = 0  # first column is Gene ID
        rows = []
        for row in rdr:
            if not row:
                continue
            row[gid_col] = row[gid_col].split("|")[0]  # clean id
            rows.append(row)
    with open(ov_out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t"); w.writerow(header); w.writerows(rows)

    # --- 2. synthesize a Prokka-style GFF laying each protein on one contig ---
    prokka_dir = os.path.join(root, "annotation", "prokka", "all", args.sample)
    os.makedirs(prokka_dir, exist_ok=True)
    gff = os.path.join(prokka_dir, f"{args.sample}.gff")
    pos = 1
    with open(gff, "w") as fo:
        fo.write("##gff-version 3\n")
        for pid in order:
            L = lengths.get(pid, 300)
            nt = L * 3
            start, end = pos, pos + nt - 1
            fo.write(f"{contig}\tdbcan4\tCDS\t{start}\t{end}\t.\t+\t0\tID={pid};locus_tag={pid}\n")
            pos = end + 1

    # --- 3. copy protein FASTA into the IPS location so sequences load (clean ids) ---
    ips_dir = os.path.join(root, "protein_annotation", "interproscan")
    os.makedirs(ips_dir, exist_ok=True)
    faa_out = os.path.join(ips_dir, f"{args.sample}_cleaned.faa")
    with open(args.faa) as fi, open(faa_out, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                fo.write(">" + line[1:].strip().split()[0].split("|")[0] + "\n")
            else:
                fo.write(line)
    print(f"[baseline] {ov_out}: {len(rows)} rows; {gff}: {len(order)} CDS on {contig}; faa->{faa_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
