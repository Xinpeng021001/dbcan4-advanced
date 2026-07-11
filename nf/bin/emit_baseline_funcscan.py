#!/usr/bin/env python3
"""Reshape a run_dbcan V5 protein run into the layout BioForge ingests.

dbCAN4 is a **fungal, protein-input** tool. run_dbcan in --mode protein takes a
protein FASTA and produces family calls — there is no genome, no gene prediction,
and **no Prokka** (Prokka is a bacterial/archaeal annotator). So this reshaper
emits only the two things the ingester actually needs in protein mode:

  <root>/cazyme/dbcan/cazyme_annotation/<sample>_overview.tsv   (family calls)
  <root>/protein_annotation/interproscan/<sample>_cleaned.faa   (protein FASTA)

BioForge's discovery anchors on the dbCAN overview + protein FASTA and builds one
gene per protein (coordinate-free). It does **not** synthesize a Prokka GFF.

A GFF is *optional*: if the user has genuine genomic coordinates for these
proteins, pass ``--gff <file>`` and it is copied to
``annotation/prokka/all/<sample>/<sample>.gff`` so the ingester runs in
coordinate mode. Absent that flag, no GFF is written.

The dbCAN_hmm column is renamed to the literal 'HMMER' that parse_dbcan.py
expects (run_dbcan V5 writes 'dbCAN_hmm').
"""
from __future__ import annotations
import argparse, csv, os, shutil


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overview", required=True, help="run_dbcan overview.tsv")
    ap.add_argument("--faa", required=True, help="input protein FASTA")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True, help="output root to populate")
    ap.add_argument("--gff", default=None,
                    help="OPTIONAL user-provided genomic GFF; if given, copied to the "
                         "prokka-style coordinate anchor so the ingester adds coordinates")
    args = ap.parse_args()

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

    # --- 2. copy protein FASTA into the IPS location so sequences load (clean ids) ---
    ips_dir = os.path.join(root, "protein_annotation", "interproscan")
    os.makedirs(ips_dir, exist_ok=True)
    faa_out = os.path.join(ips_dir, f"{args.sample}_cleaned.faa")
    with open(args.faa) as fi, open(faa_out, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                fo.write(">" + line[1:].strip().split()[0].split("|")[0] + "\n")
            else:
                fo.write(line)

    # --- 3. OPTIONAL: copy a user-provided GFF into the coordinate anchor ---
    gff_msg = "no GFF (protein mode)"
    if args.gff:
        prokka_dir = os.path.join(root, "annotation", "prokka", "all", args.sample)
        os.makedirs(prokka_dir, exist_ok=True)
        gff_out = os.path.join(prokka_dir, f"{args.sample}.gff")
        shutil.copyfile(args.gff, gff_out)
        gff_msg = f"user GFF -> {gff_out}"

    print(f"[baseline] {ov_out}: {len(rows)} rows; faa->{faa_out}; {gff_msg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
