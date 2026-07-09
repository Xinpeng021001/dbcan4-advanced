#!/usr/bin/env python
"""
Build a sequence-based ID mapping between CAZyme3D_id50 structures (accessions
like A0A3D9I5V7, QLD87080.1) and Mycocosm fungal CAZyme proteins (from
all_genome/*/*/overview.tsv + uniInput.faa), and report coverage.

Strategy: exact-sequence (MD5) matching. CAZyme3D_id50 accessions are
UniProt/RefSeq-style, NOT Mycocosm JGI protein IDs, so accession-string
matching is not viable; sequence identity is the robust bridge.

Steps:
  1. Walk all_genome/*/*/overview.tsv (dbCAN CAZyme calls: any of
     diamond/hmm/dbCAN_sub hit == candidate CAZyme) to get the set of
     "Gene ID"s per genome that are CAZyme candidates.
  2. Look up each candidate's AA sequence in the matching uniInput.faa,
     compute MD5.
  3. Extract AA sequences from ALL CAZyme3D_id50 PDB structures via
     `foldseek structureto3didescriptor` (bulk, one shot, chain-name-mode 1),
     compute MD5 per accession.
  4. Report: total unique Mycocosm CAZyme proteins (by MD5), how many have an
     exact-sequence match to a CAZyme3D_id50 structure, coverage fraction.

Usage:
    python cazyme3d_mapping.py --genome_root /array1/xinpeng/all_genome \
        --cazyme3d_dir /array1/xinpeng/cazyme3d/extracted/cazyme_id50 \
        --foldseek /usr/local/bin/foldseek \
        --out_prefix /array1/xinpeng/dbcan4-advanced/structure/cazyme3d_mapping/mapping
"""
import argparse
import glob
import hashlib
import json
import os
import subprocess
import sys
import time


def md5_of(seq):
    return hashlib.md5(seq.encode()).hexdigest()


def read_fasta_dict(path):
    """Return dict gene_id -> seq for a fasta file, keyed by the FIRST token
    after '>' (whitespace-split)."""
    d = {}
    name = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    d[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        d[name] = "".join(seq)
    return d


def collect_mycocosm_cazymes(genome_root):
    """Walk all_genome/*/*/overview.tsv + uniInput.faa. Returns:
    - md5_to_genes: dict md5 -> list of (genome, gene_id, families)
    - n_candidate_rows: total overview.tsv candidate rows (any tool hit)
    - n_genomes: number of genome dirs processed
    """
    md5_to_genes = {}
    n_candidate_rows = 0
    n_genomes = 0
    n_seq_missing = 0

    overview_files = glob.glob(os.path.join(genome_root, "*", "*", "overview.tsv"))
    print(f"[map] found {len(overview_files)} overview.tsv files", file=sys.stderr)

    for i, ov in enumerate(overview_files):
        gdir = os.path.dirname(ov)
        genome = os.path.basename(gdir)
        uni_path = os.path.join(gdir, "uniInput.faa")
        if not os.path.exists(uni_path):
            continue
        # parse overview.tsv candidates: any row with #ofTools column >=1 (i.e. present at all,
        # since overview.tsv itself is already restricted to hits)
        candidate_ids = []
        families = {}
        with open(ov) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                gid_idx = header.index("Gene ID")
            except ValueError:
                gid_idx = 0
            try:
                rec_idx = header.index("Recommend Results")
            except ValueError:
                rec_idx = None
            try:
                dia_idx = header.index("DIAMOND")
            except ValueError:
                dia_idx = None
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= gid_idx:
                    continue
                gid = parts[gid_idx]
                candidate_ids.append(gid)
                fam = None
                if rec_idx is not None and rec_idx < len(parts) and parts[rec_idx] != "-":
                    fam = parts[rec_idx]
                elif dia_idx is not None and dia_idx < len(parts) and parts[dia_idx] != "-":
                    fam = parts[dia_idx]
                families[gid] = fam

        n_candidate_rows += len(candidate_ids)
        seqs = read_fasta_dict(uni_path)
        for gid in candidate_ids:
            seq = seqs.get(gid)
            if seq is None:
                n_seq_missing += 1
                continue
            m = md5_of(seq.upper())
            md5_to_genes.setdefault(m, []).append((genome, gid, families.get(gid)))

        n_genomes += 1
        if (i + 1) % 200 == 0:
            print(f"[map] processed {i+1}/{len(overview_files)} genomes, "
                  f"{len(md5_to_genes)} unique MD5s so far", file=sys.stderr)

    print(f"[map] done: {n_genomes} genomes, {n_candidate_rows} candidate rows, "
          f"{n_seq_missing} seq-lookup misses, {len(md5_to_genes)} unique MD5s",
          file=sys.stderr)
    return md5_to_genes, n_candidate_rows, n_genomes


def extract_cazyme3d_md5s(cazyme3d_dir, foldseek_bin, tmp_tsv):
    """Bulk-run foldseek structureto3didescriptor over the whole directory of
    CAZyme3D_id50 PDBs; parse AA sequences and compute MD5 per accession."""
    cmd = [foldseek_bin, "structureto3didescriptor", "--threads", "32",
           "--chain-name-mode", "1", cazyme3d_dir, tmp_tsv]
    print(f"[map] running: {' '.join(cmd)}", file=sys.stderr)
    t0 = time.time()
    subprocess.run(cmd, check=True)
    print(f"[map] foldseek bulk extraction done in {time.time()-t0:.1f}s", file=sys.stderr)

    md5_to_acc = {}
    n = 0
    with open(tmp_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            desc, aa_seq = parts[0], parts[1]
            acc = desc.split(" ")[0]
            if acc.endswith(".pdb_A"):
                acc = acc[: -len(".pdb_A")]
            elif acc.endswith("_A"):
                acc = acc[: -2]
            m = md5_of(aa_seq.upper())
            md5_to_acc.setdefault(m, []).append(acc)
            n += 1
    print(f"[map] parsed {n} CAZyme3D_id50 structures -> {len(md5_to_acc)} unique MD5s",
          file=sys.stderr)
    return md5_to_acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome_root", default="/array1/xinpeng/all_genome")
    ap.add_argument("--cazyme3d_dir", default="/array1/xinpeng/cazyme3d/extracted/cazyme_id50")
    ap.add_argument("--foldseek", default="/usr/local/bin/foldseek")
    ap.add_argument("--out_prefix", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)

    md5_to_genes, n_candidate_rows, n_genomes = collect_mycocosm_cazymes(args.genome_root)

    tmp_tsv = args.out_prefix + "_foldseek_3di.tsv"
    md5_to_acc = extract_cazyme3d_md5s(args.cazyme3d_dir, args.foldseek, tmp_tsv)

    matched_md5s = set(md5_to_genes.keys()) & set(md5_to_acc.keys())
    n_unique_cazyme_md5 = len(md5_to_genes)
    n_matched_md5 = len(matched_md5s)
    n_matched_gene_rows = sum(len(md5_to_genes[m]) for m in matched_md5s)

    coverage_by_md5 = n_matched_md5 / n_unique_cazyme_md5 if n_unique_cazyme_md5 else 0.0
    coverage_by_rows = n_matched_gene_rows / n_candidate_rows if n_candidate_rows else 0.0

    summary = {
        "n_genomes_processed": n_genomes,
        "n_cazyme_candidate_rows_total": n_candidate_rows,
        "n_unique_cazyme_md5": n_unique_cazyme_md5,
        "n_cazyme3d_structures": sum(len(v) for v in md5_to_acc.values()),
        "n_cazyme3d_unique_md5": len(md5_to_acc),
        "n_matched_unique_md5": n_matched_md5,
        "n_matched_gene_rows": n_matched_gene_rows,
        "coverage_fraction_by_unique_md5": coverage_by_md5,
        "coverage_fraction_by_gene_rows": coverage_by_rows,
    }
    with open(args.out_prefix + "_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))

    # write the mapping table: gene rows with a CAZyme3D structure match
    with open(args.out_prefix + "_matched_genes.tsv", "w") as fh:
        fh.write("genome\tgene_id\tfamily\tseq_md5\tcazyme3d_accessions\n")
        for m in matched_md5s:
            accs = ",".join(md5_to_acc[m])
            for genome, gid, fam in md5_to_genes[m]:
                fh.write(f"{genome}\t{gid}\t{fam or '-'}\t{m}\t{accs}\n")

    print(f"[map] wrote summary + matched_genes.tsv to {args.out_prefix}_*", file=sys.stderr)


if __name__ == "__main__":
    main()
