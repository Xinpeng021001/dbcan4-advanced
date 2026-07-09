#!/usr/bin/env python
"""
Fetch AA sequences for a list of protein_ids (format: jgi-<Genome>-<...>) from
all_genome/<phylum>/<Genome>/uniInput.faa files, and write a single FASTA.
"""
import argparse
import glob
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_tsv", required=True)
    ap.add_argument("--genome_root", default="/array1/xinpeng/all_genome")
    ap.add_argument("--out_faa", required=True)
    args = ap.parse_args()

    ids = []
    with open(args.sample_tsv) as fh:
        header = fh.readline()
        for line in fh:
            pid = line.split("\t")[0]
            ids.append(pid)

    # group by genome (second field after splitting on '-')
    genome_to_ids = {}
    for pid in ids:
        parts = pid.split("-")
        genome = parts[1] if len(parts) > 1 else None
        genome_to_ids.setdefault(genome, []).append(pid)

    print(f"[fetch] {len(ids)} ids across {len(genome_to_ids)} genomes")

    # build genome -> uniInput.faa path lookup
    all_dirs = glob.glob(os.path.join(args.genome_root, "*", "*"))
    genome_path = {}
    for d in all_dirs:
        g = os.path.basename(d)
        genome_path[g] = os.path.join(d, "uniInput.faa")

    n_found = 0
    n_missing_genome = 0
    n_missing_seq = 0
    with open(args.out_faa, "w") as out:
        for genome, gids in genome_to_ids.items():
            path = genome_path.get(genome)
            if path is None or not os.path.exists(path):
                n_missing_genome += len(gids)
                continue
            wanted = set(gids)
            found_here = {}
            name = None
            seq = []
            with open(path) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    if line.startswith(">"):
                        if name is not None and name in wanted:
                            found_here[name] = "".join(seq)
                        name = line[1:].split()[0]
                        seq = []
                    else:
                        seq.append(line.strip())
                if name is not None and name in wanted:
                    found_here[name] = "".join(seq)
            for gid in gids:
                s = found_here.get(gid)
                if s:
                    out.write(f">{gid}\n{s}\n")
                    n_found += 1
                else:
                    n_missing_seq += 1

    print(f"[fetch] found={n_found} missing_genome={n_missing_genome} missing_seq={n_missing_seq}")


if __name__ == "__main__":
    main()
