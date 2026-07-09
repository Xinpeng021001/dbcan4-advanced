#!/usr/bin/env python
"""Extract a stratified validation sample from Track A's tiered_proteins.tsv.gz
for structure-evidence scoring: gray_zone (target), plus high_confidence_cazyme
and high_confidence_non_cazyme (positive/negative controls)."""
import argparse
import gzip
import random


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiered_tsv_gz", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n_gray", type=int, default=1500)
    ap.add_argument("--n_pos", type=int, default=500)
    ap.add_argument("--n_neg", type=int, default=500)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    random.seed(args.seed)
    gray, hc_cazyme, hc_noncazyme = [], [], []

    with gzip.open(args.tiered_tsv_gz, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        tier_idx = header.index("tier")
        pid_idx = header.index("protein_id")
        fam_idx = header.index("recommend_family") if "recommend_family" in header else None
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= tier_idx:
                continue
            tier = parts[tier_idx]
            pid = parts[pid_idx]
            fam = parts[fam_idx] if fam_idx is not None and fam_idx < len(parts) else "-"
            if tier == "gray_zone":
                gray.append((pid, fam))
            elif tier == "high_confidence_cazyme":
                hc_cazyme.append((pid, fam))
            elif tier == "high_confidence_non_cazyme":
                hc_noncazyme.append((pid, fam))

    print(f"totals: gray={len(gray)} hc_cazyme={len(hc_cazyme)} hc_noncazyme={len(hc_noncazyme)}")
    g_s = random.sample(gray, min(args.n_gray, len(gray)))
    c_s = random.sample(hc_cazyme, min(args.n_pos, len(hc_cazyme)))
    n_s = random.sample(hc_noncazyme, min(args.n_neg, len(hc_noncazyme)))

    with open(args.out, "w") as out:
        out.write("protein_id\ttier\tfamily\n")
        for p, fam in g_s:
            out.write(f"{p}\tgray_zone\t{fam}\n")
        for p, fam in c_s:
            out.write(f"{p}\thigh_confidence_cazyme\t{fam}\n")
        for p, fam in n_s:
            out.write(f"{p}\thigh_confidence_non_cazyme\t{fam}\n")
    print(f"wrote {len(g_s)+len(c_s)+len(n_s)} rows to {args.out}")


if __name__ == "__main__":
    main()
