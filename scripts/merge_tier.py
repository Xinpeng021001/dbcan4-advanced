#!/usr/bin/env python3
"""Merge tier/family from sample_ids.tsv into structure_evidence_scores.tsv,
producing the 9-column _final-style layout:
  protein_id, tier, family, foldseek_3di_best_hit, foldseek_3di_pident,
  foldseek_3di_evalue, foldseek_3di_bits_norm, saprot_cosine_to_cazyme_centroid,
  structure_evidence_score
"""
import sys
scores, sample_ids, out = sys.argv[1:4]
tf = {}
with open(sample_ids) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    pi, ti, fi = h.index("protein_id"), h.index("tier"), h.index("family")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) > fi:
            tf[p[pi]] = (p[ti], p[fi])
with open(scores) as fh, open(out, "w") as o:
    hdr = fh.readline().rstrip("\n").split("\t")   # protein_id, foldseek_3di_best_hit, ...
    o.write("\t".join([hdr[0], "tier", "family"] + hdr[1:]) + "\n")
    n = 0
    for line in fh:
        p = line.rstrip("\n").split("\t")
        t, f = tf.get(p[0], ("-", "-"))
        o.write("\t".join([p[0], t, f] + p[1:]) + "\n")
        n += 1
print(f"[merge] wrote {n} rows to {out}", file=sys.stderr)
