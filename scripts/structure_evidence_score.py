#!/usr/bin/env python
"""
Compute a structure-similarity-based CAZyme-likelihood score for a set of
query proteins, using ProstT5-predicted 3Di strings (AA->3Di, no folding
needed) and comparing against the CAZyme3D_id50 reference set two ways:

  (A) foldseek 3Di-sequence alignment (fast, alignment-based; treats the
      predicted 3Di string as a query sequence against a 3Di sequence DB
      built from real CAZyme3D_id50 structures) -> best-hit bit score /
      identity / e-value per query.

  (B) SaProt (westlake-repl/SaProt_650M_AF2) mean-pooled embedding cosine
      similarity to the centroid embedding of a reference set of known
      CAZyme3D_id50 structures (embeddings computed once, reused) -> a
      continuous embedding-similarity score per query.

Combines into a single "structure_evidence_score" per protein
(0-1, higher = more CAZyme-like by structure).

Usage:
    python structure_evidence_score.py \
        --prostt5_tsv all_prostt5.tsv \
        --cazyme3d_3di_tsv mapping_foldseek_3di.tsv \
        --out_dir /path/to/output_dir \
        --device cuda:0
    # writes <out_dir>/structure_evidence_scores.tsv (add --skip_saprot to
    # run only the foldseek-3Di-search component, no GPU needed)
"""
import argparse
import csv
import os
import subprocess
import sys

import numpy as np
import torch
from transformers import EsmTokenizer, EsmForMaskedLM


def read_prostt5_tsv(path):
    """protein_id, aa_seq, aa_len, di3_string, di3_len, gen_seconds"""
    d = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            d[row["protein_id"]] = (row["aa_seq"], row["di3_string"])
    return d


def build_cazyme3d_3di_fasta(cazyme3d_3di_tsv, out_fasta):
    """Parse the raw foldseek structureto3didescriptor output (desc, aa, 3di, ...)
    and write a lower-case 3Di fasta for foldseek sequence-search."""
    n = 0
    with open(cazyme3d_3di_tsv) as fh, open(out_fasta, "w") as out:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            desc, aa, di3 = parts[0], parts[1], parts[2]
            acc = desc.split(" ")[0]
            if acc.endswith(".pdb_A"):
                acc = acc[: -len(".pdb_A")]
            out.write(f">{acc}\n{di3.lower()}\n")
            n += 1
    print(f"[score] wrote {n} CAZyme3D_id50 3Di sequences to {out_fasta}", file=sys.stderr)
    return n


def foldseek_3di_search(query_3di_fasta, ref_3di_fasta, foldseek_bin, out_dir):
    """Compare ProstT5-predicted 3Di strings (no real 3D coordinates available,
    so foldseek's structure-based createdb/TM-align pipeline cannot be used --
    that requires actual PDB coordinates) against the CAZyme3D_id50 3Di
    alphabet strings (derived from real structures via foldseek
    structureto3didescriptor). We treat the 3Di alphabet (~20 letters, same
    cardinality as amino acids) as a generic sequence alphabet and align with
    mmseqs2 easy-search -- an established proxy used in the ProstT5/SaProt
    community for "structure-free" structural homology search. This is an
    approximation: mmseqs2's default substitution matrix is BLOSUM (tuned for
    amino acids, not 3Di), so identity/coverage are more informative here than
    the bit score scale."""
    os.makedirs(out_dir, exist_ok=True)
    tmp = os.path.join(out_dir, "tmp")
    m8 = os.path.join(out_dir, "hits.m8")
    cmd = (
        f"mmseqs easy-search {query_3di_fasta} {ref_3di_fasta} {m8} {tmp} "
        f"--threads 32 -s 7.5 --max-seqs 5 --alignment-mode 3 "
        f"--format-output query,target,pident,alnlen,evalue,bits"
    )
    subprocess.run(cmd, shell=True, check=True)
    return m8


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prostt5_tsv", required=True)
    ap.add_argument("--cazyme3d_3di_tsv", required=True,
                    help="raw foldseek structureto3didescriptor tsv (desc, aa, 3di, ...)")
    ap.add_argument("--foldseek", default="/usr/local/bin/foldseek")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--saprot_model", default="westlake-repl/SaProt_650M_AF2")
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--skip_saprot", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    prostt5 = read_prostt5_tsv(args.prostt5_tsv)
    print(f"[score] {len(prostt5)} query proteins with ProstT5 3Di", file=sys.stderr)

    # --- (A) foldseek 3Di sequence alignment ---
    query_3di_fasta = os.path.join(args.out_dir, "query_3di.fasta")
    with open(query_3di_fasta, "w") as out:
        for pid, (aa, di3) in prostt5.items():
            out.write(f">{pid}\n{di3}\n")

    ref_3di_fasta = os.path.join(args.out_dir, "ref_3di.fasta")
    build_cazyme3d_3di_fasta(args.cazyme3d_3di_tsv, ref_3di_fasta)

    m8 = foldseek_3di_search(query_3di_fasta, ref_3di_fasta, args.foldseek,
                              os.path.join(args.out_dir, "foldseek_3di_search"))

    best_3di_hit = {}
    with open(m8) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            q, t, pident, alnlen, evalue, bits = p[0], p[1], float(p[2]), int(p[3]), float(p[4]), float(p[5])
            if q not in best_3di_hit or bits > best_3di_hit[q][3]:
                best_3di_hit[q] = (t, pident, evalue, bits)
    print(f"[score] foldseek 3Di search: {len(best_3di_hit)}/{len(prostt5)} queries with a hit", file=sys.stderr)

    # normalize bit score to 0-1 via simple saturating transform
    max_bits = max((v[3] for v in best_3di_hit.values()), default=1.0)

    # --- (B) SaProt embedding similarity (optional, GPU) ---
    saprot_sim = {}
    if not args.skip_saprot:
        device = torch.device(args.device if torch.cuda.is_available() else "cpu")
        tokenizer = EsmTokenizer.from_pretrained(args.saprot_model)
        model = EsmForMaskedLM.from_pretrained(args.saprot_model).to(device)
        model.eval()
        print("[score] SaProt model loaded", file=sys.stderr)

        # embed CAZyme3D_id50 reference structures for the centroid: reservoir-sample
        # 500 structures UNIFORMLY across the full reference file (NOT a head-slice --
        # a fixed-size prefix would bias the centroid toward whatever ordering the
        # upstream foldseek/tar extraction happened to produce).
        import random
        random.seed(1)
        n_ref_sample = 500
        reservoir = []
        n_seen = 0
        with open(args.cazyme3d_3di_tsv) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                desc, aa, di3 = parts[0], parts[1], parts[2]
                acc = desc.split(" ")[0]
                item = (acc, aa, di3.lower())
                n_seen += 1
                if len(reservoir) < n_ref_sample:
                    reservoir.append(item)
                else:
                    j = random.randint(0, n_seen - 1)
                    if j < n_ref_sample:
                        reservoir[j] = item
        ref_sample = {acc: (aa, di3) for acc, aa, di3 in reservoir}
        print(f"[score] reservoir-sampled {len(ref_sample)} of {n_seen} CAZyme3D_id50 "
              f"structures for the SaProt reference centroid", file=sys.stderr)

        def embed_combined(aa, di3, max_len=1024):
            combined = "".join(a + b for a, b in zip(aa[:max_len], di3[:max_len]))
            inputs = tokenizer(combined, return_tensors="pt", truncation=True, max_length=max_len * 2 + 2).to(device)
            with torch.no_grad():
                out = model(**inputs, output_hidden_states=True)
            return out.hidden_states[-1][0].mean(dim=0).float().cpu().numpy()

        ref_embs = []
        for i, (acc, (aa, di3)) in enumerate(ref_sample.items()):
            ref_embs.append(embed_combined(aa, di3))
            if (i + 1) % 100 == 0:
                print(f"[score] ref embed {i+1}/{len(ref_sample)}", file=sys.stderr)
        ref_centroid = np.mean(np.vstack(ref_embs), axis=0)
        ref_centroid /= np.linalg.norm(ref_centroid)

        for i, (pid, (aa, di3)) in enumerate(prostt5.items()):
            try:
                emb = embed_combined(aa, di3)
                emb_n = emb / (np.linalg.norm(emb) + 1e-9)
                cos = float(np.dot(emb_n, ref_centroid))
                saprot_sim[pid] = cos
            except Exception as e:
                print(f"[score] WARN saprot embed failed {pid}: {e}", file=sys.stderr)
            if (i + 1) % 100 == 0:
                print(f"[score] query embed {i+1}/{len(prostt5)}", file=sys.stderr)

    # --- combine ---
    with open(os.path.join(args.out_dir, "structure_evidence_scores.tsv"), "w") as out:
        out.write("protein_id\tfoldseek_3di_best_hit\tfoldseek_3di_pident\tfoldseek_3di_evalue\t"
                   "foldseek_3di_bits_norm\tsaprot_cosine_to_cazyme_centroid\tstructure_evidence_score\n")
        for pid in prostt5:
            hit = best_3di_hit.get(pid)
            if hit:
                t, pident, evalue, bits = hit
                bits_norm = bits / max_bits if max_bits else 0.0
            else:
                t, pident, evalue, bits_norm = "-", 0.0, float("nan"), 0.0
            cos = saprot_sim.get(pid, float("nan"))
            # combine: average of foldseek-bits-norm and (cos+1)/2 rescale, ignoring missing
            parts = [bits_norm]
            if not np.isnan(cos):
                parts.append((cos + 1) / 2)
            score = float(np.mean(parts)) if parts else 0.0
            out.write(f"{pid}\t{t}\t{pident}\t{evalue}\t{bits_norm:.4f}\t{cos:.4f}\t{score:.4f}\n")

    print(f"[score] wrote structure_evidence_scores.tsv to {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
