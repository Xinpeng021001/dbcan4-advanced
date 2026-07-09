#!/usr/bin/env python
"""
Extract SaProt (westlake-repl/SaProt_650M_AF2) structure-aware embeddings for a
directory of PDB structures, using foldseek to derive per-residue 3Di tokens
(via the SaProt repo's get_struc_seq helper) and combining with AA tokens into
the 441-token structure-aware vocabulary SaProt expects.

Also supports feeding externally-supplied lower-case 3Di strings (e.g. from
ProstT5 AA->3Di prediction) instead of running foldseek, via --di3_tsv
(columns: protein_id, aa_seq, di3_string).

Usage (from PDB structures):
    python saprot_embed.py --pdb_dir /path/to/pdbs --foldseek /usr/local/bin/foldseek \
        --model_dir westlake-repl/SaProt_650M_AF2 --out embeddings.npz --tsv_out summary.tsv

Usage (from ProstT5-predicted 3Di, no structures needed):
    python saprot_embed.py --di3_tsv prostt5_out.tsv --model_dir westlake-repl/SaProt_650M_AF2 \
        --out embeddings.npz --tsv_out summary.tsv
"""
import argparse
import glob
import os
import sys

import numpy as np
import torch
from transformers import EsmTokenizer, EsmForMaskedLM

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from foldseek_util import get_struc_seq


def combined_seqs_from_pdb_dir(pdb_dir, foldseek_bin, plddt_mask=False, limit=None):
    """Return dict protein_id -> (aa_seq, di3_seq, combined_seq) using foldseek."""
    pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")))
    if limit:
        pdb_files = pdb_files[:limit]
    out = {}
    for i, pf in enumerate(pdb_files):
        pid = os.path.splitext(os.path.basename(pf))[0]
        try:
            d = get_struc_seq(foldseek_bin, pf, plddt_mask=plddt_mask)
            # take chain 'A' if present else first chain
            key = "A" if "A" in d else next(iter(d))
            aa_seq, di3_seq, combined_seq = d[key]
            out[pid] = (aa_seq, di3_seq, combined_seq)
        except Exception as e:
            print(f"[saprot] WARN: failed on {pid}: {e}", file=sys.stderr)
        if (i + 1) % 20 == 0:
            print(f"[saprot] foldseek 3Di extraction {i+1}/{len(pdb_files)}", file=sys.stderr)
    return out


def combined_seqs_from_tsv(tsv_path):
    """di3_tsv columns: protein_id, aa_len, di3_string, di3_len[, ...]. We need aa_seq
    too -- if not present in the tsv we cannot build the combined seq (need AA chars).
    Expects an optional 'aa_seq' column; if absent, raises."""
    import csv
    out = {}
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "aa_seq" not in reader.fieldnames:
            raise ValueError("di3_tsv must contain an 'aa_seq' column to build SaProt combined tokens")
        for row in reader:
            pid = row["protein_id"]
            aa_seq = row["aa_seq"]
            di3_seq = row["di3_string"]
            combined = "".join(a + b.lower() for a, b in zip(aa_seq, di3_seq))
            out[pid] = (aa_seq, di3_seq, combined)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb_dir", default=None)
    ap.add_argument("--foldseek", default="/usr/local/bin/foldseek")
    ap.add_argument("--di3_tsv", default=None)
    ap.add_argument("--model_dir", default="westlake-repl/SaProt_650M_AF2")
    ap.add_argument("--out", required=True, help="npz file: ids, embeddings (N x D)")
    ap.add_argument("--tsv_out", required=True)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--max_len", type=int, default=1024)
    args = ap.parse_args()

    if args.pdb_dir:
        seqs = combined_seqs_from_pdb_dir(args.pdb_dir, args.foldseek, limit=args.limit)
    elif args.di3_tsv:
        seqs = combined_seqs_from_tsv(args.di3_tsv)
    else:
        raise SystemExit("Must supply --pdb_dir or --di3_tsv")

    print(f"[saprot] {len(seqs)} combined AA+3Di sequences ready", file=sys.stderr)

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    tokenizer = EsmTokenizer.from_pretrained(args.model_dir)
    model = EsmForMaskedLM.from_pretrained(args.model_dir)
    model = model.to(device)
    model.eval()
    print(f"[saprot] model loaded on {device}", file=sys.stderr)

    ids = []
    embs = []
    rows = []
    with torch.no_grad():
        for i, (pid, (aa_seq, di3_seq, combined)) in enumerate(seqs.items()):
            combined_trunc = combined[: args.max_len * 2]  # each residue -> 2 chars
            inputs = tokenizer(combined_trunc, return_tensors="pt").to(device)
            out = model(**inputs, output_hidden_states=True)
            last_hidden = out.hidden_states[-1][0]  # (L, D)
            emb = last_hidden.mean(dim=0).float().cpu().numpy()
            ids.append(pid)
            embs.append(emb)
            rows.append((pid, len(aa_seq), len(di3_seq)))
            if (i + 1) % 10 == 0:
                print(f"[saprot] embedded {i+1}/{len(seqs)}", file=sys.stderr)

    embs = np.vstack(embs) if embs else np.zeros((0, 0))
    np.savez(args.out, ids=np.array(ids), embeddings=embs)
    with open(args.tsv_out, "w") as fh:
        fh.write("protein_id\taa_len\tdi3_len\n")
        for pid, alen, dlen in rows:
            fh.write(f"{pid}\t{alen}\t{dlen}\n")
    print(f"[saprot] wrote {len(ids)} embeddings (dim={embs.shape[1] if embs.size else 0}) to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
