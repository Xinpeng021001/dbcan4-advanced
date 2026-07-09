#!/usr/bin/env python
"""
Validate ProstT5 (Rostlab/ProstT5) AA -> 3Di translation on a small test set.

Usage:
    python prostt5_validate.py --fasta test.faa --out out.tsv [--device cuda:0]

Produces a TSV: protein_id, aa_len, di3_string, di3_len (lower-cased 3Di, ready
for SaProt / foldseek-style downstream consumption).
"""
import argparse
import sys
import time

import torch
from transformers import T5Tokenizer, T5ForConditionalGeneration


def read_fasta(path):
    recs = {}
    name = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        recs[name] = "".join(seq)
    return recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--max_len", type=int, default=1024, help="truncate AA sequences longer than this")
    ap.add_argument("--model_name", default="Rostlab/ProstT5")
    args = ap.parse_args()

    recs = read_fasta(args.fasta)
    print(f"[prostt5] loaded {len(recs)} sequences from {args.fasta}", file=sys.stderr)

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    print(f"[prostt5] loading tokenizer/model on {device} ...", file=sys.stderr)
    t0 = time.time()
    tokenizer = T5Tokenizer.from_pretrained(args.model_name, do_lower_case=False)
    model = T5ForConditionalGeneration.from_pretrained(args.model_name)
    model = model.to(device)
    if device.type == "cuda":
        model = model.half()
    model.eval()
    print(f"[prostt5] model loaded in {time.time()-t0:.1f}s", file=sys.stderr)

    results = []
    with torch.no_grad():
        for i, (pid, seq) in enumerate(recs.items()):
            seq = seq[: args.max_len]
            # ProstT5 expects spaced residues, uppercase, non-standard -> X, prefixed with <AA2fold>
            seq_spaced = " ".join(list(seq.upper()))
            seq_spaced = seq_spaced.replace("U", "X").replace("Z", "X").replace("O", "X").replace("B", "X")
            input_str = "<AA2fold> " + seq_spaced
            ids = tokenizer(input_str, add_special_tokens=True, return_tensors="pt").to(device)
            t1 = time.time()
            gen = model.generate(
                ids.input_ids,
                attention_mask=ids.attention_mask,
                max_length=len(seq) + 5,
                min_length=len(seq) + 1,
                do_sample=False,
                num_beams=1,
            )
            dt = time.time() - t1
            di3 = tokenizer.decode(gen[0], skip_special_tokens=True).replace(" ", "")
            di3_lower = di3.lower()
            results.append((pid, seq, len(seq), di3_lower, len(di3_lower), dt))
            print(f"[prostt5] {i+1}/{len(recs)} {pid} len={len(seq)} 3di_len={len(di3_lower)} t={dt:.2f}s", file=sys.stderr)

    with open(args.out, "w") as fh:
        fh.write("protein_id\taa_seq\taa_len\tdi3_string\tdi3_len\tgen_seconds\n")
        for pid, seq, alen, di3, dlen, dt in results:
            fh.write(f"{pid}\t{seq}\t{alen}\t{di3}\t{dlen}\t{dt:.3f}\n")
    print(f"[prostt5] wrote {len(results)} rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
