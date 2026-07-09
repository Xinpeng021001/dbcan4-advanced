#!/usr/bin/env python3
"""
Embed protein sequences with ESM-C (EvolutionaryScale) and mean-pool per sequence.

One process per GPU (set via --shard / --nshards + CUDA_VISIBLE_DEVICES in the wrapper).
Reads a FASTA whose headers are "ID|FAM[|FAM...]", processes the shard's slice,
writes <out_prefix>.shardK.npz with:
    ids       (N,)   protein IDs (str)
    fams      (N,)   comma-joined family labels (str)
    emb       (N,D)  float16 mean-pooled residue embeddings (D=1152 for esmc_600m)

Mean-pool excludes BOS/EOS. Sequences longer than --max-len are truncated.
"""
import argparse, sys, time
import numpy as np

def read_fasta(path):
    hid, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if hid is not None:
                    yield hid, "".join(buf)
                hid = line[1:].strip()   # strip ALL trailing whitespace (MMseqs rep_seq headers carry a trailing space)
                buf = []
            else:
                buf.append(line.strip())
        if hid is not None:
            yield hid, "".join(buf)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--model", default="esmc_600m")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--max-len", type=int, default=1500)
    ap.add_argument("--log-every", type=int, default=500)
    args = ap.parse_args()

    import torch
    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein, LogitsConfig

    # gather this shard's records
    recs = []
    for i, (hid, seq) in enumerate(read_fasta(args.fasta)):
        if i % args.nshards != args.shard:
            continue
        pid = hid.split("|")[0]
        fams = "|".join(hid.split("|")[1:]) if "|" in hid else ""
        recs.append((pid, fams, seq[:args.max_len]))
    print(f"[shard {args.shard}/{args.nshards}] {len(recs)} sequences", file=sys.stderr, flush=True)

    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model = ESMC.from_pretrained(args.model).to(dev).eval()

    ids, fams, embs = [], [], []
    t0 = time.time()
    with torch.no_grad():
        for k, (pid, fam, seq) in enumerate(recs):
            if not seq:
                continue
            p = ESMProtein(sequence=seq)
            enc = model.encode(p)
            out = model.logits(enc, LogitsConfig(sequence=True, return_embeddings=True))
            e = out.embeddings[0]                    # (L+2, D)
            mp = e[1:-1].mean(0).float().cpu().numpy().astype(np.float16)
            ids.append(pid); fams.append(fam); embs.append(mp)
            if (k+1) % args.log_every == 0:
                rate = (k+1)/(time.time()-t0)
                print(f"[shard {args.shard}] {k+1}/{len(recs)}  {rate:.1f} seq/s", file=sys.stderr, flush=True)

    emb = np.stack(embs).astype(np.float16)
    np.savez_compressed(f"{args.out_prefix}.shard{args.shard}.npz",
                        ids=np.array(ids), fams=np.array(fams), emb=emb)
    print(f"[shard {args.shard}] wrote {emb.shape} in {time.time()-t0:.0f}s", file=sys.stderr, flush=True)

if __name__ == "__main__":
    main()
