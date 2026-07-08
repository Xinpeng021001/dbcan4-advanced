#!/usr/bin/env python3
"""
Fold sequences to 3D structures with ESMFold (facebook/esmfold_v1 via transformers).
Shardable across GPUs (--shard/--nshards + CUDA_VISIBLE_DEVICES). Writes one PDB per
sequence to --outdir, plus a TSV of (id, length, plddt_mean, seconds, status).

Length cap (--max-len) skips proteins ESMFold can't hold in 24 GB; skipped ids recorded.
"""
import argparse, os, re, sys, time, json

_STD = set("ACDEFGHIKLMNPQRSTVWY")
def sanitize(seq):
    # map anything ESMFold's tokenizer can't handle (*, gaps, U/B/Z/O, lowercase) to X
    return re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', 'X', seq.upper())

def read_fasta(path):
    hid,buf=None,[]
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if hid is not None: yield hid,"".join(buf)
                hid=line[1:].rstrip("\n"); buf=[]
            else: buf.append(line.strip())
    if hid is not None: yield hid,"".join(buf)

def safe_id(h):  # header "414114|GH47" -> "414114"
    return h.split("|")[0].split()[0].replace("/","_")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--manifest", required=True, help="TSV log path (per shard)")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--max-len", type=int, default=1000)
    ap.add_argument("--min-len", type=int, default=20)
    args=ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding
    dev="cuda" if torch.cuda.is_available() else "cpu"
    tok=AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model=EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True).to(dev).eval()
    if dev=="cuda":
        model.esm=model.esm.half()
        try: model.trunk.set_chunk_size(64)   # lower peak memory for long seqs
        except Exception: pass

    seqs=[(safe_id(h),h,s) for h,s in read_fasta(args.fasta)]
    seqs=[x for i,x in enumerate(seqs) if i % args.nshards == args.shard]
    print(f"shard {args.shard}/{args.nshards}: {len(seqs)} sequences", file=sys.stderr, flush=True)

    with open(args.manifest,"w") as mf:
        mf.write("id\theader\tlength\tplddt_mean\tseconds\tstatus\n")
        for k,(sid,hdr,seq) in enumerate(seqs):
            L=len(seq)
            outp=os.path.join(args.outdir, f"{sid}.pdb")
            if os.path.exists(outp) and os.path.getsize(outp)>0:
                mf.write(f"{sid}\t{hdr}\t{L}\t\t\tcached\n"); continue
            if L<args.min_len or L>args.max_len:
                mf.write(f"{sid}\t{hdr}\t{L}\t\t\tskip_len\n"); continue
            try:
                ii=tok([sanitize(seq)], return_tensors="pt", add_special_tokens=False)["input_ids"].to(dev)
                t0=time.time()
                with torch.no_grad(): out=model(ii)
                if dev=="cuda": torch.cuda.synchronize()
                dt=time.time()-t0
                pdb=model.output_to_pdb(out)[0]
                open(outp,"w").write(pdb)
                pl=round(float(out["plddt"].mean()),4)
                mf.write(f"{sid}\t{hdr}\t{L}\t{pl}\t{dt:.2f}\tok\n"); mf.flush()
                if k%50==0: print(f"  {k}/{len(seqs)} {sid} L={L} {dt:.1f}s plddt={pl}", file=sys.stderr, flush=True)
            except torch.cuda.OutOfMemoryError:
                torch.cuda.empty_cache()
                mf.write(f"{sid}\t{hdr}\t{L}\t\t\toom\n"); mf.flush()
            except Exception as e:
                mf.write(f"{sid}\t{hdr}\t{L}\t\t\terr:{type(e).__name__}\n"); mf.flush()
    print(f"shard {args.shard} done", file=sys.stderr, flush=True)

if __name__=="__main__":
    main()
