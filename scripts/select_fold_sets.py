#!/usr/bin/env python3
"""
Select the sequences to fold for the structure tier (self-contained temporal design):

  EVAL side (from eval_2025.faa):
    - ALL truly-novel-base-family seqs (6 fams, no 2024 parent) -> the key test
    - ALL novel-subfamily-of-known-parent seqs
    - a capped random sample of novel-seq known-family seqs (controls)
  REFERENCE side (from reference_2024.faa):
    - up to --reps-per-family sequences per BASE family (Foldseek target DB)

Writes fold_eval.faa and fold_ref.faa; prints counts.
"""
import argparse, re, random
from collections import defaultdict

BF=re.compile(r'^((?:GH|GT|PL|CE|AA|CBM)\d+)(?:_\d+)?$')
def base(l):
    m=BF.match(l); return m.group(1) if m else None

def read_fasta(path):
    hid,buf=None,[]
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if hid is not None: yield hid,"".join(buf)
                hid=line[1:].rstrip("\n"); buf=[]
            else: buf.append(line.strip())
    if hid is not None: yield hid,"".join(buf)

def fams_of(header):
    if "|" not in header: return set()
    return set(f for f in re.split(r'[|,]', header.split("|",1)[1]) if BF.match(f))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--eval", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out-eval", required=True)
    ap.add_argument("--out-ref", required=True)
    ap.add_argument("--reps-per-family", type=int, default=5)
    ap.add_argument("--control-sample", type=int, default=500)
    ap.add_argument("--max-len", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    args=ap.parse_args()
    random.seed(args.seed)

    TRUE_NOVEL_BASE={"CBM104","CBM3","CBM8","GT109","GT119","PL29"}
    # novelty per protein
    nov={}; truth={}
    with open(args.labels) as fh:
        h=fh.readline().rstrip("\n").split("\t"); fi,ni,pi=h.index("families"),h.index("novelty"),h.index("protein_id")
        for line in fh:
            p=line.rstrip("\n").split("\t"); nov[p[pi]]=p[ni]; truth[p[pi]]=set(x for x in p[fi].split(",") if x)

    def pid_of(h): return h.split("|")[0].split()[0]
    eval_seqs=[(pid_of(h),h,s) for h,s in read_fasta(args.eval) if args.max_len>=len(s)>=20]
    truenovel=[]; subnovel=[]; controls=[]
    for pid,h,s in eval_seqs:
        nv=nov.get(pid); tb=set(base(f) for f in truth.get(pid,set()) if base(f))
        if nv=="novel_family" and (tb & TRUE_NOVEL_BASE): truenovel.append((pid,h,s))
        elif nv=="novel_family": subnovel.append((pid,h,s))
        elif nv=="novel_seq": controls.append((pid,h,s))
    random.shuffle(controls); controls=controls[:args.control_sample]
    chosen_eval = truenovel + subnovel + controls
    with open(args.out_eval,"w") as fo:
        for pid,h,s in chosen_eval: fo.write(f">{h}\n{s}\n")
    print(f"EVAL fold set: {len(chosen_eval)} "
          f"(truly-novel-base={len(truenovel)}, novel-subfam={len(subnovel)}, controls={len(controls)})")

    # reference: reps per base family
    fam_members=defaultdict(list)
    for h,s in read_fasta(args.ref):
        if not (args.max_len>=len(s)>=20): continue
        for bf in set(base(f) for f in fams_of(h) if base(f)):
            fam_members[bf].append((pid_of(h),h,s))
    chosen_ref=[]; seen=set()
    for bf,mem in fam_members.items():
        random.shuffle(mem)
        for pid,h,s in mem[:args.reps_per_family]:
            if pid in seen: continue
            seen.add(pid); chosen_ref.append((pid,h,s))
    with open(args.out_ref,"w") as fo:
        for pid,h,s in chosen_ref: fo.write(f">{h}\n{s}\n")
    print(f"REF fold set: {len(chosen_ref)} across {len(fam_members)} base families "
          f"(<= {args.reps_per_family}/family)")

if __name__=="__main__":
    main()
