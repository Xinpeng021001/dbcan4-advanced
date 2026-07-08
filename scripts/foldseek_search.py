#!/usr/bin/env python3
"""
Structure tier: Foldseek search of eval structures vs a 2024-reference structure DB.

Given ESMFold PDBs for eval and reference (built at base-family or subfamily granularity),
this:
  1. foldseek createdb on the reference PDB dir
  2. foldseek easy-search eval PDB dir vs reference DB (TM-align mode, --alignment-type 1)
  3. best structural hit per query -> predicted family (from the reference PDB's header/label map)
  4. scores exact-subfamily + parent-family recall by novelty bucket, and reports structural
     novelty signals (top TM-score, structural vs sequence agreement).

Reference/eval family labels come from a TSV mapping structure-id -> families (built from the
fold-set FASTA headers).
"""
import argparse, csv, json, os, re, subprocess, sys
from collections import defaultdict, Counter

BF=re.compile(r'^((?:GH|GT|PL|CE|AA|CBM)\d+)(?:_\d+)?$')
def base(l):
    m=BF.match(l); return m.group(1) if m else l
def base_set(fs): return set(base(f) for f in fs)

def load_label_map(fasta):
    """structure-id (header field 0) -> set(families) from a fold-set FASTA."""
    m={}
    for line in open(fasta):
        if not line.startswith(">"): continue
        h=line[1:].rstrip("\n"); sid=h.split("|")[0].split()[0]
        fams=set(f for f in re.split(r'[|,]', h.split("|",1)[1]) if BF.match(f)) if "|" in h else set()
        m[sid]=fams
    return m

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--eval-pdb", required=True)
    ap.add_argument("--ref-pdb", required=True)
    ap.add_argument("--eval-fasta", required=True)
    ap.add_argument("--ref-fasta", required=True)
    ap.add_argument("--labels", required=True, help="eval_2025_labels.tsv for novelty buckets")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out-pred", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--foldseek", default="foldseek")
    args=ap.parse_args()
    os.makedirs(args.workdir, exist_ok=True)

    ref_lab=load_label_map(args.ref_fasta)
    eval_lab=load_label_map(args.eval_fasta)

    # truth + novelty
    truth={}; nov={}
    with open(args.labels) as fh:
        h=fh.readline().rstrip("\n").split("\t"); fi,ni,pi=h.index("families"),h.index("novelty"),h.index("protein_id")
        for line in fh:
            p=line.rstrip("\n").split("\t"); truth[p[pi]]=set(x for x in p[fi].split(",") if x); nov[p[pi]]=p[ni]

    refdb=os.path.join(args.workdir,"refdb")
    aln=os.path.join(args.workdir,"aln.m8")
    tmp=os.path.join(args.workdir,"tmp")
    # 1. easy-search does createdb internally; TM-align mode (alignment-type 1) for structural sim
    cmd=[args.foldseek,"easy-search",args.eval_pdb,args.ref_pdb,aln,tmp,
         "--alignment-type","1","-e","10","--max-seqs","50","--threads",str(args.threads),
         "--format-output","query,target,alntmscore,evalue,bits,lddt"]
    print("running:"," ".join(cmd), file=sys.stderr, flush=True)
    subprocess.run(cmd, check=True)

    # 2. parse best hit per query by TM-score
    best=defaultdict(lambda:(None,-1.0,None))  # query -> (target, tmscore, evalue)
    with open(aln) as fh:
        for line in fh:
            p=line.rstrip("\n").split("\t")
            if len(p)<4: continue
            q=p[0].replace(".pdb",""); t=p[1].replace(".pdb",""); tm=float(p[2]); ev=float(p[3])
            q=q.split("|")[0]; t=t.split("|")[0]
            if tm>best[q][1]: best[q]=(t,tm,ev)

    # 3. predicted family = families of best structural target
    rows=[]
    for q,(t,tm,ev) in best.items():
        pf=ref_lab.get(t,set()) if t else set()
        rows.append((q,t,tm,ev,pf))

    # 4. score
    def scored(level,bucket):
        n=exact=0
        for q,(t,tm,ev) in best.items():
            if q not in truth: continue
            if bucket!="overall" and nov.get(q)!=bucket: continue
            n+=1; pf=ref_lab.get(t,set()); tf=truth[q]
            if level=="sub" and pf and pf==tf: exact+=1
            if level=="parent" and base_set(pf) and base_set(pf)==base_set(tf): exact+=1
        return {"n":n,"exact":round(exact/n,4) if n else 0}

    buckets=["overall","novel_seq","novel_family"]
    summary={"n_eval_with_hit":len(best),
             "subfamily":{b:scored("sub",b) for b in buckets},
             "parent_family":{b:scored("parent",b) for b in buckets}}
    with open(args.out_pred,"w") as fo:
        fo.write("query_id\tnovelty\ttrue_families\tbest_target\ttm_score\tevalue\tpred_families\n")
        for q,(t,tm,ev) in sorted(best.items()):
            fo.write(f"{q}\t{nov.get(q,'?')}\t{','.join(sorted(truth.get(q,set())))}\t{t}\t{tm:.4f}\t{ev:.2e}\t{','.join(sorted(ref_lab.get(t,set())))}\n")
    json.dump(summary, open(args.out_summary,"w"), indent=2)
    print(json.dumps(summary, indent=2))

if __name__=="__main__":
    main()
