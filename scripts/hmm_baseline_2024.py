#!/usr/bin/env python3
"""
Temporally-clean HMMER baseline: build per-base-family profile HMMs from the 2024 reference ONLY,
then hmmscan the 2025 eval set. This is the fair HMM counterpart to the DIAMOND-2024 baseline —
it answers "what does an HMM tier that has only seen 2024 do on 2025 families?" (as opposed to the
production dbCAN.hmm, whose current DB already contains the 2025 families).

Pipeline (all local tools on met: mafft, hmmbuild, hmmpress, hmmscan):
  1. group reference_2024 sequences by BASE family (GH2_13 -> GH2), matching dbCAN.hmm granularity
  2. per family with >= --min-seqs: sample up to --max-seqs, mafft --auto, hmmbuild
  3. cat -> hmmpress -> hmmscan eval_2025
  4. best (lowest full-seq e-value) profile per query -> predicted base family
Outputs a per-query TSV (query_id, pred_base_family, evalue) for downstream scoring.
"""
import argparse, os, re, subprocess, sys, random
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
    tail=header.split("|",1)[1]
    return set(f for f in re.split(r'[|,]', tail) if BF.match(f))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--eval", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out-pred", required=True)
    ap.add_argument("--min-seqs", type=int, default=5)
    ap.add_argument("--max-seqs", type=int, default=300)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--seed", type=int, default=0)
    args=ap.parse_args()
    random.seed(args.seed)
    os.makedirs(args.workdir, exist_ok=True)
    msa_dir=os.path.join(args.workdir,"msa"); hmm_dir=os.path.join(args.workdir,"hmm")
    os.makedirs(msa_dir, exist_ok=True); os.makedirs(hmm_dir, exist_ok=True)

    # 1. group by base family
    fam_seqs=defaultdict(list)
    for hid,seq in read_fasta(args.ref):
        if not seq: continue
        for bf in set(base(f) for f in fams_of(hid) if base(f)):
            fam_seqs[bf].append((hid.split("|")[0], seq))
    fams=[f for f,s in fam_seqs.items() if len(s)>=args.min_seqs]
    print(f"{len(fams)} families with >= {args.min_seqs} seqs (of {len(fam_seqs)} total)", file=sys.stderr, flush=True)

    # 2. write per-family fasta (sampled)
    for bf in fams:
        seqs=fam_seqs[bf]
        if len(seqs)>args.max_seqs: seqs=random.sample(seqs, args.max_seqs)
        with open(os.path.join(msa_dir,f"{bf}.faa"),"w") as fo:
            for i,(pid,s) in enumerate(seqs): fo.write(f">{pid}_{i}\n{s}\n")

    # 2b. mafft + hmmbuild per family, parallel via xargs
    build_sh=os.path.join(args.workdir,"build_one.sh")
    with open(build_sh,"w") as fo:
        fo.write(f"""#!/bin/bash
bf="$1"; MSA={msa_dir}; HMM={hmm_dir}
if [ ! -s "$HMM/$bf.hmm" ]; then
  mafft --auto --thread 1 --quiet "$MSA/$bf.faa" > "$MSA/$bf.aln" 2>/dev/null
  hmmbuild --amino -n "$bf" "$HMM/$bf.hmm" "$MSA/$bf.aln" > /dev/null 2>&1
fi
""")
    os.chmod(build_sh,0o755)
    fam_list=os.path.join(args.workdir,"fams.txt"); open(fam_list,"w").write("\n".join(fams)+"\n")
    print("building HMMs (mafft+hmmbuild, parallel)...", file=sys.stderr, flush=True)
    subprocess.run(f"cat {fam_list} | xargs -P {args.threads} -I{{}} {build_sh} {{}}",
                   shell=True, check=True)

    # 3. cat + press
    allhmm=os.path.join(args.workdir,"ref2024_fam.hmm")
    subprocess.run(f"cat {hmm_dir}/*.hmm > {allhmm}", shell=True, check=True)
    n_built=int(subprocess.run(f"grep -c '^NAME' {allhmm}", shell=True, capture_output=True, text=True).stdout.strip() or 0)
    print(f"built {n_built} family HMMs", file=sys.stderr, flush=True)
    for ext in (".h3f",".h3i",".h3m",".h3p"):
        p=allhmm+ext
        if os.path.exists(p): os.remove(p)
    subprocess.run(f"hmmpress {allhmm}", shell=True, check=True)

    # 4. hmmscan eval
    domtbl=os.path.join(args.workdir,"eval2025.domtbl")
    subprocess.run(f"hmmscan --cpu {args.threads} --domtblout {domtbl} -E 1e-3 {allhmm} {args.eval} > /dev/null 2>&1",
                   shell=True, check=True)

    # parse best full-seq e-value per query
    best={}
    with open(domtbl) as fh:
        for line in fh:
            if line.startswith("#"): continue
            p=line.split()
            if len(p)<13: continue
            fam=p[0]; qid=p[3].split("|")[0]; ev=float(p[6])
            if qid not in best or ev<best[qid][1]: best[qid]=(fam,ev)
    with open(args.out_pred,"w") as fo:
        fo.write("query_id\tpred_base_family\tevalue\n")
        for qid,(fam,ev) in best.items():
            fo.write(f"{qid}\t{fam}\t{ev:.2e}\n")
    print(f"scanned; {len(best)} queries with a hit -> {args.out_pred}", file=sys.stderr, flush=True)

if __name__=="__main__":
    main()
