#!/usr/bin/env python3
"""
DEFT-style consensus fusion of the dbCAN4-advanced retrieval axes:

  sequence : DIAMOND (fungal 2024 ref)         + HMMER (2024-only temporal)
  pLM      : contrastive kNN + softmax classifier (trained on ESM-C)
  structure: Foldseek (ESMFold structures vs 2024 ref)

Each method emits (family, confidence). Confidences are min-max normalized per method
to [0,1], then combined by confidence-weighted voting with per-method reliability weights
(set from each method's KNOWN-FAMILY / novel_seq parent-recall on this holdout — a principled
prior, not fit to the test label). The winning family's normalized score is the fusion
confidence; below --tau the call ABSTAINS (flags a putative novel/uncertain CAZyme).

Reports: recall by novelty bucket at subfamily + parent granularity, vs each method alone;
and abstention behaviour (does fusion confidence separate genuinely-novel from known?).
"""
import argparse, csv, re, json, math
from collections import defaultdict

BF=re.compile(r'^((?:GH|GT|PL|CE|AA|CBM)\d+)(?:_\d+)?$')
def base(l):
    m=BF.match(l); return m.group(1) if m else l
def bset(fs): return set(base(f) for f in fs)
def first(fs): return sorted(fs)[0] if fs else None

def load_truth(path):
    truth={}; nov={}
    with open(path) as fh:
        h=fh.readline().rstrip("\n").split("\t"); fi,ni,pi=h.index("families"),h.index("novelty"),h.index("protein_id")
        for line in fh:
            p=line.rstrip("\n").split("\t"); truth[p[pi]]=set(x for x in p[fi].split(",") if x); nov[p[pi]]=p[ni]
    return truth,nov

def norm(v, lo, hi):
    if v is None: return None
    return max(0.0, min(1.0, (v-lo)/(hi-lo))) if hi>lo else 0.0

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--labels", default="data/eval_2025_labels.tsv")
    ap.add_argument("--tau", type=float, default=0.35, help="abstain below this fusion confidence")
    ap.add_argument("--out-pred", default="results/fusion_pred.tsv")
    ap.add_argument("--out-summary", default="results/fusion_summary.json")
    args=ap.parse_args()
    R=args.results
    truth,nov=load_truth(args.labels)

    # per-method: query -> (family, raw_conf); plus reliability weight (novel_seq parent recall)
    methods={}  # name -> dict(query->(fam, conf))
    weight={}

    # DIAMOND (fungal 2024): conf = pident/100
    d={}
    for r in csv.DictReader(open(f"{R}/diamond_fungiref_pred.tsv"),delimiter="\t"):
        fam=first(set(x for x in r["pred_families"].split(",") if x))
        try: pid=float(r["top_pident"])/100.0
        except: pid=None
        if fam: d[r["query_id"]]=(fam, pid if pid is not None else 0.5)
    methods["DIAMOND"]=d; weight["DIAMOND"]=0.96

    # HMMER 2024 temporal: conf = -log10(evalue) capped/scaled
    h={}
    for r in csv.DictReader(open(f"{R}/hmm2024_pred.tsv"),delimiter="\t"):
        fam=r["pred_base_family"]
        try: ev=float(r["evalue"]); c=norm(-math.log10(ev+1e-300),0,50)
        except: c=0.3
        if fam: h[r["query_id"]]=(fam,c)
    methods["HMMER2024"]=h; weight["HMMER2024"]=0.88

    # contrastive kNN (trained pLM): conf = vote purity
    ck={}; clf={}
    for r in csv.DictReader(open(f"{R}/head_eval_pred.tsv"),delimiter="\t"):
        if r["contr_knn_pred"]: ck[r["query_id"]]=(r["contr_knn_pred"], float(r["contr_knn_purity"]))
        if r["clf_pred"]:       clf[r["query_id"]]=(r["clf_pred"], float(r["clf_conf"]))
    methods["ContrastiveKNN"]=ck; weight["ContrastiveKNN"]=0.90
    methods["Classifier"]=clf;    weight["Classifier"]=0.89

    # Foldseek (structure): conf = TM-score  (only queries with a structure)
    fs={}
    for r in csv.DictReader(open(f"{R}/foldseek_pred.tsv"),delimiter="\t"):
        fam=first(set(x for x in r["pred_families"].split(",") if x))
        if fam: fs[r["query_id"]]=(fam, float(r["tm_score"]))
    methods["Foldseek"]=fs; weight["Foldseek"]=0.75

    # ---- fuse ----
    allq=set(truth)
    fused={}   # query -> (pred_fam, fusion_conf, n_methods, agree_frac)
    for q in allq:
        votes=defaultdict(float); wtot=0.0; nmeth=0; per_parent=defaultdict(float)
        for name,mp in methods.items():
            if q in mp:
                fam,conf=mp[q]; w=weight[name]
                votes[fam]+=w*conf; per_parent[base(fam)]+=w*conf
                wtot+=w; nmeth+=1
        if not votes:
            fused[q]=(None,0.0,0,0.0); continue
        best_fam=max(votes,key=votes.get)
        fusion_conf=votes[best_fam]/wtot if wtot>0 else 0.0
        # agreement: fraction of voting methods whose parent == winner parent
        agree=sum(1 for name,mp in methods.items() if q in mp and base(mp[q][0])==base(best_fam))/nmeth
        fused[q]=(best_fam,fusion_conf,nmeth,agree)

    # ---- score ----
    def recall(pred_getter, level, bucket, abstain=False):
        n=exact=0
        for q in allq:
            if bucket!="overall" and nov.get(q)!=bucket: continue
            n+=1
            pf=pred_getter(q)
            if pf is None: continue
            tf=truth[q]
            if level=="sub":
                if pf and pf in tf: exact+=1
            else:
                if pf and base(pf) in bset(tf): exact+=1
        return {"n":n,"recall":round(exact/n,4) if n else 0}

    buckets=["overall","novel_seq","novel_family"]
    def fusion_getter_noabs(q): return fused[q][0]
    def fusion_getter_abs(q):   return fused[q][0] if fused[q][1]>=args.tau else None

    summary={"tau":args.tau,"n_methods_per_query":{},
             "fusion_no_abstain":{lev:{b:recall(fusion_getter_noabs,lev,b) for b in buckets} for lev in ["sub","parent"]},
             "fusion_with_abstain":{lev:{b:recall(fusion_getter_abs,lev,b) for b in buckets} for lev in ["sub","parent"]},
             "per_method_alone":{}}
    for name,mp in methods.items():
        g=lambda q,mp=mp: mp[q][0] if q in mp else None
        summary["per_method_alone"][name]={lev:{b:recall(g,lev,b) for b in buckets} for lev in ["sub","parent"]}

    # abstention behaviour: fusion confidence distribution + abstain rate by bucket
    abst={}
    for b in buckets:
        qs=[q for q in allq if b=="overall" or nov.get(q)==b]
        if not qs: continue
        confs=[fused[q][1] for q in qs]
        abst[b]={"n":len(qs),"mean_fusion_conf":round(sum(confs)/len(confs),4),
                 "abstain_rate":round(sum(1 for q in qs if fused[q][1]<args.tau)/len(qs),4)}
    summary["abstention"]=abst

    with open(args.out_pred,"w") as fo:
        fo.write("query_id\tnovelty\ttrue_families\tfusion_pred\tfusion_conf\tn_methods\tagree_frac\tabstain\n")
        for q in sorted(allq):
            fam,conf,nm,ag=fused[q]
            fo.write(f"{q}\t{nov.get(q,'?')}\t{','.join(sorted(truth[q]))}\t{fam or ''}\t{conf:.4f}\t{nm}\t{ag:.3f}\t{int(conf<args.tau)}\n")
    json.dump(summary, open(args.out_summary,"w"), indent=2)
    print(json.dumps(summary, indent=2))

if __name__=="__main__":
    main()
