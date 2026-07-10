#!/usr/bin/env python3
"""Parse hmmsearch domtblout -> per-domain sequences (envelope coords), build
independent CAZyme3D truth slice. dbCAN family HMMER convention: i-Evalue<1e-15 & hmm-coverage>0.35."""
import re, json, random, os, sys
random.seed(13)

DATA="/array1/xinpeng/dbcan4-advanced/data"
REF_FAA=f"{DATA}/reference_2024.faa"
EVAL_FAA=f"{DATA}/eval_2025.faa"
EVAL_LAB=f"{DATA}/eval_2025_labels.tsv"
CAZ3D_DIR="/array1/xinpeng/cazyme3d/extracted/cazyme_id50"
CAZ3D_FAA="/array1/xinpeng/cazyme3d/cazyme3d_id50_aa.faa"
FAM_RE=re.compile(r'^(GH|GT|PL|CE|AA|CBM)\d+')
IEVAL_MAX=1e-15
COV_MIN=0.35

def read_fasta(path):
    hid,buf=None,[]
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if hid is not None: yield hid,"".join(buf)
                hid=line[1:].split()[0]; buf=[]
            else: buf.append(line.strip())
        if hid is not None: yield hid,"".join(buf)

def fam_of(name):
    n=name[:-4] if name.endswith(".hmm") else name
    m=FAM_RE.match(n)
    return m.group(0) if m else n

def parse_domtbl(path):
    """yield dict per domain hit (env coords on the SEQUENCE=target)."""
    out=[]
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            f=line.split()
            if len(f)<22: continue
            try:
                tname=f[0]; qname=f[3]; qlen=int(f[5])
                ieval=float(f[12]); hmm_from=int(f[15]); hmm_to=int(f[16])
                env_from=int(f[19]); env_to=int(f[20])
            except (ValueError,IndexError): continue
            cov=(hmm_to-hmm_from+1)/max(qlen,1)
            if ieval<=IEVAL_MAX and cov>=COV_MIN:
                out.append(dict(pid=tname, fam=fam_of(qname), hmm_name=qname[:-4] if qname.endswith('.hmm') else qname,
                                ieval=ieval, env_from=env_from, env_to=env_to, score=float(f[13])))
    return out

def resolve_overlaps(doms):
    """greedy: sort by ieval asc; accept if <50% overlap w/ any accepted (of shorter)."""
    doms=sorted(doms, key=lambda d:(d["ieval"], -(d["env_to"]-d["env_from"])))
    kept=[]
    for d in doms:
        a0,a1=d["env_from"],d["env_to"]; la=a1-a0+1; ok=True
        for k in kept:
            b0,b1=k["env_from"],k["env_to"]
            ov=max(0,min(a1,b1)-max(a0,b0)+1)
            if ov> 0.5*min(la,(b1-b0+1)): ok=False; break
        if ok: kept.append(d)
    return sorted(kept, key=lambda d:d["env_from"])

def build_domains(domtbl, faa, tag):
    hits={}
    for h in parse_domtbl(domtbl): hits.setdefault(h["pid"],[]).append(h)
    seqs=dict(read_fasta(faa))
    rows=[]; fasta=[]; per_prot={}
    for pid,ds in hits.items():
        seq=seqs.get(pid)
        if not seq: continue
        kept=resolve_overlaps(ds); per_prot[pid]=[k["fam"] for k in kept]
        for i,d in enumerate(kept):
            s=max(1,d["env_from"]-10); e=min(len(seq),d["env_to"]+10)  # pad envelope
            sub=seq[s-1:e]
            if len(sub)<20: continue
            did=f"{pid}__{d['fam']}_{s}_{e}"
            is_cbm=d["fam"].startswith("CBM")
            rows.append((did,pid,d["fam"],d["hmm_name"],s,e,len(sub),int(is_cbm),f"{d['ieval']:.2e}"))
            fasta.append((f"{did}|{d['fam']}", sub))
    return rows, fasta, per_prot

def write_tsv(path, header, rows):
    with open(path,"w") as f:
        f.write("\t".join(header)+"\n")
        for r in rows: f.write("\t".join(str(x) for x in r)+"\n")
def write_faa(path, recs):
    with open(path,"w") as f:
        for h,s in recs: f.write(f">{h}\n{s}\n")

HDR=["domain_id","parent_protein_id","family","hmm_name","env_from","env_to","domain_len","is_cbm","i_evalue"]

# ---- TRAIN (reference_2024) ----
tr_rows, tr_faa, tr_perprot = build_domains("ref2024.domtbl", REF_FAA, "train")
write_tsv("domains_train.tsv", HDR, tr_rows); write_faa("domains_train.faa", tr_faa)

# ---- EVAL (2025) ----
ev_rows, ev_faa, ev_perprot = build_domains("eval2025.domtbl", EVAL_FAA, "eval")
write_tsv("domains_eval.tsv", HDR, ev_rows); write_faa("domains_eval.faa", ev_faa)

# multidomain eval proteins per CURATED labels (independent of my hmmsearch)
mdom=set()
with open(EVAL_LAB) as f:
    next(f)
    for line in f:
        p=line.rstrip("\n").split("\t")
        if len(p)>=3 and len([x for x in p[2].split(",") if x])>=2: mdom.add(p[0])
n_mdom_recovered=sum(1 for pid in mdom if len(ev_perprot.get(pid,[]))>=2)

# per-family domain counts (train)
from collections import Counter
fc=Counter(r[2] for r in tr_rows)
json.dump(dict(sorted(fc.items(), key=lambda kv:-kv[1])), open("domain_family_counts.json","w"), indent=2)

# ---- INDEPENDENT TRUTH SLICE from CAZyme3D structure families ----
acc2fam={}
for sub in os.listdir(CAZ3D_DIR):
    if not sub.endswith("_id50"): continue
    fam=sub[:-5]
    d=os.path.join(CAZ3D_DIR,sub)
    if not os.path.isdir(d): continue
    for fn in os.listdir(d):
        if fn.endswith(".pdb"): acc2fam[fn[:-4]]=fam
# sample ~1200 across families, cap 15/family
byfam={}
for acc,fam in acc2fam.items(): byfam.setdefault(fam,[]).append(acc)
pick=[]
for fam,accs in byfam.items():
    random.shuffle(accs); pick += [(a,fam) for a in accs[:15]]
random.shuffle(pick); pick=pick[:1200]
pickset={a for a,_ in pick}
caz_seq={}
for hid,seq in read_fasta(CAZ3D_FAA):
    if hid in pickset: caz_seq[hid]=seq
trows=[]; tfaa=[]
for acc,fam in pick:
    s=caz_seq.get(acc)
    if s and 20<=len(s)<=1500:
        trows.append((acc,fam,1,len(s),"cazyme3d_structure"))
        tfaa.append((f"{acc}|{fam}", s))
write_tsv("domain_truth_slice.tsv", ["protein_id","family","boundary_start","boundary_end","truth_source"], trows)
write_faa("domain_truth_slice.faa", tfaa)

summary=dict(n_train_domains=len(tr_rows), n_train_parent_proteins=len(tr_perprot),
             n_eval_domains=len(ev_rows), n_eval_parent_proteins=len(ev_perprot),
             n_multidomain_eval_curated=len(mdom), n_multidomain_eval_recovered_ge2dom=n_mdom_recovered,
             n_cbm_train=sum(r[7] for r in tr_rows), n_families_train=len(fc),
             top10_families={k:fc[k] for k in list(dict(sorted(fc.items(),key=lambda kv:-kv[1])))[:10]},
             truth_slice_n=len(trows), truth_slice_families=len({f for _,f in pick}))
json.dump(summary, open("domain_extraction_summary.json","w"), indent=2)
print(json.dumps(summary, indent=2))
