#!/usr/bin/env python3
"""Independent domain-truth slice: CAZyme3D_id50 structures labeled by CAZy family
via the all-kingdom 2024 CAZyDB (accession match, else MD5 sequence match).
Family labels are CAZy-curated (independent of dbCAN HMMER on eval). Exclude any
sequence whose MD5 is in reference_2024 training (no leakage)."""
import os, re, random, hashlib, csv
random.seed(13)
CAZ3D_FAA="/array1/xinpeng/cazyme3d/cazyme3d_id50_aa.faa"
CAZYDB="/array1/xinpeng/dbcan_db_2024/CAZyDB.07142024.fa"
REFLAB="/array1/xinpeng/dbcan4-advanced/data/reference_labels_2024.tsv"
FAM_RE=re.compile(r'^(GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?$')

def md5(s): return hashlib.md5(s.encode()).hexdigest()
def read_fasta(path):
    hid,buf=None,[]
    for line in open(path):
        if line.startswith(">"):
            if hid is not None: yield hid,"".join(buf)
            hid=line[1:].rstrip("\n"); buf=[]
        else: buf.append(line.strip())
    if hid is not None: yield hid,"".join(buf)

# 1) CAZyDB: accession -> set(families), and md5 -> set(families)
acc2fam={}; md52fam={}
cur_acc=None; cur_seq=[]
def flush(acc,seq):
    if not acc or not seq: return
    s="".join(seq); fams={f for f in acc2fam.get(acc,set())}
for hid,seq in read_fasta(CAZYDB):
    # header like ACC|FAM (may repeat acc across families)
    parts=hid.split("|")
    acc=parts[0]; fam=parts[1] if len(parts)>1 else None
    if fam and FAM_RE.match(fam):
        acc2fam.setdefault(acc,set()).add(fam)
        md52fam.setdefault(md5(seq),set()).add(fam)
print(f"CAZyDB accessions={len(acc2fam)} md5s={len(md52fam)}")

# 2) reference_2024 training md5 set (to exclude leakage)
train_md5=set()
with open(REFLAB) as f:
    r=csv.DictReader(f,delimiter="\t")
    for row in r: train_md5.add(row["seq_md5"])
print(f"training md5s={len(train_md5)}")

# 3) walk CAZyme3D sequences, assign family, drop leakage + multi-family ambiguity
cand=[]  # (acc, fam, seq)
n_acc_hit=n_md5_hit=n_leak=n_ambig=n_none=0
for acc,seq in read_fasta(CAZ3D_FAA):
    acc=acc.split()[0]
    if not seq or not (20<=len(seq)<=1500): continue
    m=md5(seq)
    if m in train_md5: n_leak+=1; continue          # exclude training leakage
    fams=None
    if acc in acc2fam: fams=acc2fam[acc]; src="acc"
    elif m in md52fam: fams=md52fam[m]; src="md5"
    if not fams: n_none+=1; continue
    if len(fams)!=1: n_ambig+=1; continue            # keep single-family for clean truth
    fam=next(iter(fams))
    cand.append((acc,fam,seq))
    if src=="acc": n_acc_hit+=1
    else: n_md5_hit+=1
print(f"acc_hits={n_acc_hit} md5_hits={n_md5_hit} leak_excluded={n_leak} ambiguous={n_ambig} unlabeled={n_none} candidates={len(cand)}")

# 4) stratified sample: cap 12/family, target 1500
byfam={}
for acc,fam,seq in cand: byfam.setdefault(fam,[]).append((acc,seq))
pick=[]
for fam,lst in byfam.items():
    random.shuffle(lst); pick+=[(a,fam,s) for a,s in lst[:12]]
random.shuffle(pick); pick=pick[:1500]

with open("domain_truth_slice.faa","w") as fo, open("domain_truth_slice.tsv","w") as ft:
    ft.write("protein_id\tfamily\tseq_len\ttruth_source\n")
    n=0
    for acc,fam,seq in pick:
        fo.write(f">{acc}|{fam}\n{seq}\n")
        ft.write(f"{acc}\t{fam}\t{len(seq)}\tcazyme3d_structure_cazydb_family\n"); n+=1
print(f"truth slice written: {n} structures, {len({f for _,f,_ in pick})} families")
