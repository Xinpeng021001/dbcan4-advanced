#!/usr/bin/env python3
"""Faithful full-page poster of the comprehensive 267317 CAZyme page (real data).
Native browser render is blocked in this sandbox (Chrome/WebKit/QuickLook all
killed by seccomp), so this reproduces the exact page layout + content from the
same DB the live route reads."""
import json, numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import matplotlib
matplotlib.rcParams["font.family"]="DejaVu Sans"

D=json.load(open("handoff/ui_render_267317.json"))
rec=D["comprehensive"]; calls=D["real_calls"]
L=rec["length_aa"]; cc=rec["cazyme_call"]; fn=rec["function"]; pc=rec["physicochemistry"]
st=rec["structure"]; ts=rec["topology_secretion"]["deeptmhmm"]; loc=rec["localization"]; doms=rec["domains_pfam"]["domains"]
plddt=[float(l[60:66]) for l in open("real_structures/267317.pdb").read().splitlines() if l.startswith("ATOM") and l[12:16].strip()=="CA"]
COL={"HMMER":"#64748b","dbCAN_sub":"#64748b","DIAMOND":"#64748b","ESM-C-kNN":"#2563eb","ESM-C-centroid":"#7c3aed","ESM-C-contrastive":"#c026d3","Foldseek-CAZyme3D":"#059669"}
INK="#1a2233"; MUT="#64728c"; LINE="#dbe2ee"

fig=plt.figure(figsize=(12.4,15.8)); ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")
ax.add_patch(Rectangle((0,0),1,1,fc="#f7f9fc",ec="none",zorder=-5))
def rbox(x,y,w,h,fc,ec="none",lw=0,rs=0.008,z=1):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0.002,rounding_size={rs}",fc=fc,ec=ec,lw=lw,zorder=z))
def tag(x,y,txt,col,fs=7.0):
    w=0.0068*len(txt)+0.012; rbox(x,y-0.008,w,0.016,col,rs=0.008,z=3)
    ax.text(x+w/2,y,txt,fontsize=fs,color="white",ha="center",va="center",zorder=4,fontweight="bold"); return w

# HERO
rbox(0.03,0.958,0.94,0.035,"#123a5c",rs=0.008)
ax.text(0.046,0.985,"267317  \u00b7  multi-domain glycoside hydrolase (GH28-type + GH78 \u03b1-L-rhamnosidase)",fontsize=13,fontweight="bold",color="white",va="top")
ax.text(0.046,0.969,"Comprehensive multi-tool annotation of a dbCAN4-identified fungal CAZyme candidate",fontsize=8,color="#c9d8e8",va="top")
cx=0.046
for ch in [f"{L:,} aa",f"EC {fn['ec_number']}","Extracellular",f"pLDDT {st['plddt_mean']:.0f}",f"MW {pc['molecular_weight']/1000:.0f} kDa",f"pI {pc['theoretical_pI']}"]:
    w=0.0068*len(ch)+0.014; rbox(cx,0.945,w,0.016,"#1f4e78","#3a6491",0.6,0.007); ax.text(cx+w/2,0.953,ch,fontsize=6.6,color="white",ha="center",va="center"); cx+=w+0.006

# BANNER
rbox(0.03,0.917,0.94,0.020,"#fff8ec","#f2d5ae",1,0.007)
ax.text(0.046,0.927,"Why comprehensive annotation matters:  different tools key on different domains of this multi-domain protein \u2014 no single family label is complete.",fontsize=7.8,color="#7a4d12",va="center",style="italic")

# CAZYME CALL
rbox(0.03,0.758,0.94,0.150,"white",LINE,1,0.01)
ax.text(0.046,0.900,"CAZyme family \u2014 multi-tool evidence",fontsize=10.5,fontweight="bold",color=INK,va="top")
ax.text(0.046,0.884,"Per-tool calls on this multi-domain protein; the value is in seeing all signals together",fontsize=7,color=MUT,va="top")
ax.text(0.046,0.867,"ADVANCED (comprehensive release)",fontsize=7,color=MUT,fontweight="bold",va="top")
yy=0.851
for c in [c for c in calls if c["release_id"]==3]:
    col=COL.get(c["tool"],"#64748b"); tag(0.046,yy,c["tool"],col)
    ax.text(0.285,yy,c["cazy_family"],fontsize=8.5,fontweight="bold",color=INK,va="center")
    rbox(0.35,yy-0.005,0.09,0.011,"#eef1f6",rs=0.005); rbox(0.35,yy-0.005,0.09*float(c["confidence"]),0.011,col,rs=0.005)
    ax.text(0.448,yy,f"{float(c['confidence']):.3f}",fontsize=7.2,color=INK,va="center"); yy-=0.0185
tag(0.046,yy,"baseline HMMER / dbCAN_sub / DIAMOND  \u2192  GH28","#64748b",6.6)
ax.text(0.55,0.867,"WHAT EACH SIGNAL SEES",fontsize=7,color=MUT,fontweight="bold",va="top")
sy=0.849
for name,desc in [("baseline HMMER/DIAMOND","GH28 \u2014 N-terminal Glyco_hydro_28 domain"),("ESM-C (kNN + contrastive)","GH78 \u2014 C-terminal rhamnosidase signal"),("Pfam / hmmscan","BOTH domains (PF00295 + PF17389)"),("Foldseek / CAZyme3D","GH28-like \u2014 full-length fold, N-term dominant"),("CLEAN seq\u2192EC","EC 3.2.1.40 (independent, low conf 0.11)")]:
    ax.text(0.55,sy,f"\u2022 {name}",fontsize=7.3,color=INK,va="top",fontweight="bold"); ax.text(0.562,sy-0.011,desc,fontsize=6.8,color=MUT,va="top"); sy-=0.0206
ax.text(0.046,0.766,"GH78 refs exist in the CAZyme3D set (n=2) but do not surface among 296 hits \u2014 a genuine structural observation.",fontsize=6.2,color=MUT,va="top",style="italic")

# FUNCTION (left)
rbox(0.03,0.608,0.455,0.142,"white",LINE,1,0.01)
ax.text(0.046,0.742,"Function \u00b7 EC \u00b7 substrate",fontsize=10,fontweight="bold",color=INK,va="top")
ax.text(0.046,0.727,"Family EC corroborated by independent seq\u2192EC predictor",fontsize=6.8,color=MUT,va="top")
fy=0.710
for k,v in [("Activity",fn["activity_name"]),("EC (family)",f"{fn['ec_number']}  \u00b7  dbCAN"),("EC (CLEAN)",f"{fn['ec_source_clean']['predicted_ec']}  (conf {fn['ec_source_clean']['confidence']})"),("GO (MF)",fn['go_molecular_function']['id']),("GO name",fn['go_molecular_function']['name'][:32]),("Substrate",fn["substrate"][:38]),("Reaction",fn["reaction"][:42]+"\u2026")]:
    ax.text(0.046,fy,k,fontsize=7,color=MUT,va="top"); ax.text(0.155,fy,str(v),fontsize=7.2,color=INK,va="top"); fy-=0.0165

# SECRETION (right)
rbox(0.515,0.608,0.455,0.142,"white",LINE,1,0.01)
ax.text(0.531,0.742,"Secretion & topology",fontsize=10,fontweight="bold",color=INK,va="top")
ax.text(0.531,0.727,"DeepTMHMM \u2014 real run (DTU via BioLib cloud)",fontsize=6.8,color=MUT,va="top")
sy2=0.710
for k,v in [("Prediction",f"{ts['prediction']}  (classically secreted)"),("Signal peptide",f"residues {ts['signal_peptide_span'][0]}\u2013{ts['signal_peptide_span'][1]}"),("TM helices",str(ts["n_tm_helices"])),("SignalP 6.0","not installed (license-gated)")]:
    ax.text(0.531,sy2,k,fontsize=7,color=MUT,va="top"); ax.text(0.645,sy2,str(v),fontsize=7.2,color=INK,va="top"); sy2-=0.0175
ax.text(0.531,sy2-0.006,"Cleaved N-terminal SP with 0 TM helices \u2192 a secreted enzyme,\nconsistent with a plant-cell-wall-degrading fungal CAZyme.",fontsize=6.6,color=MUT,va="top",style="italic")

# DOMAIN ARCH
rbox(0.03,0.468,0.94,0.130,"white",LINE,1,0.01)
ax.text(0.046,0.590,"Domain architecture",fontsize=10,fontweight="bold",color=INK,va="top")
ax.text(0.046,0.576,"Pfam / hmmscan (--cut_ga)  \u00b7  pLDDT confidence track",fontsize=6.8,color=MUT,va="top")
axd=fig.add_axes([0.05,0.480,0.90,0.076]); axd.set_xlim(1,L); axd.set_ylim(-1.4,1.7); axd.axis("off")
xs=np.arange(1,len(plddt)+1)
axd.fill_between(xs,[(p/100)*0.5+0.9 for p in plddt],0.9,color="#c7d2e5",lw=0)
axd.plot(xs,[(p/100)*0.5+0.9 for p in plddt],color="#2563eb",lw=0.5)
axd.add_patch(Rectangle((1,-0.12),L-1,0.24,fc="#e2e8f2",ec="none"))
DC={"Glyco_hydro_28":"#94a3b8","Bac_rhamnosid6H":"#059669"}
for d in doms:
    s,e,nm=d["seq_start"],d["seq_end"],d["name"]; col=DC.get(nm,"#7c3aed")
    axd.add_patch(Rectangle((s,-0.32),e-s,0.64,fc=col,ec="white",lw=1))
    axd.text((s+e)/2,0,nm,ha="center",va="center",fontsize=7.5,color="white",fontweight="bold")
    axd.text((s+e)/2,-0.58,f"{d['pfam_acc_base']} \u00b7 {s}\u2013{e} \u00b7 iE {d['i_evalue']:.0e}",ha="center",va="top",fontsize=6,color="#334155")
axd.add_patch(Rectangle((1,-0.32),17,0.64,fc="#ea580c",ec="white",lw=0.6))
ngly=rec["sequence_features"]["n_glycosylation_sequons"]["sites"]
for s in ngly:
    p=s["position"] if isinstance(s,dict) else s; axd.plot([p,p],[0.32,0.46],color="#0d9488",lw=0.4)
for t in range(0,L+1,200): axd.text(max(t,1),-0.95,str(t),ha="center",va="top",fontsize=6,color="#94a3b8")
lx=0.05
for lab,col in [("signal peptide 1\u201318","#ea580c"),("Glyco_hydro_28 (GH28, N-term)","#94a3b8"),("Bac_rhamnosid6H (GH78, C-term)","#059669"),(f"N-glyc sequons (n={len(ngly)})","#0d9488")]:
    rbox(lx,0.472,0.010,0.010,col,rs=0.004,z=3); ax.text(lx+0.015,0.477,lab,fontsize=6.4,color=MUT,va="center"); lx+=0.010+0.0065*len(lab)+0.026

# FOLDSEEK (left)
rbox(0.03,0.318,0.58,0.140,"white",LINE,1,0.01)
ax.text(0.046,0.450,"Structural homologs",fontsize=10,fontweight="bold",color=INK,va="top")
ax.text(0.046,0.436,f"Foldseek vs CAZyme3D \u00b7 {st['foldseek_summary']['n_hits']} hits \u00b7 high TM at low seq-id = remote homology",fontsize=6.5,color=MUT,va="top")
hx=[0.046,0.17,0.27,0.34,0.41,0.48]; 
for x,h in zip(hx,["Target","Family","TM","prob","LDDT","bits"]): ax.text(x,0.419,h,fontsize=6.4,color=MUT,fontweight="bold",va="top")
hy=0.405
for h in st["foldseek_top_hits"][:6]:
    for x,v,mono,bold in zip(hx,[h["target"][:12],h["family"],f"{h['tmscore']:.2f}",f"{h['prob']:.2f}",f"{h['lddt']:.2f}",str(h["bits"])],[1,0,0,0,0,0],[0,1,0,0,0,0]):
        ax.text(x,hy,v,fontsize=6.6,color=INK,va="top",fontweight="bold" if bold else "normal",family="monospace" if mono else "DejaVu Sans")
    hy-=0.0135
ax.text(0.046,hy-0.003,"Top neighbours are GH28 (TM 0.80\u20130.82); full-length fold is GH28-dominated.",fontsize=6.3,color=MUT,va="top",style="italic")

# LOCALIZATION (right top)
rbox(0.63,0.392,0.34,0.066,"white",LINE,1,0.01)
ax.text(0.646,0.450,"Subcellular localization",fontsize=9.5,fontweight="bold",color=INK,va="top")
for i,(k,v) in enumerate([("Call",loc["localization_call"].split(" (")[0]),("Basis","N-terminal signal peptide"),("DeepLoc-2.0","not installed (scaffold)")]):
    ax.text(0.646,0.431-i*0.014,k,fontsize=6.7,color=MUT,va="top"); ax.text(0.74,0.431-i*0.014,str(v),fontsize=6.9,color=INK,va="top")

# PHYSICOCHEM (right bottom)
rbox(0.63,0.318,0.34,0.066,"white",LINE,1,0.01)
ax.text(0.646,0.376,"Physicochemistry",fontsize=9.5,fontweight="bold",color=INK,va="top")
for i,(k,v) in enumerate([("MW",f"{pc['molecular_weight']:,.0f} Da"),("pI",str(pc['theoretical_pI'])),("Instability",f"{pc['instability_index']} ({pc['instability_class']})"),("GRAVY",str(pc['gravy']))]):
    ax.text(0.646,0.357-i*0.0125,k,fontsize=6.7,color=MUT,va="top"); ax.text(0.74,0.357-i*0.0125,str(v),fontsize=6.9,color=INK,va="top")

# PROVENANCE
import csv as _csv; from collections import Counter
prov=list(_csv.DictReader(open("annotation_267317/evidence_provenance.tsv"),delimiter="\t"))
etot=Counter(r["evidence"] for r in prov)
rbox(0.03,0.238,0.94,0.070,"white",LINE,1,0.01)
ax.text(0.046,0.300,"Evidence provenance",fontsize=10,fontweight="bold",color=INK,va="top")
ax.text(0.046,0.286,f"Every value traced to its tool \u00b7 {len(prov)} annotations",fontsize=6.8,color=MUT,va="top")
ecol={"real_run":"#12925a","reference_lookup":"#2563eb","derived":"#7c3aed","not_installed":"#b25c00"}
elab={"real_run":"real tool run","reference_lookup":"reference lookup","derived":"derived","not_installed":"not installed"}
x0=0.046; tot=len(prov); barw=0.90
for et in ["real_run","reference_lookup","derived","not_installed"]:
    n=etot.get(et,0)
    if n==0: continue
    w=barw*n/tot; rbox(x0,0.256,w,0.017,ecol[et],rs=0.004,z=3); ax.text(x0+w/2,0.2645,str(n),fontsize=7,color="white",ha="center",va="center",zorder=4,fontweight="bold"); x0+=w
lx=0.046
for et in ["real_run","reference_lookup","derived","not_installed"]:
    rbox(lx,0.246,0.010,0.010,ecol[et],rs=0.004,z=3); ax.text(lx+0.015,0.251,f"{elab[et]} ({etot.get(et,0)})",fontsize=6.4,color=MUT,va="center"); lx+=0.20

ax.text(0.5,0.226,"Real runs: ESMFold, Pfam/hmmscan, Foldseek, DeepTMHMM, CLEAN, Biopython  \u00b7  SignalP6/DeepLoc are license-gated scaffolds  \u00b7  EC/GO/substrate from curated dbCAN + QuickGO",fontsize=6.3,color=MUT,ha="center",va="top",style="italic")

fig.savefig("screenshots/comprehensive_267317_poster.png", dpi=150, bbox_inches="tight", facecolor="#f7f9fc")
print("clean poster written")
