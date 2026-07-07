#!/usr/bin/env python3
"""
Train two lightweight heads on FROZEN ESM-C embeddings for CAZy-family assignment:

  (A) Supervised-contrastive projection head (MLP -> L2-normalized proj), trained with a
      SupCon loss (Khosla et al.) so same-family embeddings cluster and different-family
      separate. Evaluated by kNN / nearest-centroid in the PROJECTED space.
  (B) Linear/MLP classifier head over the reference families, softmax cross-entropy.
      Max-softmax probability is the confidence / novelty (abstain) score.

Frozen embeddings mean this trains in minutes on one GPU. Reference (2024) embeddings are the
train set; a held-out slice of reference is the val set; eval_2025 embeddings are the test set.

Reports, per novelty bucket on eval_2025:
  - contrastive kNN & centroid exact/overlap family recall (vs off-the-shelf baseline)
  - classifier top-1 exact + max-softmax novelty-detection AUROC (novel_family vs novel_seq)
Saves head weights + a metrics JSON + projected eval embeddings for the fusion step.
"""
import argparse, glob, json, sys, time
import numpy as np

def load_shards(prefix):
    ids, fams, embs = [], [], []
    for f in sorted(glob.glob(f"{prefix}.shard*.npz")):
        d = np.load(f, allow_pickle=True)
        ids.append(d["ids"]); fams.append(d["fams"]); embs.append(d["emb"])
    return (np.concatenate(ids), np.concatenate(fams),
            np.concatenate(embs).astype(np.float32))

def primary_fam(f):
    toks = [x for x in f.replace("|",",").split(",") if x]
    return toks[0] if toks else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-prefix", required=True)
    ap.add_argument("--eval-prefix", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--proj-dim", type=int, default=256)
    ap.add_argument("--hidden", type=int, default=1024)
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--batch", type=int, default=4096)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--temp", type=float, default=0.1)
    ap.add_argument("--min-count", type=int, default=2, help="drop families with < this many ref seqs")
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    import torch, torch.nn as nn, torch.nn.functional as F
    torch.manual_seed(args.seed); np.random.seed(args.seed)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    import os; os.makedirs(args.outdir, exist_ok=True)

    # ---- data ----
    rid, rfam_raw, remb = load_shards(args.ref_prefix)
    rfam = np.array([primary_fam(f) for f in rfam_raw])
    keep = rfam != None
    rid, rfam, remb = rid[keep], rfam[keep], remb[keep]
    # family vocab with min-count
    from collections import Counter
    cnt = Counter(rfam)
    fams = sorted([f for f,c in cnt.items() if c >= args.min_count])
    f2i = {f:i for i,f in enumerate(fams)}
    m = np.array([f in f2i for f in rfam])
    rid, rfam, remb = rid[m], rfam[m], remb[m]
    y = np.array([f2i[f] for f in rfam])
    print(f"ref: {remb.shape}, families(>= {args.min_count}): {len(fams)}", file=sys.stderr, flush=True)

    # normalize inputs
    def l2(x):
        n = np.linalg.norm(x,axis=1,keepdims=True); n[n==0]=1; return x/n
    remb = l2(remb)

    # train/val split (stratified-ish: 5% val)
    rng = np.random.default_rng(args.seed)
    idx = rng.permutation(len(y)); nval = max(2000, int(0.05*len(y)))
    val_idx, tr_idx = idx[:nval], idx[nval:]
    Xtr = torch.tensor(remb[tr_idx], device=dev); ytr = torch.tensor(y[tr_idx], device=dev)
    Xval = torch.tensor(remb[val_idx], device=dev); yval = torch.tensor(y[val_idx], device=dev)

    D = remb.shape[1]; C = len(fams)

    # ---- models ----
    class Proj(nn.Module):
        def __init__(s):
            super().__init__()
            s.net = nn.Sequential(nn.Linear(D,args.hidden), nn.GELU(),
                                  nn.Linear(args.hidden,args.proj_dim))
        def forward(s,x): return F.normalize(s.net(x), dim=1)
    class Clf(nn.Module):
        def __init__(s):
            super().__init__()
            s.net = nn.Sequential(nn.Linear(D,args.hidden), nn.GELU(),
                                  nn.Dropout(0.1), nn.Linear(args.hidden,C))
        def forward(s,x): return s.net(x)

    proj = Proj().to(dev); clf = Clf().to(dev)

    def supcon(z, labels, temp):
        # z: (B,d) normalized; labels: (B,)
        sim = z @ z.T / temp
        sim = sim - sim.max(1, keepdim=True).values.detach()
        exp = torch.exp(sim)
        B = z.shape[0]
        mask_self = torch.eye(B, device=z.device).bool()
        exp = exp.masked_fill(mask_self, 0)
        pos = (labels[:,None]==labels[None,:]) & ~mask_self
        denom = exp.sum(1)
        log_prob = sim - torch.log(denom + 1e-12)
        # mean over positives
        pos_cnt = pos.sum(1)
        valid = pos_cnt > 0
        lp = (log_prob * pos).sum(1)[valid] / pos_cnt[valid]
        return -lp.mean()

    opt = torch.optim.AdamW(list(proj.parameters())+list(clf.parameters()), lr=args.lr, weight_decay=1e-4)
    ce = nn.CrossEntropyLoss()
    ntr = Xtr.shape[0]
    t0=time.time()
    for ep in range(args.epochs):
        perm = torch.randperm(ntr, device=dev)
        proj.train(); clf.train(); tot=0.0; nb=0
        for b in range(0, ntr, args.batch):
            bi = perm[b:b+args.batch]
            xb, yb = Xtr[bi], ytr[bi]
            z = proj(xb); logits = clf(xb)
            loss = supcon(z, yb, args.temp) + ce(logits, yb)
            opt.zero_grad(); loss.backward(); opt.step()
            tot+=loss.item(); nb+=1
        # val clf acc
        proj.eval(); clf.eval()
        with torch.no_grad():
            va = (clf(Xval).argmax(1)==yval).float().mean().item()
        if ep%5==0 or ep==args.epochs-1:
            print(f"epoch {ep:3d} loss {tot/nb:.4f} val_clf_acc {va:.4f} ({time.time()-t0:.0f}s)",
                  file=sys.stderr, flush=True)

    torch.save({"proj":proj.state_dict(),"clf":clf.state_dict(),"fams":fams,
                "args":vars(args)}, f"{args.outdir}/heads.pt")

    # ---- project reference + eval; evaluate ----
    truth, novelty = {}, {}
    with open(args.labels) as fh:
        h=fh.readline().rstrip("\n").split("\t"); fi,ni,pi=h.index("families"),h.index("novelty"),h.index("protein_id")
        for line in fh:
            p=line.rstrip("\n").split("\t")
            truth[p[pi]]=set(x for x in p[fi].split(",") if x); novelty[p[pi]]=p[ni]

    proj.eval(); clf.eval()
    with torch.no_grad():
        Zref = proj(torch.tensor(remb, device=dev)).cpu().numpy()
        eid, efam, eemb = load_shards(args.eval_prefix)
        eemb = l2(eemb)
        Xe = torch.tensor(eemb, device=dev)
        Ze = proj(Xe).cpu().numpy()
        logits_e = clf(Xe)
        probs_e = torch.softmax(logits_e,1).cpu().numpy()
        clf_pred_i = probs_e.argmax(1); clf_conf = probs_e.max(1)

    np.savez_compressed(f"{args.outdir}/proj_eval.npz", ids=eid, fams=efam, emb=Ze.astype(np.float16))
    np.savez_compressed(f"{args.outdir}/proj_ref.npz", ids=rid, fams=rfam, emb=Zref.astype(np.float16))

    # centroids in projected space
    from collections import defaultdict
    fam2idx=defaultdict(list)
    for i,f in enumerate(rfam): fam2idx[f].append(i)
    cf=sorted(fam2idx); cent=np.zeros((len(cf),args.proj_dim),dtype=np.float32)
    for j,f in enumerate(cf):
        v=Zref[fam2idx[f]].mean(0); nv=np.linalg.norm(v); cent[j]=v/nv if nv>0 else v

    # score eval
    def _avg_ranks(x):
        x=np.asarray(x,float); order=np.argsort(x,kind="mergesort"); xs=x[order]
        ranks=np.empty(len(x)); i=0
        while i<len(x):
            j=i
            while j+1<len(x) and xs[j+1]==xs[i]: j+=1
            ranks[order[i:j+1]]=(i+j)/2.0+1.0; i=j+1
        return ranks
    def auroc(score,pos):
        pos=np.asarray(pos,bool); npos=int(pos.sum()); nneg=int((~pos).sum())
        if npos==0 or nneg==0: return None
        r=_avg_ranks(score)
        return float((r[pos].sum()-npos*(npos+1)/2)/(npos*nneg))

    nov_arr=np.array([novelty[p] for p in eid])
    # contrastive centroid pred
    csims = Ze @ cent.T
    cen_pred = np.array([cf[i] for i in csims.argmax(1)])
    cen_conf = csims.max(1)
    cen_margin = csims.max(1) - np.partition(csims,-2,axis=1)[:,-2]
    # contrastive kNN (projected)
    ksims = Ze @ Zref.T
    from collections import Counter as Ctr
    knn_pred=[]; knn_pur=[]
    topk = np.argpartition(-ksims, args.k, axis=1)[:,:args.k]
    for r in range(Ze.shape[0]):
        idx=topk[r]; order=idx[np.argsort(-ksims[r,idx])]
        nf=rfam[order]; v=Ctr(nf); bf,bc=v.most_common(1)[0]
        knn_pred.append(bf); knn_pur.append(bc/args.k)
    knn_pred=np.array(knn_pred); knn_pur=np.array(knn_pur)
    clf_pred = np.array([fams[i] for i in clf_pred_i])

    def bucket_metrics(pred, qmask):
        n=exact=overlap=0
        for i,pid in enumerate(eid):
            if not qmask(nov_arr[i]): continue
            n+=1; tf=truth[pid]; pf={pred[i]} if pred[i] else set()
            if pf and pf==tf: exact+=1
            if tf&pf: overlap+=1
        return {"n":n,"exact":round(exact/n,4) if n else 0,"overlap":round(overlap/n,4) if n else 0}

    buckets={"overall":lambda x:True,"novel_seq":lambda x:x=="novel_seq","novel_family":lambda x:x=="novel_family"}
    m2=(nov_arr=="novel_seq")|(nov_arr=="novel_family"); pos=(nov_arr[m2]=="novel_seq")
    out={
      "n_families_trained":len(fams),
      "contrastive_knn":{b:bucket_metrics(knn_pred,f) for b,f in buckets.items()},
      "contrastive_centroid":{b:bucket_metrics(cen_pred,f) for b,f in buckets.items()},
      "classifier":{b:bucket_metrics(clf_pred,f) for b,f in buckets.items()},
      "novelty_auroc":{
         "contrastive_centroid_margin":auroc(cen_margin[m2],pos),
         "contrastive_centroid_cosine":auroc(cen_conf[m2],pos),
         "contrastive_knn_purity":auroc(knn_pur[m2],pos),
         "classifier_maxsoftmax":auroc(clf_conf[m2],pos),
      },
    }
    # persist per-query classifier confidence for fusion
    with open(f"{args.outdir}/head_eval_pred.tsv","w") as fo:
        fo.write("query_id\tnovelty\ttrue_families\tclf_pred\tclf_conf\tcontr_cent_pred\tcontr_cent_margin\tcontr_knn_pred\tcontr_knn_purity\n")
        for i,pid in enumerate(eid):
            fo.write(f"{pid}\t{nov_arr[i]}\t{','.join(sorted(truth[pid]))}\t{clf_pred[i]}\t{clf_conf[i]:.4f}\t"
                     f"{cen_pred[i]}\t{cen_margin[i]:.4f}\t{knn_pred[i]}\t{knn_pur[i]:.4f}\n")
    json.dump(out, open(f"{args.outdir}/head_metrics.json","w"), indent=2)
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
