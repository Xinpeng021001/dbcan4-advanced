#!/usr/bin/env python3
"""
ESM-C retrieval baselines for CAZy-family assignment.

Loads sharded reference (2024) and eval (2025) embeddings, L2-normalizes, then runs:
  (1) kNN         — cosine kNN over reference; predicted family = majority vote of top-k
                    neighbours; confidence = top-1 cosine similarity.
  (2) centroid    — one L2-normalized mean embedding (prototype) per family; predict
                    argmax cosine to centroids; confidence = best centroid cosine.
                    (This is the CLEAN-style nearest-prototype scheme.)

For each scheme we sweep a confidence threshold: below it the call is ABSTAIN. Abstention is
the correct behaviour on novel-family proteins (their true family is not in the 2024 reference),
so we report, per novelty bucket:
  - exact/overlap family recall among CONFIDENT calls (>= threshold)
  - abstain rate
and for novel_family specifically, abstain rate is the useful signal (high = correctly flags
"this is not a known family").

Outputs a JSON summary and a per-query TSV for the chosen operating threshold.
"""
import argparse, glob, json, sys
import numpy as np
from collections import defaultdict, Counter

def load_shards(prefix):
    ids, fams, embs = [], [], []
    files = sorted(glob.glob(f"{prefix}.shard*.npz"))
    if not files:
        sys.exit(f"no shards at {prefix}.shard*.npz")
    for f in files:
        d = np.load(f, allow_pickle=True)
        ids.append(d["ids"]); fams.append(d["fams"]); embs.append(d["emb"])
    ids = np.concatenate(ids)
    fams = np.concatenate(fams)
    emb = np.concatenate(embs).astype(np.float32)
    n = np.linalg.norm(emb, axis=1, keepdims=True); n[n==0]=1
    emb = emb / n
    return ids, fams, emb

def fam_set(fam_field):
    return set(x for x in fam_field.replace("|", ",").split(",") if x)

def primary_fam(fam_field):
    # first family token as the single-label target for centroid building
    toks = [x for x in fam_field.replace("|", ",").split(",") if x]
    return toks[0] if toks else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-prefix", required=True)
    ap.add_argument("--eval-prefix", required=True)
    ap.add_argument("--labels", required=True, help="eval labels tsv with protein_id,families,novelty")
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-pred", required=True)
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--op-threshold", type=float, default=0.5, help="confidence cut for the per-query TSV")
    ap.add_argument("--batch", type=int, default=512)
    args = ap.parse_args()

    # eval truth / novelty
    truth, novelty = {}, {}
    with open(args.labels) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        fi, ni, pi = h.index("families"), h.index("novelty"), h.index("protein_id")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            truth[p[pi]] = set(x for x in p[fi].split(",") if x)
            novelty[p[pi]] = p[ni]

    print("loading reference embeddings...", file=sys.stderr, flush=True)
    rid, rfam, remb = load_shards(args.ref_prefix)
    print(f"ref: {remb.shape}", file=sys.stderr, flush=True)
    rprimary = np.array([primary_fam(f) for f in rfam])

    print("loading eval embeddings...", file=sys.stderr, flush=True)
    eid, efam, eemb = load_shards(args.eval_prefix)
    print(f"eval: {eemb.shape}", file=sys.stderr, flush=True)

    # ---- build family centroids from reference ----
    fam2idx = defaultdict(list)
    for i, f in enumerate(rprimary):
        if f: fam2idx[f].append(i)
    cent_fams = sorted(fam2idx)
    cent = np.zeros((len(cent_fams), remb.shape[1]), dtype=np.float32)
    for j, f in enumerate(cent_fams):
        v = remb[fam2idx[f]].mean(0)
        nv = np.linalg.norm(v); cent[j] = v/nv if nv>0 else v
    print(f"centroids: {cent.shape} over {len(cent_fams)} families", file=sys.stderr, flush=True)

    # ---- scoring, batched matmul (cosine since all normalized) ----
    knn_pred, knn_conf = [], []
    knn_purity, knn_margin = [], []          # novelty scores
    cen_pred, cen_conf, cen_margin = [], [], []
    for b in range(0, eemb.shape[0], args.batch):
        q = eemb[b:b+args.batch]                       # (B,D)
        # kNN vs reference
        sims = q @ remb.T                              # (B,Nref)
        topk = np.argpartition(-sims, args.k, axis=1)[:, :args.k]
        for r in range(q.shape[0]):
            idx = topk[r]
            order = idx[np.argsort(-sims[r, idx])]
            top1 = order[0]
            neigh_fams = rprimary[order]
            votes = Counter(neigh_fams)
            best_fam, best_ct = votes.most_common(1)[0]
            knn_pred.append(best_fam)
            knn_conf.append(float(sims[r, top1]))
            # novelty scores: vote purity (fraction of k neighbours in winning family)
            knn_purity.append(best_ct / args.k)
            # similarity margin: top-1 sim minus best sim of a DIFFERENT family
            diff = [sims[r, order[j]] for j in range(len(order)) if neigh_fams[j] != best_fam]
            knn_margin.append(float(sims[r, top1] - max(diff)) if diff else 1.0)
        # centroid
        csims = q @ cent.T                             # (B,Ncent)
        for r in range(q.shape[0]):
            crow = csims[r]
            cbest = int(np.argmax(crow))
            cen_pred.append(cent_fams[cbest])
            cen_conf.append(float(crow[cbest]))
            # margin to 2nd-best centroid
            second = np.partition(crow, -2)[-2]
            cen_margin.append(float(crow[cbest] - second))
        print(f"scored {min(b+args.batch, eemb.shape[0])}/{eemb.shape[0]}", file=sys.stderr, flush=True)

    knn_pred=np.array(knn_pred); knn_conf=np.array(knn_conf)
    knn_purity=np.array(knn_purity); knn_margin=np.array(knn_margin)
    cen_pred=np.array(cen_pred); cen_conf=np.array(cen_conf); cen_margin=np.array(cen_margin)

    # novelty-detection AUROC: known (novel_seq) should score HIGH, novel_family LOW
    def auroc(score, pos):
        pos_s=score[pos]; neg_s=score[~pos]
        if len(pos_s)==0 or len(neg_s)==0: return None
        alls=np.concatenate([pos_s,neg_s]); o=np.argsort(alls)
        ranks=np.empty(len(o)); ranks[o]=np.arange(1,len(o)+1)
        return float((ranks[:len(pos_s)].sum()-len(pos_s)*(len(pos_s)+1)/2)/(len(pos_s)*len(neg_s)))
    nov_arr=np.array([novelty[p] for p in eid])
    m=(nov_arr=="novel_seq")|(nov_arr=="novel_family"); pos=(nov_arr[m]=="novel_seq")
    novelty_auroc={
        "knn_top1_cosine":auroc(knn_conf[m],pos),
        "knn_vote_purity":auroc(knn_purity[m],pos),
        "knn_margin":auroc(knn_margin[m],pos),
        "centroid_cosine":auroc(cen_conf[m],pos),
        "centroid_margin":auroc(cen_margin[m],pos),
    }
    print("novelty-detection AUROC:", json.dumps(novelty_auroc), file=sys.stderr, flush=True)

    order_ids = list(eid)
    def metrics_at(pred, conf, thr, qids_mask):
        n=0; conf_n=0; exact=0; overlap=0; abstain=0
        for i, pid in enumerate(order_ids):
            if not qids_mask(pid): continue
            n+=1
            tf = truth[pid]
            if conf[i] < thr:
                abstain+=1; continue
            conf_n+=1
            pf = {pred[i]} if pred[i] else set()
            if pf and pf==tf: exact+=1
            if tf & pf: overlap+=1
        return {"n":n, "confident":conf_n, "abstain":abstain,
                "abstain_frac": round(abstain/n,4) if n else 0,
                "exact_of_confident": round(exact/conf_n,4) if conf_n else 0,
                "overlap_of_confident": round(overlap/conf_n,4) if conf_n else 0,
                "exact_of_all": round(exact/n,4) if n else 0,
                "overlap_of_all": round(overlap/n,4) if n else 0}

    buckets = {"overall": lambda p: True}
    for nv in set(novelty.values()):
        buckets[nv] = (lambda x: (lambda p: novelty[p]==x))(nv)

    thr_grid = [0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    summary = {"k":args.k, "n_ref":int(remb.shape[0]), "n_families_ref":len(cent_fams),
               "novelty_detection_auroc": novelty_auroc, "schemes":{}}
    for name, pred, conf in [("knn",knn_pred,knn_conf), ("centroid",cen_pred,cen_conf)]:
        summary["schemes"][name] = {}
        for bname, mask in buckets.items():
            summary["schemes"][name][bname] = {f"thr={t}": metrics_at(pred,conf,t,mask) for t in thr_grid}

    with open(args.out_summary,"w") as fo:
        json.dump(summary, fo, indent=2)

    # per-query TSV at operating threshold
    with open(args.out_pred,"w") as fo:
        fo.write("query_id\tnovelty\ttrue_families\tknn_pred\tknn_conf\tknn_purity\tknn_margin\t"
                 "cent_pred\tcent_conf\tcent_margin\n")
        for i, pid in enumerate(order_ids):
            tf = ",".join(sorted(truth[pid]))
            fo.write(f"{pid}\t{novelty[pid]}\t{tf}\t{knn_pred[i]}\t{knn_conf[i]:.4f}\t"
                     f"{knn_purity[i]:.4f}\t{knn_margin[i]:.4f}\t"
                     f"{cen_pred[i]}\t{cen_conf[i]:.4f}\t{cen_margin[i]:.4f}\n")

    # console digest
    for name in ("knn","centroid"):
        print(f"\n=== {name} ===")
        for bname in ("overall","novel_seq","novel_family"):
            d = summary["schemes"][name][bname]["thr=0.0"]
            print(f"  {bname:13s} n={d['n']:5d} exact_all={d['exact_of_all']:.3f} overlap_all={d['overlap_of_all']:.3f}")

if __name__ == "__main__":
    main()
