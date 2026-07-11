#!/usr/bin/env python3
"""Label-free ESM-C inference for CAZy-family assignment (dbCAN4-advanced).

Unlike retrieval_esmc.py / train_heads.py (benchmark harnesses that REQUIRE a
ground-truth labels file to score recall), this predicts families + confidences
for an arbitrary input FASTA's embeddings against the PRECOMPUTED reference
index and the ALREADY-TRAINED heads. No truth, no novelty, no metrics — the
whole point of the product: drop in a novel sequence, get calls out.

Four schemes, one pass, each written as a normalize_predictions-ready raw TSV:
  ESM-C-kNN        cosine kNN vote over raw-space reference embeddings
  ESM-C-centroid   nearest family centroid (raw space)
  ESM-C-contrastive  trained SupCon projection: kNN over proj_ref.npz + softmax classifier
                     (we emit the classifier max-softmax call as the contrastive tool's
                      primary, with contr-kNN purity as a secondary signal)

Inputs (all precomputed on met — see nextflow.config assets/ref):
  --emb            query embeddings .npz (ids, emb[, fams])  [from embed_esmc.py]
  --ref-prefix     reference_2024.esmc.shard*.npz            [raw ESM-C ref index]
  --heads          results/heads/heads.pt                    [trained proj+clf+fams]
  --proj-ref       results/heads/proj_ref.npz                [ref in projected space]
Outputs: one raw TSV per scheme (query_id, pred, conf, purity/margin) — the
pipeline's normalize_predictions.py maps these to the §2.1 contract files.
"""
from __future__ import annotations
import argparse, glob, json, sys
import numpy as np
from collections import defaultdict, Counter


def _l2(x):
    x = x.astype(np.float32)
    n = np.linalg.norm(x, axis=1, keepdims=True); n[n == 0] = 1
    return x / n


def load_npz_index(prefix_or_file):
    """Load an embedding index from either a .npz file or a shard prefix."""
    files = [prefix_or_file] if prefix_or_file.endswith(".npz") else sorted(glob.glob(f"{prefix_or_file}.shard*.npz"))
    if not files:
        sys.exit(f"no embeddings at {prefix_or_file}")
    ids, fams, embs = [], [], []
    for f in files:
        d = np.load(f, allow_pickle=True)
        ids.append(d["ids"]); embs.append(d["emb"])
        fams.append(d["fams"] if "fams" in d.files else np.array([""] * len(d["ids"])))
    return np.concatenate(ids), np.concatenate(fams), _l2(np.concatenate(embs))


def primary_fam(field):
    toks = [x for x in str(field).replace("|", ",").split(",") if x and x != "None"]
    return toks[0] if toks else None


def knn_vote(q, ref_emb, ref_prim, k):
    """cosine kNN majority vote. Returns (pred, conf=top1 sim, purity, margin)."""
    preds, confs, purities, margins = [], [], [], []
    B = 512
    for b in range(0, q.shape[0], B):
        sims = q[b:b+B] @ ref_emb.T
        kk = min(k, ref_emb.shape[0])
        topk = np.argpartition(-sims, kk-1, axis=1)[:, :kk]
        for r in range(sims.shape[0]):
            idx = topk[r]; order = idx[np.argsort(-sims[r, idx])]
            top1 = order[0]; neigh = ref_prim[order]
            votes = Counter(neigh); best_fam, best_ct = votes.most_common(1)[0]
            preds.append(best_fam); confs.append(float(sims[r, top1]))
            purities.append(best_ct / kk)
            diff = [sims[r, order[j]] for j in range(len(order)) if neigh[j] != best_fam]
            margins.append(float(sims[r, top1] - max(diff)) if diff else 1.0)
    return preds, confs, purities, margins


def centroids(ref_emb, ref_prim):
    fam2idx = defaultdict(list)
    for i, f in enumerate(ref_prim):
        if f: fam2idx[f].append(i)
    fams = sorted(fam2idx)
    cent = np.zeros((len(fams), ref_emb.shape[1]), dtype=np.float32)
    for j, f in enumerate(fams):
        v = ref_emb[fam2idx[f]].mean(0); nv = np.linalg.norm(v)
        cent[j] = v/nv if nv > 0 else v
    return fams, cent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", required=True, help="query embeddings npz (from embed_esmc.py)")
    ap.add_argument("--ref-prefix", required=True, help="raw ESM-C reference shard prefix")
    ap.add_argument("--heads", default=None, help="trained heads.pt (enables contrastive scheme)")
    ap.add_argument("--proj-ref", default=None, help="reference embeddings in projected space (proj_ref.npz)")
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--out-knn", required=True)
    ap.add_argument("--out-centroid", required=True)
    ap.add_argument("--out-contrastive", default=None)
    args = ap.parse_args()

    qid, _, qemb = load_npz_index(args.emb)
    print(f"[infer] query: {qemb.shape}", file=sys.stderr, flush=True)
    rid, rfam, remb = load_npz_index(args.ref_prefix)
    rprim = np.array([primary_fam(f) for f in rfam])
    print(f"[infer] reference: {remb.shape}", file=sys.stderr, flush=True)

    # ---- raw-space kNN ----
    kp, kc, kpur, kmar = knn_vote(qemb, remb, rprim, args.k)
    with open(args.out_knn, "w") as fo:
        fo.write("query_id\tknn_pred\tknn_conf\tknn_purity\tknn_margin\n")
        for i, p in enumerate(qid):
            fo.write(f"{p}\t{kp[i]}\t{kc[i]:.4f}\t{kpur[i]:.4f}\t{kmar[i]:.4f}\n")
    print(f"[infer] wrote {args.out_knn}", file=sys.stderr, flush=True)

    # ---- raw-space centroid ----
    cfams, cent = centroids(remb, rprim)
    csims = qemb @ cent.T
    cbest = csims.argmax(1)
    with open(args.out_centroid, "w") as fo:
        fo.write("query_id\tcent_pred\tcent_conf\tcent_margin\n")
        for i, p in enumerate(qid):
            row = csims[i]; b = int(cbest[i])
            second = np.partition(row, -2)[-2] if len(row) > 1 else 0.0
            fo.write(f"{p}\t{cfams[b]}\t{float(row[b]):.4f}\t{float(row[b]-second):.4f}\n")
    print(f"[infer] wrote {args.out_centroid}", file=sys.stderr, flush=True)

    # ---- trained-head contrastive (optional; needs torch + heads.pt) ----
    if args.out_contrastive and args.heads:
        import torch, torch.nn as nn, torch.nn.functional as F
        ckpt = torch.load(args.heads, map_location="cpu", weights_only=False)
        fams = ckpt["fams"]; targs = ckpt.get("args", {})
        D = qemb.shape[1]; H = targs.get("hidden", 1024); PD = targs.get("proj_dim", 256); C = len(fams)

        class Proj(nn.Module):
            def __init__(s):
                super().__init__()
                s.net = nn.Sequential(nn.Linear(D, H), nn.GELU(), nn.Linear(H, PD))
            def forward(s, x): return F.normalize(s.net(x), dim=1)

        class Clf(nn.Module):
            def __init__(s):
                super().__init__()
                s.net = nn.Sequential(nn.Linear(D, H), nn.GELU(), nn.Dropout(0.1), nn.Linear(H, C))
            def forward(s, x): return s.net(x)

        proj = Proj(); proj.load_state_dict(ckpt["proj"]); proj.eval()
        clf = Clf(); clf.load_state_dict(ckpt["clf"]); clf.eval()
        with torch.no_grad():
            Xq = torch.tensor(qemb)
            zq = proj(Xq).numpy()                        # projected query
            probs = torch.softmax(clf(Xq), 1).numpy()    # classifier
        clf_pred_i = probs.argmax(1); clf_conf = probs.max(1)

        # contrastive-kNN over projected reference (if provided)
        ck_pred = ck_pur = None
        if args.proj_ref:
            pr_id, pr_fam, pr_emb = load_npz_index(args.proj_ref)
            pr_prim = np.array([primary_fam(f) for f in pr_fam])
            ck_pred, _, ck_pur, _ = knn_vote(_l2(zq), pr_emb, pr_prim, args.k)

        with open(args.out_contrastive, "w") as fo:
            fo.write("query_id\tclf_pred\tclf_conf\tcontr_knn_pred\tcontr_knn_purity\n")
            for i, p in enumerate(qid):
                ckp = ck_pred[i] if ck_pred else fams[clf_pred_i[i]]
                ckpu = ck_pur[i] if ck_pur else 0.0
                fo.write(f"{p}\t{fams[clf_pred_i[i]]}\t{clf_conf[i]:.4f}\t{ckp}\t{ckpu:.4f}\n")
        print(f"[infer] wrote {args.out_contrastive}", file=sys.stderr, flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
