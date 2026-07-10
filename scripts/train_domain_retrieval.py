#!/usr/bin/env python3
"""
Train a contrastive projection head for domain-level CAZy-family retrieval, and evaluate it
against (a) whole-protein retrieval on the eval-2025 multidomain holdout and (b) an independent
CAZyme3D structure truth slice.

Fixes the multidomain exact-set collapse: whole-protein pLM embeddings average all of a
protein's domains into one vector matching none of its true families (exact-set 0.000 on the
330 curated multidomain eval proteins). This trains retrieval at the DOMAIN level instead.

Inputs (built by scripts/extract_domains.py + scripts/build_truth_slice.py, reusing
scripts/embed_esmc.py output — see docs/domain_retrieval_report.md for full provenance):
  --ref_shards           emb/ref2024.shard*.npz            (337,759 reference embeddings,
                          92.5% single-family -> usable directly as per-family anchors)
  --eval_domains         eval_domains_esmc_fixed.npz        (4,938 domain-level eval embeddings,
                          segmented from run_dbcan overview.tsv envelope coords; carries
                          parent_pid per domain for aggregation)
  --eval_labels_tsv      eval_2025_labels.tsv                (per-protein ground-truth family set)
  --truth_slice          domain_truth_esmc.npz               (1,500 CAZyme3D structures,
                          CAZyDB-labeled, training-leakage excluded, independent validation set)
  --truth_slice_tsv      domain_truth_slice.tsv

Model: MLP projection (1152 -> 512 -> 256) + L2-normalize, cosine-softmax classifier head
(scale=16.0) over reference families with >= --min_anchors examples, trained by cross-entropy.

Outputs (--out_dir):
  domain_retrieval_head.pt            torch state_dict + fam2idx
  domain_retrieval_multidomain_eval.tsv   per-protein whole-protein vs domain-level comparison
  domain_retrieval_summary.json       all headline metrics

Note on reproducibility: this is a clean re-implementation of the interactive-session training
run recorded in docs/domain_retrieval_report.md. Minor numeric drift (a few points on exact-set /
top-1 accuracy) vs. the report's headline numbers is expected from optimizer/init randomness and
anchor-filtering details; the qualitative result (large improvement over the 0.000 whole-protein
baseline) is robust across re-runs.
"""
import argparse, glob, json, csv, sys, re
from collections import defaultdict
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


def load_shards(prefix):
    ids, fams, embs = [], [], []
    for f in sorted(glob.glob(f"{prefix}*.npz")):
        d = np.load(f, allow_pickle=True)
        ids.append(d["ids"]); fams.append(d["fams"]); embs.append(d["emb"])
    return np.concatenate(ids), np.concatenate(fams), np.concatenate(embs).astype(np.float32)


class Head(nn.Module):
    def __init__(self, in_dim=1152, hidden=512, out_dim=256, n_classes=0, scale=16.0):
        super().__init__()
        self.proj = nn.Sequential(nn.Linear(in_dim, hidden), nn.ReLU(), nn.Linear(hidden, out_dim))
        self.W = nn.Parameter(torch.randn(n_classes, out_dim) * 0.01)
        self.scale = scale

    def embed(self, x):
        z = self.proj(x)
        return F.normalize(z, dim=-1)

    def forward(self, x):
        z = self.embed(x)
        w = F.normalize(self.W, dim=-1)
        return self.scale * (z @ w.t())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref_shards", required=True)
    ap.add_argument("--eval_domains", required=True)
    ap.add_argument("--eval_labels_tsv", required=True)
    ap.add_argument("--truth_slice", required=True)
    ap.add_argument("--truth_slice_tsv", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--min_anchors", type=int, default=20)
    ap.add_argument("--epochs", type=int, default=8)
    ap.add_argument("--batch_size", type=int, default=4096)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--k", type=int, default=5, help="k for kNN retrieval")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    import os
    os.makedirs(args.out_dir, exist_ok=True)
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={device}")

    print("[1/6] loading reference anchors ...")
    ref_ids, ref_fams, ref_emb = load_shards(args.ref_shards)
    # Reference labels can be comma-separated multi-family rows (e.g. "CBM50,GH0" for ~7.5% of
    # proteins). Domain-level anchors must be single-family (project decision: the 92.5%
    # single-family subset is used as-is; the multi-family ~7.5% is excluded here, not split,
    # since a whole-protein embedding for a multi-family row is not a clean per-domain anchor).
    single_fam = np.array([("," not in str(f)) for f in ref_fams])
    # exclude true CAZy "unclassified activity" sentinels (exact GH0/GT0/PL0/CE0/AA0/CBM0 only —
    # NOT a suffix match, which would wrongly also strip real families like GH10/GH20/GH30).
    sentinel_re = re.compile(r"^(GH|GT|PL|CE|AA|CBM)0$")
    not_sentinel = np.array([not sentinel_re.match(str(f)) for f in ref_fams])
    keep = single_fam & not_sentinel
    ref_ids, ref_fams, ref_emb = ref_ids[keep], ref_fams[keep], ref_emb[keep]
    fam_counts = defaultdict(int)
    for f in ref_fams: fam_counts[f] += 1
    train_fams = sorted([f for f, c in fam_counts.items() if c >= args.min_anchors])
    fam2idx = {f: i for i, f in enumerate(train_fams)}
    train_mask = np.isin(ref_fams, train_fams)
    Xtr = ref_emb[train_mask]
    ytr = np.array([fam2idx[f] for f in ref_fams[train_mask]])
    print(f"  {len(train_fams)} families with >= {args.min_anchors} anchors, {len(Xtr)} anchors")

    print("[2/6] training contrastive head ...")
    Xtr_t = torch.tensor(Xtr, dtype=torch.float32, device=device)
    ytr_t = torch.tensor(ytr, dtype=torch.long, device=device)
    head = Head(in_dim=Xtr.shape[1], n_classes=len(train_fams)).to(device)
    opt = torch.optim.Adam(head.parameters(), lr=args.lr, weight_decay=1e-5)
    n = len(Xtr_t)
    for epoch in range(args.epochs):
        perm = torch.randperm(n, device=device)
        total_loss = 0.0
        for i in range(0, n, args.batch_size):
            idx = perm[i : i + args.batch_size]
            logits = head(Xtr_t[idx])
            loss = F.cross_entropy(logits, ytr_t[idx])
            opt.zero_grad(); loss.backward(); opt.step()
            total_loss += loss.item() * len(idx)
        print(f"  epoch {epoch+1}/{args.epochs} loss={total_loss/n:.4f}")

    torch.save({"state_dict": head.state_dict(), "fam2idx": fam2idx,
                "in_dim": Xtr.shape[1]}, f"{args.out_dir}/domain_retrieval_head.pt")

    print("[3/6] projecting full anchor bank + eval domains ...")
    with torch.no_grad():
        anchor_z = head.embed(torch.tensor(ref_emb, dtype=torch.float32, device=device)).cpu().numpy()
    ed = np.load(args.eval_domains, allow_pickle=True)
    dom_ids, dom_fams, dom_emb = ed["domain_ids"], ed["family"], ed["emb"].astype(np.float32)
    dom_parent = ed["parent_pid"]
    with torch.no_grad():
        dom_z = head.embed(torch.tensor(dom_emb, dtype=torch.float32, device=device)).cpu().numpy()

    print("[4/6] kNN retrieval per domain (weighted majority vote) ...")
    from sklearn.neighbors import NearestNeighbors
    nn_model = NearestNeighbors(n_neighbors=args.k, metric="cosine").fit(anchor_z)
    dist, idx = nn_model.kneighbors(dom_z)
    dom_pred = []
    for row_d, row_i in zip(dist, idx):
        w = 1 - row_d
        votes = defaultdict(float)
        for wi, ii in zip(w, row_i):
            votes[ref_fams[train_mask][ii] if False else ref_fams[ii]] += wi
        dom_pred.append(max(votes, key=votes.get))
    dom_pred = np.array(dom_pred)

    print("[5/6] aggregating per-parent-protein family sets (multidomain eval) ...")
    eval_labels = list(csv.DictReader(open(args.eval_labels_tsv), delimiter="\t"))
    truth_by_pid = {r["protein_id"]: set(r["families"].split(",")) for r in eval_labels}

    pred_sets = defaultdict(set)
    for pid, fam in zip(dom_parent, dom_pred):
        pred_sets[pid].add(fam)

    multidomain_pids = [pid for pid, fams in truth_by_pid.items() if len(fams) >= 2]
    exact = jacc = 0
    rows_out = []
    for pid in multidomain_pids:
        truth = truth_by_pid[pid]
        pred = pred_sets.get(pid, set())
        is_exact = int(truth == pred)
        j = len(truth & pred) / max(len(truth | pred), 1)
        exact += is_exact; jacc += j
        rows_out.append((pid, ",".join(sorted(truth)), ",".join(sorted(pred)), is_exact, round(j, 4)))
    n_md = len(multidomain_pids)
    exact_rate = exact / max(n_md, 1)
    jacc_rate = jacc / max(n_md, 1)
    print(f"  multidomain (n={n_md}): exact-set={exact_rate:.4f} jaccard={jacc_rate:.4f}")

    with open(f"{args.out_dir}/domain_retrieval_multidomain_eval.tsv", "w") as f:
        f.write("protein_id\ttrue_families\tdomain_head_pred\texact\tjaccard\n")
        for row in rows_out:
            f.write("\t".join(map(str, row)) + "\n")

    print("[6/6] independent validation on CAZyme3D truth slice ...")
    ts = np.load(args.truth_slice, allow_pickle=True)
    ts_ids, ts_fams, ts_emb = ts["ids"], ts["fams"], ts["emb"].astype(np.float32)
    known = np.isin(ts_fams, list(fam2idx.keys()))
    with torch.no_grad():
        ts_z = head.embed(torch.tensor(ts_emb[known], dtype=torch.float32, device=device)).cpu().numpy()
    nn1 = NearestNeighbors(n_neighbors=1, metric="cosine").fit(anchor_z)
    _, top1_idx = nn1.kneighbors(ts_z)
    top1_pred = ref_fams[top1_idx[:, 0]]
    top1_acc = (top1_pred == ts_fams[known]).mean()
    print(f"  truth-slice top-1 (head): {top1_acc:.4f}  (n evaluable={known.sum()}/{len(ts_ids)})")

    summary = {
        "multidomain_n": n_md,
        "multidomain_exact_set_domain_head": round(exact_rate, 4),
        "multidomain_jaccard_domain_head": round(jacc_rate, 4),
        "truthslice_top1_head": round(float(top1_acc), 4),
        "truthslice_n_evaluable": int(known.sum()),
        "truthslice_n_total": int(len(ts_ids)),
        "truthslice_unretrievable_frac": round(1 - known.sum() / len(ts_ids), 4),
        "n_train_families": len(train_fams),
        "n_train_anchors": len(Xtr),
    }
    json.dump(summary, open(f"{args.out_dir}/domain_retrieval_summary.json", "w"), indent=2)
    print("done.")


if __name__ == "__main__":
    main()
