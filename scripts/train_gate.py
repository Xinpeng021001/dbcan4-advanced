#!/usr/bin/env python3
"""
Train the CAZyme-vs-non-CAZyme precision gate on frozen ESM-C embeddings.

Data (built by scripts/sample_negatives.py + reused reference embeddings):
  - reference_2024 CAZyme positives:  emb/ref2024.shard*.npz         (337,759 proteins)
  - decoy/negative set:               negatives/negatives_esmc.npz    (112,818: 97,818 natural
                                       non-CAZyme + 15,000 shuffled-domain decoys)
  - held-out realistic-imbalance slice (NEVER used for training):
                                       negatives/realistic_slice_esmc.npz (32,320 proteins,
                                       2.07% CAZyme, 5 whole genomes / 5 tax classes, zero
                                       genome overlap with training positives or negatives)

Model: L2-normalize embeddings -> LogisticRegression(C=1.0, max_iter=2000).

Reports (written to --out_dir):
  - cazyme_gate_model.joblib          trained sklearn model
  - gate_pr_curve.npz                 PR curve arrays for the realistic slice + clean cazyme/non
  - gate_operating_points.tsv         precision/recall/F1 at a grid of thresholds
  - gate_calibration.json             t_lo / t_hi abstention thresholds + tier separation stats
  - gate_grayzone_adjudication.json   gray-zone triage counts at the calibrated thresholds
  - gate_eval_realistic_slice.tsv     per-protein score + truth for the realistic slice

Usage:
  python train_gate.py \
    --ref_shards emb/ref2024.shard \
    --negatives negatives/negatives_esmc.npz \
    --realistic_slice negatives/realistic_slice_esmc.npz \
    --realistic_slice_tsv negatives/realistic_imbalance_slice.tsv \
    --out_dir gate_out/
"""
import argparse, glob, json, sys
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve
import joblib


def load_shards(prefix):
    ids, fams, embs = [], [], []
    for f in sorted(glob.glob(f"{prefix}*.npz")):
        d = np.load(f, allow_pickle=True)
        ids.append(d["ids"]); fams.append(d["fams"]); embs.append(d["emb"])
    return np.concatenate(ids), np.concatenate(fams), np.concatenate(embs).astype(np.float32)


def l2norm(x):
    n = np.linalg.norm(x, axis=1, keepdims=True)
    n[n == 0] = 1.0
    return x / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref_shards", required=True, help="prefix for reference_2024 embedding shards")
    ap.add_argument("--negatives", required=True, help="negatives_esmc.npz (natural + shuffled decoys)")
    ap.add_argument("--realistic_slice", required=True, help="realistic_slice_esmc.npz (held-out eval)")
    ap.add_argument("--realistic_slice_tsv", required=True, help="realistic_imbalance_slice.tsv (truth column)")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--n_train_pos", type=int, default=90000)
    ap.add_argument("--n_train_neg", type=int, default=90000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    import os
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    print("[1/5] loading reference (positive) embeddings ...")
    pos_ids, pos_fams, pos_emb = load_shards(args.ref_shards)
    pos_emb = l2norm(pos_emb)
    print(f"  {len(pos_ids)} reference CAZyme positives")

    print("[2/5] loading negative/decoy embeddings ...")
    neg_d = np.load(args.negatives, allow_pickle=True)
    neg_ids, neg_emb = neg_d["ids"], l2norm(neg_d["emb"].astype(np.float32))
    print(f"  {len(neg_ids)} negatives (natural non-CAZyme + shuffled decoys)")

    print("[3/5] train/held-out split (no overlap with realistic slice by construction)")
    pos_idx = rng.permutation(len(pos_ids))
    neg_idx = rng.permutation(len(neg_ids))
    train_pos = pos_idx[: args.n_train_pos]; heldout_pos = pos_idx[args.n_train_pos :]
    train_neg = neg_idx[: args.n_train_neg]; heldout_neg = neg_idx[args.n_train_neg :]

    Xtr = np.concatenate([pos_emb[train_pos], neg_emb[train_neg]])
    ytr = np.concatenate([np.ones(len(train_pos)), np.zeros(len(train_neg))])
    Xho = np.concatenate([pos_emb[heldout_pos], neg_emb[heldout_neg]])
    yho = np.concatenate([np.ones(len(heldout_pos)), np.zeros(len(heldout_neg))])

    print("[4/5] training LogisticRegression gate ...")
    clf = LogisticRegression(C=1.0, max_iter=2000)
    clf.fit(Xtr, ytr)
    p_ho = clf.predict_proba(Xho)[:, 1]
    auroc_bal = roc_auc_score(yho, p_ho)
    auprc_bal = average_precision_score(yho, p_ho)
    print(f"  balanced held-out AUROC={auroc_bal:.4f} AUPRC={auprc_bal:.4f}")
    joblib.dump(clf, f"{args.out_dir}/cazyme_gate_model.joblib")

    print("[5/5] scoring the realistic-imbalance slice (never seen in training) ...")
    import csv
    truth_rows = list(csv.DictReader(open(args.realistic_slice_tsv), delimiter="\t"))
    truth_by_id = {r["protein_id"]: r["truth"] for r in truth_rows}

    rs_d = np.load(args.realistic_slice, allow_pickle=True)
    rs_ids, rs_emb = rs_d["ids"], l2norm(rs_d["emb"].astype(np.float32))
    rs_truth = np.array([truth_by_id.get(i, "unknown") for i in rs_ids])
    p_slice = clf.predict_proba(rs_emb)[:, 1]

    # clean cazyme-vs-non subset (exclude gray, no ground truth there)
    clean_mask = np.isin(rs_truth, ["cazyme", "non_cazyme"])
    y_clean = (rs_truth[clean_mask] == "cazyme").astype(int)
    p_clean = p_slice[clean_mask]
    auroc_real = roc_auc_score(y_clean, p_clean)
    auprc_real = average_precision_score(y_clean, p_clean)
    base_rate = y_clean.mean()
    print(f"  realistic slice (n={clean_mask.sum()}, base_rate={base_rate:.4f}): "
          f"AUROC={auroc_real:.4f} AUPRC={auprc_real:.4f}")

    precision, recall, thresholds = precision_recall_curve(y_clean, p_clean)
    np.savez(f"{args.out_dir}/gate_pr_curve.npz",
             precision=precision, recall=recall, thresholds=thresholds,
             p_slice=p_slice, rs_truth=rs_truth, y_clean=y_clean, p_clean=p_clean)

    # operating points
    ops = []
    for t in [0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.994]:
        flagged = p_clean >= t
        tp = (flagged & (y_clean == 1)).sum(); fp = (flagged & (y_clean == 0)).sum()
        fn = ((~flagged) & (y_clean == 1)).sum()
        prec = tp / max(tp + fp, 1); rec = tp / max(tp + fn, 1)
        f1 = 2 * prec * rec / max(prec + rec, 1e-9)
        ops.append((t, int(flagged.sum()), prec, rec, f1))
    with open(f"{args.out_dir}/gate_operating_points.tsv", "w") as f:
        f.write("threshold\tn_flagged\tprecision\trecall\tf1\n")
        for t, n, p, r, f1v in ops:
            f.write(f"{t}\t{n}\t{p:.4f}\t{r:.4f}\t{f1v:.4f}\n")

    # single consistent abstention threshold: retire non-CAZymes at >=99.99% purity
    t_lo = None
    for t in np.linspace(0.01, 0.5, 200):
        below = p_clean < t
        if below.sum() == 0:
            continue
        purity = (y_clean[below] == 0).mean()
        if purity >= 0.9999:
            t_lo = t
    if t_lo is None:
        t_lo = 0.10
    below = p_clean < t_lo
    coverage = (y_clean[below] == 0).sum() / max((y_clean == 0).sum(), 1)
    purity = (y_clean[below] == 0).mean() if below.sum() else float("nan")
    lost = int((y_clean[below] == 1).sum())

    gray_mask = rs_truth == "gray"
    gray_scores = p_slice[gray_mask]
    gray_below = (gray_scores < t_lo).sum()

    calib = {
        "t_lo": float(t_lo),
        "coverage_non_cazyme_retired": float(coverage),
        "purity_of_retired": float(purity),
        "true_cazyme_lost": lost,
        "auroc_balanced": float(auroc_bal), "auprc_balanced": float(auprc_bal),
        "auroc_realistic": float(auroc_real), "auprc_realistic": float(auprc_real),
        "realistic_base_rate": float(base_rate),
        "tier_means": {
            "cazyme": float(p_slice[rs_truth == "cazyme"].mean()),
            "gray": float(p_slice[rs_truth == "gray"].mean()),
            "non_cazyme": float(p_slice[rs_truth == "non_cazyme"].mean()),
        },
    }
    json.dump(calib, open(f"{args.out_dir}/gate_calibration.json", "w"), indent=2)

    adjud = {
        "n_gray": int(gray_mask.sum()),
        "gray_retired_as_non_at_t_lo": int(gray_below),
        "gray_retired_frac": float(gray_below / max(gray_mask.sum(), 1)),
        "t_lo": float(t_lo),
        "note": "gray-zone proteins have no ground truth; this is a triage suggestion, not a validated call.",
    }
    json.dump(adjud, open(f"{args.out_dir}/gate_grayzone_adjudication.json", "w"), indent=2)

    with open(f"{args.out_dir}/gate_eval_realistic_slice.tsv", "w") as f:
        f.write("protein_id\ttruth\tgate_score\n")
        for pid, truth, score in zip(rs_ids, rs_truth, p_slice):
            f.write(f"{pid}\t{truth}\t{score:.6f}\n")

    print(f"done. wrote outputs to {args.out_dir}/")


if __name__ == "__main__":
    main()
