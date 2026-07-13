#!/usr/bin/env python
"""
run_example.py -- reproduce the dbCAN4-advanced reasoning-track prototype outputs.

This script regenerates the six numeric/tabular deliverables produced during the reasoning-track
prototype build, using ONLY the staged input files listed below, plus the three library modules
`triage_rule.py`, `conformal.py`, `ood_novelty.py` (must be importable from the same directory or
on PYTHONPATH).

INPUTS (all read from --assets DIR, default '.'):
  Eval slice (n=4726 fungal CAZymes, used for triage/conformal/OOD calibration + evaluation):
    - head_eval_pred.tsv        classifier + trained contrastive-centroid/kNN head predictions
    - esmc_retrieval_pred.tsv   off-the-shelf ESM-C kNN + centroid retrieval predictions
    - diamond_eval2025_pred.tsv DIAMOND sequence-homology baseline predictions
  3 worked-example proteins (267317 / 602276 / 169208):
    - raw_knn.tsv, raw_centroid.tsv, raw_contrastive.tsv   per-head predictions
    - fusion_raw.tsv                                       fusion-layer final calls + vote agreement

  NOT used by this script (present in the staged bundle but not needed for these six outputs):
    - eval2025_overview.tsv, head_metrics.json, esmc_retrieval_summary.json: these carry the
      dbCAN3-baseline overview table and the PROJECT'S OWN previously reported baseline metrics
      (novelty_detection_auroc.centroid_margin = 0.6549 etc.). This script independently
      RE-COMPUTES that 0.6549 baseline number directly from esmc_retrieval_pred.tsv's cent_margin
      column rather than reading it out of esmc_retrieval_summary.json, and prints both for
      cross-check.
    - domains.tsv, ec_prediction.tsv, deeptmhmm.tsv, localization.tsv, physicochem.tsv,
      structures.tsv: these fed the GROUNDED PER-PROTEIN REASONING REPORTS (report_267317.md /
      report_602276.md / report_169208.md), which were drafted via an LLM call (host.llm) under a
      strict grounding constraint. That step is NOT reproduced by this script -- it is not a
      deterministic computation, and re-running it would call an LLM and could legitimately return
      different prose on a re-run even with unchanged inputs. This script reproduces only the
      DETERMINISTIC, code-only pieces: the triage rule, conformal sets, and OOD/novelty scoring.

OUTPUTS (written to --outdir, default '.'):
  - triage_eval_slice.csv              per-protein triage score/tier on the full eval slice
  - conformal_calibration_results.csv  coverage/set-size at 4 target alphas, calibration/test split
  - conformal_demo_predictions.json    conformal prediction sets for the 3 example proteins
  - ood_novelty_results.csv            baseline vs. energy vs. logistic-regression AUROC summary
  - ood_eval_slice_scored.csv          per-protein energy + logistic novelty scores on eval slice
  - ood_demo_scores.json               energy novelty scores for the 3 example proteins

Run:
    python run_example.py --assets /path/to/staged_bundle --outdir ./out
"""
import argparse
import json
import os
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import MinMaxScaler, StandardScaler

from triage_rule import compute_triage_score
from conformal import split_conformal_calibrate, build_prediction_set
from ood_novelty import (
    ENERGY_DIRECTIONS, LOGISTIC_FEATURES,
    fit_energy_scaler, energy_novelty_score,
)

PROTEINS = ['267317', '602276', '169208']
TRUE_FAMILY = {
    '267317': 'GH78 (multidomain GH28+GH78)',
    '602276': 'GH11',
    '169208': 'GH183',
}
RANDOM_SEED = 42
HEADS_CFG = [
    ('classifier', 'clf_pred', 'clf_conf'),
    ('ESM-C-kNN', 'knn_pred', 'knn_conf'),
    ('ESM-C-centroid', 'cent_pred', 'cent_conf'),
    ('ESM-C-contrastive-kNN', 'contr_knn_pred', 'contr_knn_purity'),
    ('DIAMOND', 'diamond_pred', 'diamond_conf'),
]


def load_eval_slice(assets_dir):
    """Load + merge the 3 eval-slice prediction tables into one n=4726 dataframe."""
    head_eval = pd.read_csv(os.path.join(assets_dir, 'head_eval_pred.tsv'), sep='\t')
    esmc_eval = pd.read_csv(os.path.join(assets_dir, 'esmc_retrieval_pred.tsv'), sep='\t')
    diamond_eval = pd.read_csv(os.path.join(assets_dir, 'diamond_eval2025_pred.tsv'), sep='\t')

    merged = head_eval.merge(esmc_eval, on=['query_id', 'novelty', 'true_families'],
                              suffixes=('_clf', '_esmc'))
    full = merged.merge(diamond_eval[['query_id', 'pred_families', 'top_pident']], on='query_id')
    full['diamond_conf'] = full['top_pident'].fillna(0) / 100.0
    full['diamond_pred'] = full['pred_families'].fillna('NO_HIT').apply(lambda s: s.split(',')[0])
    return full


def true_label_set(true_families_str):
    return set(true_families_str.split(','))


def build_menu(row, heads_cfg=HEADS_CFG):
    menu = {}
    for _, pred_col, conf_col in heads_cfg:
        fam, conf = row[pred_col], row[conf_col]
        if fam not in menu or conf > menu[fam]:
            menu[fam] = conf
    return menu


def nonconformity_for_label(row, y, heads_cfg=HEADS_CFG):
    best_conf, found = 0.0, False
    for _, pred_col, conf_col in heads_cfg:
        if row[pred_col] == y:
            found = True
            best_conf = max(best_conf, row[conf_col])
    return 1 - best_conf if found else 1.0


def true_label_nc_score(row):
    trues = true_label_set(row['true_families'])
    return min(nonconformity_for_label(row, y) for y in trues)


def add_menu_and_conformal_columns(full):
    full = full.copy()
    full['menu'] = full.apply(build_menu, axis=1)
    full['menu_size'] = full['menu'].apply(len)
    full['true_in_menu'] = full.apply(
        lambda r: len(true_label_set(r['true_families']) & set(r['menu'].keys())) > 0, axis=1)
    full['true_nc_score'] = full.apply(true_label_nc_score, axis=1)
    return full


# ---------------------------------------------------------------------------
# 1. Disagreement-detection triage, evaluated on the eval slice
# ---------------------------------------------------------------------------

def eval_triage_row(row):
    calls = {'clf': row['clf_pred'], 'cent': row['contr_cent_pred'], 'knn': row['contr_knn_pred']}
    confs = {'clf': row['clf_conf'], 'cent': row['contr_cent_margin'], 'knn': row['contr_knn_purity']}
    vote_counts = Counter(calls.values())
    _, majority_n = vote_counts.most_common(1)[0]
    agreement_frac = majority_n / 3.0
    max_head_conf = max(confs.values())
    fusion_like_conf = row['clf_conf']
    conf_collapse = max(0.0, max_head_conf - fusion_like_conf)
    score = 0.40 * (1 - agreement_frac) + 0.35 * (1 - fusion_like_conf) + 0.25 * min(conf_collapse, 1.0)
    tier = "FLAG" if score > 0.35 else ("WATCH" if score > 0.15 else "ACCEPT")
    return pd.Series({'triage_score_v2': score, 'triage_tier_v2': tier})


def run_triage(full_eval, outdir):
    scored = full_eval.join(full_eval.apply(eval_triage_row, axis=1))
    scored['exact_clf'] = scored.apply(lambda r: r['clf_pred'] in r['true_families'].split(','), axis=1)
    # Write only the original triage-relevant columns (the menu/conformal columns computed upstream in
    # the shared full_eval pipeline are a different deliverable and are dropped here for a byte-faithful
    # reproduction of triage_eval_slice.csv).
    triage_cols = ['query_id', 'novelty', 'true_families', 'clf_pred', 'clf_conf', 'contr_cent_pred',
                   'contr_cent_margin', 'contr_knn_pred', 'contr_knn_purity', 'knn_pred', 'knn_conf',
                   'knn_purity', 'knn_margin', 'cent_pred', 'cent_conf', 'cent_margin',
                   'triage_score_v2', 'triage_tier_v2', 'exact_clf']
    scored[triage_cols].to_csv(os.path.join(outdir, 'triage_eval_slice.csv'), index=False)

    y_wrong = (~scored['exact_clf']).astype(int)
    auc = roc_auc_score(y_wrong, scored['triage_score_v2'])
    summary = scored.groupby('triage_tier_v2').agg(n=('exact_clf', 'size'),
                                                     accuracy=('exact_clf', 'mean')).reset_index()
    return scored, auc, summary


def demo_triage(evidence_bundles):
    rows = []
    for pid in PROTEINS:
        b = evidence_bundles[pid]
        k, c, ct, fu = b['knn'], b['centroid'], b['contrastive'], b['fusion']
        head_confs = {'kNN': k['knn_conf'], 'centroid': c['cent_conf'], 'contrastive_clf': ct['clf_conf']}
        result = compute_triage_score(head_confs, fusion_conf=fu['confidence'],
                                       agreement=fu['agreement'], n_total_votes=4)
        result['protein_id'] = pid
        result['true_family'] = TRUE_FAMILY[pid]
        result['fusion_call'] = fu['family']
        result['fusion_conf'] = fu['confidence']
        rows.append(result)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 2. Split-conformal prediction sets
# ---------------------------------------------------------------------------

def run_conformal(full_eval, outdir):
    n = len(full_eval)
    np.random.seed(RANDOM_SEED)
    idx = np.random.permutation(n)
    n_calib = n // 2
    calib = full_eval.iloc[idx[:n_calib]].reset_index(drop=True)
    test = full_eval.iloc[idx[n_calib:]].reset_index(drop=True)

    rows = []
    for alpha, label in [(0.30, '70%'), (0.20, '80%'),
                          (0.10, '90% (target, see ceiling note)'),
                          (0.05, '95% (target, see ceiling note)')]:
        q_hat = split_conformal_calibrate(calib, alpha)
        pred_sets = test.apply(lambda r: build_prediction_set(r, q_hat), axis=1)
        covered = [len(true_label_set(t) & set(ps)) > 0
                   for t, ps in zip(test['true_families'], pred_sets)]
        sizes = pred_sets.apply(len)
        rows.append({
            'alpha': alpha, 'target_coverage_label': label, 'q_hat': round(q_hat, 4),
            'empirical_coverage_test': round(float(np.mean(covered)), 4),
            'mean_set_size': round(sizes.mean(), 3),
            'pct_size_0_abstain': round((sizes == 0).mean(), 4),
            'pct_size_1': round((sizes == 1).mean(), 4),
            'pct_size_ge2': round((sizes >= 2).mean(), 4),
        })
    report_table = pd.DataFrame(rows)
    report_table.to_csv(os.path.join(outdir, 'conformal_calibration_results.csv'), index=False)

    ceiling_test = test['true_in_menu'].mean()
    # q_hat values re-derived on the FULL eval slice (used for out-of-sample demo-protein application)
    q_hat_70 = split_conformal_calibrate(full_eval, 0.30)
    q_hat_80 = split_conformal_calibrate(full_eval, 0.20)
    q_hat_90 = split_conformal_calibrate(full_eval, 0.10)
    return report_table, ceiling_test, (q_hat_70, q_hat_80, q_hat_90)


def demo_menu(evidence_bundles, pid):
    b = evidence_bundles[pid]
    k, c, ct = b['knn'], b['centroid'], b['contrastive']
    menu = {}
    for fam, conf in [(k['knn_pred'], k['knn_conf']), (c['cent_pred'], c['cent_conf']),
                      (ct['clf_pred'], ct['clf_conf']), (ct['contr_knn_pred'], ct['contr_knn_purity'])]:
        if fam not in menu or conf > menu[fam]:
            menu[fam] = conf
    return menu


def demo_conformal(evidence_bundles, q_hats, outdir):
    q_hat_70, q_hat_80, q_hat_90 = q_hats
    export = []
    for pid in PROTEINS:
        menu = demo_menu(evidence_bundles, pid)
        set70 = [f for f, c in menu.items() if (1 - c) <= q_hat_70]
        set80 = [f for f, c in menu.items() if (1 - c) <= q_hat_80]
        set90 = [f for f, c in menu.items() if (1 - c) <= q_hat_90]
        export.append({
            'protein_id': pid, 'true_family': TRUE_FAMILY[pid], 'menu': menu,
            'pred_set_70pct': set70, 'pred_set_80pct': set80, 'pred_set_90pct_saturated': set90,
        })
    with open(os.path.join(outdir, 'conformal_demo_predictions.json'), 'w') as f:
        json.dump(export, f, indent=2)
    return export


# ---------------------------------------------------------------------------
# 3. OOD / novelty scoring
# ---------------------------------------------------------------------------

def run_ood(full_eval, outdir):
    y = (full_eval['novelty'] == 'novel_family').astype(int).values

    # baseline: off-the-shelf centroid_margin (lower margin => more novel)
    auc_baseline = roc_auc_score(y, -full_eval['cent_margin'].values)

    # unsupervised energy score (calibration-only scaler fit, using the same 50/50 split as conformal)
    np.random.seed(RANDOM_SEED)
    n = len(full_eval)
    idx = np.random.permutation(n)
    n_calib = n // 2
    calib_idx, test_idx = idx[:n_calib], idx[n_calib:]

    energy_scaler = fit_energy_scaler(full_eval.iloc[calib_idx])
    energy_scores_full = energy_novelty_score(full_eval, energy_scaler)
    auc_energy_full = roc_auc_score(y, energy_scores_full)

    # supervised logistic score, evaluated by 5-fold CV
    X = full_eval[LOGISTIC_FEATURES].values
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    cv_aucs, baseline_cv = [], []
    for tr_idx, te_idx in skf.split(X, y):
        scaler = StandardScaler().fit(X[tr_idx])
        clf = LogisticRegression(max_iter=1000, class_weight='balanced')
        clf.fit(scaler.transform(X[tr_idx]), y[tr_idx])
        p = clf.predict_proba(scaler.transform(X[te_idx]))[:, 1]
        cv_aucs.append(roc_auc_score(y[te_idx], p))
        baseline_cv.append(roc_auc_score(y[te_idx], -full_eval['cent_margin'].values[te_idx]))

    # full-data fit, used only to SCORE every protein for the per-protein CSV (not for the CV AUROC claim)
    scaler_full = StandardScaler().fit(X)
    clf_full = LogisticRegression(max_iter=1000, class_weight='balanced').fit(scaler_full.transform(X), y)
    logistic_scores_full = clf_full.predict_proba(scaler_full.transform(X))[:, 1]

    results = pd.DataFrame([
        {'method': 'Baseline: off-the-shelf centroid_margin (current pipeline)',
         'auroc_reported': 0.6549, 'auroc_reproduced_here': round(auc_baseline, 4),
         'evaluation': 'full eval slice n={} (as in esmc_retrieval_summary.json)'.format(len(full_eval))},
        {'method': 'Energy score (unsupervised, 5-signal combination)',
         'auroc_reported': None, 'auroc_reproduced_here': round(auc_energy_full, 4),
         'evaluation': 'full eval slice n={}, calibration-only scaler fit'.format(len(full_eval))},
        {'method': 'Logistic regression (supervised, 10-signal combination)',
         'auroc_reported': None,
         'auroc_reproduced_here': round(float(np.mean(cv_aucs)), 4),
         'evaluation': '5-fold CV mean +/- {:.4f}, n={}'.format(np.std(cv_aucs), len(full_eval))},
    ])
    results.to_csv(os.path.join(outdir, 'ood_novelty_results.csv'), index=False)

    scored = full_eval[['query_id', 'novelty', 'true_families']].copy()
    scored['energy_novelty_score'] = energy_scores_full
    scored['logistic_novelty_score_infold'] = logistic_scores_full
    scored.to_csv(os.path.join(outdir, 'ood_eval_slice_scored.csv'), index=False)

    return results, auc_baseline, auc_energy_full, cv_aucs, baseline_cv


def run_demo_ood(evidence_bundles, full_eval, outdir):
    demo_directions = {'cent_margin': -1, 'cent_conf': +1, 'clf_conf': -1}
    np.random.seed(RANDOM_SEED)
    n = len(full_eval)
    idx = np.random.permutation(n)
    calib_idx = idx[:n // 2]
    mm = MinMaxScaler().fit(full_eval[list(demo_directions.keys())].values[calib_idx])

    energy_scaler_full = fit_energy_scaler(full_eval.iloc[calib_idx])
    energy_full_scores = energy_novelty_score(full_eval, energy_scaler_full)

    export = []
    for pid in PROTEINS:
        b = evidence_bundles[pid]
        vals = np.array([[b['centroid']['cent_margin'], b['centroid']['cent_conf'],
                           b['contrastive']['clf_conf']]])
        scaled = mm.transform(vals)[0]
        e = 0.0
        for i, (sig, d) in enumerate(demo_directions.items()):
            e += scaled[i] if d == 1 else (1 - scaled[i])
        e /= len(demo_directions)
        pct = float((energy_full_scores < e).mean() * 100)
        export.append({'protein_id': pid, 'true_family': TRUE_FAMILY[pid],
                        'energy_novelty_score': round(float(e), 4),
                        'percentile_vs_eval_slice': round(pct, 1)})
    with open(os.path.join(outdir, 'ood_demo_scores.json'), 'w') as f:
        json.dump(export, f, indent=2)
    return export


# ---------------------------------------------------------------------------
# Evidence bundles for the 3 demo proteins (triage + conformal + OOD demo application)
# ---------------------------------------------------------------------------

def load_evidence_bundles(assets_dir):
    def get_row(df, pid, col='query_id'):
        r = df[df[col].astype(str) == str(pid)]
        return r.iloc[0].to_dict() if len(r) else None

    knn = pd.read_csv(os.path.join(assets_dir, 'raw_knn.tsv'), sep='\t')
    cent = pd.read_csv(os.path.join(assets_dir, 'raw_centroid.tsv'), sep='\t')
    contr = pd.read_csv(os.path.join(assets_dir, 'raw_contrastive.tsv'), sep='\t')
    fusion = pd.read_csv(os.path.join(assets_dir, 'fusion_raw.tsv'), sep='\t')

    bundles = {}
    for pid in PROTEINS:
        pid_int = int(pid)
        bundles[pid] = {
            'knn': get_row(knn, pid_int),
            'centroid': get_row(cent, pid_int),
            'contrastive': get_row(contr, pid_int),
            'fusion': get_row(fusion, pid_int, col='protein_id'),
        }
    return bundles


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--assets', default='.', help='Directory containing the staged input files')
    parser.add_argument('--outdir', default='.', help='Directory to write the regenerated CSV/JSON files')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("=" * 78)
    print("dbCAN4-advanced reasoning-track prototype: reproduction run")
    print("=" * 78)

    full_eval = load_eval_slice(args.assets)
    full_eval = add_menu_and_conformal_columns(full_eval)
    evidence_bundles = load_evidence_bundles(args.assets)
    print(f"\nLoaded eval slice: n={len(full_eval)} "
          f"({(full_eval['novelty']=='novel_seq').sum()} novel_seq / "
          f"{(full_eval['novelty']=='novel_family').sum()} novel_family)")

    # --- 1. Triage ---
    scored, triage_auc, triage_summary = run_triage(full_eval, args.outdir)
    demo_triage_df = demo_triage(evidence_bundles)
    print("\n--- Disagreement-detection triage ---")
    print(f"AUROC (triage_score for flagging wrong classifier calls, n={len(scored)}): {triage_auc:.4f}")
    print(triage_summary.to_string(index=False))
    print("\nDemo proteins:")
    print(demo_triage_df[['protein_id', 'true_family', 'fusion_call', 'fusion_conf',
                           'triage_score', 'triage_tier']].to_string(index=False))

    # --- 2. Conformal ---
    conformal_table, ceiling_test, q_hats = run_conformal(full_eval, args.outdir)
    demo_conformal_export = demo_conformal(evidence_bundles, q_hats, args.outdir)
    print("\n--- Conformal prediction sets ---")
    print(conformal_table.to_string(index=False))
    print(f"Coverage ceiling on test split (true label reachable by some head): {ceiling_test:.4f}")
    print("\nDemo proteins (menu / 70% / 80% / 90%-saturated sets):")
    for row in demo_conformal_export:
        print(f"  {row['protein_id']} (true {row['true_family']}): menu={row['menu']}")
        print(f"    70%: {row['pred_set_70pct']}  80%: {row['pred_set_80pct']}  "
              f"90%(sat): {row['pred_set_90pct_saturated']}")

    # --- 3. OOD / novelty ---
    ood_results, auc_baseline, auc_energy, cv_aucs, baseline_cv = run_ood(full_eval, args.outdir)
    demo_ood_export = run_demo_ood(evidence_bundles, full_eval, args.outdir)
    print("\n--- OOD / novelty scoring ---")
    print(ood_results.to_string(index=False))
    print(f"\n5-fold CV logistic AUROC per fold: {[round(a,4) for a in cv_aucs]}")
    print(f"5-fold CV baseline (cent_margin) AUROC per fold: {[round(a,4) for a in baseline_cv]}")
    print("\nDemo proteins (energy novelty score + percentile vs. eval slice):")
    for row in demo_ood_export:
        print(f"  {row['protein_id']} (true {row['true_family']}): "
              f"score={row['energy_novelty_score']:.4f} "
              f"({row['percentile_vs_eval_slice']:.1f}th percentile)")

    print("\n" + "=" * 78)
    print("Outputs written to:", os.path.abspath(args.outdir))
    for fname in ['triage_eval_slice.csv', 'conformal_calibration_results.csv',
                  'conformal_demo_predictions.json', 'ood_novelty_results.csv',
                  'ood_eval_slice_scored.csv', 'ood_demo_scores.json']:
        print(" -", fname)
    print("=" * 78)


if __name__ == '__main__':
    main()
