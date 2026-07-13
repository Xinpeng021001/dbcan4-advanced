"""
dbCAN4-advanced split-conformal prediction sets over the fusion-tier candidate menu.

Design
------
For each protein, we build a candidate "menu" = union of family predictions across the available
per-method heads (classifier, ESM-C-kNN, ESM-C-centroid, ESM-C-contrastive-kNN, DIAMOND), each
carrying its own reported confidence in [0,1] (softmax/purity/margin-derived; see the individual
head TSVs). For a candidate family y with confidence c(y) in a protein's menu, the "nonconformity
score" is s(x,y) = 1 - c(y); if y is not in the menu at all, s(x,y) = 1 (max nonconformity).

Split-conformal calibration (standard Vovk et al. procedure):
  1. On a held-out calibration set, compute the nonconformity score of the TRUE label for every
     calibration example: s_i = min over true label(s) of s(x_i, y_true).
  2. q_hat = the ceil((n+1)(1-alpha))/n empirical quantile of {s_i}, capped at 1.0.
  3. For a new/test protein, the prediction SET = { y in menu(x) : s(x,y) <= q_hat }.

This gives the standard split-conformal marginal coverage guarantee **PROVIDED the true label is
representable at all** -- i.e. provided the true family appears in SOME head's menu. If it does not
(entirely novel family relative to what any head proposes), no achievable q_hat can cover that
example, and the guaranteed coverage saturates below 100% at the "menu coverage ceiling"
P(true label in menu). This is a real, reported limitation below, not smoothed over.

Usage
-----
    from conformal import split_conformal_calibrate, build_prediction_set

    q_hat = split_conformal_calibrate(calibration_df, alpha=0.20)   # target 80% coverage
    pred_set = build_prediction_set(row, q_hat)   # row must carry a 'menu' dict {family: confidence}
"""
import numpy as np

def split_conformal_calibrate(calib_df, alpha):
    """calib_df must have a 'true_nc_score' column: nonconformity of the true label per example."""
    n_c = len(calib_df)
    q_level = min(np.ceil((n_c + 1) * (1 - alpha)) / n_c, 1.0)
    return np.quantile(calib_df['true_nc_score'], q_level, method='higher')

def build_prediction_set(row, q_hat):
    """row must have a 'menu' column: {family: confidence} across all available per-method heads."""
    return [fam for fam, conf in row['menu'].items() if (1 - conf) <= q_hat]
