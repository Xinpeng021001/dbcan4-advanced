"""
build_visualizations.py — dbCAN4-advanced viz-track deliverables.

Regenerates the three interactive HTML explorers and two static PNG companions
built for the ESM-C protein-language-model tier of dbCAN4-advanced, from the
pre-staged viz_track_assets bundle. No network access and no repo/live-compute
dependency — every input is read from local files.

Input files (all expected under --assets DIR, default '.'):
    embedding_umap_coords.json   - reference + query UMAP coords (raw ESM-C and
                                    trained-head projection spaces)
    head_metrics.json            - trained-head eval metrics (exact/overlap
                                    accuracy per scheme, novelty AUROC)
    esmc_retrieval_summary.json  - off-the-shelf (untrained) ESM-C retrieval
                                    baseline metrics
    train_heads.log              - real per-epoch training log for the
                                    contrastive/classifier head (7 logged
                                    epochs: 0, 5, 10, 15, 20, 25, 29)
    head_eval_pred.tsv           - per-protein eval predictions vs. truth
                                    (4,726 rows): query_id, novelty,
                                    true_families, clf_pred, clf_conf,
                                    contr_cent_pred, contr_cent_margin,
                                    contr_knn_pred, contr_knn_purity

Outputs (written under --outdir DIR, default '.'):
    embedding_explorer.html               - Deliverable 1: reference embedding
                                             space, raw-ESM-C vs trained-head
                                             UMAP toggle, family coloring +
                                             query-protein overlay
    training_dashboard.html               - Deliverable 2: loss/val-accuracy
                                             curve, per-family accuracy,
                                             margin/purity histograms
    training_summary.png                  - static companion to deliverable 2
                                             (figure-style rules applied)
    calibration_confusion_explorer.html   - Deliverable 3: reliability diagram
                                             + clickable family x family
                                             confusion matrix
    reliability_diagram.png               - static companion to deliverable 3

Usage:
    python build_visualizations.py --assets /path/to/viz_track_assets --outdir ./out

Requirements: see requirements.txt (plotly, pandas, numpy, matplotlib).
"""
import argparse
import json
import os
import re

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.offline as pyo

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# Static-figure styling (publication-grade rcParams; standalone reimplementation
# of the `figure-style` skill's apply_figure_style(), inlined here so this
# script has no dependency outside the packages in requirements.txt).
# ----------------------------------------------------------------------------

def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    """Set matplotlib rcParams for publication-grade static output.

    frame : 'open' (bottom+left spines, default) | 'boxed' (all four) | 'none'
    font  : sans-serif family name; None = system default sans-serif
    sizes : (base, secondary, tick) - titles/axis-labels, legend/annotation, ticks
    grid  : whether to draw axes.grid (default False)
    """
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    matplotlib.rcParams.update(rc)


# ----------------------------------------------------------------------------
# Shared loaders
# ----------------------------------------------------------------------------

def load_bundle(assets_dir):
    """Load all five staged input files into memory."""
    with open(os.path.join(assets_dir, "embedding_umap_coords.json")) as f:
        emb = json.load(f)
    with open(os.path.join(assets_dir, "head_metrics.json")) as f:
        hm = json.load(f)
    with open(os.path.join(assets_dir, "esmc_retrieval_summary.json")) as f:
        rs = json.load(f)
    with open(os.path.join(assets_dir, "train_heads.log")) as f:
        log_lines = f.readlines()
    pred = pd.read_csv(os.path.join(assets_dir, "head_eval_pred.tsv"), sep="\t")
    return emb, hm, rs, log_lines, pred


def parse_train_log(log_lines):
    """Parse the real per-epoch training log. Returns (epochs, losses, accs).

    The log contains exactly 7 logged epochs (0, 5, 10, 15, 20, 25, 29) out of
    30 total training epochs -- this is the only per-epoch record that exists;
    do not interpolate a finer curve.
    """
    epochs, losses, accs = [], [], []
    pat = re.compile(r"epoch\s+(\d+)\s+loss\s+([\d.]+)\s+val_clf_acc\s+([\d.]+)")
    for line in log_lines:
        m = pat.search(line)
        if m:
            epochs.append(int(m.group(1)))
            losses.append(float(m.group(2)))
            accs.append(float(m.group(3)))
    return epochs, losses, accs


def add_correctness_columns(pred):
    """Derive overlap-match correctness booleans and the primary true-family
    token (first item of a possibly multi-label, comma-separated string)."""
    pred = pred.copy()

    def token_set(s):
        return set(s.split(","))

    pred["true_set"] = pred["true_families"].apply(token_set)
    pred["true_primary"] = pred["true_families"].apply(lambda s: s.split(",")[0])
    pred["clf_correct_overlap"] = pred.apply(lambda r: r["clf_pred"] in r["true_set"], axis=1)
    pred["clf_correct_exact"] = pred["clf_pred"] == pred["true_families"]
    pred["knn_correct_overlap"] = pred.apply(lambda r: r["contr_knn_pred"] in r["true_set"], axis=1)
    pred["cent_correct_overlap"] = pred.apply(lambda r: r["contr_cent_pred"] in r["true_set"], axis=1)
    return pred


# ----------------------------------------------------------------------------
# Deliverable 1: interactive embedding explorer
# ----------------------------------------------------------------------------

def build_embedding_explorer(emb, outpath, top_n=24):
    """Reference-space scatter (9,815 proteins / 814 families), colored by
    family (top-N distinct, rest grouped as 'other'), with a button toggle
    between the raw-ESM-C UMAP and the trained-head UMAP, and the 3 named
    query proteins overlaid as star markers.
    """
    df = pd.DataFrame({
        "id": emb["ref_ids"],
        "fam": emb["ref_fams"],
        "raw_x": emb["raw_umap_x"],
        "raw_y": emb["raw_umap_y"],
        "trained_x": emb["trained_umap_x"],
        "trained_y": emb["trained_umap_y"],
    })

    fam_counts = df["fam"].value_counts()
    top_fams = fam_counts.head(top_n).index.tolist()
    df["fam_group"] = df["fam"].apply(lambda f: f if f in top_fams else "other")

    palette = px.colors.qualitative.Dark24
    color_map = {f: palette[i % len(palette)] for i, f in enumerate(top_fams)}
    color_map["other"] = "#cccccc"

    plot_order = ["other"] + top_fams

    fig = go.Figure()
    for view, xcol, ycol in [("raw", "raw_x", "raw_y"), ("trained", "trained_x", "trained_y")]:
        for fam in plot_order:
            sub = df[df["fam_group"] == fam]
            is_other = (fam == "other")
            fig.add_trace(go.Scattergl(
                x=sub[xcol], y=sub[ycol],
                mode="markers",
                name=fam,
                legendgroup=fam,
                marker=dict(size=4 if is_other else 6, color=color_map[fam],
                            opacity=0.35 if is_other else 0.85, line=dict(width=0)),
                text=[f"id: {i}<br>family: {f}" for i, f in zip(sub["id"], sub["fam"])],
                hoverinfo="text",
                visible=(view == "raw"),
                showlegend=(view == "raw"),
            ))

    n_fam_traces_per_view = len(plot_order)

    query_ids = emb["query_ids"]
    query_fams = emb["query_fams"]
    for view, xkey, ykey in [
        ("raw", "query_raw_umap_x", "query_raw_umap_y"),
        ("trained", "query_trained_umap_x", "query_trained_umap_y"),
    ]:
        qx, qy = emb[xkey], emb[ykey]
        fig.add_trace(go.Scatter(
            x=qx, y=qy,
            mode="markers+text",
            name="query proteins",
            legendgroup="query",
            marker=dict(size=16, color="black", symbol="star", line=dict(width=2, color="white")),
            text=[f"Q{i + 1}" for i in range(len(query_ids))],
            textposition="top center",
            hovertext=[f"QUERY id: {qi}<br>family: {qf}" for qi, qf in zip(query_ids, query_fams)],
            hoverinfo="text",
            visible=(view == "raw"),
            showlegend=(view == "raw"),
        ))

    raw_visible = [True] * n_fam_traces_per_view + [False] * n_fam_traces_per_view + [True, False]
    trained_visible = [False] * n_fam_traces_per_view + [True] * n_fam_traces_per_view + [False, True]
    assert len(raw_visible) == len(fig.data)

    updatemenus = [dict(
        type="buttons", direction="right", x=0.02, y=1.12, xanchor="left", yanchor="top",
        showactive=True,
        buttons=[
            dict(label="Raw ESM-C UMAP", method="update",
                 args=[{"visible": raw_visible},
                       {"title.text": "Reference embedding space — RAW ESM-C (cosine UMAP)"}]),
            dict(label="Trained-head UMAP", method="update",
                 args=[{"visible": trained_visible},
                       {"title.text": "Reference embedding space — TRAINED HEAD (256-dim, L2-norm, cosine UMAP)"}]),
        ],
    )]

    fig.update_layout(
        title=dict(text="Reference embedding space — RAW ESM-C (cosine UMAP)", x=0.5),
        updatemenus=updatemenus,
        width=1150, height=850,
        legend=dict(title=f"Family (top {top_n} shown; rest = 'other')", itemsizing="constant",
                    font=dict(size=9), tracegroupgap=0),
        xaxis_title="UMAP-1", yaxis_title="UMAP-2",
        margin=dict(t=100),
        plot_bgcolor="white",
    )
    fig.update_xaxes(showgrid=True, gridcolor="#eee", zeroline=False)
    fig.update_yaxes(showgrid=True, gridcolor="#eee", zeroline=False)

    meta = emb["meta"]
    annotation_text = (
        f"n={meta['n_subsample']} proteins stratified from {meta['n_families_ref_total']} "
        f"reference families (cap {meta['per_fam_cap']}/family). "
        f"Trained-head proj: {meta['proj_arch']['in_dim']}\u2192{meta['proj_arch']['hidden']}\u2192"
        f"{meta['proj_arch']['out_dim']} dims, L2-normalized, cosine metric."
    )
    fig.add_annotation(text=annotation_text, xref="paper", yref="paper", x=0, y=-0.08, showarrow=False,
                        font=dict(size=10, color="#555"), align="left")

    fig.write_html(outpath, include_plotlyjs=True, full_html=True)


# ----------------------------------------------------------------------------
# Deliverable 2: training dashboard (interactive HTML + static PNG)
# ----------------------------------------------------------------------------

def build_training_dashboard(epochs, losses, accs, pred, outpath, min_support=5):
    """4-panel interactive dashboard: loss/val-acc curve, per-family accuracy
    (families with >= min_support eval proteins), and margin/purity
    histograms split by correctness."""
    from plotly.subplots import make_subplots

    fam_support = pred["true_primary"].value_counts()
    eligible_fams = fam_support[fam_support >= min_support].index
    fam_acc = pred[pred["true_primary"].isin(eligible_fams)].groupby("true_primary").agg(
        n=("clf_correct_overlap", "size"),
        clf_acc_overlap=("clf_correct_overlap", "mean"),
        clf_acc_exact=("clf_correct_exact", "mean"),
        novelty=("novelty", lambda s: s.mode()[0] if s.nunique() == 1 else "mixed"),
    ).reset_index().sort_values("clf_acc_overlap")

    fig2 = make_subplots(
        rows=2, cols=2,
        specs=[[{"secondary_y": True}, {"type": "xy"}],
               [{"type": "xy"}, {"type": "xy"}]],
        subplot_titles=(
            "Training curve: loss & val classifier accuracy (7 logged epochs)",
            f"Per-family classifier accuracy (overlap-match), families with \u2265{min_support} eval support (n={len(fam_acc)})",
            "Contrastive-centroid margin distribution: correct vs incorrect",
            "Contrastive-kNN purity distribution: correct vs incorrect",
        ),
        column_widths=[0.45, 0.55], row_heights=[0.5, 0.5],
        vertical_spacing=0.12, horizontal_spacing=0.09,
    )

    fig2.add_trace(go.Scatter(x=epochs, y=losses, mode="lines+markers", name="train loss",
                               line=dict(color="#d62728", width=2), marker=dict(size=8)),
                   row=1, col=1, secondary_y=False)
    fig2.add_trace(go.Scatter(x=epochs, y=accs, mode="lines+markers", name="val clf accuracy",
                               line=dict(color="#1f77b4", width=2), marker=dict(size=8)),
                   row=1, col=1, secondary_y=True)
    fig2.update_xaxes(title_text="epoch", row=1, col=1)
    fig2.update_yaxes(title_text="loss", row=1, col=1, secondary_y=False, color="#d62728")
    fig2.update_yaxes(title_text="val clf accuracy", row=1, col=1, secondary_y=True, color="#1f77b4", range=[0, 1])

    novelty_colors = {"novel_seq": "#2ca02c", "novel_family": "#ff7f0e", "mixed": "#7f7f7f"}
    for nb in ["novel_seq", "novel_family", "mixed"]:
        sub = fam_acc[fam_acc["novelty"] == nb]
        if len(sub) == 0:
            continue
        fig2.add_trace(go.Bar(
            x=sub["true_primary"], y=sub["clf_acc_overlap"], name=f"novelty: {nb}",
            marker=dict(color=novelty_colors[nb]),
            text=sub["n"], hovertemplate="family=%{x}<br>acc=%{y:.3f}<br>n=%{text}<extra></extra>",
        ), row=1, col=2)
    fig2.update_xaxes(title_text="family (true, primary token)", tickfont=dict(size=6), row=1, col=2)
    fig2.update_yaxes(title_text="clf accuracy (overlap)", range=[0, 1.02], row=1, col=2)

    for correct, label, color in [(True, "correct", "#1f77b4"), (False, "incorrect", "#d62728")]:
        sub = pred[pred["cent_correct_overlap"] == correct]
        fig2.add_trace(go.Histogram(x=sub["contr_cent_margin"], name=f"centroid margin ({label})",
                                     marker=dict(color=color), opacity=0.6, nbinsx=40,
                                     legendgroup=label), row=2, col=1)
    fig2.update_xaxes(title_text="contrastive-centroid margin", row=2, col=1)
    fig2.update_yaxes(title_text="count", row=2, col=1)

    for correct, label, color in [(True, "correct", "#1f77b4"), (False, "incorrect", "#d62728")]:
        sub = pred[pred["knn_correct_overlap"] == correct]
        fig2.add_trace(go.Histogram(x=sub["contr_knn_purity"], name=f"kNN purity ({label})",
                                     marker=dict(color=color), opacity=0.6, nbinsx=40,
                                     legendgroup=label, showlegend=False), row=2, col=2)
    fig2.update_xaxes(title_text="contrastive-kNN purity (k=15 vote fraction)", row=2, col=2)
    fig2.update_yaxes(title_text="count", row=2, col=2)

    fig2.update_layout(
        title=dict(text="dbCAN4-advanced — head training & evaluation dashboard", x=0.5),
        barmode="overlay", width=1500, height=1000,
        legend=dict(font=dict(size=9)), plot_bgcolor="white",
    )
    fig2.update_xaxes(showgrid=True, gridcolor="#eee")
    fig2.update_yaxes(showgrid=True, gridcolor="#eee")

    fig2.write_html(outpath, include_plotlyjs=True, full_html=True)
    return fam_acc


def build_training_summary_png(epochs, losses, accs, hm, rs, outpath):
    """Static figure-style companion: training curve beside an overall
    exact/overlap accuracy comparison (untrained ESM-C retrieval baseline vs
    the three trained-head schemes)."""
    apply_figure_style()

    schemes = ["untrained\nkNN retrieval", "contrastive\nkNN (trained)",
               "contrastive\ncentroid (trained)", "classifier\n(trained)"]
    exact_vals = [rs["schemes"]["knn"]["overall"]["thr=0.0"]["exact_of_all"],
                  hm["contrastive_knn"]["overall"]["exact"],
                  hm["contrastive_centroid"]["overall"]["exact"],
                  hm["classifier"]["overall"]["exact"]]
    overlap_vals = [rs["schemes"]["knn"]["overall"]["thr=0.0"]["overlap_of_all"],
                    hm["contrastive_knn"]["overall"]["overlap"],
                    hm["contrastive_centroid"]["overall"]["overlap"],
                    hm["classifier"]["overall"]["overlap"]]

    x = np.arange(len(schemes))
    w = 0.35

    fig3, axes = plt.subplots(1, 2, figsize=(12, 3.8), constrained_layout=True)

    ax = axes[0]
    ax2 = ax.twinx()
    ax.plot(epochs, losses, marker="o", color="#b03a2e", lw=1.8, label="train loss")
    ax2.plot(epochs, accs, marker="s", color="#1f6f9c", lw=1.8, label="val classifier accuracy")
    ax.set_xlabel("epoch")
    ax.set_ylabel("train loss", color="#b03a2e")
    ax2.set_ylabel("val classifier accuracy", color="#1f6f9c")
    ax2.set_ylim(0, 1.08)
    ax.tick_params(axis="y", colors="#b03a2e")
    ax2.tick_params(axis="y", colors="#1f6f9c")
    ax.set_title("Head training converges by epoch 29 (7 logged points)")
    ax.margins(x=0.08)
    ax.set_xticks(epochs)

    ax_b = axes[1]
    ax_b.bar(x - w / 2, exact_vals, width=w, color="#4c72b0", label="exact-family accuracy")
    ax_b.bar(x + w / 2, overlap_vals, width=w, color="#95b8dd", label="overlap accuracy")
    for xi, v in zip(x - w / 2, exact_vals):
        ax_b.text(xi, v + 0.03, f"{v:.2f}", ha="center", va="bottom", fontsize=6)
    for xi, v in zip(x + w / 2, overlap_vals):
        ax_b.text(xi, v + 0.03, f"{v:.2f}", ha="center", va="bottom", fontsize=6)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(schemes, fontsize=6)
    ax_b.set_ylabel("accuracy (overall, n=4726)")
    ax_b.set_ylim(0, 1.15)
    ax_b.set_xlim(-0.7, 3.7)
    ax_b.legend(frameon=False, loc="upper left", fontsize=6)
    ax_b.set_title("Trained heads beat off-the-shelf ESM-C retrieval overall")

    fig3.savefig(outpath, dpi=300)
    plt.close(fig3)


# ----------------------------------------------------------------------------
# Deliverable 3: calibration + confusion explorer (interactive HTML + static PNG)
# ----------------------------------------------------------------------------

def build_calibration_confusion_explorer(pred, outpath, min_support=15):
    """Reliability diagram (clf_conf vs overlap-match empirical accuracy, 10
    fixed-width bins) plus a clickable family x family confusion matrix
    restricted to true families with >= min_support eval proteins; predictions
    to families below that threshold are bucketed as 'other'. Clicking a cell
    lists the (capped-at-50) misclassified protein ids in a side panel via a
    plain vanilla-JS click handler (no framework)."""
    bins = np.linspace(0, 1, 11)
    pred = pred.copy()
    pred["conf_bin"] = pd.cut(pred["clf_conf"], bins=bins, include_lowest=True)
    rel = pred.groupby("conf_bin", observed=True).agg(
        n=("clf_correct_overlap", "size"),
        mean_conf=("clf_conf", "mean"),
        emp_acc=("clf_correct_overlap", "mean"),
    ).reset_index()
    ece = (rel["n"] * (rel["mean_conf"] - rel["emp_acc"]).abs()).sum() / len(pred)

    fam_support = pred["true_primary"].value_counts()
    top_conf_fams = fam_support[fam_support >= min_support].index.tolist()
    sub = pred[pred["true_primary"].isin(top_conf_fams)].copy()
    sub["pred_bucket"] = sub["clf_pred"].apply(lambda p: p if p in top_conf_fams else "other")

    fam_order_conf = fam_support.loc[top_conf_fams].sort_values(ascending=False).index.tolist()
    col_order = fam_order_conf + ["other"]

    conf_matrix = pd.crosstab(sub["true_primary"], sub["pred_bucket"])
    conf_matrix = conf_matrix.reindex(index=fam_order_conf, columns=col_order, fill_value=0)

    cell_ids = {}
    for true_fam in fam_order_conf:
        for pred_fam in col_order:
            rows = sub[(sub["true_primary"] == true_fam) & (sub["pred_bucket"] == pred_fam)]
            if len(rows) > 0:
                cell_ids[f"{true_fam}|||{pred_fam}"] = rows["query_id"].astype(str).tolist()[:50]

    z = conf_matrix.values.astype(float)
    z_log = np.log10(z + 1)

    # customdata mirrors the (true, pred) key used by the JS click handler below;
    # kept for parity with the original build even though the handler reads
    # pt.x / pt.y directly rather than pt.customdata.
    customdata = np.empty(z.shape, dtype=object)
    for i, true_fam in enumerate(fam_order_conf):
        for j, pred_fam in enumerate(col_order):
            customdata[i, j] = f"{true_fam}|||{pred_fam}"

    hovertext = np.empty(z.shape, dtype=object)
    for i, true_fam in enumerate(fam_order_conf):
        for j, pred_fam in enumerate(col_order):
            hovertext[i, j] = f"true: {true_fam}<br>predicted: {pred_fam}<br>n = {int(z[i, j])}"

    fig4 = go.Figure(data=go.Heatmap(
        z=z_log, x=col_order, y=fam_order_conf,
        customdata=customdata,
        text=hovertext, hoverinfo="text",
        colorscale="Blues", showscale=True,
        colorbar=dict(title="log10(n+1)"),
        xgap=0.3, ygap=0.3,
    ))
    fig4.update_layout(
        title=dict(text=(f"Classifier confusion matrix — families with \u2265{min_support} eval support "
                          f"(n={len(fam_order_conf)} families, {sub.shape[0]} proteins). "
                          f"Click a cell to list protein ids."), x=0.5, font=dict(size=13)),
        xaxis=dict(title="predicted family (clf_pred)", tickfont=dict(size=6), tickangle=90),
        yaxis=dict(title="true family (primary token)", tickfont=dict(size=6), autorange="reversed"),
        width=1400, height=1300, margin=dict(t=100, b=150, l=150),
    )

    fig5 = go.Figure()
    fig5.add_trace(go.Bar(x=rel["mean_conf"], y=rel["emp_acc"], width=0.08,
                           marker=dict(color="#4c72b0"), name="empirical accuracy",
                           text=rel["n"],
                           hovertemplate="mean conf=%{x:.3f}<br>emp acc=%{y:.3f}<br>n=%{text}<extra></extra>"))
    fig5.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(color="black", dash="dash"),
                               name="perfect calibration"))
    fig5.update_layout(
        title=dict(text=f"Reliability diagram: clf_conf vs empirical accuracy (overlap-match) — 10 bins, ECE={ece:.3f}",
                   x=0.5, font=dict(size=13)),
        xaxis=dict(title="mean predicted confidence (clf_conf) in bin", range=[0, 1]),
        yaxis=dict(title="empirical accuracy (overlap-match)", range=[0, 1]),
        width=750, height=650, plot_bgcolor="white", legend=dict(x=0.02, y=0.98),
    )
    fig5.update_xaxes(showgrid=True, gridcolor="#eee")
    fig5.update_yaxes(showgrid=True, gridcolor="#eee")

    rel_div = fig5.to_html(full_html=False, include_plotlyjs=False, div_id="rel-div")
    conf_div = fig4.to_html(full_html=False, include_plotlyjs=False, div_id="conf-div")

    id_to_true = dict(zip(pred["query_id"].astype(str), pred["true_families"]))
    id_to_novelty = dict(zip(pred["query_id"].astype(str), pred["novelty"]))
    cell_ids_json = json.dumps(cell_ids)
    id_to_true_json = json.dumps(id_to_true)
    id_to_novelty_json = json.dumps(id_to_novelty)
    plotly_js_script = pyo.get_plotlyjs()

    html_page = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8"/>
<title>dbCAN4-advanced — Calibration & Confusion Explorer</title>
<script type="text/javascript">{plotly_js_script}</script>
<style>
  body {{ font-family: -apple-system, Arial, sans-serif; margin: 20px; color: #222; }}
  h1 {{ font-size: 20px; }}
  .container {{ display: flex; flex-wrap: wrap; gap: 30px; align-items: flex-start; }}
  #detail-panel {{
    border: 1px solid #ccc; border-radius: 6px; padding: 12px; width: 420px;
    max-height: 700px; overflow-y: auto; font-size: 12px; background: #fafafa;
  }}
  #detail-panel h3 {{ margin-top: 0; font-size: 14px; }}
  .idlist {{ list-style: none; padding: 0; margin: 0; }}
  .idlist li {{ padding: 2px 0; border-bottom: 1px solid #eee; }}
  .novel_family {{ color: #d9534f; font-weight: bold; }}
  .novel_seq {{ color: #333; }}
  .caveat {{ background: #fff3cd; border: 1px solid #ffe08a; border-radius: 4px; padding: 8px; font-size: 12px; margin-bottom: 16px; }}
</style>
</head>
<body>
<h1>dbCAN4-advanced — Calibration &amp; Confusion Explorer</h1>
<div class="caveat">
Source: <code>head_eval_pred.tsv</code> (4,726 proteins). Confusion matrix restricted to families with
&ge;{min_support} eval support ({len(fam_order_conf)} families, {sub.shape[0]} proteins covered; predictions
to families below this threshold are bucketed as "other"). Reliability diagram uses <code>clf_conf</code>
against overlap-match correctness (<code>clf_pred</code> in the true family set) across all 4,726 proteins, 10 fixed-width bins.
Click any confusion-matrix cell to list the protein ids in the panel on the right.
</div>
<div class="container">
  <div>{conf_div}</div>
  <div id="detail-panel"><h3>Click a confusion-matrix cell</h3><div id="detail-content">No cell selected yet.</div></div>
</div>
<hr/>
<div>{rel_div}</div>

<script type="text/javascript">
const CELL_IDS = {cell_ids_json};
const ID_TO_TRUE = {id_to_true_json};
const ID_TO_NOVELTY = {id_to_novelty_json};

function renderCell(trueFam, predFam) {{
  const key = trueFam + "|||" + predFam;
  const ids = CELL_IDS[key] || [];
  const panel = document.getElementById('detail-content');
  if (ids.length === 0) {{
    panel.innerHTML = "<p>No proteins in cell <b>" + trueFam + " &rarr; " + predFam + "</b>.</p>";
    return;
  }}
  let html = "<p><b>true:</b> " + trueFam + " &rarr; <b>predicted:</b> " + predFam +
             " (" + ids.length + (ids.length === 50 ? "+ shown, capped at 50" : "") + ")</p>";
  html += "<ul class='idlist'>";
  for (const id of ids) {{
    const novelty = ID_TO_NOVELTY[id] || '?';
    const trueFull = ID_TO_TRUE[id] || '?';
    html += "<li class='" + novelty + "'>" + id + " &nbsp;<span style='color:#888'>(true: " + trueFull + ", " + novelty + ")</span></li>";
  }}
  html += "</ul>";
  panel.innerHTML = html;
}}

document.addEventListener('DOMContentLoaded', function() {{
  const gd = document.getElementById('conf-div');
  gd.on('plotly_click', function(data) {{
    const pt = data.points[0];
    const trueFam = pt.y;
    const predFam = pt.x;
    renderCell(trueFam, predFam);
  }});
}});
</script>
</body>
</html>
"""
    with open(outpath, "w") as f:
        f.write(html_page)

    return rel, ece


def build_reliability_diagram_png(rel, ece, n_total, outpath):
    """Static figure-style companion to the reliability diagram."""
    apply_figure_style()
    fig6, ax6 = plt.subplots(figsize=(5.2, 4.8))
    ax6.bar(rel["mean_conf"], rel["emp_acc"], width=0.07, color="#4c72b0", edgecolor="white",
            label="empirical accuracy (this bin)")
    ax6.plot([0, 1], [0, 1], linestyle="--", color="black", lw=1.2, label="perfect calibration")
    ax6.set_xlim(0, 1.02)
    ax6.set_ylim(0, 1.05)
    ax6.set_xlabel("mean predicted confidence (clf_conf)")
    ax6.set_ylabel("empirical accuracy (overlap-match)")
    ax6.set_title(f"Overconfident in the dominant >0.9-confidence bin (90% of n={n_total}, ECE={ece:.2f})")
    ax6.legend(frameon=False, loc="upper left", fontsize=7)
    ax6.margins(0.03)
    fig6.tight_layout()
    fig6.savefig(outpath, dpi=300)
    plt.close(fig6)


# ----------------------------------------------------------------------------
# Runner
# ----------------------------------------------------------------------------

def main(assets_dir, outdir):
    os.makedirs(outdir, exist_ok=True)
    emb, hm, rs, log_lines, pred_raw = load_bundle(assets_dir)
    epochs, losses, accs = parse_train_log(log_lines)
    pred = add_correctness_columns(pred_raw)

    build_embedding_explorer(emb, os.path.join(outdir, "embedding_explorer.html"))

    build_training_dashboard(epochs, losses, accs, pred, os.path.join(outdir, "training_dashboard.html"))
    build_training_summary_png(epochs, losses, accs, hm, rs, os.path.join(outdir, "training_summary.png"))

    rel, ece = build_calibration_confusion_explorer(
        pred, os.path.join(outdir, "calibration_confusion_explorer.html"))
    build_reliability_diagram_png(rel, ece, len(pred), os.path.join(outdir, "reliability_diagram.png"))

    print(f"Wrote 3 HTML + 2 PNG deliverables to {outdir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assets", default=".", help="Directory containing the staged input files")
    parser.add_argument("--outdir", default=".", help="Directory to write the output deliverables")
    args = parser.parse_args()
    main(args.assets, args.outdir)
