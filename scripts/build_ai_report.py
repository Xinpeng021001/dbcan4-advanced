#!/usr/bin/env python3
"""
build_ai_report.py  --  dbCAN4-advanced per-protein "AI report" generator
=========================================================================

Purpose
-------
For each protein processed by the dbCAN4-advanced fungal CAZyme annotation
pipeline, assemble a single, self-contained, downloadable JSON "AI report"
in the **prompt-pack** format chosen by the project owner. A user pastes /
uploads the report into any LLM (Claude, ChatGPT, ...) and the LLM can then
describe the protein and answer questions about it -- grounded STRICTLY in
the evidence carried inside the report.

A report contains three parts and nothing else:

  1. system_prompt      -- a grounded instruction block telling the receiving
                           LLM to act as a CAZyme-annotation expert who uses
                           ONLY the evidence in the report, never invents a
                           family / EC / number / citation, explains tool
                           disagreements instead of silently picking, and says
                           "not determinable from the provided evidence" when
                           asked something the evidence does not cover.
  2. evidence           -- the FULL structured multi-tool evidence: per-head
                           and fusion CAZyme calls (with confidence /
                           agreement), sequence baselines, Pfam domains, CLEAN
                           EC, ESMFold structure, DeepTMHMM topology/secretion,
                           localization, physicochem, and tool provenance.
  3. suggested_questions-- starter questions the user can ask.

The report deliberately carries NO pre-written prose description of the
protein and NO baked-in answer -- the owner chose the pure prompt-pack, so
the receiving LLM generates the description at use time from the evidence.
(`annotation_summary` is the pipeline's own STRUCTURED verdict -- family,
confidence, agreement, review flag -- which is itself evidence, not prose.)

Grounding contract
------------------
Every value in the evidence block is read from the staged evidence files in
the --assets directory; nothing is hardcoded from prior knowledge. Where a
file is known to be stale (the 267317 ESMFold record in structures.tsv), the
value is cross-checked against the authoritative comprehensive JSON and the
correction is recorded transparently in a provenance note.

CLI
---
  # one protein
  python build_ai_report.py --assets DIR --protein 267317 --out 267317_ai_report.json

  # every protein found in the bundle
  python build_ai_report.py --assets DIR --all --outdir reports/

Evidence files consumed (all under --assets)
-------------------------------------------
  raw_knn.tsv            ESM-C kNN head          -> evidence.cazyme_calls.esm_c_heads.knn
  raw_centroid.tsv       ESM-C centroid head     -> evidence.cazyme_calls.esm_c_heads.centroid
  raw_contrastive.tsv    ESM-C contrastive head  -> evidence.cazyme_calls.esm_c_heads.contrastive
  head_eval_pred.tsv     contrastive sub-signals + curated true_families (enrichment)
  fusion_raw.tsv         fusion layer final call -> evidence.cazyme_calls.fusion
  diamond_eval2025_pred.tsv  DIAMOND vs 2025 ref -> evidence.cazyme_calls.sequence_baselines.diamond_2025ref
  real3_baseline_overview.tsv dbCAN3 standalone  -> evidence.cazyme_calls.sequence_baselines.dbcan3_standalone
  domains.tsv            Pfam hmmscan domains    -> evidence.pfam_domains
  ec_prediction.tsv      CLEAN seq->EC           -> evidence.ec_prediction
  structures.tsv         ESMFold pLDDT/length    -> evidence.structure
  deeptmhmm.tsv          topology / signal pep.  -> evidence.topology_secretion
  localization.tsv       secreted/extracellular  -> evidence.localization
  physicochem.tsv        MW/pI/instability/...   -> evidence.physicochem
  manifest.json          pipeline + tool versions-> evidence.tool_provenance
  <pid>_comprehensive.json  authoritative cross-check / foldseek enrichment (267317)
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import OrderedDict

SCHEMA_VERSION = "1.0.0"
GENERATED_BY = "dbcan4-advanced build_ai_report.py"

# Fusion abstain threshold used by the project fusion/consensus layer.
TAU_ABSTAIN = 0.35
# Confidence below which an otherwise-consistent call still warrants a second look.
HIGH_CONF_BAND = 0.90

# --------------------------------------------------------------------------
# Low-level file readers
# --------------------------------------------------------------------------

def read_tsv(path):
    """Read a tab-separated file into a list of OrderedDict rows.

    Handles the bundle's quoting convention (JSON blobs live in an `extra`
    column, quoted with doubled inner quotes) via the stdlib csv module.
    """
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [OrderedDict(row) for row in reader]


def index_by(rows, key):
    """Index rows (or lists of rows) by a key column."""
    out = {}
    for r in rows:
        out.setdefault(r[key], []).append(r)
    return out


def parse_extra(row):
    """Parse the JSON `extra` column if present; return {} otherwise."""
    raw = (row or {}).get("extra")
    if not raw:
        return {}
    try:
        return json.loads(raw)
    except (ValueError, TypeError):
        return {}


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def _num(x):
    """Best-effort numeric coercion; leaves non-numeric strings untouched."""
    if x is None or x == "" or x == "-":
        return None
    try:
        f = float(x)
        return int(f) if f.is_integer() else f
    except (ValueError, TypeError):
        return x


# --------------------------------------------------------------------------
# Assets container
# --------------------------------------------------------------------------

class Assets:
    """Lazy-ish loader that pulls every evidence file from an assets dir."""

    CORE_ID_FILES = ["raw_knn.tsv", "raw_centroid.tsv", "raw_contrastive.tsv",
                     "fusion_raw.tsv"]

    def __init__(self, assets_dir):
        self.dir = assets_dir
        self._cache = {}

    def path(self, name):
        return os.path.join(self.dir, name)

    def has(self, name):
        return os.path.exists(self.path(name))

    def tsv(self, name, key):
        """Return {id: [rows]} for a TSV, or {} if the file is absent."""
        ck = ("tsv", name, key)
        if ck not in self._cache:
            self._cache[ck] = index_by(read_tsv(self.path(name)), key) if self.has(name) else {}
        return self._cache[ck]

    def one(self, name, key, pid):
        """First row for a protein id in a TSV, or None."""
        rows = self.tsv(name, key).get(pid)
        return rows[0] if rows else None

    def comprehensive(self, pid):
        name = f"{pid}_comprehensive.json"
        return load_json(self.path(name)) if self.has(name) else None

    def manifest(self):
        return load_json(self.path("manifest.json")) if self.has("manifest.json") else {}

    def protein_ids(self):
        """Union of protein ids across the core prediction files, sorted."""
        ids = set()
        specs = [("raw_knn.tsv", "query_id"), ("raw_centroid.tsv", "query_id"),
                 ("raw_contrastive.tsv", "query_id"), ("fusion_raw.tsv", "protein_id")]
        for name, key in specs:
            ids.update(self.tsv(name, key).keys())
        return sorted(ids)


# --------------------------------------------------------------------------
# Evidence assembly (one function per evidence block)
# --------------------------------------------------------------------------

def build_curated_reference(A, pid):
    """Benchmark ground-truth label for this example protein (eval slice)."""
    row = (A.one("head_eval_pred.tsv", "query_id", pid)
           or A.one("esmc_retrieval_pred.tsv", "query_id", pid)
           or A.one("diamond_eval2025_pred.tsv", "query_id", pid))
    if not row:
        return None
    fams = row.get("true_families")
    return OrderedDict([
        ("curated_families", fams.split(",") if fams else []),
        ("novelty_bucket", row.get("novelty")),
        ("source_file", "head_eval_pred.tsv (true_families column)"),
        ("note", "Benchmark ground-truth label for this example/eval protein. "
                 "This is NOT available for a novel query protein in deployment; "
                 "it is included here so tool agreement/disagreement can be judged."),
    ])


def build_sequence_baselines(A, pid):
    out = OrderedDict()

    # dbCAN3 standalone (shipped baseline: dbCAN_hmm + dbCAN_sub + DIAMOND consensus)
    b = A.one("real3_baseline_overview.tsv", "Gene ID", pid)
    if b:
        ntools = _num(b.get("#ofTools"))
        recommend = b.get("Recommend Results")
        called = recommend not in (None, "", "-")
        out["dbcan3_standalone"] = OrderedDict([
            ("call", recommend if called else "no CAZyme call"),
            ("n_tools_supporting", ntools),
            ("dbCAN_hmm", b.get("dbCAN_hmm")),
            ("dbCAN_sub", b.get("dbCAN_sub")),
            ("DIAMOND", b.get("DIAMOND")),
            ("baseline_missed", (not called)),
            ("source_file", "real3_baseline_overview.tsv"),
            ("note", "dbCAN3 standalone (shipped baseline) returned no CAZyme call "
                     "for this protein -- a baseline-missed case that the advanced "
                     "ESM-C / structure tier is designed to recover."
                     if not called else
                     "dbCAN3 standalone consensus call."),
        ])

    # DIAMOND vs 2025 reference (research comparator; carries curated match flags)
    d = A.one("diamond_eval2025_pred.tsv", "query_id", pid)
    if d:
        out["diamond_2025ref"] = OrderedDict([
            ("call", d.get("pred_families")),
            ("top_percent_identity", _num(d.get("top_pident"))),
            ("top_evalue", d.get("top_evalue")),
            ("exact_match_vs_curated", _num(d.get("exact")) == 1),
            ("overlap_match_vs_curated", _num(d.get("overlap")) == 1),
            ("source_file", "diamond_eval2025_pred.tsv"),
            ("note", "DIAMOND against the 2025 CAZy reference (research comparator, "
                     "distinct from the dbCAN3 standalone run above)."),
        ])
    return out


def build_esm_c_heads(A, pid):
    out = OrderedDict()

    knn = A.one("raw_knn.tsv", "query_id", pid)
    if knn:
        out["knn"] = OrderedDict([
            ("predicted_family", knn.get("knn_pred")),
            ("confidence", _num(knn.get("knn_conf"))),
            ("neighborhood_purity", _num(knn.get("knn_purity"))),
            ("margin", _num(knn.get("knn_margin"))),
            ("source_file", "raw_knn.tsv"),
            ("description", "k-nearest-neighbour retrieval over the ESM-C embedding "
                            "reference; purity = fraction of the k neighbours sharing "
                            "the winning family."),
        ])

    cent = A.one("raw_centroid.tsv", "query_id", pid)
    if cent:
        out["centroid"] = OrderedDict([
            ("predicted_family", cent.get("cent_pred")),
            ("confidence", _num(cent.get("cent_conf"))),
            ("margin", _num(cent.get("cent_margin"))),
            ("source_file", "raw_centroid.tsv"),
            ("description", "nearest per-family centroid (CLEAN-style prototype) in "
                            "ESM-C embedding space; small margin = the top two family "
                            "prototypes are nearly equidistant."),
        ])

    contr = A.one("raw_contrastive.tsv", "query_id", pid)
    heval = A.one("head_eval_pred.tsv", "query_id", pid)
    if contr:
        block = OrderedDict([
            ("predicted_family", contr.get("clf_pred")),
            ("confidence", _num(contr.get("clf_conf"))),
            ("primary_signal", "classifier"),
            ("source_file", "raw_contrastive.tsv"),
        ])
        # Sub-signals inside the trained contrastive module.
        sub = OrderedDict()
        if contr.get("contr_knn_pred"):
            sub["contrastive_knn_pred"] = contr.get("contr_knn_pred")
            sub["contrastive_knn_purity"] = _num(contr.get("contr_knn_purity"))
        if heval and heval.get("contr_cent_pred"):
            sub["contrastive_centroid_pred"] = heval.get("contr_cent_pred")
            sub["contrastive_centroid_margin"] = _num(heval.get("contr_cent_margin"))
        if sub:
            sub["source_file"] = "raw_contrastive.tsv + head_eval_pred.tsv"
            block["retrieval_sub_signals"] = sub
        block["description"] = (
            "trained contrastive head. Its reported prediction is the softmax "
            "classifier; the retrieval sub-signals (contrastive kNN / centroid over "
            "the trained embedding) are shown separately and may disagree with the "
            "classifier.")
        out["contrastive"] = block
    return out


def _head_preds(esm_heads):
    """Top predicted family per primary head (for agreement analysis)."""
    return OrderedDict([
        ("knn", esm_heads.get("knn", {}).get("predicted_family")),
        ("centroid", esm_heads.get("centroid", {}).get("predicted_family")),
        ("contrastive", esm_heads.get("contrastive", {}).get("predicted_family")),
    ])


def build_fusion(A, pid):
    f = A.one("fusion_raw.tsv", "protein_id", pid)
    if not f:
        return None
    votes = {}
    try:
        votes = json.loads(f.get("votes") or "{}")
    except (ValueError, TypeError):
        votes = {}
    signals = []
    try:
        signals = json.loads(f.get("signals") or "[]")
    except (ValueError, TypeError):
        signals = []
    conf = _num(f.get("confidence"))
    agreement = _num(f.get("agreement"))
    fams = f.get("all_families")
    return OrderedDict([
        ("final_family", f.get("family")),
        ("confidence", conf),
        ("abstain", (conf is not None and conf < TAU_ABSTAIN)),
        ("abstain_threshold", TAU_ABSTAIN),
        ("agreement", f"{agreement}/4" if agreement is not None else None),
        ("agreement_count", agreement),
        ("candidate_families", fams.split(",") if fams else []),
        ("per_method_votes", votes),
        ("signals", signals),
        ("source_file", "fusion_raw.tsv"),
        ("description", "confidence-weighted vote across DIAMOND / HMMER / the three "
                        "ESM-C heads; calls below the abstain threshold are flagged as "
                        "putative-novel / uncertain."),
    ])


def build_pfam_domains(A, pid):
    rows = A.tsv("domains.tsv", "protein_id").get(pid, [])
    if not rows:
        return None
    domains = []
    for r in sorted(rows, key=lambda x: _num(x.get("start")) or 0):
        ex = parse_extra(r)
        domains.append(OrderedDict([
            ("pfam_accession", r.get("acc")),
            ("name", r.get("name")),
            ("seq_start", _num(r.get("start"))),
            ("seq_end", _num(r.get("end"))),
            ("evalue", r.get("evalue")),
            ("bitscore", _num(r.get("score"))),
            ("hmm_coverage", ex.get("hmm_coverage")),
        ]))
    architecture = " - ".join(d["name"] for d in domains)
    return OrderedDict([
        ("architecture", architecture),
        ("n_domains", len(domains)),
        ("domains", domains),
        ("source_file", "domains.tsv"),
        ("tool", "hmmscan (HMMER3) vs Pfam-A"),
    ])


def _ec_band(conf):
    if conf is None:
        return None
    if conf >= 0.50:
        return "HIGH"
    if conf >= 0.05:
        return "LOW"
    return "VERY_LOW"


def build_ec(A, pid):
    e = A.one("ec_prediction.tsv", "protein_id", pid)
    if not e:
        return None
    ex = parse_extra(e)
    conf = _num(e.get("confidence"))
    return OrderedDict([
        ("ec_number", e.get("ec_number")),
        ("confidence", conf),
        ("confidence_band", _ec_band(conf)),
        ("confidence_type", ex.get("confidence_type")),
        ("rank", _num(e.get("rank"))),
        ("tool", e.get("tool")),
        ("source_file", "ec_prediction.tsv"),
        ("description", "CLEAN sequence->EC prediction (contrastive-learning EC "
                        "assignment), independent of the family-inherited EC. Low "
                        "confidence means the sequence is far from CLEAN's EC "
                        "training space and the EC should be treated as weak."),
    ])


def build_structure(A, pid):
    s = A.one("structures.tsv", "protein_id", pid)
    if not s:
        return None
    ex = parse_extra(s)
    plddt = _num(s.get("plddt"))
    length = _num(s.get("length"))
    source = s.get("source")
    model = ex.get("model")
    note = None

    # Authoritative cross-check against the comprehensive JSON, if present.
    comp = A.comprehensive(pid)
    if comp and isinstance(comp.get("structure"), dict):
        c_plddt = comp["structure"].get("plddt_mean")
        c_len = comp["structure"].get("n_residues")
        if c_plddt is not None and plddt is not None and abs(c_plddt - plddt) > 0.05:
            note = (f"structures.tsv listed mean pLDDT {plddt} over {length} residues; "
                    f"corrected to {c_plddt} over {c_len} residues from "
                    f"{pid}_comprehensive.json (the structures.tsv record is stale).")
            plddt, length = c_plddt, (c_len if c_len is not None else length)

    block = OrderedDict([
        ("tool", source or "ESMFold"),
        ("model", model),
        ("mean_plddt", plddt),
        ("plddt_scale", "0-100"),
        ("length_aa", length),
        ("source_file", "structures.tsv"),
    ])
    if note:
        block["provenance_note"] = note

    # Foldseek structural-homology enrichment (only where the comprehensive JSON carries it).
    if comp and isinstance(comp.get("structure"), dict):
        cs = comp["structure"]
        if cs.get("fold"):
            block["fold"] = cs.get("fold")
        if cs.get("foldseek_top_hits"):
            fs = OrderedDict([
                ("reference", (cs.get("foldseek_summary") or {}).get("reference", "CAZyme3D")),
                ("tool", (cs.get("foldseek_summary") or {}).get("tool", "Foldseek")),
                ("n_hits", (cs.get("foldseek_summary") or {}).get("n_hits")),
                ("top_hits", [
                    OrderedDict([
                        ("target", h.get("target")),
                        ("family", h.get("family")),
                        ("tmscore", h.get("tmscore")),
                        ("prob", h.get("prob")),
                        ("evalue", h.get("evalue")),
                    ]) for h in cs["foldseek_top_hits"][:5]
                ]),
                ("source_file", f"{pid}_comprehensive.json"),
            ])
            block["foldseek"] = fs
    else:
        block["foldseek"] = "not available in this bundle for this protein"
    return block


def build_topology_secretion(A, pid):
    t = A.one("deeptmhmm.tsv", "protein_id", pid)
    if not t:
        return None
    ex = parse_extra(t)
    topo = t.get("topology") or ""
    sp_span = None
    for seg in topo.split("|"):
        if seg.startswith("signal:"):
            rng = seg.split(":", 1)[1]
            try:
                a, b = rng.split("-")
                sp_span = [int(a), int(b)]
            except ValueError:
                sp_span = None
    has_sp = ex.get("has_signal_peptide")
    return OrderedDict([
        ("tool", ex.get("tool", "DeepTMHMM")),
        ("prediction", t.get("prediction")),
        ("n_tm_helices", _num(t.get("n_tm"))),
        ("has_signal_peptide", has_sp),
        ("signal_peptide_span", sp_span),
        ("topology", topo),
        ("secretion_call", "classically secreted (signal peptide, 0 TM helices)"
            if has_sp else "no signal peptide (globular / non-secretory-pathway)"),
        ("source_file", "deeptmhmm.tsv"),
    ])


def build_localization(A, pid):
    l = A.one("localization.tsv", "protein_id", pid)
    if not l:
        return None
    ex = parse_extra(l)
    return OrderedDict([
        ("localization", l.get("localization")),
        ("confidence", l.get("confidence") or None),
        ("method", l.get("method")),
        ("basis", ex.get("basis")),
        ("sp_prediction", ex.get("sp_prediction")),
        ("source_file", "localization.tsv"),
        ("note", "Derived from the signal-peptide call (labelled derived-from-SP), "
                 "not a DeepLoc run."),
    ])


def _instability_class(v):
    if v is None:
        return None
    return "stable" if v < 40 else "unstable"


def build_physicochem(A, pid):
    p = A.one("physicochem.tsv", "protein_id", pid)
    if not p:
        return None
    ex = parse_extra(p)
    inst = _num(p.get("instability"))
    return OrderedDict([
        ("molecular_weight_da", _num(p.get("mw_da"))),
        ("theoretical_pi", _num(p.get("pi"))),
        ("instability_index", inst),
        ("instability_class", _instability_class(inst)),
        ("gravy", _num(p.get("gravy"))),
        ("aromaticity", _num(p.get("aromaticity"))),
        ("n_glycosylation_sequons", ex.get("n_glyc_sequons")),
        ("source_file", "physicochem.tsv"),
        ("tool", "Biopython ProtParam"),
    ])


def build_tool_provenance(A):
    m = A.manifest()
    return OrderedDict([
        ("pipeline", m.get("pipeline", "dbcan4-advanced")),
        ("pipeline_version", m.get("pipeline_version")),
        ("contract_version", m.get("contract_version")),
        ("release_label", m.get("release_label")),
        ("tool_versions", m.get("tool_versions", {})),
        ("assets_manifest", "manifest.json"),
    ])


# --------------------------------------------------------------------------
# Review-flag logic (head agreement + fusion-vs-best-head mismatch)
# --------------------------------------------------------------------------

def compute_review_flag(esm_heads, fusion, pfam, baselines):
    """Return a 3-level review verdict computed from tool agreement.

    Levels: 'clean' < 'attention' < 'review'. `flagged` is True only at
    'review'. Mirrors the project triage rule: low head agreement, or fusion
    following a wrong-but-confident head, or a sub-threshold fusion confidence.
    """
    heads = _head_preds(esm_heads)
    head_vals = [v for v in heads.values() if v]
    distinct = sorted(set(head_vals))
    n_heads = len(head_vals)

    # Majority family among the primary heads.
    majority_family, majority_count = None, 0
    for fam in distinct:
        c = head_vals.count(fam)
        if c > majority_count:
            majority_family, majority_count = fam, c
    has_majority = majority_count >= 2

    fusion_family = (fusion or {}).get("final_family")
    fusion_conf = (fusion or {}).get("confidence")
    fusion_agree = (fusion or {}).get("agreement_count")

    review_reasons, attention_reasons = [], []

    # ---- REVIEW-level triggers ----
    if fusion_conf is not None and fusion_conf < TAU_ABSTAIN:
        review_reasons.append(
            f"fusion confidence {fusion_conf:.4f} is below the abstain threshold "
            f"{TAU_ABSTAIN:.2f} (fusion is effectively abstaining)")
    if n_heads >= 3 and len(distinct) == n_heads:
        review_reasons.append(
            f"all {n_heads} ESM-C heads disagree ({', '.join(heads[k] for k in heads if heads[k])}) "
            f"-- no majority family")
    if fusion_agree is not None and fusion_agree <= 2:
        review_reasons.append(f"fusion agreement is only {fusion_agree}/4")
    if has_majority and fusion_family and fusion_family != majority_family:
        review_reasons.append(
            f"fusion picked {fusion_family} but the head majority is {majority_family} "
            f"(fusion may be following a confident-but-outvoted head)")

    # ---- ATTENTION-level triggers (only matter if not already REVIEW) ----
    if len(distinct) > 1:
        attention_reasons.append(
            f"ESM-C heads are not unanimous ({', '.join(f'{k}={heads[k]}' for k in heads if heads[k])})")
    if fusion_conf is not None and fusion_conf < HIGH_CONF_BAND:
        attention_reasons.append(
            f"fusion confidence {fusion_conf:.4f} is below the high-confidence band "
            f"({HIGH_CONF_BAND:.2f})")
    if pfam and pfam.get("n_domains", 0) > 1:
        attention_reasons.append(
            f"multidomain architecture ({pfam.get('architecture')}) -- a single family "
            f"label may not capture the protein")
    diamond = (baselines or {}).get("diamond_2025ref")
    if diamond and fusion_family and diamond.get("call") and diamond.get("call") != fusion_family:
        attention_reasons.append(
            f"sequence baseline (DIAMOND {diamond.get('call')}) disagrees with the "
            f"advanced fusion call ({fusion_family})")

    if review_reasons:
        level = "review"
        reasons = review_reasons + attention_reasons
    elif attention_reasons:
        level = "attention"
        reasons = attention_reasons
    else:
        level = "clean"
        reasons = ["all heads agree, fusion is high-confidence and unanimous"]

    return OrderedDict([
        ("review_level", level),
        ("flagged", level == "review"),
        ("head_predictions", heads),
        ("head_majority_family", majority_family),
        ("head_majority_count", f"{majority_count}/{n_heads}" if n_heads else None),
        ("reasons", reasons),
    ])


def build_annotation_summary(esm_heads, fusion, review):
    """Structured pipeline verdict (NOT a prose description)."""
    level = review["review_level"]
    if level == "review":
        interp = ("flagged for expert review: the tools do not converge and the "
                  "fusion call should not be trusted without checking the individual "
                  "heads and baselines.")
    elif level == "attention":
        interp = ("correct/consistent top call but with genuine complexity (e.g. "
                  "multidomain architecture or a diverging head) -- worth a look, "
                  "not a hard review block.")
    else:
        interp = "clean unanimous call; no disagreement signal."
    return OrderedDict([
        ("final_family", (fusion or {}).get("final_family")),
        ("final_confidence", (fusion or {}).get("confidence")),
        ("fusion_agreement", (fusion or {}).get("agreement")),
        ("tools_in_agreement", review["head_majority_count"]),
        ("review_flag", OrderedDict([
            ("review_level", review["review_level"]),
            ("flagged", review["flagged"]),
            ("reasons", review["reasons"]),
        ])),
        ("verdict", interp),
    ])


# --------------------------------------------------------------------------
# System prompt + suggested questions
# --------------------------------------------------------------------------

def build_system_prompt(pid):
    return (
        "You are a carbohydrate-active-enzyme (CAZyme) annotation expert helping a "
        "researcher interpret an automated annotation for a single fungal protein "
        f"(protein_id {pid}) produced by the dbCAN4-advanced pipeline.\n\n"
        "This message contains a self-contained evidence report. You must follow "
        "these rules WITHOUT EXCEPTION:\n\n"
        "1. GROUND EVERYTHING IN THE REPORT. Use ONLY the evidence provided in the "
        "`evidence` and `annotation_summary` blocks of this report. Do not use outside "
        "knowledge to assert facts about THIS protein.\n"
        "2. NEVER INVENT. Do not invent or guess a CAZy family, EC number, domain, "
        "numeric value (confidence, pLDDT, coverage, MW, etc.), organism, or literature "
        "citation. Every family/EC/number you state must appear verbatim in the report.\n"
        "3. EXPLAIN DISAGREEMENT, DO NOT SILENTLY PICK. When the tools disagree (the "
        "kNN, centroid, and contrastive ESM-C heads; the fusion call; the sequence "
        "baselines), lay out what each method said and why they differ. If the report "
        "carries a review flag, foreground it and explain the reason rather than "
        "presenting a single confident answer.\n"
        "4. RESPECT CONFIDENCE. Report confidence values and their bands as given. Treat "
        "low-confidence or abstaining calls (fusion below its abstain threshold, a "
        "VERY_LOW/LOW EC) as weak evidence and say so.\n"
        "5. SAY WHEN YOU CANNOT ANSWER. If the user asks something the report does not "
        "cover, respond exactly: \"not determinable from the provided evidence\" (and "
        "state what additional evidence would be needed). Do not fill the gap with "
        "general knowledge presented as fact about this protein.\n"
        "6. GENERAL CAZyme BACKGROUND is allowed ONLY when clearly framed as general "
        "context (e.g. what a GH78 family does in general), never as a measured property "
        "of this specific protein, and never to override or add to the report's evidence.\n\n"
        "When you answer, prefer to: (a) give the pipeline's final call and its "
        "confidence, (b) summarize the supporting and dissenting evidence, (c) surface "
        "any review flag and its cause, and (d) state the limitations. You are an "
        "interpreter of the evidence in this report, not an independent predictor."
    )


def build_about():
    """Compact, grounded method glossary so the report is self-contained."""
    return OrderedDict([
        ("pipeline", "dbCAN4-advanced fungal CAZyme annotation (standalone prototype)"),
        ("family_call_methods", OrderedDict([
            ("sequence_baselines", "dbCAN3 standalone (HMMER + dbCAN_sub + DIAMOND "
                                   "consensus) and DIAMOND vs a 2025 CAZy reference."),
            ("esm_c_heads", "Three retrieval heads over ESM-C protein-language-model "
                            "embeddings: kNN (majority vote of nearest neighbours), "
                            "centroid (nearest per-family prototype), and a trained "
                            "contrastive head (softmax classifier + contrastive "
                            "kNN/centroid sub-signals)."),
            ("fusion", "A confidence-weighted vote across the baselines and the ESM-C "
                       "heads; produces the final family + confidence + agreement, and "
                       "abstains below its confidence threshold."),
        ])),
        ("supporting_evidence", "Pfam domain architecture (hmmscan), CLEAN sequence->EC, "
                                "ESMFold structure (mean pLDDT) with optional Foldseek "
                                "structural homology, DeepTMHMM topology/signal-peptide, "
                                "localization, and physicochemistry."),
        ("confidence_semantics", OrderedDict([
            ("fusion_abstain_threshold", TAU_ABSTAIN),
            ("ec_bands", "HIGH >= 0.50, LOW 0.05-0.50, VERY_LOW < 0.05"),
            ("plddt_scale", "0-100 (per-residue confidence, mean reported)"),
        ])),
    ])


def build_suggested_questions(pid, review, fusion, pfam, curated):
    fam = (fusion or {}).get("final_family")
    qs = [
        "Describe this protein: what is the pipeline's best call, how confident is it, "
        "and what evidence supports it?",
        "Do the different tools agree on the CAZy family? Walk me through what each "
        "method (baselines, kNN, centroid, contrastive, fusion) predicted.",
        "Is this protein secreted, and what does the structure/topology evidence say?",
        "How reliable is the predicted EC number, and does it fit the family call?",
        "What are the main limitations or caveats of this annotation?",
    ]
    # Protein-specific starters driven by the data.
    if review["review_level"] == "review":
        qs.insert(1, "This annotation is flagged for review -- why, and which "
                     "prediction (if any) should I trust?")
        if curated and curated.get("curated_families"):
            qs.insert(2, "Which tool, if any, recovered the curated family "
                         f"({','.join(curated['curated_families'])}), and which tools missed it?")
    elif review["review_level"] == "attention":
        if pfam and pfam.get("n_domains", 0) > 1:
            qs.insert(1, f"This protein has multiple Pfam domains ({pfam.get('architecture')}) "
                         "-- how should I interpret a single family label, and why do the "
                         "sequence baseline and the ESM-C call differ?")
    if fam:
        qs.append(f"What substrate class and reaction would a {fam} enzyme typically act "
                  "on (general context only)?")
    return qs


# --------------------------------------------------------------------------
# Top-level report assembly
# --------------------------------------------------------------------------

def build_report(A, pid):
    esm_heads = build_esm_c_heads(A, pid)
    fusion = build_fusion(A, pid)
    pfam = build_pfam_domains(A, pid)
    baselines = build_sequence_baselines(A, pid)
    curated = build_curated_reference(A, pid)

    review = compute_review_flag(esm_heads, fusion, pfam, baselines)

    evidence = OrderedDict()
    if curated:
        evidence["curated_reference"] = curated
    evidence["cazyme_calls"] = OrderedDict([
        ("sequence_baselines", baselines),
        ("esm_c_heads", esm_heads),
        ("fusion", fusion),
    ])
    evidence["pfam_domains"] = pfam
    evidence["ec_prediction"] = build_ec(A, pid)
    evidence["structure"] = build_structure(A, pid)
    evidence["topology_secretion"] = build_topology_secretion(A, pid)
    evidence["localization"] = build_localization(A, pid)
    evidence["physicochem"] = build_physicochem(A, pid)
    evidence["tool_provenance"] = build_tool_provenance(A)

    report = OrderedDict([
        ("schema_version", SCHEMA_VERSION),
        ("generated_by", GENERATED_BY),
        ("report_type", "prompt_pack"),
        ("protein_id", pid),
        ("about", build_about()),
        ("system_prompt", build_system_prompt(pid)),
        ("evidence", evidence),
        ("annotation_summary", build_annotation_summary(esm_heads, fusion, review)),
        ("suggested_questions", build_suggested_questions(pid, review, fusion, pfam, curated)),
        ("usage", "Paste or upload this entire JSON into an LLM. The `system_prompt` "
                  "field is the instruction block; the rest is the grounded evidence. "
                  "Then ask any of the suggested_questions or your own question."),
    ])
    return report


def write_report(report, out_path):
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")
    return out_path


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Generate dbCAN4-advanced per-protein AI reports (prompt-pack JSON).")
    ap.add_argument("--assets", required=True, help="directory with the staged evidence files")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--protein", help="single protein id to report on")
    g.add_argument("--all", action="store_true", help="one report per protein found in the bundle")
    ap.add_argument("--out", help="output path for --protein (default: <id>_ai_report.json)")
    ap.add_argument("--outdir", default="reports", help="output directory for --all (default: reports/)")
    args = ap.parse_args(argv)

    if not os.path.isdir(args.assets):
        ap.error(f"--assets directory not found: {args.assets}")

    A = Assets(args.assets)

    if args.protein:
        pid = args.protein
        if pid not in A.protein_ids():
            ap.error(f"protein {pid} not found in {args.assets} "
                     f"(available: {', '.join(A.protein_ids())})")
        out = args.out or f"{pid}_ai_report.json"
        report = build_report(A, pid)
        write_report(report, out)
        print(f"wrote {out}  (review_level="
              f"{report['annotation_summary']['review_flag']['review_level']})")
    else:
        ids = A.protein_ids()
        if not ids:
            ap.error(f"no proteins found in {args.assets}")
        for pid in ids:
            out = os.path.join(args.outdir, f"{pid}_ai_report.json")
            report = build_report(A, pid)
            write_report(report, out)
            print(f"wrote {out}  (review_level="
                  f"{report['annotation_summary']['review_flag']['review_level']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
