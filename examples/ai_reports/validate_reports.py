#!/usr/bin/env python3
"""
validate_reports.py -- ground-truth validator for dbCAN4-advanced AI reports.

For each report it checks three things:
  (A) the file is well-formed JSON with the required top-level structure;
  (B) every family / EC / key numeric value in the evidence block traces back
      verbatim to the source TSV/JSON in the assets bundle (anti-fabrication);
  (C) the review flag fires for 169208 and does NOT fire for 602276.

It re-reads the raw source files independently of build_ai_report.py, so it is
a genuine cross-check, not a restatement of the generator's own logic.
"""
import csv, json, os, sys
from collections import OrderedDict

ASSETS = sys.argv[1] if len(sys.argv) > 1 else "assets"
REPORT_DIR = sys.argv[2] if len(sys.argv) > 2 else "."

def rd(name):
    with open(os.path.join(ASSETS, name), newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def by(rows, key):
    d = {}
    for r in rows:
        d[r[key]] = r  # last wins; fine, ids are unique per file here
    return d

# ---- independently load the source truth ----
KNN   = by(rd("raw_knn.tsv"), "query_id")
CENT  = by(rd("raw_centroid.tsv"), "query_id")
CONTR = by(rd("raw_contrastive.tsv"), "query_id")
FUS   = by(rd("fusion_raw.tsv"), "protein_id")
EC    = by(rd("ec_prediction.tsv"), "protein_id")
PHYS  = by(rd("physicochem.tsv"), "protein_id")
STRUCT= by(rd("structures.tsv"), "protein_id")
DIA   = by(rd("diamond_eval2025_pred.tsv"), "query_id")
HEVAL = by(rd("head_eval_pred.tsv"), "query_id")
DOM = {}
for r in rd("domains.tsv"):
    DOM.setdefault(r["protein_id"], []).append(r)
COMP = {}
for pid in ("267317", "602276", "169208"):
    p = os.path.join(ASSETS, f"{pid}_comprehensive.json")
    if os.path.exists(p):
        COMP[pid] = json.load(open(p))

def approx(a, b, tol=1e-6):
    return a is not None and b is not None and abs(float(a) - float(b)) <= tol

results = []
def check(pid, name, ok, detail=""):
    results.append((pid, name, bool(ok), detail))

REQUIRED_TOP = ["schema_version", "generated_by", "protein_id", "system_prompt",
                "evidence", "annotation_summary", "suggested_questions"]
REQUIRED_EVID = ["cazyme_calls", "pfam_domains", "ec_prediction", "structure",
                 "topology_secretion", "localization", "physicochem", "tool_provenance"]

for pid in ("267317", "602276", "169208"):
    path = os.path.join(REPORT_DIR, f"{pid}_ai_report.json")

    # (A) well-formed JSON + structure
    try:
        R = json.load(open(path))
        check(pid, "well_formed_json", True)
    except Exception as e:
        check(pid, "well_formed_json", False, str(e))
        continue
    check(pid, "top_level_keys", all(k in R for k in REQUIRED_TOP),
          "missing: " + ",".join(k for k in REQUIRED_TOP if k not in R))
    ev = R.get("evidence", {})
    check(pid, "evidence_keys", all(k in ev for k in REQUIRED_EVID),
          "missing: " + ",".join(k for k in REQUIRED_EVID if k not in ev))
    # prompt-pack purity: no baked-in prose description/answer field
    check(pid, "no_baked_description",
          not any(k in R for k in ("description", "answer", "summary_text", "narrative")),
          "prompt pack must carry evidence+prompt+questions only")

    heads = ev["cazyme_calls"]["esm_c_heads"]
    fus = ev["cazyme_calls"]["fusion"]

    # (B) trace every head family + confidence to source
    check(pid, "knn_family", heads["knn"]["predicted_family"] == KNN[pid]["knn_pred"])
    check(pid, "knn_conf", approx(heads["knn"]["confidence"], KNN[pid]["knn_conf"]))
    check(pid, "knn_purity", approx(heads["knn"]["neighborhood_purity"], KNN[pid]["knn_purity"]))
    check(pid, "centroid_family", heads["centroid"]["predicted_family"] == CENT[pid]["cent_pred"])
    check(pid, "centroid_conf", approx(heads["centroid"]["confidence"], CENT[pid]["cent_conf"]))
    check(pid, "centroid_margin", approx(heads["centroid"]["margin"], CENT[pid]["cent_margin"]))
    check(pid, "contrastive_family", heads["contrastive"]["predicted_family"] == CONTR[pid]["clf_pred"])
    check(pid, "contrastive_conf", approx(heads["contrastive"]["confidence"], CONTR[pid]["clf_conf"]))

    # fusion final call + confidence + agreement
    check(pid, "fusion_family", fus["final_family"] == FUS[pid]["family"])
    check(pid, "fusion_conf", approx(fus["confidence"], FUS[pid]["confidence"]))
    check(pid, "fusion_agreement", str(fus["agreement_count"]) == FUS[pid]["agreement"])
    src_fams = FUS[pid]["all_families"].split(",")
    check(pid, "fusion_candidates", fus["candidate_families"] == src_fams)

    # EC
    eco = ev["ec_prediction"]
    check(pid, "ec_number", eco["ec_number"] == EC[pid]["ec_number"])
    check(pid, "ec_conf", approx(eco["confidence"], EC[pid]["confidence"]))

    # Pfam: every accession + coordinates traced
    dom_ok = True
    src_doms = {d["acc"]: d for d in DOM[pid]}
    for d in ev["pfam_domains"]["domains"]:
        s = src_doms.get(d["pfam_accession"])
        if not s or d["name"] != s["name"] or str(d["seq_start"]) != s["start"] or str(d["seq_end"]) != s["end"]:
            dom_ok = False
    check(pid, "pfam_domains_traced", dom_ok and len(ev["pfam_domains"]["domains"]) == len(DOM[pid]))

    # physicochem
    ph = ev["physicochem"]
    check(pid, "mw", approx(ph["molecular_weight_da"], PHYS[pid]["mw_da"]))
    check(pid, "pi", approx(ph["theoretical_pi"], PHYS[pid]["pi"]))
    check(pid, "gravy", approx(ph["gravy"], PHYS[pid]["gravy"]))

    # structure: 267317 must use the corrected comprehensive values, others match TSV
    st = ev["structure"]
    if pid == "267317":
        c = COMP[pid]["structure"]
        check(pid, "struct_plddt_corrected", approx(st["mean_plddt"], c["plddt_mean"]),
              f"expected {c['plddt_mean']} got {st['mean_plddt']}")
        check(pid, "struct_len_corrected", st["length_aa"] == c["n_residues"],
              f"expected {c['n_residues']} got {st['length_aa']}")
        check(pid, "struct_stale_noted", "provenance_note" in st)
    else:
        check(pid, "struct_plddt", approx(st["mean_plddt"], STRUCT[pid]["plddt"]))
        check(pid, "struct_len", st["length_aa"] == int(STRUCT[pid]["length"]))

    # DIAMOND baseline traced (incl. the 169208 GH183-recovery flag)
    db = ev["cazyme_calls"]["sequence_baselines"].get("diamond_2025ref")
    if db:
        check(pid, "diamond_call", db["call"] == DIA[pid]["pred_families"])
        check(pid, "diamond_exact", db["exact_match_vs_curated"] == (DIA[pid]["exact"] == "1"))

    # curated true family traced
    cr = ev.get("curated_reference")
    if cr:
        check(pid, "curated_family", cr["curated_families"] == HEVAL[pid]["true_families"].split(","))

    # (C) review-flag expectations
    flagged = R["annotation_summary"]["review_flag"]["flagged"]
    if pid == "169208":
        check(pid, "REVIEW_FLAG_fires", flagged is True, "169208 must be flagged")
    elif pid == "602276":
        check(pid, "REVIEW_FLAG_absent", flagged is False, "602276 must NOT be flagged")
    else:  # 267317 intermediate: not a hard flag, but level == attention
        lvl = R["annotation_summary"]["review_flag"]["review_level"]
        check(pid, "267317_intermediate", (flagged is False and lvl == "attention"),
              f"expected attention/not-flagged, got {lvl}/{flagged}")

# ---- report ----
n_pass = sum(1 for *_, ok, _ in ((r[0], r[1], r[2], r[3]) for r in results) if ok)
by_pid = {}
for pid, name, ok, detail in results:
    by_pid.setdefault(pid, []).append((name, ok, detail))

print("=" * 70)
print("dbCAN4-advanced AI-report validation")
print("=" * 70)
all_ok = True
for pid in ("267317", "602276", "169208"):
    rows = by_pid.get(pid, [])
    fails = [(n, d) for n, ok, d in rows if not ok]
    status = "PASS" if not fails else "FAIL"
    if fails:
        all_ok = False
    print(f"\n[{pid}]  {len([1 for _,ok,_ in rows if ok])}/{len(rows)} checks  -> {status}")
    for n, d in fails:
        print(f"    FAIL: {n}  {d}")

print("\n" + "-" * 70)
total = len(results)
passed = sum(1 for *_, ok, _ in [(r[0], r[1], r[2], r[3]) for r in results] if ok)
print(f"TOTAL: {passed}/{total} checks passed across 3 reports")
f169 = json.load(open(os.path.join(REPORT_DIR, "169208_ai_report.json")))["annotation_summary"]["review_flag"]
f602 = json.load(open(os.path.join(REPORT_DIR, "602276_ai_report.json")))["annotation_summary"]["review_flag"]
f267 = json.load(open(os.path.join(REPORT_DIR, "267317_ai_report.json")))["annotation_summary"]["review_flag"]
print()
print(f"  169208: flagged={f169['flagged']} (level={f169['review_level']})  [expect True/review]")
print(f"  602276: flagged={f602['flagged']} (level={f602['review_level']})  [expect False/clean]")
print(f"  267317: flagged={f267['flagged']} (level={f267['review_level']})  [expect False/attention]")
print("-" * 70)
print("OVERALL:", "ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
sys.exit(0 if all_ok else 1)
