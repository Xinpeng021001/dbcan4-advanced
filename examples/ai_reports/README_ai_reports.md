# dbCAN4-advanced — per-protein "AI report" generator

A prototype that turns the dbCAN4-advanced multi-tool evidence for one fungal
protein into a single, self-contained **prompt-pack JSON**. A user pastes or
uploads that JSON into any LLM (Claude, ChatGPT, …) and the LLM can then
describe the protein and answer questions about it — grounded **strictly** in
the evidence carried inside the report.

This is a standalone review prototype. Nothing here touches a git repo or gets
pushed anywhere.

---

## 1. What a report is (the "prompt pack")

Each report is a JSON object with three functional parts and nothing else:

1. **`system_prompt`** — a grounded instruction block for the receiving LLM. It
   tells the LLM to act as a CAZyme-annotation expert who uses **only** the
   evidence in the report, **never invents** a family / EC / number / citation,
   **explains disagreement** instead of silently picking a family, respects the
   reported confidences, and answers *"not determinable from the provided
   evidence"* when asked something the report does not cover.
2. **`evidence`** — the full structured multi-tool evidence (see the field map
   below). This is the ground truth the LLM must reason from.
3. **`suggested_questions`** — starter questions, partly tailored to the
   protein (e.g. a flagged protein gets a "why is this flagged?" starter).

Supporting fields: `about` (a compact, grounded glossary of the pipeline's
methods and confidence semantics so the pack is self-contained), an
`annotation_summary` (the pipeline's **structured** verdict — final family,
confidence, agreement, and the review flag), and a `usage` note.

**The report deliberately carries no pre-written prose description and no
baked-in answer.** The owner chose the pure prompt-pack: the receiving LLM
writes the description at use time from the evidence. `annotation_summary` is a
structured verdict (values, not prose), so it counts as evidence, not a
description. (The validator enforces this: it rejects any top-level
`description`/`answer`/`narrative` field.)

---

## 2. Report schema (top level)

```
schema_version        "1.0.0"
generated_by          "dbcan4-advanced build_ai_report.py"
report_type           "prompt_pack"
protein_id            e.g. "267317"
about                 { pipeline, family_call_methods, supporting_evidence,
                        confidence_semantics }        # grounded glossary
system_prompt         <the grounded instruction block, verbatim>
evidence              { … see field map below … }
annotation_summary    { final_family, final_confidence, fusion_agreement,
                        tools_in_agreement, review_flag{…}, verdict }
suggested_questions   [ "…", … ]
usage                 <how to paste it into an LLM>
```

### Evidence field → source-file map

Every evidence value is read from the staged bundle; nothing is hardcoded.
Each evidence block also carries its own `source_file` so provenance is visible
inside the report itself.

| `evidence.` field | Source file in `--assets` | Notes |
|---|---|---|
| `curated_reference` | `head_eval_pred.tsv` (`true_families`) | Benchmark ground-truth label for these example/eval proteins; **not** available for a novel deployment query. Included so agreement/disagreement can be judged. |
| `cazyme_calls.sequence_baselines.dbcan3_standalone` | `real3_baseline_overview.tsv` | dbCAN3 standalone (HMMER + dbCAN_sub + DIAMOND consensus). For all 3 examples this is a **baseline-miss** (0 tools). |
| `cazyme_calls.sequence_baselines.diamond_2025ref` | `diamond_eval2025_pred.tsv` | DIAMOND vs a 2025 CAZy reference (research comparator); carries curated exact/overlap flags. |
| `cazyme_calls.esm_c_heads.knn` | `raw_knn.tsv` | family, confidence, neighborhood purity, margin. |
| `cazyme_calls.esm_c_heads.centroid` | `raw_centroid.tsv` | family, confidence, margin. |
| `cazyme_calls.esm_c_heads.contrastive` | `raw_contrastive.tsv` (+ `head_eval_pred.tsv`) | classifier family + confidence, plus the contrastive-kNN / contrastive-centroid **sub-signals** (which can disagree with the classifier). |
| `cazyme_calls.fusion` | `fusion_raw.tsv` | final family, confidence, abstain flag, agreement (x/4), candidate families, per-method votes. |
| `pfam_domains` | `domains.tsv` | accession, name, coordinates, e-value, bitscore, HMM coverage; domains sorted by start. |
| `ec_prediction` | `ec_prediction.tsv` | CLEAN seq→EC, confidence + band (HIGH/LOW/VERY_LOW). |
| `structure` | `structures.tsv` (+ `<pid>_comprehensive.json`) | ESMFold mean pLDDT + length; Foldseek enrichment where the comprehensive JSON has it. **See the correction rule below.** |
| `topology_secretion` | `deeptmhmm.tsv` | signal peptide, TM helices, topology string, secretion call. |
| `localization` | `localization.tsv` | derived-from-SP call (labelled; not a DeepLoc run). |
| `physicochem` | `physicochem.tsv` | MW, pI, instability (+class), GRAVY, aromaticity, N-glyc sequons. |
| `tool_provenance` | `manifest.json` | pipeline + tool versions. |

### Structure-record correction (267317)

`structures.tsv` lists 267317 as mean pLDDT **69.2 over 1089 residues**, which
is stale. The generator cross-checks against `267317_comprehensive.json`
(authoritative: **76.6 over 1088 residues**) and, when they differ, uses the
comprehensive value and records the substitution in
`evidence.structure.provenance_note`. The check is generic (any protein with a
comprehensive JSON), not a 267317 special case.

---

## 3. Review-flag logic

Computed from tool agreement in `compute_review_flag()`. Three levels:

- **`clean`** — all three ESM-C heads agree, fusion is unanimous (4/4) and
  high-confidence. → *not flagged.*
- **`attention`** — a correct/consistent top call but with real complexity:
  heads not unanimous, fusion confidence below the 0.90 high-confidence band,
  a multidomain Pfam architecture, or a sequence-baseline vs fusion mismatch.
  → *not flagged*, but surfaced.
- **`review`** (**flagged**) — any of: fusion confidence below the abstain
  threshold **τ = 0.35**; all heads mutually disagree (no majority); fusion
  agreement ≤ 2/4; or fusion picking a family that contradicts the head
  majority (fusion following a confident-but-outvoted head).

This mirrors the project triage rule (low head agreement, or fusion following a
wrong-but-confident head, ⇒ flag). Outcomes for the three examples:

| Protein | True family | Heads (kNN / centroid / contrastive) | Fusion | Level | Flagged |
|---|---|---|---|---|---|
| **602276** | GH11 | GH11 / GH11 / GH11 | GH11 @ 0.9823, 4/4 | `clean` | **No** |
| **267317** | GH28,GH78 | GH78 / GH92 / GH78 | GH78 @ 0.6671, 3/4 | `attention` | **No** (multidomain complexity) |
| **169208** | GH183 | GH43_6 / GH183 / PL42 | GH43_6 @ 0.3082, 2/4 | `review` | **Yes** |

169208 is the honest hard case: **only the centroid head recovered the true
GH183**; fusion followed the confident-but-wrong kNN (GH43_6) and abstains
(0.3082 < τ). The report says this plainly and shows that the DIAMOND-2025
baseline *did* recover GH183 while dbCAN3 standalone missed the protein
entirely.

---

## 4. Validation

`validate_reports.py` re-reads the raw source files **independently** of the
generator and checks, per report:

- **(A)** well-formed JSON with all required top-level and evidence keys, and
  prompt-pack purity (no baked-in description field).
- **(B)** every head/fusion family, every confidence/purity/margin, EC, Pfam
  accession + coordinates, MW/pI/GRAVY, and DIAMOND/curated call traces back
  **verbatim** to the source file — a genuine anti-fabrication check. For 267317
  it also asserts the pLDDT/length correction fired and is noted.
- **(C)** the review flag fires for **169208**, is absent for **602276**, and
  267317 lands at **attention / not-flagged**.

Result: **85/85 checks pass** across the three reports.

---

## 5. Grounding self-test

`267317_grounding_selftest.txt` records a real paste-test: the 267317 report's
`system_prompt` is used as the system message and its evidence + a user
question ("Describe this protein and explain why the tools disagree… also, what
organism is it from?") as the user turn, answered via `host.llm` (no external
network LLM). Automated checks on the answer confirm it: names GH78 vs GH28 vs
GH92 from the evidence, cites both Pfam domains, uses the corrected pLDDT
(76.6, never 69.2), reports the exact fusion confidence, invents **no** new
family, and refuses the unanswerable organism question with the exact required
phrase.

**Known limitation (noted honestly):** the evidence carries family/EC *values*
but not human-readable *names* for every family/EC, so the LLM supplies names
from general knowledge and can slip (it called EC 3.2.1.40 "exopolygalacturonase"
whereas it is α-L-rhamnosidase). No evidence *value* was fabricated. The clean
fix — carrying canonical family/EC display names in the evidence (e.g. from
dbCAN `fam-substrate-mapping`) — is a future enhancement.

---

## 6. Regeneration commands

From a directory containing `build_ai_report.py` and the unpacked `assets/`
bundle:

```bash
# one report per protein found in the bundle -> reports/
python build_ai_report.py --assets assets --all --outdir reports/

# or one protein at a time
python build_ai_report.py --assets assets --protein 267317 --out 267317_ai_report.json
python build_ai_report.py --assets assets --protein 602276 --out 602276_ai_report.json
python build_ai_report.py --assets assets --protein 169208 --out 169208_ai_report.json

# validate every report against the raw sources
python validate_reports.py assets .
```

`build_ai_report.py` uses only the Python standard library. `validate_reports.py`
likewise. The grounding self-test uses `host.llm` inside the analysis kernel.

## 7. Files

- `build_ai_report.py` — the generator (stdlib only, documented, reusable).
- `validate_reports.py` — independent anti-fabrication + review-flag validator.
- `267317_ai_report.json`, `602276_ai_report.json`, `169208_ai_report.json` — example reports.
- `267317_grounding_selftest.txt` — the host.llm grounding self-test transcript.
- `README_ai_reports.md` — this note.
