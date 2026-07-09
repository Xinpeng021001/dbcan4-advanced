---
name: deeploc
description: >-
  Predict protein subcellular localization with DeepLoc-2.0 (DTU Health Tech).
  Use this skill whenever the task is to predict where a protein localizes in the
  cell, to ask whether a protein is secreted / extracellular / membrane-bound /
  cytoplasmic / nuclear, or to run DeepLoc on a FASTA. DeepLoc-2.0 is license-gated
  and not pip/conda-installable, so this skill wraps the CLI when present and
  otherwise reports its real install path — and, as a transparent fallback,
  derives a clearly-sourced localization call from a signal-peptide flag plus GO
  cellular-component terms (labelled as derived, never attributed to DeepLoc).
  Pairs with the protein-function skill, which supplies the signal region and GO.
---

# deeploc

Subcellular localization for proteins. DeepLoc-2.0 predicts one of 10
compartments (Cytoplasm, Nucleus, Extracellular, Cell membrane, Mitochondrion,
Plastid, ER, Lysosome/Vacuole, Golgi, Peroxisome) from sequence using
protein-language-model embeddings.

**This is a scaffold skill.** DeepLoc-2.0 is distributed by DTU Health Tech under
an academic license and is **not** on PyPI or conda, so it cannot be installed
non-interactively. The `kernel.py` sidecar wraps the CLI when it is on PATH and,
when it is not, returns an honest `not_installed` status with the real install
instructions — it never fabricates a localization.

## Workflow

### 1. Try DeepLoc directly

```python
run = run_deeploc("proteins.fasta", out_dir="deeploc_out", model="Accurate")
run["status"]   # 'success' | 'not_installed' | 'failed'
```

- `success` → `run["predictions"]` holds the parsed per-protein rows from
  DeepLoc's `results_*.csv` (localization + per-class probabilities).
- `not_installed` → `run["install_doc"]` has the DTU download/install steps
  (register at the DTU DeepLoc-2.0 service, `pip install` the downloaded tarball,
  verify with `deeploc2 --help`; the CLI is `deeploc2 -f in.fasta -o out -m Accurate`).
- `failed` → `run["error"]` carries the real stderr.

Report the status truthfully. Do not present a `not_installed` result as if a
prediction had been made.

### 2. Derived fallback (transparent, not DeepLoc)

When DeepLoc is unavailable you can still make a *sourced* localization call from
evidence you already have — an N-terminal signal peptide and GO
cellular-component terms:

```python
loc = derive_localization_from_evidence(
    signal_region=sig,     # from protein-function's signal_region()
    go_terms=go,           # normalised GO dicts; only cellular_component used
    cazy_family="GH78",
)
loc["localization_call"], loc["confidence"], loc["evidence"]
```

Logic (fully transparent, in `loc["evidence"]`): GO cellular-component terms decide
when present; otherwise a predicted signal peptide implies secretory routing
(Extracellular for a secreted CAZyme). The result is stamped
`method: "derived from signal-peptide + GO-CC (NOT DeepLoc)"` — always keep that
label so the provenance is unambiguous. This is a heuristic prior, weaker than a
DeepLoc run, and should be reported as such.

### 3. Save

Save predictions/derived call as JSON and `save_artifacts`. In your summary,
state clearly which route produced the call (DeepLoc vs derived) and its confidence.

## Honesty

The whole point of the scaffold is that the *absence* of a tool is reported, not
papered over. A `not_installed` status plus a clearly-labelled derived call is a
valid, honest deliverable. An invented DeepLoc probability is not.
