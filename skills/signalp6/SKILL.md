---
name: signalp6
description: Predict signal peptides and their cleavage sites with the licensed SignalP-6.0 (`signalp6` command-line tool). Use whenever you need to classify a protein\'s N-terminal signal peptide type - Sec/SPI (standard secretory, "SP"), Sec/SPII (lipoprotein, "LIPO"), Tat/SPI ("TAT") - versus no signal peptide, and to locate the cleavage site; for example confirming that a CAZyme or other protein is classically secreted. Reach for this skill when the user mentions "SignalP", "signal peptide cleavage site", "SPI/SPII/TAT", "lipoprotein signal", or asks "is this secreted" and wants the dedicated SignalP model. IMPORTANT: SignalP-6.0 is license-gated (DTU academic download) and is NOT pip-installable without the license file; if the `signalp6` binary is not on PATH this skill returns status `not_installed` with install guidance and never fabricates a prediction.
---

# SignalP-6.0: signal peptide prediction + cleavage site

SignalP-6.0 (Teufel et al. 2022) detects all five signal-peptide types across all
domains of life and predicts the cleavage site. It is the reference tool for
answering "does this protein have a signal peptide, of what type, and where is it
cleaved?".

## Licensing and installation (scaffold - read before use)

SignalP-6.0 is **license-gated**: the weights and CLI are distributed only via an
academic-use download from DTU, not from PyPI/conda. This skill therefore wraps
the `signalp6` CLI **only when it is already present on PATH**. When it is not, the
helpers return `status="not_installed"` with the install steps below - they never
invent a prediction.

To install (one-time, requires accepting the academic license):
1. Request the download at
   `https://services.healthtech.dtu.dk/services/SignalP-6.0/` (choose the "fast"
   or "slow" model; "fast" is sufficient for most screens).
2. Unpack and `pip install signalp-6-package/` (the bundle ships its own installer),
   then run `signalp6-register <path-to-model-weights>` to place the weights.
3. Confirm `signalp6 --version` works and the binary is on PATH.

## Requirements

- `signalp6` on PATH (verify with the `find_signalp6()` helper).
- No network access needed - SignalP-6.0 runs fully locally once installed.

The kernel sidecar (`kernel.py`) defines three helpers loaded automatically with
this skill: `find_signalp6()`, `run_signalp6(...)`, and `parse_signalp6(...)`.

## Workflow

### 1. Check availability first

```python
binpath = find_signalp6()   # -> "/path/to/signalp6" or None
```

If `None`, report to the user that SignalP-6.0 is not installed (status
`not_installed`) and point them at the install steps above. Do NOT substitute
another tool\'s output as if it were SignalP.

### 2. Run (only if installed)

```python
res = run_signalp6("myprotein.fasta", out_dir="signalp6_out",
                   organism="eukarya", mode="fast", fmt="txt")
# res: {"status": "ok"|"not_installed"|"failed", "binary": ..., "command": ...,
#       "returncode": ..., "stdout_tail": ..., "stderr_tail": ...}
```

`organism` is one of `eukarya`, `other` (bacteria/archaea both SP types), `arch`,
or `gram+`/`gram-` depending on your SignalP build; `mode` is `fast` or `slow`.

### 3. Parse the prediction TSV

```python
parsed = parse_signalp6("signalp6_out")
```

SignalP-6.0 writes `prediction_results.txt` (a TSV). The parser extracts, per
protein id, the predicted type, the probability, and the cleavage site.

## Output schema

`parse_signalp6` returns:

```
{
  "status": "ok" | "no_output_found",
  "n_proteins": int,
  "proteins": {
     "<id>": {
        "prediction":      "SP" | "LIPO" | "TAT" | "TATLIPO" | "PILIN" | "NO_SP",
        "prediction_raw":  "<raw SignalP label; OTHER = no signal peptide>",
        "sp_probability":  float,          # probability of the winning signal-peptide class
        "cleavage_site":   "<pos>" | null, # e.g. "23-24" (between residues 23 and 24)
        "cleavage_prob":   float | null,
        "all_probs":       {<class>: float, ...}   # every column SignalP reported
     }, ...
  },
  "source_files": {"tsv": [...]}
}
```

Prediction mapping: SignalP\'s `OTHER` class means **no signal peptide** and is
normalized to `NO_SP`. `SP` = Sec/SPI (standard secretory), `LIPO` = Sec/SPII
(lipoprotein), `TAT` = Tat/SPI, `TATLIPO` = Tat/SPII, `PILIN` = Sec/SPIII.

## Honesty contract

This is a scaffold for a license-gated tool. If `signalp6` is absent, the only
correct output is `status="not_installed"` plus the install guidance - never a
fabricated SP type, probability, or cleavage coordinate. The parser format above
follows the documented SignalP-6.0 `prediction_results.txt` layout; verify it
against the real file the first time the tool is actually run in your environment,
since column sets vary by `organism` and SignalP build.
