---
name: deeptmhmm
description: Predict transmembrane topology AND signal peptides for any protein sequence with DeepTMHMM (DTU/DeepTMHMM), run via pybiolib on the BioLib cloud. Use this whenever you need to know whether a protein is membrane-embedded, how many transmembrane helices it has, its inside/outside membrane topology, or whether it carries an N-terminal signal peptide - for example deciding whether a CAZyme or other protein is secreted, membrane-anchored, or cytosolic. Make sure to reach for this skill whenever the user mentions "transmembrane", "TM helices", "membrane topology", "DeepTMHMM", asks "is this protein secreted / membrane-bound / cytosolic", or hands you a FASTA together with any membrane-topology question, even if they do not name DeepTMHMM explicitly.
---

# DeepTMHMM: transmembrane topology + signal peptide prediction

DeepTMHMM (Hallgren et al. 2022) is a deep-learning predictor that labels every
residue of a protein as signal peptide, transmembrane helix (or beta strand),
inside (cytoplasmic), or outside (extracellular/luminal). It classifies each
protein into one of six categories and returns a full per-residue topology, which
is exactly what you need to decide whether a protein is secreted, membrane-bound,
or globular/cytosolic.

This skill runs DeepTMHMM through **pybiolib**, which dispatches the job to the
DeepTMHMM app on the BioLib cloud - so you do not install the model or its
weights locally.

## Requirements

- `pybiolib` installed in the active environment: `pip install pybiolib`.
- Outbound network access to `biolib.com` (the job runs on BioLib's servers).
- Public DeepTMHMM runs work anonymously in most cases. If BioLib requires a
  signed-in account for compute, the run raises an auth error - the helper
  captures the real error text and returns `status="failed"`; do NOT fabricate a
  result. A token can be supplied out-of-band via `biolib login` / the
  `BIOLIB_TOKEN` environment variable.

The kernel sidecar (`kernel.py`) defines two helpers that load automatically when
this skill is loaded: `run_deeptmhmm(...)` and `parse_deeptmhmm(...)`.

## Workflow

### 1. Run the prediction

```python
res = run_deeptmhmm("myprotein.fasta", workdir="deeptmhmm_out")
# res: {"status": "ok"|"not_installed"|"failed", "out_dir": ..., "exit_code": ...,
#       "pybiolib_version": ..., "error"/"traceback" on failure}
```

`run_deeptmhmm` accepts any FASTA (one or many sequences) and saves DeepTMHMM's
output files into `workdir`. Check `res["status"]` before parsing:
- `"ok"` - files saved, proceed to parse.
- `"not_installed"` - pybiolib is not importable; `res["error"]` says how to install.
- `"failed"` - the run raised (e.g. auth/network); `res["error"]`/`res["traceback"]`
  carry the real message. Report it truthfully.

### 2. Parse the output

```python
parsed = parse_deeptmhmm("deeptmhmm_out")
```

DeepTMHMM writes `predicted_topologies.3line` (per-residue labels) and
`TMRs.gff3` (region spans + a "Number of predicted TMRs" count). The parser reads
both and returns a dict keyed by protein id.

## Output schema

`parse_deeptmhmm` returns:

```
{
  "status": "ok" | "no_3line_found",
  "n_proteins": int,
  "proteins": {
     "<id>": {
        "prediction":        "TM" | "SP+TM" | "SP" | "Globular" | "TM_beta" | "SP+TM_beta",
        "prediction_raw":    "<raw DeepTMHMM label: TM/SP+TM/SP/GLOB/BETA/SP+BETA>",
        "n_tm_helices":      int,          # authoritative count (from GFF3 if present)
        "topology_string":   "SSSS...IIII...MMMM...OOOO",
        "has_signal_peptide": bool,
        "signal_peptide_span": [start, end] | null,   # 1-based, inclusive
        "length":            int,
        "regions": [ {"type": "signal_peptide"|"inside"|"outside"|"TMhelix"|"beta_strand"|"periplasm",
                      "start": int, "end": int}, ... ]   # 1-based, inclusive, contiguous
     }, ...
  },
  "source_files": {"three_line": [...], "gff3": [...]}
}
```

Per-residue topology label meanings: `S` signal peptide, `I` inside
(cytoplasmic), `O` outside, `M` alpha-helical TM, `B` beta-strand TM, `P`
periplasm (beta-barrel loops).

## Interpreting the prediction categories

- **Globular (`GLOB`)** - soluble protein, no TM segments, no signal peptide.
  Cytosolic unless secreted by a non-classical route.
- **SP** - has an N-terminal signal peptide but no TM segment: typically a
  **secreted** protein (the mature chain after cleavage is extracellular/luminal).
- **TM** - one or more transmembrane segments, no signal peptide: an integral
  membrane protein; the topology string gives inside/outside orientation.
- **SP+TM** - signal peptide plus TM segment(s): membrane protein routed through
  the secretory pathway, or a signal-anchored membrane protein.
- **BETA / SP+BETA** - beta-barrel membrane protein (outer-membrane type).

For CAZymes and other secreted enzymes, `SP` (or `has_signal_peptide=True` with
`n_tm_helices=0`) is the signature of a classically secreted enzyme; `n_tm_helices>=1`
indicates membrane association (e.g. a C-terminal anchor).

## Notes

- Keep the skill general: it works for any protein FASTA, single or multi-record.
- DeepTMHMM predicts topology; it does not annotate catalytic function. Combine
  with family/domain evidence for a full picture.
- If you need only signal-peptide type + cleavage site with the specialized
  SignalP model, see the `signalp6` skill (licensed, local CLI).
