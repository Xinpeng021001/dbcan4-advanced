---
name: esmfold-fold
description: Predict a protein 3D structure from sequence with ESMFold (facebook/esmfold_v1 via HuggingFace transformers). Use when you need a fast single-sequence predicted fold (PDB + per-residue pLDDT) for a protein or a CAZyme, e.g. to feed a structural-homology search (Foldseek) or a web-UI 3D viewer. Runs on a CUDA GPU. Not for complexes or nucleic acids (use Boltz/Chai for those).
---

# esmfold-fold

Single-sequence structure prediction with **ESMFold** (`facebook/esmfold_v1`).
Emits a PDB whose CA B-factor column holds per-residue **pLDDT** (0-100).

## When to use
- You have a protein sequence and want a predicted 3D fold fast (seconds-minutes, no MSA).
- Downstream: Foldseek structural search, 3D viewer, pocket/active-site inspection.

## Requirements
- CUDA GPU (~16 GB+; a 24 GB card folds up to ~1000-1400 aa with chunking).
- `transformers`, `torch`, `accelerate` in the env. First run downloads ~2.8 GB weights to `HF_HOME`.

## Gotchas (verified on met RTX A5500, 24 GB)
1. **Sanitize the sequence**: ESMFold's tokenizer throws on non-standard residues; map anything outside the 20 canonical aa (especially the `*` stop symbol) to `X`. `sanitize_sequence()` does this.
2. **One model per GPU**: two ESMFold models will not fit one 24 GB card — serialize folds, `torch.cuda.empty_cache()` between them.
3. **Long sequences**: `model.trunk.set_chunk_size(64)` and half-precision the language model (`model.esm = model.esm.half()`) to fit long fungal CAZymes; cost is O(L^2) (~1 s for 80 aa, up to ~90 s for 600-900 aa).
4. Load time ~90 s (weights download on first run).

## Workflow
`kernel.py` provides:
- `sanitize_sequence(seq)` -> uppercase, strip `*`, non-canonical -> `X`.
- `esmfold_script(fasta_name, out_dir)` -> returns a self-contained Python script (string) that loads ESMFold once and folds every record in `fasta_name`, writing `<out_dir>/<id>.pdb` + `fold_summary.json`. Submit this as a GPU job.
- `plddt_from_pdb(pdb_path)` -> `(n_residues, mean_pLDDT, min, max)` from CA B-factors.

Typical remote run (SSH host, no scheduler), from the `repl` tool:
```python
from pathlib import Path
Path("fold.py").write_text(esmfold_script("in.fasta", "pdb_out"))
c = host.compute.create("ssh:<host>")
job = c.submit_job(
    intent="ESMFold: fold N proteins on 1 GPU",
    command=("export HF_HOME=/path/to/hf_cache\n"
             "export CUDA_VISIBLE_DEVICES=$(nvidia-smi --query-gpu=index,memory.free "
             "--format=csv,noheader,nounits | sort -t, -k2 -nr | head -1 | cut -d, -f1)\n"
             "<venv>/bin/python fold.py\ncp pdb_out/*.pdb ./ 2>/dev/null || true\nls -lh *.pdb"),
    inputs=[{"src": "in.fasta", "dst_filename": "in.fasta"},
            {"src": "fold.py", "dst_filename": "fold.py"}],
    outputs=[{"glob": "*.pdb", "visibility": "featured"},
             {"glob": "fold_summary.json", "visibility": "featured"}],
    timeout_seconds=3600)
```
Then park on `wait_for_notification`; `save_artifacts(payload["featured_files"])`.

## Interpreting pLDDT
Per-residue confidence in the B-factor column: >90 very high, 70-90 confident,
50-70 low, <50 very low (often disordered/linker). Colour the viewer by it.
ESMFold is single-sequence, so pLDDT is typically a bit below AlphaFold2-with-MSA
on the same target; treat <50 regions as unreliable geometry.
