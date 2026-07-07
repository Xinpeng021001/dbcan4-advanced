# Compute environment — met.unl.edu

**Host:** met.unl.edu — 128 CPU, 8× RTX A5500 (24 GB each), 386 GB RAM, CUDA driver 580.126.09.
**Scheduler:** none. Jobs run **directly** on the host (no SLURM/sbatch). Scratch: `/array1/xinpeng/scratch`.

## GPU environment
Path: `/array1/xinpeng/dbcan4-advanced/venv` (Python 3.11, built with `uv`).
- `torch` 2.12.1 (CUDA) — verified `torch.cuda.is_available()` across all 8 GPUs.
- EvolutionaryScale **`esm` 3.2.1** — provides **ESM-C** (`esmc_600m` → 1152-dim embeddings).
- `faiss-cpu`, `scikit-learn`, `biopython`, `h5py`, `pandas`, `numpy`.
- `HF_HOME=/array1/xinpeng/dbcan4-advanced/hf_cache`.

**Gotcha:** do **not** also install `fair-esm` (legacy ESM-2). It clashes with EvolutionaryScale
`esm` on the `esm/` import namespace and breaks `from esm.models.esmc import ESMC`. If it sneaks in
as a transitive dep, `pip uninstall fair-esm` and `pip install --force-reinstall --no-deps esm`.

System binaries already on met: `foldseek`, `diamond`, `TMalign`, `hmmscan` (structure/baseline tiers).

Frozen package list: [`met_env_freeze.txt`](../met_env_freeze.txt).

## Data
`/array1/xinpeng/fungi-cazyme-project/` — CAZy 2024 (`CAZyDB.07142024.fungi.faa`) and
2025 (`CAZyDB.07242025.fungi.faa`) fungal CAZymes; family label is in each FASTA header.
Splits built into `/array1/xinpeng/dbcan4-advanced/data/` by `scripts/build_reference.py`.
