# Installation

## Requirements

### Hardware

- **GPU** (NVIDIA, ≥16 GB) for the ESM-C embedding and ESMFold structure steps. The reference
  deployment uses 8× RTX A5500 (24 GB). CPU-only works for the stub DAG and baseline dbCAN, but
  not for the real advanced/structure tiers.
- ~30 GB disk for model weights + databases.

### Software

| Component | Version | Notes |
|---|---|---|
| Python | ≥3.10 (3.11 on the reference host) | |
| PyTorch | ≥2.0 (CUDA build) | `torch.cuda.is_available()` must be `True` |
| EvolutionaryScale `esm` | 3.2.1 | provides **ESM-C** (`esmc_600m` → 1152-dim) |
| `run_dbcan` | V5 (dev 5.0.7 on the reference host) | baseline HMMER/dbCAN_sub/DIAMOND |
| Nextflow | ≥24 (needs Java 17+) | only for `dbcan4 run` |
| foldseek, diamond, hmmscan | on `PATH` | structure/baseline tiers |
| FastAPI/uvicorn/SQLAlchemy/Alembic | vendored in-repo (`src/bioforge`) | web stack |

!!! danger "Do NOT install `fair-esm`"
    `fair-esm` (legacy ESM-2) clashes with EvolutionaryScale `esm` on the `esm/` import
    namespace and breaks `from esm.models.esmc import ESMC`. If it sneaks in as a transitive
    dependency: `pip uninstall fair-esm` then `pip install --force-reinstall --no-deps esm`.

License-gated (optional, handled as honest fallbacks if absent): **SignalP-6.0**, **DeepLoc-2.0**
(DTU academic download). **CLEAN** and **DeepTMHMM** run in isolated environments.

## Data assets

These are **not** in the git repo (`.gitignore` excludes all `*.npz/*.pt/*.hmm/*.dmnd` and
`data/`). On a fresh machine you must **copy them from the reference host or rebuild them**:

| Asset | Path (reference host) | How to obtain |
|---|---|---|
| ESM-C reference index | `emb/ref2024.shard{0..7}.npz` | copy, or rebuild: `build_reference.py` → `embed_esmc.py` |
| Trained heads | `results/heads/{heads.pt,proj_ref.npz}` | copy, or rebuild: `train_heads.py` |
| dbCAN database (~7.4 GB) | `/array1/xinpeng/dbcan_db` | `run_dbcan database --db_dir dbcan_db` |
| Pfam-A (~2.2 GB, pressed) | `/array1/xinpeng/pfam/Pfam-A.hmm` | download + `hmmpress` |
| ESM-C / ESMFold weights (~16 GB) | `hf_cache/` | auto-download to `HF_HOME` on first use |

## On the reference host (met) — nothing to install

All assets live under `$REPO=/array1/xinpeng/dbcan4-advanced`. Verify with:

```bash
source /array1/xinpeng/scratch/biodb_venv/bin/activate   # or any venv where you ran 'pip install -e .'
dbcan4 info      # prints resolved pipeline / reference-index / heads paths
```

Expected output resolves every path (`heads.pt … (exists)`). If it does, skip to the
[Quick start](quickstart.md).

## From scratch on a new GPU machine

```bash
# 1. code
git clone https://github.com/Xinpeng021001/dbcan4-advanced.git
cd dbcan4-advanced

# 2. engine venv (torch + EvolutionaryScale esm 3.2.1 = ESM-C). Do NOT install fair-esm.
python -m venv venv && source venv/bin/activate
pip install torch --index-url https://download.pytorch.org/whl/cu121   # match your CUDA
pip install "esm==3.2.1" faiss-cpu scikit-learn biopython h5py pandas numpy
pip install -e .            # installs the `dbcan4` console script + the vendored bioforge
                            # web layer (bioforge-ingest, bioforge-ingest-advanced, web app)

# 3. data assets — copy from the reference host (fastest) or rebuild
#    scp -r met:/array1/xinpeng/dbcan4-advanced/emb ./emb
#    scp -r met:/array1/xinpeng/dbcan4-advanced/results/heads ./results/heads
#    run_dbcan database --db_dir dbcan_db          # ~7.4 GB

# 4. verify
dbcan4 info
```

## Overriding asset locations

The engine resolves assets in this order: `--assets` flag / `DBCAN4_ASSETS` → bundled `nf/` →
repo checkout → reference-host defaults. Heavy data paths are overridable via environment
variables:

| Variable | Overrides |
|---|---|
| `DBCAN4_ASSETS` | asset root (dir with `nf/` + `results/`) |
| `DBCAN4_REF_EMB` | reference embedding index prefix |
| `DBCAN4_HEADS` | trained heads `.pt` |
| `DBCAN4_PROJ_REF` | projected reference `.npz` |
| `DBCAN4_ENGINE_PYTHON` | the Python that runs `embed`/`infer` (must have torch + esm) |
