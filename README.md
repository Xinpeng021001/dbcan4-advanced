# dbCAN4-advanced

**Protein-language-model + structure-similarity CAZyme annotation for fungal proteins — beyond HMMER/DIAMOND.**

[![Docs](https://img.shields.io/badge/docs-mkdocs--material-blue)](docs/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/nextflow-DSL2-brightgreen)](https://www.nextflow.io/)

> **One command turns a protein FASTA into a browsable CAZyme-annotation database.**
>
> ```bash
> bash dbcan4_workup.sh proteins.faa --serve      # → http://127.0.0.1:8000
> ```

Current dbCAN (`run_dbcan` / dbCAN3) assigns CAZy families by **sequence similarity**
(HMMER, dbCAN_sub, DIAMOND). This misses **remote-homolog CAZymes** — enzymes that share
fold, mechanism, and active-site geometry with known families but have drifted below the
sequence-identity detection threshold. dbCAN4-advanced adds an orthogonal
**protein-language-model (ESM-C) + structure (ESMFold/Foldseek/CAZyme3D)** tier, a
consensus **fusion** layer with an abstain option, and a full per-protein functional
workup — then ingests everything into a versioned database with a per-gene deep-dive web UI.

![Product architecture](docs/architecture_product.png)

---

## Table of contents

1. [What you get](#what-you-get)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Quick start](#quick-start)
5. [Detailed usage](#detailed-usage)
6. [The CLI](#the-cli)
7. [Understanding the output](#understanding-the-output)
8. [The web UI](#the-web-ui)
9. [Results & benchmarks](#results--benchmarks)
10. [Worked example — hero protein 267317](#worked-example--hero-protein-267317)
11. [Troubleshooting](#troubleshooting)
12. [Prototypes](#prototypes)
13. [Documentation](#documentation)
14. [Citing the methods](#citing-the-methods)
15. [License](#license)

---

## What you get

dbCAN4-advanced is **three parts that compose into one product**:

| Part | What it is | Where |
|---|---|---|
| **`dbcan4` Python package + CLI** | the annotation engine (`embed` / `infer` / `annotate` / `run`) | `src/dbcan4_advanced/`, `pyproject.toml` |
| **Nextflow pipeline** | baseline dbCAN + advanced ESM-C/structure tiers + 8 feature tracks → a standard v1.1 output contract | `nf/` |
| **BioForge database + web UI** | versioned SQLite schema + FastAPI web app that ingests the contract and serves per-gene deep-dive pages | [`github.com/Xinpeng021001/biodb`](https://github.com/Xinpeng021001/biodb) (branch `feature/advanced-cazyme-integration`) |

**Per-protein functional workup** (8 feature tracks): Pfam domains (hmmscan), EC number
(CLEAN), TM topology + signal peptide (DeepTMHMM), 3D structure (ESMFold), subcellular
localization, physicochemistry (Biopython), plus structural-homology hits (Foldseek vs CAZyme3D).

### Honest positioning

On a 2024→2025 fungal temporal holdout the method is **competitive with, not superior to,
sequence baselines**, and the docs say so: contrastive-kNN on frozen ESM-C **ties DIAMOND at
family level** (overlap 0.973 vs 0.981) and beats dbCAN-sub on subfamily. The eval mass sits
at high identity (median ~81% to the 2024 reference), so it measures near-term annotation, not
twilight-zone remote-homolog recovery. **The value is orthogonality, calibrated abstention,
subfamily resolution, and fungal calibration** — plus a rigorous DB-vintage leakage control
most tool comparisons skip. Full detail: [`docs/benchmark_report.md`](docs/benchmark_report.md).

---

## Requirements

### Hardware
- **GPU** (NVIDIA, ≥16 GB) for the ESM-C embedding and ESMFold structure steps. The reference
  deployment uses 8× RTX A5500 (24 GB). CPU-only works for the stub DAG and for baseline dbCAN,
  but not for the real advanced/structure tiers.
- ~30 GB disk for model weights + databases (see [data assets](#3-data-assets)).

### Software
| Component | Version | Notes |
|---|---|---|
| Python | ≥3.10 (3.11 on the reference host) | |
| PyTorch | ≥2.0 (CUDA build) | `torch.cuda.is_available()` must be `True` |
| EvolutionaryScale `esm` | 3.2.1 | provides **ESM-C** (`esmc_600m` → 1152-dim). **Do NOT also install `fair-esm`** — it clashes on the `esm/` namespace |
| `run_dbcan` | V5 (dev 5.0.7 on the reference host) | baseline HMMER/dbCAN_sub/DIAMOND |
| Nextflow | ≥24 (needs Java 17+) | only for `dbcan4 run` |
| foldseek, diamond, hmmscan | on `PATH` | structure/baseline tiers |
| FastAPI/uvicorn/SQLAlchemy/Alembic | via `biodb` | web stack |

License-gated (optional, handled as honest fallbacks if absent): **SignalP-6.0**, **DeepLoc-2.0**
(DTU academic download). **CLEAN** and **DeepTMHMM** run in isolated environments (see below).

### 3. Data assets

These are **not** in the git repo (`.gitignore` excludes all `*.npz/*.pt/*.hmm/*.dmnd` and `data/`).
On a fresh machine you must **copy them from the reference host or rebuild them**:

| Asset | Path (reference host) | How to obtain |
|---|---|---|
| ESM-C reference index | `emb/ref2024.shard{0..7}.npz` | copy, or rebuild: `build_reference.py` → `embed_esmc.py` |
| Trained heads | `results/heads/{heads.pt,proj_ref.npz}` | copy, or rebuild: `train_heads.py` |
| dbCAN database (~7.4 GB) | `/array1/xinpeng/dbcan_db` | `run_dbcan database --db_dir dbcan_db` |
| Pfam-A (~2.2 GB, pressed) | `/array1/xinpeng/pfam/Pfam-A.hmm` | download + `hmmpress` |
| ESM-C / ESMFold weights (~16 GB) | `hf_cache/` | auto-download to `HF_HOME` on first use |

> **On the reference host (`met.unl.edu`) everything above is already installed** — skip straight to [Quick start](#quick-start).

---

## Installation

### On the reference host (met) — nothing to install

All assets live under `$REPO=/array1/xinpeng/dbcan4-advanced`. Verify with:

```bash
source /array1/xinpeng/scratch/biodb_venv/bin/activate
dbcan4 info      # prints resolved pipeline / reference-index / heads paths
```

### From scratch on a new GPU machine

```bash
# 1. code
git clone https://github.com/Xinpeng021001/dbcan4-advanced.git
cd dbcan4-advanced
git clone -b feature/advanced-cazyme-integration https://github.com/Xinpeng021001/biodb.git

# 2. engine venv (torch + EvolutionaryScale esm 3.2.1 = ESM-C). Do NOT install fair-esm.
python -m venv venv && source venv/bin/activate
pip install torch --index-url https://download.pytorch.org/whl/cu121   # match your CUDA
pip install "esm==3.2.1" faiss-cpu scikit-learn biopython h5py pandas numpy
pip install -e .            # installs the `dbcan4` console script
pip install -e biodb        # installs bioforge-ingest, bioforge-ingest-advanced, web app

# 3. data assets — copy from the reference host (fastest) or rebuild (see table above)
#    scp -r met:/array1/xinpeng/dbcan4-advanced/emb ./emb
#    scp -r met:/array1/xinpeng/dbcan4-advanced/results/heads ./results/heads
#    run_dbcan database --db_dir dbcan_db          # ~7.4 GB

# 4. verify
dbcan4 info
```

Point the engine at non-default asset locations with `DBCAN4_REF_EMB`, `DBCAN4_HEADS`,
`DBCAN4_PROJ_REF`, `DBCAN4_ENGINE_PYTHON`, or `--assets`.

---

## Quick start

Three paths, fastest first. All use the shipped 3-protein example `examples/real3.faa`.

### A. Label-free family calls (real ESM-C engine, ~1–2 min on GPU)

```bash
source /array1/xinpeng/scratch/biodb_venv/bin/activate      # reference host
CUDA_VISIBLE_DEVICES=0 dbcan4 annotate examples/real3.faa --outdir calls_out
# → calls_out/ESM-C-{kNN,centroid,contrastive}.raw.tsv
```

### B. Full DAG with no GPU/tools (stub — proves the pipeline anywhere, ~1 min)

```bash
source /array1/xinpeng/scratch/bin/nxf_env.sh               # Nextflow + Java
dbcan4 run --fasta examples/real3.faa --sample smoke --outdir stub_out --profile stub --stub
# → stub_out/cazyme_advanced/manifest.json  (contract v1.1: 6 methods + 8 feature types)
```

### C. The whole product — real workup + web UI (~10–15 min for 3 proteins)

```bash
bash dbcan4_workup.sh examples/real3.faa --serve --gpu 0
# → runs baseline + advanced + all 8 feature tracks, ingests, serves http://127.0.0.1:8000
```

To view the UI from your laptop (uvicorn binds `127.0.0.1` on the server):

```bash
ssh -L 8000:127.0.0.1:8000 xinpeng@met.unl.edu     # then open http://localhost:8000
```

---

## Detailed usage

### The one-command workup: `dbcan4_workup.sh`

This is the **verified real path** — it runs each proven tool directly (handling every
tool-specific gotcha internally) and is the recommended way to annotate real sequences.

```bash
bash dbcan4_workup.sh <proteins.faa> [options]
```

It runs, in order:

1. **Baseline dbCAN** (`run_dbcan` V5 → HMMER · dbCAN_sub · DIAMOND)
2. **Advanced ESM-C** embed + label-free infer (kNN / centroid / contrastive)
3. **Fusion** consensus
4. **Pfam domains** (hmmscan vs Pfam-A, `--cut_ga`)
5. **Physicochemistry** (Biopython: MW, pI, GRAVY, N-glyc sequons)
6. **DeepTMHMM** TM topology + signal peptide (+ derived localization)
7. **CLEAN** EC number
   7b. **ESMFold** 3D structure (local GPU)
8. **v1.1 manifest** → (with `--serve`) ingest + launch the web UI

| Option | Meaning |
|---|---|
| `--outdir DIR` | output directory (default `./<sample>_workup`) |
| `--sample NAME` | sample key (default: FASTA basename) |
| `--serve` | ingest into BioForge SQLite + launch the web UI |
| `--port N` | web UI port (default 8000) |
| `--gpu N` | CUDA device index (default 0) |
| `--gff FILE` | optional genomic GFF (adds a genome track; **omit for protein-input mode**) |
| `--no-deeptmhmm` | skip DeepTMHMM (e.g. offline — it uses the BioLib cloud) |
| `--no-clean` | skip CLEAN EC (CPU-slow, ~2 min/protein) |
| `--no-structure` | skip ESMFold folding (e.g. no GPU) |

> **dbCAN4 is fungal + protein-input**: genes are built straight from the protein FASTA — **no
> genome, no Prokka, no GFF**. Each protein becomes its own gene with residue coordinates `1–L`.
> Pass `--gff` only when you genuinely have genomic coordinates.

**Fast smoke variant** (skips the slow/network steps — proves the whole chain except those 3 tracks):

```bash
bash dbcan4_workup.sh examples/real3.faa --serve --gpu 0 --no-deeptmhmm --no-clean --no-structure
```

### Run on your own data

Just point it at any fungal protein FASTA:

```bash
bash dbcan4_workup.sh /path/to/my_proteome.faa --sample my_fungus --serve --gpu 0
```

---

## The CLI

```
dbcan4 embed     FASTA -> ESM-C embeddings (.npz)                 [GPU]
dbcan4 infer     embeddings -> label-free family calls (TSVs)     [CPU/GPU]
dbcan4 annotate  FASTA -> family calls in one step (embed+infer)  [GPU]
dbcan4 run       FASTA -> full Nextflow pipeline (+ --serve to ingest & launch web UI)
dbcan4 info      show resolved asset paths + versions
```

- `dbcan4 annotate` = `embed` + `infer`; emits `ESM-C-{kNN,centroid,contrastive}.raw.tsv`.
  Label-free — needs only the precomputed reference index + trained heads, **no ground-truth labels**.
- `dbcan4 run` wraps the Nextflow pipeline. `-profile stub -stub-run` (`--profile stub --stub`)
  proves the whole DAG with no GPU/tools; `--serve` chains `alembic upgrade → bioforge-ingest →
  bioforge-ingest-advanced → uvicorn`.
- The heavy GPU steps shell out to the engine venv (torch + esm) via `DBCAN4_ENGINE_PYTHON`, so
  the CLI itself can be installed in a lightweight venv.

---

## Understanding the output

The pipeline publishes a **standardized, versioned layout** — one `manifest.json` plus a
funcscan-style tree. Downstream consumers (the BioForge ingester) read only this contract, never
a tool's raw output. Full schema: [`nf/OUTPUT_CONTRACT.md`](nf/OUTPUT_CONTRACT.md).

```
<outdir>/cazyme_advanced/
  manifest.json                      # the contract entry point
  predictions/<sample>/
    ESM-C-kNN.tsv  ESM-C-centroid.tsv  ESM-C-contrastive.tsv
    Foldseek-CAZyme3D.tsv  SaProt.tsv  fusion.tsv
  features/<sample>/
    domains.tsv  ec_prediction.tsv  deeptmhmm.tsv  signalp6.tsv
    localization.tsv  physicochem.tsv  structures.tsv
    structures/<protein_id>.pdb
```

Every prediction TSV shares the same normalized schema
(`protein_id`, `family`, `confidence`, `ec`, `all_families`, `extra`), and each `manifest.json`
becomes one **versioned BioForge release** loaded additively next to the baseline — so
"advanced vs baseline" is a query across releases, never a mutation.

![Feature coverage](docs/fig_feature_coverage.png)

---

## The web UI

Each gene page is a genuine deep-dive. Every card is populated from real data.

**Dashboard & browse:**

![Dashboard](docs/ui/ui_dashboard.png)
![Browse](docs/ui/ui_browse.png)

**Per-gene deep-dive** (hero protein 267317) — advanced-vs-baseline family calls, a per-residue
DeepTMHMM topology track, the ESMFold model in an interactive **3Dmol** viewer coloured by pLDDT,
plus EC / localization / physicochemistry / Pfam-domain / GO cards:

![Hero gene page 267317](docs/ui/ui_hero_267317.png)

> The screenshots are **true-browser captures** (headless Chrome against the live server — real
> CSS + JavaScript + the 3Dmol WebGL viewer), produced by `capture_ui.sh`.

---

## Results & benchmarks

Fungal 2024→2025 temporal holdout. Full report: [`docs/benchmark_report.md`](docs/benchmark_report.md).

![Master benchmark](docs/figures/master_benchmark_figure.png)

| Method | Known (subfamily) | Known (parent) | Genuinely novel (parent) |
|---|---|---|---|
| DIAMOND (seq, 2024) | 0.98 | 0.98 | 0.33 |
| Contrastive kNN (pLM) | 0.97 | 0.98 | 0.33 |
| Foldseek (structure)* | 0.72 | 0.75 | 0.00 |
| **FUSION (consensus)** | **0.98** | **0.98** | **0.33** |

\*on its 347-structure subset. Known-family recall is effectively tied across DIAMOND, the
trained pLM head, and fusion; fusion's differentiator is **novelty handling** — its abstention
flags genuinely-novel CAZymes at 66.7% vs 0.3% for new-to-fungi families.

**DB-vintage leakage control** (the standout result) — DIAMOND against the *current* dbCAN DB is
contaminated by 2025 eval sequences (novel-to-fungi subfamily recall 0.001 fair-2024 → 0.992
current DB); pLM/structure/fusion tiers train on 2024 only and are already clean:

![Fair benchmark vs 2024 DB](docs/figures/fair_benchmark_2024db.png)

**Domain-level retrieval** fixes the multidomain collapse (whole-protein exact-set 0.006 → 0.412):

![Domain retrieval](docs/figures/domain_retrieval_multidomain.png)

---

## Worked example — hero protein 267317

A 1,088-residue multi-domain fungal α-L-rhamnosidase (truth **GH78**, with an N-terminal GH28
domain) — the ideal case for why comprehensive multi-tool annotation matters:

| head | call | note |
|---|---|---|
| ESM-C-kNN | **GH78** (0.995) ✓ | keys on the rhamnosidase domain |
| ESM-C-centroid | GH92 (0.977) ✗ | dissents |
| ESM-C-contrastive | **GH78** (0.950) ✓ | |
| **fusion** | **GH78** ✓ | majority-correct |
| baseline HMMER/DIAMOND | GH28 | keys on the strong N-terminal Glyco_hydro_28 domain |

Real Pfam architecture: Glyco_hydro_28 (PF00295) + Bac_rhamnosid6H (PF17389); CLEAN EC 3.2.1.40;
DeepTMHMM secreted (SP cleave@18); ESMFold mean pLDDT 76.6. Full record in
[`examples/267317_comprehensive/`](examples/267317_comprehensive/).

The comprehensive poster below combines **all** signals in one view — the multi-tool CAZyme
evidence and function/secretion/domain/localization/physicochemistry cards, the real ESMFold **3D
structure** (3Dmol cartoon coloured by per-residue pLDDT, confident blue core → low-confidence red
termini — the two-domain GH28 + GH78 architecture is clearly visible), **and** the full 1,088-aa
protein sequence:

![Hero 267317 comprehensive poster](examples/267317_comprehensive/comprehensive_267317_full.png)

Regenerate it (renders the structure via 3Dmol in headless Chrome, then captures the full page):

```bash
cd examples/267317_comprehensive
python build_comprehensive_poster.py --html . --outdir out
# → out/comprehensive_267317_full.png  (cards + 3D structure + sequence)
#   out/structure_267317.png           (standalone ESMFold cartoon, pLDDT-coloured)
```

### Structure colouring schemes

`render_structure_views.py` renders the ESMFold model in **four colouring schemes** — and includes a
reusable **colour-by-domain** mode driven by Pfam/CAZyme boundaries:

| mode | what it shows |
|---|---|
| `plddt` | per-residue ESMFold confidence (blue = confident → red = low) |
| `spectrum` | sequence position (N-terminus blue → C-terminus red) |
| `sstruc` | secondary structure (helix / sheet / coil) |
| `domain` | **Pfam / CAZyme domains** — each domain a distinct colour, linker grey |

![267317 structure colouring schemes](examples/267317_comprehensive/structure_views_267317.png)

The **domain** view maps the annotation onto the fold: the N-terminal **GH28** (Glyco_hydro_28,
PF00295, res 62–412) is the right-handed β-helix; the C-terminal **GH78** (Bac_rhamnosid6H, PF17389,
res 667–886) is the α-helical rhamnosidase domain — the two domains that make the baseline (GH28) and
advanced ESM-C (GH78) calls disagree, now spatially distinct on one structure:

![267317 coloured by Pfam/CAZyme domain](examples/267317_comprehensive/structure_domain_267317.png)

```bash
# all four modes + composite panels
python render_structure_views.py --pdb 267317.pdb --outdir out
# colour any structure by your own domains (JSON: [{name,family,start,end,color}, ...])
python render_structure_views.py --pdb my.pdb --domains 267317_domains.json --modes plddt,domain
```

The reusable primitive is `domain_cartoon_style_js(domains)` — hand it Pfam hmmscan hits (each mapped
to a CAZy family + a colour) and it emits the 3Dmol `setStyle()` calls that light up each domain on
the cartoon, so the same function colours any protein by its CAZyme/Pfam architecture.

**Interactive viewer** — the static panels above are a preview; for exploring the structure yourself,
`build_structure_viewer.py` emits a single self-contained HTML file (3Dmol.js embedded inline, works
offline) with a **colour-scheme drop-down**. Open it in any browser, drag to rotate, and switch
between the four schemes live — no need to display every option at once:

```bash
python build_structure_viewer.py --pdb 267317.pdb --out structure_viewer_267317.html
# open structure_viewer_267317.html in a browser → dropdown: domain / pLDDT / spectrum / sstruc
```

The same generator colours any structure by your own domain spec:
`python build_structure_viewer.py --pdb my.pdb --domains my_domains.json`.

---

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `dbcan4 info` shows `MISSING` for heads | data assets not present — copy from the reference host or set `DBCAN4_HEADS`/`DBCAN4_REF_EMB` |
| `from esm.models.esmc import ESMC` fails | `fair-esm` is installed and clashing — `pip uninstall fair-esm` then `pip install --force-reinstall --no-deps esm` |
| DeepTMHMM `FileNotFoundError: hash/…` | biolib stages by **basename** — the workup `cd`s into the FASTA dir and passes a bare filename; use `--no-deeptmhmm` offline |
| CLEAN errors on a `*` char / very long protein | ESM-1b has no `*` token and caps at 1022 aa — `run_clean.sh` strips stops + truncates; use `--no-clean` to skip |
| `dbcan4 run -profile met` (real Nextflow) not fully working | **use `dbcan4_workup.sh` for real runs** — the full-Nextflow *real* path is not proven end-to-end; only `-profile stub` is validated |
| Nextflow "requires Java 17+" | system Java is 11 — `source /array1/xinpeng/scratch/bin/nxf_env.sh` (portable Temurin 21) |
| Web UI blank / no 3D structure | pass `--structures-dir` on ingest (the workup does this); the 3Dmol viewer needs the served PDB |

More gotchas: [`REPRODUCE_PRODUCT.md`](REPRODUCE_PRODUCT.md), [`nf/TOOLS.md`](nf/TOOLS.md).

---

## Prototypes

Two standalone, review-stage prototypes live in [`prototypes/`](prototypes/README.md). Each runs
end-to-end from a small checked-in `example_data/` bundle (2.1 MB) — no GPU, no reference index, no
live compute, no network:

- **Visualization** ([`prototypes/viz/`](prototypes/viz/)) — an interactive embedding explorer
  (raw-ESM-C ⇄ trained-head UMAP toggle over 9,815 reference proteins / 814 families), a training
  dashboard from the real per-epoch head-training log, and a calibration + confusion explorer. One
  command regenerates all of it: `python build_visualizations.py --assets ../example_data`.
- **Reasoning** ([`prototypes/reasoning/`](prototypes/reasoning/)) — grounded per-protein reasoning
  reports (incl. the 169208 disagreement case), a disagreement-detection / review-triage rule
  (AUROC 0.61 for flagging wrong calls), split-conformal prediction sets, and an OOD/novelty score.
  The supervised novelty score reaches **AUROC 0.786** (5-fold CV) vs. the current pipeline's 0.655
  baseline — a **+0.13** gain using only signals the pipeline already computes. Reproduce with
  `PYTHONPATH="$PWD" python run_example.py --assets ../example_data`.

See [`prototypes/SYNTHESIS_REPORT.md`](prototypes/SYNTHESIS_REPORT.md) for a cross-track review and
the recommendation of what to integrate into the product. These are prototypes, not yet wired into
the main pipeline or web UI.

---

## Documentation

Full documentation is written for **MkDocs Material** and is ready to publish on ReadTheDocs —
the build config (`mkdocs.yml`, `.readthedocs.yaml`, `docs/requirements.txt`) ships in this repo.

**Build & preview locally** (works today, no external service):

```bash
pip install mkdocs-material pymdown-extensions
mkdocs serve      # → http://127.0.0.1:8000
```

**Publish on ReadTheDocs** (one-time, requires repo-owner action — the site is *not* live until
you do this):

1. Sign in at [readthedocs.org](https://readthedocs.org) with the GitHub account that owns the repo.
2. **Import** the `dbcan4-advanced` repo — this installs the build webhook.
3. The first build runs automatically from the committed `.readthedocs.yaml`. The site then goes
   live at `https://dbcan4-advanced.readthedocs.io` (or the slug RTD assigns if that name is taken).

Key pages: [Design & architecture](docs/design_dbcan4_advanced.md) ·
[Benchmark report](docs/benchmark_report.md) · [Output contract](nf/OUTPUT_CONTRACT.md) ·
[Feature tools](nf/TOOLS.md).

---

## Citing the methods

dbCAN4-advanced builds on: **ESM-C** (EvolutionaryScale), **CLEAN** (Yu et al., *Science* 2023),
**Foldseek** (van Kempen et al., *Nat. Biotechnol.* 2023), **CAZyme3D** (Zheng, Yin et al.,
bioRxiv 2024), **SaProt** (*ICLR* 2024), **ProstT5** (Heinzinger et al., 2024), **ESMFold**,
**DeepTMHMM**, **Pfam/HMMER**, and **CAZyLingua** (Thurimella et al., *BMC Bioinformatics* 2025)
as the closest prior pLM-CAZyme tool. Full reference list with DOIs and provenance:
[`docs/design_dbcan4_advanced.md`](docs/design_dbcan4_advanced.md).

## License

MIT — see [LICENSE](LICENSE).
