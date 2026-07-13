# Reproduce the product end-to-end

Two paths: **(A) stub — provable anywhere** (no GPU, no tools, no data; proves the whole
DAG + output contract + ingest + web in ~1 min), and **(B) real — on met** (genuine
run_dbcan + ESM-C + features on real sequences).

Every command below was run exactly as shown. `$REPO = /array1/xinpeng/dbcan4-advanced`.

---

## 0. Toolchain (once)

Nextflow needs Java 17+. met ships Java 11, so a portable Temurin JDK 21 + Nextflow live
in scratch:

```bash
source /array1/xinpeng/scratch/bin/nxf_env.sh   # sets JAVA_HOME (Temurin 21), PATH (nextflow 26.04.4),
                                                 # NXF_HOME, CAPSULE_LOG=none
nextflow -version                                # 26.04.4 build 12445
```

Install the two packages:

```bash
# engine venv (has torch 2.12.1+cu130, esm 3.2.1 ESM-C) — the project venv
source $REPO/venv/bin/activate

# a venv that has the web stack (FastAPI/uvicorn/sqlalchemy/alembic) + the dbcan4 CLI
source /array1/xinpeng/scratch/biodb_venv/bin/activate   # or any venv with `pip install -e .`
pip install -e $REPO/repo_clone      # -> `dbcan4` CLI + the vendored bioforge web layer
                                     #    (bioforge-ingest, bioforge-ingest-advanced, web app)

dbcan4 info                          # resolves pipeline dir, ref index, heads.pt
```

---

## A. Stub — provable anywhere (no GPU / tools / data)

```bash
source /array1/xinpeng/scratch/bin/nxf_env.sh
source /array1/xinpeng/scratch/biodb_venv/bin/activate   # or any venv with `pip install -e .`
cd /array1/xinpeng/scratch/nf_stub

# one command: build samplesheet from a FASTA, run all 16 processes, publish v1.1 manifest
dbcan4 run --fasta demo_fungal.faa --sample demo --outdir cli_results --profile stub --stub

# result: cli_results/cazyme_advanced/manifest.json  (contract 1.1, 6 preds + 8 feats)
#         cli_results/funcscan/{annotation,cazyme,protein_annotation}/...

# ingest + verify the web serves (stub uses self-consistent demo_p01/p02/p03 ids)
DB=$PWD/cli_serve.db; rm -f $DB; export DATABASE_URL=sqlite:///$DB
( cd $REPO/repo_clone && alembic upgrade head )   # alembic.ini + db/alembic are vendored at the checkout root
bioforge-ingest          cli_results/funcscan
bioforge-ingest-advanced cli_results/cazyme_advanced/manifest.json
uvicorn bioforge.api.main:app --host 127.0.0.1 --port 8000   # browse http://127.0.0.1:8000
```

Stub run → **3 genes, 7 baseline + 16 advanced cazymes, 14 features across all 8 types**.

---

## B. Real — on met (genuine annotation of real sequences)

The entire real workup is one command — **`dbcan4_workup.sh`** at the repo root. Give it a
protein FASTA; it runs baseline dbCAN + advanced ESM-C/fusion + all comprehensive feature
tracks, writes the v1.1 manifest, and (with `--serve`) ingests + launches the web UI:

```bash
source /array1/xinpeng/scratch/bin/nxf_env.sh      # (only needed if you also use `dbcan4 run`)
# ONE command: FASTA -> full CAZyme workup -> browsable web UI
bash $REPO/repo_clone/dbcan4_workup.sh my_proteins.faa --serve
#   → http://127.0.0.1:8000
```

That single call runs, in order: **(1)** baseline dbCAN (run_dbcan V5) → **(2)** ESM-C embed
+ label-free infer (kNN/centroid/contrastive) → **(3)** fusion → **(4)** Pfam domains
(hmmscan) → **(5)** physicochem (Biopython) → **(6)** DeepTMHMM topology + signal peptide +
derived localization → **(7)** CLEAN EC → **(7b)** ESMFold 3D structure (local GPU) →
**(8)** v1.1 manifest, then ingest + serve. It handles all the real-tool gotchas internally
(biolib bare-filename staging, CLEAN stop-codon strip + truncation, flat manifest staging).

dbCAN4 is **fungal + protein-input**: genes come straight from the FASTA — **no genome, no
Prokka, no GFF**. Pass `--gff FILE` only if you genuinely have genomic coordinates.

```
Options:  --outdir DIR   --sample NAME   --gff FILE   --serve   --port N   --gpu N
          --no-deeptmhmm    (skip the BioLib-cloud step, e.g. offline)
          --no-clean        (skip CLEAN EC)
          --no-structure    (skip ESMFold folding, e.g. no GPU)
```

Verified end-to-end on the 3 held-out proteins (`real3.faa`): the one command produced
**3 genes, 12 advanced calls, 22 protein features across 7 real tracks (domains,
ec_prediction, localization, physicochem, signal_peptide, tm_topology, structure),
0 unmatched**, and served `/genes/1` (267317), `/genes/2` (602276), `/genes/3` (169208)
with **every feature card populated (no "not analysed")** — the 3Dmol viewer renders the
ESMFold model and each page shows the per-residue DeepTMHMM topology track. The funcscan
tree carries **no `annotation/prokka/` directory** (protein-input mode).

**True-browser screenshots:** `capture_ui.sh <base_url> <out_dir> [gene_ids...]` drives
headless Google Chrome (real CSS + JS + WebGL) against the live server, so the 3Dmol viewer
and topology tracks are captured as a user actually sees them — not a print-render
substitute.

For just label-free family calls on a novel sequence (no features, no web), the CLI is a
one-liner too:

```bash
dbcan4 annotate my_proteins.faa --outdir calls_out
#   267317 -> GH78 (kNN 0.9953, contrastive 0.9495)      [truth GH78,GH28]
#   602276 -> GH11 (kNN 0.9905, contrastive 0.9986)      [truth GH11]
#   169208 -> centroid GH183 ✓, kNN GH43_6, fusion GH43_6 [truth GH183 — hard case]
```

---

## Notes / gotchas

- **Engine python**: `dbcan4 embed/infer/annotate` shell out to `$REPO/venv/bin/python`
  (torch + esm), resolved as `Assets.engine_python` (override `DBCAN4_ENGINE_PYTHON`).
  The CLI itself can live in any venv.
- **write_manifest.py** enumerates prediction/feature TSVs **flat in `--stage-dir`** (as
  Nextflow stages them in the task CWD) and records contract-relative
  `predictions/<sample>/` + `features/<sample>/` paths. For a manual run, symlink the
  nested TSVs flat into a stage dir first (see `run_real_demo.sh` step 5).
- **DATABASE_URL** drives everything (alembic, both ingesters, the web app) — set it once,
  all four agree (`bioforge.config.database_url()`).
- `alembic upgrade head` applies 3 revisions: initial → v2 (sequences/CGC/ARG/GO) → v3
  (advanced methods + protein features).
- **DeepTMHMM (biolib)**: the biolib venv is `/array1/xinpeng/scratch/venv_biolib`
  (biolib 1.4.262). `biolib run ... --fasta X` stages the arg by **basename** on the
  cloud, so `cd` into the FASTA's directory and pass a **bare filename** (an absolute path
  fails cloud-side). One run yields both `deeptmhmm.tsv` (tm_topology) and, via
  `deeptmhmm_to_tsv.py --out-sp`, the signal-peptide TSV.
- **CLEAN EC**: run via `run_clean.sh` (in `nf/bin/`), which encodes three gotchas —
  (1) strip the trailing stop codon `*` from every sequence (ESM-1b alphabet has no `*`
  token and errors the whole batch), (2) pass the CLEAN data name **bare** (`-d real3`;
  `CLEAN_infer_fasta.py` prepends `inputs/` itself), (3) `--truncation_seq_length 1022`
  for the ESM-1b 1022-token limit (267317 is 1089 aa).
- **SignalP-6.0 / DeepLoc-2.0 are license-gated and not installed.** The pipeline handles
  both as honest fallbacks: `SIGNALP6` uses the real binary if it's on `PATH`, else copies
  the DeepTMHMM-derived signal peptide; localization is a transparent derived rule
  (secreted if SP + no TM), labelled "not DeepLoc". Neither fabricates a probability.
- **ESMFold structure**: `scripts/fold_esmfold.py` folds each protein on the local GPU via
  `facebook/esmfold_v1` (transformers) from the **engine venv** (torch + esm — same stack
  as ESM-C, so no separate env). `structures_to_tsv.py` then emits `structures.tsv` and
  **rescales the per-residue pLDDT B-factors ×100 (0–1 → 0–100)** so the 3Dmol viewer's
  colour ramp (min 50, max 90) is correct. Weights (~16 GB) cache under `hf_cache/`; the
  hero (1089 aa) folds in ~4 min, the two ~200-aa proteins in seconds.
- **Tool environments**: each Nextflow module carries its own `conda`/`container` directive
  so envs never conflict. Pinned conda recipes are in `nf/envs/`; `nf/TOOLS.md` covers all
  three install routes (conda YAML, BioContainers image, license-gated manual) and how to
  add a new feature tool.
