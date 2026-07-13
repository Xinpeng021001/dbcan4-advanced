# CLI reference

```
dbcan4 embed     FASTA -> ESM-C embeddings (.npz)                 [GPU]
dbcan4 infer     embeddings -> label-free family calls (TSVs)     [CPU/GPU]
dbcan4 annotate  FASTA -> family calls in one step (embed+infer)  [GPU]
dbcan4 run       FASTA -> full Nextflow pipeline (+ --serve)
dbcan4 info      show resolved asset paths + versions
```

The heavy GPU steps shell out to the engine venv (torch + esm) via `DBCAN4_ENGINE_PYTHON`, so the
CLI itself can be installed in a lightweight venv.

## `dbcan4 info`

Prints the resolved asset paths and tool versions — run it first to confirm your install:

```bash
dbcan4 info
```

```text
  pipeline    : .../nf/main.nf
  ref index   : /array1/xinpeng/dbcan4-advanced/emb/ref2024
  heads.pt    : /array1/xinpeng/dbcan4-advanced/results/heads/heads.pt  (exists)
  nextflow    : /array1/xinpeng/scratch/bin/nextflow
```

## `dbcan4 embed`

FASTA → ESM-C embeddings. Mean-pooled last-layer `esmc_600m` representation (1152-dim).

```bash
dbcan4 embed proteins.faa --out emb.npz
```

## `dbcan4 infer`

Embeddings → label-free family calls. Runs three heads against the reference index:

- **kNN** — k=15 nearest reference embeddings; reports purity + margin.
- **centroid** — nearest family centroid.
- **contrastive** — the trained projection head (`heads.pt` + `proj_ref.npz`).

```bash
dbcan4 infer --emb emb.npz --out-knn knn.tsv --out-centroid cent.tsv --out-contrastive contr.tsv
```

## `dbcan4 annotate`

`embed` + `infer` in one step. **Label-free** — needs only the precomputed reference index and
trained heads, no ground-truth labels.

```bash
CUDA_VISIBLE_DEVICES=0 dbcan4 annotate proteins.faa --outdir calls_out
# → calls_out/ESM-C-{kNN,centroid,contrastive}.raw.tsv
```

## `dbcan4 run`

Wraps the Nextflow pipeline (`nf/main.nf`).

```bash
# stub: proves the whole DAG with no GPU/tools
dbcan4 run --fasta proteins.faa --sample s1 --outdir out --profile stub --stub

# stub + ingest + serve
dbcan4 run --fasta proteins.faa --sample s1 --outdir out \
    --profile stub --stub --serve --port 8000
```

`--serve` chains `alembic upgrade → bioforge-ingest → bioforge-ingest-advanced → uvicorn`.

!!! warning
    `-profile met` (the real Nextflow run) is **not proven end-to-end** — for real annotation use
    [`dbcan4_workup.sh`](usage.md).
