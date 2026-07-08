# dbCAN4-advanced annotation pipeline (Nextflow DSL2)

A reproducible, one-command pipeline that runs the **advanced** fungal-protein
CAZyme annotation tiers and publishes the **standard output contract**
(`OUTPUT_CONTRACT.md`) that BioForge ingests as a new versioned release.

```
proteins.faa ─┬─ ESM-C embed ─┬─ kNN + centroid ─┐
              │               └─ contrastive     │
              ├─ ESMFold ─┬─ Foldseek/CAZyme3D ───┼─ FUSION ─┐
              │           └─ SaProt ──────────────┘          │
              ├─ SignalP6 ───────────(feature)               ├─ manifest.json
              └─ DeepTMHMM ──────────(feature)      structures┘
                                                              ↓
                             bioforge-ingest-advanced  →  SQL DB + web UI
```

## Layout

| path | role |
|------|------|
| `main.nf` | workflow (the DAG above) |
| `modules/esmc.nf` | ESM-C embed, kNN+centroid retrieval, contrastive head |
| `modules/structure.nf` | ESMFold, Foldseek vs CAZyme3D, SaProt |
| `modules/features.nf` | SignalP6, DeepTMHMM, fusion, manifest collation |
| `bin/normalize_predictions.py` | project any tool's raw TSV → standard §2.1 schema |
| `bin/fuse_predictions.py` | weighted-vote fusion (real logic, GPU-free) |
| `bin/write_manifest.py` | assemble `cazyme_advanced/manifest.json` |
| `nextflow.config` | params + profiles (`stub`, `met`, `conda`, `docker`, `singularity`) |
| `OUTPUT_CONTRACT.md` | the standardized output layout + TSV schemas + manifest spec |
| `assets/stub/` | tiny contract-shaped fixtures for `-stub-run` |

## Run

**Prove the DAG + contract anywhere (no GPU, no tools installed):**
```bash
nextflow run nf/main.nf -stub-run -profile stub \
    --input samplesheet.csv --outdir results
```
`samplesheet.csv` has columns `sample,faa`. Output lands in
`results/cazyme_advanced/` per the contract; feed it straight to
`bioforge-ingest-advanced results/cazyme_advanced/manifest.json`.

**Real execution on met (8× A5500, conda envs, CAZyme3D DB):**
```bash
nextflow run nf/main.nf -profile met \
    --input samplesheet.csv --outdir results \
    --esmc_env /array1/xinpeng/dbcan4-advanced/venv \
    --cazyme3d_db /array1/xinpeng/dbcan4-advanced/cazyme3d/db
```
met has no scheduler, so the `met` profile uses Nextflow's `local` executor
(jobs run directly, matching the project convention). Each process declares its
own `conda` env + a `stub:` block, so swapping stub→real is just dropping
`-stub-run`.

## Adding / swapping a method

1. Add a process to the relevant `modules/*.nf` (real command + `stub:` block).
2. Add its raw→standard column map to `bin/normalize_predictions.py::METHOD_MAP`.
3. Register the tool in `bioforge/methods.py` (family/kind/colour/display).

No change to the database or web layer is needed — they read the contract, not
the tool.
