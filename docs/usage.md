# Running the workup

`dbcan4_workup.sh` is the **verified real path** — it runs each proven tool directly (handling
every tool-specific gotcha internally) and is the recommended way to annotate real sequences.

```bash
bash dbcan4_workup.sh <proteins.faa> [options]
```

## What it does, step by step

| # | Step | Tool | Notes |
|---|---|---|---|
| 1 | Baseline CAZyme calls | `run_dbcan` V5 (HMMER · dbCAN_sub · DIAMOND) | the sequence-similarity tier |
| 2 | Advanced family calls | ESM-C embed + label-free infer (kNN / centroid / contrastive) | GPU |
| 3 | Fusion consensus | rule-based vote across tiers, with abstain | |
| 4 | Pfam domains | `hmmscan` vs Pfam-A (`--cut_ga`) | domain architecture |
| 5 | Physicochemistry | Biopython | MW, pI, GRAVY, instability, N-glyc sequons |
| 6 | TM topology + signal peptide | DeepTMHMM (BioLib cloud) | derives localization |
| 7 | EC number | CLEAN | independent sequence→EC |
| 7b | 3D structure | ESMFold | local GPU; pLDDT B-factors |
| 8 | Manifest + (optional) serve | `write_manifest.py` → BioForge ingest → uvicorn | v1.1 contract |

## Options

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

!!! note "dbCAN4 is fungal + protein-input"
    Genes are built straight from the protein FASTA — **no genome, no Prokka, no GFF**. Each
    protein becomes its own gene with residue coordinates `1–L`. Pass `--gff` only when you
    genuinely have genomic coordinates.

## Timing

For the 3-protein example, DeepTMHMM (cloud), CLEAN (~2 min/protein CPU), and ESMFold (the hero
1089-aa protein ~4 min on GPU) dominate — expect ~10–15 min total.

**Fast smoke variant** — skips the slow/network steps, proving the whole chain except those 3 tracks:

```bash
bash dbcan4_workup.sh examples/real3.faa --serve --gpu 0 \
    --no-deeptmhmm --no-clean --no-structure
```

## Run on your own data

Point it at any fungal protein FASTA:

```bash
bash dbcan4_workup.sh /path/to/my_proteome.faa --sample my_fungus --serve --gpu 0
```

!!! warning "Scale"
    The full functional workup is verified at the 3-protein scale. At proteome scale (thousands),
    ESMFold folding is the bottleneck — consider `--no-structure` for a first pass, or fold a
    shortlist afterward.

## The two other entry points

- **`dbcan4 annotate` / `dbcan4 run --stub`** — see the [Quick start](quickstart.md) and the
  [CLI reference](cli.md).
- **`dbcan4 run -profile met`** (full Nextflow *real* run) — **not proven end-to-end**; use
  `dbcan4_workup.sh` for real runs and `dbcan4 run --stub` for the DAG proof. See
  [Troubleshooting](troubleshooting.md).
