# Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `dbcan4 info` shows `MISSING` for heads | data assets not present — copy from the reference host or set `DBCAN4_HEADS` / `DBCAN4_REF_EMB` |
| `from esm.models.esmc import ESMC` fails | `fair-esm` is installed and clashing — `pip uninstall fair-esm` then `pip install --force-reinstall --no-deps esm` |
| DeepTMHMM `FileNotFoundError: hash/…` | biolib stages by **basename** — the workup `cd`s into the FASTA dir and passes a bare filename; use `--no-deeptmhmm` offline |
| CLEAN errors on a `*` char / very long protein | ESM-1b has no `*` token and caps at 1022 aa — `run_clean.sh` strips stops + truncates; use `--no-clean` to skip |
| `dbcan4 run -profile met` (real Nextflow) not fully working | **use `dbcan4_workup.sh` for real runs** — the full-Nextflow *real* path is not proven end-to-end; only `-profile stub` is validated |
| Nextflow "requires Java 17+" | system Java is 11 — `source /array1/xinpeng/scratch/bin/nxf_env.sh` (portable Temurin 21) |
| Web UI blank / no 3D structure | pass `--structures-dir` on ingest (the workup does this); the 3Dmol viewer needs the served PDB |
| GPU out of memory during embed | lower the batch size, or set `CUDA_VISIBLE_DEVICES` to a free device; ESM-C 600M fits comfortably on a 24 GB card |

## Known limitations

- **`dbcan4 run -profile met` is not proven end-to-end.** Some `-profile met` conda directives
  point at env prefixes that are not all real conda envs. Use `dbcan4_workup.sh` for real runs
  and `dbcan4 run --stub` for the DAG proof.
- **DeepTMHMM needs outbound internet** (BioLib cloud). Use `--no-deeptmhmm` offline.
- **CLEAN is CPU-slow** (~2 min/protein) and truncates >1022 aa (ESM-1b limit). Use `--no-clean`
  to skip.
- **Fusion consensus can be led astray by a confident wrong call** — e.g. protein 169208, where
  the correct GH183 is recovered by the centroid head but fusion lands on the high-confidence kNN
  call (GH43_6). This is surfaced, not hidden.
- **Novelty / abstention detection is weak** (novelty AUROC ~0.63); the novel-family discovery
  pipeline is a triage step, not a certification pipeline.
- **Scale**: the full functional workup is verified at the 3-protein scale; ESMFold folding is
  the bottleneck at proteome scale.

More detail: [`REPRODUCE_PRODUCT.md`](https://github.com/Xinpeng021001/dbcan4-advanced/blob/main/REPRODUCE_PRODUCT.md)
and the [Feature tools](tools.md) page.
