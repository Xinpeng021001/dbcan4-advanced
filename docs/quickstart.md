# Quick start

Three paths, fastest first. All use the shipped 3-protein example `examples/real3.faa`
(the held-out fungal CAZymes 267317 GH78, 602276 GH11, 169208 GH183).

=== "A · Label-free calls (GPU, ~1–2 min)"

    The real ESM-C engine only — no features, no web. Needs the reference index + trained heads.

    ```bash
    source /array1/xinpeng/scratch/biodb_venv/bin/activate      # reference host
    CUDA_VISIBLE_DEVICES=0 dbcan4 annotate examples/real3.faa --outdir calls_out
    ```

    Output — `calls_out/ESM-C-{kNN,centroid,contrastive}.raw.tsv`:

    | query | truth | kNN | centroid | contrastive |
    |---|---|---|---|---|
    | 267317 | GH78 | **GH78** (0.995) | GH92 (0.977) | **GH78** (0.950) |
    | 602276 | GH11 | **GH11** (0.990) | **GH11** (0.955) | **GH11** (0.999) |
    | 169208 | GH183 | GH43_6 (0.986) | **GH183** (0.984) | PL42 (0.297) |

    169208 is the honest hard case — only the centroid head recovers GH183.

=== "B · Stub DAG (no GPU/tools, ~1 min)"

    Proves the full Nextflow DAG + the v1.1 output contract on any machine.

    ```bash
    source /array1/xinpeng/scratch/bin/nxf_env.sh               # Nextflow + Java
    dbcan4 run --fasta examples/real3.faa --sample smoke \
        --outdir stub_out --profile stub --stub
    ```

    Output — `stub_out/cazyme_advanced/manifest.json` (contract v1.1: 6 prediction methods +
    8 feature types).

=== "C · The whole product (GPU, ~10–15 min)"

    Baseline + advanced + all 8 feature tracks, ingested and served.

    ```bash
    bash dbcan4_workup.sh examples/real3.faa --serve --gpu 0
    # → http://127.0.0.1:8000
    ```

    View the UI from your laptop (uvicorn binds `127.0.0.1` on the server):

    ```bash
    ssh -L 8000:127.0.0.1:8000 xinpeng@met.unl.edu     # then open http://localhost:8000
    ```

!!! success "Verified"
    Paths **A** (`dbcan4 annotate`) and **B** (`dbcan4 run --stub`) were re-verified from a clean
    working directory on the reference host on 2026-07-12 (the family calls and manifest shown
    above are from that run). Path **C** (the full `dbcan4_workup.sh --serve` product workup —
    DeepTMHMM/CLEAN/ESMFold/ingest/web UI) was verified end-to-end earlier, on 2026-07-10, on the
    same three held-out proteins. See [Running the workup](usage.md) for the full option set and
    how to run on your own data.
