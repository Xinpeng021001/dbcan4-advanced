# dbCAN4-advanced

**Protein-language-model + structure-similarity CAZyme annotation for fungal proteins — beyond HMMER/DIAMOND.**

!!! tip "One command turns a protein FASTA into a browsable CAZyme-annotation database"
    ```bash
    bash dbcan4_workup.sh proteins.faa --serve      # → http://127.0.0.1:8000
    ```

Current dbCAN (`run_dbcan` / dbCAN3) assigns CAZy families by **sequence similarity**
(HMMER, dbCAN_sub, DIAMOND). This misses **remote-homolog CAZymes** — enzymes that share
fold, mechanism, and active-site geometry with known families but have drifted below the
sequence-identity detection threshold. dbCAN4-advanced adds an orthogonal
**protein-language-model (ESM-C) + structure (ESMFold/Foldseek/CAZyme3D)** tier, a consensus
**fusion** layer with an abstain option, and a full per-protein functional workup — then
ingests everything into a versioned database with a per-gene deep-dive web UI.

![Product architecture](architecture_product.png)

## Where to go next

<div class="grid cards" markdown>

-   :material-download: **[Installation](installation.md)**

    Requirements, the from-scratch install, and the data assets you need.

-   :material-rocket-launch: **[Quick start](quickstart.md)**

    Three paths — label-free calls, the stub DAG, and the full real workup.

-   :material-book-open-variant: **[Running the workup](usage.md)**

    Every step of `dbcan4_workup.sh`, all the options, and running on your own data.

-   :material-console: **[CLI reference](cli.md)**

    `embed` / `infer` / `annotate` / `run` / `info`.

-   :material-file-tree: **[Output contract](output-contract.md)**

    The standardized v1.1 manifest + funcscan tree the pipeline publishes.

-   :material-chart-box: **[Benchmarks](benchmark_report.md)**

    The 2024→2025 temporal holdout, honestly reported.

</div>

## What you get

dbCAN4-advanced is **three parts that compose into one product**:

| Part | What it is |
|---|---|
| **`dbcan4` Python package + CLI** | the annotation engine (`embed` / `infer` / `annotate` / `run`) |
| **Nextflow pipeline** | baseline dbCAN + advanced ESM-C/structure tiers + 8 feature tracks → a standard v1.1 output contract |
| **BioForge database + web UI** | versioned SQLite schema + FastAPI web app that ingests the contract and serves per-gene deep-dive pages |

**Per-protein functional workup** (8 feature tracks): Pfam domains (hmmscan), EC number (CLEAN),
TM topology + signal peptide (DeepTMHMM), 3D structure (ESMFold), subcellular localization,
physicochemistry (Biopython), plus structural-homology hits (Foldseek vs CAZyme3D).

## Honest positioning

On a 2024→2025 fungal temporal holdout the method is **competitive with, not superior to,
sequence baselines**, and these docs say so. Contrastive-kNN on frozen ESM-C **ties DIAMOND at
family level** (overlap 0.973 vs 0.981) and beats dbCAN-sub on subfamily. The eval mass sits at
high identity (median ~81% to the 2024 reference), so it measures near-term annotation, not
twilight-zone remote-homolog recovery. **The value is orthogonality, calibrated abstention,
subfamily resolution, and fungal calibration** — plus a rigorous DB-vintage leakage control most
tool comparisons skip. See the [benchmark report](benchmark_report.md).
