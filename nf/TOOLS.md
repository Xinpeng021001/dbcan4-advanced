# Feature tools — what runs, why it matters, and how to install it

dbCAN4-advanced is built as **one Nextflow module per tool tier**
(`nf/modules/*.nf`). Every process carries its **own** `conda` / `container`
directive, so each tool resolves its own environment and **environments never
conflict** — DeepTMHMM's pinned client can sit next to the CUDA folding stack
next to CLEAN's Torch 1.11 without a single version clash. You pick how each
env is provided (conda prefix, conda YAML, or a container) per run via a
Nextflow profile; the pipeline code doesn't change.

```
nf/modules/
  baseline.nf        BASELINE_DBCAN                          run_dbcan V5 (HMMER+dbCAN-sub+DIAMOND)
  esmc.nf            ESMC_EMBED / _RETRIEVAL / _CONTRASTIVE   ESM-C family calls (advanced tier)
  features.nf        DEEPTMHMM / SIGNALP6 / FUSION            topology + signal peptide + consensus
  features_extra.nf  PFAM_DOMAINS / STRUCTURE_HITS /          domains, Foldseek, localization,
                     LOCALIZATION / PHYSICOCHEM / CLEAN_EC    physicochem, EC number
  structure.nf       ESMFOLD / FOLDSEEK_CAZYME3D / SAPROT     3D structure + structural homology
```

---

## Which feature tools matter for a CAZyme / enzyme (and their status here)

Ranked by how much they add to a fungal CAZyme annotation:

| Tool | Feature it adds | Why it matters for a CAZyme | Status |
|------|-----------------|------------------------------|--------|
| **run_dbCAN V5** | CAZy family (HMMER + dbCAN-sub + DIAMOND) | The core call — which GH/GT/PL/CE/AA/CBM family | **integrated** (conda) |
| **ESM-C heads** | family by embedding (kNN/centroid/contrastive→fusion) | Recovers novel/divergent fungal CAZymes the DBs miss | **integrated** (conda) |
| **DeepTMHMM** | per-residue TM topology + signal peptide | Is the enzyme **secreted**? Secreted CAZymes act on extracellular polysaccharides | **integrated** (pybiolib cloud) |
| **CLEAN** | EC number | Catalytic activity (e.g. EC 3.2.1.x glycosidase) independent of family | **integrated** (conda, research license) |
| **ESMFold** | 3D structure + per-residue pLDDT | Fold confirms family when sequence is ambiguous; feeds Foldseek | **integrated** (local GPU) |
| **Foldseek vs CAZyme3D** | structural homology to known CAZyme folds | Structure-level family evidence for the hardest cases | **wired** (needs CAZyme3D_id50 DB) |
| **Pfam / hmmscan** | domain architecture | Multi-domain CAZymes (e.g. GH + CBM) — domain map | **integrated** (conda) |
| **Biopython** | physicochemistry (MW, pI, GRAVY) | Basic biochemical context | **integrated** (conda) |
| **SignalP 6.0** | secretion signal (higher precision than DeepTMHMM's SP) | Sharper secreted/not call; eukaryote-tuned | **gated** — honest fallback used |
| **DeepLoc 2.0** | subcellular localization (10 compartments) | Where the enzyme acts (wall / vacuole / secreted) | **gated** — derived fallback used |

"Gated" = free for academics but **license-gated at DTU Health Tech**; not
pip/conda-installable and not redistributable. Until you register and drop the
tarball in, the pipeline uses an **honest fallback** (localization derived from
the real DeepTMHMM signal-peptide + topology, labelled as such — never
fabricated) and clearly marks the source on every row and in the web UI.

Worth considering later (not yet wired): **InterProScan** (families+GO+domains,
one big Java run), **DeepFRI** (structure→GO function), **eCAMI/CGC** substrate
context (dbCAN already emits CGC + dbCAN-PUL substrate).

---

## Three ways to provide a tool's environment

### 1. Conda YAML (recommended for the open tools)
Pinned recipes live in `nf/envs/`. Build once, point the pipeline at the prefix:

```bash
mamba env create -f nf/envs/engine.yml      # engine tier (GPU pLM + fold + baseline)
                                            # ESMFold rides this env: same torch/
                                            # transformers/esm stack as ESM-C, verified.
mamba env create -f nf/envs/deeptmhmm.yml    # pybiolib client
mamba env create -f nf/envs/foldseek.yml     # foldseek + hmmer

nextflow run nf/main.nf -profile met \
    --esmc_env    $(conda info --base)/envs/dbcan4-engine \
    --dbcan_env   $(conda info --base)/envs/dbcan4-engine \
    --esmfold_env $(conda info --base)/envs/dbcan4-engine \
    --deeptmhmm_env $(conda info --base)/envs/dbcan4-deeptmhmm \
    --foldseek_env  $(conda info --base)/envs/dbcan4-foldseek
```

Nextflow can also build these for you automatically: with `-profile conda`
(or `met`, which sets `conda.enabled=true`) point a process's `conda` directive
at the YAML and Nextflow creates the env on first run and caches it.

### 2. Container (best for reproducibility / HPC / no conda)
Each process can carry a `container` directive instead of `conda`. Turn it on
with `-profile docker` (or `singularity` / `apptainer` on HPC). Use community
**BioContainers** for the open tools — no image-building needed. Every
bioconda package has an auto-built image at
`quay.io/biocontainers/<tool>:<version>--<build>`; look up the exact
`<version>--<build>` tag for a package on <https://biocontainers.pro> or the
quay.io tag list (the build suffix is hashed, so don't guess it):

```groovy
// shape of the directive — resolve each real tag from the registry
process BASELINE_DBCAN    { container 'quay.io/biocontainers/dbcan:<ver>--<build>' }
process PFAM_DOMAINS      { container 'quay.io/biocontainers/hmmer:<ver>--<build>' }
process FOLDSEEK_CAZYME3D { container 'quay.io/biocontainers/foldseek:<ver>--<build>' }
process ESMFOLD           { container 'nvcr.io/nvidia/pytorch:<tag>' }  // + pip esm/transformers, or build once
```
```bash
nextflow run nf/main.nf -profile met,singularity     # profiles compose
```
Pin to a **digest** (`@sha256:…`) for a byte-identical rerun. CLEAN ships its own
`Dockerfile` (in the CLEAN repo) — build that image and set it as CLEAN_EC's
`container`.

### 3. License-gated manual install (SignalP6, DeepLoc, CLEAN)
These require accepting a license at the vendor, so they can't be automated:

- **SignalP 6.0** — register at <https://services.healthtech.dtu.dk/>, download
  the academic tarball, `pip install signalp-6-package/` into a dedicated env,
  and point `--signalp_env` at it. The SIGNALP6 process auto-detects the binary
  on PATH and **falls back honestly** to the DeepTMHMM signal-peptide call if
  absent (never fabricated; the fallback source is labelled on every row).
- **DeepLoc 2.0** — same DTU registration + tarball. Not yet a wired env param:
  LOCALIZATION currently **derives** the compartment from the real DeepTMHMM
  signal-peptide + topology (labelled "derived", not "DeepLoc"). To use real
  DeepLoc, add a `deeploc_env` param + a DeepLoc branch in the LOCALIZATION
  process (see "Adding a new feature tool" below) — the honest-fallback pattern
  is already there to copy.
- **CLEAN** — clone the CLEAN repo, follow its README
  (`conda create -n clean python==3.10.4 && pip install -r requirements.txt`,
  Torch 1.11), and set `--clean_env`. Non-exclusive research-use license.
  Already wired here via the CLEAN_EC process.

---

## Adding a *new* feature tool (the pattern)

1. Add a process to the right `nf/modules/*.nf` with a `script:` (real command +
   a converter that writes the v1.1 TSV) **and** a `stub:` (cp a canned TSV, so
   the DAG runs tool-free in CI).
2. Give it a `conda "${params.<tool>_env}"` (and/or a `container`) directive and
   add the env param to `nextflow.config`.
3. Write a converter in `nf/bin/` that emits the feature's TSV in the shape
   `parse_advanced.py` expects (see `OUTPUT_CONTRACT.md`).
4. Register the TSV → (feature_type, tool) in `nf/bin/write_manifest.py`.
5. Add a card to `gene_detail.html` if it should show in the web UI.

That's the whole contract — the DB/web layers pick it up from the manifest with
no further changes.
