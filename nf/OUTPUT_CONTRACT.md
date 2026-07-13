# dbCAN4-advanced — Standard Output Contract (v1.1)

> **v1.1** adds five per-protein feature types produced by the comprehensive
> annotation stack: `domains` (§2.5, Pfam/hmmscan), `structure_hits` (§2.6, Foldseek),
> `localization` (§2.7), `physicochem` (§2.8), and `ec_prediction` (§2.9, CLEAN).
> All map onto the existing generic `ProteinFeature` table — **no schema migration**.

The advanced fungal-protein annotation pipeline (`nf/main.nf`) publishes a
**standardized, versioned output layout**. Downstream consumers — first of all
the BioForge ingester (`bioforge-ingest-advanced`) — read *only* this contract,
never a tool's raw output. That decoupling is the whole point: add or swap an
annotation method by changing one Nextflow process + one manifest entry, and the
database/web layer keeps working unchanged.

This contract deliberately **extends the nf-core/funcscan layout** that BioForge
already ingests, rather than replacing it. A funcscan run gives genes, baseline
dbCAN CAZyme calls, CGCs, InterPro/GO, sequences; this pipeline adds an
`cazyme_advanced/` tree with pLM/structure CAZyme calls and per-protein features.

---

## 1. Directory layout

```
<outdir>/
  cazyme_advanced/
    manifest.json                     # ← the contract entry point (see §3)
    predictions/
      <sample>/
        ESM-C-kNN.tsv                 # one normalized file per method (§2.1)
        ESM-C-centroid.tsv
        ESM-C-contrastive.tsv
        Foldseek-CAZyme3D.tsv
        SaProt.tsv
        fusion.tsv
    features/
      <sample>/
        signalp6.tsv                  # §2.2
        deeptmhmm.tsv                 # §2.3
        structures.tsv                # §2.4 (index)
        domains.tsv                   # §2.5 (Pfam/hmmscan)
        structure_hits.tsv            # §2.6 (Foldseek)
        localization.tsv              # §2.7 (DeepLoc/derived)
        physicochem.tsv               # §2.8 (Biopython)
        ec_prediction.tsv             # §2.9 (CLEAN)
        structures/
          <protein_id>.pdb            # served 3D structures
  funcscan/                           # baseline tree (bioforge-ingest)
    protein_annotation/
      interproscan/
        <sample>_interproscan_faa.tsv # §2.10 InterPro domains + GO terms
  pipeline_info/
    dbcan4_advanced_software_versions.yml
```

`<sample>` is the same `sample_key` BioForge already uses (from the funcscan
samplesheet / Prokka dir). `protein_id` in every TSV **must** match the gene
`ID`/`locus_tag` in that sample's Prokka GFF and `_cleaned.faa`, so the ingester
can join advanced calls onto the genes funcscan already loaded.

---

## 2. File schemas (TSV, tab-separated, `-` = null/abstain)

### 2.1 CAZyme prediction — `predictions/<sample>/<TOOL>.tsv`

Every method emits the **same normalized schema** (one row per protein per
call), regardless of how the underlying tool reports its result:

| column         | type   | meaning                                                       |
|----------------|--------|---------------------------------------------------------------|
| `protein_id`   | str    | joins to the gene key in the sample GFF / faa                 |
| `family`       | str    | predicted CAZy family/subfamily (e.g. `GH5`, `GH13_31`); `-` = abstain |
| `confidence`   | float  | calibrated score in [0,1]; `-` if the method gives none       |
| `ec`           | str    | EC number if the method assigns one, else `-`                 |
| `all_families` | str    | comma-separated alternatives / multi-domain families, else `-`|
| `extra`        | json   | method-specific detail (see below), one-line JSON object      |

`extra` carries the evidence a reviewer wants but the schema shouldn't hard-code:
- kNN: `{"k":15,"purity":1.0,"margin":0.42,"neighbors":["GH5","GH5",...]}`
- centroid/contrastive: `{"margin":0.13,"runner_up":"GH9"}`
- Foldseek-CAZyme3D: `{"target":"AF-...","tmscore":0.91,"lddt":0.88,"prob":0.99}`
- SaProt: `{"nn_id":"...","cosine":0.83}`
- fusion: `{"votes":{"ESM-C-kNN":"GH5",...},"agreement":3,"signals":["sequence","structure"]}`

The **`tool` key, its display name, method family (`baseline`/`advanced`),
method kind (`sequence-plm`/`structure`/`fusion`) and colour are NOT in the TSV**
— they live in the manifest (§3) and in `bioforge.methods.REGISTRY`, the single
source of truth shared by pipeline and web app.

### 2.2 Signal peptide — `features/<sample>/signalp6.tsv`  (SignalP 6.0)

| `protein_id` | `prediction` | `sp_prob` | `cs_position` | `extra` |
|---|---|---|---|---|
| str | `SP`/`NO_SP`/`LIPO`/`TAT` | float [0,1] | int (cleavage site aa) or `-` | json: per-class probs, cleavage motif |

### 2.3 Membrane topology — `features/<sample>/deeptmhmm.tsv`  (DeepTMHMM)

| `protein_id` | `prediction` | `n_tm` | `topology` | `extra` |
|---|---|---|---|---|
| str | `TM`/`SP+TM`/`SP`/`Globular` | int | topology string (e.g. `i12-34o...`) | json: region spans |

### 2.4 Structure index — `features/<sample>/structures.tsv`

| `protein_id` | `source` | `plddt` | `length` | `path` | `extra` |
|---|---|---|---|---|---|
| str | `ESMFold`/`AlphaFold`/`PDB` | float mean pLDDT | int aa | rel. path to `structures/<id>.pdb` | json: model, db accession |

### 2.5 Protein domains — `features/<sample>/domains.tsv`  (Pfam / hmmscan)

One row **per domain occurrence** (a protein with N domains has N rows), ordered N→C.

| `protein_id` | `acc` | `name` | `start` | `end` | `evalue` | `score` | `extra` |
|---|---|---|---|---|---|---|---|
| str | Pfam acc (e.g. `PF00295`) | domain name | int seq_start | int seq_end | float i-Evalue | float bitscore | json: hmm_coverage, clan, thresholding |

→ ingested as `feature_type="domain"`, `tool="Pfam/hmmscan"`, `label=<name>`, `score=<bitscore>`, `start/end`, `attributes={acc,evalue,hmm_coverage,clan}`.

### 2.6 Structural-homology hits — `features/<sample>/structure_hits.tsv`  (Foldseek)

Top structural homologs of the predicted structure against a reference set (e.g. CAZyme3D).
One row per hit (keep the top-K, ranked by bitscore).

| `protein_id` | `target` | `target_family` | `tmscore` | `prob` | `lddt` | `evalue` | `extra` |
|---|---|---|---|---|---|---|---|
| str | reference id | CAZy family of target (`-` if unknown) | float [0,1] | float [0,1] | float [0,1] | float | json: bits, fident, alnlen, reference_db |

→ ingested as `feature_type="structure_hit"`, `tool="Foldseek-CAZyme3D"`, `label=<target_family>`, `score=<tmscore>`, `attributes={target,prob,lddt,evalue,rank}`.

### 2.7 Subcellular localization — `features/<sample>/localization.tsv`  (DeepLoc / derived)

| `protein_id` | `localization` | `confidence` | `method` | `extra` |
|---|---|---|---|---|
| str | e.g. `Extracellular` | float or qualitative | `DeepLoc-2.0`/`derived-from-SP+GO-CC` | json: signals, go_cc_terms |

→ `feature_type="localization"`, `tool=<method>`, `label=<localization>`, `score=<confidence if numeric>`, `attributes={method,basis,signals}`.

### 2.8 Physicochemistry — `features/<sample>/physicochem.tsv`  (Biopython)

One row per protein (summary features).

| `protein_id` | `mw_da` | `pi` | `instability` | `gravy` | `aromaticity` | `extra` |
|---|---|---|---|---|---|---|
| str | float | float | float | float | float | json: aa_composition, n_glyc_sequons |

→ `feature_type="physicochem"`, `tool="Biopython"`, `score=<mw_da>`, `attributes={pi,instability,gravy,aromaticity,aa_composition,n_glyc_sequons}`.

### 2.9 EC-number prediction (sequence-based) — `features/<sample>/ec_prediction.tsv`  (CLEAN)

Independent, sequence-based EC prediction (orthogonal to family-inherited EC).
One row per predicted EC (keep the top ranks).

| `protein_id` | `ec_number` | `confidence` | `rank` | `tool` | `extra` |
|---|---|---|---|---|---|
| str | e.g. `3.2.1.40` | float [0,1] | int | `CLEAN` | json: confidence_type, model, agreement_with_family |

→ `feature_type="ec_prediction"`, `tool="CLEAN"`, `label=<ec_number>`, `score=<confidence>`, `attributes={rank,confidence_type,model}`.

> **Note on `ec_number`**: the CAZyme call's own EC (family-inherited) continues to live on
> `CazymeAnnotation.ec_number` (§2.1). §2.9 is the *independent predictor's* EC and is stored as a
> `ProteinFeature` so both lines of evidence coexist and can be compared in the UI.

### 2.10 InterPro domains + GO terms — `funcscan/protein_annotation/interproscan/<sample>_interproscan_faa.tsv`

Unlike §2.1–2.9 (which live under `cazyme_advanced/` and are read from the manifest),
InterPro/GO ride on the **baseline funcscan tree** and are loaded by `bioforge-ingest`
(not `bioforge-ingest-advanced`), populating the gene page's **Gene Ontology** card and
**InterPro domains** table. The file is the standard **headerless, positional InterProScan v5
TSV** (`bioforge.ingest.parse_interpro`): cols `1 protein  2 md5  3 length  4 analysis
5 signature_acc  6 signature_desc  7 start  8 stop  9 score/evalue  10 status  11 date
12 interpro_acc  13 interpro_desc  14 GO(|-sep)  15 pathways`.

Produced by the `INTERPROSCAN` process (`nf/modules/interproscan.nf`) / step 4b of
`dbcan4_workup.sh`, which uses **real InterProScan** when `params.interproscan_sh` /
`INTERPROSCAN_SH` points at an install, otherwise derives Analysis=Pfam signatures + GO
from the §2.5 Pfam domains via the bundled `pfam2go` map (offline). `interpro_acc`/`desc`
are filled only when a real InterProScan run or an optional `--pfam2interpro` map supplies them.

---

## 3. `manifest.json` — the contract descriptor

The ingester reads this **first** and follows it; it never globs blindly. It
declares which files exist, which registry `tool` each maps to, and (optionally)
column overrides so even a *non-normalized* legacy/raw TSV can be ingested
without new code.

```jsonc
{
  "contract_version": "1.0",
  "pipeline": "dbcan4-advanced",
  "pipeline_version": "0.1.0",
  "created": "2026-07-08T00:00:00Z",
  "release_label": "advanced-2026-07-08",       // becomes the BioForge Release label
  "release_notes": "ESM-C + structure advanced CAZyme calls",
  "tool_versions": { "esm": "3.2.1", "foldseek": "9.427df8a", "signalp": "6.0", ... },
  "samples": [
    {
      "sample_key": "demo_fungal",
      "cazyme_predictions": [
        { "tool": "ESM-C-kNN", "path": "predictions/demo_fungal/ESM-C-kNN.tsv" },
        { "tool": "fusion",    "path": "predictions/demo_fungal/fusion.tsv" },
        // Optional raw-TSV mode — map arbitrary columns to the standard fields:
        { "tool": "ESM-C-centroid",
          "path": "raw/esmc_retrieval_pred.tsv",
          "id_col": "query_id", "family_col": "cent_pred", "confidence_col": "cent_conf" }
      ],
      "protein_features": [
        { "feature_type": "signal_peptide", "tool": "SignalP6",         "path": "features/demo_fungal/signalp6.tsv" },
        { "feature_type": "tm_topology",    "tool": "DeepTMHMM",         "path": "features/demo_fungal/deeptmhmm.tsv" },
        { "feature_type": "structure",      "tool": "ESMFold",           "path": "features/demo_fungal/structures.tsv" },
        { "feature_type": "domain",         "tool": "Pfam/hmmscan",      "path": "features/demo_fungal/domains.tsv" },
        { "feature_type": "structure_hit",  "tool": "Foldseek-CAZyme3D", "path": "features/demo_fungal/structure_hits.tsv" },
        { "feature_type": "localization",   "tool": "DeepLoc",           "path": "features/demo_fungal/localization.tsv" },
        { "feature_type": "physicochem",    "tool": "Biopython",         "path": "features/demo_fungal/physicochem.tsv" },
        { "feature_type": "ec_prediction",  "tool": "CLEAN",             "path": "features/demo_fungal/ec_prediction.tsv" }
      ]
    }
  ]
}
```

A `tool` value **must** exist in `bioforge.methods.REGISTRY` (else the ingester
errors loudly, exactly like the baseline dbCAN parser does on an unknown column).
When `family_col`/`confidence_col` are omitted the ingester assumes the
normalized §2.1 schema (`family`,`confidence`).

---

## 4. Provenance & versioning

- Each `manifest.json` → one BioForge **Release** (`release_label`), loaded
  additively alongside the baseline release. Advanced-vs-baseline is then a
  *query across releases*, never a mutation.
- The ingester records per-file sha256 + `tool_versions` as `Provenance` rows,
  so re-ingesting an identical manifest is an idempotent no-op (same mechanism
  as the baseline loader).
- `contract_version` is bumped on any breaking schema change; the ingester
  checks it.
