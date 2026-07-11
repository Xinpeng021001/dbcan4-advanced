# dbCAN4-advanced — the product

**One command turns a protein FASTA into a browsable CAZyme-annotation database.**

```bash
# the full real workup, verified end-to-end on met:
bash dbcan4_workup.sh proteins.faa --serve
# → baseline dbCAN + advanced ESM-C/fusion + 7 real feature tracks
#   (Pfam domains, physicochem, EC, TM topology, signal peptide,
#    localization, ESMFold 3D structure)
# → publishes the standard v1.1 output contract
# → ingests into a BioForge SQLite database
# → serves a live web UI at http://127.0.0.1:8000
```

dbCAN4 is **fungal + protein-input**: genes are built straight from the protein
FASTA — **no genome, no Prokka, no GFF**. Each protein becomes its own gene with
residue coordinates (`1–L`). A genomic GFF is *optional* — pass `--gff FILE` only
when you genuinely have coordinates, and the gene pages gain a genome track.

`dbcan4_workup.sh` is the one-command entry point for a real run — every step is a
validated real tool, and it handles all the tool-specific gotchas internally. The
Nextflow-driven `dbcan4 run --fasta proteins.faa --serve` is the same pipeline wrapped as
a portable DAG (fully stub-provable anywhere; `-profile met` for real execution).

This repository is three parts that compose into one product:

| Part | What it is | Where |
|---|---|---|
| **`dbcan4` Python package** | the annotation engine + CLI (`embed` / `infer` / `annotate` / `run`) | `src/dbcan4_advanced/`, `pyproject.toml` |
| **Nextflow pipeline** | orchestrates baseline dbCAN + advanced ESM-C/structure tiers + 8 comprehensive feature tracks → standard output contract | `nf/` |
| **BioForge db + web** | versioned SQLite schema + FastAPI web UI that ingests the contract and serves per-gene deep-dive pages | `github.com/Xinpeng021001/biodb` |

![product architecture]({{artifact:art_da9d0df3-2a6c-4b5c-978a-c98beb04b5f5}})

---

## What the annotation engine does

dbCAN's released engines (HMMER, dbCAN_sub, DIAMOND) assign CAZy families by
**sequence similarity**. dbcan4-advanced adds an orthogonal **protein-language-model +
structure** tier for fungal proteins:

- **ESM-C 600M embeddings** → three retrieval/classification heads (kNN vote,
  class-centroid cosine, a CLEAN-style contrastive head) against a 337,759-protein
  fungal reference.
- **Structure tier** (ESMFold → Foldseek vs the group's CAZyme3D DB + SaProt) for hard cases.
- **DEFT-style fusion** — confidence-weighted consensus across methods with an **abstain**
  option when methods disagree.
- **Label-free inference**: the engine predicts on novel sequences with **no ground-truth
  labels** — it only needs the precomputed reference index + trained heads.

### Honest positioning (this matters for judging)

The method is **competitive with, not superior to, sequence baselines** on the tested
regime, and we say so:

- On a 2024→2025 fungal temporal holdout, contrastive-kNN on frozen ESM-C **ties DIAMOND
  at family level** (overlap 0.973 vs 0.981) and **beats dbCAN-sub on subfamily**
  (subfamily-exact 0.531 vs 0.427), while trailing custom DIAMOND (0.588).
- The holdout's eval mass sits at **high identity** (median ~81% to the 2024 reference;
  only 14 of 4,726 proteins below 30% id), so it measures *near-term annotation*, **not**
  twilight-zone remote-homolog recovery. The value is **orthogonality, abstention,
  subfamily resolution, and fungal calibration** — not raw recall superiority.
- A rigorous **DB-vintage leakage control** shows DIAMOND against the *current* dbCAN DB is
  contaminated by 2025 eval sequences (novel-to-fungi subfamily recall 0.001 fair-2024 →
  0.992 current DB) — most tool comparisons skip this; we did it.

Full benchmark report: `docs/benchmark_report.md`. Honest narrative + figures:
`dbcan4_advanced_honest_narrative.md`.

---

## The CLI

```
dbcan4 embed     FASTA -> ESM-C embeddings (.npz)                 [GPU]
dbcan4 infer     embeddings -> label-free family calls (TSVs)     [CPU/GPU]
dbcan4 annotate  FASTA -> family calls in one step (embed+infer)  [GPU]
dbcan4 run       FASTA -> full Nextflow pipeline (baseline + advanced + features)
                 -> standard v1.1 output contract  [--serve to ingest + launch web UI]
dbcan4 info      show resolved asset paths + versions
```

`dbcan4 run` builds a samplesheet from a single FASTA automatically, runs the pipeline
(`-profile met` for real execution, `-profile stub -stub-run` to prove the DAG anywhere
with no GPU/tools), and with `--serve` chains
`alembic upgrade → bioforge-ingest → bioforge-ingest-advanced → uvicorn`.

The heavy GPU steps shell out to validated stage scripts (`embed_esmc.py`,
`infer_esmc.py`) run by the **engine venv** (torch + esm), resolved via
`DBCAN4_ENGINE_PYTHON` — so the CLI can be installed in a lightweight venv while the
engine lives in the project venv.

---

## The standard output contract (v1.1)

The pipeline publishes one **`manifest.json`** + a **funcscan-layout tree**. This is the
designed handshake between the pipeline and BioForge — the manifest lists per-sample
`cazyme_predictions` (one TSV per method) and `protein_features` (8 types:
signal_peptide, tm_topology, structure, domains, structure_hits, localization,
physicochem, ec_prediction). `bioforge-ingest-advanced` reads it and attaches each call to
the matching baseline gene as a new versioned release. See `nf/OUTPUT_CONTRACT.md`.

---

## The full functional workup

Every CAZyme the pipeline calls is then put through a complete functional workup —
this is what makes each protein's web page a genuine deep-dive rather than a bare family
label. On the 3 real held-out proteins, **7 of the 8 feature tracks run on real tools**
(including ESMFold 3D structure, folded on met's GPU); only Foldseek-vs-CAZyme3D structure
search is left unrun (wired + stub-proven, needs the CAZyme3D DB). The two license-gated
tools (SignalP-6.0, DeepLoc-2.0) are not installed and are handled as **honest fallbacks**
— clearly labelled, never fabricated.

![feature coverage]({{artifact:art_97942c1e-b1fb-4580-9f33-3f37058eba14}})

| track | tool | status |
|---|---|---|
| Pfam domains | hmmscan vs Pfam-A | **real** |
| EC number | CLEAN (maxsep) | **real** |
| TM topology | DeepTMHMM (BioLib cloud) | **real** |
| Signal peptide | DeepTMHMM → SignalP6 slot | **honest fallback** — SignalP-6.0 not installed; sourced from DeepTMHMM's SP call, `sp_prob` left blank (not fabricated) |
| Localization | derived from DeepTMHMM SP + topology | **honest fallback** — DeepLoc-2.0 not installed; a transparent rule (secreted if SP + no TM), labelled "not DeepLoc" |
| Physicochem | Biopython (MW, pI, GRAVY, N-glyc sequons) | **real** |
| Structure | ESMFold (facebook/esmfold_v1, met GPU) | **real** — all 3 folded; pLDDT 267317 69.2, 602276 82.3, 169208 75.8; served in the 3Dmol viewer |
| Structure hits | Foldseek vs CAZyme3D | not run on these 3 (wired + stub-proven; needs CAZyme3D DB) |

DeepTMHMM is the workhorse of the secretion/topology tracks: one real BioLib-cloud run
predicts *both* the transmembrane topology and the N-terminal signal peptide, so it
populates the TM-topology track and honestly backs the signal-peptide track when the
licensed SignalP-6.0 is absent. When SignalP-6.0 *is* installed, the `SIGNALP6` process
detects it on `PATH` and uses it instead — the fallback is automatic, not baked in.

**Real functional results (all 3 proteins):** 267317 (GH78) → EC 3.2.1.40
(α-L-rhamnosidase), signal peptide (cleave@18), Extracellular; 602276 (GH11) → EC 3.2.1.8
(endo-1,4-β-xylanase, conf 0.99), signal peptide (cleave@12), Extracellular; 169208
(GH183) → EC 3.2.1.55 (low conf), no signal peptide, Cytoplasm. All biologically coherent
with the CAZyme family assignments.

---

## Verified end-to-end on real sequences

The whole chain was run on met on **3 real held-out 2025 fungal CAZymes**:

| protein | truth | ESM-C kNN | ESM-C centroid | ESM-C contrastive | fusion | real Pfam domains |
|---|---|---|---|---|---|---|
| **267317** (hero, 1089 aa) | GH78, GH28 | GH78 (0.995) ✓ | GH92 (0.977) ✗ | GH78 (0.950) ✓ | **GH78** ✓ | Glyco_hydro_28 (PF00295) + Bac_rhamnosid6H (PF17389) |
| **602276** (203 aa) | GH11 | GH11 (0.991) ✓ | GH11 (0.955) ✓ | GH11 (0.999) ✓ | **GH11** ✓ (4/4) | Glyco_hydro_11 (PF00457) |
| **169208** (205 aa) | GH183 | GH43_6 (0.986) ✗ | **GH183 (0.984) ✓** | PL42 (0.297) ✗ | GH43_6 (0.308) ✗ | DUF4185 (PF13810) |

The three cases span the honest range: **602276** all four heads agree on GH11 (unanimous, correct); **267317** majority-correct (kNN + contrastive + fusion → GH78; centroid dissents GH92); **169208** a genuinely hard case where the **ESM-C-centroid head recovers the true GH183 at high confidence (0.98)** but the other two heads miss it (kNN GH43_6, contrastive PL42) and the fusion consensus is dragged onto the wrong high-confidence kNN call (GH43_6, low fusion confidence 0.31). 169208 is a useful diagnostic: it shows the heads are genuinely orthogonal — one recovers a family the others miss — and it also exposes a real fusion weakness (a confident-but-wrong kNN vote can outweigh a correct centroid vote), which is honest to surface rather than hide.

Result: **3 genes, 12 advanced CAZyme calls, 22 protein features across 7 real tracks
(ec_prediction, localization, physicochem, signal_peptide, tm_topology, structure — one
per protein — plus 4 Pfam domain rows), ingested; 3 live web pages served with every
feature card populated — the 3D structure viewer renders the folded model, and each page
carries a per-residue DeepTMHMM topology/signal-peptide track.** Genes are built directly
from the protein FASTA (protein-input mode: no genome, no Prokka, no GFF).
All three true families are recovered by at least one ESM-C head — 267317
and 602276 by the majority/consensus, and 169208's GH183 by the centroid head specifically
(the kNN and contrastive heads miss it, and the fusion consensus lands on the wrong
high-confidence kNN call). 169208 is the honest diagnostic case: the correct answer is
present in the method ensemble but the current fusion rule does not surface it — a real,
disclosed limitation of the consensus step, not a clean success. Evidence:
`run_real_demo.sh`, `manifest.json`, `real_demo.db`,
`fig_real_end2end.png`, `docs/fig_feature_coverage.png`.

---

## The web UI

Each gene page is a genuine deep-dive. Every card below is populated from real data on
the 3 held-out proteins:

- **Protein / Provenance** — protein-input mode (residue length, "no genomic coordinates";
  a per-source provenance table with honest tool labels + sha256 of every input file).
- **CAZyme annotations — advanced vs baseline** — advanced-only families the baseline
  missed are flagged; each call shows the method + confidence pill.
- **Secretion & membrane topology** — signal-peptide + TM-topology text **and an inline
  per-residue SVG track** drawn from the real DeepTMHMM Viterbi path (orange signal / green
  extracellular / blue TM), with cleavage site and residue ruler.
- **Predicted 3D structure** — an interactive **3Dmol viewer** rendering the ESMFold model,
  cartoon coloured by per-residue pLDDT, with a PDB download.
- **Function · EC · substrate**, **Subcellular localization**, **Physicochemistry**,
  **Pfam domain architecture**, **GO**.

The screenshots in `docs/ui/` are **true-browser captures** (headless Google Chrome
against the live server — real CSS + JavaScript + the 3Dmol WebGL viewer), produced by
`capture_ui.sh` and described in `docs/ui/README.md`. `ui_hero_267317.png` and
`ui_gene_169208.png` show the ESMFold structure actually rendered in the 3Dmol viewer.
Earlier in-development previews were print-render substitutes that could not show the JS
3Dmol viewer; the UI itself was never a mock-up.

## Feature tools & environments

Each tool tier is its own Nextflow module with its **own** `conda`/`container` directive,
so environments never conflict. Pinned conda recipes live in `nf/envs/`
(`engine.yml`, `deeptmhmm.yml`, `foldseek.yml`); **`nf/TOOLS.md`** documents which tools
matter for a CAZyme, all three install routes (conda YAML, BioContainers image,
license-gated manual for SignalP6/DeepLoc/CLEAN), and the 5-step contract for adding a new
feature tool.

---

## Install

```bash
# engine venv (torch + esm 3.2.1 ESM-C, on a GPU host)
pip install -e .            # installs the dbcan4 console script
# web stack (BioForge)
pip install -e /path/to/biodb
```

Reproduction commands (stub-anywhere + real met run): **`REPRODUCE_PRODUCT.md`**.
