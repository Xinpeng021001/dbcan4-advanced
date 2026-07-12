# Feature tools

The comprehensive per-protein workup runs eight feature tracks, each a **reusable annotation-tool
skill** (`skills/`) with a documented install path and honest fallback. Five ran for real on the
reference host; two are license-gated scaffolds (SignalP-6.0, DeepLoc-2.0).

| Skill | Tool | Kind | What it does |
|---|---|---|---|
| `esmfold-fold` | ESMFold | real | 3D structure from sequence; pLDDT B-factors |
| `pfam-scan` | HMMER hmmscan vs Pfam-A | real | All Pfam domains with `--cut_ga`; domain architecture |
| `deeptmhmm` | DeepTMHMM (BioLib) | real | Signal peptide + transmembrane topology |
| `protein-function` | InterPro + dbCAN + Biopython | real | GO/EC/substrate + physicochemistry + N-glyc sequons |
| `clean-ec` | CLEAN (Yu et al. 2023) | real | Independent sequence→EC prediction |
| `signalp6` | SignalP-6.0 (DTU) | scaffold | Signal-peptide prediction (license-gated) |
| `deeploc` | DeepLoc-2.0 (DTU) | scaffold | Subcellular localization (license-gated) |

The canonical `nf/TOOLS.md` — what each tier runs, why it matters, and how to install it — is
reproduced below.

--8<-- "nf/TOOLS.md"
