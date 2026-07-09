# Annotation tool skills

One reusable skill per annotation tool in the comprehensive CAZyme-annotation
stack. Each skill is a self-contained `SKILL.md` (usage + gotchas) plus a
`kernel.py` sidecar (helper functions auto-loaded into the analysis kernel).
All were authored and verified against real runs on met.unl.edu.

| Skill | Tool | Kind | What it does |
|-------|------|------|--------------|
| `esmfold-fold` | ESMFold (ESM-C venv / NIM) | real | Predict 3D structure from sequence; pLDDT B-factors |
| `pfam-scan` | HMMER hmmscan vs Pfam-A | real | All Pfam domains with `--cut_ga`; domain architecture |
| `deeptmhmm` | DeepTMHMM (DTU, BioLib) | real | Signal peptide + transmembrane topology |
| `protein-function` | InterPro + dbCAN + Biopython | real | GO/EC/substrate + physicochemistry + N-glyc sequons |
| `clean-ec` | CLEAN (Yu et al. 2023) | real | Independent sequence→EC prediction with confidence |
| `signalp6` | SignalP-6.0 (DTU) | scaffold | Signal-peptide prediction (license-gated; install stub) |
| `deeploc` | DeepLoc-2.0 (DTU) | scaffold | Subcellular localization (license-gated; install stub) |

**real** = ran for real on met with real outputs in this project's example.
**scaffold** = license-gated (DTU academic download); skill documents the
install path and provides the honest derived fallback.

See `examples/267317_comprehensive/` for a full real annotation record
(multi-domain GH28/GH78 α-L-rhamnosidase) produced by running this stack.
