# BioForge web UI — true-browser captures

These are **real headless-Google-Chrome screenshots** of the live BioForge web
server (`uvicorn bioforge.api.main:app`), taken by [`capture_ui.sh`](../../capture_ui.sh)
on the 3 held-out fungal proteins. They render the actual CSS + JavaScript + the
3Dmol WebGL structure viewer — i.e. what a user sees in a browser, not a static
figure.

| file | page | shows |
|---|---|---|
| `ui_dashboard.png` | `/` | dataset overview: 3 protein-mode genes, `dbcan4-protein` tool, 12 CAZyme calls, 6 families |
| `ui_browse.png` | `/browse` | gene table with protein coordinates (`protein:1–1,089`) and advanced-only flags |
| `ui_hero_267317.png` | `/genes/1` | hero (GH78, 1089 aa): **rendered ESMFold 3D structure** (cartoon by pLDDT), per-residue DeepTMHMM topology/SP track, EC, localization, physicochem, honest provenance |
| `ui_gene_169208.png` | `/genes/3` | hard case (GH183, 205 aa): **rendered ESMFold structure**, all-green globular topology (no SP) |

Notes:
- The 3Dmol viewer needs WebGL; headless Chrome has no GPU, so `capture_ui.sh`
  runs ANGLE's SwiftShader software-WebGL backend with a long virtual-time budget.
  Even so the software render is timing-sensitive — occasionally a small protein's
  canvas doesn't finish before the screenshot; that is a capture artifact, not a
  serving bug (every PDB serves `HTTP 200 chemical/x-pdb` and every other card
  populates). In a real browser all three load normally.
- Earlier in-development previews were WeasyPrint PDF→PNG renders, which flatten
  CSS and drop all JavaScript (no 3Dmol viewer, no SVG tracks). Those are **not**
  in this repo; only these true-browser captures are.
