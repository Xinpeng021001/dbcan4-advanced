# Vendored BioForge web layer

This `bioforge` package is **vendored** (copied) into dbCAN4-advanced so the whole
product — pipeline **and** the web UI/ingest layer — runs from a single
`git clone` with no second repository to fetch.

- Upstream: https://github.com/Xinpeng021001/biodb  (branch `feature/advanced-cazyme-integration`)
- Snapshot: commit `1f5d2cd` — "Protein-input ingest: genes from FASTA, GFF optional; topology track + structure viewer" (2026-07-10)

To refresh from upstream, re-copy `src/bioforge/`, `db/alembic/`, and `alembic.ini`
from the biodb checkout (see `scripts/vendor_bioforge.sh`).
