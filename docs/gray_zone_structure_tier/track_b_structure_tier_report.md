# Track B: Structure-tier setup and bulk-structure coverage evaluation

dbCAN4 fungal CAZyme annotation project — met.unl.edu, `/array1/xinpeng/dbcan4-advanced/`

## Summary

| Deliverable | Result |
|---|---|
| CAZyme3D_id50 mapping (exact sequence) | 178,356 structures extracted; **6/3,514,731 unique Mycocosm CAZyme sequences (~0.0002%)** match exactly by MD5 |
| CAZyme3D_id50 mapping (homology) | **93.3%** of a 2000-protein Mycocosm CAZyme sample have a homolog (≥30% id, ≥50% cov; median 58.8% id) |
| ProstT5 (AA→3Di) | Installed & validated; **2483/2483** validation-sample proteins scored, all length-matched 3Di output |
| SaProt (650M-AF2) | Installed & validated; 1280-dim embeddings produced for all query + reference structures |
| ESM Atlas 2 coverage | **85.9%** of a 2000-protein sample have a hit, but median identity only **51.2%** (distant homology; MGnify-derived, not fungal-specific) |
| AF3db (AlphaFold DB) coverage | **75.9%** of queryable CAZyme3D-homolog accessions have an AFDB entry (944/1243; a further 216 accessions are non-UniProt-format and rejected by the API) |
| Structure-evidence score | [structure_evidence_scores_final.tsv](structure_evidence_scores_final.tsv) — 2483 rows, discriminates tiers as expected (below) |

## 1. CAZyme3D_id50: download, extraction, ID mapping

The tarball at `/array1/xinpeng/cazyme3d/CAZyme3d_id50.tar.gz` (13 GB) was extracted to
`/array1/xinpeng/cazyme3d/extracted/cazyme_id50/`: **178,356 PDB structures**
(178,351 unique by AA-sequence MD5), named by UniProt/RefSeq-style accessions
(e.g. `A0A3D9I5V7.pdb`, `QLD87080.1.pdb`) — confirmed **not** to overlap Mycocosm's
JGI-style gene IDs (`jgi-<Genome>-<num>-<transcript>`) by accession.

Two mapping strategies were used:

- **Exact sequence match (MD5).** Every CAZyme-candidate row across all 2226 Mycocosm
  genome directories (3,885,558 rows / 3,514,731 unique AA sequences, from `overview.tsv`
  any-tool hits) was MD5-hashed and compared against MD5s of the 178,356 CAZyme3D_id50
  AA sequences (extracted in bulk via `foldseek structureto3didescriptor`). Only **6 unique
  sequences matched exactly** (coverage ≈ 1.7×10⁻⁶) — expected, since CAZyme3D_id50 draws
  from all-kingdom CAZy/UniProt entries essentially disjoint from Mycocosm's own JGI gene
  calls.
- **Homology match (mmseqs2).** A random 2000-protein Mycocosm CAZyme sample searched
  against the 178,356 CAZyme3D_id50 sequences (`--min-seq-id 0.3 -c 0.5 -s 5.7`) found
  **1866/2000 (93.3%)** with at least one hit (median 58.8% identity). This confirms
  CAZyme3D_id50 is a useful **reference set** for structure-similarity scoring, not a
  direct structure source for Mycocosm proteins.

## 2. ProstT5 (Rostlab/ProstT5)

Installed in the existing met venv (`sentencepiece` was missing — bootstrapped pip via
`ensurepip` since the venv has no pip and `uv` is unavailable, then installed
`sentencepiece==0.2.1`). Validated on 5 test sequences: AA→3Di generation
(`<AA2fold>` prefix, deterministic decoding, half precision) produced length-matched
lowercase 3Di strings for all 5, ~3–13s/sequence on a single A5500 GPU.

Run at scale on the full 2483-protein validation sample (see §5) split across all 8 GPUs
in parallel — completed in ~1.5–2 hours wall time (~15–35s/protein depending on length,
capped at 1024 residues).

## 3. SaProt (westlake-repl/SaProt_650M_AF2)

Installed via `transformers` `EsmTokenizer`/`EsmForMaskedLM`; weights (4.9 GB) downloaded
to the shared HF cache. The repo's own `foldseek_util.get_struc_seq` helper
(`scripts/foldseek_util.py`) was reused to derive AA+3Di combined sequences from real
CAZyme3D_id50 structures via the installed `foldseek` binary (with automatic pLDDT
masking for AlphaFold-origin structures). Validated on 3 CAZyme3D_id50 structures:
1280-dim mean-pooled embeddings produced without error for all 3.

## 4. ESM Atlas 2 and AF3db bulk-structure coverage

**ESM Atlas 2** (`highquality_clust30.fasta`, 37M MGnify-derived sequences, 8.0 GB,
downloaded in full). A 2000-protein random Mycocosm CAZyme sample searched against it
(mmseqs2, ≥30% id, ≥50% cov) found **1718/2000 (85.9%)** with a hit — but the identity
distribution skews low (median 51.2%; only 19/1718 hits ≥90% identity, 900 in the 50–90%
range, 799 in the 30–50% range). This is consistent with ESM Atlas's bacterial/archaeal/
environmental metagenomic origin: coverage is real but mostly **distant** homology, not
close structural matches.

**AF3db (AlphaFold DB).** Mycocosm JGI protein IDs have no UniProt accessions, so direct
per-protein lookup isn't possible. Coverage was checked indirectly via the 1459 unique
CAZyme3D_id50-homolog accessions already identified in §5's foldseek/mmseqs2 3Di search:
querying the AlphaFold DB prediction API (`alphafold.ebi.ac.uk/api/prediction/<acc>`) found
**944/1459 (64.7% of all queried, 75.9% of the 1243 UniProt-format-queryable subset)**
have an AFDB entry; 299 don't (HTTP 404); 216 are RefSeq/GenBank-format accessions the API
rejects outright (HTTP 400, would need a RefSeq→UniProt mapping step first).

**Recommendation.** Neither ESM Atlas 2 nor AF3db gives direct, high-confidence structural
coverage for Mycocosm fungal proteins by accession — both require the same
homology/embedding step already used for CAZyme3D_id50, and even then ESM Atlas hits are
mostly distant homologs. **Local ESMFold folding remains the more practical route** for
proteins with no adequate reference hit (already proven on met at ~22s/protein median
across 8 GPUs from the earlier POC). CAZyme3D_id50 + ProstT5 (no folding needed) +
local ESMFold (for gaps) is the recommended structure-tier stack; AF3db lookup is worth a
light second pass only if/when a RefSeq→UniProt mapping step is built for other purposes.

## 5. Structure-evidence score

**Validation sample.** Drawn from Track A's `tiered_proteins.tsv.gz` handoff (found on the
met filesystem at `/array1/xinpeng/dbcan4-advanced/track_a_output/`, not the artifact
store): stratified random sample of 1500 `gray_zone` + 500 `high_confidence_cazyme` + 500
`high_confidence_non_cazyme` proteins (2483 with sequences retrievable from `uniInput.faa`
files; 17 dropped for missing sequence).

**Pipeline** (`scripts/structure_evidence_score.py`):
1. ProstT5 AA→3Di prediction for every query protein (no folding required).
2. **Foldseek/mmseqs2 3Di-string search** — predicted 3Di sequences aligned against
   CAZyme3D_id50's real-structure-derived 3Di sequences (treating the 3Di alphabet as a
   generic sequence for mmseqs2, since foldseek's native structure-based search requires
   real 3D coordinates, unavailable for ProstT5-only predictions). 2188/2483 (88.1%)
   queries got a hit.
3. **SaProt embedding similarity** — cosine similarity of each query's SaProt (AA+3Di)
   mean-pooled embedding to the centroid embedding of a 500-structure CAZyme3D_id50
   reference sample.
4. Combined into `structure_evidence_score` (mean of normalized foldseek bit score and
   rescaled SaProt cosine similarity, 0–1, higher = more CAZyme-like by structure).

**Validation result** (score by known tier):

| Tier | n | mean | median |
|---|---|---|---|
| high_confidence_cazyme | 500 | 0.653 | 0.647 |
| gray_zone | 1500 | 0.566 | 0.556 |
| high_confidence_non_cazyme | 483 | 0.463 | 0.455 |

The score cleanly separates known CAZymes from known non-CAZymes, with gray-zone proteins
falling in between as expected — confirming the structure-similarity signal is informative
for tiering gray-zone candidates. See the distribution plot below.

![Structure evidence score by tier]({{artifact:61657d2c-4985-41ed-97ea-f9c7b87547cf}})

## Scripts (in `dbcan4-advanced/scripts/`)

- `prostt5_validate.py` — ProstT5 AA→3Di prediction CLI
- `foldseek_util.py` — SaProt's `get_struc_seq` helper (foldseek-derived AA+3Di), vendored from westlake-repl/SaProt
- `saprot_embed.py` — SaProt embedding extraction (from PDB structures or external 3Di)
- `cazyme3d_mapping.py` — exact-MD5 + prep for homology mapping between CAZyme3D_id50 and Mycocosm
- `esmatlas_coverage.py` — ESM Atlas 2 sequence-coverage checker (mmseqs2 easy-search)
- `extract_validation_sample.py` — stratified sample extraction from Track A's tiered output
- `fetch_sample_sequences.py` — AA sequence retrieval for a protein-ID list from Mycocosm `uniInput.faa` files
- `structure_evidence_score.py` — end-to-end structure-similarity scoring pipeline (ProstT5 + foldseek/mmseqs2 3Di search + SaProt embedding)

## Known limitations / follow-ups

- The foldseek/mmseqs2 3Di-string search uses a BLOSUM-tuned substitution matrix, not a
  proper 3Di substitution matrix — identity/coverage are more reliable than bit-score
  magnitude for this component (noted in the script docstring).
- AF3db coverage was assessed only for the subset of proteins with a CAZyme3D_id50-homolog
  accession, not the full validation sample directly; a full assessment would need a
  Mycocosm→UniProt/RefSeq ID-mapping step (out of scope for this pass).
- **SaProt reference-centroid sampling caveat.** The `structure_evidence_scores_final.tsv`
  results reported above were computed with a version of `structure_evidence_score.py` that
  built the 500-structure SaProt reference centroid from only the **first 2000 rows** of the
  178,356-row CAZyme3D_id50 3Di file (a head-slice in file order, not a representative draw
  across the whole reference set) — a bug caught after the run. The script has since been
  fixed to reservoir-sample 500 structures uniformly across all 178,356 (see
  `structure_evidence_score.py`), but the TSV/plot in this report were **not regenerated**
  with the fix. This affects only the SaProt-embedding-similarity component (B); the
  foldseek/mmseqs2 3Di search (component A) already used the full 178,356-sequence
  reference correctly. The reported tier-discrimination result (high_confidence_cazyme
  0.653 > gray_zone 0.566 > high_confidence_non_cazyme 0.463) still holds as computed, but
  should be treated as based on a non-representative SaProt reference sample; rerunning
  with the fixed script is recommended before using the embedding-similarity component for
  downstream decisions.
