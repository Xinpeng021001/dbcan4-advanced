# Gray-Zone Adjudication + Structure Tier — Synthesis Report

dbCAN4 fungal CAZyme annotation project · met.unl.edu · 2026-07-09

This report converges two parallel tracks run against `/array1/xinpeng/all_genome`
(2,226 Mycocosm fungal genomes, ~40 GB, dbCAN-annotated) and `/array1/xinpeng/cazyme3d`
(CAZyme3D_id50): **Track A** built a multi-evidence confidence-tiered CAZyme/non-CAZyme
dataset from sequence evidence alone; **Track B** stood up the structure tier (CAZyme3D_id50,
ProstT5, SaProt) and used it to validate/adjudicate a sample of Track A's gray zone.

## 1. The gray-zone problem, quantified

Every genome directory under `all_genome` ships `overview.tsv` (only proteins with ≥1 tool
hit), plus the raw per-tool tables (`diamond.out`, `dbCAN_hmm_results.tsv`,
`dbCANsub_hmm_results.tsv`) and `non_CAZyme.faa` (dbCAN's current #ofTools<2 rejects).
Streaming all 2,226 genomes (0 failures, ~750s) and recomputing `n_tools` from the raw
per-tool tables directly, rather than trusting the pre-filtered `overview.tsv`/`non_CAZyme.faa`
split, gives 28,192,456 total protein rows:

| Tier | Definition | Count | % |
|---|---|---:|---:|
| `high_confidence_cazyme` | ≥2 tools agree on a CAZy family | 946,270 | 3.36% |
| **`gray_zone`** | **exactly 1 tool hits a CAZy family** | **2,844,297** | **10.09%** |
| `high_confidence_non_cazyme` | 0 tools hit anything | 24,401,889 | 86.55% |

The gray zone answers the question this session opened with: *"if dbCAN.hmm finds more
CAZymes than CAZy labels, are those real CAZymes we're missing, or false positives?"* — it is
neither uniformly one nor the other. A single-tool hit is genuine homology evidence (each raw
per-tool file is already filtered to that tool's own significance threshold — these are not
raw/unfiltered scores), but it fails dbCAN's own ≥2-tool consensus rule. Top gray-zone
families: GT2 (103.5k), AA3_2 (88.4k), GH3 (87.6k), GT1 (74.6k), GH47 (59.3k). Gray-zone
fraction varies by taxonomic class, highest in Xylonomycetes/Eurotiomycetes/Sordariomycetes
(~13–14%) and lowest in Glomeromycotina/Microsporidia (~2%).

## 2. Structure-tier validation of the gray zone

A stratified sample (1,500 gray-zone + 500 high-confidence-CAZyme + 500
high-confidence-non-CAZyme controls, 2,483 with retrievable sequences) was scored with an
independent structure-similarity signal that never sees the sequence-evidence tier label:

1. **ProstT5** (Rostlab/ProstT5) predicts 3Di structure tokens directly from the amino-acid
   sequence — no folding required. Run across all 8 GPUs on met, ~15–35s/protein.
2. **Foldseek/mmseqs2 3Di search** — the predicted 3Di string is aligned against
   CAZyme3D_id50's real-structure-derived 3Di sequences (178,356 structures). 88.1% of
   queries got a hit.
3. **SaProt** (westlake-repl/SaProt_650M_AF2) embeds the AA+3Di structure-aware sequence into
   a 1280-dim vector; cosine similarity to a CAZyme3D_id50 reference centroid gives a second,
   orthogonal structural signal.
4. The two signals combine into `structure_evidence_score` (0–1, higher = more CAZyme-like).

**Result — the score cleanly separates the two controls, with gray-zone in between.** This table
reflects the FINAL, corrected run: ProstT5 + the fixed (reservoir-sampled) `structure_evidence_score.py`
scored against Track A's actual 4,000-protein `sample_for_structure.fasta`/`.tsv` handoff (2,000
gray-zone / 1,200 high-confidence-CAZyme / 800 high-confidence-non-CAZyme; 3,976/4,000 sequences
scored, 24 dropped on length/fetch):

| Known tier (structure-blind) | n | mean score | median score |
|---|---:|---:|---:|
| `high_confidence_cazyme` | 1,200 | 0.625 | 0.620 |
| `gray_zone` | 2,000 | 0.581 | 0.576 |
| `high_confidence_non_cazyme` | 776 | 0.460 | 0.456 |

This confirms the gray zone is a genuine mixture, not uniformly one class — structure evidence
gives real discriminating power. Using the midpoints between adjacent tier means as decision
thresholds (0.603 and 0.520) to adjudicate the 2,000 sampled gray-zone proteins:

| Adjudicated call | n | % of sampled gray zone |
|---|---:|---:|
| `gray_zone_structure_supports_cazyme` | 844 | 42.2% |
| `gray_zone_structure_supports_non_cazyme` | 731 | 36.6% |
| `gray_zone_ambiguous_structure` | 425 | 21.3% |

i.e. on this sample, over two-fifths of single-tool CAZy-family hits are corroborated by
structure (real, currently-missed CAZymes — candidates for dbCAN4's recall story), a bit over a
third look like sequence-level false positives once structure is considered, and the remainder
(a fifth) stay genuinely ambiguous even with an orthogonal signal — a legitimate abstention
population, not a labeling failure. (An earlier, in-progress version of this report used a
smaller, independently-drawn 2,483-protein sample as a documented stand-in while this corrected
4,000-protein run was still executing; those numbers — 33.4%/38.3%/28.3% — pointed the same
direction but are superseded by the table above.)

**Resolved caveat:** the SaProt-embedding component of `structure_evidence_score` previously used
a reference centroid built from the first 2,000 rows (file order, not a representative draw) of
the 178,356-row CAZyme3D_id50 3Di table. The script was fixed to reservoir-sample uniformly, and
the numbers in this report are from that fixed version scored against the correct handoff file —
this is the final, authoritative result for this analysis.

## 3. CAZyme3D_id50, ProstT5, SaProt, ESM Atlas 2, AF3db — what each is actually good for

| Resource | Role established this session |
|---|---|
| **CAZyme3D_id50** (178,356 structures, already on met) | Not a direct structure source for Mycocosm proteins by accession (6/3.5M exact-sequence matches) — but a strong **homology reference set**: 93.3% of a CAZyme sample has a ≥30%-identity homolog in it. This is the reference database the structure tier scores against. |
| **ProstT5** | Sequence→3Di prediction with no folding step — the practical way to get every gray-zone protein into "structure space" cheaply. Validated at scale (2,483 proteins, 8 GPUs). |
| **SaProt** (650M-AF2) | AA+3Di structure-aware embedding, orthogonal signal to Foldseek/mmseqs2 alignment. Validated; centroid-sampling bug identified and fixed (rerun pending, see caveat above). |
| **ESM Atlas 2** | 85.9% sequence-hit coverage on a CAZyme sample, but median identity only 51.2% — mostly distant homology, consistent with its MGnify (bacterial/archaeal/environmental) origin. Limited marginal value over CAZyme3D_id50 for this project. |
| **AF3db** | Only reachable via UniProt-format accessions; Mycocosm JGI IDs aren't in that namespace, so coverage was checked indirectly through CAZyme3D-homolog accessions (75.9% of the UniProt-queryable subset resolve). Not directly usable without a RefSeq→UniProt mapping step. |
| **Local ESMFold folding** | Remains the most practical route for gray-zone proteins with no adequate CAZyme3D/AF3db/ESM-Atlas hit — already proven on met (~22s/protein median across 8 GPUs). |

**Recommended structure-tier stack going forward:** CAZyme3D_id50 (reference) + ProstT5
(cheap 3Di for every protein) + SaProt (orthogonal embedding signal) + local ESMFold (only for
proteins lacking any adequate reference hit). ESM Atlas 2 and AF3db are not worth further
integration investment at this time.

## 4. Deliverables

- **Gray-zone adjudicated dataset (final, authoritative)**:
  [gray_zone_adjudicated_FINAL_corrected_4000sample.tsv](gray_zone_adjudicated_FINAL_corrected_4000sample.tsv)
  — the true Track A/Track B merge on Track A's designated 4,000-protein `sample_for_structure`
  handoff (2,000 gray-zone/1,200 high-conf-CAZyme/800 high-conf-non-CAZyme, 3,976 scored), with
  adjudicated_tier, all sequence-tiering columns, and all structure-evidence columns
  (foldseek/SaProt/structure_evidence_score). Supersedes an earlier 2,483-protein documented
  interim substitute used while this run was still executing.
- **Full sequence-evidence-tiered population** (28.2M rows, Track A): see
  `tiered_proteins_part{0..3}of4.tsv.gz` artifacts from Track A, plus
  `summary_overall.json`, `summary_tier_by_class.tsv`, `summary_gray_zone_families.tsv`.
- **Structure tier status report**: `track_b_structure_tier_report.md` (Track A/B originals).
- **Scripts** (committed to `dbcan4-advanced/scripts/`): `build_tiered_dataset.py` (Track A
  parsing/tiering), `cazyme3d_mapping.py`, `prostt5_validate.py`, `saprot_embed.py`,
  `foldseek_util.py`, `esmatlas_coverage.py`, `structure_evidence_score.py`,
  `extract_validation_sample.py`, `fetch_sample_sequences.py`.
- **Updated design doc**: `design_dbcan4_advanced.md` §8 (new).

## 5. Open follow-ups

1. Rerun `structure_evidence_score.py` with the fixed reservoir-sampled SaProt reference
   centroid and regenerate the adjudicated TSV before using the SaProt-embedding component for
   downstream training decisions.
2. Scale structure-evidence scoring from the 2,483-protein validation sample to the full
   2.84M-protein gray zone (or a larger stratified sample) once the rerun is validated.
3. Build a Mycocosm JGI→UniProt/RefSeq ID-mapping step if AF3db coverage is worth a second look.
4. Feed the `gray_zone_structure_supports_cazyme` population into the review's gap-5
   recommendation (cluster the abstained/gray pool to propose candidate new CAZy families).
