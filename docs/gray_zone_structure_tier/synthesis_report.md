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

**Result — the score cleanly separates the two controls, with gray-zone in between:**

| Known tier (structure-blind) | n | mean score | median score |
|---|---:|---:|---:|
| `high_confidence_cazyme` | 500 | 0.653 | 0.647 |
| `gray_zone` | 1,500 | 0.566 | 0.556 |
| `high_confidence_non_cazyme` | 483 | 0.463 | 0.455 |

This confirms the gray zone is a genuine mixture, not uniformly one class — structure evidence
gives real discriminating power. Using the midpoints between adjacent tier means as decision
thresholds (0.609 and 0.514) to adjudicate the 1,500 sampled gray-zone proteins:

| Adjudicated call | n | % of sampled gray zone |
|---|---:|---:|
| `gray_zone_structure_supports_non_cazyme` | 575 | 38.3% |
| `gray_zone_structure_supports_cazyme` | 501 | 33.4% |
| `gray_zone_ambiguous_structure` | 424 | 28.3% |

i.e. on this sample, roughly a third of single-tool CAZy-family hits are corroborated by
structure (real, currently-missed CAZymes — candidates for dbCAN4's recall story), roughly
two-fifths look like sequence-level false positives once structure is considered, and the
remainder stay genuinely ambiguous even with an orthogonal signal — a legitimate abstention
population, not a labeling failure.

**Caveat (carried over from Track B, not yet resolved):** the SaProt-embedding component of
`structure_evidence_score` was computed with a reference centroid built from the first 2,000
rows (file order, not a representative draw) of the 178,356-row CAZyme3D_id50 3Di table. The
script has been fixed to reservoir-sample uniformly, but the numbers above have **not** been
regenerated with the fix. The foldseek/mmseqs2 3Di-search component (component A) already used
the full reference set correctly, and the reported tier-discrimination direction/ranking is
unlikely to reverse, but the exact adjudication counts above should be treated as provisional
until rerun.

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

- **Gray-zone adjudicated dataset**: [gray_zone_adjudicated_structure_validated.tsv](gray_zone_adjudicated_structure_validated.tsv)
  — 2,483 structure-validated proteins (protein_id, sequence tier, adjudicated tier, family,
  foldseek/SaProt evidence columns, structure_evidence_score).
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
