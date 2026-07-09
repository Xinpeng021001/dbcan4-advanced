# dbCAN4 Track A — Multi-Evidence CAZyme/non-CAZyme Confidence Tiering

## Purpose

dbCAN's current release logic calls a protein a CAZyme only when **≥2 of 3 tools**
(DIAMOND-vs-CAZy, dbCAN HMM, dbCAN-sub HMM) hit it (the `#ofTools` / `Recommend
Results` columns in `overview.tsv`). That binary cutoff throws away information:
a protein with a single strong DIAMOND hit at 40% identity to a genuine CAZy
enzyme is treated identically to a protein with *zero* hits from any tool, even
though the two have very different prior probabilities of really being CAZymes.
This is the "gray zone" this track quantifies and stratifies, so it can be
handed to structure-based validation (Track B, ESMFold + Foldseek) and
eventually re-integrated as training/eval data for dbCAN4.

## Data basis (confirmed by direct inspection on met, 2026-07-09)

- 2,226 genome directories across 28 fungal taxonomic classes under
  `/array1/xinpeng/all_genome/<TaxClass>/<GenomeID>/`, ~40 GB, ~28.2M total
  proteins (`uniInput.faa`, excluding the redundant `cgc_cazyme_only/uniInput.faa`
  duplicate).
- `overview.tsv` **only lists proteins with ≥1 hit from ≥1 tool** — in the
  sampled genome (Abobi1, Agaricomycetes), 1,476 of 11,987 proteins appear in
  `overview.tsv`; the other 10,511 never appear there at all and are only
  recoverable by set-differencing `uniInput.faa` IDs against `overview.tsv` IDs.
  These are the "zero-hit" proteins.
- The **raw per-tool output files** (`diamond.out`, `dbCAN_hmm_results.tsv`,
  `dbCANsub_hmm_results.tsv`) are already filtered to each tool's own internal
  significance threshold (every row is a real hit by that tool's own standard —
  verified: every Gene ID with a DIAMOND row in `diamond.out` matches exactly the
  set of Gene IDs with `DIAMOND != "-"` in `overview.tsv`, and same for the two
  HMM tools). So a **single-tool hit is genuine homology evidence** — it simply
  fails dbCAN's own ≥2-tool consensus bar. This is the direct empirical basis
  for treating "1-tool-hit" proteins as an intermediate confidence class rather
  than merging them into either extreme.
- Gene ID separator differs by file: FASTA headers (`uniInput.faa`,
  `non_CAZyme.faa`) and `diamond.out`'s Gene ID column use `|`
  (`jgi|Abobi1|105364|CE105363_8952`); `overview.tsv`, `dbCAN_hmm_results.tsv`,
  `dbCANsub_hmm_results.tsv` use `-` (`jgi-Abobi1-105364-CE105363_8952`). All IDs
  are canonicalized to the dash form in the aggregation script.
- `diamond.out`'s `CAZy ID` field is `<CAZy_seq_id>|<source_fasta>|<Family>` —
  family is the token after the **last** `|`.
- A protein can have multiple rows in a raw per-tool file (multiple domain hits);
  the aggregation script keeps the single best-scoring row per protein per tool
  (lowest i-Evalue for the two HMM-based tools, lowest E-value for DIAMOND).

## Per-tool evidence components

For each tool we compute a bounded `[0, 1]` component score that combines
statistical significance and alignment coverage/identity, so tools with
different score scales (E-value vs i-Evalue vs %identity) become comparable:

```
strength(evalue) = clip(-log10(evalue), 0, 20) / 20     # saturates at evalue <= 1e-20

hmm_component      = strength(dbCAN_hmm i-Evalue)   * dbCAN_hmm coverage
sub_component      = strength(dbCAN-sub i-Evalue)   * dbCAN-sub coverage
diamond_component  = strength(DIAMOND E-value)      * (DIAMOND %identity / 100)

evidence_score = max(hmm_component, sub_component, diamond_component)
```

`coverage` (HMM-based tools) and `%identity` (DIAMOND) act as a penalty for
short/partial-domain or low-identity hits that pass the E-value bar only by
alignment-length inflation. `evidence_score` is the strongest single signal
across tools — it is what actually distinguishes a "confident 1-tool hit" from
a "marginal 1-tool hit" inside the gray zone (see subtiers below).

## Tier definitions

Tiering is driven primarily by `n_tools` (recomputed independently from the
per-tool best-hit tables — not merely copied from `overview.tsv`'s `#ofTools`,
though the two agree by construction since `overview.tsv` also uses per-tool
`!= "-"` columns as its source), with `evidence_score` used only to split the
1-tool bucket into subtiers:

| Tier | Condition | Subtier | Meaning |
|---|---|---|---|
| **high_confidence_cazyme** | `n_tools >= 2` | `1A_all_tools_agree` (n_tools=3) | All three tools (DIAMOND, dbCAN HMM, dbCAN-sub) agree — dbCAN's own strongest consensus class. |
| | | `1B_two_tools_agree` (n_tools=2) | Two of three tools agree — dbCAN's current release-quality CAZyme call. |
| **gray_zone** | `n_tools == 1` | `2A_gray_high` (`evidence_score >= 0.5`) | Single strong tool hit (high identity/significance + good coverage) that only just misses the 2-tool consensus bar — most likely true CAZymes, prime candidates for structural confirmation. |
| | | `2B_gray_low` (`evidence_score < 0.5`) | Single weak/marginal tool hit — plausible remote homologs, pseudo-enzymes, or noise; the least certain gray-zone proteins. |
| **high_confidence_non_cazyme** | `n_tools == 0` | `3_zero_hits` | No hit from any of the three tools at any threshold — the population dbCAN would currently call "not a CAZyme" with no ambiguity. |

A **reserved column**, `structure_evidence_score`, is carried through the whole
schema (initialized to null/blank) so that Track B's ESMFold+Foldseek
structure-similarity results can be merged back in by `protein_id` without any
schema change. It is intentionally *not* used in the current tier computation —
tiering here is sequence/profile-evidence only; a later phase can re-tier using
`structure_evidence_score` as a fourth input (e.g. promote a `2B_gray_low`
protein to high-confidence if it structurally aligns well to a CAZy fold, or
demote a `1B_two_tools_agree` call that fails to fold into any known CAZyme
architecture).

## Output schema (`tiered_proteins.tsv.gz`, one row per protein)

```
protein_id, genome_id, tax_class, seq_length,
ec_numbers,
dbcan_hmm_family, dbcan_hmm_ievalue, dbcan_hmm_coverage,
dbcan_sub_family, dbcan_sub_ievalue, dbcan_sub_coverage, dbcan_sub_substrate,
diamond_family, diamond_pident, diamond_evalue, diamond_bitscore, diamond_cazy_hit_id,
n_tools, recommend_family,
hmm_component, sub_component, diamond_component, evidence_score,
tier, subtier,
structure_evidence_score   # placeholder, null until Track B merges in
```

## Sampling strategy for structure validation

Per-(tax_class, tier) reservoir sampling (algorithm R, cap configurable, default
400/stratum) keeps memory bounded while giving every taxonomic class fair
representation, then a final stratified draw targets a fixed total sample size
(default 4,000) with tier weights **gray_zone 50% / high_confidence_cazyme 30%
/ high_confidence_non_cazyme 20%** — gray-zone proteins are the primary object
of interest for structural confirmation, but confident positives and confident
negatives are retained in the sample as calibration/negative controls for
Track B's structure-similarity scoring.

## Known limitations

- Tiering uses only the tool evidence already computed by the existing dbCAN
  run (no new HMM/DIAMOND search was re-run); it is a re-scoring/re-thresholding
  of existing outputs, not a re-annotation.
- `evidence_score`'s 0.5 split point for `2A`/`2B` is a heuristic threshold, not
  fit against a labeled gold-standard gray-zone set — it should be revisited
  once Track B's structural evidence is available to check calibration (e.g.
  does `2A_gray_high` actually structurally resemble CAZy folds more often than
  `2B_gray_low`?).
- `recommend_family` / per-tool family columns can disagree with each other in
  edge cases (e.g. DIAMOND calls a different family from dbCAN-sub); the gray
  zone family distribution summary uses `recommend_family` when present, else
  falls back to `dbcan_hmm_family` → `dbcan_sub_family` → `diamond_family` in
  that order, and only the family before any `+`-joined multi-family suffix.
