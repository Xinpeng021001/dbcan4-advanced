# dbCAN4 Domain-Level Contrastive Retrieval: fixing the multidomain failure

## The problem (review priority #1)
Single-label pLM/fusion assigned ONE family to each whole protein. On multidomain
proteins (≥2 distinct CAZyme families) this scored **0.006 exact-set match**: a whole-protein
ESM-C embedding averages its domains into a vector that matches neither. Example — protein
2625 (true CBM91,GH43_14) → whole-protein kNN predicts a single GT1 (purity 0.20).

## The fix
1. **Segment** each eval protein into domains using run_dbcan `overview.tsv` envelope
   coordinates (dbCAN_hmm ∪ dbCAN_sub) — no hmmsearch rerun. 4,938 domains from 4,550 proteins.
2. **Retrieve per domain** against a 307,299-anchor bank of single-family reference-2024
   CAZymes (already domain-level; 92.5% of the reference is single-family).
3. **Train a contrastive projection head** (1152→512→256, cosine-softmax over 193 families
   with ≥20 anchors, 8 epochs) to sharpen family discrimination.
4. **Aggregate** the set of per-domain family calls as the protein's prediction.

## Results (parent-family level, eval-2025 temporal holdout)

### Multidomain proteins (n=330)
| approach | exact-set | Jaccard | per-domain family acc |
|----------|-----------|---------|-----------------------|
| whole-protein kNN (baseline) | **0.000** | 0.477 | — |
| domain-level raw ESM-C kNN | 0.179 (59) | 0.479 | 0.673 |
| **domain-level + trained head** | **0.412 (136)** | **0.685** | **0.897** |

### No regression on the easy cases
- Single-domain proteins (n=4,220): exact-set **0.971**
- **Overall eval exact-set: 0.897 → 0.931** (whole-protein → domain+head) — the gain
  comes entirely from multidomain, with single-domain essentially unchanged.

### Independent validation — CAZyme3D structure truth slice
1,500 CAZyme3D_id50 structures labeled by all-kingdom CAZyDB (accession/MD5), with any
training-MD5 sequence excluded (no leakage):
- top-1 family accuracy **0.779** (head) vs **0.269** (raw ESM-C) — 2.9× gain, confirming
  the head learned real family structure, not eval overfitting.
- top-5 family recall 0.834.

## Honest limits
- **Retrieval ceiling from reference coverage:** 607/1,500 (40.5%) truth-slice structures
  belong to families with NO fungal anchor — structurally unretrievable. This re-confirms
  the main benchmark's "reference must be comprehensive" finding, now for structures.
- **CBM modules are the residual failure:** short accessory CBM domains embed less
  distinctly; most remaining multidomain errors are a correctly-segmented CBM retrieved to
  the wrong family, or a short CBM below dbCAN's HMMER coverage threshold (not segmented at all).
- **Ceiling on segmentation:** 218/332 curated-multidomain proteins have ≥2 domains with
  coordinates in the overview; the other ~114 are single-domain in both HMMER layers (2nd
  family from coordinate-less DIAMOND), so coordinate-based segmentation cannot recover them.

## Design conclusion for dbCAN4
Domain-level contrastive retrieval **replaces whole-protein pLM labeling** for multidomain
architecture and lifts exact-set from ~0 to 0.41 on multidomain / 0.93 overall, with no
single-domain regression. It should be the pLM branch's operating mode in the fusion stack.

## Artifacts
- `domain_retrieval_head.pt` — trained contrastive head (torch state + fam2idx)
- `domain_retrieval_multidomain_eval.tsv` — per-protein WP-vs-domain comparison (330 multidomain)
- `domain_retrieval_summary.json` — all metrics
- inputs: `eval_domains.tsv/.faa`, `eval_domains_esmc_fixed.npz`, `domain_truth_slice.tsv/.faa`,
  `domain_truth_esmc.npz` (saved earlier)
