
import json, csv
import numpy as np

with open('novel_family_work/ranked_clusters_raw.json') as f:
    results = json.load(f)

# combined score: reward (a) low pct_vs_novelfam (distant from known -> more novel),
# (b) high structural coherence (mean_3di_pident_intra),
# (c) high frac_agree_top_0fam_class (cross-validation support),
# (d) cluster size (more evidence)
# normalize each component 0-1 across clusters, weighted sum
def minmax(vals):
    vals = np.array(vals, dtype=float)
    lo, hi = vals.min(), vals.max()
    if hi - lo < 1e-9:
        return np.zeros_like(vals)
    return (vals - lo) / (hi - lo)

novelty_raw = np.array([100 - r['pctile_vs_eval2025_novel_family_calib'] for r in results])  # higher = more distant from known fam
coherence_raw = np.array([r['mean_3di_pident_intra'] for r in results])
crossval_raw = np.array([r['frac_agree_top_0fam_class'] for r in results])
size_raw = np.array([np.log1p(r['n']) for r in results])

novelty_n = minmax(novelty_raw)
coherence_n = minmax(coherence_raw)
crossval_n = minmax(crossval_raw)
size_n = minmax(size_raw)

W_NOVELTY, W_COHERENCE, W_CROSSVAL, W_SIZE = 0.35, 0.30, 0.25, 0.10
combined = W_NOVELTY*novelty_n + W_COHERENCE*coherence_n + W_CROSSVAL*crossval_n + W_SIZE*size_n

for i, r in enumerate(results):
    r['novelty_score_0to100'] = float(novelty_raw[i])
    r['combined_evidence_score'] = float(combined[i])

results.sort(key=lambda r: -r['combined_evidence_score'])

# write final ranked tsv
cols = ['rank','cluster','n','combined_evidence_score','novelty_score_0to100',
        'mean_cosine_to_nearest_family','nearest_known_family_top','nearest_known_family_top_frac',
        'mean_3di_pident_intra','n_intra_3di_hits',
        'mean_structure_evidence_score',
        'top_cazy_0fam_class','frac_agree_top_0fam_class','mean_cosine_to_cazy_0fam',
        'n_genomes','n_tax_classes','top_tax_class','top_tax_class_frac',
        'representative_protein']
with open('novel_family_work/ranked_candidate_clusters.tsv','w') as out:
    out.write('\t'.join(cols)+'\n')
    for i, r in enumerate(results):
        r['rank'] = i+1
        out.write('\t'.join(str(r.get(c,'')) for c in cols)+'\n')

with open('novel_family_work/ranked_candidate_clusters_full.json','w') as f:
    json.dump(results, f, indent=2, default=str)

print('Top 15 ranked clusters:')
print('rank\tcluster\tn\tcombined\tnovelty\tsim_fam\tcoherence\t0fam_class(frac)\ttax(frac)')
for i, r in enumerate(results[:15]):
    print(f"{i+1}\t{r['cluster']}\t{r['n']}\t{r['combined_evidence_score']:.3f}\t{r['novelty_score_0to100']:.1f}\t"
          f"{r['mean_cosine_to_nearest_family']:.3f}\t{r['mean_3di_pident_intra']:.1f}\t"
          f"{r['top_cazy_0fam_class']}({r['frac_agree_top_0fam_class']:.2f})\t{r['top_tax_class']}({r['top_tax_class_frac']:.2f})")

