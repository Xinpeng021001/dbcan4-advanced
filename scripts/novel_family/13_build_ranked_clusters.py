
import numpy as np, csv, json
from collections import defaultdict, Counter

# load calibration
dcal = np.load('novel_family_work/calibration_novelty_sims.npz')
novel_family_sims = dcal['novel_family_sims']
novel_seq_sims = dcal['novel_seq_sims']

def percentile_rank(x, arr):
    return float((arr < x).mean() * 100)

# load merged table
rows = []
with open('novel_family_work/candidate_v2_merged.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        rows.append(row)
print('total rows', len(rows))

by_cluster = defaultdict(list)
for row in rows:
    by_cluster[int(row['cluster'])].append(row)

# structural coherence
with open('novel_family_work/cluster_structural_coherence.json') as f:
    coherence = json.load(f)

results = []
for lab, members in by_cluster.items():
    if lab == -1:
        continue
    n = len(members)
    sims_fam = [float(m['cosine_sim_to_nearest_family']) for m in members if m['cosine_sim_to_nearest_family']]
    sims_0fam = [float(m['cosine_sim_to_0fam']) for m in members if m['cosine_sim_to_0fam']]
    struct_scores = [float(m['structure_evidence_score']) for m in members if m['structure_evidence_score']]
    taxa = Counter(m['tax_class'] for m in members)
    genomes = Counter(m['genome_id'] for m in members)
    nearest_fams = Counter(m['nearest_known_family_esmc'] for m in members)
    fam0_classes = Counter(m['best_0fam_class'] for m in members)
    diamond_fams = Counter(m['diamond_family'] for m in members if m['diamond_family'] not in ('-',''))

    mean_sim_fam = np.mean(sims_fam) if sims_fam else None
    mean_sim_0fam = np.mean(sims_0fam) if sims_0fam else None
    mean_struct = np.mean(struct_scores) if struct_scores else None

    pct_vs_novelfam = percentile_rank(mean_sim_fam, novel_family_sims) if mean_sim_fam else None
    pct_vs_novelseq = percentile_rank(mean_sim_fam, novel_seq_sims) if mean_sim_fam else None

    coh = coherence.get(str(lab), {})

    top_fam0, top_fam0_n = fam0_classes.most_common(1)[0]
    frac_fam0_agree = top_fam0_n / n

    results.append({
        'cluster': lab, 'n': n,
        'n_genomes': len(genomes), 'n_tax_classes': len(taxa),
        'top_tax_class': taxa.most_common(1)[0][0], 'top_tax_class_frac': taxa.most_common(1)[0][1]/n,
        'mean_structure_evidence_score': mean_struct,
        'mean_3di_pident_intra': coh.get('mean_3di_pident', 0),
        'median_3di_pident_intra': coh.get('median_3di_pident', 0),
        'n_intra_3di_hits': coh.get('n_intra_hits', 0),
        'nearest_known_family_top': nearest_fams.most_common(1)[0][0],
        'nearest_known_family_top_frac': nearest_fams.most_common(1)[0][1]/n,
        'mean_cosine_to_nearest_family': mean_sim_fam,
        'pctile_vs_eval2025_novel_family_calib': pct_vs_novelfam,
        'pctile_vs_eval2025_novel_seq_calib': pct_vs_novelseq,
        'mean_cosine_to_cazy_0fam': mean_sim_0fam,
        'top_cazy_0fam_class': top_fam0, 'frac_agree_top_0fam_class': frac_fam0_agree,
        'diamond_family_hits': dict(diamond_fams),
        'representative_protein': members[0]['protein_id'],
    })

results.sort(key=lambda r: -r['n'])
with open('novel_family_work/ranked_clusters_raw.json','w') as f:
    json.dump(results, f, indent=2, default=str)
print('n clusters', len(results))
for r_ in results[:10]:
    print(r_['cluster'], r_['n'], 'coh=%.1f'%r_['mean_3di_pident_intra'], 'sim_fam=%.3f'%r_['mean_cosine_to_nearest_family'],
          'pct_vs_novelfam=%.1f'%r_['pctile_vs_eval2025_novel_family_calib'], 'top_fam0=%s(%.2f)'%(r_['top_cazy_0fam_class'], r_['frac_agree_top_0fam_class']),
          'tax=%s(%.2f)'%(r_['top_tax_class'], r_['top_tax_class_frac']))

