
import numpy as np, csv
from collections import defaultdict

d = np.load('novel_family_work/candidate_v2_clusters.npz', allow_pickle=True)
ids = d['ids']; labels = d['labels']

de = np.load('novel_family_work/emb_candidate_v2/cand.shard0.npz', allow_pickle=True)
cand_ids = de['ids']; cand_emb = de['emb'].astype(np.float32)
norm = np.linalg.norm(cand_emb, axis=1, keepdims=True); norm[norm==0]=1
cand_embn = cand_emb / norm

# match order (should already match since same source)
assert list(ids) == list(cand_ids), "id order mismatch"

dc = np.load('novel_family_work/known_family_centroids.npz', allow_pickle=True)
fam_names = dc['fam_names']; cent = dc['cent'].astype(np.float32)

# cosine sim to each known family centroid -> max = nearest known family
sims = cand_embn @ cent.T   # (N, F)
best_idx = np.argmax(sims, axis=1)
best_sim = sims[np.arange(len(cand_ids)), best_idx]
best_fam = fam_names[best_idx]

# per-protein table
with open('novel_family_work/candidate_v2_nearest_family.tsv','w') as out:
    out.write('protein_id\tcluster\tnearest_known_family\tcosine_sim_to_nearest_family\n')
    for i in range(len(cand_ids)):
        out.write(f'{cand_ids[i]}\t{labels[i]}\t{best_fam[i]}\t{best_sim[i]:.4f}\n')

# per-cluster summary
by_cluster = defaultdict(list)
for i in range(len(cand_ids)):
    by_cluster[labels[i]].append(i)

print('cluster\tn\tmean_sim_to_nearest_fam\tmax_sim\tmost_common_nearest_fam')
rows = []
for lab, idxs in sorted(by_cluster.items(), key=lambda x: -len(x[1])):
    if lab == -1:
        continue
    sims_c = best_sim[idxs]
    fams_c = best_fam[idxs]
    from collections import Counter
    top_fam, top_count = Counter(fams_c).most_common(1)[0]
    rows.append((lab, len(idxs), sims_c.mean(), sims_c.max(), top_fam, top_count))
    print(f"{lab}\t{len(idxs)}\t{sims_c.mean():.4f}\t{sims_c.max():.4f}\t{top_fam}({top_count}/{len(idxs)})")

import json
with open('novel_family_work/cluster_nearest_family_summary.json','w') as f:
    json.dump([{'cluster':int(r[0]),'n':int(r[1]),'mean_sim':float(r[2]),'max_sim':float(r[3]),
                'top_fam':str(r[4]),'top_fam_count':int(r[5])} for r in rows], f, indent=2)
print('done')

