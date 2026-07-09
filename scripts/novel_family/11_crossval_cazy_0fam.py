
import numpy as np, glob, csv, re
from collections import defaultdict, Counter

# clean family map from tsv (single column already GH0/GT0/etc)
fam_map = {}
with open('structure/novel_family/cazydb_0fam_entries.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        fam_map[row['protein_id']] = row['family']

fam0_ids, fam0_emb = [], []
for f in sorted(glob.glob('novel_family_work/emb_cazy0fam/c0.shard*.npz')):
    d = np.load(f, allow_pickle=True)
    fam0_ids.append(d['ids']); fam0_emb.append(d['emb'])
fam0_ids = np.concatenate(fam0_ids)
fam0_emb = np.concatenate(fam0_emb).astype(np.float32)
norm = np.linalg.norm(fam0_emb, axis=1, keepdims=True); norm[norm==0]=1
fam0_embn = fam0_emb / norm

fam0_class = np.array([fam_map.get(str(pid), 'UNK') for pid in fam0_ids])
print('0fam total', fam0_embn.shape)
print('class counts:', Counter(fam0_class))

d = np.load('novel_family_work/emb_candidate_v2/cand.shard0.npz', allow_pickle=True)
cand_ids = d['ids']; cand_emb = d['emb'].astype(np.float32)
norm = np.linalg.norm(cand_emb, axis=1, keepdims=True); norm[norm==0]=1
cand_embn = cand_emb / norm

dcl = np.load('novel_family_work/candidate_v2_clusters.npz', allow_pickle=True)
assert list(dcl['ids']) == list(cand_ids)
labels = dcl['labels']

sims = cand_embn @ fam0_embn.T
best_idx = np.argmax(sims, axis=1)
best_sim = sims[np.arange(len(cand_ids)), best_idx]
best_0fam_id = fam0_ids[best_idx]
best_0fam_class = fam0_class[best_idx]

with open('novel_family_work/candidate_v2_0fam_nn.tsv','w') as out:
    out.write('protein_id\tcluster\tbest_0fam_hit\tbest_0fam_class\tcosine_sim\n')
    for i in range(len(cand_ids)):
        out.write(f'{cand_ids[i]}\t{labels[i]}\t{best_0fam_id[i]}\t{best_0fam_class[i]}\t{best_sim[i]:.4f}\n')

by_cluster = defaultdict(list)
for i in range(len(cand_ids)):
    by_cluster[labels[i]].append(i)

print()
print('cluster\tn\tmean_sim_to_0fam\tmax_sim\ttop_0fam_class')
rows=[]
for lab, idxs in sorted(by_cluster.items(), key=lambda x:-len(x[1])):
    if lab==-1: continue
    s = best_sim[idxs]
    classes = best_0fam_class[idxs]
    top_class, top_count = Counter(classes).most_common(1)[0]
    rows.append((lab,len(idxs),s.mean(),s.max(),top_class,top_count))
    print(f"{lab}\t{len(idxs)}\t{s.mean():.4f}\t{s.max():.4f}\t{top_class}({top_count}/{len(idxs)})")

import json
with open('novel_family_work/cluster_0fam_crossval.json','w') as f:
    json.dump([{'cluster':int(r[0]),'n':int(r[1]),'mean_sim_0fam':float(r[2]),'max_sim_0fam':float(r[3]),
                'top_0fam_class':str(r[4]),'top_0fam_class_count':int(r[5])} for r in rows], f, indent=2)
print('done')

