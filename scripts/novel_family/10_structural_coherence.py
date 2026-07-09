
import numpy as np, csv, json
from collections import defaultdict

dcl = np.load('novel_family_work/candidate_v2_clusters.npz', allow_pickle=True)
ids = dcl['ids']; labels = dcl['labels']
id_to_cluster = dict(zip(ids, labels))
cluster_members = defaultdict(set)
for pid, lab in zip(ids, labels):
    cluster_members[lab].add(pid)

# parse self.m8 hits, keep only within-cluster pairs (excl self)
intra_pident = defaultdict(list)
with open('novel_family_work/foldseek_selfsearch/self.m8') as f:
    for line in f:
        q, t, pident, alnlen, evalue, bits = line.rstrip('\n').split('\t')
        if q == t:
            continue
        pident = float(pident)
        cq = id_to_cluster.get(q)
        ct = id_to_cluster.get(t)
        if cq is not None and cq == ct and cq != -1:
            intra_pident[cq].append(pident)

print('cluster\tn\tn_intra_hits\tmean_3di_pident\tmedian_3di_pident')
coherence = {}
for lab in sorted(set(labels)):
    if lab == -1: continue
    n = (labels==lab).sum()
    hits = intra_pident.get(lab, [])
    if hits:
        mean_p = np.mean(hits); med_p = np.median(hits)
    else:
        mean_p = 0.0; med_p = 0.0
    coherence[int(lab)] = {'n': int(n), 'n_intra_hits': len(hits), 'mean_3di_pident': float(mean_p), 'median_3di_pident': float(med_p)}
    print(f"{lab}\t{n}\t{len(hits)}\t{mean_p:.1f}\t{med_p:.1f}")

with open('novel_family_work/cluster_structural_coherence.json','w') as f:
    json.dump(coherence, f, indent=2)
print('saved')

