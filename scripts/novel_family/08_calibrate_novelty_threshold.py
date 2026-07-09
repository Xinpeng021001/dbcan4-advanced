
import numpy as np, csv, glob
from collections import defaultdict

# load eval2025 embeddings
eval_ids, eval_fams, eval_emb = [], [], []
for f in sorted(glob.glob('emb/eval2025.shard*.npz')):
    d = np.load(f, allow_pickle=True)
    eval_ids.append(d['ids']); eval_fams.append(d['fams']); eval_emb.append(d['emb'])
eval_ids = np.concatenate(eval_ids)
eval_emb = np.concatenate(eval_emb).astype(np.float32)
norm = np.linalg.norm(eval_emb, axis=1, keepdims=True); norm[norm==0]=1
eval_emb = eval_emb / norm
print('eval total', eval_emb.shape)

# load novelty labels
novelty = {}
with open('data/eval_2025_labels.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        novelty[row['protein_id']] = row['novelty']

dc = np.load('novel_family_work/known_family_centroids.npz', allow_pickle=True)
fam_names = dc['fam_names']; cent = dc['cent'].astype(np.float32)

sims = eval_emb @ cent.T
best_sim = sims.max(axis=1)

novel_family_sims = []
novel_seq_sims = []
for i, pid in enumerate(eval_ids):
    nv = novelty.get(str(pid))
    if nv == 'novel_family':
        novel_family_sims.append(best_sim[i])
    elif nv == 'novel_seq':
        novel_seq_sims.append(best_sim[i])

novel_family_sims = np.array(novel_family_sims)
novel_seq_sims = np.array(novel_seq_sims)
print('novel_family (true new family, ground truth): n=%d mean=%.4f median=%.4f p75=%.4f p90=%.4f max=%.4f'%(
    len(novel_family_sims), novel_family_sims.mean(), np.median(novel_family_sims),
    np.percentile(novel_family_sims,75), np.percentile(novel_family_sims,90), novel_family_sims.max()))
print('novel_seq (known family, new sequence): n=%d mean=%.4f median=%.4f p10=%.4f p25=%.4f min=%.4f'%(
    len(novel_seq_sims), novel_seq_sims.mean(), np.median(novel_seq_sims),
    np.percentile(novel_seq_sims,10), np.percentile(novel_seq_sims,25), novel_seq_sims.min()))

np.savez('novel_family_work/calibration_novelty_sims.npz',
         novel_family_sims=novel_family_sims, novel_seq_sims=novel_seq_sims)

