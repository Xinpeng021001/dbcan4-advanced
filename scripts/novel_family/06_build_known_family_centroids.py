
import numpy as np, glob
from collections import defaultdict

# load reference 2024 embeddings + family labels
ref_ids, ref_fams, ref_embs = [], [], []
for f in sorted(glob.glob('emb/ref2024.shard*.npz')):
    d = np.load(f, allow_pickle=True)
    ref_ids.append(d['ids']); ref_fams.append(d['fams']); ref_embs.append(d['emb'])
ref_ids = np.concatenate(ref_ids)
ref_fams = np.concatenate(ref_fams)
ref_emb = np.concatenate(ref_embs).astype(np.float32)
norm = np.linalg.norm(ref_emb, axis=1, keepdims=True); norm[norm==0]=1
ref_emb = ref_emb / norm
print('ref total', ref_emb.shape)

# build per-family centroid (primary family = first token before comma)
fam_to_idx = defaultdict(list)
for i, fam_field in enumerate(ref_fams):
    fam_field = str(fam_field)
    toks = [x for x in fam_field.replace('|', ',').split(',') if x]
    if not toks:
        continue
    primary = toks[0]
    fam_to_idx[primary].append(i)

centroids = {}
for fam, idxs in fam_to_idx.items():
    if len(idxs) < 2:
        continue
    v = ref_emb[idxs].mean(axis=0)
    v = v / (np.linalg.norm(v)+1e-9)
    centroids[fam] = v

print('n known families with centroid:', len(centroids))
fam_names = list(centroids.keys())
cent_mat = np.vstack([centroids[f] for f in fam_names])
np.savez('novel_family_work/known_family_centroids.npz', fam_names=np.array(fam_names), cent=cent_mat)
print('saved centroids')

