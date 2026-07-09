
import numpy as np, csv, json

d = np.load('novel_family_work/emb_candidate_v2/cand.shard0.npz', allow_pickle=True)
ids = d['ids']; emb = d['emb'].astype(np.float32)
print('candidate emb', emb.shape)

# L2 normalize
norm = np.linalg.norm(emb, axis=1, keepdims=True); norm[norm==0]=1
embn = emb / norm

import umap
import hdbscan

reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=10, metric='cosine', random_state=42)
emb_umap = reducer.fit_transform(embn)
print('umap done', emb_umap.shape)

clusterer = hdbscan.HDBSCAN(min_cluster_size=4, min_samples=2, metric='euclidean', cluster_selection_method='eom')
labels = clusterer.fit_predict(emb_umap)
print('n clusters (excl noise):', len(set(labels)) - (1 if -1 in labels else 0))
print('n noise:', (labels==-1).sum(), '/', len(labels))

# also 2D umap for plotting
reducer2d = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, metric='cosine', random_state=42)
emb_umap2d = reducer2d.fit_transform(embn)

np.savez('novel_family_work/candidate_v2_clusters.npz', ids=ids, labels=labels,
         umap2d=emb_umap2d, umap10d=emb_umap, prob=clusterer.probabilities_)
print('saved')

