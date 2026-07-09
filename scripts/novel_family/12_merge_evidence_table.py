
import csv
from collections import defaultdict

# load all pieces
meta = {}
with open('novel_family_work/candidate_pool_v2_meta.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        meta[row['protein_id']] = row

scores = {}
with open('structure/validation_sample/scores/structure_evidence_scores_final.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        scores[row['protein_id']] = row

nearest_fam = {}
with open('novel_family_work/candidate_v2_nearest_family.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        nearest_fam[row['protein_id']] = row

nn0fam = {}
with open('novel_family_work/candidate_v2_0fam_nn.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        nn0fam[row['protein_id']] = row

out = open('novel_family_work/candidate_v2_merged.tsv','w')
cols = ['protein_id','cluster','tier','genome_id','tax_class','seq_length',
        'dbcan_hmm_family','dbcan_sub_family','diamond_family',
        'structure_evidence_score','foldseek_3di_best_hit','foldseek_3di_pident','saprot_cosine_to_cazyme_centroid',
        'nearest_known_family_esmc','cosine_sim_to_nearest_family',
        'best_0fam_hit','best_0fam_class','cosine_sim_to_0fam']
out.write('\t'.join(cols)+'\n')
n=0
for pid, m in meta.items():
    s = scores.get(pid, {})
    nf = nearest_fam.get(pid, {})
    n0 = nn0fam.get(pid, {})
    row = [pid, nf.get('cluster',''), m.get('tier',''), m.get('genome_id',''), m.get('tax_class',''), m.get('seq_length',''),
           m.get('dbcan_hmm_family','-'), m.get('dbcan_sub_family','-'), m.get('diamond_family','-'),
           s.get('structure_evidence_score',''), s.get('foldseek_3di_best_hit',''), s.get('foldseek_3di_pident',''),
           s.get('saprot_cosine_to_cazyme_centroid',''),
           nf.get('nearest_known_family',''), nf.get('cosine_sim_to_nearest_family',''),
           n0.get('best_0fam_hit',''), n0.get('best_0fam_class',''), n0.get('cosine_sim',''),
           ]
    out.write('\t'.join(str(x) for x in row)+'\n')
    n+=1
out.close()
print('merged', n)

