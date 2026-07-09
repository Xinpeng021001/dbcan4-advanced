
import csv, gzip

want = set()
with open('novel_family_work/candidate_pool_v2_ids.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        want.add(row['protein_id'])
print('want', len(want))

out = open('novel_family_work/candidate_pool_v2_meta.tsv','w')
cols = ['protein_id','genome_id','tax_class','seq_length','n_tools','recommend_family',
        'dbcan_hmm_family','dbcan_sub_family','diamond_family','evidence_score','tier','subtier']
out.write('\t'.join(cols)+'\n')
found = 0
with gzip.open('track_a_output/tiered_proteins.tsv.gz','rt') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        if row['protein_id'] in want:
            out.write('\t'.join(row.get(c,'') for c in cols)+'\n')
            found += 1
            if found == len(want):
                break
out.close()
print('found', found)

