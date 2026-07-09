
import csv

want = set()
with open('novel_family_work/candidate_pool_v2_ids.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        want.add(row['protein_id'])

out = open('novel_family_work/candidate_v2_3di.fasta','w')
n=0
with open('structure/validation_sample/all_prostt5.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        if row['protein_id'] in want:
            di3 = row['di3_string']
            out.write(f">{row['protein_id']}\n{di3}\n")
            n+=1
out.close()
print('wrote', n)

