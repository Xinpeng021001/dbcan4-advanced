
import csv

scores = {}
with open('structure/validation_sample/scores/structure_evidence_scores_final.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        try:
            s = float(row['structure_evidence_score'])
        except:
            continue
        scores[row['protein_id']] = (row['tier'], s, row.get('foldseek_3di_best_hit','-'), row.get('foldseek_3di_pident','-'))

print('total scored:', len(scores))

THRESH = 0.60
gray_high = [(pid,s) for pid,(t,s,_,_) in scores.items() if t=='gray_zone' and s>=THRESH]
noncaz_high = [(pid,s) for pid,(t,s,_,_) in scores.items() if t=='high_confidence_non_cazyme' and s>=THRESH]
print('gray_zone >= 0.60:', len(gray_high))
print('high_confidence_non_cazyme >= 0.60:', len(noncaz_high))

cand_ids = set(pid for pid,_ in gray_high) | set(pid for pid,_ in noncaz_high)
print('total candidate pool:', len(cand_ids))

with open('novel_family_work/candidate_pool_v2_ids.tsv','w') as out:
    out.write('protein_id\ttier\tstructure_evidence_score\tfoldseek_best_hit\tfoldseek_pident\n')
    for pid,(t,s,hit,pid2) in scores.items():
        if pid in cand_ids:
            out.write(f'{pid}\t{t}\t{s}\t{hit}\t{pid2}\n')
print('written')

