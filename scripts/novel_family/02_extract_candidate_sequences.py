
import csv

want = set()
with open('novel_family_work/candidate_pool_v2_ids.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        want.add(row['protein_id'])
print('want', len(want))

got = 0
out = open('novel_family_work/candidate_pool_v2.fasta','w')
pid=None; buf=[]
def flush(pid, buf):
    global got
    if pid in want:
        seq = ''.join(buf)
        out.write(f'>{pid}\n{seq}\n')
        got += 1
with open('structure/validation_sample/sample.faa') as f:
    for line in f:
        if line.startswith('>'):
            flush(pid, buf)
            pid = line[1:].strip().split()[0]
            buf=[]
        else:
            buf.append(line.strip())
    flush(pid, buf)
out.close()
print('extracted', got)

