
import re
def clean(seq):
    seq = seq.upper().rstrip('*')
    return re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', 'X', seq)

pid=None; buf=[]
out = open('novel_family_work/candidate_pool_v2.clean.fasta','w')
def flush(pid, buf):
    if pid is None: return
    seq = clean(''.join(buf))
    if seq:
        out.write(f'>{pid}\n{seq[:1500]}\n')
with open('novel_family_work/candidate_pool_v2.fasta') as f:
    for line in f:
        if line.startswith('>'):
            flush(pid, buf)
            pid = line[1:].strip().split()[0]
            buf=[]
        else:
            buf.append(line.strip())
    flush(pid, buf)
out.close()
print('done')

