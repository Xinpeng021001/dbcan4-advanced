#!/usr/bin/env python3
"""
Build a labeled CAZyme reference FASTA for ESM-C embedding, from any CAZy release
FASTA (fungal-only or all-kingdom). Applies a temporal-cleanliness guard and exact
duplicate removal so the reference is comparable across taxonomic scope.

Header layouts handled (same as build_reference.py):
    NCBI:  >ACCESSION|FAM[|FAM...][|EC]                 id = field[0]
    JGI24: >NUMID|Genome...aa.fasta|FAM[|FAM...]        id = field before .fasta
    JGI25: >FAM[|FAM...]|NUMID|Genome...aa.fasta        id = field before .fasta
A family is any '|'-field matching ^(GH|GT|PL|CE|AA|CBM)\\d+(_\\d+)?$.

Steps:
  1. Read --exclude-faa (the fixed eval_2025 set); collect md5 of each sequence
     (uppercased, trailing '*' stripped) -> EVAL_MD5.
  2. Stream --in-faa. For each record with >=1 family:
       - compute md5; if md5 in EVAL_MD5  -> SKIP  (temporal / self-match leak guard)
       - if md5 already emitted           -> SKIP  (exact-duplicate removal)
       - else emit  >ID|fam1,fam2,...  and the sequence.
  3. Write <out-prefix>.faa and <out-prefix>_labels.tsv (protein_id, seq_md5, families).

The eval guard matters for the all-kingdom reference: a 2025 fungal eval sequence can
appear verbatim in a 2024 non-fungal genome; keeping it would let kNN retrieve a
cosine=1.0 self-match (memorization), the same leak we criticized for DIAMOND. We
exclude exact matches so the fungi-vs-all-kingdom contrast measures homology-based
retrieval, not exact-duplicate recall. Applied identically to both references.
"""
import argparse, hashlib, re, sys, time

FAM_RE = re.compile(r'^(GH|GT|PL|CE|AA|CBM)\d+(_\d+)?$')

def header_id(fields):
    for i, f in enumerate(fields):
        if ".fasta" in f:
            return fields[i-1] if i > 0 else fields[0]
    return fields[0]

def header_families(fields):
    seen, out = set(), []
    for f in fields:
        if FAM_RE.match(f) and f not in seen:
            seen.add(f); out.append(f)
    return out

def iter_fasta(path):
    hid, buf = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if hid is not None:
                    yield hid, "".join(buf)
                hid = line[1:].rstrip("\n"); buf = []
            else:
                buf.append(line.strip())
        if hid is not None:
            yield hid, "".join(buf)

def md5_seq(seq):
    return hashlib.md5(seq.upper().rstrip("*").encode()).hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-faa", required=True)
    ap.add_argument("--exclude-faa", required=True, help="eval FASTA whose sequences must NOT appear in the reference")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-len", type=int, default=30)
    ap.add_argument("--max-len", type=int, default=1500)
    args = ap.parse_args()

    t0 = time.time()
    eval_md5 = set()
    for _, seq in iter_fasta(args.exclude_faa):
        if seq:
            eval_md5.add(md5_seq(seq))
    print(f"[guard] {len(eval_md5):,} eval sequence md5s to exclude", file=sys.stderr, flush=True)

    n_in = n_lab = n_excl_eval = n_dup = n_len = n_out = 0
    emitted = set()
    with open(args.out_prefix + ".faa", "w") as fo, \
         open(args.out_prefix + "_labels.tsv", "w") as fl:
        fl.write("protein_id\tseq_md5\tfamilies\n")
        for hid, seq in iter_fasta(args.in_faa):
            n_in += 1
            if n_in % 500000 == 0:
                print(f"  ...{n_in:,} read, {n_out:,} kept ({time.time()-t0:.0f}s)", file=sys.stderr, flush=True)
            fields = hid.split("|")
            fams = header_families(fields)
            if not fams:
                continue
            n_lab += 1
            if not seq:
                continue
            L = len(seq.rstrip("*"))
            if L < args.min_len or L > args.max_len:
                n_len += 1; continue
            h = md5_seq(seq)
            if h in eval_md5:
                n_excl_eval += 1; continue
            if h in emitted:
                n_dup += 1; continue
            emitted.add(h)
            pid = header_id(fields)
            fo.write(f">{pid}|{','.join(fams)}\n{seq}\n")
            fl.write(f"{pid}\t{h}\t{','.join(fams)}\n")
            n_out += 1

    print(f"[done] read={n_in:,} labeled={n_lab:,} kept={n_out:,} "
          f"| excluded: eval_leak={n_excl_eval:,} exact_dup={n_dup:,} len_filter={n_len:,} "
          f"({time.time()-t0:.0f}s)", file=sys.stderr, flush=True)
    # machine-readable tail
    print(f'{{"kept":{n_out},"labeled":{n_lab},"excl_eval":{n_excl_eval},"excl_dup":{n_dup},"excl_len":{n_len}}}')

if __name__ == "__main__":
    main()
