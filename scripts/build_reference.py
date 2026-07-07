#!/usr/bin/env python3
"""
Build the labeled reference (2024) + temporal-holdout test (2025) for dbCAN4-advanced.

Input  (on met, /array1/xinpeng/fungi-cazyme-project/):
    CAZyDB.07142024.fungi.faa   fungal CAZymes, 2024 CAZy release  (reference/train)
    CAZyDB.07242025.fungi.faa   fungal CAZymes, 2025 CAZy release  (test)

CAZy family label lives in the FASTA header. Header layouts handled:
    NCBI:  >ACCESSION|FAM[|FAM...][|EC]                 id = field[0]
    JGI24: >NUMID|Genome...aa.fasta|FAM[|FAM...]          id = field before .fasta
    JGI25: >FAM[|FAM...]|NUMID|Genome...aa.fasta          id = field before .fasta
A family is any '|'-field matching ^(GH|GT|PL|CE|AA|CBM)\\d+(_\\d+)?$ ; this excludes
EC numbers, genome-file fields, and numeric IDs. Families may repeat/multi.

Outputs (written to --out-dir):
    reference_labels_2024.tsv   protein_id \\t seq_len \\t families(csv) \\t classes(csv)
    test_2025_labels.tsv        protein_id \\t seq_len \\t families(csv) \\t classes(csv) \\t novelty
    novelty_summary.json        counts per novelty bucket, family-space diff
    eval_2025.faa               bounded, stratified evaluation FASTA (subset of 2025-novel)
    eval_2025_labels.tsv        labels for eval_2025.faa
    reference_2024.faa          reference FASTA (optionally subsampled per family)

novelty buckets (accession-based, first cut; sequence-identity split done separately):
    carried_over  : accession present in BOTH 2024 and 2025      (memorizable -> excluded from eval)
    novel_seq     : accession only in 2025, ALL families exist in 2024
    novel_family  : accession only in 2025, >=1 family absent from 2024
"""
import argparse, json, re, random, sys
from collections import defaultdict

FAM_RE = re.compile(r'^(GH|GT|PL|CE|AA|CBM)\d+(_\d+)?$')
CLASS_RE = re.compile(r'^(GH|GT|PL|CE|AA|CBM)')

def header_id(fields):
    for i, f in enumerate(fields):
        if ".fasta" in f:
            return fields[i-1] if i > 0 else fields[0]
    return fields[0]

def header_families(fields):
    fams = [f for f in fields if FAM_RE.match(f)]
    # dedupe preserving order
    seen, out = set(), []
    for f in fams:
        if f not in seen:
            seen.add(f); out.append(f)
    return out

def parse_faa(path):
    """Yield (protein_id, families[list], seq_len, seq_md5). Streaming.

    seq_md5 is the md5 of the uppercased sequence with any trailing '*' stripped —
    an accession-scheme-independent identity key (JGI numeric IDs collide across
    genomes and JGI genome-file names carry dates that differ between DB builds,
    so accession matching across years is unreliable; exact-sequence identity is not).
    """
    import hashlib
    pid, fams, buf = None, [], []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if pid is not None:
                    seq = "".join(buf).upper().rstrip("*")
                    yield pid, fams, len(seq), hashlib.md5(seq.encode()).hexdigest()
                fields = line[1:].rstrip("\n").split("|")
                pid = header_id(fields)
                fams = header_families(fields)
                buf = []
            else:
                buf.append(line.strip())
        if pid is not None:
            seq = "".join(buf).upper().rstrip("*")
            yield pid, fams, len(seq), hashlib.md5(seq.encode()).hexdigest()

def load_seqs(path, wanted):
    """Return {pid: seq} for pids in `wanted`."""
    out, pid, keep, buf = {}, None, False, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if pid is not None and keep:
                    out[pid] = "".join(buf)
                fields = line[1:].rstrip("\n").split("|")
                pid = header_id(fields)
                keep = pid in wanted
                buf = []
            elif keep:
                buf.append(line.strip())
        if pid is not None and keep:
            out[pid] = "".join(buf)
    return out

def fam_class(fam):
    m = CLASS_RE.match(fam)
    return m.group(1) if m else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/array1/xinpeng/fungi-cazyme-project")
    ap.add_argument("--out-dir", default="/array1/xinpeng/dbcan4-advanced/data")
    ap.add_argument("--faa-2024", default="CAZyDB.07142024.fungi.faa")
    ap.add_argument("--faa-2025", default="CAZyDB.07242025.fungi.faa")
    ap.add_argument("--eval-per-bucket", type=int, default=4000,
                    help="max eval sequences sampled per novelty bucket (novel_seq, novel_family)")
    ap.add_argument("--ref-per-family", type=int, default=0,
                    help="if >0, cap reference sequences per family (0 = keep all)")
    ap.add_argument("--min-len", type=int, default=30)
    ap.add_argument("--max-len", type=int, default=1500)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    random.seed(args.seed)

    import os
    os.makedirs(args.out_dir, exist_ok=True)
    p24 = os.path.join(args.data_dir, args.faa_2024)
    p25 = os.path.join(args.data_dir, args.faa_2025)

    # ---- pass 1: parse labels ----
    print("[1/4] parsing 2024 ...", file=sys.stderr)
    lab24 = {}       # pid -> (families, seq_len, md5)
    hashes24 = set() # md5 of every 2024 sequence (identity key)
    for pid, fams, slen, h in parse_faa(p24):
        if fams:
            lab24[pid] = (fams, slen, h)
            hashes24.add(h)
    print(f"      2024: {len(lab24):,} labeled proteins, {len(hashes24):,} unique sequences",
          file=sys.stderr)

    print("[2/4] parsing 2025 ...", file=sys.stderr)
    lab25 = {}
    for pid, fams, slen, h in parse_faa(p25):
        if fams:
            lab25[pid] = (fams, slen, h)
    print(f"      2025: {len(lab25):,} labeled proteins", file=sys.stderr)

    fams24 = set()
    for fams, _, _ in lab24.values():
        fams24.update(fams)
    fams25 = set()
    for fams, _, _ in lab25.values():
        fams25.update(fams)
    new_fams = sorted(fams25 - fams24)

    # ---- write reference labels (2024) ----
    ref_path = os.path.join(args.out_dir, "reference_labels_2024.tsv")
    with open(ref_path, "w") as fo:
        fo.write("protein_id\tseq_len\tseq_md5\tfamilies\tclasses\n")
        for pid, (fams, slen, h) in lab24.items():
            cls = sorted({fam_class(f) for f in fams if fam_class(f)})
            fo.write(f"{pid}\t{slen}\t{h}\t{','.join(fams)}\t{','.join(cls)}\n")

    # ---- classify 2025 novelty (exact-sequence identity vs 2024) ----
    print("[3/4] classifying 2025 novelty (sequence-hash) ...", file=sys.stderr)
    buckets = defaultdict(list)   # bucket -> [pid]
    test_path = os.path.join(args.out_dir, "test_2025_labels.tsv")
    with open(test_path, "w") as fo:
        fo.write("protein_id\tseq_len\tseq_md5\tfamilies\tclasses\tnovelty\n")
        for pid, (fams, slen, h) in lab25.items():
            if h in hashes24:
                nov = "carried_over"          # exact sequence already in 2024
            elif any(f not in fams24 for f in fams):
                nov = "novel_family"          # >=1 family absent from 2024
            else:
                nov = "novel_seq"             # new sequence, family exists in 2024
            cls = sorted({fam_class(f) for f in fams if fam_class(f)})
            fo.write(f"{pid}\t{slen}\t{h}\t{','.join(fams)}\t{','.join(cls)}\t{nov}\n")
            buckets[nov].append(pid)

    summary = {
        "novelty_definition": "exact-sequence-md5 vs 2024; carried_over = seq present in 2024",
        "n_2024_labeled": len(lab24),
        "n_2024_unique_seqs": len(hashes24),
        "n_2025_labeled": len(lab25),
        "n_families_2024": len(fams24),
        "n_families_2025": len(fams25),
        "n_new_families_2025": len(new_fams),
        "new_families_all": new_fams,
        "buckets": {k: len(v) for k, v in buckets.items()},
    }
    with open(os.path.join(args.out_dir, "novelty_summary.json"), "w") as fo:
        json.dump(summary, fo, indent=2)
    print("      " + json.dumps(summary["buckets"]), file=sys.stderr)
    print(f"      new families in 2025: {len(new_fams)}", file=sys.stderr)

    # ---- build bounded, stratified eval set from novel_seq + novel_family ----
    print("[4/4] building eval set ...", file=sys.stderr)
    eval_pids = []
    for bucket in ("novel_seq", "novel_family"):
        pids = [p for p in buckets[bucket]
                if args.min_len <= lab25[p][1] <= args.max_len]
        random.shuffle(pids)
        eval_pids.extend(pids[:args.eval_per_bucket])
    eval_set = set(eval_pids)
    seqs25 = load_seqs(p25, eval_set)
    with open(os.path.join(args.out_dir, "eval_2025.faa"), "w") as fo, \
         open(os.path.join(args.out_dir, "eval_2025_labels.tsv"), "w") as fl:
        fl.write("protein_id\tseq_len\tfamilies\tclasses\tnovelty\n")
        for pid in eval_pids:
            if pid not in seqs25:
                continue
            fams, slen, _ = lab25[pid]
            nov = "novel_family" if any(f not in fams24 for f in fams) else "novel_seq"
            cls = sorted({fam_class(f) for f in fams if fam_class(f)})
            fo.write(f">{pid}|{','.join(fams)}\n{seqs25[pid]}\n")
            fl.write(f"{pid}\t{slen}\t{','.join(fams)}\t{','.join(cls)}\t{nov}\n")

    # ---- reference FASTA (all or capped per family) ----
    ref_wanted = set(lab24)
    if args.ref_per_family > 0:
        by_fam = defaultdict(list)
        for pid, (fams, slen, _) in lab24.items():
            if args.min_len <= slen <= args.max_len:
                for f in fams:
                    by_fam[f].append(pid)
        ref_wanted = set()
        for f, pids in by_fam.items():
            random.shuffle(pids)
            ref_wanted.update(pids[:args.ref_per_family])
    seqs24 = load_seqs(p24, ref_wanted)
    with open(os.path.join(args.out_dir, "reference_2024.faa"), "w") as fo:
        for pid in ref_wanted:
            if pid in seqs24:
                fams, _, _ = lab24[pid]
                fo.write(f">{pid}|{','.join(fams)}\n{seqs24[pid]}\n")

    print(f"DONE. eval={len(eval_pids)} ref_fasta={len(ref_wanted)} out={args.out_dir}", file=sys.stderr)

if __name__ == "__main__":
    sys.exit(main())
