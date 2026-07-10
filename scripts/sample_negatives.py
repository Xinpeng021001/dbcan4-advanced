#!/usr/bin/env python3
"""
Build the shared dbCAN4 NEGATIVE / DECOY set from Track A tiered output.

Three products:
  (1) natural non-CAZyme negatives  -- stratified by taxonomic class to match the
      high_confidence_cazyme (positive) class distribution, ~100k, for TRAINING.
      Drawn from tier==high_confidence_non_cazyme, 30..1500 aa, per-class reservoir.
  (2) shuffled-domain decoys        -- residue-permuted known CAZymes (preserve length
      + AA composition, destroy fold/order), ~N_SHUFFLED, hard negatives.
  (3) realistic-imbalance slice     -- ALL proteins of a handful of whole genomes
      (>90% non-CAZyme), DISJOINT genomes from the training negatives, with a
      per-protein truth column (cazyme/non_cazyme/gray) from the Track A tier.

Deterministic (fixed seeds). Sequences pulled from each genome's uniInput.faa
(headers exactly match tiered_proteins protein_id; hyphen form).
"""
import gzip, os, sys, glob, random, argparse
from collections import defaultdict

TIERED    = "/array1/xinpeng/dbcan4-advanced/track_a_output/tiered_proteins.tsv.gz"
ALLGENOME = "/array1/xinpeng/all_genome"
REF2024   = "/array1/xinpeng/dbcan4-advanced/data/reference_2024.faa"

TRUTH_MAP = {
    "high_confidence_cazyme": "cazyme",
    "high_confidence_non_cazyme": "non_cazyme",
    "gray_zone": "gray",
}

def uniinput_path(cl, g):
    p = os.path.join(ALLGENOME, cl, g, "uniInput.faa")
    if os.path.exists(p):
        return p
    hits = glob.glob(os.path.join(ALLGENOME, "*", g, "uniInput.faa"))
    return hits[0] if hits else None

def iter_fasta(path):
    hid, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if hid is not None:
                    yield hid, "".join(buf)
                hid = line[1:].split()[0].rstrip("\n")   # first token
                buf = []
            else:
                buf.append(line.strip())
        if hid is not None:
            yield hid, "".join(buf)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=".")
    ap.add_argument("--n_neg", type=int, default=100000)
    ap.add_argument("--n_shuffled", type=int, default=15000)
    ap.add_argument("--n_slice_genomes", type=int, default=5)
    ap.add_argument("--min_len", type=int, default=30)
    ap.add_argument("--max_len", type=int, default=1500)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    random.seed(args.seed)
    os.makedirs(args.out, exist_ok=True)

    # ---------- PASS 1: per-genome tier counts + per-class positive counts ----------
    genome_tier  = defaultdict(lambda: defaultdict(int))
    genome_class = {}
    class_pos    = defaultdict(int)   # high_confidence_cazyme per class
    with gzip.open(TIERED, "rt") as f:
        header = f.readline().rstrip("\n").rstrip("\r").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        ti, pi, gi, ci = idx["tier"], idx["protein_id"], idx["genome_id"], idx["tax_class"]
        for line in f:
            p = line.rstrip("\n").rstrip("\r").split("\t")
            if len(p) <= ti:
                continue
            tier, g, cl = p[ti], p[gi], p[ci]
            genome_tier[g][tier] += 1
            genome_class[g] = cl
            if tier == "high_confidence_cazyme":
                class_pos[cl] += 1
    total_pos = sum(class_pos.values())
    print(f"[neg] {len(genome_tier)} genomes; total high_confidence_cazyme positives={total_pos}",
          file=sys.stderr, flush=True)

    # ---------- pick realistic-imbalance slice genomes (disjoint, >90% non-CAZyme) ----------
    cands = []
    for g, tc in genome_tier.items():
        tot = sum(tc.values())
        ncaz = tc.get("high_confidence_non_cazyme", 0)
        caz  = tc.get("high_confidence_cazyme", 0)
        if tot == 0:
            continue
        frac_non = ncaz / tot
        if 3000 <= tot <= 9000 and frac_non > 0.90 and caz >= 20:
            cands.append((genome_class[g], g, tot, frac_non, caz))
    cands.sort(key=lambda x: (x[0], x[1]))     # deterministic; prefer class diversity
    slice_genomes, seen_class = [], set()
    for cl, g, tot, fn, caz in cands:
        if cl in seen_class:
            continue
        slice_genomes.append(g); seen_class.add(cl)
        if len(slice_genomes) >= args.n_slice_genomes:
            break
    slice_set = set(slice_genomes)
    print(f"[neg] realistic-slice genomes ({len(slice_set)}): "
          + ", ".join(f"{g}[{genome_class[g]}]" for g in slice_genomes), file=sys.stderr, flush=True)

    # ---------- per-class negative targets proportional to positive distribution ----------
    targets = {}
    for cl, npos in class_pos.items():
        targets[cl] = int(round(args.n_neg * npos / total_pos))
    print(f"[neg] {len([t for t in targets.values() if t>0])} classes with a target; "
          f"sum={sum(targets.values())}", file=sys.stderr, flush=True)

    # ---------- PASS 2: reservoir-sample negatives per class + collect slice rows ----------
    reservoirs = {cl: [] for cl in targets}
    seen_cnt   = {cl: 0  for cl in targets}
    slice_rows = []   # (pid, genome, class, len, tier)
    with gzip.open(TIERED, "rt") as f:
        header = f.readline().rstrip("\n").rstrip("\r").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        ti, pi, gi, ci, li = idx["tier"], idx["protein_id"], idx["genome_id"], idx["tax_class"], idx["seq_length"]
        for line in f:
            p = line.rstrip("\n").rstrip("\r").split("\t")
            if len(p) <= li:
                continue
            tier, pid, g, cl = p[ti], p[pi], p[gi], p[ci]
            try:
                L = int(p[li])
            except ValueError:
                continue
            if g in slice_set:
                slice_rows.append((pid, g, cl, L, tier))
                continue
            if tier != "high_confidence_non_cazyme":
                continue
            if not (args.min_len <= L <= args.max_len):
                continue
            k = targets.get(cl, 0)
            if k <= 0:
                continue
            seen_cnt[cl] += 1
            res = reservoirs[cl]
            if len(res) < k:
                res.append((pid, g, cl, L))
            else:
                j = random.randint(0, seen_cnt[cl] - 1)
                if j < k:
                    res[j] = (pid, g, cl, L)

    negatives = [row for cl in reservoirs for row in reservoirs[cl]]
    print(f"[neg] sampled {len(negatives)} natural negatives; slice rows={len(slice_rows)}",
          file=sys.stderr, flush=True)

    # ---------- extract sequences from uniInput.faa (group by genome) ----------
    need = defaultdict(dict)   # genome -> {pid: class}
    for pid, g, cl, L in negatives:
        need[g][pid] = cl
    for pid, g, cl, L, tier in slice_rows:
        need[g][pid] = cl
    seqs = {}
    missing_genomes = 0
    for gi_n, (g, ids) in enumerate(need.items()):
        cl = genome_class[g]
        path = uniinput_path(cl, g)
        if path is None:
            missing_genomes += 1
            continue
        want = set(ids)
        for hid, seq in iter_fasta(path):
            if hid in want:
                seqs[hid] = seq
        if (gi_n + 1) % 200 == 0:
            print(f"[neg] extracted sequences from {gi_n+1}/{len(need)} genomes", file=sys.stderr, flush=True)
    print(f"[neg] sequences found for {len(seqs)} proteins; missing genome dirs={missing_genomes}",
          file=sys.stderr, flush=True)

    # ---------- shuffled-domain decoys from reference_2024 CAZymes ----------
    shuf_reservoir = []   # (pid, fam, seq)
    n_seen = 0
    for hid, seq in iter_fasta(REF2024):
        L = len(seq)
        if not (args.min_len <= L <= args.max_len):
            continue
        toks = hid.split("|")
        opid = toks[0]
        ofam = toks[1] if len(toks) > 1 else "-"
        n_seen += 1
        if len(shuf_reservoir) < args.n_shuffled:
            shuf_reservoir.append((opid, ofam, seq))
        else:
            j = random.randint(0, n_seen - 1)
            if j < args.n_shuffled:
                shuf_reservoir[j] = (opid, ofam, seq)
    shuffled = []
    for opid, ofam, seq in shuf_reservoir:
        chars = list(seq)
        random.shuffle(chars)                       # permutation -> preserve length + AA comp
        shuffled.append((f"shuffled_{opid}", ofam, "".join(chars), len(seq)))
    print(f"[neg] built {len(shuffled)} shuffled-domain decoys from {n_seen} candidate CAZymes",
          file=sys.stderr, flush=True)

    # ---------- write decoy_set.tsv / .faa ----------
    n_natural_written = 0
    with open(os.path.join(args.out, "decoy_set.tsv"), "w") as t, \
         open(os.path.join(args.out, "decoy_set.faa"), "w") as fa:
        t.write("protein_id\tdecoy_type\tsource_tier\ttax_class\tgenome_id\tseq_len\n")
        for pid, g, cl, L in negatives:
            if pid not in seqs:
                continue
            t.write(f"{pid}\tnatural_non_cazyme\thigh_confidence_non_cazyme\t{cl}\t{g}\t{L}\n")
            fa.write(f">{pid}\n{seqs[pid]}\n")
            n_natural_written += 1
        for dpid, ofam, dseq, L in shuffled:
            # source_tier holds the CAZy family the decoy was scrambled from (provenance)
            t.write(f"{dpid}\tshuffled_domain\t{ofam}\t-\t-\t{L}\n")
            fa.write(f">{dpid}\n{dseq}\n")

    # ---------- write realistic_imbalance_slice.tsv / .faa ----------
    n_slice_written = 0
    slice_truth_counts = defaultdict(int)
    with open(os.path.join(args.out, "realistic_imbalance_slice.tsv"), "w") as t, \
         open(os.path.join(args.out, "realistic_imbalance_slice.faa"), "w") as fa:
        t.write("protein_id\ttruth\ttier\ttax_class\tgenome_id\tseq_len\n")
        for pid, g, cl, L, tier in slice_rows:
            if pid not in seqs:
                continue
            truth = TRUTH_MAP.get(tier, "unknown")
            t.write(f"{pid}\t{truth}\t{tier}\t{cl}\t{g}\t{L}\n")
            fa.write(f">{pid}\n{seqs[pid]}\n")
            n_slice_written += 1
            slice_truth_counts[truth] += 1

    # ---------- summary ----------
    import json
    summary = {
        "n_natural_negatives": n_natural_written,
        "n_shuffled_decoys": len(shuffled),
        "n_realistic_imbalance_slice": n_slice_written,
        "slice_genomes": slice_genomes,
        "slice_genome_classes": {g: genome_class[g] for g in slice_genomes},
        "slice_truth_counts": dict(slice_truth_counts),
        "realistic_slice_cazyme_frac": round(slice_truth_counts.get("cazyme", 0) / n_slice_written, 6) if n_slice_written else 0.0,
        "n_classes_sampled": len([t for t in targets.values() if t > 0]),
        "seed": args.seed,
    }
    with open(os.path.join(args.out, "decoy_sampling_summary.json"), "w") as j:
        json.dump(summary, j, indent=2)
    print("[neg] SUMMARY " + json.dumps(summary), file=sys.stderr, flush=True)

if __name__ == "__main__":
    main()
