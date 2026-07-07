#!/usr/bin/env python3
"""
Parse DIAMOND BLASTp output for the sequence-similarity baseline.

DIAMOND is run as (in the job wrapper):
    diamond makedb --in CAZyDB.07142024.fa -d cazy2024
    diamond blastp -q eval_2025.faa -d cazy2024 -o hits.tsv \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --max-target-seqs 5 --evalue 1e-15 -p N

Both query and subject headers are "ID|FAM[|FAM...]" (no spaces), so qseqid and
sseqid carry the family label directly — the same FAM_RE used in build_reference.py
extracts it. Best hit (top bitscore) per query defines the predicted family set.

Outputs: predictions TSV (query_id, true_families, pred_families, top_pident,
top_evalue, correct) + a JSON summary of recall/precision at several e-value cuts.
"""
import argparse, json, re, sys
from collections import defaultdict

FAM_RE = re.compile(r'^(GH|GT|PL|CE|AA|CBM)\d+(_\d+)?$')

def fams_from_token(tok):
    # subject headers vary: "ID|FAM|FAM" (full CAZy) or "ID|FAM,FAM" (fungal reference)
    return [f for f in re.split(r'[|,]', tok) if FAM_RE.match(f)]

def load_truth(labels_tsv):
    truth, novelty = {}, {}
    with open(labels_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        fi = header.index("families"); ni = header.index("novelty"); pi = header.index("protein_id")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            truth[p[pi]] = set(x for x in p[fi].split(",") if x)
            novelty[p[pi]] = p[ni]
    return truth, novelty

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out-pred", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--evalues", default="1e-15,1e-10,1e-5,1e-3")
    args = ap.parse_args()

    truth, novelty = load_truth(args.labels)

    # collect all hits per query: (evalue, bitscore, pident, subject_families)
    hits = defaultdict(list)
    with open(args.hits) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            q, s, pid, length, ev, bit = f[0], f[1], f[2], f[3], f[4], f[5]
            qid = q.split("|")[0]
            sfams = fams_from_token(s)
            hits[qid].append((float(ev), float(bit), float(pid), sfams))

    ev_cuts = [float(x) for x in args.evalues.split(",")]

    # best-hit prediction at the loosest cut for the per-query table
    loosest = max(ev_cuts)
    rows = []
    for qid, tfams in truth.items():
        qhits = [h for h in hits.get(qid, []) if h[0] <= loosest]
        if qhits:
            qhits.sort(key=lambda x: (-x[1], x[0]))  # top bitscore
            best = qhits[0]
            pred = set(best[3])
            rows.append((qid, tfams, pred, best[2], best[0], novelty[qid]))
        else:
            rows.append((qid, tfams, set(), 0.0, None, novelty[qid]))

    with open(args.out_pred, "w") as fo:
        fo.write("query_id\tnovelty\ttrue_families\tpred_families\ttop_pident\ttop_evalue\texact\toverlap\n")
        for qid, tf, pf, pident, ev, nov in rows:
            exact = int(tf == pf and len(pf) > 0)
            overlap = int(len(tf & pf) > 0)
            fo.write(f"{qid}\t{nov}\t{','.join(sorted(tf))}\t{','.join(sorted(pf))}\t"
                     f"{pident:.1f}\t{'' if ev is None else ev}\t{exact}\t{overlap}\n")

    # metrics at each e-value cut, overall and per novelty bucket
    def metrics(subset_qids):
        out = {}
        for cut in ev_cuts:
            n = len(subset_qids); hit=0; exact=0; overlap=0
            for qid in subset_qids:
                tf = truth[qid]
                qhits = [h for h in hits.get(qid, []) if h[0] <= cut]
                if not qhits:
                    continue
                qhits.sort(key=lambda x: (-x[1], x[0]))
                pf = set(qhits[0][3])
                if pf: hit += 1
                if pf == tf and pf: exact += 1
                if tf & pf: overlap += 1
            out[f"evalue<={cut:g}"] = {
                "n": n,
                "any_hit": hit, "any_hit_frac": round(hit/n, 4) if n else 0,
                "exact_family": exact, "exact_frac": round(exact/n, 4) if n else 0,
                "overlap_family": overlap, "overlap_frac": round(overlap/n, 4) if n else 0,
            }
        return out

    all_q = list(truth)
    by_nov = defaultdict(list)
    for qid in all_q:
        by_nov[novelty[qid]].append(qid)

    summary = {
        "n_eval": len(all_q),
        "overall": metrics(all_q),
        "by_novelty": {k: metrics(v) for k, v in by_nov.items()},
    }
    with open(args.out_summary, "w") as fo:
        json.dump(summary, fo, indent=2)
    print(json.dumps({"overall": summary["overall"],
                      "by_novelty_n": {k: len(v) for k, v in by_nov.items()}}, indent=2))

if __name__ == "__main__":
    sys.exit(main())
