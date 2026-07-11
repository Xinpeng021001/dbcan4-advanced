#!/usr/bin/env python3
"""Convert raw tool output into the dbCAN4-advanced v1.1 feature-contract TSVs.

Subcommands, one per contract feature type (§2.5-2.9):
  domains        hmmscan --domtblout  -> features/<s>/domains.tsv        (§2.5)
  structure_hits foldseek aln tsv     -> features/<s>/structure_hits.tsv (§2.6)
  localization   SP + topology        -> features/<s>/localization.tsv   (§2.7, derived)
  physicochem    protein FASTA        -> features/<s>/physicochem.tsv     (§2.8, Biopython)
  ec_prediction  CLEAN output         -> features/<s>/ec_prediction.tsv   (§2.9)

Each writes the exact columns bioforge.ingest.parse_advanced expects.
"""
from __future__ import annotations
import argparse, csv, json, os, re, sys


# ---------- §2.5 domains (hmmscan --domtblout) ----------
def domains(a):
    """Parse hmmscan --domtblout into one row per domain occurrence, N->C."""
    rows = []
    with open(a.domtbl) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            # hmmscan domtblout columns (Pfam): 0 tname 1 tacc ... 3 qname ...
            # 6 full-eval 7 full-score ... 11 i-Evalue 13 domain-score
            # 15 hmm_from 16 hmm_to 17 ali_from 18 ali_to ... 21 hmm_len? use tlen(2)
            tname, tacc = f[0], f[1]
            qname = f[3]
            tlen = int(f[2])
            i_eval = float(f[12]); dom_score = float(f[13])
            hmm_from, hmm_to = int(f[15]), int(f[16])
            ali_from, ali_to = int(f[17]), int(f[18])
            hmm_cov = (hmm_to - hmm_from + 1) / tlen if tlen else 0.0
            acc = tacc.split(".")[0] if tacc and tacc != "-" else tname
            rows.append((qname, acc, tname, ali_from, ali_to, i_eval, dom_score, hmm_cov))
    rows.sort(key=lambda r: (r[0], r[3]))  # per protein, N->C
    with open(a.out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["protein_id", "acc", "name", "start", "end", "evalue", "score", "extra"])
        for qname, acc, name, s, e, ev, sc, cov in rows:
            w.writerow([qname, acc, name, s, e, f"{ev:.2e}", f"{sc:.1f}",
                        json.dumps({"hmm_coverage": round(cov, 4)})])
    print(f"[domains] {a.out}: {len(rows)} domain rows")


# ---------- §2.6 structure_hits (foldseek aln tsv) ----------
def structure_hits(a):
    """foldseek easy-search tsv -> top-K structural homologs per protein."""
    # foldseek fmt: query,target,fident,alnlen,...,evalue,bits,...  we read a
    # normalized tsv with header if present; else assume standard 'query target
    # fident alnlen mismatch gapopen qstart qend tstart tend evalue bits'
    by_q = {}
    with open(a.aln) as fh:
        first = fh.readline().rstrip("\n").split("\t")
        has_hdr = "query" in [c.lower() for c in first]
        def rowget(cols, names, idx):
            if has_hdr:
                h = {c.lower(): i for i, c in enumerate(first)}
                for n in names:
                    if n in h and h[n] < len(cols):
                        return cols[h[n]]
                return None
            return cols[idx] if idx < len(cols) else None
        lines = fh if has_hdr else [("\t".join(first))] + list(fh)
        for line in lines:
            cols = line.rstrip("\n").split("\t") if isinstance(line, str) else line
            if not cols or not cols[0]:
                continue
            q = rowget(cols, ["query"], 0)
            t = rowget(cols, ["target"], 1)
            fident = float(rowget(cols, ["fident"], 2) or 0)
            evalue = float(rowget(cols, ["evalue"], 10) or 0)
            bits = float(rowget(cols, ["bits"], 11) or 0)
            tm = rowget(cols, ["alntmscore", "tmscore"], -1)
            lddt = rowget(cols, ["lddt"], -1)
            # target_family: parse a trailing _FAM tag or leave '-'
            m = re.search(r"([A-Z]{2,3}\d+)", str(t) or "")
            fam = m.group(1) if m else "-"
            by_q.setdefault(q, []).append(
                {"target": t, "target_family": fam,
                 "tmscore": float(tm) if tm not in (None, "") else "",
                 "lddt": float(lddt) if lddt not in (None, "") else "",
                 "evalue": evalue, "bits": bits, "fident": fident})
    with open(a.out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["protein_id", "target", "target_family", "tmscore", "prob", "lddt", "evalue", "extra"])
        n = 0
        for q, hits in by_q.items():
            hits.sort(key=lambda h: -h["bits"])
            for rank, h in enumerate(hits[:a.topk], 1):
                w.writerow([q, h["target"], h["target_family"],
                            (f"{h['tmscore']:.4f}" if h["tmscore"] != "" else ""), "",
                            (f"{h['lddt']:.4f}" if h["lddt"] != "" else ""),
                            f"{h['evalue']:.2e}",
                            json.dumps({"bits": h["bits"], "fident": h["fident"], "rank": rank})])
                n += 1
    print(f"[structure_hits] {a.out}: {n} hit rows")


# ---------- §2.7 localization (derived from SP) ----------
def localization(a):
    """Derive localization from SignalP call: SP present -> Extracellular (secreted)."""
    # Read all SP rows and write one localization row per protein.
    out_rows = []
    if a.signalp and os.path.exists(a.signalp):
        with open(a.signalp) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                pid = row.get("protein_id")
                pred = (row.get("prediction") or "").upper()
                prob = row.get("sp_prob") or ""
                is_sp = ("SP" in pred) or (pred in ("SIGNAL", "SEC/SPI"))
                lz = "Extracellular" if is_sp else "Intracellular"
                out_rows.append((pid, lz, prob, "derived-from-SP",
                                 json.dumps({"basis": "signalp", "sp_prediction": pred})))
    with open(a.out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["protein_id", "localization", "confidence", "method", "extra"])
        w.writerows(out_rows)
    print(f"[localization] {a.out}: {len(out_rows)} rows")


# ---------- §2.8 physicochem (Biopython) ----------
def physicochem(a):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    rows = []
    pid, seq = None, []
    def flush(pid, seq):
        s = "".join(seq).upper().replace("*", "").replace("X", "")
        if not pid or not s:
            return None
        pa = ProteinAnalysis(s)
        try:
            mw = pa.molecular_weight()
        except Exception:
            return None
        n_glyc = len(re.findall(r"N[^P][ST]", "".join(seq).upper()))
        return (pid, mw, pa.isoelectric_point(), pa.instability_index(),
                pa.gravy(), pa.aromaticity(), n_glyc)
    for line in open(a.faa):
        if line.startswith(">"):
            r = flush(pid, seq)
            if r:
                rows.append(r)
            pid = line[1:].strip().split()[0].split("|")[0]; seq = []
        else:
            seq.append(line.strip())
    r = flush(pid, seq)
    if r:
        rows.append(r)
    with open(a.out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["protein_id", "mw_da", "pi", "instability", "gravy", "aromaticity", "extra"])
        for pid, mw, pi, inst, grav, arom, nglyc in rows:
            w.writerow([pid, f"{mw:.2f}", f"{pi:.2f}", f"{inst:.2f}", f"{grav:.4f}",
                        f"{arom:.4f}", json.dumps({"n_glyc_sequons": nglyc})])
    print(f"[physicochem] {a.out}: {len(rows)} rows")


# ---------- §2.9 ec_prediction (CLEAN) ----------
def ec_prediction(a):
    """Parse CLEAN maxsep output (protein_id, EC:conf, EC:conf, ...) -> top-K rows."""
    rows = []
    with open(a.clean) as fh:
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 2:
                continue
            pid = parts[0]
            rank = 0
            for p in parts[1:]:
                m = re.search(r"EC:([\d.\-]+)/?([\d.]+)?", p) or re.match(r"([\d.]+):([\d.]+)", p.strip())
                if not m:
                    # format "EC:3.2.1.40/0.9953"
                    mm = re.match(r"EC:([\d.\-]+)/([\d.]+)", p.strip())
                    if mm:
                        ec, conf = mm.group(1), float(mm.group(2))
                    else:
                        continue
                else:
                    ec = m.group(1); conf = float(m.group(2)) if m.group(2) else 0.0
                rank += 1
                rows.append((pid, ec, conf, rank))
    with open(a.out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["protein_id", "ec_number", "confidence", "rank", "tool", "extra"])
        for pid, ec, conf, rank in rows:
            w.writerow([pid, ec, f"{conf:.4f}", rank, "CLEAN",
                        json.dumps({"confidence_type": "maxsep"})])
    print(f"[ec_prediction] {a.out}: {len(rows)} rows")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("domains"); p.add_argument("--domtbl", required=True); p.add_argument("--out", required=True); p.set_defaults(fn=domains)
    p = sub.add_parser("structure_hits"); p.add_argument("--aln", required=True); p.add_argument("--out", required=True); p.add_argument("--topk", type=int, default=10); p.set_defaults(fn=structure_hits)
    p = sub.add_parser("localization"); p.add_argument("--signalp", required=True); p.add_argument("--out", required=True); p.set_defaults(fn=localization)
    p = sub.add_parser("physicochem"); p.add_argument("--faa", required=True); p.add_argument("--out", required=True); p.set_defaults(fn=physicochem)
    p = sub.add_parser("ec_prediction"); p.add_argument("--clean", required=True); p.add_argument("--out", required=True); p.set_defaults(fn=ec_prediction)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
