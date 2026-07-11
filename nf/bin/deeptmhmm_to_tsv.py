#!/usr/bin/env python3
"""Convert DeepTMHMM output into the v1.1 contract feature TSVs.

DeepTMHMM (run via `biolib run DTU/DeepTMHMM`) writes biolib_results/TMRs.gff3
(region spans + "Number of predicted TMRs") and predicted_topologies.3line
(per-residue S/I/O/M/B labels + a per-protein category header ">id | SP").

This emits TWO contract TSVs from that one real run:
  --out-tm      deeptmhmm.tsv  (feature_type tm_topology): protein_id, prediction, n_tm, topology, extra
  --out-sp      signalp6.tsv   (feature_type signal_peptide): protein_id, prediction, sp_prob, cs_position, extra
                 -- honestly sourced from DeepTMHMM's SP call when SignalP-6.0 is not installed;
                    sp_prob is left blank (DeepTMHMM emits no calibrated SP probability -> never fabricated).
"""
from __future__ import annotations
import argparse, csv, json, os, re

def parse_gff3(path):
    """Return {pid: {'length':int,'n_tm':int,'regions':[(type,start,end)]}}."""
    prot = {}
    cur = None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("#"):
            m = re.match(r"#\s+(\S+)\s+Length:\s+(\d+)", line)
            if m:
                cur = m.group(1)
                prot.setdefault(cur, {"length": int(m.group(2)), "n_tm": 0, "regions": []})
                continue
            m = re.match(r"#\s+(\S+)\s+Number of predicted TMRs:\s+(\d+)", line)
            if m:
                prot.setdefault(m.group(1), {"length": 0, "n_tm": 0, "regions": []})
                prot[m.group(1)]["n_tm"] = int(m.group(2))
                continue
        elif line and not line.startswith("//"):
            f = line.split("\t")
            if len(f) >= 4:
                pid, rtype, start, end = f[0], f[1], f[2], f[3]
                prot.setdefault(pid, {"length": 0, "n_tm": 0, "regions": []})
                try:
                    prot[pid]["regions"].append((rtype, int(start), int(end)))
                except ValueError:
                    pass
    return prot

def parse_3line(path):
    """Return {pid: {'category':'SP'|'GLOB'|...,'topology':'SSOO...'}}."""
    out = {}
    if not os.path.exists(path):
        return out
    lines = [l.rstrip("\n") for l in open(path)]
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            hdr = lines[i][1:].strip()
            parts = [p.strip() for p in hdr.split("|")]
            pid = parts[0].split()[0]
            cat = parts[1] if len(parts) > 1 else ""
            seq = lines[i+1] if i+1 < len(lines) else ""
            topo = lines[i+2] if i+2 < len(lines) else ""
            out[pid] = {"category": cat, "topology": topo}
            i += 3
        else:
            i += 1
    return out

CAT_MAP = {"SP": "SP", "SP+TM": "SP+TM", "TM": "TM", "GLOB": "Globular",
           "BETA": "TM_beta", "SP+BETA": "SP+TM_beta"}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff3", required=True, help="biolib_results/TMRs.gff3")
    ap.add_argument("--three-line", default=None, help="biolib_results/predicted_topologies.3line")
    ap.add_argument("--out-tm", required=True, help="deeptmhmm.tsv (tm_topology)")
    ap.add_argument("--out-sp", default=None, help="signalp6.tsv (signal_peptide, derived from DeepTMHMM)")
    a = ap.parse_args()

    gff = parse_gff3(a.gff3)
    tl = parse_3line(a.three_line) if a.three_line else {}

    tm_rows, sp_rows = [], []
    for pid, d in gff.items():
        regions = d["regions"]
        sp_region = next((r for r in regions if r[0] == "signal"), None)
        has_sp = sp_region is not None
        cat_raw = tl.get(pid, {}).get("category", "")
        prediction = CAT_MAP.get(cat_raw, cat_raw or ("SP" if has_sp else "Globular"))
        # Compact, lossless region encoding (e.g. "signal:1-18|outside:19-1089")
        # instead of the raw 1-char-per-residue string, which is unreadable and
        # spills over many display pages for long proteins.
        topo = "|".join(f"{rtype}:{s}-{e}" for rtype, s, e in regions) if regions else ""
        full_topo = tl.get(pid, {}).get("topology", "")   # keep raw string in extra
        tm_rows.append((pid, prediction, d["n_tm"], topo,
                        json.dumps({"tool": "DeepTMHMM", "has_signal_peptide": has_sp,
                                    "residue_topology": full_topo})))
        # signal-peptide TSV derived from DeepTMHMM SP region
        sp_pred = "SP" if has_sp else "NO_SP"
        cs = sp_region[2] if has_sp else ""   # cleavage after last SP residue
        sp_rows.append((pid, sp_pred, "", cs,
                        json.dumps({"source": "DeepTMHMM (SignalP-6.0 not installed)",
                                    "sp_span": [sp_region[1], sp_region[2]] if has_sp else None})))

    with open(a.out_tm, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["protein_id", "prediction", "n_tm", "topology", "extra"])
        w.writerows(tm_rows)
    print(f"[deeptmhmm_to_tsv] {a.out_tm}: {len(tm_rows)} rows")

    if a.out_sp:
        with open(a.out_sp, "w", newline="") as fo:
            w = csv.writer(fo, delimiter="\t")
            w.writerow(["protein_id", "prediction", "sp_prob", "cs_position", "extra"])
            w.writerows(sp_rows)
        print(f"[deeptmhmm_to_tsv] {a.out_sp}: {len(sp_rows)} rows (signal_peptide derived from DeepTMHMM)")

if __name__ == "__main__":
    main()
