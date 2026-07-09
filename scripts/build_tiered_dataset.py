#!/usr/bin/env python3
"""
build_tiered_dataset.py -- dbCAN4 Track A: multi-evidence CAZyme/non-CAZyme tiering.

Walks /array1/xinpeng/all_genome/<TaxClass>/<GenomeID>/ (2,226 genome directories),
aggregates the per-genome dbCAN annotation outputs (overview.tsv, diamond.out,
dbCAN_hmm_results.tsv, dbCANsub_hmm_results.tsv, non_CAZyme.faa, uniInput.faa) into
one streamed per-protein TSV, computes a multi-evidence confidence tier for every
protein in every proteome, and (memory-boundedly, via per-stratum reservoir
sampling) builds a stratified representative sample for downstream structure-based
validation (Track B).

Design notes / non-obvious facts this script relies on (confirmed by direct
inspection on met, 2026-07-09):
  - overview.tsv only lists proteins that got >=1 hit from >=1 tool. Proteins with
    ZERO hits from all three tools (diamond, dbCAN_hmm, dbCAN_sub) never appear in
    overview.tsv at all -- they are the "clean negative" candidates and are only
    enumerated via uniInput.faa (full proteome) minus overview.tsv's Gene IDs.
  - The raw per-tool files (diamond.out, dbCAN_hmm_results.tsv,
    dbCANsub_hmm_results.tsv) are PRE-FILTERED to each tool's own significance
    threshold already (i.e. every row in them is a "real" hit by that tool's
    standard) -- a single-tool hit is genuine homology evidence, it just fails
    dbCAN's own >=2-tool consensus rule that decides the "Recommend Results"
    column. This is exactly the "gray zone" this project cares about.
  - Gene ID separator differs by file: uniInput.faa / non_CAZyme.faa headers use
    "|" (e.g. "jgi|Abobi1|105364|CE105363_8952"); overview.tsv / dbCAN_hmm_results
    .tsv / dbCANsub_hmm_results.tsv use "-" (e.g. "jgi-Abobi1-105364-CE105363_8952").
    diamond.out's "Gene ID" column uses "|" like the FASTA headers. All are
    normalized to the dash form here (`canon_id`).
  - diamond.out's "CAZy ID" field is "<CAZy_seq_id>|<source_fasta>|<Family>" --
    family is the substring after the LAST "|".
  - A Gene ID can appear multiple times in dbCAN_hmm_results.tsv /
    dbCANsub_hmm_results.tsv / diamond.out (multiple HMM/domain hits); we keep the
    single best-scoring row per protein per tool (lowest i-Evalue for HMM-based
    tools, lowest E-value / highest bit score for DIAMOND).

Usage (on met, inside the project venv which already has pandas):
    /array1/xinpeng/dbcan4-advanced/venv/bin/python scripts/build_tiered_dataset.py \
        --data-root /array1/xinpeng/all_genome \
        --out-dir /array1/xinpeng/dbcan4-advanced/track_a_output \
        --sample-size 4000 \
        --reservoir-per-stratum 300

Outputs (written to --out-dir):
    tiered_proteins.tsv.gz      full per-protein tiered table, one row per protein
                                 in uniInput.faa across all genomes (streamed, gzip)
    summary_tier_by_class.tsv   counts per (tax_class, tier, subtier)
    summary_gray_zone_families.tsv  family distribution restricted to gray-zone rows
    summary_overall.json        headline totals + run metadata
    sample_for_structure.fasta  stratified representative sample (protein sequences)
    sample_for_structure.tsv    matching per-protein evidence/tier table for the sample
    failed_genomes.tsv          any genome dirs that errored out during parsing (id + reason)
"""
import argparse
import csv
import gzip
import json
import math
import os
import random
import sys
import time
from collections import defaultdict, Counter

FIELDS = [
    "protein_id", "genome_id", "tax_class", "seq_length",
    "ec_numbers",
    "dbcan_hmm_family", "dbcan_hmm_ievalue", "dbcan_hmm_coverage",
    "dbcan_sub_family", "dbcan_sub_ievalue", "dbcan_sub_coverage", "dbcan_sub_substrate",
    "diamond_family", "diamond_pident", "diamond_evalue", "diamond_bitscore", "diamond_cazy_hit_id",
    "n_tools", "recommend_family",
    "hmm_component", "sub_component", "diamond_component", "evidence_score",
    "tier", "subtier",
    "structure_evidence_score",
]


def canon_id(raw_id):
    """Normalize a Gene ID to the dash-separated form used by overview.tsv/HMM files."""
    return raw_id.replace("|", "-")


def parse_fasta_lengths(path):
    """Return {canon_id: seq_length} for every record in a FASTA file, streamed."""
    lengths = {}
    gid = None
    ln = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if gid is not None:
                    lengths[gid] = ln
                gid = canon_id(line[1:].strip())
                ln = 0
            else:
                ln += len(line.strip())
        if gid is not None:
            lengths[gid] = ln
    return lengths


def parse_fasta_seqs(path):
    """Return {canon_id: sequence} for every record in a FASTA file (used only for
    genomes that contain at least one reservoir-sampled protein -- not called for
    every genome, to keep memory bounded)."""
    seqs = {}
    gid = None
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if gid is not None:
                    seqs[gid] = "".join(buf)
                gid = canon_id(line[1:].strip())
                buf = []
            else:
                buf.append(line.strip())
        if gid is not None:
            seqs[gid] = "".join(buf)
    return seqs


def safe_float(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def last_field(s, sep="|"):
    return s.split(sep)[-1] if s and s != "-" else "-"


def log_component(ievalue, coverage, sat_log=20.0):
    """Map (i-Evalue, coverage) -> [0,1] component score.
    -log10(ievalue) saturates at `sat_log` (i.e. ievalue <= 1e-20 -> full strength),
    scaled by fractional coverage of the profile/alignment."""
    if ievalue is None or coverage is None:
        return 0.0
    ievalue = max(ievalue, 1e-300)
    strength = min(max(-math.log10(ievalue), 0.0), sat_log) / sat_log
    return max(0.0, min(1.0, strength)) * max(0.0, min(1.0, coverage))


def diamond_component_score(evalue, pident, sat_log=20.0):
    if evalue is None or pident is None:
        return 0.0
    evalue = max(evalue, 1e-300)
    strength = min(max(-math.log10(evalue), 0.0), sat_log) / sat_log
    ident_frac = max(0.0, min(1.0, pident / 100.0))
    return max(0.0, min(1.0, strength)) * ident_frac


def parse_overview(path):
    """Return dict: canon_id -> row dict (only proteins with >=1 tool hit)."""
    rows = {}
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            if len(row) < 7:
                continue
            gid, ec, hmm_fam, sub_fam, diamond_fam, n_tools, recommend = row[:7]
            rows[gid] = {
                "ec_numbers": ec,
                "dbcan_hmm_family": hmm_fam,
                "dbcan_sub_family": sub_fam,
                "diamond_family": diamond_fam,
                "n_tools": int(n_tools) if n_tools.strip().isdigit() else 0,
                "recommend_family": recommend,
            }
    return rows


def parse_dbcan_hmm(path):
    """Return dict: canon_id -> best (lowest i-Evalue) {ievalue, coverage} row."""
    best = {}
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            if len(row) < 10:
                continue
            target = row[2]
            ievalue = safe_float(row[4])
            coverage = safe_float(row[9])
            if ievalue is None:
                continue
            cur = best.get(target)
            if cur is None or ievalue < cur["ievalue"]:
                best[target] = {"ievalue": ievalue, "coverage": coverage}
    return best


def parse_dbcansub_hmm(path):
    """Return dict: canon_id -> best (lowest i-Evalue) {ievalue, coverage, substrate} row."""
    best = {}
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            if len(row) < 13:
                continue
            substrate = row[3]
            target = row[5]
            ievalue = safe_float(row[7])
            coverage = safe_float(row[12])
            if ievalue is None:
                continue
            cur = best.get(target)
            if cur is None or ievalue < cur["ievalue"]:
                best[target] = {"ievalue": ievalue, "coverage": coverage, "substrate": substrate}
    return best


def parse_diamond(path):
    """Return dict: canon_id -> best (lowest E-value) {pident, evalue, bitscore, cazy_id, family} row."""
    best = {}
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            if len(row) < 12:
                continue
            gid = canon_id(row[0])
            cazy_id_field = row[1]
            pident = safe_float(row[2])
            evalue = safe_float(row[10])
            bitscore = safe_float(row[11])
            if evalue is None:
                continue
            cur = best.get(gid)
            if cur is None or evalue < cur["evalue"]:
                best[gid] = {
                    "pident": pident, "evalue": evalue, "bitscore": bitscore,
                    "cazy_id": cazy_id_field.split("|")[0] if cazy_id_field != "-" else "-",
                    "family": last_field(cazy_id_field),
                }
    return best


def tier_for_row(n_tools, hmm_c, sub_c, dia_c):
    evidence_score = max(hmm_c, sub_c, dia_c)
    if n_tools >= 2:
        tier = "high_confidence_cazyme"
        subtier = "1A_all_tools_agree" if n_tools >= 3 else "1B_two_tools_agree"
    elif n_tools == 1:
        tier = "gray_zone"
        subtier = "2A_gray_high" if evidence_score >= 0.5 else "2B_gray_low"
    else:
        tier = "high_confidence_non_cazyme"
        subtier = "3_zero_hits"
    return tier, subtier, evidence_score


class ReservoirSampler:
    """Simple per-key reservoir sampler with algorithm R, capped at `cap` items/key."""

    def __init__(self, cap):
        self.cap = cap
        self.reservoirs = defaultdict(list)
        self.seen = Counter()

    def offer(self, key, item):
        self.seen[key] += 1
        res = self.reservoirs[key]
        n = self.seen[key]
        if len(res) < self.cap:
            res.append(item)
        else:
            j = random.randint(1, n)
            if j <= self.cap:
                res[j - 1] = item

    def all_items(self):
        out = []
        for items in self.reservoirs.values():
            out.extend(items)
        return out


def process_genome(genome_dir, tax_class, genome_id, writer, tier_counts, gray_family_counts,
                    reservoir, seq_cache_needed):
    overview_path = os.path.join(genome_dir, "overview.tsv")
    hmm_path = os.path.join(genome_dir, "dbCAN_hmm_results.tsv")
    sub_path = os.path.join(genome_dir, "dbCANsub_hmm_results.tsv")
    diamond_path = os.path.join(genome_dir, "diamond.out")
    uniinput_path = os.path.join(genome_dir, "uniInput.faa")

    for p in (overview_path, hmm_path, sub_path, diamond_path, uniinput_path):
        if not os.path.isfile(p):
            raise FileNotFoundError(p)

    overview = parse_overview(overview_path)
    hmm_best = parse_dbcan_hmm(hmm_path)
    sub_best = parse_dbcansub_hmm(sub_path)
    dia_best = parse_diamond(diamond_path)
    lengths = parse_fasta_lengths(uniinput_path)

    n_rows = 0

    for gid, seq_len in lengths.items():
        ov = overview.get(gid)
        hb = hmm_best.get(gid)
        sb = sub_best.get(gid)
        db = dia_best.get(gid)

        n_tools = ov["n_tools"] if ov else 0
        ec_numbers = ov["ec_numbers"] if ov else "-"
        recommend_family = ov["recommend_family"] if ov else "-"
        dbcan_hmm_family = ov["dbcan_hmm_family"] if ov else "-"
        dbcan_sub_family = ov["dbcan_sub_family"] if ov else "-"
        diamond_family = ov["diamond_family"] if ov else "-"

        hmm_ievalue = hb["ievalue"] if hb else None
        hmm_coverage = hb["coverage"] if hb else None
        sub_ievalue = sb["ievalue"] if sb else None
        sub_coverage = sb["coverage"] if sb else None
        sub_substrate = sb["substrate"] if sb else "-"
        dia_pident = db["pident"] if db else None
        dia_evalue = db["evalue"] if db else None
        dia_bitscore = db["bitscore"] if db else None
        dia_cazy_id = db["cazy_id"] if db else "-"

        hmm_c = log_component(hmm_ievalue, hmm_coverage)
        sub_c = log_component(sub_ievalue, sub_coverage)
        dia_c = diamond_component_score(dia_evalue, dia_pident)

        tier, subtier, evidence_score = tier_for_row(n_tools, hmm_c, sub_c, dia_c)

        row = {
            "protein_id": gid, "genome_id": genome_id, "tax_class": tax_class, "seq_length": seq_len,
            "ec_numbers": ec_numbers,
            "dbcan_hmm_family": dbcan_hmm_family,
            "dbcan_hmm_ievalue": "" if hmm_ievalue is None else f"{hmm_ievalue:.3g}",
            "dbcan_hmm_coverage": "" if hmm_coverage is None else f"{hmm_coverage:.4f}",
            "dbcan_sub_family": dbcan_sub_family,
            "dbcan_sub_ievalue": "" if sub_ievalue is None else f"{sub_ievalue:.3g}",
            "dbcan_sub_coverage": "" if sub_coverage is None else f"{sub_coverage:.4f}",
            "dbcan_sub_substrate": sub_substrate,
            "diamond_family": diamond_family,
            "diamond_pident": "" if dia_pident is None else f"{dia_pident:.2f}",
            "diamond_evalue": "" if dia_evalue is None else f"{dia_evalue:.3g}",
            "diamond_bitscore": "" if dia_bitscore is None else f"{dia_bitscore:.1f}",
            "diamond_cazy_hit_id": dia_cazy_id,
            "n_tools": n_tools,
            "recommend_family": recommend_family,
            "hmm_component": f"{hmm_c:.4f}",
            "sub_component": f"{sub_c:.4f}",
            "diamond_component": f"{dia_c:.4f}",
            "evidence_score": f"{evidence_score:.4f}",
            "tier": tier,
            "subtier": subtier,
            "structure_evidence_score": "",
        }
        writer.writerow(row)
        n_rows += 1

        tier_counts[(tax_class, tier, subtier)] += 1
        if tier == "gray_zone":
            fam = recommend_family if recommend_family != "-" else (
                dbcan_hmm_family if dbcan_hmm_family != "-" else (
                    dbcan_sub_family if dbcan_sub_family != "-" else diamond_family))
            fam = fam.split("+")[0] if fam and fam != "-" else "UNKNOWN"
            gray_family_counts[fam] += 1

        key = (tax_class, tier)
        reservoir.offer(key, {
            "protein_id": gid, "genome_id": genome_id, "tax_class": tax_class, "tier": tier,
            "subtier": subtier, "evidence_score": evidence_score, "n_tools": n_tools,
            "recommend_family": recommend_family, "dbcan_hmm_family": dbcan_hmm_family,
            "dbcan_sub_family": dbcan_sub_family, "diamond_family": diamond_family,
            "seq_length": seq_len,
        })
    return n_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="/array1/xinpeng/all_genome")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--sample-size", type=int, default=4000)
    ap.add_argument("--reservoir-per-stratum", type=int, default=300,
                     help="max proteins kept per (tax_class, tier) stratum before final downselect")
    ap.add_argument("--limit-genomes", type=int, default=0, help="debug: only process first N genome dirs")
    ap.add_argument("--seed", type=int, default=13)
    args = ap.parse_args()

    random.seed(args.seed)
    os.makedirs(args.out_dir, exist_ok=True)

    tier_counts = Counter()
    gray_family_counts = Counter()
    reservoir = ReservoirSampler(cap=args.reservoir_per_stratum)
    failed = []

    out_tsv = os.path.join(args.out_dir, "tiered_proteins.tsv.gz")
    t0 = time.time()
    n_genomes = 0
    n_rows_total = 0

    with gzip.open(out_tsv, "wt", newline="") as gz:
        writer = csv.DictWriter(gz, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()

        tax_classes = sorted(d for d in os.listdir(args.data_root)
                              if os.path.isdir(os.path.join(args.data_root, d)))
        genome_dirs = []
        for tc in tax_classes:
            tc_path = os.path.join(args.data_root, tc)
            for gid in sorted(os.listdir(tc_path)):
                gpath = os.path.join(tc_path, gid)
                if os.path.isdir(gpath):
                    genome_dirs.append((tc, gid, gpath))

        if args.limit_genomes:
            genome_dirs = genome_dirs[:args.limit_genomes]

        print(f"[build_tiered_dataset] {len(genome_dirs)} genome dirs to process", file=sys.stderr)

        for i, (tc, gid, gpath) in enumerate(genome_dirs):
            try:
                n_rows = process_genome(gpath, tc, gid, writer, tier_counts, gray_family_counts,
                                         reservoir, seq_cache_needed=None)
                n_rows_total += n_rows
                n_genomes += 1
            except Exception as e:
                failed.append((tc, gid, repr(e)))
                print(f"[build_tiered_dataset] FAILED {tc}/{gid}: {e!r}", file=sys.stderr)
                continue

            if (i + 1) % 100 == 0:
                elapsed = time.time() - t0
                print(f"[build_tiered_dataset] {i+1}/{len(genome_dirs)} genomes, "
                      f"{n_rows_total} rows, {elapsed:.0f}s elapsed", file=sys.stderr)

    print(f"[build_tiered_dataset] done: {n_genomes} genomes ok, {len(failed)} failed, "
          f"{n_rows_total} total protein rows, {time.time()-t0:.0f}s", file=sys.stderr)

    # --- summary: tier x class counts ---
    with open(os.path.join(args.out_dir, "summary_tier_by_class.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["tax_class", "tier", "subtier", "count"])
        for (tc, tier, subtier), cnt in sorted(tier_counts.items()):
            w.writerow([tc, tier, subtier, cnt])

    # --- summary: gray zone family distribution ---
    with open(os.path.join(args.out_dir, "summary_gray_zone_families.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family", "gray_zone_count"])
        for fam, cnt in gray_family_counts.most_common():
            w.writerow([fam, cnt])

    # --- failed genomes ---
    with open(os.path.join(args.out_dir, "failed_genomes.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["tax_class", "genome_id", "error"])
        for row in failed:
            w.writerow(row)

    tier_totals = Counter()
    for (tc, tier, subtier), cnt in tier_counts.items():
        tier_totals[tier] += cnt

    summary = {
        "data_root": args.data_root,
        "n_genomes_attempted": len(genome_dirs) if not args.limit_genomes else args.limit_genomes,
        "n_genomes_ok": n_genomes,
        "n_genomes_failed": len(failed),
        "n_protein_rows_total": n_rows_total,
        "tier_totals": dict(tier_totals),
        "elapsed_seconds": time.time() - t0,
        "reservoir_per_stratum": args.reservoir_per_stratum,
        "seed": args.seed,
    }
    with open(os.path.join(args.out_dir, "summary_overall.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))

    # --- final stratified sample: downselect from reservoir pool, then fetch sequences ---
    pool = reservoir.all_items()
    print(f"[build_tiered_dataset] reservoir pool size: {len(pool)}", file=sys.stderr)

    # Target allocation across tiers: weight gray_zone highest (primary interest),
    # but keep meaningful counts of the other two tiers for calibration/negative controls.
    tier_weight = {"gray_zone": 0.5, "high_confidence_cazyme": 0.3, "high_confidence_non_cazyme": 0.2}
    by_tier = defaultdict(list)
    for item in pool:
        by_tier[item["tier"]].append(item)

    target = args.sample_size
    selected = []
    for tier, w in tier_weight.items():
        items = by_tier.get(tier, [])
        random.shuffle(items)
        n_take = min(len(items), round(target * w))
        selected.extend(items[:n_take])

    # top up from whatever remains if under target (small pools in some tiers)
    if len(selected) < target:
        remaining = [it for tier_items in by_tier.values() for it in tier_items if it not in selected]
        random.shuffle(remaining)
        selected.extend(remaining[:target - len(selected)])

    print(f"[build_tiered_dataset] final sample size: {len(selected)}", file=sys.stderr)

    # fetch sequences: group selected by (tax_class, genome_id) so each uniInput.faa is opened once
    by_genome = defaultdict(list)
    for item in selected:
        by_genome[(item["tax_class"], item["genome_id"])].append(item)

    sample_tsv_path = os.path.join(args.out_dir, "sample_for_structure.tsv")
    sample_fasta_path = os.path.join(args.out_dir, "sample_for_structure.fasta")
    sample_fields = ["protein_id", "genome_id", "tax_class", "tier", "subtier", "evidence_score",
                      "n_tools", "recommend_family", "dbcan_hmm_family", "dbcan_sub_family",
                      "diamond_family", "seq_length", "structure_evidence_score"]

    with open(sample_tsv_path, "w", newline="") as tsv_fh, open(sample_fasta_path, "w") as fasta_fh:
        tsv_w = csv.DictWriter(tsv_fh, fieldnames=sample_fields, delimiter="\t")
        tsv_w.writeheader()
        for (tc, gid), items in by_genome.items():
            gpath = os.path.join(args.data_root, tc, gid, "uniInput.faa")
            seqs = parse_fasta_seqs(gpath)
            for item in items:
                seq = seqs.get(item["protein_id"], "")
                item_out = {k: item.get(k, "") for k in sample_fields}
                item_out["structure_evidence_score"] = ""
                tsv_w.writerow(item_out)
                fasta_fh.write(f">{item['protein_id']} tier={item['tier']} class={tc}\n")
                for k in range(0, len(seq), 80):
                    fasta_fh.write(seq[k:k+80] + "\n")

    print(f"[build_tiered_dataset] wrote sample TSV -> {sample_tsv_path}", file=sys.stderr)
    print(f"[build_tiered_dataset] wrote sample FASTA -> {sample_fasta_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
