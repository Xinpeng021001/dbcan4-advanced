#!/usr/bin/env python3
"""Turn an ESMFold run (PDBs + fold_manifest.tsv) into the BioForge `structures.tsv`.

fold_esmfold.py writes one PDB per protein with per-residue pLDDT in the B-factor
column on a **0-1** scale, plus a manifest TSV (id, header, length, plddt_mean,
seconds, status). The BioForge advanced ingester (parse_advanced.py, ft=="structure")
expects a `structures.tsv` with columns:

    protein_id   plddt   path   length   source   extra

and the served 3Dmol viewer colours by B-factor on a **0-100** pLDDT scale
(min 50, max 90). So this converter:

  1. copies each PDB into <out-dir>/structures/<id>.pdb with B-factors rescaled
     ×100 (0-1 -> 0-100) so the viewer's colour ramp is correct;
  2. writes structures.tsv next to the other feature TSVs, with `plddt` = mean
     pLDDT ×100 and `path` = "structures/<id>.pdb" (relative to the TSV's parent,
     which is what the loader resolves against).

Honest: pLDDT is ESMFold's own confidence; nothing is fabricated. Proteins the
fold step skipped (too long / failed) simply get no row.
"""
from __future__ import annotations
import argparse, csv, json, os, shutil


def _rescale_bfactors(src_pdb: str, dst_pdb: str) -> None:
    """Copy a PDB, multiplying the B-factor column (cols 61-66) by 100 if 0-1."""
    # decide scale from the max B-factor seen
    mx = 0.0
    with open(src_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    mx = max(mx, float(line[60:66]))
                except ValueError:
                    pass
    scale = 100.0 if mx <= 1.5 else 1.0
    with open(src_pdb) as fi, open(dst_pdb, "w") as fo:
        for line in fi:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
                try:
                    b = float(line[60:66]) * scale
                    line = line[:60] + f"{min(b,999.99):6.2f}" + line[66:]
                except ValueError:
                    pass
            fo.write(line)


def _mean_bfactor(pdb: str) -> float:
    """Mean CA-atom pLDDT recovered from the PDB B-factor column.

    ESMFold stores per-residue pLDDT in the B-factor column, so the mean over
    CA atoms is one value per residue -- matching fold_esmfold.py's own
    ``out["plddt"].mean()`` (the manifest ``plddt_mean``) and the canonical
    per-residue pLDDT reported everywhere else in the product. Averaging over
    ALL atoms instead would weight residues by atom count and drift a few points
    low (e.g. 73.1 vs the true 76.6 for the multidomain example), so this is
    deliberately CA-only. Falls back to all-atom only if a PDB somehow has no CA
    records, so a value is still recovered rather than 0."""
    ca, allatom = [], []
    with open(pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
                try:
                    b = float(line[60:66])
                except ValueError:
                    continue
                allatom.append(b)
                if line[12:16].strip() == "CA":
                    ca.append(b)
    vals = ca or allatom
    return sum(vals) / len(vals) if vals else 0.0


def _plddt_0_100(plddt_mean, src_pdb: str) -> float:
    """pLDDT on a 0-100 scale. Prefer the manifest's plddt_mean; when it is
    blank/unparseable — which is exactly what fold_esmfold.py writes for a
    `cached` row on any re-run — fall back to the source PDB's B-factor column.
    Either input scale (0-1 or 0-100) is normalized to 0-100."""
    v = None
    if plddt_mean not in (None, ""):
        try:
            v = float(plddt_mean)
        except (TypeError, ValueError):
            v = None
    if v is None:
        v = _mean_bfactor(src_pdb)
    return round(v * 100, 1) if v <= 1.5 else round(v, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--esmfold-dir", required=True, help="dir with <id>.pdb + fold_manifest.tsv")
    ap.add_argument("--manifest", default=None, help="fold manifest TSV (default <esmfold-dir>/fold_manifest.tsv)")
    ap.add_argument("--out-dir", required=True, help="feature dir to write structures.tsv + structures/ into")
    ap.add_argument("--source", default="ESMFold", help="structure source label")
    args = ap.parse_args()

    man = args.manifest or os.path.join(args.esmfold_dir, "fold_manifest.tsv")
    struct_out = os.path.join(args.out_dir, "structures")
    os.makedirs(struct_out, exist_ok=True)
    tsv_out = os.path.join(args.out_dir, "structures.tsv")

    rows = []
    with open(man) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pid = r.get("id")
            if not pid:
                continue
            # A structure exists iff there is a real (non-empty) PDB on disk.
            # Do NOT gate on the manifest `status` string: a re-run logs
            # status="cached" (with a blank plddt_mean) for every PDB folded on
            # a previous pass, while genuine non-folds (skip_len / oom / err:*)
            # simply never wrote a PDB. The old `status != "ok"` skip therefore
            # dropped all cached structures, leaving structures.tsv header-only
            # and removing the entire 3D-structure card (and its colour dropdown)
            # from the gene page on any second run. File existence is the
            # invariant fold_esmfold.py itself uses for the cached decision.
            src_pdb = os.path.join(args.esmfold_dir, f"{pid}.pdb")
            if not (os.path.exists(src_pdb) and os.path.getsize(src_pdb) > 0):
                continue
            _rescale_bfactors(src_pdb, os.path.join(struct_out, f"{pid}.pdb"))
            plddt = _plddt_0_100(r.get("plddt_mean"), src_pdb)
            rows.append({
                "protein_id": pid,
                "plddt": plddt,
                "path": f"structures/{pid}.pdb",
                "length": r.get("length", ""),
                "source": args.source,
                "extra": json.dumps({"model": "facebook/esmfold_v1",
                                     "note": f"ESMFold prediction (mean pLDDT {plddt}/100)"}),
            })

    with open(tsv_out, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=["protein_id", "plddt", "path", "length", "source", "extra"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"[structures] {tsv_out}: {len(rows)} structure row(s); PDBs -> {struct_out} (B-factors 0-100)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
