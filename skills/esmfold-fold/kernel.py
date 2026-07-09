"""esmfold-fold kernel helpers."""

CANONICAL = "ACDEFGHIKLMNPQRSTVWY"


def sanitize_sequence(seq):
    """Uppercase, strip a trailing stop symbol, map non-canonical residues to X."""
    s = "".join(seq.split()).upper().rstrip("*")
    return "".join(c if c in CANONICAL else "X" for c in s)


def plddt_from_pdb(pdb_path):
    """Return (n_residues, mean_pLDDT, min_pLDDT, max_pLDDT) from CA B-factors."""
    vals = []
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    pass
    if not vals:
        return (0, 0.0, 0.0, 0.0)
    return (len(vals), sum(vals) / len(vals), min(vals), max(vals))


def esmfold_script(fasta_name="in.fasta", out_dir="pdb_out"):
    """Return a self-contained ESMFold folding script (string) to run as a GPU job.

    Loads facebook/esmfold_v1 once, folds every record in `fasta_name`
    shortest-first, writes <out_dir>/<id>.pdb (CA B-factor = pLDDT) and
    fold_summary.json. Sanitizes sequences; serializes folds for one GPU.
    """
    tmpl = r"""
import os, time, json, torch
from transformers import AutoTokenizer, EsmForProteinFolding

CANON = set("ACDEFGHIKLMNPQRSTVWY")
def sanitize(s):
    s = "".join(s.split()).upper().rstrip("*")
    return "".join(c if c in CANON else "X" for c in s)

def read_fasta(p):
    d, k, b = {}, None, []
    for ln in open(p):
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            if k: d[k] = "".join(b)
            k = ln[1:].split()[0]; b = []
        else:
            b.append(ln)
    if k: d[k] = "".join(b)
    return d

seqs = read_fasta("__FASTA__")
print("loaded", len(seqs), "sequences", flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
model = model.cuda().eval()
model.esm = model.esm.half()
model.trunk.set_chunk_size(64)
print("model loaded", flush=True)

os.makedirs("__OUT__", exist_ok=True)
summary = []
for pid, raw in sorted(seqs.items(), key=lambda kv: len(kv[1])):
    s = sanitize(raw); L = len(s); t0 = time.time()
    try:
        with torch.no_grad():
            pdb = model.infer_pdb(s)
        open("__OUT__/" + pid + ".pdb", "w").write(pdb)
        pls = [float(l[60:66]) for l in pdb.splitlines()
               if l.startswith("ATOM") and l[12:16].strip() == "CA"]
        mp = sum(pls) / len(pls) if pls else 0.0
        summary.append({"id": pid, "len": L, "mean_plddt": round(mp, 1),
                        "sec": round(time.time() - t0, 1)})
        print("OK", pid, "L=%d pLDDT=%.1f" % (L, mp), flush=True)
    except Exception as e:
        summary.append({"id": pid, "len": L, "mean_plddt": -1, "error": str(e)[:120]})
        print("FAIL", pid, type(e).__name__, str(e)[:120], flush=True)
    torch.cuda.empty_cache()
json.dump(summary, open("fold_summary.json", "w"))
print("DONE", len(summary), flush=True)
"""
    return tmpl.replace("__FASTA__", fasta_name).replace("__OUT__", out_dir)
