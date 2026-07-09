"""Kernel helpers for the clean-ec skill: run CLEAN (Yu et al., Science 2023)
for EC-number prediction from protein sequence, and parse its results.

CLEAN predicts Enzyme Commission numbers directly from sequence via contrastive
learning on ESM-1b embeddings, retrieving against EC-cluster centroids. It is an
INDEPENDENT sequence-based EC predictor -- orthogonal to inheriting an EC from a
CAZy/Pfam family assignment.

All helpers are stdlib-only; CLEAN itself runs as a subprocess in its own app
dir + isolated interpreter (see SKILL.md "Isolated environment" -- do NOT put
fair-esm in an env that also has EvolutionaryScale ESM-C).
"""
import os
import re
import csv
import shutil
import subprocess

# Canonical, tested install facts (from the CLEAN Dockerfile / README).
CLEAN_REPO = "https://github.com/tttianhao/CLEAN"
CLEAN_WEIGHTS_GDRIVE_ID = "1gsxjSf2CtXzgW1XsennTr-TcvSoTSDtk"  # pretrained.zip
ESM_REPO = "https://github.com/facebookresearch/esm"           # provides scripts/extract.py
# Pins that keep ESM-1b checkpoint loading + the sklearn GMM pickle working.
# torch<2.6 so torch.load defaults to weights_only=False (ESM-1b .pt needs it);
# sklearn 1.2.0 matches the GMM confidence pickle shipped in pretrained.zip.
CLEAN_TORCH_PIN = "torch==2.2.2"
VALID_AA = "ACDEFGHIKLMNPQRSTVWY"


def sanitize_fasta(in_path, out_path, entry_id=None):
    """Write a CLEAN-safe FASTA: single-token header, non-standard residues->X.

    CLEAN maps the header (up to first whitespace) to the ESM embedding filename
    and to the query id, so a header with spaces silently breaks the id match.
    ESM-1b also rejects non-standard residues (e.g. '*'). One record in/out.
    """
    raw = open(in_path).read().strip().splitlines()
    hdr = raw[0][1:].strip() if raw and raw[0].startswith(">") else "query"
    if entry_id is None:
        entry_id = hdr.split()[0] if hdr.split() else "query"
    seq = "".join(l.strip() for l in raw if not l.startswith(">"))
    seq = re.sub("[^" + VALID_AA + "]", "X", seq.upper())
    with open(out_path, "w") as fh:
        fh.write(">" + entry_id + "\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    return {"entry_id": entry_id, "length": len(seq)}


def clean_infer_command(fasta_name, model="split100", python_bin="python"):
    """Return the shell one-liner that runs max-separation inference inside a
    CLEAN app dir on data/inputs/<fasta_name>.fasta.

    Generalizes CLEAN_infer_fasta.py to any pretrained split (split100/split70)
    while keeping GMM confidence. Reusable both locally (subprocess) and in a
    remote submit_job `command`. Assumes cwd == <CLEAN>/app and that
    data/inputs/<fasta_name>.fasta already exists.
    """
    td = "inputs/" + fasta_name
    py = (
        "import sys; sys.path.insert(0,'src');"
        "from CLEAN.utils import prepare_infer_fasta;"
        "from CLEAN.infer import infer_maxsep;"
        "prepare_infer_fasta('%s');"
        "infer_maxsep('%s','%s',report_metrics=False,pretrained=True,"
        "gmm='./data/pretrained/gmm_ensumble.pkl')" % (td, model, td)
    )
    return "mkdir -p data/inputs data/esm_data results/inputs && %s -c \"%s\"" % (python_bin, py)


def run_clean(fasta, workdir=None, clean_app_dir=None, python_bin=None, model="split100"):
    """Run CLEAN locally (subprocess) on a single-record FASTA and return where
    the result CSV landed. Requires an already-installed CLEAN (see SKILL.md
    "Isolated environment" / clean_setup_script). For a REMOTE CLEAN install,
    don't call this -- build the command with clean_infer_command() and pass it
    to host.compute submit_job instead.

    Returns {results_csv, entry_id, returncode, cmd, stdout_tail}.
    """
    if clean_app_dir is None:
        raise ValueError(
            "clean_app_dir is required: path to <CLEAN>/app of an installed CLEAN. "
            "Install it first (see SKILL.md); or for a remote host use clean_infer_command()+submit_job.")
    if python_bin is None:
        python_bin = "python"
    app = os.path.abspath(clean_app_dir)
    name = os.path.splitext(os.path.basename(fasta))[0]
    dst = os.path.join(app, "data", "inputs", name + ".fasta")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    os.makedirs(os.path.join(app, "results", "inputs"), exist_ok=True)
    info = sanitize_fasta(fasta, dst, entry_id=name)
    cmd = clean_infer_command(name, model=model, python_bin=python_bin)
    env = dict(os.environ)
    env["PYTHONPATH"] = os.path.join(app, "src") + os.pathsep + env.get("PYTHONPATH", "")
    proc = subprocess.run(cmd, shell=True, cwd=app, env=env,
                          capture_output=True, text=True)
    out_csv = os.path.join(app, "results", "inputs", name + "_maxsep.csv")
    return {
        "results_csv": out_csv if os.path.exists(out_csv) else None,
        "entry_id": info["entry_id"],
        "returncode": proc.returncode,
        "cmd": cmd,
        "stdout_tail": (proc.stdout or "")[-1500:] + (proc.stderr or "")[-800:],
    }


def parse_clean_results(path):
    """Parse a CLEAN *_maxsep.csv into a list of per-query dicts.

    CLEAN row format: query, EC:<ec>/<score>, EC:<ec>/<score>, ...
    <score> is a GMM confidence in [0,1] (higher=more confident) when inference
    used gmm=..., else a Euclidean distance (lower=closer). We report the raw
    number as `score` and flag which regime it looks like.

    Returns [{query, predictions:[{rank, ec_number, score, raw}], top_ec,
    top_score, score_type}].
    """
    rows = []
    with open(path, newline="") as fh:
        for parts in csv.reader(fh):
            if not parts or not parts[0].strip():
                continue
            query = parts[0].strip()
            preds = []
            for rank, tok in enumerate(parts[1:], start=1):
                tok = tok.strip()
                if not tok:
                    continue
                m = re.match(r"EC:([0-9n.\-]+)/([-0-9.eE]+)", tok)
                if m:
                    ec, sc = m.group(1), float(m.group(2))
                else:  # be lenient about format drift
                    ec = tok.replace("EC:", "").split("/")[0]
                    try:
                        sc = float(tok.split("/")[-1])
                    except ValueError:
                        sc = float("nan")
                preds.append({"rank": rank, "ec_number": ec, "score": sc, "raw": tok})
            scores = [p["score"] for p in preds if p["score"] == p["score"]]
            score_type = "gmm_confidence" if scores and max(scores) <= 1.0 else "distance"
            rows.append({
                "query": query,
                "predictions": preds,
                "top_ec": preds[0]["ec_number"] if preds else None,
                "top_score": preds[0]["score"] if preds else None,
                "score_type": score_type,
            })
    return rows


def compare_ec(predicted_ec, family_ec):
    """Compare two EC numbers at increasing specificity.

    Returns {level} in {exact, class3 (first 3 digits), class2, class1,
    different, incomparable}. Family EC is a reference for agreement checking,
    not ground truth.
    """
    if not predicted_ec or not family_ec:
        return {"level": "incomparable"}
    p = predicted_ec.replace("EC:", "").split(".")
    f = family_ec.replace("EC:", "").split(".")
    if p == f:
        return {"level": "exact"}
    for n, name in ((3, "class3"), (2, "class2"), (1, "class1")):
        if p[:n] == f[:n]:
            return {"level": name}
    return {"level": "different"}
