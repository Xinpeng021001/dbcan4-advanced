---
name: clean-ec
description: >-
  Predict Enzyme Commission (EC) numbers directly from a protein amino-acid
  sequence using CLEAN (Yu et al., Science 2023) — contrastive learning on
  ESM-1b embeddings retrieved against EC-cluster centroids. Use this skill
  whenever you need a sequence-based EC-number prediction with a confidence
  score for one or many proteins from a FASTA file: enzyme function
  annotation, checking or supplementing a family-derived EC (e.g. an EC
  inherited from a CAZy/Pfam family), screening putative enzymes, or asking
  "what reaction does this protein catalyze?". CLEAN gives an INDEPENDENT
  prediction orthogonal to family/homology assignment, so reach for it even
  when a family-based EC already exists and you want a second, model-based
  opinion. Triggers on: "predict EC number", "what enzyme class", "EC from
  sequence", "run CLEAN", "enzyme function prediction", "annotate enzyme
  commission".
---

# clean-ec: EC-number prediction from sequence with CLEAN

CLEAN (**C**ontrastive **Lea**rning enabled Enzyme **AN**notation; Yu et al.,
*Science* 2023, doi:10.1126/science.adf2465; repo
https://github.com/tttianhao/CLEAN) assigns Enzyme Commission (EC) numbers to a
protein from its amino-acid sequence. It embeds each sequence with ESM-1b, maps
it through a contrastively-trained projection, and predicts the EC number(s)
whose training-set cluster centroids are closest, using a greedy
**max-separation** rule. A GMM converts the retrieval distance into a
confidence in [0, 1].

This is a **sequence-based, model-driven** EC predictor. It does not need a
family assignment, a structure, or a BLAST hit — so it is a genuinely
independent line of evidence you can set against a family-inherited EC.

## When to use it

Any time you have a protein sequence (a FASTA) and want its EC number(s) plus a
confidence. Common cases: annotating hypothetical enzymes, giving a second
opinion on an EC that was inherited from a CAZy/Pfam family, or flagging
proteins whose predicted function disagrees with their family.

## Output format (what CLEAN produces)

CLEAN writes `results/inputs/<name>_maxsep.csv`. One row per query:

```
267317,EC:3.2.1.40/0.9876,EC:3.2.1.20/0.4013
```

- First field: the query id (the FASTA header up to the first space).
- Each further field: `EC:<number>/<score>`. Multiple fields = multiple
  predicted functions (max-separation can call several).
- `<score>` semantics depend on how inference was run:
  - **With `gmm=...` (the default here):** a GMM confidence in **[0, 1]**,
    higher = more confident.
  - **Without gmm:** a **Euclidean distance**, lower = closer/better.

`parse_clean_results()` reports the raw number as `score` and sets
`score_type` to `gmm_confidence` or `distance` by inspecting the range, so
you never misread a distance as a confidence.

## Kernel helpers (auto-loaded with this skill)

- `sanitize_fasta(in_path, out_path, entry_id=None)` — write a CLEAN-safe FASTA
  (single-token header, non-standard residues → `X`). Do this before every run:
  CLEAN maps the header to the ESM embedding filename, so a header with spaces
  silently breaks the id match, and ESM-1b rejects non-standard residues.
- `clean_infer_command(fasta_name, model="split100", python_bin="python")` —
  return the exact shell command for max-separation inference on
  `data/inputs/<fasta_name>.fasta`. Use its string both locally and as a remote
  `submit_job` command. Generalizes the repo's `CLEAN_infer_fasta.py` to any
  pretrained split while keeping GMM confidence.
- `run_clean(fasta, clean_app_dir=..., python_bin=..., model="split100")` — run
  an already-installed **local** CLEAN as a subprocess; returns
  `{results_csv, entry_id, returncode, cmd, stdout_tail}`.
- `parse_clean_results(path) -> list[dict]` — parse the maxsep CSV. Each query:
  `{query, predictions:[{rank, ec_number, score, raw}], top_ec, top_score,
  score_type}`.
- `compare_ec(predicted_ec, family_ec) -> {level}` — agreement at increasing
  specificity: `exact` / `class3` (first 3 EC digits) / `class2` / `class1` /
  `different`. Handy for "does CLEAN agree with the family EC?".

`model` picks the pretrained split: `"split100"` (100 %-identity clustering of
SwissProt, the default/recommended) or `"split70"`.

## Isolated environment — read this before installing

CLEAN depends on **`fair-esm`** (Facebook ESM-1b). `fair-esm` installs a
top-level `esm` package that **collides with EvolutionaryScale's `esm`
(ESM-C / ESM-3)**. If both land in one environment, the ESM-C import breaks.

So **always install CLEAN in a dedicated interpreter** — a fresh
`python -m venv` or a fresh conda env — never into an environment that has (or
will have) EvolutionaryScale ESM. Point HF/torch caches at scratch; the ESM-1b
checkpoint alone is ~7.8 GB.

Version pins that matter (baked into the setup script below):
- `torch==2.2.2` (CPU or CUDA). torch ≥ 2.6 flips `torch.load` to
  `weights_only=True` by default, which **fails to load the ESM-1b
  checkpoint**. Staying < 2.6 avoids it.
- The repo's pinned `requirements.txt` (numpy 1.22.3, pandas 1.4.2,
  scikit-learn 1.2.0, fair-esm 2.0.0). scikit-learn 1.2.0 matches the shipped
  GMM confidence pickle — a much newer sklearn can fail to unpickle it.

## Install recipe (tested, idempotent)

Run inside a dedicated env. `$BASE` is a scratch dir with ~10 GB free. Every
step is safe to re-run, so a mid-download timeout just resumes.

```bash
set -o pipefail
BASE=/path/to/scratch                 # >= 10 GB free
CLEAN_DIR="$BASE/CLEAN"; VENV="$BASE/venv_clean"
export TORCH_HOME="$BASE/clean_cache/torch" HF_HOME="$BASE/clean_cache/hf"
mkdir -p "$BASE/clean_cache"

# 1. clone
[ -d "$CLEAN_DIR/.git" ] || git clone --depth 1 https://github.com/tttianhao/CLEAN "$CLEAN_DIR"
cd "$CLEAN_DIR/app"

# 2. DEDICATED venv (never share with EvolutionaryScale ESM)
[ -x "$VENV/bin/python" ] || python3 -m venv "$VENV"
source "$VENV/bin/activate"
python -m pip install -q --upgrade pip setuptools wheel

# 3. pinned deps + torch<2.6 CPU wheel (use the CUDA index for GPU)
pip install -q -r requirements.txt
pip install -q "torch==2.2.2" --extra-index-url https://download.pytorch.org/whl/cpu
pip install -q "gdown>=5.1.0"

# 4. ESM repo (provides scripts/extract.py) + dirs
[ -d esm ] || git clone --depth 1 https://github.com/facebookresearch/esm
mkdir -p data/esm_data data/pretrained data/inputs results/inputs

# 5. pretrained weights (split100/split70 .pth + gmm_ensumble.pkl + 100.pt/70.pt)
if [ ! -f data/pretrained/split100.pth ]; then
  gdown "https://drive.google.com/uc?id=1gsxjSf2CtXzgW1XsennTr-TcvSoTSDtk" -O pretrained.zip
  unzip -o -q pretrained.zip -d data/pretrained
  for d in data/pretrained/CLEAN_pretrained*/; do [ -d "$d" ] && mv "$d"* data/pretrained/ || true; done
fi

# 6. build the CLEAN package (PYTHONPATH=src also works if build is skipped)
python build.py install 2>&1 | tail -3 || true
export PYTHONPATH="$CLEAN_DIR/app/src:$PYTHONPATH"
python -c "import CLEAN.infer, CLEAN.utils; print('CLEAN OK')"
```

The first ESM-1b embedding call downloads the ~7.8 GB checkpoint into
`TORCH_HOME`; later runs reuse it.

## Running inference

Once installed, per query: stage a sanitized FASTA into `data/inputs/`, then run
the max-separation command from `<CLEAN>/app`.

```python
# CLEAN installed locally:
sanitize_fasta("my_protein.fasta", "staged.fasta")   # header hygiene
r = run_clean("staged.fasta", clean_app_dir="<BASE>/CLEAN/app",
              python_bin="<BASE>/venv_clean/bin/python", model="split100")
preds = parse_clean_results(r["results_csv"])
```

### Running on a remote host (recommended for the heavy ESM-1b step)

CLEAN's cost is the ESM-1b embedding. On a cluster, do install + inference in
**one** `submit_job` so the big download and the run share a workdir. Build the
command with `clean_infer_command()` and append it to the install recipe:

```python
# repl tool
c = host.compute.create("ssh:<host>")
cmd = clean_infer_command("267317", model="split100",
                          python_bin="$VENV/bin/python")   # from <CLEAN>/app
job = c.submit_job(command=SETUP_RECIPE + "\ncd $CLEAN_DIR/app\n" + cmd +
                   "\ncp results/inputs/267317_maxsep.csv $SUBMIT_DIR/",
                   intent="CLEAN EC prediction on 267317",
                   inputs=[{"src": "267317.fasta", "dst_filename": "267317.fasta"}],
                   outputs=[{"glob": "*_maxsep.csv", "visibility": "featured"}],
                   timeout_seconds=3000)
```

Then parse the harvested CSV with `parse_clean_results()`.

## Interpreting agreement with a family-based EC

When a protein already has an EC inherited from a family assignment, use
`compare_ec(clean_top_ec, family_ec)`. Report the agreement level honestly —
`exact`, `class3` (same first three EC digits = same reaction class, different
substrate specificity), or `different`. Disagreement is informative, not an
error: CLEAN and family assignment are independent, and either can be right.

## Honesty

CLEAN's confidence can be low or its call wrong, especially for sequences far
from SwissProt. Always report the real predicted EC(s), the real score, and the
`score_type`. If setup or a run fails, report the actual error — never invent an
EC, a score, or an agreement result.
