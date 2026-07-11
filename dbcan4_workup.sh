#!/usr/bin/env bash
# =============================================================================
# dbcan4_workup.sh — ONE command: protein FASTA -> full CAZyme workup -> web UI
# =============================================================================
# Runs the complete, validated real workup on met and (optionally) serves the
# browsable BioForge web UI:
#     baseline dbCAN (run_dbcan V5)
#   + advanced ESM-C family calls (kNN / centroid / contrastive -> fusion)
#   + Pfam domains (hmmscan)          + physicochemistry (Biopython)
#   + EC number (CLEAN)               + TM topology & signal peptide (DeepTMHMM)
#   + subcellular localization (derived from DeepTMHMM; honest fallback)
#   -> v1.1 output contract -> BioForge SQLite -> FastAPI web UI
#
# Usage:
#     dbcan4_workup.sh <proteins.faa> [--outdir DIR] [--sample NAME]
#                      [--serve] [--port N] [--gpu N] [--no-deeptmhmm] [--no-clean]
#
# Example (the whole product in one line):
#     dbcan4_workup.sh my_proteins.faa --serve
# =============================================================================
set -euo pipefail

# ---- defaults (met paths; all overridable) ----
REPO=/array1/xinpeng/dbcan4-advanced
NF=$REPO/repo_clone/nf
VENV=$REPO/venv                                   # engine venv (torch + esm)
BIODB_VENV=/array1/xinpeng/scratch/biodb_venv     # web stack + dbcan4 CLI
BIODB_SRC=/array1/xinpeng/biodb
BIOLIB_BIN=/array1/xinpeng/scratch/venv_biolib/bin/biolib
DBCAN_DB=/array1/xinpeng/dbcan_db
PFAM=/array1/xinpeng/pfam/Pfam-A.hmm

FAA=""; OUTDIR=""; SAMPLE=""; SERVE=0; PORT=8000; GPU=0
DO_DEEPTMHMM=1; DO_CLEAN=1
HOST=127.0.0.1

# ---- arg parse ----
while [ $# -gt 0 ]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --serve) SERVE=1; shift;;
    --port) PORT="$2"; shift 2;;
    --gpu) GPU="$2"; shift 2;;
    --host) HOST="$2"; shift 2;;
    --no-deeptmhmm) DO_DEEPTMHMM=0; shift;;
    --no-clean) DO_CLEAN=0; shift;;
    -h|--help) sed -n '2,25p' "$0"; exit 0;;
    -*) echo "unknown option: $1" >&2; exit 2;;
    *) FAA="$1"; shift;;
  esac
done
[ -n "$FAA" ] || { echo "ERROR: provide a protein FASTA. See --help." >&2; exit 2; }
[ -f "$FAA" ] || { echo "ERROR: FASTA not found: $FAA" >&2; exit 2; }
FAA=$(cd "$(dirname "$FAA")" && pwd)/$(basename "$FAA")   # absolutize
SAMPLE="${SAMPLE:-$(basename "${FAA%.*}")}"
OUTDIR="${OUTDIR:-$(pwd)/${SAMPLE}_workup}"
mkdir -p "$OUTDIR"; cd "$OUTDIR"
OUT="$OUTDIR/results"
FEAT="$OUT/cazyme_advanced/features/$SAMPLE"
PRED="$OUT/cazyme_advanced/predictions/$SAMPLE"
mkdir -p "$FEAT" "$PRED"

export PATH="$NF/bin:$REPO/scripts:$PATH"
export HF_HOME=$REPO/hf_cache
export CUDA_VISIBLE_DEVICES=$GPU
NPROT=$(grep -c '^>' "$FAA" || echo 0)
echo "### dbcan4 full workup | sample=$SAMPLE | $NPROT proteins | outdir=$OUTDIR"

echo "=== 1/8 baseline dbCAN (run_dbcan V5) ==="
"$VENV/bin/run_dbcan" CAZyme_annotation --mode protein \
    --input_raw_data "$FAA" --db_dir "$DBCAN_DB" \
    --output_dir rundbcan --methods diamond,hmm,dbCANsub
"$VENV/bin/python" "$NF/bin/emit_baseline_funcscan.py" \
    --overview rundbcan/overview.tsv --faa "$FAA" --sample "$SAMPLE" --outdir "$OUT/funcscan"

echo "=== 2/8 advanced ESM-C embed + label-free infer ==="
"$VENV/bin/python" "$REPO/scripts/embed_esmc.py" --fasta "$FAA" \
    --out-prefix query.esmc --model esmc_600m --nshards 1 --shard 0
[ -f query.esmc.shard0.npz ] && mv -f query.esmc.shard0.npz query.esmc.npz || true
"$VENV/bin/python" "$REPO/scripts/infer_esmc.py" --emb query.esmc.npz \
    --ref-prefix "$REPO/emb/ref2024" --heads "$REPO/results/heads/heads.pt" \
    --proj-ref "$REPO/results/heads/proj_ref.npz" --k 15 \
    --out-knn raw_knn.tsv --out-centroid raw_centroid.tsv --out-contrastive raw_contrastive.tsv
for pair in "ESM-C-kNN:raw_knn.tsv" "ESM-C-centroid:raw_centroid.tsv" "ESM-C-contrastive:raw_contrastive.tsv"; do
  "$VENV/bin/python" "$NF/bin/normalize_predictions.py" --tool "${pair%%:*}" --in "${pair##*:}" --out "$PRED/${pair%%:*}.tsv"
done

echo "=== 3/8 fusion consensus ==="
"$VENV/bin/python" "$NF/bin/fuse_predictions.py" --inputs "$PRED/"*.tsv \
    --out fusion_raw.tsv --weights '{}' --min-confidence 0.0
"$VENV/bin/python" "$NF/bin/normalize_predictions.py" --tool fusion --in fusion_raw.tsv --out "$PRED/fusion.tsv"

echo "=== 4/8 Pfam domains (hmmscan) ==="
hmmscan --domtblout domains.domtbl --cut_ga -o /dev/null "$PFAM" "$FAA"
"$VENV/bin/python" "$NF/bin/feature_converters.py" domains --domtbl domains.domtbl --out "$FEAT/domains.tsv"

echo "=== 5/8 physicochemistry (Biopython) ==="
"$VENV/bin/python" "$NF/bin/feature_converters.py" physicochem --faa "$FAA" --out "$FEAT/physicochem.tsv"

if [ "$DO_DEEPTMHMM" = 1 ]; then
  echo "=== 6/8 DeepTMHMM (TM topology + signal peptide, BioLib cloud) ==="
  mkdir -p dtm && cp -f "$FAA" dtm/query.faa
  ( cd dtm && "$BIOLIB_BIN" run DTU/DeepTMHMM --fasta query.faa )
  "$NF/bin/deeptmhmm_to_tsv.py" --gff3 dtm/biolib_results/TMRs.gff3 \
      --three-line dtm/biolib_results/predicted_topologies.3line \
      --out-tm "$FEAT/deeptmhmm.tsv" --out-sp "$FEAT/signalp6.tsv"
  echo "=== 6b/8 localization (derived from DeepTMHMM SP; honest fallback, not DeepLoc) ==="
  "$VENV/bin/python" "$NF/bin/feature_converters.py" localization --signalp "$FEAT/signalp6.tsv" --out "$FEAT/localization.tsv"
else
  echo "=== 6/8 DeepTMHMM SKIPPED (--no-deeptmhmm) ==="
fi

if [ "$DO_CLEAN" = 1 ]; then
  echo "=== 7/8 CLEAN EC number ==="
  "$NF/bin/run_clean.sh" "$FAA" "$OUTDIR/clean_out.csv"
  "$VENV/bin/python" "$NF/bin/feature_converters.py" ec_prediction --clean "$OUTDIR/clean_out.csv" --out "$FEAT/ec_prediction.tsv"
else
  echo "=== 7/8 CLEAN SKIPPED (--no-clean) ==="
fi

echo "=== 8/8 write v1.1 manifest ==="
STAGE="$OUTDIR/stage"; rm -rf "$STAGE"; mkdir -p "$STAGE"
cp "$PRED/"*.tsv "$STAGE/" 2>/dev/null || true
cp "$FEAT/"*.tsv "$STAGE/" 2>/dev/null || true
( cd "$STAGE" && "$VENV/bin/python" "$NF/bin/write_manifest.py" \
    --outdir "$OUT" --sample "$SAMPLE" --stage-dir . \
    --release-label "workup-$(date +%Y-%m-%d)" --pipeline-version '0.1.0' \
    --tool-versions '{"esm":"3.2.1","deeptmhmm":"1.0.24","clean":"maxsep","hmmscan":"3.4","biopython":"1.85"}' )
echo "manifest -> $OUT/cazyme_advanced/manifest.json"

if [ "$SERVE" = 1 ]; then
  echo "=== ingest + serve ==="
  source "$BIODB_VENV/bin/activate"
  DB="$OUTDIR/${SAMPLE}.db"; rm -f "$DB"; export DATABASE_URL="sqlite:///$DB"
  export BIOFORGE_TRACKS_DIR="$OUTDIR/web_static/tracks"
  ( cd "$BIODB_SRC" && alembic upgrade head >/dev/null 2>&1 )
  bioforge-ingest          "$OUT/funcscan"
  bioforge-ingest-advanced "$OUT/cazyme_advanced/manifest.json"
  echo "### serving http://$HOST:$PORT  (Ctrl-C to stop)"
  exec uvicorn bioforge.api.main:app --host "$HOST" --port "$PORT"
else
  echo "### done. Ingest + serve with:  --serve   (or dbcan4 run --serve on the manifest)"
fi
