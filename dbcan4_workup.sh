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
#   + 3D structure (ESMFold, local GPU) + subcellular localization (derived)
#   + InterPro domains + GO terms (real InterProScan if configured, else Pfam->GO)
#   -> v1.1 output contract -> BioForge SQLite -> FastAPI web UI
#
# dbCAN4 is fungal + protein-input: no genome, no Prokka, no GFF. Genes are built
# straight from the protein FASTA. A GFF is OPTIONAL — pass --gff only if you have
# genuine genomic coordinates for these proteins.
#
# The web UI (BioForge) is VENDORED into this repo (src/bioforge) — no second repo
# to clone. Every path below is overridable via an env var; see the README section
# "Environment variables & portability".
#
# Usage:
#     dbcan4_workup.sh <proteins.faa> [--outdir DIR] [--sample NAME] [--gff FILE]
#                      [--serve] [--port N] [--gpu N] [--host H]
#                      [--no-deeptmhmm] [--no-clean] [--no-structure]
#                      [--no-ai-report] [--no-interpro]
#
# Example (the whole product in one line):
#     dbcan4_workup.sh my_proteins.faa --serve
# =============================================================================
set -euo pipefail

# =============================================================================
# Paths — every one is overridable via an environment variable, so this script
# runs unchanged after a fresh `git clone` on any machine. Defaults point at the
# reference host (met). See README "Environment variables & portability".
# =============================================================================
# SELF_DIR = the checkout this script lives in (portable; via BASH_SOURCE).
# Resolve symlinks first so it works even when the script is symlinked onto PATH.
_SELF_SRC="${BASH_SOURCE[0]}"
_SELF_SRC="$(readlink -f "$_SELF_SRC" 2>/dev/null || echo "$_SELF_SRC")"
SELF_DIR="$(cd "$(dirname "$_SELF_SRC")" && pwd)"
# DBCAN4_ROOT = where the heavy data assets live (engine venv, ESM-C reference
# index, trained heads, HF cache). On met these sit one level ABOVE the checkout;
# set DBCAN4_ROOT to wherever you unpacked/downloaded them on another machine.
DBCAN4_ROOT="${DBCAN4_ROOT:-$(cd "$SELF_DIR/.." && pwd)}"

REPO="$DBCAN4_ROOT"                                          # asset root (emb/, results/heads/, hf_cache/)
NF="$SELF_DIR/nf"                                            # Nextflow pipeline + bin/ (in the checkout)
SCRIPTS="$SELF_DIR/scripts"                                  # stage scripts (in the checkout)
VENV="${DBCAN4_ENGINE_VENV:-$DBCAN4_ROOT/venv}"             # engine venv (torch + esm + run_dbcan)
# Web layer: bioforge is VENDORED into this checkout (src/bioforge), so the web
# venv just needs `pip install -e .` in the checkout. Defaults to the engine venv.
BIODB_VENV="${BIODB_VENV:-$VENV}"                           # venv where `bioforge` is importable
BIODB_SRC="${BIODB_SRC:-$SELF_DIR}"                         # dir with alembic.ini + db/alembic (this checkout)
BIOLIB_BIN="${BIOLIB_BIN:-/array1/xinpeng/scratch/venv_biolib/bin/biolib}"   # DeepTMHMM / BioLib CLI
DBCAN_DB="${DBCAN_DB:-/array1/xinpeng/dbcan_db}"           # run_dbcan V5 database (~7.4 GB)
PFAM="${PFAM_HMM:-/array1/xinpeng/pfam/Pfam-A.hmm}"        # Pfam-A HMM (§2.5 domains + GO fallback)
# InterProScan / GO. Real InterProScan is optional (tens of GB); if INTERPROSCAN_SH
# is unset or not executable we derive GO + InterPro from the Pfam domains offline.
INTERPROSCAN_SH="${INTERPROSCAN_SH:-}"
PFAM2GO="${PFAM2GO:-$NF/assets/mappings/pfam2go.txt}"
PFAM2INTERPRO="${PFAM2INTERPRO:-}"                          # optional PF->IPR map (fills InterPro accession col)

FAA=""; OUTDIR=""; SAMPLE=""; SERVE=0; PORT=8000; GPU=0
DO_DEEPTMHMM=1; DO_CLEAN=1; DO_STRUCTURE=1; DO_AI_REPORT=1; DO_INTERPRO=1
GFF=""                                            # OPTIONAL user-provided genomic GFF
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
    --no-structure) DO_STRUCTURE=0; shift;;
    --no-ai-report) DO_AI_REPORT=0; shift;;
    --no-interpro) DO_INTERPRO=0; shift;;
    --gff) GFF="$2"; shift 2;;
    -h|--help) sed -n '2,30p' "$0"; exit 0;;
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

export PATH="$NF/bin:$SCRIPTS:$PATH"
export HF_HOME=$REPO/hf_cache
export CUDA_VISIBLE_DEVICES=$GPU
NPROT=$(grep -c '^>' "$FAA" || echo 0)
echo "### dbcan4 full workup | sample=$SAMPLE | $NPROT proteins | outdir=$OUTDIR"

echo "=== 1/8 baseline dbCAN (run_dbcan V5) ==="
"$VENV/bin/run_dbcan" CAZyme_annotation --mode protein \
    --input_raw_data "$FAA" --db_dir "$DBCAN_DB" \
    --output_dir rundbcan --methods diamond,hmm,dbCANsub
GFF_ARG=""; [ -n "$GFF" ] && GFF_ARG="--gff $GFF"
"$VENV/bin/python" "$NF/bin/emit_baseline_funcscan.py" \
    --overview rundbcan/overview.tsv --faa "$FAA" --sample "$SAMPLE" --outdir "$OUT/funcscan" $GFF_ARG

echo "=== 2/8 advanced ESM-C embed + label-free infer ==="
"$VENV/bin/python" "$SCRIPTS/embed_esmc.py" --fasta "$FAA" \
    --out-prefix query.esmc --model esmc_600m --nshards 1 --shard 0
[ -f query.esmc.shard0.npz ] && mv -f query.esmc.shard0.npz query.esmc.npz || true
"$VENV/bin/python" "$SCRIPTS/infer_esmc.py" --emb query.esmc.npz \
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

if [ "$DO_INTERPRO" = 1 ]; then
  echo "=== 4b/8 InterPro domains + GO terms ==="
  # Lands in the funcscan tree so `bioforge-ingest` populates the Gene Ontology
  # card + InterPro-domains table on the gene page. Real InterProScan if configured,
  # else GO derived from the Pfam domains above (offline, self-contained).
  IPR_DIR="$OUT/funcscan/protein_annotation/interproscan"; mkdir -p "$IPR_DIR"
  IPR_OUT="$IPR_DIR/${SAMPLE}_interproscan_faa.tsv"
  P2IPR_ARG=""; [ -n "$PFAM2INTERPRO" ] && P2IPR_ARG="--pfam2interpro $PFAM2INTERPRO"
  if [ -n "$INTERPROSCAN_SH" ] && [ -x "$INTERPROSCAN_SH" ]; then
    echo "  real InterProScan: $INTERPROSCAN_SH"
    "$INTERPROSCAN_SH" -i "$FAA" -f tsv -goterms -pa -o "$IPR_OUT" -T "$OUTDIR/ips_tmp" \
      || { echo "  real InterProScan failed -> Pfam->GO fallback" >&2
           "$VENV/bin/python" "$NF/bin/pfam_to_interproscan.py" --domains-tsv "$FEAT/domains.tsv" \
               --pfam2go "$PFAM2GO" $P2IPR_ARG --faa "$FAA" --out "$IPR_OUT"; }
  else
    echo "  interproscan.sh not set -> deriving GO from Pfam domains (offline)"
    "$VENV/bin/python" "$NF/bin/pfam_to_interproscan.py" --domains-tsv "$FEAT/domains.tsv" \
        --pfam2go "$PFAM2GO" $P2IPR_ARG --faa "$FAA" --out "$IPR_OUT"
  fi
  echo "  interpro/GO -> $IPR_OUT"
else
  echo "=== 4b/8 InterPro/GO SKIPPED (--no-interpro) ==="
fi

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

if [ "$DO_STRUCTURE" = 1 ]; then
  echo "=== 7b/8 ESMFold 3D structure (local GPU) ==="
  # Folds every protein on the local GPU (facebook/esmfold_v1 via transformers).
  # No hosted service, no length cap beyond GPU memory — runs anywhere with CUDA,
  # so users can reproduce on their own server. Writes one PDB + a pLDDT manifest.
  "$VENV/bin/python" "$SCRIPTS/fold_esmfold.py" \
      --fasta "$FAA" --outdir esmfold_out --manifest esmfold_out/fold_manifest.tsv \
      --max-len 1100 --min-len 20
  # -> structures.tsv (+ structures/<id>.pdb, B-factors rescaled to 0-100 pLDDT)
  "$VENV/bin/python" "$NF/bin/structures_to_tsv.py" \
      --esmfold-dir esmfold_out --out-dir "$FEAT" --source ESMFold
else
  echo "=== 7b/8 ESMFold SKIPPED (--no-structure) ==="
fi

echo "=== 8/8 write v1.1 manifest ==="
STAGE="$OUTDIR/stage"; rm -rf "$STAGE"; mkdir -p "$STAGE"
cp "$PRED/"*.tsv "$STAGE/" 2>/dev/null || true
cp "$FEAT/"*.tsv "$STAGE/" 2>/dev/null || true
# structure PDBs must ride next to structures.tsv so the ingester resolves them
[ -d "$FEAT/structures" ] && cp -r "$FEAT/structures" "$STAGE/structures"
( cd "$STAGE" && "$VENV/bin/python" "$NF/bin/write_manifest.py" \
    --outdir "$OUT" --sample "$SAMPLE" --stage-dir . \
    --release-label "workup-$(date +%Y-%m-%d)" --pipeline-version '0.1.0' \
    --tool-versions '{"esm":"3.2.1","deeptmhmm":"1.0.24","clean":"maxsep","hmmscan":"3.4","biopython":"1.85"}' )
echo "manifest -> $OUT/cazyme_advanced/manifest.json"

if [ "$DO_AI_REPORT" = 1 ]; then
  echo "=== 8b/8 AI-ready per-protein JSON reports ==="
  # Assemble the evidence the report generator reads: the raw ESM-C head +
  # fusion calls (written in CWD=$OUTDIR), the per-protein feature TSVs, and the
  # manifest. Then emit one grounded "prompt pack" JSON per protein.
  AIA="$OUTDIR/ai_report_assets"; rm -rf "$AIA"; mkdir -p "$AIA"
  for f in raw_knn.tsv raw_centroid.tsv raw_contrastive.tsv fusion_raw.tsv; do
    [ -f "$OUTDIR/$f" ] && cp "$OUTDIR/$f" "$AIA/"
  done
  cp "$FEAT/"*.tsv "$AIA/" 2>/dev/null || true
  cp "$OUT/cazyme_advanced/manifest.json" "$AIA/" 2>/dev/null || true
  AIOUT="$OUT/cazyme_advanced/ai_reports"
  "$VENV/bin/python" "$SCRIPTS/build_ai_report.py" --assets "$AIA" --all --outdir "$AIOUT" \
    && echo "AI reports -> $AIOUT/<protein_id>_ai_report.json" \
    || echo "WARN: AI-report generation failed (non-fatal); continuing."
else
  echo "=== 8b/8 AI reports SKIPPED (--no-ai-report) ==="
fi

if [ "$SERVE" = 1 ]; then
  echo "=== ingest + serve ==="
  source "$BIODB_VENV/bin/activate"
  # Preflight: the vendored web layer must be installed in this venv.
  if ! command -v bioforge-ingest >/dev/null 2>&1; then
    echo "ERROR: 'bioforge-ingest' not found in \$BIODB_VENV ($BIODB_VENV)." >&2
    echo "       The web UI is vendored in this repo — install it once with:" >&2
    echo "         pip install -e \"$SELF_DIR\"" >&2
    echo "       or point BIODB_VENV at a venv that already has it." >&2
    exit 3
  fi
  DB="$OUTDIR/${SAMPLE}.db"; rm -f "$DB"; export DATABASE_URL="sqlite:///$DB"
  export BIOFORGE_TRACKS_DIR="$OUTDIR/web_static/tracks"
  ( cd "$BIODB_SRC" && alembic upgrade head )
  bioforge-ingest          "$OUT/funcscan"
  bioforge-ingest-advanced "$OUT/cazyme_advanced/manifest.json" \
      --structures-dir "$OUTDIR/web_static/structures"
  echo "### serving http://$HOST:$PORT  (Ctrl-C to stop)"
  exec uvicorn bioforge.api.main:app --host "$HOST" --port "$PORT"
else
  echo "### done. Ingest + serve with:  --serve   (or dbcan4 run --serve on the manifest)"
fi
