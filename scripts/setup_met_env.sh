#!/usr/bin/env bash
# Build the dbcan4-advanced GPU env on met at a PERSISTENT path (job workdir is ephemeral).
set -uo pipefail
ENVDIR=/array1/xinpeng/dbcan4-advanced/venv
export HF_HOME=/array1/xinpeng/dbcan4-advanced/hf_cache
mkdir -p /array1/xinpeng/dbcan4-advanced "$HF_HOME"

echo "=== egress check ==="
curl -sSI https://pypi.org/simple/ 2>&1 | head -1
curl -sSI https://huggingface.co 2>&1 | head -1

echo "=== install uv ==="
python3 -m pip install --user -q uv 2>&1 | tail -2 || pip install --user -q uv 2>&1 | tail -2
UV=$(command -v uv || echo "$HOME/.local/bin/uv")
"$UV" --version || { echo "UV_MISSING"; exit 2; }

echo "=== create venv (py3.11) ==="
"$UV" venv --python 3.11 "$ENVDIR" 2>&1 | tail -3
source "$ENVDIR/bin/activate"

echo "=== install torch (CUDA 12.1) + esm + retrieval stack ==="
"$UV" pip install --no-cache-dir "torch==2.4.1" --index-url https://download.pytorch.org/whl/cu121 2>&1 | tail -3
"$UV" pip install --no-cache-dir esm fair-esm scikit-learn faiss-cpu biopython pandas numpy h5py tqdm 2>&1 | tail -5

echo "=== smoke test: torch + cuda ==="
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), 'ngpu', torch.cuda.device_count())"
echo "=== smoke test: esm import ==="
python -c "from esm.models.esmc import ESMC; print('esm ESMC import OK')" 2>&1 | tail -3
echo "=== freeze ==="
"$UV" pip freeze > ./met_env_freeze.txt 2>&1
wc -l ./met_env_freeze.txt
echo "ENV_SETUP_DONE at $ENVDIR"
