#!/usr/bin/env bash
# Reproducible GPU setup for AlphaFold2-Multimer (ColabFold) in WSL.
# Fixes the "installer ships CPU-only JAX" problem: localColabFold installs a
# CPU jaxlib, so the GPU is never used until a CUDA build is added.
#
# Run inside WSL as root:  wsl -u root bash scripts/setup_gpu.sh
# Idempotent: skips steps already done; ends by verifying the GPU is visible.
set -e

CF_DIR=/root/localcolabfold
PY=$CF_DIR/colabfold-conda/bin/python
PIP=$CF_DIR/colabfold-conda/bin/pip

# 1. Install localColabFold if absent (uses MMseqs2 web API for MSAs, so no
#    multi-hundred-GB sequence database is needed locally).
if [ ! -x "$CF_DIR/colabfold-conda/bin/colabfold_batch" ]; then
  echo "[setup] installing localColabFold ..."
  cd /root
  wget -q https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
  bash install_colabbatch_linux.sh
fi

# 2. The critical fix: replace CPU jax with a CUDA build matched to the
#    installed jax version. jax[cuda12] pulls its own CUDA libs via pip, so no
#    system CUDA toolkit is required beyond the NVIDIA driver + WSL passthrough.
JAXVER=$("$PY" -c "import jax;print(jax.__version__)" 2>/dev/null || echo "")
if "$PY" -c "import jax,sys;sys.exit(0 if any('cuda' in str(d).lower() for d in jax.devices()) else 1)" 2>/dev/null; then
  echo "[setup] CUDA JAX already active."
else
  echo "[setup] installing CUDA JAX (jax[cuda12]==$JAXVER) ..."
  "$PIP" install --upgrade "jax[cuda12]==${JAXVER}"
fi

# 3. Verify.
echo "[setup] verifying GPU visibility:"
"$PY" -c "import jax; d=jax.devices(); print('  devices:', d); \
print('  GPU_OK', any('cuda' in str(x).lower() for x in d))"
