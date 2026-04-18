#!/usr/bin/env bash
#SBATCH --job-name=psb_swapdata_test
#SBATCH --partition=boost_usr_prod
#SBATCH --time=02:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --gpus-per-node=4
#SBATCH --export=ALL
#SBATCH -A CNHPC_1736213
#SBATCH --output=psb_comm_swapdata_%j.out
#SBATCH --error=psb_comm_swapdata_%j.err

set -euo pipefail

# Environment tuned like the main comm test script.
export UCX_VFS_ENABLE=n
export UCX_VFS_USE_FUSE=n
export UCX_STATS_DEST=none
export UCX_LOG_LEVEL=error
export UCX_TLS=dc,sm,self
export UCX_NET_DEVICES=mlx5_0:1
export OMPI_MCA_coll=^hcoll,han
export OMPI_MCA_coll_hcoll_enable=0
export UCX_MEMTYPE_CACHE=n
export UCX_CLOSE_TIMEOUT=10s

EXEC=./test/comm/swapdata/runs/psb_comm_test
DIM_START=${DIM_START:-20}
DIM_END=${DIM_END:-300}
DIM_STEP=${DIM_STEP:-20}
ITERS=${ITERS:-10}
MODE=${MODE:-both}
DEBUG=${DEBUG:-0}

EXTRA_ARGS=()
if [[ "$DEBUG" != "0" ]]; then
  EXTRA_ARGS+=(--debug)
fi

if [[ ! -x "$EXEC" ]]; then
  echo "Executable not found: $EXEC" >&2
  echo "Build the test first with 'make' in test/comm/swapdata/" >&2
  exit 1
fi

echo "Running swapdata comm test"
echo "  EXEC=$EXEC"
echo "  DIM_START=$DIM_START"
echo "  DIM_END=$DIM_END"
echo "  DIM_STEP=$DIM_STEP"
echo "  ITERS=$ITERS"
echo "  MODE=$MODE"
echo "  DEBUG=$DEBUG"

for DIM in $(seq "$DIM_START" "$DIM_STEP" "$DIM_END"); do
  echo ""
  echo "=== Running DIM=$DIM ==="
  srun --exclusive -N2 -n64 "$EXEC" --dim "$DIM" --iters "$ITERS" --mode "$MODE" "${EXTRA_ARGS[@]}"
done
