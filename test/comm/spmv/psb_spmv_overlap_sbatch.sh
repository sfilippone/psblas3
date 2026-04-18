#!/usr/bin/env bash
#SBATCH --job-name=psb_spmv_overlap
#SBATCH --partition=boost_usr_prod
#SBATCH --time=01:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --gpus-per-node=4
#SBATCH --export=ALL
#SBATCH -A CNHPC_1736213
#SBATCH --output=psb_spmv_overlap_%j.out
#SBATCH --error=psb_spmv_overlap_%j.err

set -euo pipefail

# Environment tuned like the existing comm test script.
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

EXEC=./test/comm/spmv/runs/spmv_overlap
IDIM_LIST=${IDIM_LIST:-"20 40 60 80 100 140 180 220 260 300"}
TIMES=${TIMES:-100}
BUILD_IF_MISSING=${BUILD_IF_MISSING:-1}

if [[ ! -x "$EXEC" ]]; then
  echo "Executable not found: $EXEC" >&2
  exit 1
fi

echo "Running SpMV overlap comm test"
echo "  EXEC=$EXEC"
echo "  IDIM_LIST=$IDIM_LIST"
echo "  TIMES=$TIMES"
echo "  BUILD_IF_MISSING=$BUILD_IF_MISSING"

for IDIM in $IDIM_LIST; do
  echo ""
  echo "=== Running IDIM=$IDIM TIMES=$TIMES ==="
  IDIM="$IDIM" TIMES="$TIMES" srun --exclusive -N2 -n64 "$EXEC"
done
