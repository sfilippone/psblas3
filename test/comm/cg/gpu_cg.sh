#!/bin/bash
# ============================================================================
# gpu_cg.sh   GPU Conjugate Gradient across comm backends -- "sparse ranks"
#
# Final test of the comm-scheme study: a full CG solve (not just spmv) on GPU,
# with few MPI ranks (1 rank = 1 GPU) spread over MANY nodes so the halo
# exchange is dominated by inter-node traffic -- the regime where RMA may pay
# off. The test sweeps ALL backends (P2P, NEIGHBOR, PNEIGHBOR, MPI_GET,
# MPI_PUT) x 3 preconditioners internally, so there is no per-mode loop here.
#
# It also self-verifies, for every scheme, that an internally-allocated work
# vector (like CG's hidden r/p/q/z) inherits the scheme from desc_a%comm_type
# -- look for "[OK] internal-style work vector inherited scheme" in run.out,
# and abort on "INTERNAL-VECTOR SCHEME MISMATCH".
# ============================================================================
#SBATCH --job-name=cg_gpu
#SBATCH --output=cg_gpu_%j.out
#SBATCH --error=cg_gpu_%j.err
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4         # MN5-ACC: 4x H100 per node
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=20          # 80 cores / 4 ranks
#SBATCH --time=00:30:00
#SBATCH --qos=acc_debug
#SBATCH --account=YOUR_ACCOUNT_HERE

# ============================================================================
# USER CONFIGURATION
# ============================================================================
NREP=8                              # CG solves per (scheme,prec) for statistics
NWARM=1                             # warm-up solves (discarded)
ITMAX=500                           # CG iteration cap
IDIM=200                            # ignored when --matrix is given

# --- "sparse" knob: GPUs (= ranks) per node.
#     RANKS_PER_NODE=1 -> 1 GPU/node, maximally spread, pure inter-node comm.
RANKS_PER_NODE=1
RANK_POINTS="2 4 8"                 # total ranks per scale point

EXE=$HOME/Desktop/scorep/psblas3/test/comm/cg/runs/psb_comm_cg_test
MATRIX=$HOME/Desktop/scorep/psblas3/test/comm/data/Geo_1438.mtx

# ============================================================================
# ENVIRONMENT
# ============================================================================
module purge
module load bsc/1.0
module load nvidia-hpc-sdk/25.3
module load gcc/12.3.0
module load ucx/1.16.0-gcc
module load openmpi/5.0.5-gcc
module load openblas/0.3.27-gcc
export PATH="$HOME/scorep-cuda/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/scorep-cuda/lib:$LD_LIBRARY_PATH"

export OMPI_MCA_coll_hcoll_enable=0
export SCOREP_CUDA_ENABLE=yes
export SCOREP_CUDA_BUFFER=48M

RESDIR=$SLURM_SUBMIT_DIR/results_cg_gpu_${SLURM_JOB_ID}
mkdir -p $RESDIR

echo "=== CG GPU (sparse-ranks study) ==="
echo "  nrep=$NREP nwarm=$NWARM itmax=$ITMAX ranks_per_node=$RANKS_PER_NODE"
echo "  reserved_nodes=$SLURM_NNODES rank_points=[$RANK_POINTS]"
echo "==================================="

for NRANKS in $RANK_POINTS; do
    NNODES=$(( (NRANKS + RANKS_PER_NODE - 1) / RANKS_PER_NODE ))
    STEP_DIR=$RESDIR/${NRANKS}ranks
    mkdir -p $STEP_DIR

    echo ""
    echo ">>> CG GPU point: $NRANKS ranks on $NNODES nodes ($RANKS_PER_NODE GPU/node)"

    srun -N $NNODES -n $NRANKS --ntasks-per-node=$RANKS_PER_NODE \
         --gres=gpu:$RANKS_PER_NODE --gpu-bind=single:1 --cpus-per-task=20 \
    $EXE $IDIM $NREP $NWARM $ITMAX --gpu=TRUE --matrix=$MATRIX --fmt=MM \
    > $STEP_DIR/run.out 2>&1

    echo ">>> exit=$? output=$STEP_DIR/run.out"
    # quick propagation sanity at submit time:
    grep -q "INTERNAL-VECTOR SCHEME MISMATCH" $STEP_DIR/run.out \
        && echo "    !! SCHEME PROPAGATION FAILED at $NRANKS ranks -- check run.out"
done

echo ""
echo "=== CG GPU DONE. Results: $RESDIR ==="
ls -R $RESDIR
