#!/bin/bash
# ============================================================================
# gpu_spmv.sh   GPU SpMV across comm backends -- "sparse ranks" study
#
# Goal: few MPI ranks spread over MANY nodes (1 rank = 1 GPU), so that the
# halo exchange is dominated by INTER-node traffic. This is where RMA
# (MPI_Get/MPI_Put) may finally pay off: with few, fat messages over the
# network the one-time window-creation cost can amortize (break-even n*),
# unlike the dense CPU runs where RMA was dominated.
#
# The spmv test sweeps ALL backends internally (P2P, NEIGHBOR, PNEIGHBOR,
# MPI_GET, MPI_PUT) in one run -> no per-mode loop here.
# ============================================================================
#SBATCH --job-name=spmv_gpu
#SBATCH --output=spmv_gpu_%j.out
#SBATCH --error=spmv_gpu_%j.err
#SBATCH --nodes=8                   # reserve up front; we use a subset per point
#SBATCH --ntasks-per-node=4         # MN5-ACC has 4x H100 per node
#SBATCH --gres=gpu:4                # request all 4 GPUs on each reserved node
#SBATCH --cpus-per-task=20          # 80 cores / 4 ranks = 20 cores per rank
#SBATCH --time=00:30:00
#SBATCH --qos=acc_debug
#SBATCH --account=YOUR_ACCOUNT_HERE

# ============================================================================
# USER CONFIGURATION
# ============================================================================
TIMES=50                            # SpMV repetitions (timed)
GPU_FMT=HLG                         # GPU matrix format (HLG|ELG|CSRG|HDIAG)
PROFILE=true                        # Score-P profiling on/off

# --- "sparse" knob: how many ranks (= GPUs) to put on each node.
#     RANKS_PER_NODE=1  -> 1 GPU/node, maximally spread, pure inter-node comm
#                          (this is the configuration that isolates RMA-over-network)
#     RANKS_PER_NODE=4  -> pack all 4 GPUs/node (denser, more intra-node)
RANKS_PER_NODE=1

# total ranks per scale point. With RANKS_PER_NODE=1 these are also node counts,
# so "2 4 8" means 2, 4, 8 nodes each holding a single GPU.
RANK_POINTS="2 4 8"

EXE=$HOME/Desktop/scorep/psblas3/test/comm/spmv/runs/psb_spmv_kernel
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
export SCOREP_CUDA_ENABLE=yes           # GPU run: turn the CUDA adapter ON
export SCOREP_CUDA_BUFFER=48M           # per-stream CUDA record buffer

RESDIR=$SLURM_SUBMIT_DIR/results_spmv_gpu_${SLURM_JOB_ID}
mkdir -p $RESDIR

# --- Score-P region filter (same as CPU: drop tiny high-frequency USR funcs) ---
FILTER=$RESDIR/scorep.filt
cat > $FILTER <<'EOF'
SCOREP_REGION_NAMES_BEGIN
  EXCLUDE
    psb_indx_map_mod::*
    psb_desc_mod::*
    psb_error_mod::*
    psb_gen_block_map_mod::*
    psi_penv_mod::*
    psb_hash_mod::*
SCOREP_REGION_NAMES_END
EOF

echo "=== SPMV GPU (sparse-ranks study) ==="
echo "  times=$TIMES gpu_fmt=$GPU_FMT ranks_per_node=$RANKS_PER_NODE"
echo "  reserved_nodes=$SLURM_NNODES rank_points=[$RANK_POINTS] profile=$PROFILE"
echo "====================================="

for NRANKS in $RANK_POINTS; do
    NNODES=$(( (NRANKS + RANKS_PER_NODE - 1) / RANKS_PER_NODE ))
    STEP_DIR=$RESDIR/${NRANKS}ranks
    mkdir -p $STEP_DIR

    if [ "$PROFILE" = "true" ]; then
        export SCOREP_ENABLE_PROFILING=true
        export SCOREP_ENABLE_TRACING=false
        export SCOREP_TOTAL_MEMORY=256M     # bumped: CUDA adapter needs more
        export SCOREP_FILTERING_FILE=$FILTER
        export SCOREP_EXPERIMENT_DIRECTORY=$STEP_DIR/scorep_profile
    else
        export SCOREP_ENABLE_PROFILING=false
        export SCOREP_ENABLE_TRACING=false
    fi

    echo ""
    echo ">>> GPU point: $NRANKS ranks on $NNODES nodes ($RANKS_PER_NODE GPU/node)"

    srun -N $NNODES -n $NRANKS --ntasks-per-node=$RANKS_PER_NODE \
         --gres=gpu:$RANKS_PER_NODE --gpu-bind=single:1 --cpus-per-task=20 \
    $EXE --gpu=TRUE --gpu_fmt=$GPU_FMT --matrix=$MATRIX --fmt=MM --times=$TIMES \
    > $STEP_DIR/run.out 2>&1

    echo ">>> exit=$? output=$STEP_DIR/run.out"
done

echo ""
echo "=== SPMV GPU DONE. Results: $RESDIR ==="
ls -R $RESDIR
