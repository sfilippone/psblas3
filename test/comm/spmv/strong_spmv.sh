#!/bin/bash
# ============================================================================
# strong_spmv.sbatch  STRONG scaling, CPU-only SpMV across comm backends
#
# Fixed total problem size; MPI ranks grow across scale points. The spmv test
# sweeps ALL comm backends internally (P2P, NEIGHBOR, PNEIGHBOR, MPI_GET,
# MPI_PUT) in a single run, so there is no per-mode loop.
# All scale points run inside ONE reserved allocation (same hardware).
# ============================================================================
#SBATCH --job-name=spmv_strong
#SBATCH --output=spmv_strong_%j.out
#SBATCH --error=spmv_strong_%j.err
#SBATCH --nodes=8                   # reserve the maximum up front
#SBATCH --ntasks-per-node=80        # 80 cores per ACC node -> 80 ranks/node
#SBATCH --cpus-per-task=1           # CPU-only: one core per rank
#SBATCH --time=00:30:00
#SBATCH --qos=acc_debug
#SBATCH --account=YOUR_ACCOUNT_HERE

# ============================================================================
# USER CONFIGURATION
# ============================================================================
DIM=280                             # FIXED problem size (idim^3 unknowns)
TIMES=50                            # SpMV repetitions (timed)
CPU_FMT=CSR                         # CPU matrix storage format
PROFILE=true                        # Score-P profiling on/off
RANK_POINTS="80 160 320 640"        # total ranks per scale point (multiples of 80)

EXE=$HOME/Desktop/scorep/psblas3/test/comm/spmv/runs/psb_spmv_kernel

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
export SCOREP_CUDA_ENABLE=no        # CPU-only run: silence Score-P CUDA adapter

RESDIR=$SLURM_SUBMIT_DIR/results_spmv_strong_${SLURM_JOB_ID}
mkdir -p $RESDIR

# --- Score-P region filter: exclude the high-frequency tiny USR functions so
#     the profile reflects real MPI/comm time, not instrumentation overhead. ---
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

echo "=== SPMV STRONG scaling (CPU-only) ==="
echo "  fixed_dim=$DIM times=$TIMES cpu_fmt=$CPU_FMT"
echo "  reserved_nodes=$SLURM_NNODES rank_points=[$RANK_POINTS] profile=$PROFILE"
echo "======================================"

for NRANKS in $RANK_POINTS; do
    NNODES=$(( (NRANKS + 79) / 80 ))     # 80 ranks per node
    STEP_DIR=$RESDIR/${NRANKS}ranks
    mkdir -p $STEP_DIR

    if [ "$PROFILE" = "true" ]; then
        export SCOREP_ENABLE_PROFILING=true
        export SCOREP_ENABLE_TRACING=false
        export SCOREP_TOTAL_MEMORY=128M
        export SCOREP_FILTERING_FILE=$FILTER
        export SCOREP_EXPERIMENT_DIRECTORY=$STEP_DIR/scorep_profile
    else
        export SCOREP_ENABLE_PROFILING=false
        export SCOREP_ENABLE_TRACING=false
    fi

    echo ""
    echo ">>> STRONG point: $NRANKS ranks ($NNODES nodes), fixed dim=$DIM"
    srun -N $NNODES -n $NRANKS --ntasks-per-node=80 --cpus-per-task=1 \
        $EXE --gpu=FALSE --dim=$DIM --times=$TIMES --cpu_fmt=$CPU_FMT \
        > $STEP_DIR/run.out 2>&1
    echo ">>> exit=$? output=$STEP_DIR/run.out"
done

echo ""
echo "=== SPMV STRONG DONE. Results: $RESDIR ==="
ls -R $RESDIR
