#!/bin/bash
# ============================================================================
# extract_scorep.sh  -- auto-extract comm-scheme measurements from a
#                       results_*_<jobid> directory (swapdata/spmv style).
#
# For every <N>ranks/ subdir it collects:
#   1) run.out wall-clock avg time per backend  (the REAL headline number)
#   2) cube_stat -r per swap routine from scorep_profile/profile.cubex:
#        - INCL (total time in the routine, summed over call paths)
#        - one-time setup  (Win_create/Win_allocate/topology_init/
#                            ini_memory_buffer_layout/alltoallv_init)
#        - steady = INCL - one-time
#        - every MPI_* component (to scorep_components.csv)
#
# Output (written into RESULTS_DIR):
#   runout.csv            ranks,backend,total_s,avg_s
#   scorep_summary.csv    ranks,scheme,incl_s,onetime_s,steady_s
#   scorep_components.csv ranks,scheme,component,time_s
# plus pretty pivot tables on stdout.
#
# Usage:
#   ./extract_scorep.sh                 # process ALL <N>ranks/ in current dir
#   ./extract_scorep.sh 80ranks         # process ONE rank-point dir
#   ./extract_scorep.sh path/to/results # process ALL <N>ranks/ under that dir
#
# CSVs ACCUMULATE (header written only if absent), so you can run one
# rank-point at a time and the pivots grow. To start over: rm *.csv
# ============================================================================
set -u
ARG="${1:-.}"
OUTDIR="${OUTDIR:-.}"     # where the CSVs go (default: current dir)

command -v cube_stat >/dev/null 2>&1 || { echo "ERROR: cube_stat not in PATH"; exit 1; }

# scheme label : mangled routine name used by cube_stat -r
SCHEMES=(
  "baseline:psi_d_comm_v_mod.psi_d_swapdata_impl::psi_dswap_baseline_vect"
  "neighbor:psi_d_comm_v_mod.psi_d_swapdata_impl::psi_dswap_neighbor_topology_vect"
  "persistent:psi_d_comm_v_mod.psi_d_swapdata_impl::psi_dswap_neighbor_persistent_topology_vect"
  "rma_pull:psi_d_comm_v_mod.psi_d_swapdata_impl::psi_dswap_rma_pull_vect"
  "rma_push:psi_d_comm_v_mod.psi_d_swapdata_impl::psi_dswap_rma_push_vect"
)

RUNOUT="$OUTDIR/runout.csv"
SUMMARY="$OUTDIR/scorep_summary.csv"
COMPONENTS="$OUTDIR/scorep_components.csv"
# write header only if the file does not exist yet (append/accumulate mode)
[ -f "$RUNOUT" ]     || echo "ranks,backend,total_s,avg_s,setup_s"     > "$RUNOUT"
[ -f "$SUMMARY" ]    || echo "ranks,scheme,incl_s,setup_s,gpuwait_s,pack_s,mpimove_s,mpiwait_s,netmpi_s" > "$SUMMARY"
[ -f "$COMPONENTS" ] || echo "ranks,scheme,bucket,component,time_s"     > "$COMPONENTS"

# Select dirs: a single rank-point dir (has scorep_profile/run.out) or, if ARG
# is a results dir, all its <N>ranks/ subdirs.
if [ -f "$ARG/scorep_profile/profile.cubex" ] || [ -f "$ARG/run.out" ]; then
  DIRS=( "$ARG" )
else
  DIRS=( $(ls -d "$ARG"/*ranks 2>/dev/null | sort -t/ -k99 -V) )
fi
[ ${#DIRS[@]} -eq 0 ] && { echo "No rank-point dir(s) found under '$ARG'"; exit 1; }

for d in "${DIRS[@]}"; do
  ranks=$(basename "$d" | sed 's/[^0-9]//g')
  [ -z "$ranks" ] && continue

  # ---- 1) run.out wall-clock avg + setup per backend ----
  # One run.out holds several "comm backend" blocks in sequence; we accumulate
  # tot/avg/setup for the current block and flush it when the next block starts
  # (or at EOF). "setup time" is optional: 0 if the benchmark does not emit it
  # (i.e. until the warm-up-timing patch lands), so the break-even falls back to
  # pure steady-state comparison.
  if [ -f "$d/run.out" ]; then
    awk -v r="$ranks" '
      function flush(){ if (be!=""){ print r","be","tot","avg","setup } be="";tot=0;avg=0;setup=0 }
      /comm backend/  { flush(); split($0,a,":"); be=a[2]; gsub(/ /,"",be) }
      /total time/    { split($0,a,":"); tot=a[2]+0 }
      /avg time/      { split($0,a,":"); avg=a[2]+0 }
      /setup time/    { split($0,a,":"); setup=a[2]+0 }
      END { flush() }
    ' "$d/run.out" >> "$RUNOUT"
  fi

  # ---- 2) scorep cube_stat per scheme ----
  cubex="$d/scorep_profile/profile.cubex"
  [ -f "$cubex" ] || { echo "WARN: no profile.cubex in $d" >&2; continue; }

  for entry in "${SCHEMES[@]}"; do
    label="${entry%%:*}"; routine="${entry#*:}"
    out=$(cube_stat -r "$routine" "$cubex" 2>/dev/null)
    [ -z "$out" ] && continue
    printf '%s\n' "$out" | awk -v ranks="$ranks" -v scheme="$label" -v comp="$COMPONENTS" '
      # Bucketize every child component of the swap routine. cube_stat -r reports,
      # for each child, its inclusive time within the routine subtree, and the
      # children partition INCL -- so per-bucket sums are meaningful.
      #
      # Why this matters on GPU: the single largest component in EVERY scheme is
      # d_cuda_device_wait (~0.036s, the GPU stream sync), which is NOT network
      # time and is ~equal across schemes. The old "steady = INCL - onetime" left
      # it inside, masking the real differentiator -- the MPI wait/sync time.
      # Buckets (each component lands in exactly one, order matters):
      #   setup   : one-time  Win_create / topology_init / alltoallv_init /
      #                       ini_memory_buffer_layout
      #   gpuwait : GPU stream sync (device_wait, cuda_sync) -- exclude from MPI
      #   mpiwait : blocking MPI sync (MPI_Wait, Win_wait/post/start/complete)
      #             <-- the discriminating metric ("guarda le wait")
      #   mpimove : actual transfer (Isend/Irecv/Put/Get/Start/[I]neighbor_alltoallv)
      #   pack    : gather/scatter/buffer (gthz*, sctb*, new_buffer, realloc)
      #   netmpi  = mpimove + mpiwait  (honest per-iteration network cost)
      {
        k=split($0,a,","); v=a[k]+0          # value = last comma-field
        nm=$0; sub(/,[^,]*$/,"",nm)           # name  = everything before it
        if (nm ~ /^INCL\(/) { incl+=v; next }
        if (nm ~ /^EXCL\(/) { next }
        bkt="other"
        if      (nm ~ /Win_create|Win_allocate|topology_init|ini_memory_buffer_layout|alltoallv_init/) { setup+=v;   bkt="setup" }
        else if (nm ~ /device_wait|cuda_sync/)                                                          { gpuwait+=v; bkt="gpuwait" }
        else if (nm ~ /MPI_Wait|Win_wait|Win_post|Win_start|Win_complete/)                              { mpiwait+=v; bkt="mpiwait" }
        else if (nm ~ /MPI_Isend|MPI_Irecv|MPI_Send|MPI_Recv|MPI_Put|MPI_Get|MPI_Start|neighbor_alltoallv|Neighbor_alltoallv/) { mpimove+=v; bkt="mpimove" }
        else if (nm ~ /gthz|sctb|new_buffer|realloc/)                                                   { pack+=v;    bkt="pack" }
        cn=nm; sub(/.*::/,"",cn); gsub(/"/,"",cn)
        print ranks","scheme","bkt","cn","v >> comp
      }
      END { printf "%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", \
                   ranks, scheme, incl, setup, gpuwait, pack, mpimove, mpiwait, mpimove+mpiwait }
    ' >> "$SUMMARY"
  done
done

# ---- pretty pivots ----
pivot() {  # $1=csv  $2=value-col-index(1-based, after ranks,key)  $3=title
  awk -F, -v vc="$2" -v title="$3" '
    NR==1 { next }
    { r=$1; k=$2; val=$(2+vc); v[r","k]=val; rs[r]=1; ks[k]=1 }
    END {
      n=0; for (k in ks) cols[++n]=k
      # stable-ish column order
      order="baseline neighbor persistent rma_pull rma_push P2P NEIGHBOR PNEIGHBOR MPI_GET MPI_PUT"
      m=split(order,oc," "); delete cols; n=0
      for (i=1;i<=m;i++) if (oc[i] in ks) cols[++n]=oc[i]
      printf "\n=== %s ===\n", title
      printf "%-8s", "ranks"; for (i=1;i<=n;i++) printf "%14s", cols[i]; printf "\n"
      nr=0; for (r in rs) rr[++nr]=r
      # numeric sort of ranks
      for (i=1;i<=nr;i++) for (j=i+1;j<=nr;j++) if (rr[i]+0>rr[j]+0){t=rr[i];rr[i]=rr[j];rr[j]=t}
      for (i=1;i<=nr;i++){ printf "%-8s", rr[i]; for (c=1;c<=n;c++){key=rr[i]","cols[c]; printf "%14s", (key in v)?v[key]:"-"} printf "\n" }
    }' "$1"
}

# ---- break-even pivot: iterations after which scheme X beats P2P ----
# Model for an iterative solver doing n_it spmv (halo exchanges):
#     T_X(n) = setup_X + n * steady_X     (steady = run.out avg, setup excluded
#                                           by the warm-up; setup = run.out "setup time")
# X beats the P2P reference when T_X(n) < T_P2P(n). Solving for the crossover:
#     n* = (setup_X - setup_P2P) / (steady_P2P - steady_X)
# Output legend (per cell):
#   <number> : break-even n*  (X wins for n_it > n*)   -- more setup, less steady
#   always   : X wins from iteration 1 (cheaper setup AND cheaper steady) -> dominates P2P
#   never    : P2P wins for all n_it (more setup AND slower steady) -> X is dominated
#   <n*      : X wins ONLY for n_it < n*  (cheaper setup but slower steady)
#              ^ this is the GPU/few-ranks RMA case: if RMA steady drops below P2P
#                you'll instead see a number or "always" here.
breakeven() {
  awk -F, '
    NR==1 { next }
    { r=$1; be=$2; steady[r","be]=$4+0; setp[r","be]=$5+0; rs[r]=1; bes[be]=1 }
    END {
      order="NEIGHBOR PNEIGHBOR MPI_GET MPI_PUT"; m=split(order,oc," ")
      printf "\n=== break-even n* vs P2P  (iterations to amortize setup; needs run.out '\''setup time'\'') ===\n"
      printf "%-8s", "ranks"; for (i=1;i<=m;i++) if (oc[i] in bes) printf "%14s", oc[i]; printf "\n"
      nr=0; for (r in rs) rr[++nr]=r
      for (i=1;i<=nr;i++) for (j=i+1;j<=nr;j++) if (rr[i]+0>rr[j]+0){t=rr[i];rr[i]=rr[j];rr[j]=t}
      anysetup=0
      for (i=1;i<=nr;i++){
        r=rr[i]; printf "%-8s", r
        sref=steady[r",P2P"]; pref=setp[r",P2P"]
        for (c=1;c<=m;c++){ be=oc[c]; if (!(be in bes)) continue
          key=r","be
          if (!(key in steady) || !(r",P2P" in steady)) { printf "%14s","-"; continue }
          dsteady = sref - steady[key]      # >0 : X faster per-iter than P2P
          dsetup  = setp[key] - pref        # >0 : X pays more setup than P2P
          if (setp[key]>0 || pref>0) anysetup=1
          if (dsteady > 0)      printf "%14s", (dsetup<=0 ? "always" : sprintf("%.0f", dsetup/dsteady))
          else if (dsteady < 0) printf "%14s", (dsetup>=0 ? "never"  : sprintf("<%.0f", dsetup/dsteady))
          else                  printf "%14s", (dsetup<0  ? "always" : "never")
        }
        printf "\n"
      }
      if (!anysetup) print "\n  NOTE: all setup_s == 0 -> run.out has no \"setup time\" yet (apply the\n        warm-up-timing patch to the benchmark). Columns reflect steady-state only."
    }' "$RUNOUT"
}

echo "Wrote: $RUNOUT  $SUMMARY  $COMPONENTS"
pivot "$RUNOUT"  2 "run.out  AVG time per backend [s]  (real wall-clock = steady-state per iter)"
pivot "$RUNOUT"  3 "run.out  SETUP time per backend [s]  (one-time, paid in warm-up)"
breakeven
pivot "$SUMMARY" 6 "scorep  MPI WAIT per scheme [s]  <== THE discriminator (sync/imbalance; 'guarda le wait')"
pivot "$SUMMARY" 7 "scorep  NET MPI per scheme [s]  (mpimove + mpiwait; excl GPU sync, pack, one-time setup)"
pivot "$SUMMARY" 5 "scorep  MPI MOVE per scheme [s]  (pure transfer: Isend/Irecv/Put/Get/a2av)"
pivot "$SUMMARY" 2 "scorep  SETUP per scheme [s]  (one-time: Win_create/topology_init/alltoallv_init/buf_layout)"
pivot "$SUMMARY" 3 "scorep  GPU WAIT per scheme [s]  (d_cuda_device_wait; sanity: ~equal across schemes)"
pivot "$SUMMARY" 1 "scorep  INCL per scheme [s]  (everything; dominated by GPU wait -> misleading alone)"

echo
echo "Tip: full per-component breakdown (now bucketed) is in $COMPONENTS"
echo "     (ranks,scheme,bucket,component,time).  e.g.:"
echo "       grep ',rma_pull,' $COMPONENTS | sort -t, -k5 -gr   # rma_pull hot spots"
echo "       awk -F, '\$3==\"mpiwait\"' $COMPONENTS              # only the waits"
