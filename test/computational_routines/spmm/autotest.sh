#!/bin/bash

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/../common/testlib.sh"

# Variables definition
num_procs=$(get_num_procs)
log_dir="logs"
log_file_name="${log_dir}/psblas_spmm_test.log"

ensure_dirs runs serial parallel vectors matrix "${log_dir}"
log_header "${log_file_name}" "psb_spmm"

# Ensure required directories exist
# Ensure a default matrix is available
if [ ! -f "matrix/1138_bus.mtx" ]; then
  warn "No default matrix found. Place a MatrixMarket file at matrix/1138_bus.mtx"
fi

# Check if the executable ELF file exists
if [ ! -f "./runs/psb_spmm_test" ]; then
  warn "Executable not found. Running make..."
  make
  if [ ! -f "./runs/psb_spmm_test" ]; then
    err "Failed to create executable. Check make command."
    exit 1
  fi
else
  info "The executable already exists. Skipping the make process."
fi

# Execute tests and save results
info "Running the PSBLAS psb_spmm test..."
echo "" >> "${log_file_name}"
run_mpi 1 ./runs/psb_spmm_test "single process computation"
run_mpi "$num_procs" ./runs/psb_spmm_test "${num_procs} processes computation"

# Compare serial and parallel outputs if present.
# spmm results are single precision and the distributed sparse mat-vec sums the
# row contributions in a different order than the serial run, so the outputs only
# agree to single-precision roundoff. The proper bound is the finite-precision
# error gamma_n = n*u/(1-n*u): EPS_MODE=gamma_n computes it (u = single-precision
# unit roundoff 1.19e-7), COMPARE=relative applies it as |fl-exact| <= gamma_n*|value|
# (an absolute tolerance is unsatisfiable for values of magnitude ~1e2).
# N=1000 generously bounds the per-row accumulation length (max 11 nnz/row, ~22
# for symmetric storage) to absorb cancellation in y=alpha*A*x+beta*y and the
# scatter/gather roundoff, while staying tight enough to catch real errors
# (observed max relative diff ~1.9e-5, ~6x below this gamma_n ~1.2e-4).
PSBLAS_TEST_EPS_MODE=gamma_n PSBLAS_TEST_N=1000 PSBLAS_TEST_UNIT_ROUNDOFF=1.19e-7 PSBLAS_TEST_COMPARE=relative \
    compare_dirs "serial" "parallel" "${log_file_name}"
compare_status=$?

info "PSBLAS psb_spmm test completed."
exit $compare_status
