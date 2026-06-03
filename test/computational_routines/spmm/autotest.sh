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

# Compare serial and parallel outputs if present
compare_dirs "serial" "parallel" "${log_file_name}"

info "PSBLAS psb_spmm test completed."
