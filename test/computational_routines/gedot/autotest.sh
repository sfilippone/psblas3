#!/bin/bash

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/../common/testlib.sh"

# Variables definition
dir1="serial"
dir2="parallel"
num_procs=$(get_num_procs)
log_dir="logs"
log_file_name="${log_dir}/psblas_gedot_test.log"

ensure_dirs runs serial parallel vectors "${log_dir}"
log_header "${log_file_name}" "psb_gedot"


# Check if the executable ELF file exists
if [ ! -f "./runs/psb_gedot_test" ]; then
    warn "Executable not found. Running make..."
    make
    if [ ! -f "./runs/psb_gedot_test" ]; then
        err "Failed to create executable. Check make command."
    fi
else
    info "The executable already exists. Skipping the make process."
fi


# Excecute tests and save results
info "Running the PSBLAS psb_gedot test..."
echo "" >> "${log_file_name}"
run_mpi 1 ./runs/psb_gedot_test "single process computation"
run_mpi "$num_procs" ./runs/psb_gedot_test "${num_procs} processes computation"


# Iterate through files in the first directory
compare_dirs "$dir1" "$dir2" "${log_file_name}"

info "PSBLAS psb_gedot test successfully completed."