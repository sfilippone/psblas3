#!/bin/bash

# Common helpers for computational_routines autotest.sh scripts.

# Color codes
GREEN="\033[0;32m"
RED="\033[0;31m"
BLUE="\033[0;34m"
YELLOW="\033[33m"
RESET="\033[0m"

info() {
  echo -e "${BLUE}[INFO]\t  $*${RESET}"
}

warn() {
  echo -e "${YELLOW}[WARNING] $*${RESET}"
}

err() {
  echo -e "${RED}[ERROR]\t$*${RESET}"
}

ensure_dirs() {
  for d in "$@"; do
    mkdir -p "$d"
  done
}

log_header() {
  local logfile="$1"
  local title="$2"
  {
    echo "[RUN] ${title}"
    echo "[DATE] $(date +"%Y-%m-%d %H:%M:%S")"
  } >> "$logfile"
}

run_mpi() {
  local np="$1"
  local exe="$2"
  local label="$3"
  info "Starting ${label}"
  mpirun -np "$np" "$exe"
}

compare_dirs() {
  local dir1="$1"
  local dir2="$2"
  local logfile="$3"
  local eps

  eps=${PSBLAS_TEST_TOL:-1e-6}

  if [[ ! -d "$dir1" || ! -d "$dir2" ]]; then
    warn "Missing directories for comparison (${dir1}, ${dir2})."
    return 0
  fi

  for file1 in "$dir1"/*; do
    if [[ ! -f "$file1" ]]; then
      continue
    fi
    local filename
    filename=$(basename "$file1")
    local file2="$dir2/$filename"
    if [[ -f "$file2" ]]; then
      local diff_count
      if [[ "${filename}" == *.mtx ]]; then
        diff_count=$(awk -v f1="$file1" -v f2="$file2" -v eps="$eps" '
          function readvals(fname, vals,   line, n, header_seen) {
            n = 0
            while ((getline line < fname) > 0) {
              if (line ~ /^%/) continue
              if (line ~ /^[[:space:]]*$/) continue
              if (!header_seen) { header_seen = 1; continue }
              n++; vals[n] = line + 0
            }
            close(fname)
            return n
          }
          BEGIN {
            n1 = readvals(f1, a)
            n2 = readvals(f2, b)
            n = (n1 < n2 ? n1 : n2)
            diff = 0
            for (i = 1; i <= n; i++) {
              d = a[i] - b[i]
              if (d < 0) d = -d
              if (d > eps) diff++
            }
            if (n1 != n2) diff += (n1 > n2 ? n1 - n2 : n2 - n1)
            print diff
          }')
      else
        diff_count=$(diff -U 0 "$file1" "$file2" | grep -E '^[+-]' | grep -v '^[+-]{3}' | wc -l)
      fi
      echo "[DIFF] ${file1} vs ${file2}: ${diff_count} differences" >> "$logfile"
    else
      err "File ${filename} does not exist in ${dir2}"
    fi
  done
}

get_num_procs() {
  local detected
  local max
  detected=${PSBLAS_TEST_NP:-$(nproc)}
  max=${PSBLAS_TEST_MAX_NP:-4}

  if [[ "$detected" -gt "$max" ]]; then
    detected="$max"
  fi

  echo "$detected"
}
