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
  local eps_mode
  local cmp_mode
  local n_gamma
  local u_round
  local total_files=0
  local total_fail=0
  local total_diff=0

  # EPS_MODE selects how the tolerance value is computed (absolute | gamma_n);
  # COMPARE selects how it is applied (absolute | relative). They are
  # orthogonal: e.g. gamma_n + relative gives the standard finite-precision
  # error bound |fl - exact| <= gamma_n * |value| used for spmm.
  eps_mode=${PSBLAS_TEST_EPS_MODE:-absolute}
  cmp_mode=${PSBLAS_TEST_COMPARE:-absolute}
  eps=${PSBLAS_TEST_TOL:-1e-6}

  if [[ "$eps_mode" == "gamma_n" ]]; then
    n_gamma=${PSBLAS_TEST_N:-0}
    u_round=${PSBLAS_TEST_UNIT_ROUNDOFF:-1.19e-7}
    if [[ "$n_gamma" -gt 0 ]]; then
      eps=$(awk -v n="$n_gamma" -v u="$u_round" 'BEGIN { g = (n*u)/(1.0 - n*u); printf "%.12g", g }')
    fi
  fi

  if [[ ! -d "$dir1" || ! -d "$dir2" ]]; then
    warn "Missing directories for comparison (${dir1}, ${dir2})."
    return 0
  fi

  for file1 in "$dir1"/*; do
    if [[ ! -f "$file1" ]]; then
      continue
    fi
    total_files=$((total_files + 1))
    local filename
    filename=$(basename "$file1")
    local file2="$dir2/$filename"
    if [[ -f "$file2" ]]; then
      local diff_count
      if [[ "${filename}" == *.mtx ]]; then
        diff_count=$(awk -v f1="$file1" -v f2="$file2" -v eps="$eps" -v mode="$cmp_mode" '
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
              tol = eps
              if (mode == "relative") {
                # mixed relative/absolute tolerance: scale eps by the
                # magnitude of the operands (floor of 1 keeps it absolute
                # near zero). Required for finite-precision results whose
                # summation order differs between serial and parallel runs.
                m = (a[i] < 0 ? -a[i] : a[i])
                p = (b[i] < 0 ? -b[i] : b[i])
                scale = (m > p ? m : p)
                if (scale < 1) scale = 1
                tol = eps * scale
              }
              if (d > tol) diff++
            }
            if (n1 != n2) diff += (n1 > n2 ? n1 - n2 : n2 - n1)
            print diff
          }')
      else
        diff_count=$(diff -U 0 "$file1" "$file2" | grep -E '^[+-]' | grep -v '^[+-]{3}' | wc -l)
      fi
      echo "[DIFF] ${file1} vs ${file2}: ${diff_count} differences" >> "$logfile"
      total_diff=$((total_diff + diff_count))
      if [[ "$diff_count" -gt 0 ]]; then
        total_fail=$((total_fail + 1))
      fi
    else
      err "File ${filename} does not exist in ${dir2}"
      total_fail=$((total_fail + 1))
    fi
  done

  if [[ "$total_files" -eq 0 ]]; then
    warn "No files to compare in ${dir1}."
    return 0
  fi

  if [[ "$total_fail" -eq 0 ]]; then
    echo -e "${GREEN}[PASS]\t  ${dir1} vs ${dir2}: ${total_files}/${total_files} tests passed${RESET}"
    return 0
  elif [[ "$total_fail" -eq "$total_files" ]]; then
    echo -e "${RED}[FAIL]\t  ${dir1} vs ${dir2}: 0/${total_files} tests passed${RESET}"
    return 2
  else
    local passed=$((total_files - total_fail))
    echo -e "${YELLOW}[WARN]\t  ${dir1} vs ${dir2}: ${passed}/${total_files} tests passed, ${total_fail} failed${RESET}"
    return 1
  fi
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
