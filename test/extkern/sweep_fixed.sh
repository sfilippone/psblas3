#!/bin/bash
# Serial sweep after the DIA/HDIA fixes: CSR, HLL (six hack sizes), HDIA.
set -u
HERE=/home/stack/Desktop/PSBLAS/hll/psblas3/test/extkern
BIN=$HERE/runs/dpdegenmv
OUT=${OUT:-$HERE/results-sweep}
CORE=${CORE:-6}
IDIMS=${IDIMS:-"120 160 200 240 280 320 360 400 440 480 500"}
SPECS=${SPECS:-"CSR:32 HLL:4 HLL:8 HLL:16 HLL:32 HLL:64 HLL:128 HDIA:32"}
mkdir -p "$OUT"
# one sweep at a time: two instances would share the pinned core and spoil the timings
exec 9>"$OUT/.lock"
flock -n 9 || { echo "another sweep is already running, aborting"; exit 1; }
CSV=$OUT/sweep.csv
[ -f "$CSV" ] || echo "fmt,hks,idim,nrows,nnz,mem,ntests,time,mflops" > "$CSV"
for d in $IDIMS; do
  # memory guard: the biggest sizes need tens of GB, skip if the machine is busy
  avail=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)
  need=$(( d*d*d/1000000*3/10 + 2 ))
  if [ "$avail" -lt "$((need+8))" ]; then
    echo "[$(date +%H:%M:%S)] idim=$d skipped: ${avail}GB available, ~${need}GB needed"; continue
  fi
  for spec in $SPECS; do
    f=${spec%:*}; h=${spec#*:}
    key=$f,$( [ "$f" = HLL ] && echo "$h" || echo 0 ),$d,
    grep -q "^$key" "$CSV" && continue
    out=$OUT/run.$f.$h.$d.txt
    printf "%s\n%s\n%s\nF\n" "$f" "$d" "$h" | numactl -C "$CORE" -m 0 "$BIN" > "$out" 2>&1
    awk -v f="$f" -v h="$( [ "$f" = HLL ] && echo "$h" || echo 0 )" -v d="$d" \
      '/Size of matrix/{nr=$NF} /Number of nonzeros/{nz=$NF} /Memory occupation/{m=$NF} \
       /Time for/{nt=$3;t=$NF} /^MFLOPS/{printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n",f,h,d,nr,nz,m,nt,t,$NF}' \
      "$out" >> "$CSV"
    echo "[$(date +%H:%M:%S)] $f h=$h idim=$d  ->  $(tail -1 "$CSV")"
  done
done
echo "sweep finished"
