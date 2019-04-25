#!/bin/bash

# https://sparse.tamu.edu/RB/CPM/cz148.tar.gz
# input is filename.rb

rb=$1
dir=${rb%.rb}
inp=${dir}.inp
mkdir ${dir}
mv ${rb} ${dir}
cd    ${dir}

# --- write out input file ---
echo "11                  Number of inputs"   >  ${inp}
echo "${rb}"                                  >> ${inp}
echo "NONE"                                   >> ${inp}
echo "HB                  HB: Harwell-Boeing" >> ${inp}
echo "BiCGSTAB            Iterative method"   >> ${inp}
echo "BJAC                Preconditioner "    >> ${inp}
echo "CSR                 Storage format "    >> ${inp}
echo "BLOCK               Partition method"   >> ${inp}
echo "2                   ISTOPC"             >> ${inp}
echo "00500               ITMAX"              >> ${inp}
echo "-1                  ITRACE"             >> ${inp}
echo "002                 IRST "              >> ${inp}
echo "1.d-7               EPS"                >> ${inp}
