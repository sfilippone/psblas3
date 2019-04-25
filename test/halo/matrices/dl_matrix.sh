#!/bin/bash

# https://sparse.tamu.edu/RB/CPM/cz148.tar.gz

if [ `expr match "$1" "https://sparse.tamu.edu/"` == 24 ]; then
    f=$(basename $1)
    dir=${f%.tar.gz}
    rb=${dir}.rb
    inp=${dir}.inp
    wget $1
    tar zxf ${f}
    rm ${f}
    cd ${dir}/
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
else
    echo "Please submit a .tar.gz file from SuiteSparse Matrix Collection"
fi
#                      123456789012345678901234


# declare -a matrix_size=( "148" "308" "628" "1268" "2548" "5108" "10228" "20468" "40948" )

# --- dl matrices ---
# for i in "${matrix_size[@]}"
# do
#     n=$i
#     wget "https://sparse.tamu.edu/RB/CPM/cz${n}.tar.gz"
# done


# --- untar and mv and rm ---
# for n in "${matrix_size[@]}"
# do
#     tar zxf cz${n}.tar.gz
#     mv cz${n}/cz${n}.rb ./
#     rm -r cz${n}
# done
# rm *.tar.gz

# --- create input files ---
