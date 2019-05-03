#!/bin/bash
matrix_dir=../matrices
bin_dir=../bin
declare -a matrices=(
    "cz148"
    "cz308"
    "cz628"
    "cz1268"
    "cz2548"
    "cz5108"
    "cz10228"
    "cz20468"
    "cz40948"
)

declare -a iterations=(
  "1" "10" "25" "50" "100" "150" "200" "300" "600" "1000" "1500" "2000" "3000"
)

declare -a numprocs=(
  "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16"
)


declare -a swapchoices=(
    "psb_swap_persistent_" "psb_swap_nonpersistent_"
)

# -- for testing --
# declare -a iterations=(
#   "1" "10"
# )
declare -a numprocs=(
  "16"
)

declare -a matrices=(
    "cz20468"
    "cz40948"
)
# -----------------

# ----------- BUILD -----------
build=true
if [ "$build" = true ]; then
for swap in "${swapchoices[@]}"; do
for iter in "${iterations[@]}";  do
    make clean
    make ITERATIONS=${iter} SWAPCHOICE=${swap}
    mv halo_${swap}${iter}_iter bin/
done
done
fi

# ------------ RUN ------------
run=false
if [ "$run" = true ]; then
for matrix in "${matrices[@]}"; do
    rm run/*
    cd run
for swap in "${swapchoices[@]}"; do
for iter in "${iterations[@]}"; do
    # echo $matrix $iter
    bin=halo_${swap}${iter}_iter

    cp ${matrix_dir}/${matrix}/* ./
    cp ${bin_dir}/${bin} ./

    for np in "${numprocs[@]}"; do
    SCOREP_TIMER=gettimeofday \
    		SCOREP_EXPERIMENT_DIRECTORY=../results/scorep/${matrix}${swap}${iter} \
    		SCOREP_WRAPPER_COMPILER_FLAGS="-O2" \
    		SCOREP_MPI_ENABLE_GROUP=all \
		mpirun --map-by ppr:${np}:node -machinefile $PBS_NODEFILE ./${bin} ${matrix}.inp

    # mpirun -np ${np} ./${bin} ${matrix}.inp

    done
done
done
    write_to=../results/${matrix}.output
    mv halo_output.txt ${write_to}
    cd ..
done
fi
