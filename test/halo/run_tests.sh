#!/bin/bash
matrix_dir=../matrices
bin_dir=../bin
    # "cz148"
    # "cz308"
declare -a matrices=(
    "cz628"
    "cz1268"
    "cz2548"
    "cz5108"
    "cz10228"
    "cz20468"
    "cz40948"
)

declare -a iterations=(
    "10" "300" "600" "900" "1200" "1500" "1800" "2100" "2400" "2700" "3000"
)
# declare -a iterations=(
#     "10" "100"
# )


# ----------- BUILD -----------
# for iter in "${iterations[@]}"; do
#     make clean
#     make ITERATIONS=${iter}
#     mv halo_persistant_${iter}_iter bin/
# done


# ------------ RUN ------------



# go to directory to run
for matrix in "${matrices[@]}"; do
for iter in "${iterations[@]}"; do
    # echo $matrix $iter
    rm run/*
    cd run
    bin=halo_persistant_${iter}_iter
    # matrix=cz308
    cp ${matrix_dir}/${matrix}/* ./
    cp ${bin_dir}/${bin} ./

    mkdir ../results/halo-persistant/${matrix}_matrix
    write_to=../results/halo-persistant/${matrix}_matrix/${iter}_iter_${matrix}_matrix
    SCOREP_TIMER=gettimeofday \
    		SCOREP_EXPERIMENT_DIRECTORY=${write_to} \
    		SCOREP_WRAPPER_COMPILER_FLAGS="-O2" \
    		SCOREP_MPI_ENABLE_GROUP=all \
    		mpirun --map-by ppr:16:node -machinefile $PBS_NODEFILE ./${bin} ${matrix}.inp
    mv halo_persistant_output.txt ${write_to}
    cd ..
done
done
