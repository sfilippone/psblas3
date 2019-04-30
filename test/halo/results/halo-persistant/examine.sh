#!/bin/bash

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
    "10" "300" "600" "900" "1200" "1500" "1800" "2100" "2400" "2700" "3000"
)



for matrix in "${matrices[@]}"; do
    echo "${matrix}"
    eval "scorep-score -r cz628_matrix/10_iter_cz628_matrix/profile.cubex | grep -e max_buf"

    # for iter in "${iterations[@]}"; do
    # 	cmd="scorep-score -r ${matrix}_matrix/${iter}_iter_${matrix}_matrix/profile.cubex | grep -e MPI_Wait"
    # 	eval $cmd
    # done

    # for iter in "${iterations[@]}"; do
    # 	cmd="scorep-score -r ${matrix}_matrix/${iter}_iter_${matrix}_matrix/profile.cubex | grep -e MPI_Start"
    # 	eval $cmd
    # done

    # for iter in "${iterations[@]}"; do
    # 	cmd="scorep-score -r ${matrix}_matrix/${iter}_iter_${matrix}_matrix/profile.cubex | grep -e MPI_Dist_graph_create_adjacent"
    # 	eval $cmd
    # done

    for iter in "${iterations[@]}"; do
    	cmd="scorep-score -r ${matrix}_matrix/${iter}_iter_${matrix}_matrix/profile.cubex | grep -e MPI_Dist_graph_neighbors"
    	eval $cmd
    done

    # for iter in "${iterations[@]}"; do
    # 	cmd="scorep-score -r ${matrix}_matrix/${iter}_iter_${matrix}_matrix/profile.cubex | grep -e MPI_Dist_graph_neighbors_count"
    # 	eval $cmd
    # done
done
