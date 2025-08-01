#!/bin/bash

# Variables definition
num_procs=$(nproc)

# Define color codes
GREEN="\033[0;32m"
RED="\033[0;31m"
BLUE="\033[0;34m"
YELLOW="\033[33m" 
RESET="\033[0m"


# Check if the executable ELF file exists
if [ ! -f "./runs/psb_gedots_test" ]; then
    echo -e "${YELLOW}[WARNING] Executable not found. Running make...${RESET}"
    make
    if [ ! -f "./runs/psb_geaxpbys_test" ]; then
        echo -e "${RED}[ERROR]   Failed to create executable. Check make command.${RESET}"
    fi
else
    echo -e "${BLUE}[INFO]\t  The executable already exists. Skipping the make process.${RESET}"
fi


# Excecute tests and save results
echo -e "${BLUE}[INFO]\t  Running the PSBLAS psb_gedots test...${RESET}"
echo ""
echo -e "${BLUE}[INFO]\t  Starting single process computation${RESET}"
mpirun -np 1 ./runs/psb_gedots_test
echo -e "${BLUE}[INFO]\t  Single process computation terminated correctly${RESET}"
echo ""
echo -e "${BLUE}[INFO]\t  Starting $num_procs processes computation${RESET}"
mpirun -np $num_procs ./runs/psb_gedots_test
echo -e "${BLUE}[INFO]\t  Multiple processes computation terminated correctly${RESET}"

rm -f results/*
rm -f -r results/

echo -e "${GREEN}[INFO]\t  PSBLAS psb_gedots test succesfully completed.${RESET}"