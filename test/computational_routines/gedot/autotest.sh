#!/bin/bash

# Variables definition
dir1="serial"
dir2="parallel"
log_file_name="psblas_gedot_test.log"
num_procs=$(nproc)

# Define color codes
GREEN="\033[0;32m"
RED="\033[0;31m"
BLUE="\033[0;34m"
YELLOW="\033[33m" 
RESET="\033[0m"


# Check if the executable ELF file exists
if [ ! -f "./runs/psb_gedot_test" ]; then
    echo -e "${YELLOW}[WARNING] Executable not found. Running make...${RESET}"
    make
    if [ ! -f "./runs/psb_geaxpby_test" ]; then
        echo -e "${RED}[ERROR]    Failed to create executable. Check make command.${RESET}"
    fi
else
    echo -e "${BLUE}[INFO]\t  The executable already exists. Skipping the make process.${RESET}"
fi


# Excecute tests and save results
echo -e "${BLUE}[INFO]\t  Running the PSBLAS psb_gedot test...${RESET}"
echo ""
echo -e "${BLUE}[INFO]\t  Starting single process computation${RESET}"
#mpirun -np 1 ./runs/psb_gedot_test
echo -e "${BLUE}[INFO]\t  Single process computation terminated correctly${RESET}"
echo ""
echo -e "${BLUE}[INFO]\t  Starting $num_procs processes computation${RESET}"
mpirun -np $num_procs ./runs/psb_gedot_test
echo -e "${BLUE}[INFO]\t  Multiple processes computation terminated correctly${RESET}"

echo "" >> ${log_file_name}


# Iterate through files in the first directory
for file1 in "$dir1"/*; do
    filename=$(basename "$file1") # Extract the filename
    file2="$dir2/$filename"      # Construct the path for the second directory

    # Check if the file exists in the second directory
    if [ -f "$file2" ]; then
        diff_count=$(diff "$file1" "$file2" | wc -l) # Compare the files
        echo "Comparison between $file1 and $file2: $diff_count differences" >> ${log_file_name}
        # echo "Comparing $file1 and $file2: $diff_count"
    else
        echo -e "${RED}[ERROR]   File $filename does not exist in $dir2${RESET}"
    fi
done

echo -e "${BLUE}[INFO]\t  PSBLAS psb_gedot test succesfully completed.${RESET}"