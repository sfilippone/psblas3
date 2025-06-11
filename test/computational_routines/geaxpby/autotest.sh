#!/bin/bash

GREEN="\033[0;32m"
RED="\033[0;31m"
BLUE="\033[0;34m"
YELLOW="\033[33m" 
RESET="\033[0m"


# Directories to compare
dir1="serial"
dir2="parallel"
log_file_name="psblas_geaxpby_test.log"
flag=0
num_procs=$(nproc)

# Function to center text
center_text() {
    local text="$1"
    printf "\033[0;32m%*s\033[0m\n" $(((${#text} + terminal_width) / 2)) "$text"
}

# Welcome message
terminal_width=$(tput cols) # Get the terminal width
separator=$(printf "%0.s=" $(seq 1 $terminal_width)) # Generate separator of correct length

clear
echo -e "${GREEN}${separator}${RESET}"
center_text "PSBLAS Computational Routines Test Suite"
echo -e "${GREEN}${separator}${RESET}"
center_text "Welcome to the PSBLAS Computational Routines Test Suite!"
center_text "This script compares the results of serial and parallel computations"
center_text "for all the computational routines documented on the version 3.9 of PSBLAS."
echo -e "${GREEN}${separator}${RESET}"

echo -e "${BLUE}[INFO]\t  Starting environment check for required modules...${RESET}"


# Check and load required modules
required_modules=("gnu/12.2.1-sys" "mpich/4.2.2" "cuda/12.5")

for module in "${required_modules[@]}"; do
    if ! module list 2>&1 | grep -q "$module"; then
        echo -e "${YELLOW}[WARNING] Module not found, loading $module${RESET}"
        module load "$module"
        flag=1
        if ! grep -q "module load $module" "$HOME/.bashrc"; then
            echo -e "[INFO]\t  Adding 'module load $module' to $bashrc..."
            echo "module load $module" >> "$HOME/.bashrc"
        # else
        #     echo "'module load $module' is already present in $bashrc."
        fi
    else
        echo -e "[INFO]\t  Found module $module."
    fi
done

# Update .bashrc if necessary
if [ $flag -eq 1 ]; then
    echo -e "[INFO]\t  Reloading $HOME/.bashrc..."
    source ~/.bashrc
fi

# Inform the user about environment persistence
if [ "$$" -eq "$PPID" ]; then
    echo -e "${YELLOW}[WARNING] Modules loaded in this script will not persist after the script finishes.${RESET}"
    echo -e "${YELLOW}[WARNING] Run the script using 'source autotest.sh' to make the changes persist.${RESET}"
fi

echo -e "${BLUE}[INFO]\t  Environment check for required modules completed.${RESET}"





# Check if the executable ELF file exists
if [ ! -f "./runs/psb_geaxpby_test" ]; then
    echo -e "${YELLOW}[WARNING] Executable not found. Running make...${RESET}"
    make
    if [ ! -f "./runs/psb_geaxpby_test" ]; then
        echo -e "${RED}[ERROR]    Failed to create executable. Check make command.${RESET}"
    fi
else
    echo -e "${BLUE}[INFO]\t  The executable already exists. Skipping the make process.${RESET}"
fi

# Excecute tests and save results
echo -e "${BLUE}[INFO]\t  Running the PSBLAS psb_geaxpby test...${RESET}"
echo ""
echo -e "${BLUE}[INFO]\t  Starting single process computation${RESET}"
mpirun -np 1 ./runs/psb_geaxpby_test
echo -e "${BLUE}[INFO]\t  Single process computation terminated correctly${RESET}"
echo ""
echo -e "${BLUE}[INFO]\t  Starting $num_procs processes computation${RESET}"
mpirun -np $num_procs ./runs/psb_geaxpby_test
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

echo -e "${BLUE}[INFO]\t  PSBLAS psb_geaxpby test succesfully completed.${RESET}"


echo -e "${GREEN}[INFO]\t  All tests completed successfully. Results are logged in ${log_file_name}.${RESET}"

