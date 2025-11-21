#!/bin/bash

# Variables definition
terminal_width=$(tput cols) # Get the terminal width
separator=$(printf "%0.s=" $(seq 1 $terminal_width)) # Generate separator of correct length
flag=0
log_file_name="psblas_test_results.log"
base_dir=$(pwd)

# Define color codes
GREEN="\033[0;32m"
RED="\033[0;31m"
BLUE="\033[0;34m"
YELLOW="\033[33m" 
RESET="\033[0m"

# Function to center text
center_text() {
    local text="$1"
    printf "\033[0;32m%*s\033[0m\n" $(((${#text} + terminal_width) / 2)) "$text"
}




# Script start
clear

# Welcome message
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
echo ""


# Iterate through first-layer subdirectories
for dir in "$base_dir"/*/; do    
    # Skip the current directory itself
    if [ "$dir" = "." ]; then
        continue
    fi
    
    echo -e "${BLUE}${separator}${RESET}"
    echo -e "${BLUE}[INFO]\t  Entering directory: $(pwd)/$(basename "$dir")${RESET}"
    ( # excecute script in a subshell, otherwise the dir search will stop
        cd "$dir"

        # Check if autotest.sh exists before executing it
        if [ -f autotest.sh ]; then
            chmod +x autotest.sh
            ./autotest.sh 
        else
            echo -e "${YELLOW}[WARNING] autotest.sh not found in $(pwd). Skipping $(basename "$dir") kernel${RESET}"
        fi

        # Append contents of any .log file in the subdirectory to the main log file
        log_files=$(find . -maxdepth 1 -type f -name "*.log")
        if [ -n "$log_files" ]; then
            for log_file in $log_files; do
                cat "$log_file" >> "../${log_file_name}"
                echo ' ' >> "../${log_file_name}"
            done
        else
            echo -e "${YELLOW}[WARNING] No .log files found in $(pwd). Skipping log append.${RESET}"
        fi
    )

    # Return to the parent directory
    echo -e "${BLUE}[INFO]\t  Leaving directory: $(pwd)/$(basename "$dir")${RESET}"
    echo -e "${BLUE}${separator}${RESET}"
    echo ""
done

echo -e "${BLUE}[INFO]\t  Finished processing all subdirectories.${RESET}"


echo -e "${GREEN}[INFO]\t  All tests completed successfully. Results are logged in ${log_file_name}.${RESET}"