#!/bin/bash

# File paths for the two CSV files
read -p "Enter the path to file1: " file1
read -p "Enter the path to file2: " file2

# Read the first column of file1 into an array
IFS=$'\n' read -d '' -r -a file1_values < <(cut -d ',' -f1 "$file1")

# Read the first column of file2 into an array
IFS=$'\n' read -d '' -r -a file2_values < <(cut -d ',' -f1 "$file2")

# Find matching values
matching_values=()
for value in "${file1_values[@]}"; do
    if [[ " ${file2_values[*]} " == *" $value "* ]]; then
        matching_values+=("$value")
    fi
done

# Print the matching values
echo "Matching values:"
for value in "${matching_values[@]}"; do
    echo "$value"
done
