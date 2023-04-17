#!/bin/bash

# Remove contents of keywords.txt file
echo -n "" > keywords.txt

# Loop through matching file names
for file in latest/Output/expe_data/expe_saccharification_*_glc.txt; do
    # Extract wild card entry from file name
    entry=$(echo $file | awk -F'_' '{print $4}')
    # Write entry to keywords.txt file
    echo $entry >> keywords.txt
done

