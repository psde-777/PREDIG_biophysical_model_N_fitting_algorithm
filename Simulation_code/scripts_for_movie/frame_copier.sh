#!/bin/bash

# find the highest number in the files
highest_num=$(ls Output/3D/visualisation_total_*.txt | awk -F"_" '{print $NF}' | awk -F"." '{print $1}' | sort -n | tail -n1)



# make 30 copies of the file with the correct numbering
for i in $(seq 1 100); do
    new_num=$(expr $highest_num + $i)
    cp "Output/3D/visualisation_total_1_$highest_num.txt" "Output/3D/visualisation_total_1_$new_num.txt"
    cp "Output/3D/conversion_1_$highest_num.txt" "Output/3D/conversion_1_$new_num.txt"
done
