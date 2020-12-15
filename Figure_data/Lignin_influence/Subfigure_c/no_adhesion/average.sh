readarray -t directories < directories.txt
readarray -t percentages < percentages.txt


for j in $(seq 0 11);do
	python3 calc_mean_interpolate.py $1 enzymes_50/${percentages[j]}_percent_lignin/Output/saccharification/saccharification_ mean_saccharification_${percentages[j]}.txt soos &
done



