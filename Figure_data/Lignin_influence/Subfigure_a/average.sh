readarray -t directories < directories.txt
readarray -t percentages < percentages.txt

for i in $(seq 0 3);do
	echo $i
	cd ${directories[$i]}
	for j in $(seq 1 5);do
		python3 calc_mean_interpolate.py 100 ${percentages[j]}_percent_lignin/Output/saccharification/saccharification_ mean_saccharification_${percentages[j]}.txt soos &
	done
	python3 calc_mean_interpolate.py 100 0_percent_lignin/Output/saccharification/saccharification_ mean_saccharification_0_percent_lignin.txt soos
	cd ..
done

