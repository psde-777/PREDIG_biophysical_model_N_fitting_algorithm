readarray -t directories < directories.txt
readarray -t percentages < percentages.txt

for directory in ${directories[*]};do

    echo ${directory}
	for percentage in ${percentages[*]};do
		python3 calc_mean_interpolate.py 100 ${directory}/${percentage}_percent_lignin/Output/saccharification/saccharification_ ${directory}/mean_saccharification_${percentage}.txt foo &
	done
	cd ..
done

