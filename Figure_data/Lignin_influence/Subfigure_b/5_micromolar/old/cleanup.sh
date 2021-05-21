readarray -t directories < directories.txt
readarray -t percentages < percentages.txt

cd ${directories[0]}
for j in $(seq 0 7);do
	cd ${percentages[j]}_percent_lignin/
	rm Params/initial_configuration_parameters_*
	rm -r Params/low_med_high*
	rm -r Params/optimized_inits
	cd ..
done
cd ..


