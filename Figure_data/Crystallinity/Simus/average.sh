COUNT=1
PIDCOUNT=1
echo "Running code";
readarray -t foldernames < foldernames.txt
readarray -t percentages < percentages.txt
for foldername in ${foldernames[*]};do

	for percentage in ${percentages[*]}; do

	    python3 calc_mean_interpolate.py 100 ${foldername}/${percentage}_percent_cry/Output/saccharification/saccharification_ ${foldername}/${percentage}_percent_cry/Output/mean_saccharification.txt ~/Documents/cornwall/Code/latest/Output/expe_data/expe_saccharification_low_&

		pids[${PIDCOUNT}]=$!	
		PIDCOUNT=$((${PIDCOUNT} + 1));
		
		joblist=($(jobs -p))



		while [ ${#joblist[*]} -ge 5 ]
		do
			sleep 0.1
			joblist=($(jobs -p))
		done




	done

done

for pid in ${pids[*]}; do
    wait $pid
done
    
    
echo "Done"