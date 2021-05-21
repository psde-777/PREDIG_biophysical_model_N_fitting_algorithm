
PIDCOUNT=1
echo "Averaging";

readarray -t directories < directories.txt
readarray -t mus < mus.txt
readarray -t sigmas < sigmas.txt
COUNT=0
for directory in ${directories[*]};do
    echo ${directory}
	for mus in ${mus[*]};do
		for sigma in ${sigmas[*]};do
    	    let COUNT++
            python3 calc_mean_interpolate.py $1 ${directory}/folder_${COUNT}/Output/saccharification/saccharification_ ${directory}/mean_saccharification_${COUNT}.txt foo &
		    pids[${PIDCOUNT}]=$!	
		    PIDCOUNT=$((${PIDCOUNT} + 1));
		    
		    joblist=($(jobs -p))
		


		    while [ ${#joblist[*]} -ge $2 ]
		    do
			    sleep 0.1
			    joblist=($(jobs -p))
		    done
	   	 done
	done
done

for pid in ${pids[*]}; do
    wait $pid
done 
    
echo "Done"











