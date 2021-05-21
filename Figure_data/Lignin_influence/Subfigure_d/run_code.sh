
PIDCOUNT=1
echo "Running code";

readarray -t directories < directories.txt
readarray -t mus < mus.txt
readarray -t sigmas < sigmas.txt
COUNT=0
for directory in ${directories[*]};do
    cd ${directory}
    echo ${directory}
	for mus in ${mus[*]};do
    		for sigma in ${sigmas[*]};do
        	    let COUNT++
		    cd folder_${COUNT}
		    ./code_4 &
		    pids[${PIDCOUNT}]=$!	
		    PIDCOUNT=$((${PIDCOUNT} + 1));
		    
		    joblist=($(jobs -p))
		


		    while [ ${#joblist[*]} -ge $1 ]
		    do
			    sleep 0.1
			    joblist=($(jobs -p))
		    done
		    cd ..
	   	 done
	done
done

for pid in ${pids[*]}; do
    wait $pid
done 
    
echo "Done"
