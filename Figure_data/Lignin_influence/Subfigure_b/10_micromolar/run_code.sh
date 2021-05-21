COUNT=1
PIDCOUNT=1
echo "Running code";
readarray -t foldernames < directories.txt
readarray -t percentages < percentages.txt
for foldername in ${foldernames[*]};do
    echo ${foldername}
    cd ${foldername}
	for percentage in ${percentages[*]}; do
	    cd ${percentage}_percent_lignin
	    ./code_4 &

		pids[${PIDCOUNT}]=$!	
		PIDCOUNT=$((${PIDCOUNT} + 1));
		
		joblist=($(jobs -p))



		while [ ${#joblist[*]} -ge 5 ]
		do
			sleep 0.1
			joblist=($(jobs -p))
		done



	    cd ..
	done
	cd ..
done

for pid in ${pids[*]}; do
    wait $pid
    
    
echo "Done"