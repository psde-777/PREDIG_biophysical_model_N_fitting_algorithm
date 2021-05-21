COUNT=1
PIDCOUNT=1
echo "Running code";
readarray -t foldernames < foldernames.txt
readarray -t percentages < percentages.txt
for foldername in ${foldernames[*]};do
    cd ${foldername}
	for percentage in ${percentages[*]}; do
	    cd ${percentage}_percent_cry
	    ./code_4 &

		pids[${PIDCOUNT}]=$!	
		PIDCOUNT=$((${PIDCOUNT} + 1));
		
		joblist=($(jobs -p))



		while [ ${#joblist[*]} -ge 6 ]
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