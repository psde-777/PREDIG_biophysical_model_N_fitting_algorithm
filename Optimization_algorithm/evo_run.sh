# $1 = family
# $2 = Generation number
# $3 = Number of simus in generation
# $4 = Maximum number of cores available
COUNT=1
PIDCOUNT=1
echo "Running code";
readarray -t names < keywords.txt
for i in $(seq 1 $3)
do
	cd family_$1/Generation_$2/Run_${COUNT}
	for name in ${names[*]}; do
		./code_4 ${name} >/dev/null &
		pids[${PIDCOUNT}]=$!	
		PIDCOUNT=$((${PIDCOUNT} + 1));
		
		joblist=($(jobs -p))



		while [ ${#joblist[*]} -ge $4 ]
		do
			sleep 0.1
			joblist=($(jobs -p))
		done

	done
	cd ../../../
	let COUNT++
done

for pid in ${pids[*]}; do
    wait $pid
done
echo "Done running code";
wait $PROCESSID


