# $1 = family
# $2 = Generation number
# $3 = Number of simus in generation
# $4 = Maximum number of cores available
COUNT=1
PIDCOUNT=1
PIDRUNNINGCOUNT=0
PIDADDEDCOUNT=0
echo "Running code";
readarray -t names < keywords.txt
for i in $(seq 1 $3)
do
	cd family_$1/Generation_$2/Run_${COUNT}
	for name in ${names[*]}; do
		./code_4 ${name} >/dev/null &
		pids[${PIDADDEDCOUNT}]=$!
		PIDCOUNT=$((${PIDCOUNT} + 1));
		PIDADDEDCOUNT=$((${PIDADDEDCOUNT} + 1));
		if [ $PIDADDEDCOUNT -ge $4  ]
		then
			PIDRUNNINGCOUNT=0
#			echo "Before: $PIDRUNNINGCOUNT"
			for pid in ${pids[*]}; do
				if ps -p $pid > /dev/null
				then
					PIDRUNNINGCOUNT=$((${PIDRUNNINGCOUNT} + 1));
				fi
			done
#			echo "After: $PIDRUNNINGCOUNT"
			if [ $PIDRUNNINGCOUNT -ge $4 ]
			then
				while [ $PIDRUNNINGCOUNT -gt 0 ]
				do
			#		echo "Sleeping 1 second"
					sleep 0.1
					PIDRUNNINGCOUNT=0
					for pid in ${pids[*]}; do
						if ps -p $pid > /dev/null
						then
							PIDRUNNINGCOUNT=$((${PIDRUNNINGCOUNT} + 1));
						fi
					done
				done
			fi
			PIDADDEDCOUNT=0
		fi

	done
	cd ../../../
	let COUNT++
done

for pid in ${pids[*]}; do
    wait $pid
done
echo "Done running code";
wait $PROCESSID


