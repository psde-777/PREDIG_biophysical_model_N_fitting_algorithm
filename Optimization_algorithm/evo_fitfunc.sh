# $1 : family
# $2 : Generation
# $3 : Runs in generation
# $4 : Max number of CPU cores
PIDCOUNT=1
PIDRUNNINGCOUNT=0
PIDADDEDCOUNT=0
COUNT=$2
readarray -t names < keywords.txt
for i in $(seq 1 $3)
do
	echo $i
	cd family_$1/Generation_${COUNT}/Run_${i}/Output/
	for name in ${names[*]}; do
		python3 fitfunc.py ${name} &
		pids[${PIDADDEDCOUNT}]=$!
		PIDCOUNT=$((${PIDCOUNT} + 1));	
		PIDADDEDCOUNT=$((${PIDADDEDCOUNT} + 1));
		if [ $PIDADDEDCOUNT -ge $4  ]
		then
			PIDRUNNINGCOUNT=0
			echo "Before: $PIDRUNNINGCOUNT"
			for pid in ${pids[*]}; do
				if ps -p $pid > /dev/null
				then
					PIDRUNNINGCOUNT=$((${PIDRUNNINGCOUNT} + 1));
				fi
			done
			echo "After: $PIDRUNNINGCOUNT"
			if [ $PIDRUNNINGCOUNT -ge $4 ]
			then
				while [ $PIDRUNNINGCOUNT -gt 0 ]
				do
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
#	python3 fitfunc.py low &
#	pids[${i}]=$!
#	python3 fitfunc.py medium &
#	pids[${i}]=$!
#	python3 fitfunc.py high &
#	pids[${i}]=$!
#	PROCESSID=$!
	cd ../../../../
done

for pid in ${pids[*]}; do
	wait $pid
done
echo "Done with fitfuncs";
wait $PROCESSID
