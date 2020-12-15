# $1 : Family
# $2 : Generation
# $3 : Runs in generation
# $4 : Max number of CPUs
COUNT=$2
FAM=$1
GEN=$2
RUNS=$3
#RUNS=$((${RUNS}+1))
echo "Runs = ${RUNS}"
PIDCOUNT=1
PIDRUNNINGCOUNT=0
PIDADDEDCOUNT=0
#shift
#shift
#shift
readarray -t names < keywords.txt
#for i in $(seq 1 $#)
#do
#	names[$i] = $1
#	echo "$1 and ${${names}[$i]}"
#	shift
#done

for i in $(seq 1 ${RUNS})
do
	echo "in evo_average: i = $i"
	for name in ${names[*]}; do
#		echo "Hi, my name is ${name}"
		python3 calc_mean_interpolate.py 20 family_${FAM}/Generation_${COUNT}/Run_${i}/Output/saccharification/saccharification_${name}_ family_${FAM}/Generation_${COUNT}/Run_${i}/Output/mean_saccharification_${name}.txt family_${FAM}/Generation_${COUNT}/Run_${i}/Output/expe_data/expe_saccharification_${name}_ >/dev/null &
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
#	PROCESSID=$!

done

for pid in ${pids[*]}; do
    wait $pid
done
wait $PROCESSID


