# $1 : Family
# $2 : Generation
# $3 : Runs in generation
# $4 and following: keywords of saccharification files, such as low, medium, high, etc
COUNT=$2
FAM=$1
GEN=$2
RUNS=$3
PIDCOUNT=1

echo "$#"
shift
shift
shift

echo "$#"
for i in $(seq 1 ${RUNS})
do
	for j in $(seq 1 $#)
	do
		echo "sees"
		python3 calc_mean_interpolate.py 20 family_${FAM}/Generation_${COUNT}/Run_${i}/Output/saccharification/saccharification_$1_ family_${FAM}/Generation_${COUNT}/Run_${i}/Output/mean_saccharification_$1.txt
		pids[${PIDCOUNT}]=$!
		shift
		PIDCOUNT=$((${PIDCOUNT} + 1));	
	done
#	PROCESSID=$!

done

for pid in ${pids[*]}; do
    wait $pid
done
wait $PROCESSID


