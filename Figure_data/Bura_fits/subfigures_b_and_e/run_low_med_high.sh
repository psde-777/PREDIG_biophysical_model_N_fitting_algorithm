#$1: number of runs

# This scripts runs the simulation code for low, medium and high pre-treatment severity at the same time
COUNT=0
./code_4 low >/dev/null &
pids[${COUNT}]=$!
let COUNT++
./code_4 medium >/dev/null &
pids[${COUNT}]=$!
let COUNT++
./code_4 high >/dev/null &
pids[${COUNT}]=$!
echo "Running"
for pid in ${pids[*]}; do
    wait $pid
done
echo "Done"
