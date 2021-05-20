#$1 : number of families
#$2 : number of generations
#$3 : number of folders in each generation
#$4 : DELTA
#$5 : Max number of cores to be used
#(echo "99999999">smallest_var.txt)
(echo "0">no_reduction_count.txt)
(echo "Generation 1")
(python3 reset_gradient.py)
(./evo_create_folders.sh $1 1 $3)
PROCESSID=$!
wait $PROCESSID
(python3 rando_specs.py 1 $3 $4 $1)
PROCESSID=$!
wait $PROCESSID
(echo "Running")
(./evo_run.sh $1 1 $3 $5)
PROCESSID=$!
wait $PROCESSID
(echo "Averaging")
(./evo_average_fit_and_print_variance.sh $1 1 $3 $5)
PROCESSID=$!
wait $PROCESSID
(rm vars.txt)
PROCESSID=$!
wait $PROCESSID
(./evo_pick_min.sh $1 1 $3)
PROCESSID=$!
wait $PROCESSID
(python3 find_min_vars.py $1 1)
PROCESSID=$!
wait $PROCESSID
echo $(<best_candidate_specs.txt) || echo "No best candidate file!"

cp vars.txt family_$1/vars_1.txt
./evo_clear_gen.sh $1 1 $(<best_candidate.txt)
PROCESSID=$!
wait $PROCESSID
#while (true)
#do
#	echo blabla
#done
#(./evo_create_folders.sh $1 $2) && (python3 rando_kin_specs.py 0.1 1 $(<run_specs.txt)) && (./evo_run.sh 1 $2) && (while(true) do echo gaga; done) &&(./evo_average.sh 1 $2) && (./evo_fitfunc.sh 1 $2) && (./evo_print_variance.sh 1 $2) && (rm vars.txt) && (rm best_candidate*) && (./evo_pick_min.sh 1 $2) && (python3 find_min_vars.py 1) &&

(for i in $(seq 2 $2)
do
#	(python3 rando_kin_specs.py 0.1 $i $(<run_specs.txt)) && (./evo_run.sh $i $2) && (./evo_average.sh $i $2) && (./evo_fitfunc.sh $i $2) && (./evo_print_variance.sh $i $2) && (rm vars.txt) && (rm best_candidate*) && (./evo_pick_min.sh $i $2) && (python3 find_min_vars.py $i)

	#python3 rando_kin_specs.py 0.1 $i $(<best_candidate_specs.txt) && ./evo_run.sh $i $2 && ./evo_average.sh $i $2 &&./evo_fitfunc.sh $i $2 && ./evo_print_variance.sh $i $2 && rm vars.txt && rm best_candidate* && ./evo_pick_min.sh $i $2 && python3 find_min_vars.pi $i
    (echo "Generation $i")
	./evo_create_folders.sh $1 ${i} $3
	PROCESSID=$!
	wait $PROCESSID
	(python3 rando_specs.py $i $3 $4 $1)
	PROCESSID=$!
	wait $PROCESSID
    (echo "Running")
	./evo_run.sh $1 $i $3 $5
	PROCESSID=$!
	wait $PROCESSID
    (echo "Averaging")
	./evo_average_fit_and_print_variance.sh $1 $i $3 $5
	PROCESSID=$!
	wait $PROCESSID
	rm vars.txt
	./evo_pick_min.sh $1 $i $3 
	PROCESSID=$!
	wait $PROCESSID
    python3 find_min_vars.py $1 $i
	PROCESSID=$!
	wait $PROCESSID
	cp vars.txt family_$1/vars_$i.txt
	PROCESSID=$!
	wait $PROCESSID
	./evo_clear_gen.sh $1 $i $(<best_candidate.txt)
	PROCESSID=$!
	wait $PROCESSID

done
)

