#!/bin/bash

#param1 : number of families (currently fixed to 1)
#param2 : number of Gens in family
#param3 : number of subfolders in a Gen
#param4 : DELTA
#param5 : number of cores to be used
#param6 : Resume interrupted fitting from this generation number


# read 5 tab-separated input parameters from fit_params.txt
read -r param1 param2 param3 param4 param5 < fit_settings.txt

# read a interruption status from generations.log
#read -r param6 < statoos.txt
last_gen=$(tail -1 "generations.log" | awk '{print $2}')

param6=`expr $last_gen + 1`
#param6=$(cut -f2 generations.log)

# print the values of the input parameters
echo "--------------------------------------------"
echo "number of families: $param1"
echo "Generations per family: $param2"
echo "Sub-folders per generation: $param3"
echo "DELTA: $param4"
echo "cores to be used: $param5"
echo "--------------------------------------------"


if [[ $last_gen -gt 0 ]]
then
	echo "Resume fitting from generation: $param6"
	./evo_gen_stat.sh $param2
	./evo_all_in_one.sh $param1 $param2 $param3 $param4 $param5 $param6
else
	./evo_all_in_one.sh $param1 $param2 $param3 $param4 $param5 1
fi


./evo_janitor.sh


echo "ALL DONE"
