# $1 family
# $2 delta
# $3 generation number
# $4,5,6,7,8,9,10,11: Initial values to be changed by +-delta
COUNT=1
for i in -$2 0 $2
do
	for j in -$2 0 $2
	do
		for k in -$2 0 $2
		do
			for l in -$2 0 $2
			do
				for m in -$2 0 $2
				do
					#for n in -$1 0 $1
					#do
					#	for o in -$1 0 $1
					#	do
					#		for p in -$1 0 $1
					#		do
								#cp init_param_changer family_$1/Generation_$3/Run_${COUNT}/
								cp kin_param_changer family_$1/Generation_$3/Run_${COUNT}/
								#cp kin_param_changer.py family_$1/Generation_$3/Run_${COUNT}/
								rm family_$1/Generation_$3/Run_${COUNT}/Params/kinetic_parameters.txt
								#rm family_$1/Generation_$3/Run_${COUNT}/Params/initial_configuration_parameters.txt
								cd family_$1/Generation_$3/Run_${COUNT}
								#python3 kin_param_changer.py $4 $5 $6 $7 $8 $9 ${10} ${11} $i $j $k $l $m 0 0 0
								./kin_param_changer $4 $5 $6 $7 $8 $9 ${10} ${11} $i $j $k $l $m 0 0 0
								#./init_param_changer 100 100 100 100 100 1 0.095 0.28 0 0.23 0.25
								cd ../../../
								#echo $3 $4 $5 $6 $7 $8 $9 ${10} ${11}
								echo $(<family_$1/Generation_$3/Run_${COUNT}/Params/kinetic_parameters.txt)
								let COUNT++
					#		done
					#	done
					#done					
				done
			done
		done
	done
done
