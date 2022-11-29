#$1 : number of families
#$2 : number of generations


mkdir TOP_5_Gens

touch TOP_5_Gens/var_vs_gen.txt



i=1
while [[ $i -le $2 ]]
do

#	echo "$i 	 "$(sort -g -k 2,2 family_1/vars_"$i".txt | head -n 1 |  cut -d ' ' -f2)"" >> TOP_5_Gens/var_vs_gen.txt
	echo "$i 	 "$(sort -g -k 2,2 family_1/vars_"$i".txt | head -n 1)"" >> TOP_5_Gens/var_vs_gen.txt

	((i = i + 1))

done


touch TOP_5_Gens/lowest_vars_gen.txt

sort -g -k 3,3 TOP_5_Gens/var_vs_gen.txt | head -n 5 >> TOP_5_Gens/lowest_vars_gen.txt





echo "Copying best 5 Gens"


while IFS=" " read gen_no run_no gen_var
do
	cp -r family_1/Generation_"$gen_no" TOP_5_Gens
    
done < "TOP_5_Gens/lowest_vars_gen.txt"
