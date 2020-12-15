readarray -t directories < directories.txt
readarray -t percentages < percentages.txt

for i in $(seq 0 3);do
	echo $i
	cd ${directories[$i]}
	for j in $(seq 1 5);do
		cd ${percentages[j]}_percent_lignin/
		./code_4 >/dev/null &
		cd ..
	done
	cd 0_percent_lignin/
	./code_4 >/dev/null
	cd ..
	cd ..
done

