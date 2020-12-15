readarray -t directories < directories.txt
readarray -t percentages < percentages.txt

cd ${directories[0]}
for j in $(seq 0 11);do
	cd ${percentages[j]}_percent_lignin/
	./code_4 >/dev/null &
	cd ..
done
cd ..


