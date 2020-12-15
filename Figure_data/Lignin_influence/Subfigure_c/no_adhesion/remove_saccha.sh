readarray -t directories < directories.txt
readarray -t percentages < percentages.txt

cd ${directories[0]}
for j in $(seq 0 7);do
	cd ${percentages[j]}_percent_lignin/
	rm Output/saccharification/*
	cd ..
done
cd ..


