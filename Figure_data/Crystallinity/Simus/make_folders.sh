COUNT=1
PIDCOUNT=1
readarray -t foldernames < foldernames.txt
readarray -t percentages < percentages.txt
for foldername in ${foldernames[*]};do
    echo ${foldername}
    cd ${foldername}
    cd template
    ./compile.sh
    cd ..
	for percentage in ${percentages[*]}; do
	    cp -r template ${percentage}_percent_cry
	done
	cd ..
done