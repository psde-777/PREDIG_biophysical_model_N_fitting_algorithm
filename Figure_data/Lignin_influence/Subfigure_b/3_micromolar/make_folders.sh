readarray -t directories < directories.txt
readarray -t percentages < percentages.txt

for directory in ${directories[*]};do

    echo ${directory}
	for percentage in ${percentages[*]};do
        cp -r template ${directory}/${percentage}_percent_lignin
	done
done

