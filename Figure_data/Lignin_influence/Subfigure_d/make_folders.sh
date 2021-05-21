readarray -t directories < directories.txt
readarray -t mus < mus.txt
readarray -t sigmas < sigmas.txt
COUNT=0
for directory in ${directories[*]};do
    echo ${directory}
	for mus in ${mus[*]};do
    	for sigma in ${sigmas[*]};do
            let COUNT++
        	cp -r template ${directory}/folder_${COUNT}
    	done
	done
done

