readarray -t directories < directories.txt
readarray -t percentages < percentages.txt

for i in $(seq 0 3);do
	echo $i
	cp plot_saccha.gnu ${directories[$i]}
	cd ${directories[$i]}
	gnuplot plot_saccha.gnu
	cd ..
done

