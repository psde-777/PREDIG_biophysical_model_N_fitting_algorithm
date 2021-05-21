
set xtics font ", 20"
set ytics font ", 20"

set rmargin at screen 0.8
#set lmargin at screen 0.15
#set tmargin at screen 0.88


set palette defined (0.1 "black",0.79 'blue',0.86 "red",0.95 "orange",0.99 "yellow")
#set palette defined (0.1 "black",0.29 'blue',0.5 "red",0.71 "orange",0.9 "yellow")

set terminal postscript eps enhanced color

set key font ", 20"


scale = 1


set key top left
unset key


set xlabel "Mean polymer covering fraction (%) \n" font ", 24"
set ylabel "Standard deviation (%) \n" font ", 24"

set cblabel offset 3 "Glucan to glucose conversion (%)" font ", 24"
set cbtics font ", 24"
set cbtics 0,20,100 font ", 24"

set xr[2:100]
set yr[2:100]
set cbrange[0:100]
directory="data"


set output "Figure_5_f.eps"
p "final_saccharification_data.txt" u($1*100):($2*100):3 w image