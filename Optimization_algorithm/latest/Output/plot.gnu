
set xtics font ", 14"
set ytics font ", 14"


#set grid ytics lt 0 lw 2 lc rgb "#555555"
#set grid xtics lt 0 lw 2 lc rgb "#555555" 

set terminal postscript eps enhanced color





set key bottom right

set xlabel "Time[hours]"
set ylabel "Percentage glucan to glucose

set xr[0:75]
set yr[0:105]
set output "Bura_comparison_glucan.eps"
set title "Glucan to glucose conversion"
p  "~/Documents/Gitstuff/cornwall/Expe-Data/Bura/Glucose_release/CornStover_Low.txt" title "Bura low" pt 5 ps 1 lt rgb "#00000000" ,\
	"~/Documents/Gitstuff/cornwall/Expe-Data/Bura/Glucose_release/CornStover_Medium.txt" title "Bura medium" pt 1 ps 1 lt rgb "#00000000" ,\
	"~/Documents/Gitstuff/cornwall/Expe-Data/Bura/Glucose_release/CornStover_High.txt" title "Bura high" pt 3 ps 1 lt rgb "#00000000" ,\
	"mean_saccharification_low.txt"  u 1:2 w l title "Simu low" lw 3 lt rgb "#00000000" ,\
	"mean_saccharification_medium.txt"  u 1:2 w l title "Simu medium" lw 3 lt rgb "#00000000" dt ".",\
	"mean_saccharification_high.txt"  u 1:2 w l title "Simu high" lw 3 lt rgb "#00000000" dt "-"


set title "Xylan to xylose conversion"
set output "Bura_comparison_xylan.eps"

p  "~/Documents/Gitstuff/cornwall/Expe-Data/Bura/Xylan_release/CornStover_Low.txt" title "Bura low" pt 5 ps 1 lt rgb "#00000000" ,\
	"~/Documents/Gitstuff/cornwall/Expe-Data/Bura/Xylan_release/CornStover_Medium.txt" title "Bura medium" pt 1 ps 1 lt rgb "#00000000" ,\
	"~/Documents/Gitstuff/cornwall/Expe-Data/Bura/Xylan_release/CornStover_High.txt" title "Bura high" pt 3 ps 1 lt rgb "#00000000" ,\
	"mean_saccharification_low.txt"  u 1:4 w l title "Simu low" lw 3 lt rgb "#00000000" ,\
	"mean_saccharification_medium.txt"  u 1:4 w l title "Simu medium" lw 3 lt rgb "#00000000" dt ".",\
	"mean_saccharification_high.txt"  u 1:4 w l title "Simu high" lw 3 lt rgb "#00000000" dt "-"
#set title "Enzyme activity"
#p  "mean_enzyme_activity.txt"  u 1:2 w l title "EG" lw 3 lt rgb "#00000000",\
#	"mean_enzyme_activity.txt"  u 1:3 w l title "CBH" lw 3 lt rgb "#00000000"dt "." ,\
#	"mean_enzyme_activity.txt"  u 1:4 w l title "BGL" lw 3 lt rgb "#00000000" dt "-",\
#	"mean_enzyme_activity.txt"  u 1:5 w l title "XYL" lw 3 lt rgb "#88888888" 
#
#set output "mean_enzyme_fraction.eps"
#set title "Enzyme fraction"
#p  "mean_enzyme_fraction.txt"  u 1:2 w l title "EG" lw 3 lt rgb "#00000000",\
#	"mean_enzyme_fraction.txt"  u 1:3 w l title "CBH" lw 3 lt rgb "#00000000"dt "." ,\
#	"mean_enzyme_fraction.txt"  u 1:4 w l title "BGL" lw 3 lt rgb "#00000000" dt "-",\
#	"mean_enzyme_fraction.txt"  u 1:5 w l title "XYL" lw 3 lt rgb "#88888888" 
