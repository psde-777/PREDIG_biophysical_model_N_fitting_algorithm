
set xtics font ", 20"
set ytics font ", 20"


#set grid ytics lt 0 lw 2 lc rgb "#555555"
#set grid xtics lt 0 lw 2 lc rgb "#555555" 

set terminal postscript eps enhanced color

mydir = "~/Documents

set key font ", 20"

set key bottom right

set xlabel "Time (hours)" font ", 24"
set ylabel "Glucan to glucose conversion (%)" font ", 24"

set xr[0:75]
set yr[0:105]
set output "Bura_comparison_glucan.eps"
#set title "Glucan to glucose conversion (%)" font ", 26"
p  "Output/expe_data/expe_saccharification_low_glc.txt" w lp title "Low" pt 4 ps 1.5 lw 2 lt rgb "#00000000" dt "." ,\
	"Output/expe_data/expe_saccharification_medium_glc.txt" w lp title "Medium" pt 6 ps 1.5 lw 2 lt rgb "#00000000" dt ".",\
	"Output/expe_data/expe_saccharification_high_glc.txt" w lp title "High" pt 8 ps 1.5 lw 2 lt rgb "#00000000" dt ".",\
	"mean_saccharification_low.txt"  u 1:2 w l notitle lw 8 lt rgb "#00000000" ,\
	"mean_saccharification_medium.txt"  u 1:2 w l notitle lw 8 lt rgb "#66666666" ,\
	"mean_saccharification_high.txt"  u 1:2 w l notitle lw 8 lt rgb "#AAAAAAAA" ,\


#set title "Xylan to xylose conversion" font ", 26"
set output "Bura_comparison_xylan.eps"
set ylabel "Xylan to xylose conversion (%)" font ", 24"

p  "Output/expe_data/expe_saccharification_low_xyl.txt" w lp title "Low" pt 4 ps 1.5 lw 2 lt rgb "#00000000" dt ".",\
	"Output/expe_data/expe_saccharification_medium_xyl.txt" w lp title "Medium" pt 6 ps 1.5 lw 2 lt rgb "#00000000" dt "." ,\
	"Output/expe_data/expe_saccharification_high_xyl.txt" w lp title "High" pt 8 ps 1.5 lw 2 lt rgb "#00000000" dt ".",\
	"mean_saccharification_low.txt"  u 1:4 w l notitle lw 8 lt rgb "#00000000" ,\
	"mean_saccharification_medium.txt"  u 1:4 w l notitle lw 8 lt rgb "#66666666" ,\
	"mean_saccharification_high.txt"  u 1:4 w l notitle lw 8 lt rgb "#AAAAAAAA" 


set ylabel "Glucan to glucose conversion (%)" font ", 24"
set output "Bura_comparison_glucan_only_expe_data.eps"
#set title "Glucan to glucose conversion (%)" font ", 26"
p  "Output/expe_data/expe_saccharification_low_glc.txt" w lp title "Low" pt 4 ps 1.5 lw 2 lt rgb "#00000000" dt "." ,\
	"Output/expe_data/expe_saccharification_medium_glc.txt" w lp title "Medium" pt 6 ps 1.5 lw 2 lt rgb "#00000000" dt ".",\
	"Output/expe_data/expe_saccharification_high_glc.txt" w lp title "High" pt 8 ps 1.5 lw 2 lt rgb "#00000000" dt "."


#set title "Xylan to xylose conversion" font ", 26"
set output "Bura_comparison_xylan_only_expe_data.eps"
set ylabel "Xylan to xylose conversion (%)" font ", 24"

p  "Output/expe_data/expe_saccharification_low_xyl.txt" w lp title "Low" pt 4 ps 1.5 lw 2 lt rgb "#00000000" dt ".",\
	"Output/expe_data/expe_saccharification_medium_xyl.txt" w lp title "Medium" pt 6 ps 1.5 lw 2 lt rgb "#00000000" dt "." ,\
	"Output/expe_data/expe_saccharification_high_xyl.txt" w lp title "High" pt 8 ps 1.5 lw 2 lt rgb "#00000000" dt "." 
	
