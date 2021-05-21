
set x2label "Degree of polymerization" font "\bold:, 26"
set ylabel offset -2.5 "Time [arbitrary units]" font "\bold:, 26"
set cblabel offset 12 "\nFraction of overall glucose" offset 1.5 font "\bold:, 26"

#set xtics font ", 28"
set ytics font ", 24"
 
#set grid ytics lt 0 lw 2 lc rgb "#555555"
#set grid xtics lt 0 lw 2 lc rgb "#555555"


set rmargin at screen 0.78
set lmargin at screen 0.15
set tmargin at screen 0.88


set terminal postscript eps enhanced

set palette defined (0.1 "black",0.29 'blue',0.5 "red",0.71 "orange",0.9 "yellow")

set key center right 

#fibrilLength = 100

set key font ", 24"

set xr[-0.5:200.5]
set x2tics 0,40,200 font ", 26"



set logscale cb
set cbrange[0.0001:1]
set cbtics font ", 26"



unset xtics







#set format y "%2.0t{/Symbol \264}10^{%L}"
set yr[72:0]reverse
set output "DP_distrib_only_EG.eps"
p "mean_DP_distrib_only_EG.txt" u 2:1:3 w image notitle

set yr[72:0]reverse
set output "DP_distrib_only_CBH.eps"
p "mean_DP_distrib_only_CBH.txt" u 2:1:3 w image notitle

set yr[72:0]reverse
set output "DP_distrib_all_cellulases.eps"
p "mean_DP_distrib_all_cellulases.txt" u 2:1:3 w image notitle