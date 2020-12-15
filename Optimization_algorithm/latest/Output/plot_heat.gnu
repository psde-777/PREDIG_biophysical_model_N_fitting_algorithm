
set xlabel "Degree of polymerization" font ", 14"

set cblabel "Average number in system" font "\:bold, 10"

set xtics font ", 14"
set ytics font ", 14"
 
#set grid ytics lt 0 lw 2 lc rgb "#555555"
#set grid xtics lt 0 lw 2 lc rgb "#555555"

#set term png 
#set term pngcairo size 1280, 768
#set term pngcairo size 800, 480
set terminal postscript eps enhanced
#set output "multiplot_DP_distrib.eps"

set palette defined (0.1 "blue",0.29 'cyan',0.36 '#27ad81',0.5 "yellow",0.64 "pink",0.71 "magenta",0.9 "red")

set key center right 

set logscale cb
set xr[0:1000]


#set multiplot layout 3, 1 title "Polymer length distribution for different enzyme combinations" font "\:bold, 18"
set cbrange[0.01:100000]

set xtics format ""
set x2tics

set ylabel "Time[arbitrary units]" font ", 14"

set xr[0:1000]
set yr[400:0]reverse
set title "DP distrib" font ", 16"

set output "DP_distrib.eps"
p "DP_distrib.txt" u 3:1:4 w image notitle# axes x2y1



#set xr[950:1000]
#set ylabel "Time[arbitrary units]" font ", 14"
#set title "CBH only (DP 950 to 1000)" font ", 16"
#p "DP_distrib_only_CBH.txt" u 3:1:4 w image notitle# axes x2y1

#unset ylabel
#set title "BGL missing" font ", 16"
#p "DP_distrib_no_BGL.txt" u 3:1:4 w image notitle #axes x2y1
