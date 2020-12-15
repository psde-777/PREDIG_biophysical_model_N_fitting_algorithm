
set xtics font ", 14"
set ytics font ", 14"


set terminal postscript eps enhanced color




scale = 100
set key top left


set xlabel "Time [arbitrary units]
set ylabel "{/Symbol m}g Glc per mg AIR"

set xr[0:30]

set title "Comparison between simulated saccharification yields for various lignin percentages"
set output "lignin_contents.eps"
p "mean_saccharification_6p25_perc_lignin.txt" u ($1*scale):2 w l title "6.25% lignin"lw 3 lt rgb "#00000000" ,\
"mean_saccharification_12p5_perc_lignin.txt" u ($1*scale):2 w l title "12.5% lignin"lw 3 lt rgb "#00000000" dt "." ,\
 "mean_saccharification_31p25_perc_lignin.txt" u ($1*scale):2 w l title "31.25% lignin"lw 3 lt rgb "#00000000" dt "-",\
 "mean_saccharification_37p5_perc_lignin.txt" u ($1*scale):2 w l title "37% lignin"lw 3 lt rgb "#00000000" dt "_",\
 "mean_saccharification_50_perc_lignin.txt" u ($1*scale):2 w l title "50% lignin"lw 3 lt rgb "#00000000" dt "._",\

set output "final_yield_vs_lignin.eps"
set title "Simulated final saccharification yield depending on lignin percentage
set xlabel "Lignin percentage"
set xr[0:55]
set yr[0:100]
set ylabel "Final yield [{/Symbol m}g Glc per mg AIR]"
f(x) = a*x + b
fit f(x) "lignin_vs_final_yield.txt" u 1:2 via a,b
p "lignin_vs_final_yield.txt" u 1:2 pt 7 ps 1 lt rgb "#000000" notitle, f(x) lw 3 lt rgb "#00000000" dt "-" title "f(x) = ".a."x + ".b
	


