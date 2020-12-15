
set xtics font ", 20"
set ytics font ", 20"


set terminal postscript eps enhanced color

set key font ", 20"


scale = 1
set key top left


set xlabel "Time [arbitrary units]\n" font ", 24"
set ylabel "Glucan to glucose conversion (%)" font ", 26"

set xr[0:0.07]
set yr[-1:105]
#set title "Comparison between simulated saccharification yields for various lignin percentages"
set output "lignin_contents.eps"
p "mean_saccharification_0.txt" u ($1*scale):2 w l title "No lignin"lw 4 lt rgb "#00000000" ,\
"mean_saccharification_10.txt" u ($1*scale):2 w l title "10% lignin"lw 4 lt rgb "#00000000" dt "." ,\
 "mean_saccharification_15.txt" u ($1*scale):2 w l title "15% lignin"lw 4 lt rgb "#00000000" dt "-",\
 "mean_saccharification_20.txt" u ($1*scale):2 w l title "20% lignin"lw 4 lt rgb "#00000000" dt "_",\
  "mean_saccharification_25.txt" u ($1*scale):2 w l title "25% lignin"lw 4 lt rgb "#00000000" dt "._",\
 "mean_saccharification_30.txt" u ($1*scale):2 w l title "30% lignin"lw 4 lt rgb "#00000000" dt ".-",\
 "mean_saccharification_40.txt" u ($1*scale):2 w l title "40% lignin"lw 4 lt rgb "#00000000" dt "._",\

# "mean_saccharification_35_percent_lignin.txt" u ($1*scale):2 w l title "35% lignin"lw 4 lt rgb "#88888888" dt "_",\
# "mean_saccharification_45_percent_lignin.txt" u ($1*scale):2 w l title "45% lignin"lw 4 lt rgb "#88888888" dt "_",\


set key top right
set output "final_yield_vs_lignin.eps"
set xlabel "Lignin percentage\n"
set xr[0:55]

#unset yr
set ylabel "Final glucose yield (%)"

p "lignin_vs_final_yield.txt" u 1:2 pt 13 ps 2 lt rgb "#000000" notitle#, f(x) lw 3 lt rgb "#00000000" dt "." title "Simulations",\


# [:51.4]untreat(x) lw 3 lt rgb "#00000000" dt "-" title "Untreated"
# [:26.7] treat(x) lw 3 lt rgb "#0000000" dt "_" title "Pretreated"



#sprintf("Y(p_{lign}) = %f p_{lign} + %f",a, b)
	


