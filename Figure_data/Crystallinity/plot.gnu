
set xtics font ", 20"
set ytics font ", 20"


#set grid ytics lt 0 lw 2 lc rgb "#555555"
#set grid xtics lt 0 lw 2 lc rgb "#555555" 

set terminal postscript eps enhanced color

sc_dot001(x) = a*x + b
sc_dot1(x) = c*x + d
sc_0(x) = e*x + f
sc_dot023(x) = l*x + m

f_cui(x) = -1.053*(x) + 124.34 #Parameters provided by Cui et al (2014)
f_pena(x) = g*x + h #Parameters are fit to data provided by Pena et al (2019)

#fit sc(x) "saccharification_vs_crystallinity.txt" u (100-$1*100):3 via a, b


fit sc_dot001(x) "final_saccharification_low_digestibility.txt" u ($1):2 via a, b
fit sc_dot1(x) "final_saccharification_high_digestibility.txt" via c, d
fit sc_0(x) "final_saccharification_no_digestibility.txt" via e, f
fit sc_dot023(x) "final_saccharification_fitted_digestibility.txt" u ($1):2 via l, m
fit f_pena(x) "data_by_pena_et_al.txt" via g, h

set key font ", 18"

set key bottom left

set xlabel "Cellulose crystallinity fraction (%) \n" font ", 24" offset 0,-0.5
set ylabel "Final glucan to glucose conversion (%)" font ", 24"

set xr[0:100]
set yr[0:105]
set output "Bura_comparison_crystallinity.eps"
#set title "Final saccharification yield vs digestibility" font ", 26"


set key spacing 1.2

p  	"final_saccharification_high_digestibility.txt" w p pt 5 ps 1.5 lt rgb "#00000000" title "r_{c,a} =      10^{-1}    ",\
	"final_saccharification_fitted_digestibility.txt" w p pt 9 ps 2 lt rgb "#00000000" title "r_{c,a} = 3.0’{\267} 10^{-2}",\
	"final_saccharification_low_digestibility.txt" w p  pt 7 ps 1.8 lt rgb "#00000000" title "r_{c,a} =      10^{-3}    " ,\
	"data_by_Cui_et_al.txt" w p pt 4 ps 1.5 lt rgb "#00000000" title "Cui et al. (2014)",\
	"data_by_pena_et_al.txt" w p pt 6 ps 1.8 lt rgb "#00000000" title "Pena et al. (2019)",\
	[0:90] sc_dot1(x) notitle lw 5 lt rgb "#00000000",\
	[0:90] sc_dot023(x) notitle lw 5 lt rgb "#00000000",\
	[0:90] sc_0(x) notitle lw 5 lt rgb "#00000000",\
	[20:75] f_cui(x) notitle lw 5 lt rgb "#BBBBBBBB" ,\
	[10:50] f_pena(x) notitle lw 5 lt rgb "#BBBBBBBB" ,\



unset arrow 1
unset arrow 2
#unset arrow 3


