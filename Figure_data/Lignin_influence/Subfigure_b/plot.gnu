
set xtics font ", 20"
set ytics font ", 20"


set terminal postscript eps enhanced color

set key font ", 20"


scale = 1


set key top right



set xlabel "Lignin content (%)\n" font ", 24"
set ylabel "Final released glucose (%)\n" font ", 26"

set xr[0:60]
set yr[-1:105]


f_50(x) = a*x + b
f_100(x) = c*x + d
f_150(x) = e*x + f


logistic_50(x) = 100/(1+exp(0.2*(x-5)))
logistic_100(x) = 98/(1+exp(0.2*(x-28)))
logistic_150(x) = 100/(1+exp(0.35*(x-44)))

fit f_50(x) "lignin_vs_final_yield_enzymes_50.txt" every ::0::3 via a,b
fit f_100(x) "lignin_vs_final_yield_enzymes_100.txt" every ::2::5 via c,d
fit f_150(x) "lignin_vs_final_yield_enzymes_150.txt" every ::4::7 via e,f



#set title "Comparison between simulated saccharification yields for various lignin percentages"
set output "final_yield_vs_lignin_enzymes.eps"

p [0:60]logistic_50(x) lw 5 lt rgb "#DDDDDDDD" notitle,\
[0:60]logistic_100(x) lw 5 lt rgb "#DDDDDDDD" notitle,\
[0:60]logistic_150(x) lw 5 lt rgb "#DDDDDDDD" notitle,\
"lignin_vs_final_yield_enzymes_50.txt" u ($1*scale):2 w p pt 5 ps 1.5 lt rgb "#00000000" title "[E] = 50",\
"lignin_vs_final_yield_enzymes_100.txt" u ($1*scale):2 w p pt 7 ps 1.7 lt rgb "#00000000" title "[E] = 100",\
"lignin_vs_final_yield_enzymes_150.txt" u ($1*scale):2 w p pt 9 ps 1.9 lt rgb "#00000000" title "[E] = 150",\
[0:53]f_50(x) lw 5 lt rgb "#00000000" notitle,\
[20:50]f_100(x) lw 5 lt rgb "#00000000" notitle,\
[40:60]f_150(x) lw 5 lt rgb "#00000000" notitle,\



