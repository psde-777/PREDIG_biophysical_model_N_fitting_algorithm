
set xtics font ", 20"
set ytics font ", 20"


set terminal postscript eps enhanced color

scale = 1


f_studer_high_s_to_g_untreated(x) = m*x + n
f_studer_low_s_to_g_untreated(x) = o*x + p
f_studer_high_s_to_g_pretreated(x) = q*x + r
f_studer_low_s_to_g_pretreated(x) = s*x + t
f_chen_untreated(x) = -0.4062*x + 208.86
f_chen_pretreated(x) = -2.604*x + 695.39

fit f_studer_high_s_to_g_untreated(x)"Studer_et_al_data_high_S_to_G_ratio_untreated.csv" u 1:2 via m,n
fit f_studer_low_s_to_g_untreated(x) "Studer_et_al_data_low_S_to_G_ratio_untreated.csv" u 1:2 via o,p
fit f_studer_high_s_to_g_pretreated(x) "Studer_et_al_data_high_S_to_G_ratio_180_degrees.csv" u 1:2 via q,r
fit f_studer_low_s_to_g_pretreated(x) "Studer_et_al_data_low_S_to_G_ratio_180_degrees.csv" u 1:2 via s,t



set xlabel "Lignin content (%)\n" font ", 24"
set ylabel "Glucan to glucose conversion (%)\n" font ", 26"

set xr[0:50]
set yr[-1:110]


f_1(x) = a*x + b
f_2(x) = c*x + d
f_3(x) = e*x + f
f_5(x) = g*x + h
f_10(x) = i*x + j


logistic_1(x) = 100/(1+exp(0.4*(x-c2)))
logistic_2(x) = 90/(1+exp(0.35*(x-c4)))
logistic_3(x) = 96/(1+exp(0.25*(x-c6)))
logistic_5(x) = 100/(1+exp(0.2*(x-c8)))
logistic_10(x) = 100/(1+exp(0.1*(x-c10)))

fit f_1(x) "final_saccharification_data_1_micromolar.txt" every ::0::1 via a,b
fit f_2(x) "final_saccharification_data_2_micromolar.txt" every ::1::3 via c,d
fit f_3(x) "final_saccharification_data_3_micromolar.txt" every ::2::4 via e,f
fit f_5(x) "final_saccharification_data_5_micromolar.txt" every ::5::7 via g,h
fit f_10(x) "final_saccharification_data_10_micromolar.txt" every ::4::7 via i,j
fit logistic_1(x) "final_saccharification_data_1_micromolar.txt" via c2
fit logistic_2(x) "final_saccharification_data_2_micromolar.txt" via c4
fit logistic_3(x) "final_saccharification_data_3_micromolar.txt" via c6
fit logistic_5(x) "final_saccharification_data_5_micromolar.txt" via c8
fit logistic_10(x) "final_saccharification_data_10_micromolar.txt" via c10



#set title "Comparison between simulated saccharification yields for various lignin percentages"
set output "Figure_5_b.eps"

#set multiplot

set label "1 {/Symbol m}M" at 1,70 rotate by 0 font ", 20"
set label "2 {/Symbol m}M" at 10,80 rotate by 0 font ", 20"
set label "3 {/Symbol m}M" at 19,85 rotate by 0 font ", 20"
set label "5 {/Symbol m}M" at 43,86 rotate by 0 font ", 20"
set label "10 {/Symbol m}M" at 52,104 rotate by 0 font ", 20"

p [0:60]logistic_1(x) lw 5 lt rgb "#DDDDDDDD" notitle,\
[0:60]logistic_2(x) lw 5 lt rgb "#DDDDDDDD" notitle,\
[0:60]logistic_3(x) lw 5 lt rgb "#DDDDDDDD" notitle,\
[0:60]logistic_5(x) lw 5 lt rgb "#DDDDDDDD" notitle,\
[0:60]logistic_10(x) lw 5 lt rgb "#DDDDDDDD" notitle,\
"final_saccharification_data_1_micromolar.txt" u ($1*scale):2 w p pt 5 ps 1.5 lt rgb "#00000000" notitle,\
"final_saccharification_data_2_micromolar.txt" u ($1*scale):2 w p pt 7 ps 1.7 lt rgb "#00000000" notitle,\
"final_saccharification_data_3_micromolar.txt" u ($1*scale):2 w p pt 11 ps 1.7 lt rgb "#00000000" notitle,\
"final_saccharification_data_5_micromolar.txt" u ($1*scale):2 w p pt 13 ps 1.7 lt rgb "#00000000" notitle,\
"final_saccharification_data_10_micromolar.txt" u ($1*scale):2 w p pt 9 ps 1.9 lt rgb "#00000000" notitle,\
[0:10]f_1(x) lw 5 lt rgb "#00000000" notitle,\
[8:19]f_2(x) lw 5 lt rgb "#00000000" notitle,\
[17:29]f_3(x) lw 5 lt rgb "#00000000" notitle,\
[39:53]f_5(x) lw 5 lt rgb "#00000000" notitle


#set origin .535, .27
#set size .23,.26
#clear


#set key top right font "bold, 8"
#unset key
#set xtics font ", 8"
#set ytics font ", 8"
#unset label
#set label "Untreated" at 10,240 rotate by -6 font "bold, 12"

#set label "Pretreated" at 35,545 rotate by -38 font "bold, 12"

#set yr[0:600]
#set xr[0:250]
#set xlabel "Lignin content \n (mg/g CWR)" font "bold, 12"
#set ylabel offset 3 "Total sugar released \n (mg/g CWR)" font "bold, 12"
#set bmargin 1
#set tmargin 1
#set lmargin 3
#set rmargin 1
#set title "Chen and Dixon (2007)" font "bold, 12"

#p[0:250] f_chen_untreated(x) lw 2 lt rgb "#BBBBBBBB" notitle,\
#[0:250] f_chen_pretreated(x) lw 2 lt rgb "#BBBBBBBB" notitle,\
#"Chen_dixon_mutants_untreated.csv" u 1:2 w p pt 4 ps 0.8 lt rgb "#000000" title "Untreated",\
#"Chen_dixon_mutants_pretreated.csv" u 1:2 w p pt 6 ps 0.8 lt rgb "#000000" title "Pretreated"
#unset multiplot
