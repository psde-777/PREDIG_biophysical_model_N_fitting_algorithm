
set xtics font ", 20"
set ytics font ", 20"


set terminal postscript eps enhanced color

set key font ", 20"


scale = 1
set key top left



set xlabel "Time [arbitrary units]\n" font ", 24"
set ylabel "Glucan to glucose conversion (%)" font ", 26"

set xr[0:0.05]
#set yr[-1:105]
directory_for_lign_plot=100

#set title "Comparison between simulated saccharification yields for various lignin percentages"



set key top right font ", 10"
set xlabel "Lignin content (%)\n"
unset xr
#set xr[0:55]

#unset yr

set xr[-1:60]
set yr[-1:60]


set ylabel "Final released glucose (%)" font ", 22"

#set key center right font ", 15"

scale_length_100=4848
scale_length_200=9648
scale_length_300=14448

f1(x) = a*x + b
f2(x) = c*x + d
f3(x) = e*x + f
f4(x) = g*x + h
f5(x) = i*x + j
f6(x) = k*x + l

f_studer_high_s_to_g_untreated(x) = m*x + n
f_studer_low_s_to_g_untreated(x) = o*x + p
f_studer_high_s_to_g_pretreated(x) = q*x + r
f_studer_low_s_to_g_pretreated(x) = s*x + t
f_chen_untreated(x) = -0.4062*x + 208.86
f_chen_pretreated(x) = -2.604*x + 695.39

f_linear(x) = 100-x


fit f1(x) "lignin_vs_final_yield_no_adhesion_no_structure.txt" u 1:2  via a,b
fit f2(x) "lignin_vs_final_yield_adhesion_no_structure.txt" u 1:2 every ::0::4 via c,d
fit f3(x) "lignin_vs_final_yield_no_adhesion_structure.txt" u 1:2 every ::6::11 via e,f
fit f4(x) "lignin_vs_final_yield_both_effects.txt" u 1:2 every ::0::3 via g,h

#fit f4(x) "lignin_vs_final_yield_s_g_less_than_2.txt" u 1:2 every ::0::3 via g,h



fit f_studer_high_s_to_g_untreated(x)"Studer_et_al_data_high_S_to_G_ratio_untreated.csv" u 1:2 via m,n
fit f_studer_low_s_to_g_untreated(x) "Studer_et_al_data_low_S_to_G_ratio_untreated.csv" u 1:2 via o,p
fit f_studer_high_s_to_g_pretreated(x) "Studer_et_al_data_high_S_to_G_ratio_180_degrees.csv" u 1:2 via q,r
fit f_studer_low_s_to_g_pretreated(x) "Studer_et_al_data_low_S_to_G_ratio_180_degrees.csv" u 1:2 via s,t


#set xr[15:30]




unset yr

set yr [-1:105]

scale = 1

set key top right font ", 20"


set output "final_yield_vs_lignin_chen_dixon_and_simus_with_hemi.eps"

set multiplot
unset key
set xr[-1:55]
set ylabel "Final released glucose (%)" font ", 22"
p "lignin_vs_final_yield_both_effects.txt" u 1:2 w p pt 5 ps 1.6 lt rgb "#000000" title "A on, S on",\
"lignin_vs_final_yield_adhesion_no_structure.txt" u 1:2 w p pt 3 ps 2 lt rgb "#000000" title "A on, S off",\
"lignin_vs_final_yield_no_adhesion_structure.txt" u 1:2 w p pt 1 ps 2.5 lt rgb "#000000" title "A off, S on",\
"lignin_vs_final_yield_no_adhesion_no_structure.txt" u 1:2 w p pt 7 ps 2 lt rgb "#000000" title "A off, S off",\
[0:53] f1(x) lw 5 lt rgb "#00000000" notitle,\
[0:42] f2(x) lw 5 lt rgb "#00000000" notitle,\
[38:52] f3(x) lw 5 lt rgb "#00000000" notitle,\
[0:35] f4(x) lw 5 lt rgb "#00000000" notitle
#f_linear(x) lw 5 lt rgb "#00000000" dt "." title "diagonal"


set origin .435, .25
set size .288,.288
clear


#set key top right font "bold, 8"
unset key
set xtics font ", 8"
set ytics font ", 8"

set label "Untreated" at 10,170 rotate by -6 font "bold, 12"

set label "Pretreated" at 40,550 rotate by -35 font "bold, 12"

set yr[0:600]
set xr[0:250]
set xlabel "Lignin content \n (mg/g CWR)" font "bold, 12"
set ylabel offset 3 "Total sugar released \n (mg/g CWR)" font "bold, 12"
set bmargin 1
set tmargin 1
set lmargin 3
set rmargin 1
set title "Chen et al (2007)" font "bold, 12"

p[0:250] f_chen_untreated(x) lw 2 lt rgb "#BBBBBBBB" notitle,\
[0:250] f_chen_pretreated(x) lw 2 lt rgb "#BBBBBBBB" notitle,\
"Chen_dixon_mutants_untreated.csv" u 1:2 w p pt 4 ps 0.8 lt rgb "#000000" title "Untreated",\
"Chen_dixon_mutants_pretreated.csv" u 1:2 w p pt 6 ps 0.8 lt rgb "#000000" title "Pretreated"
unset multiplot

