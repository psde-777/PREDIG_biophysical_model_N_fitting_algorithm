
set xtics font ", 20"
set ytics font ", 20"


set terminal postscript eps enhanced color

set key font ", 20"


scale = 1


set key top left
unset key


set xlabel "Time [arbitrary units]\n" font ", 24"
set ylabel "Glucan to glucose conversion (%)" font ", 26"

set xr[0:0.05]
set yr[-1:105]
directory_for_lign_plot=100

#set title "Comparison between simulated saccharification yields for various lignin percentages"
set output "lignin_saccharification_length_200.eps"


set label "No lignin" at 0.02,63 rotate by 44 font "bold, 20"
set label "10%" at 0.025,60 rotate by 45 font "bold, 20"
set label "20%" at 0.0295,48 rotate by 41 font "bold, 20"
set label "30%" at 0.035,25.5 rotate by 30 font "bold, 20"
set label "35%" at 0.038,12.5 rotate by 18 font "bold, 20"
set label "50%" at 0.04,2.5 rotate by 0 font "bold, 20"

p "length_200/mean_saccharification_0.txt" u ($1*scale):2 w l title "No lignin"lw 4 lt rgb "#00000000" ,\
"length_200/mean_saccharification_10.txt" u ($1*scale):2 w l title "10% lignin"lw 4 lt rgb "#00000000" dt ".",\
"length_200/mean_saccharification_20.txt" u ($1*scale):2 w l title "20% lignin"lw 4 lt rgb "#00000000" dt "-",\
"length_200/mean_saccharification_30.txt" u ($1*scale):2 w l title "30% lignin"lw 4 lt rgb "#00000000" dt "_",\
"length_200/mean_saccharification_35.txt" u ($1*scale):2 w l title "40% lignin"lw 4 lt rgb "#00000000" dt ".-",\
"length_200/mean_saccharification_50.txt" u ($1*scale):2 w l title "50% lignin"lw 4 lt rgb "#00000000" dt "._"


