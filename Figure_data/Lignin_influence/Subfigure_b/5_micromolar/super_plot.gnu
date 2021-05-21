
set xtics font ", 20"
set ytics font ", 20"


set terminal postscript eps enhanced color

set key font ", 20"


scale = 1


set key top left
unset key


set xlabel "Time [arbitrary units]\n" font ", 24"
set ylabel "Glucan to glucose conversion (%)" font ", 26"

set xr[0:40]
set yr[-1:110]
directory="data"

#set title "Comparison between simulated saccharification yields for various lignin percentages"
set output "lignin_saccharification_length_200.eps"


set label "No lignin" at 6,72 rotate by 44 font ", 20"
set label "20%" at 10,70 rotate by 45 font ", 20"
set label "30%" at 12,65 rotate by 35 font ", 20"
set label "40%" at 15,56 rotate by 27 font ", 20"
set label "50%" at 20,24 rotate by 15 font ", 20"

p directory."/mean_saccharification_0.txt" u ($1*scale):2 w l title "No lignin"lw 4 lt rgb "#00000000" ,\
directory."/mean_saccharification_20.txt" u ($1*scale):2 w l title "20% lignin"lw 4 lt rgb "#00000000" dt ".",\
directory."/mean_saccharification_30.txt" u ($1*scale):2 w l title "30% lignin"lw 4 lt rgb "#00000000" dt "-",\
directory."/mean_saccharification_40.txt" u ($1*scale):2 w l title "40% lignin"lw 4 lt rgb "#00000000" dt "_",\
directory."/mean_saccharification_50.txt" u ($1*scale):2 w l title "50% lignin"lw 4 lt rgb "#00000000" dt ".-"

#directory."/mean_saccharification_35.txt" u ($1*scale):2 w l title "35% lignin"lw 4 lt rgb "#00000000" dt "_",\