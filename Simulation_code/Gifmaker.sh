#rm -r frames
mkdir frames

#unlink tempgnuplot.sh
touch tempgnuplot.sh

fmax=$(ls -d *Output/3D/visualisation_lignin* | wc -l)
fmax=$((fmax-1))


if [[ $fmax -lt 400 ]]; then
	skip=1
else
	skip=$((fmax/400))
fi




echo "Total number of snapshots : ${fmax}"

echo "For optimal viewing experience frame is taken every: ${skip} snapshot"

echo "set terminal pngcairo size 600, 1200" >> tempgnuplot.sh

echo "unset border" >> tempgnuplot.sh

echo "unset xtics" >> tempgnuplot.sh

echo "unset ytics" >> tempgnuplot.sh
echo "unset ztics" >> tempgnuplot.sh

echo "set xr[-5:9]" >> tempgnuplot.sh
echo "set yr[-4:8]" >> tempgnuplot.sh

echo "set key above" >> tempgnuplot.sh
echo "set key font ',24'" >> tempgnuplot.sh
echo "set key spacing 2" >> tempgnuplot.sh

echo "set zr[0:200]" >> tempgnuplot.sh

echo "j=0" >> tempgnuplot.sh

echo "do for [i=0:${fmax}:${skip}] { " >> tempgnuplot.sh

echo "j=j+1" >> tempgnuplot.sh

echo 'set output sprintf("frames/frame%05.0f.png",j)

  splot "Output/3D/visualisation_cellu_1_".i.".txt" pt 7 ps 1.5  lc rgb "green" title "cellulose",\
	"Output/3D/visualisation_hemi_1_".i.".txt" pt 5 ps 1.5 lc rgb "red" title "hemicellulose",\
	"Output/3D/visualisation_lignin_1_".i.".txt" pt 9 ps 1.5 lc rgb "blue" title "lignin"
}' >> tempgnuplot.sh 

echo "Gnuplot is creating frames"

gnuplot tempgnuplot.sh >/dev/null 2>&1

unlink tempgnuplot.sh

#convert -loop 0 frames/frame*.png digestion.gif

#ffmpeg -i frames/frame%5d.png -r 1 output.mp4

echo "Creating final animation"

ffmpeg -y -i frames/frame%5d.png -vf "format=rgb24,scale=120:-1,tile=10x10:color=black" -frames:v 1 tmp.png > /dev/null 2>&1
ffmpeg -y -i tmp.png -vf "palettegen=max_colors=64:reserve_transparent=0:stats_mode=single" palette.png > /dev/null 2>&1
ffmpeg -y -framerate 30 -i frames/frame%5d.png -i palette.png -filter_complex "[0]format=rgb24[b];[b][1]paletteuse=new=1" digestion.gif > /dev/null 2>&1

unlink tmp.png
unlink palette.png

rm -r frames