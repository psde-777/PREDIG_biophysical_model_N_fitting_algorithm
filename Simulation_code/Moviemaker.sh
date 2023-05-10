#!/bin/bash

# User selection mode 1 or mode 2
while true; do
	echo -e '\033[1m\033[32m Please select mode \033[0m'
	echo -e '\033[1m\033[35m Mode 1: \033[0m Make animation for a theoretical simulation \033[1m\033[31mNOT\033[0m comparing to any experimental data.'
	echo -e '\033[1m\033[35m Mode 2: \033[0m Make animation for a simulation \033[1m\033[31mWITH\033[0m comparison experimental data.'
	echo -e 'Enter selection [\033[1m\033[35m 1 \033[0m or \033[1m\033[35m 2 \033[0m].'
	echo ""
	read num
	echo ""

  if [ $num -eq 1 ]; then
  	echo -e 'You have chosen\033[1m\033[35m 1\033[0m : The Model simulation will run until all substrate has been digested (NO comparision).'	
    break
  elif [ $num -eq 2 ]; then
  	if [ -e Output/expe_data/expe_saccharification_movie_glc.txt ]; then
    	echo -e 'You have chosen\033[1m\033[35m 2\033[0m : The Model simulation will run till the maximum time for expermental data (WITH comparision).'
    else
    	echo -e '\033[1m\033[31m You have chosen 2. But, no experimental data is detected. Switching to Mode 1. \033[0m'
    	num=1
    fi	
    break
  else
  	echo -e '\033[1m\033[31m Error: invalid input. Try again. \033[0m'
    echo "........."
  fi
done


# copying the scripts from the folder to main directory
cp scripts_for_movie/frame_*.py .



# backing up simulation_parameters if mode 2 is selected
cp Params/simulation_parameters.txt  Params/simulation_parameters_backup.txt

# Changing simulation parameters.
if [ $num -eq 1 ]; then
	python3 frame_timelimit_1.py
elif [ $num -eq 2 ]; then	 
	python3 frame_timelimit_2.py
fi


# Running the simulation

rm -r Output/3D/* >> /dev/null

pyfiglet "Running your simulation"

./code_4 movie -vid >> /dev/null


# Making the frames

python3 frame_maker_1.py

if [ $num -eq 1 ]; then
	python3 frame_maker_2a.py
elif [ $num -eq 2 ]; then
	python3 frame_maker_2b.py
fi

python3 frame_joiner.py


# Set the input file pattern and the output video file name
input_pattern="snapshots/snap_*.png"
output_video="animation.mp4"

# Use FFmpeg to create the video
ffmpeg -framerate 30 -pattern_type glob -i "$input_pattern" -c:v libx264 -preset medium -crf 23 -pix_fmt yuv420p -y "$output_video"


# Clean up

pyfiglet "Cleaning up . . ."

unlink Params/simulation_parameters.txt
mv Params/simulation_parameters_backup.txt Params/simulation_parameters.txt
rm -r frames_py
rm -r snapshots
rm -r Output/3D/*
touch Output/3D/placeholder.txt
rm -r Output/enzyme_activity/*movie*.txt
rm -r Output/enzyme_concentration/*movie*.txt
rm -r Output/enzyme_fraction/enzyme_fraction_movie_*.txt
rm -r Output/saccharification/saccharification_movie_*.txt
rm -r frame_*.py
