#!/bin/bash

rm -r Output/3D/* >> /dev/null

pyfiglet "Running your simulation"

./code_4 movie -vid >> /dev/null



python3 frame_maker.py

# Set the input file pattern and the output video file name
input_pattern="frames_py/frame_*.png"
output_file="animation.mp4"

# Use FFmpeg to create the video
ffmpeg -framerate 30 -pattern_type glob -i "$input_pattern" -c:v libx264 -preset medium -crf 23 -pix_fmt yuv420p -y "$output_file"

pyfiglet "Cleaning up . . ."

rm -r frames_py
rm -r Output/3D/*
touch Output/3D/placeholder.txt
rm -r Output/enzyme_activity/*movie*.txt
rm -r Output/enzyme_concentration/*movie*.txt
rm -r Output/enzyme_fraction/enzyme_fraction_movie_*.txt
rm -r Output/saccharification/saccharification_movie_*.txt
