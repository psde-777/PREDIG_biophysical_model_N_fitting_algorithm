import glob
import os

# Find the latest file matching the pattern 'latest/Output/expe_data/expe_saccharification_*.txt'
file_list = glob.glob('latest/Output/expe_data/expe_saccharification_*.txt')
latest_file = max(file_list, key=os.path.getctime)

# Read the last element of the first column of the latest file and convert it to an integer
with open(latest_file, 'r') as f:
    T_int = int(f.readlines()[-1].split()[0])

# Modify 'Params/simulation_parameters.txt' file by replacing number in fourth column with (T_int+1)
with open('latest/Params/simulation_parameters.txt', 'r') as f:
    lines = f.readlines()

lines[0] = lines[0].split()
lines[0][3] = str(T_int+2)
lines[0] = '    '.join(lines[0]) + '\n'

with open('latest/Params/simulation_parameters.txt', 'w') as f:
    f.writelines(lines)
