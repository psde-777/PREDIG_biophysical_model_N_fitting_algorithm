# Read experimental data for movie

with open('Output/expe_data/expe_saccharification_movie_glc.txt', 'r') as f:
    lines = f.readlines()
T_int = float(lines[-1].split()[0])



# Modify 'Params/simulation_parameters.txt' file by replacing number in 3rd column with 16, 5th column with (T_int+2), 12th column with 1
with open('Params/simulation_parameters.txt', 'r') as f:
    lines = f.readlines()

lines[0] = lines[0].split()
lines[0][3] = str(T_int+1)
lines[0][4] = "16"
lines[0][11] = "1"
lines[0] = '    '.join(lines[0]) + '\n'

with open('Params/simulation_parameters.txt', 'w') as f:
    f.writelines(lines)
