# Modify 'Params/simulation_parameters.txt' file by replacing number in 12th column with 1
with open('Params/simulation_parameters.txt', 'r') as f:
    lines = f.readlines()

lines[0] = lines[0].split()
lines[0][11] = "1"
lines[0][4] = "32"
lines[0] = '    '.join(lines[0]) + '\n'

with open('Params/simulation_parameters.txt', 'w') as f:
    f.writelines(lines)
