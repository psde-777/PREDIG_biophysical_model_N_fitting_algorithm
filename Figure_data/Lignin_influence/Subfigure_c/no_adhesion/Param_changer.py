import numpy as np
import sys
percentages = ["0","10","15","20","25","30","40","43","47","48","49","50"]

percentage_dir_list = [percentages[i] + "_percent_lignin/" for i in range(len(percentages))]

kin_params = ["1","1","1","1",str(sys.argv[1]),"0.009","0.09","0.09","1","1"]

simu_params = [200000, 200000, 400000, 100000, 10, 1, 1, 2, 2, -1, 1, 1, -1, int(sys.argv[2]), 2.9]

for j,item in enumerate(percentage_dir_list):
    
    with open("enzymes_50/" + item + "Params/kinetic_parameters.txt", "w") as f:
        for param in kin_params:
            f.write(str(param) + "\t")
    with open("enzymes_50/" + item + "Params/simulation_parameters.txt", "w") as f:
        for param in simu_params:
            f.write(str(param) + "\t")