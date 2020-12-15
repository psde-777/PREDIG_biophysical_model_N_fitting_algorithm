import numpy as np
import sys

names = ["0_percent_cry","10_percent_cry","30_percent_cry","50_percent_cry","80_percent_cry"]

crystallinities = [0,10,30,50,80]




with open("final_saccharification_" + sys.argv[1] + "_digestibility.txt", "w") as f:

    for i,name in enumerate(names):
        data = np.loadtxt(name + "/Output/mean_saccharification_medium.txt")
        f.write(str(crystallinities[i]) + "\t" + str(data[-1,1]) + "\n")
        
