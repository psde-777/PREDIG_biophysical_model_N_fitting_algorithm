import numpy as np
import sys


#foldernames = np.loadtxt("foldernames.txt",dtype =str)

foldernames = np.loadtxt("foldernames.txt",dtype = str)
crystallinities = np.loadtxt("percentages.txt",dtype = str)
#names = ["0_percent_cry","10_percent_cry","30_percent_cry","50_percent_cry","80_percent_cry"]

#crystallinities = [0,10,30,50,80]



for foldername in foldernames:
    with open("final_saccharification_" + foldername + ".txt", "w") as f:

        for i,name in enumerate(crystallinities):
            data = np.loadtxt(foldername + "/" + name + "_percent_cry/Output/mean_saccharification.txt")
            f.write(str(crystallinities[i]) + "\t" + str(data[-1,1]) + "\n")
        
