import numpy as np
import sys


#foldernames = np.loadtxt("foldernames.txt",dtype =str)

foldernames = np.loadtxt("directories.txt",dtype = str)
crystallinities = np.loadtxt("percentages.txt",dtype = str)
#names = ["0_percent_cry","10_percent_cry","30_percent_cry","50_percent_cry","80_percent_cry"]

#crystallinities = [0,10,30,50,80]

final_time = float(sys.argv[1])


for foldername in foldernames:
    with open("final_saccharification_" + foldername + ".txt", "w") as f:

        for i,name in enumerate(crystallinities):
            data = np.loadtxt(foldername + "/mean_saccharification_" + name + ".txt")
            written = False
            for j in range(data[:,0].size):
                if data[j,0] >= final_time:
                    f.write(str(crystallinities[i]) + "\t" + str(data[j,1]) + "\n")
                    written = True
                    break
            if written == False:
                f.write(str(crystallinities[i]) + "\t" + str(data[-1,1]) + "\n")
