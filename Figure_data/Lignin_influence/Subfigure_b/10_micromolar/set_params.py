import numpy as np
import sys
import os

foldernames = np.loadtxt("directories.txt",dtype =str)
percentages = np.loadtxt("percentages.txt",dtype = str)
percentages_float = np.loadtxt("percentages.txt",dtype = float)

print(foldernames.shape)
print(percentages.shape)

for k, foldername in enumerate(foldernames):
    for i,name in enumerate(percentages):
        init_params = np.loadtxt(foldername + "/" + name + "_percent_lignin/Params/initial_configuration_parameters.txt")
        with open(foldername + "/" + name + "_percent_lignin/Params/initial_configuration_parameters.txt", "w") as f:
            for j,param in enumerate(init_params):
                if j == 6:#hemi fraction
                    if 0.5-(percentages_float[i]/100) >= 0:
                        f.write("%1.8f" % (0.5 - (percentages_float[i]/100)) + "\t")
                    else:
                        f.write("%d" % 0 + "\t")
                elif j == 7:#lign fraction
                    f.write("%1.8f" % (percentages_float[i]/100) + "\t")
                else:
                    if j < 5:
                        f.write("%d" % param + "\t")
                    else:
                        f.write("%1.8f" % param + "\t")
                
