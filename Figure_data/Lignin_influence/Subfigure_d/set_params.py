import numpy as np
import sys
import os

foldernames = np.loadtxt("directories.txt",dtype =str)
mus = np.loadtxt("mus.txt",dtype = float)
sigmas = np.loadtxt("sigmas.txt",dtype = float)

print(foldernames.shape)
print(mus.shape)
print(sigmas.shape)
count = 0
for k, foldername in enumerate(foldernames):
    for i,mu in enumerate(mus):
        for l,sigma in enumerate(sigmas):
            count += 1
            print("count = " + str(count))
            if os.path.isfile(foldername + "/folder_" + str(count) + "/Params/simulation_parameters.txt"):
                simu_params = np.loadtxt(foldername + "/folder_" + str(count) + "/Params/simulation_parameters.txt")
                with open(foldername + "/folder_" + str(count) + "/Params/simulation_parameters.txt", "w") as f:
                    for j,param in enumerate(simu_params):
                        if j == simu_params.size-2:#mu
                            f.write("%1.8f" % (mu) + "\t")
                        elif j == simu_params.size-1:#sigma
                            f.write("%1.8f" % (sigma) + "\t")
                        else:
                            if j < simu_params.size-3:
                                f.write("%d" % param + "\t")
                            else:
                                f.write("%1.8f" % param + "\t")
            else:
                exit()