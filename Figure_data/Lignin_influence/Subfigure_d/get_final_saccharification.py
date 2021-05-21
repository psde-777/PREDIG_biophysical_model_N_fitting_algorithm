import numpy as np
import sys

final_time = float(sys.argv[1])

foldernames = np.loadtxt("directories.txt",dtype =str)
mus = np.loadtxt("mus.txt",dtype = float)
sigmas = np.loadtxt("sigmas.txt",dtype = float)

print(foldernames.shape)
print(mus.shape)
print(sigmas.shape)
count = 0
for k, foldername in enumerate(foldernames):
    with open("final_saccharification_" + foldername + ".txt", "w") as f:
        for i,mu in enumerate(mus):
            for l,sigma in enumerate(sigmas):
                count += 1
                print(count)
                data = np.loadtxt(foldername + "/mean_saccharification_" + str(count) + ".txt")
                written = False
                for j in range(data[:,0].size):
                    if data[j,0] >= final_time:
                        f.write(str(mu) + "\t" + str(sigma) + "\t" + str(data[j,1]) + "\n")
                        written = True
                        break
                if written == False:
                    f.write(str(mu) + "\t" + str(sigma) + "\t" + str(data[-1,1]) + "\n")
