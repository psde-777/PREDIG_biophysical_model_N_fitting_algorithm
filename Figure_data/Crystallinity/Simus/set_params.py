import numpy as np
import sys


foldernames = np.loadtxt("foldernames.txt",dtype =str)
percentages = np.loadtxt("percentages.txt",dtype = str)
crystallinities = [0.,10.,20.,30.,40.,50.,60.,70.,80.,90.]




for foldername in foldernames:
    for i,name in enumerate(percentages):
        init_params = np.loadtxt(foldername + "/" + name + "_percent_cry/Params/initial_configuration_parameters.txt")
        with open(foldername + "/" + name + "_percent_cry/Params/initial_configuration_parameters.txt", "w") as f:
            for j,param in enumerate(init_params):
                if j == init_params.size - 2 or j == init_params.size - 3:
                    f.write("%1.8f" % (crystallinities[i]/100.) + "\t")
                else:
                    if j < 5:
                        f.write("%d" % param + "\t")
                    else:
                        f.write("%1.8f" % param + "\t")
                
