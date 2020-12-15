import numpy as np
import sys
import random

#argv[1] = delta around which to randomly change values
delta = 0#Example: if this is set to 0.1 it means that each value will be varied within a range of +- 10 percent of its current value
rundata = np.loadtxt("Params/simulation_parameters.txt")
Nvals = rundata.size
outdata = np.full((1,Nvals),0.)


formatvals = []

for i in range(Nvals):
    formatvals.append("%d")
    if i == Nvals-2 or i == Nvals-1:
        pass
        #rundata[i] = 1.
        #if rundata[i] < 0:
        #    rundata[i] = 0
    outdata[0,i] = rundata[i]

formatvals[Nvals-2] = "%1.4f"
formatvals[Nvals-1] = "%1.4f"



np.savetxt("Params/simulation_parameters.txt", outdata, delimiter = "  ", fmt = formatvals)



