import numpy as np
import sys
import random

#argv[1] = delta around which to randomly change values
delta = 0.25#Example: if This is set to 0.1 it means that each value will be varied within a range of +- 10 percent of its current value
rundata = np.loadtxt("run_specs.txt")
Nvals = rundata.size
outdata = np.full((1,Nvals),0.)
for i in range(Nvals):
    if i < 5:
        rundata[i] += (random.random()-0.5)*2.*delta*rundata[i]
        if rundata[i] < 0:
            rundata[i] = 0
    outdata[0,i] = rundata[i]


np.savetxt("run_specs.txt", outdata, delimiter = "  ", fmt = "%1.9f")



