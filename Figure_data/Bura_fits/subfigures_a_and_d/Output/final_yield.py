import numpy as np


N = 7

data6p25 = np.loadtxt("mean_saccharification_6p25_perc_lignin.txt")
data12p5 = np.loadtxt("mean_saccharification_12p5_perc_lignin.txt")
data18p75 = np.loadtxt("mean_saccharification_18p75_perc_lignin.txt")
data24 = np.loadtxt("mean_saccharification_24_perc_lignin.txt")
data31p25 = np.loadtxt("mean_saccharification_31p25_perc_lignin.txt")
data37p5 = np.loadtxt("mean_saccharification_37p5_perc_lignin.txt")
data50 = np.loadtxt("mean_saccharification_50_perc_lignin.txt")

final_yields = np.full(N,0.)



final_yields[0] = data6p25[-1,1]
final_yields[1] = data12p5[-1,1]
final_yields[2] = data18p75[-1,1]
final_yields[3] = data24[-1,1]
final_yields[4] = data31p25[-1,1]
final_yields[5] = data37p5[-1,1]
final_yields[6] = data50[-1,1]


np.savetxt("final_yields.txt",final_yields, delimiter = "	") 
