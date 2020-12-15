import numpy as np
import sys
from scipy import optimize

def test_func(x, a, b):
    res = a*x/(x+b);
    return res;

#print("expe_data/expe_saccharification_" + str(sys.argv[1]) + "_glc.txt")

try:
    data = np.loadtxt("expe_data/expe_saccharification_" + str(sys.argv[1]) + "_glc.txt");
    simu_data = np.loadtxt("mean_saccharification_" + str(sys.argv[1]) + ".txt");
    N_points = simu_data[:,0].size;
    params_glc, params_covariance_glc = optimize.curve_fit(test_func, data[:,0], data[:,1]);
    fitdata_glc = np.full((N_points,2),0.);
    fitdata_glc[:,0] = simu_data[:,0];
    fitdata_glc[:,1] = test_func(simu_data[:,0], params_glc[0], params_glc[1]);
    #print("glc fitdata size = " + str(fitdata_glc[:,0].size) + "; simu_data size = " + str(simu_data[:,0].size));
    np.savetxt("fitdata_" + str(sys.argv[1]) + "_glc.txt", fitdata_glc, delimiter = "\t");
except IOError:
    pass;    
    #np.savetxt("fitdata_" + str(sys.argv[1] + "_glc.txt", np.full((10,2),99999999), delimiter = "\t");

try:
    data = np.loadtxt("expe_data/expe_saccharification_" + str(sys.argv[1]) + "_xyl.txt");
    simu_data = np.loadtxt("mean_saccharification_" + str(sys.argv[1]) + ".txt");
    N_points = simu_data[:,0].size;
    params_xyl, params_covariance_xyl = optimize.curve_fit(test_func, data[:,0], data[:,1]);
    fitdata_xyl = np.full((N_points,2),0.);
    fitdata_xyl[:,0] = simu_data[:,0];
    fitdata_xyl[:,1] = test_func(simu_data[:,0], params_xyl[0], params_xyl[1]);
    #print("xyl fitdata size = " + str(fitdata_xyl[:,0].size) + "; simu_data size = " + str(simu_data[:,0].size));
    np.savetxt("fitdata_" + str(sys.argv[1]) + "_xyl.txt", fitdata_xyl, delimiter = "\t");
except IOError:
    pass;
    #np.savetxt("fitdata_" + str(sys.argv[1] + "_xyl.txt", np.full((10,2),99999999), delimiter = "\t");
#simu_data_low = np.loadtxt("mean_saccharification_low.txt")
#simu_data_medium = np.loadtxt("mean_saccharification_medium.txt")
#simu_data_high = np.loadtxt("mean_saccharification_high.txt")
#N_low = simu_data_low[:,0].size
#N_medium = simu_data_medium[:,0].size
#N_high = simu_data_high[:,0].size

#data_glc_low = np.full((N_low,2),0.)
#data_xyl_low = np.full((N_low,2),0.)
#data_glc_medium = np.full((N_medium,2),0.)
#data_xyl_medium = np.full((N_medium,2),0.)
#data_glc_high = np.full((N_high,2),0.)
#data_xyl_high = np.full((N_high,2),0.)
#tmin = 0
#tmax = 70
#data_glc_low[:,0] = np.linspace(tmin, tmax, N_low)
#data_xyl_low[:,0] = np.linspace(tmin, tmax, N_low)
#data_glc_medium[:,0] = np.linspace(tmin, tmax, N_medium)
#data_xyl_medium[:,0] = np.linspace(tmin, tmax, N_medium)
#data_glc_high[:,0] = np.linspace(tmin, tmax, N_high)
#data_xyl_high[:,0] = np.linspace(tmin, tmax, N_high)

#for i in range(N_low):
#    data_glc_low[i,1] = glc_low(data_glc_low[i,0])
#    data_xyl_low[i,1] = xyl_low(data_xyl_low[i,0])
#for i in range(N_medium):
#    data_glc_medium[i,1] = glc_medium(data_glc_medium[i,0])
#    data_xyl_medium[i,1] = xyl_medium(data_xyl_medium[i,0])
#for i in range(N_high):
#    data_glc_high[i,1] = glc_high(data_glc_high[i,0])
#    data_xyl_high[i,1] = xyl_high(data_xyl_high[i,0])


#np.savetxt("fitdataCornStover_L_glc.txt", data_glc_low, delimiter = "   ")
#np.savetxt("fitdataCornStover_L_xyl.txt", data_xyl_low, delimiter = "   ")
#np.savetxt("fitdataCornStover_M_glc.txt", data_glc_medium, delimiter = "   ")
#np.savetxt("fitdataCornStover_M_xyl.txt", data_xyl_medium, delimiter = "   ")
#np.savetxt("fitdataCornStover_H_glc.txt", data_glc_high, delimiter = "   ")
#np.savetxt("fitdataCornStover_H_xyl.txt", data_xyl_high, delimiter = "   ")



