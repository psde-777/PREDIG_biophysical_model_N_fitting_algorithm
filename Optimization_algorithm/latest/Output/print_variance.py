import numpy as np
import sys



var = 0.;
var_glc = 0.;
var_xyl = 0.;
count_xyl_hit = 0;
count_glc_hit = 0;
if len(sys.argv) > 1:
    for arg in sys.argv[1:]:
        simu_data = np.loadtxt("mean_saccharification_" + str(arg)+ ".txt");#the time values for the simu and expe data should be the same due to the way the expe data are created in fitfunc.py
        try:
            fitted_data_glc = np.loadtxt("fitdata_" + str(arg) + "_glc.txt");
            #print("Difference in length = " + str(fitted_data_glc[:,0].size - simu_data[:,1].size))
            count_glc_hit += 1
            var_glc += np.sum((fitted_data_glc[:,1] - simu_data[:,1]) * (fitted_data_glc[:,1] - simu_data[:,1]))/simu_data[:,1].size;
        except IOError:
            var_glc += 0
        try:
            fitted_data_xyl = np.loadtxt("fitdata_" + str(arg) + "_xyl.txt");
            count_xyl_hit += 1
            var_xyl += np.sum((fitted_data_xyl[:,1] - simu_data[:,3])*(fitted_data_xyl[:,1] - simu_data[:,3]))/simu_data[:,3].size;
        except IOError:
            var_xyl += 0

if count_glc_hit > 0 or count_xyl_hit > 0:
    if count_glc_hit > 0 and count_xyl_hit > 0:
        var_glc /= count_glc_hit;
        var_xyl /= count_xyl_hit;
        var = 0.5*(var_glc + var_xyl);
    else:
        if count_glc_hit > 0:
            var = var_glc/count_glc_hit;
        elif count_xyl_hit > 0:
            var = var_xyl/count_xyl_hit;
else:
    var = 9999999999


outdata = np.full(1,var)
np.savetxt("var.txt", outdata, delimiter = "    ")



#############################################################################################################################





