import numpy as np

exp_data_fit_glc_low = np.loadtxt("fitdataCornStover_L_glc.txt")
exp_data_fit_xyl_low = np.loadtxt("fitdataCornStover_L_xyl.txt")
simu_data_low = np.loadtxt("mean_saccharification_low.txt")

exp_data_fit_glc_medium = np.loadtxt("fitdataCornStover_M_glc.txt")
exp_data_fit_xyl_medium = np.loadtxt("fitdataCornStover_M_xyl.txt")
simu_data_medium = np.loadtxt("mean_saccharification_medium.txt")

exp_data_fit_glc_high = np.loadtxt("fitdataCornStover_H_glc.txt")
exp_data_fit_xyl_high = np.loadtxt("fitdataCornStover_H_xyl.txt")
simu_data_high = np.loadtxt("mean_saccharification_high.txt")

simu_size_low = simu_data_low[0,:].size
simu_size_medium = simu_data_medium[0,:].size
simu_size_high = simu_data_high[0,:].size


exp_data_fit_glc_low[:,0] /= exp_data_fit_glc_low[-1,0]
exp_data_fit_xyl_low[:,0] /= exp_data_fit_xyl_low[-1,0]
exp_data_fit_glc_medium[:,0] /= exp_data_fit_glc_medium[-1,0]
exp_data_fit_xyl_medium[:,0] /= exp_data_fit_xyl_medium[-1,0]
exp_data_fit_glc_high[:,0] /= exp_data_fit_glc_high[-1,0]
exp_data_fit_xyl_high[:,0] /= exp_data_fit_xyl_high[-1,0]

#exp_data_fit[:,1] /= exp_data_fit[-1,1]

simu_data_low[:,0] /= simu_data_low[-1,0]
simu_data_medium[:,0] /= simu_data_medium[-1,0]
simu_data_high[:,0] /= simu_data_high[-1,0]
#simu_data[:,1] /= simu_data[-1,1]

if simu_size_low > 0:
    var_data_glc_low = (exp_data_fit_glc_low[:,1] - simu_data_low[:,1])*(exp_data_fit_glc_low[:,1] - simu_data_low[:,1])

if simu_size_low == 4:
    var_data_xyl_low = (exp_data_fit_xyl_low[:,1] - simu_data_low[:,3])*(exp_data_fit_xyl_low[:,1] - simu_data_low[:,3])
else:
    var_data_glc_low = np.full(exp_data_fit_glc_low[:,0].size,999999999)
    var_data_xyl_low = np.full(exp_data_fit_xyl_low[:,0].size,999999999)
if simu_size_medium > 0:
    var_data_glc_medium = (exp_data_fit_glc_medium[:,1] - simu_data_medium[:,1])*(exp_data_fit_glc_medium[:,1] - simu_data_medium[:,1])

if simu_size_medium == 4:
    var_data_xyl_medium = (exp_data_fit_xyl_medium[:,1] - simu_data_medium[:,3])*(exp_data_fit_xyl_medium[:,1] - simu_data_medium[:,3])
else:
    var_data_glc_medium = np.full(exp_data_fit_glc_medium[:,0].size,999999999)
    var_data_xyl_medium = np.full(exp_data_fit_xyl_medium[:,0].size,999999999)
if simu_size_high > 0:
    var_data_glc_high = (exp_data_fit_glc_high[:,1] - simu_data_high[:,1])*(exp_data_fit_glc_high[:,1] - simu_data_high[:,1])

if simu_size_high == 4:
    var_data_xyl_high = (exp_data_fit_xyl_high[:,1] - simu_data_high[:,3])*(exp_data_fit_xyl_high[:,1] - simu_data_high[:,3])
else:
    var_data_glc_high = np.full(exp_data_fit_glc_high[:,0].size,999999999)
    var_data_xyl_high = np.full(exp_data_fit_xyl_high[:,0].size,999999999)
#var_data = np.sqrt((exp_data_fit[:,:] - simu_data[:,:2])*(exp_data_fit[:,:] - simu_data[:,:2]))


#var_x = np.sum(var_data[:,0])/var_data[:,0].size
#var_y = np.sum(var_data[:,1])/var_data[:,1].size
#var_x = np.sum(var_data[:,0])
#var_y = np.sum(var_data[:,1])
#var = 0.5*(var_x + var_y)
var_glc = (sum(var_data_glc_low[:])/var_data_glc_low.size +  sum(var_data_glc_medium[:])/var_data_glc_medium.size +  sum(var_data_glc_high[:])/var_data_glc_high.size)/3
var_xyl = (sum(var_data_xyl_low[:])/var_data_xyl_low.size +  sum(var_data_xyl_medium[:])/var_data_xyl_medium.size +  sum(var_data_xyl_high[:])/var_data_xyl_high.size)/3
var = 0.5*(var_glc + var_xyl)
print("var_glc = " + str(var_glc) + "; var_xyl = " + str(var_xyl))
print("var = " + str(var))



outdata = np.full(1,var)
np.savetxt("var.txt", outdata, delimiter = "    ")
