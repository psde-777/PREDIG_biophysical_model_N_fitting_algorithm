import numpy as np

exp_data_fit_glc = np.loadtxt("fitdataCornStover_M_glc.txt")
exp_data_fit_xyl = np.loadtxt("fitdataCornStover_M_xyl.txt")
simu_data = np.loadtxt("mean_saccharification.txt")

simu_size = simu_data[0,:].size

exp_data_fit_glc[:,0] /= exp_data_fit_glc[-1,0]
exp_data_fit_xyl[:,0] /= exp_data_fit_xyl[-1,0]
#exp_data_fit[:,1] /= exp_data_fit[-1,1]
simu_data[:,0] /= simu_data[-1,0]
#simu_data[:,1] /= simu_data[-1,1]

if simu_size > 0:
    var_data_glc = (exp_data_fit_glc[:,1] - simu_data[:,1])*(exp_data_fit_glc[:,1] - simu_data[:,1])

if simu_size == 4:
    var_data_xyl = (exp_data_fit_xyl[:,1] - simu_data[:,3])*(exp_data_fit_xyl[:,1] - simu_data[:,3])
else:
    var_data_glc = np.full(exp_data_fit_glc[:,0].size,999999999)
    var_data_xyl = np.full(exp_data_fit_xyl[:,0].size,999999999)
#var_data = np.sqrt((exp_data_fit[:,:] - simu_data[:,:2])*(exp_data_fit[:,:] - simu_data[:,:2]))


#var_x = np.sum(var_data[:,0])/var_data[:,0].size
#var_y = np.sum(var_data[:,1])/var_data[:,1].size
#var_x = np.sum(var_data[:,0])
#var_y = np.sum(var_data[:,1])
#var = 0.5*(var_x + var_y)
var_glc = sum(var_data_glc[:])/var_data_glc.size
var_xyl = sum(var_data_xyl[:])/var_data_xyl.size


print("Optimizing only for xyl; var_xyl = " + str(var_xyl) + "(; var_glc = " + str(var_glc) + ")" )
var = var_xyl


outdata = np.full(1,var)
np.savetxt("var.txt", outdata, delimiter = "    ")
