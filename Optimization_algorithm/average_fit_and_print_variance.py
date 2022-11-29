##### Added initial guess for test function parameters :partho
import numpy as np
import sys
import os.path
from scipy import optimize
#argv[1] = file path to output
#argv[2] = Name of mean file
#argv[3] = name of file that contains the keywords (e.g. low,medium,high)

def test_func(x, a, b):
    res = a*x/(x+b);
    return res;


i = 1;
end = False;




count = 0
countFound = 0
lastCount = 0#number of lines in previous file
max_size = 0
simu_params = np.loadtxt("latest/Params/simulation_parameters.txt");
N_files = int(simu_params[-4])
print("N_files = " + str(N_files))

#N_files = int(sys.argv[1]);#Number of files
#N_files = 20
var_glc_dict = {}
var_xyl_dict = {}
keywords = np.loadtxt("keywords.txt",dtype=str)

for keyword in keywords:
    countFound = 0
    expe_glc_file = sys.argv[1] + "expe_data/expe_saccharification_" + keyword + "_glc.txt";
    print(expe_glc_file)
    expe_xyl_file = sys.argv[1] + "expe_data/expe_saccharification_" + keyword + "_xyl.txt";
    print(expe_xyl_file)

    my_file = sys.argv[1] + "saccharification/saccharification_" + keyword + "_1" + ".txt";
    #if my_file.is_file(my_file):
    if os.path.isfile(my_file):
        countFound += 1
        data = np.loadtxt(my_file);
        Ncolums = data.shape[1]
        #data = np.loadtxt(sys.argv[2]+str(1)+".txt");
        max_size = data[:,0].size;
        min_size = max_size;
        final_time = data[-1,0]
        for i in range(2,N_files+1):
            my_file = sys.argv[1] + "saccharification/saccharification_" + keyword + "_" + str(i) + ".txt";
            if os.path.isfile(my_file):
                countFound += 1
                count = 0;
                data = np.loadtxt(my_file);
                count = data[:,0].size;
                max_size = max(count,max_size);
                min_size = min(count,min_size);
                final_time += data[-1,0]
                print(count)
                count = lastCount



        print("Max_size: ", max_size)
        final_time /= countFound
        print("final time: " + str(final_time))

        fitdata_glc = np.full((max_size,2),9999.)
        fitdata_xyl = np.full((max_size,2),9999.)
        expe_data_glc = np.full((10,2),9999.)
        expe_data_xyl = np.full((10,2),9999.)    
        if os.path.isfile(expe_glc_file):
            expe_data_glc = np.loadtxt(expe_glc_file)
            final_time = expe_data_glc[-1,0];
        else:
            print("file doesn't exist") 

        if os.path.isfile(expe_xyl_file):
            expe_data_xyl = np.loadtxt(expe_xyl_file)
    #        final_time = expe_data_xyl[-1,0];
        #else:
        #    print("No experimental data found!!!!")
        #    final_time = -1;
        #    exit();
    #    elif final_time > 100:
    #        final_time = 100;
        print("final time after if: " + str(final_time))
        outFile = np.full((max_size,Ncolums),0.)


        x = np.linspace(0,final_time,max_size)
        fitdata_glc[:,0] = x[:]
        fitdata_xyl[:,0] = x[:]
 

        params_glc, params_covariance_glc = optimize.curve_fit(test_func, expe_data_glc[:,0], expe_data_glc[:,1]);  #magic_numbers for initial guess ;)
        params_xyl, params_covariance_xyl = optimize.curve_fit(test_func, expe_data_xyl[:,0], expe_data_xyl[:,1]);  #magic_numbers for initial guess ;)
        fitdata_glc[:,1] = test_func(fitdata_glc[:,0],params_glc[0],params_glc[1])
        fitdata_xyl[:,1] = test_func(fitdata_xyl[:,0],params_xyl[0],params_xyl[1])

        np.savetxt(sys.argv[1]+'myFile-'+keyword+'.txt', fitdata_glc)  #### edit partho


        print(x[-1])

    #    outFile[:,0] = x[:];
        for i in range(N_files):
            my_file = sys.argv[1] + "saccharification/saccharification_" + keyword + "_" + str(i+1) + ".txt";
            if os.path.isfile(my_file):
                data = np.loadtxt(my_file);
                outFile[:,0] += x[:];
                for j in range(1,data[1,:].size):
                    outFile[:,j] += np.interp(x,data[:,0],data[:,j])
#                    np.savetxt(sys.argv[1]+'interpol-'+keyword+'.txt', outFile)  #### edit partho


        #        outFile[1:,:]/=N_files


        outFile[:,:]/=countFound
        np.savetxt(sys.argv[1]+'interpol-'+keyword+'.txt', outFile)  #### edit partho



        var_glc_dict[keyword] = np.sum((fitdata_glc[:,1] - outFile[:,1]) * (fitdata_glc[:,1] - outFile[:,1]))/max_size;
        var_xyl_dict[keyword] = np.sum((fitdata_xyl[:,1] - outFile[:,3]) * (fitdata_xyl[:,1] - outFile[:,3]))/max_size;
        


    #    np.savetxt(sys.argv[1] + sys.argv[2] + keyword + ".txt", outFile[:,:], delimiter="    ")
    else:
        with open(sys.argv[1] + "var.txt","w") as f:
            f.write(str(999999999999999))
            exit()
    #    np.savetxt(sys.argv[1] + sys.argv[2] + keyword + ".txt", outFile[:,:], delimiter="    ") 

var_glc = 0.
var_xyl = 0.
for keyword in keywords:
    var_glc += var_glc_dict[keyword]
    var_xyl += var_xyl_dict[keyword]

var_glc/=len(keywords)
var_xyl/=len(keywords)
var = 0.5*(var_glc+var_xyl)    ##comment if no Xyl but MLG
#var=var_glc     ##uncomment if only glucose relevant MLG or Xyl
    
print("var = " + str(var))
with open(sys.argv[1] + "var.txt","w") as f:
    f.write(str(var))
    exit()
