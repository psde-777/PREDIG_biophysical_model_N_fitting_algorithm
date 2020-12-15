import numpy as np
import sys
dir_list = ["enzymes_50"]
#enzyme_conc_list = [50, 200, 100, 100]
#fibril_length_list = [200, 200, 100, 300]

#init_params = [100, 100, 100, 100, 200, 1, 0.05, 0.45, 0, 0.5, 0.5, 0.3]

percentage_dir_list = [str(10*i) for i in range(6)]

percentage_dir_list_1 = [0,10,15,20,25,30,40,50]
#percentage_dir_list_2 = [0,10,20,30,40,45,48,50]
#percentage_dir_list_3 = [0,10,20,30,40,45,48,50]
#percentage_dir_list_4 = [0,10,20,30,40,45,50]
#percentage_dir_list_5 = [0,10,20,25,30,35,40,50]
#percentage_dir_list_6 = [0,10,20,30,35,40,45,50]

percentage_dir_lists = [percentage_dir_list_1]#,percentage_dir_list_2,percentage_dir_list_3,percentage_dir_list_4,percentage_dir_list_5,percentage_dir_list_6]

time_limit = float(sys.argv[1])


for i,directory in enumerate(dir_list):
    with open("lignin_vs_final_yield_" + directory + ".txt","w") as f:
        for j,item in enumerate(percentage_dir_lists[i]):
            data = np.loadtxt("mean_saccharification_" + str(item) + ".txt")
            for k in range(data[:,0].size):
                if data[k,0] >= time_limit:
                    f.write(str(item) + "\t" + str(data[k,1]) + "\n")
                    break