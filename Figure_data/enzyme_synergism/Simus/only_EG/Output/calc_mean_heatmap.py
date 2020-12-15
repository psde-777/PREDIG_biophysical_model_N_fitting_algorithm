import numpy as np
#import matplotlib.pyplot as plt
#from qtbplot import plotting_context
#from scipy import optimize
import sys
import os.path
import random
#from matplotlib.colors import LogNorm
#from matplotlib import ticker, cm
#import pandas as pd


if len(sys.argv) < 6:
    print(len(sys.argv));
    print("Please provide the number of files, as well as the name of the input and output file, the number of time bins and, finally, the fibril length.");
    exit();

i = 1

print("Number of files: " + str(sys.argv[1]));
print("Name of first file: " + sys.argv[2]+str(i)+".txt");
print("Name of output file: " + sys.argv[3]);
print("Number of time bins: " + sys.argv[4]);
print("Fibril length: " + sys.argv[5]);



#20 DP_distrib/DP_distrib_low_ mean_DP_distrib_low.txt N_bins


############## Read in raw data ##############
all_data = []
print("Reading in data")
for i in range(1,int(sys.argv[1])+1):
    print(i)
    heatmap_data = np.loadtxt(sys.argv[2] + str(i) + ".txt")
    all_data.append(heatmap_data)
#    heatmap_plot(heatmap_data, 100)
print("Done reading")


####### Define time points to be considered
N_time_points = int(sys.argv[4])
max_time = all_data[0][-1,0];
print("max_time: " + str(max_time))
for data in all_data:
    print("max_time: " + str(data[-1,0]))
    max_time += data[-1,0]    
max_time /= float(sys.argv[1])

split_between_halfs = 0.5*max_time #Due to the mechanism of the gillespie algorithm, there are more points at the beginning than at the end (in time)

print("max time = " + str(max_time))
print("split between halfs = " + str(split_between_halfs))

time_points_first_half = np.linspace(0,split_between_halfs,int(0.5*N_time_points))
time_points_second_half = np.linspace(split_between_halfs,max_time,int(0.5*N_time_points))
time_points = np.full(N_time_points,0.);
time_points[:int(0.5*N_time_points)] = time_points_first_half[:]
time_points[int(0.5*N_time_points):] = time_points_second_half[:]



fibril_length = int(sys.argv[5]);
print("Building averaged file");
averaged_data = np.full((time_points.size*(fibril_length+1),4), 0.)# columns: time, DP, sum(data), Nbr of points in bin
for i in range(time_points.size):
    averaged_data[i*(fibril_length+1):(i+1)*(fibril_length+2),0] = time_points[i]
#    print(averaged_data[i*(fibril_length+1):(i+1)*(fibril_length+1),1].size)
#    print(np.linspace(0,fibril_length,fibril_length+1).size)
    averaged_data[i*(fibril_length+1):(i+1)*(fibril_length+1),1] = np.array(range(0,fibril_length+1))


set_count = 0
for data in all_data:
    set_count += 1
    print("File number" + str(set_count))
    for i in range(1,time_points.size):
#        print(i)
        current_bin_max = time_points[i]
        current_bin_min = time_points[i-1]
#        print(int(data[:,0].size/(fibril_length+1)))
#        print(averaged_data[:,0].size/(fibril_length+1))
        for j in range(int(data[:,0].size/(fibril_length+1))):
            current_time_point = data[j*(fibril_length),0]
            if current_time_point >= current_bin_min and current_time_point < current_bin_max:
#                print("sees")
#                print(averaged_data[j*(fibril_length+1):(j+1)*(fibril_length+1),2].shape)
#                print(data[j*(fibril_length+1):(j+1)*(fibril_length+1),2].shape)
#                print(.shape)
                averaged_data[i*(fibril_length+1):(i+1)*(fibril_length+1),2] += data[j*(fibril_length+1):(j+1)*(fibril_length+1),2];
                averaged_data[i*(fibril_length+1):(i+1)*(fibril_length+1),3] += 1;
for i in range(averaged_data[:,2].size):
    if averaged_data[i,3] != 0:
        averaged_data[i,2] /= averaged_data[i,3]

#for i in range(averaged_data[:,2].size):
#    if i > 0 and averaged_data[i,2] == 0:
            
bool_done = True
print("Saving new file");

np.savetxt(sys.argv[3], averaged_data[:,:3], delimiter = "\t", fmt="%1.8f")