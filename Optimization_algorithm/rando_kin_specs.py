from numpy import full,savetxt
import sys
import numpy as np
import sys
import random


def rnum():#return random number between -1 and 1
    return 2.*(random.random()-0.5)

#argv[1] = generation
#argv[2] = N_runs
#argv[3] = delta
#argv[4] = params_to_change
#argv[5] = parentfolder


generation=sys.argv[1]
N_runs=int(sys.argv[2])
delta=float(sys.argv[3])
current_params =sys.argv[4]
parent_folder = sys.argv[5]

N_unsuccessfull_tries = np.loadtxt("no_reduction_count.txt")
if N_unsuccessfull_tries > 2:
    delta += delta * N_unsuccessfull_tries/100.

#if delta > 1:
#    delta = 1

print("Delta = " + str(delta))
#if N_unsuccessfull_tries > 10:
#    delta *= 2
#if N_unsuccessfull_tries > 20:
#    delta *= 2
smallest_var = np.loadtxt("smallest_var.txt")
gradient_data = np.loadtxt("current_gradient.txt")
rundata = np.loadtxt(current_params)
Nvals = rundata.size
N_to_randomize = 5
max_val = 100.;

print("In rando_kinspecs.py: N_unsuccessful_tries = " + str(N_unsuccessfull_tries) + "; sum(gradient_data) = " + str(np.sum(gradient_data[:])))

if N_unsuccessfull_tries < 3 and np.sum(gradient_data[:]) != 0:
    print("Gradient param update");
    pct_old_new_mini = np.loadtxt("mini_difference_pct.txt");
    if pct_old_new_mini < 0.1:
        pct_old_new_mini = 1.
    if smallest_var < 1:
        pct_old_new_mini = 1
    outdata = np.full((1,Nvals),0.)
    formatvals = []
    for i in range(Nvals):
        formatvals.append("%1.8f")

    for j in range(int(N_runs)):
        grad_size = rnum();
        for i in range(Nvals):
            if i < N_to_randomize or i == Nvals -1 or i == Nvals -2:
                value_to_set = rundata[i] + grad_size*gradient_data[i]/pct_old_new_mini;
                if value_to_set > 0 and value_to_set <= max_val:
                    outdata[0,i] = value_to_set
                elif value_to_set == 0:
                    outdata[0,i] = delta * grad_size
                else:
                    print("Not changing param")
                    outdata[0,i] = rundata[i]
            else:
                outdata[0,i] = rundata[i]



        np.savetxt("family_" + parent_folder + "/Generation_" + str(generation) + "/Run_" + str(j+1) + "/Params/kinetic_parameters.txt", outdata, delimiter = "  ", fmt = formatvals)

else:

    print("Random param update (no gradient used)");
    outdata = np.full((1,Nvals),0.)

    formatvals = []
    for i in range(Nvals):
        formatvals.append("%1.8f")

    for j in range(int(N_runs)):

        for i in range(Nvals):
            if i < N_to_randomize or i == Nvals -1 or i == Nvals -2:
                if rundata[i] > 0:
                    value_to_set = rundata[i] + delta*rnum()*rundata[i]
                else:
                    value_to_set = delta*(rnum())
                if value_to_set > 0 and value_to_set <= max_val:
                    outdata[0,i] = value_to_set
                elif value_to_set == 0:
                    outdata[0,i] = delta * random.random()
                else:
                    outdata[0,i] = 0
            else:
                outdata[0,i] = rundata[i]



        np.savetxt("family_" + parent_folder + "/Generation_" + str(generation) + "/Run_" + str(j+1) + "/Params/kinetic_parameters.txt", outdata, delimiter = "  ", fmt = formatvals)





