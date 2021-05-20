from numpy import full,savetxt
import sys
import numpy as np
import sys
import random


def rnum():#return random number between -1 and 1
    val = random.gauss(0,0.25)
    while val < -1 or val > 1:
        val = random.gauss(0,0.25)
    return val

#argv[1] = generation
#argv[2] = N_runs
#argv[3] = delta
#argv[4] = params_to_change
#argv[5] = parentfolder


def gradient():


    keywords = np.loadtxt("keywords.txt",dtype=str)



    #Generation N-2
    old_var = np.loadtxt("old_smallest_var.txt");
    old_kin_data = np.loadtxt("old_best_kin_specs.txt")#//
    old_init_data = [np.loadtxt("old_best_init_specs_" + keyword + ".txt") for keyword in keywords];

    #Generation N-1
    new_var = np.loadtxt("smallest_var.txt");
    new_kin_data = np.loadtxt("best_kin_specs.txt")
    new_init_data = [np.loadtxt("best_init_specs_" + keyword + ".txt") for keyword in keywords];

    #Maximum values of each parameter
    max_vals_kin = np.loadtxt("max_kin_vals.txt")
    max_vals_init = np.loadtxt("max_init_vals.txt")

    #Denotes which parameters should be varied
    variation_bools_kin = np.loadtxt("kin_params_to_randomize.txt",dtype = int)#Contains boolean values for which parameters are varied
    variation_bools_init = np.loadtxt("init_params_to_randomize.txt",dtype = int)

    N_kin_vals = new_kin_data.size
    N_init_vals = new_init_data[0].size


    kin_gradient = [0. for i in range(N_kin_vals)]
    for i in range(N_kin_vals):
        if new_kin_data[i] - old_kin_data[i] != 0:
            kin_gradient[i] = (new_var- old_var) / (new_kin_data[i]-old_kin_data[i])

    init_gradients = [np.full(N_init_vals,0.) for keyword in keywords]
    for i,keyword in enumerate(keywords):
        for j in range(N_init_vals):
            if new_init_data[i][j] - old_init_data[i][j] != 0:
                init_gradients[i][j] = ((new_var - old_var) / (new_init_data[i][j]-old_init_data[i][j]))



    out_kin_data = np.full(N_kin_vals,0.)
    out_init_data = [np.full(N_init_vals,0.) for keyword in keywords]


    for j in range(int(N_runs)):
        step_size = delta;
        for i in range(N_kin_vals):
            if variation_bools_kin[i] == 1:
                value_to_set = new_kin_data[i] - step_size * kin_gradient[i];
                if value_to_set > 0 and value_to_set <= max_vals_kin[i]:
                    out_kin_data[i] = value_to_set
                else:
                    print("Not changing kin param " + str(i) + " from " + str(new_kin_data[i])+ " to " + str(value_to_set))
                    out_kin_data[i] = new_kin_data[i]
            else:
                out_kin_data[i] = new_kin_data[i]
        with open("family_" + parent_folder + "/Generation_" + str(generation) + "/Run_" + str(j+1) + "/Params/kinetic_parameters.txt","w") as f:
            for k,data_point in enumerate(out_kin_data):
                if k == 4:
                    f.write("%d\t" % data_point)
                else:
                    f.write("%1.8f\t" % data_point)

        for i,keyword in enumerate(keywords):
            for k in range(N_init_vals):
                if variation_bools_init[k] == 1:
                    value_to_set = new_init_data[i][k] - step_size * init_gradients[i][k];
                    if value_to_set > 0 and value_to_set <= max_vals_init[k]:
                        out_init_data[i][k] = value_to_set
                    else:
                        print("Not changing init param " + str(k) + " from " + str(new_init_data[i][k])+ " to " + str(value_to_set))
                        out_init_data[i][k] = new_init_data[i][k]
                else:
                    out_init_data[i][k] = new_init_data[i][k]
            with open("family_" + parent_folder + "/Generation_" + str(generation) + "/Run_" + str(j+1) + "/Params/initial_configuration_parameters_" + keyword + ".txt","w") as f:
                for k,data_point in enumerate(out_init_data[i]):
                    if k < 6:
                        f.write("%d\t" % data_point)
                    else:
                        f.write("%1.8f\t" % data_point)
################################################################################################################


def random_sampling():


    keywords = np.loadtxt("keywords.txt",dtype=str)

    #Generation N-1
    new_var = np.loadtxt("smallest_var.txt");
    new_kin_data = np.loadtxt("best_kin_specs.txt")
    new_init_data = [np.loadtxt("best_init_specs_" + keyword + ".txt") for keyword in keywords];

    #Maximum values of each parameter
    max_vals_kin = np.loadtxt("max_kin_vals.txt")
    max_vals_init = np.loadtxt("max_init_vals.txt")

    #Denotes which parameters should be varied
    variation_bools_kin = np.loadtxt("kin_params_to_randomize.txt",dtype = int)#Contains boolean values for which parameters are varied
    variation_bools_init = np.loadtxt("init_params_to_randomize.txt",dtype = int)

    N_kin_vals = new_kin_data.size
    N_init_vals = new_init_data[0].size

    out_kin_data = np.full(N_kin_vals,0.)
    out_init_data = [np.full(N_init_vals,0.) for keyword in keywords]




    for j in range(int(N_runs)):
        for i in range(N_kin_vals):
            step_size = rnum()*delta;
            if variation_bools_kin[i] == 1:
                value_to_set = new_kin_data[i] + step_size * new_kin_data[i];
                if value_to_set > 0 and value_to_set <= max_vals_kin[i]:
                    out_kin_data[i] = value_to_set
                else:
                    print("Not changing kin param " + str(i) + " from " + str(new_kin_data[i])+ " to " + str(value_to_set))
                    out_kin_data[i] = new_kin_data[i]
            else:
                out_kin_data[i] = new_kin_data[i]
        with open("family_" + parent_folder + "/Generation_" + str(generation) + "/Run_" + str(j+1) + "/Params/kinetic_parameters.txt","w") as f:
            for k,data_point in enumerate(out_kin_data):
                if k == 4:
                    f.write("%d\t" % data_point)
                else:
                    f.write("%1.8f\t" % data_point)

        for i,keyword in enumerate(keywords):
            for k in range(N_init_vals):
                step_size = rnum()*delta;
                if variation_bools_init[k] == 1:
                    value_to_set = new_init_data[i][k] + step_size * new_init_data[i][k];
                    if value_to_set > 0 and value_to_set <= max_vals_init[k]:
                        out_init_data[i][k] = value_to_set
                    else:
                        print("Not changing init param " + str(k) + " from " + str(new_init_data[i][k])+ " to " + str(value_to_set))
                        out_init_data[i][k] = new_init_data[i][k]
                else:
                    out_init_data[i][k] = new_init_data[i][k]
            with open("family_" + parent_folder + "/Generation_" + str(generation) + "/Run_" + str(j+1) + "/Params/initial_configuration_parameters_" + keyword + ".txt","w") as f:
                for k,data_point in enumerate(out_init_data[i]):
                    if k < 6:
                        f.write("%d\t" % data_point)
                    else:
                        f.write("%1.8f\t" % data_point)
################################################################################################################


generation=sys.argv[1]
N_runs=int(sys.argv[2])
delta=float(sys.argv[3])
#current_params =sys.argv[4]
parent_folder = sys.argv[4]

N_unsuccessful_tries = np.loadtxt("no_reduction_count.txt")
if N_unsuccessful_tries > 2:
    delta += delta * N_unsuccessful_tries/10.


print("Delta = " + str(delta))


print("In rando_specs.py: N_unsuccessful_tries = " + str(N_unsuccessful_tries))

print("Generation = " + str(generation))
if N_unsuccessful_tries < 1 and int(generation) != 1:
    print("Gradient param update");
    gradient()


else:

    print("Random param update (no gradient used)");
    random_sampling()





