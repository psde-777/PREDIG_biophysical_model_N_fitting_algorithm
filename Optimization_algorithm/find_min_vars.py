import numpy as np
import sys

def findmin(x,y):
    if x <= y:
        return x
    else:
        return y

#oldspecs = np.loadtxt("evo_specs.txt")
#newspecs = np.loadtxt("evo_specs.txt")

#argv[1] = family
#argv[2] = generation

vardata = np.loadtxt("vars.txt")
keywords = np.loadtxt("keywords.txt",dtype=str)
test_shape = np.full((1,2),0.)

#print("VARDATA SHAPE: " + str(vardata.shape))
noise_range = 0.;
if vardata.size == 2:
    mini_posi = np.full(1,1)
    np.savetxt("best_candidate.txt", mini_posi, delimiter = "\t", fmt = "%1i")    
else:
    mini_nr = vardata[0,0]
    mini = vardata[0,1]
    oldmini = mini

    for i in range(1,vardata[:,0].size):
        oldmini = mini
        mini = findmin(mini, vardata[i,1])
        if mini<oldmini:
            mini_nr = vardata[i,0]
    print(mini_nr)        


    previous_mini = np.loadtxt("smallest_var.txt")#This is the minimum var from the generation before
    count_unsuccessfull = np.loadtxt("no_reduction_count.txt")
    count_unsuccessfull_out = np.full(1,count_unsuccessfull)
    print("Old min variance: " + str(previous_mini))
    if mini < previous_mini:
        print ("Old mini = " + str(previous_mini) + "; new mini = " + str(mini))
        old_kin_data = np.loadtxt("best_kin_specs.txt");
        with open("old_best_kin_specs.txt", "w") as f:
            for item in old_kin_data:
                f.write("%1.8f\t" % item)
#        np.savetxt("old_best_kin_specs.txt", old_kin_data, delimiter = "  ", fmt = "%1.8f")
        with open("old_smallest_var.txt", "w") as f:
            f.write("%1.8f\t" % previous_mini)


        old_init_data = [np.loadtxt("best_init_specs_" + keyword + ".txt") for keyword in keywords]
        for j,keyword in enumerate(keywords):
            with open("old_best_init_specs_" + keyword + ".txt", "w") as output:
                for i, item in enumerate(old_init_data[j]):
                    if i < 5: 
                        output.write("%d\t" % item)
                    else:
                        output.write("%1.8f\t" % item)


        new_kin_data = np.loadtxt("family_" + str(sys.argv[1]) + "/Generation_" + str(sys.argv[2]) + "/Run_" + str(int(mini_nr)) + "/Params/kinetic_parameters.txt")
        new_kin_data_out = np.full((1,new_kin_data[:].size),0.)
        new_kin_data_out[0,:] = new_kin_data[:]
        np.savetxt("best_kin_specs.txt", new_kin_data_out, delimiter = "\t", fmt = "%1.8f")

        new_init_data = [np.loadtxt("family_" + str(sys.argv[1]) + "/Generation_" + str(sys.argv[2]) + "/Run_" + str(int(mini_nr)) + "/Params/initial_configuration_parameters_" + keyword + ".txt") for keyword in keywords]
        for j,keyword in enumerate(keywords):
            with open("best_init_specs_" + keyword + ".txt", "w") as output:
                for i, item in enumerate(new_init_data[j]):
                    if i < 5: 
                        output.write("%d\t" % item)
                    else:
                        output.write("%1.8f\t" % item)


        mini_data_out = np.full(1,mini)
        np.savetxt("smallest_var.txt", mini_data_out, delimiter = "\t", fmt = "%1.8f")
        mini_data_out[0] = int(sys.argv[2])
        np.savetxt("gen_of_smallest_var.txt", mini_data_out, delimiter = "\t", fmt = "%d")
        print( "mini_difference pct = " + str((previous_mini-mini)/previous_mini));
        if (previous_mini-mini)/previous_mini > noise_range:
            count_unsuccessfull = 0
    #        np.savetxt("no_reduction_count.txt",0)
        else:
            count_unsuccessfull += 1
    #        np.savetxt("no_reduction_count.txt",count_unsuccessfull)
    else:
        count_unsuccessfull +=1

    count_unsuccessfull_out[0] = count_unsuccessfull    
    np.savetxt("no_reduction_count.txt",count_unsuccessfull_out)
    #for i in range(1,vardata[:,0].size):
    #    if i != mini_nr:
    #        if abs(mini - vardata[i,1]) <= delta_check:
    #            count_around_delta += 1
    #            new_kin_data[:] += np.loadtxt("/mnt/data_scratch/eric/Generation_" + str(sys.argv[1]) + "/Run_" + str(int(i)) + "/Params/kinetic_parameters.txt")[:]
    #
    #new_kin_data /= (count_around_delta + 1)



    #print("======")
    #print(new_kin_data.size)
    #print(new_kin_data_out.size)
    #print("======")
    mini_posi_and_var = np.full(1,mini_nr)
    np.savetxt("best_candidate.txt", mini_posi_and_var, delimiter = "\t", fmt = "%1i")




