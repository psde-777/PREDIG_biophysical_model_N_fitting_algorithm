import numpy as np

N_data_points = np.loadtxt("latest/Params/kinetic_parameters.txt").size;
out_file = np.full((1,N_data_points),0.);
np.savetxt("current_gradient.txt", out_file, delimiter = "\t", fmt = "%1.8f")


