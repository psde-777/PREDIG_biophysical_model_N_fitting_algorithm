 

## Instructions for how to use the paramater optimization algorithm


### Setup directory and data to be fitted

- Create a directory in which you want to run the simulations (this will be referred to as topdir)

- copy the content of this directory into topdir

- If you have changed the simulation code, copy the new code into topdir/latest/, and run compile.sh via the terminal command

./compile.sh

Check that the "Params" directory, and the "Output" directory exist within the directory "latest"

- create an empty directory called "family_1" in topdir. If the directory already exists, make sure that it is empty

- Open the file "best_candidate_specs.txt" and enter those kinetic parameters from which the optimization algorithm should start. This file will be updated continuously during the algorithm

- IMPORTANT: Within topdir, there is a file "keywords.txt" (if not, you need to create it). From this file the algorithm takes
the names of the experimental data sets, and therefore it needs to contain AT LEAST ONE keyword.

- Make sure that for each keyword there exists a file "best_init_specs_keyword.txt" in topdir

- In the Output directory, create a directory "expe_data" if it does not exist yet. Within this directory,
 the data to be used for the fitting procedure needs to be placed.
 These data need to be contained in .txt files,
 in which the first column is time and the second column is the sugar release percentage.
 If you have data for glucose and xylose, these need to be contained in separate files.
 The end of the filename specifies whether it is a glucose or xylose data file, i.e. expe_saccharification_keyword_glc.txt or expe_saccharification_keyword_xyl_.txt.

- You need to have at least one experimental data file (either for glucose or xylose) for each keyword in keywords.txt. For example, if you hava set "data_1"
 and a set "data_2" for which you have glucose and xylose release data, then you need to have the following files in "expe_data":
     - "expe_saccharification_data_1_glc.txt" 
     - "expe_saccharification_data_1_xyl.txt"
     - "expe_saccharification_data_2_glc.txt"
     - "expe_saccharification_data_2_xyl.txt", 

- Open the files "kin_params_to_randomize.txt" and "init_params_to_randomize.txt". Here, you can choose, which parameters should be fitted by setting the respective value to 0 or 1

- Open the files "max_kin_vals" and "max_init_vals.txt". Here, you can specify the maximum value for each parameter.

### Running the algorithm

- Once the setup is complete, you may run the algorithm. For this, you need to specify the following five parameters:

    - the number of families **N_family** (currently to be kept at a value of 1)
    - the number of generations per family **N_gens** (choose an integer greater than or equal to 1)
    - the number of subsets per generation **N_subsets** (choose an integer greater than or equal to 1)
    - the percentage around which to vary the parameters **delta** (choose a number between 0 and 1)
    - the number of CPU cores you have available **N_cores** (choosing a number higher than the actual number of cores available will reduce the efficiency of the algorithm considerably)

- Now, to start the algorithm, execute the shell script "evo_all_in_one.sh" via the following command:

        ./evo_all_in_one.sh N_family N_gens N_subsets delta N_cores

- While the algorithm runs, the files "best_candidate_specs.txt" and "best_init_specs_keyword.txt" for each keyword in keywords.txt are updated if the difference between simulated and experimental data is reduced compared to the previous set of parameters in the file. Furthermore, the current lowest error is contained in the file "smallest_var.txt"

- If you want to visually compare the simulated and experimental data, run the code in a separate folder and use the parameters found in "best_candidate_specs.txt" in "kinetic_parameters.txt" (see Code README at the top of this repository for more detail)