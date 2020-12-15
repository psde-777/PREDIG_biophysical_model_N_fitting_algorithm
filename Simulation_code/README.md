# Simulated saccharification of a lignocellulose microfibril

This README provides instructions on how to run the simulation code, and a brief documentation of its structure.

## TLDR

You can run the code by executing code_4 via the terminal command

        ./code_4

## Compiling, setup and first run

(You can skip the first step, if you already have a file named "code_4")

To run the code, first do the following:

- Compile the source code via the terminal command 

        g++ Code/code_4.cpp Code/functions.cpp -o code_4

  and copy the code executable "code_4" to the directory from which you want to run it

- Alternatively, copy the already existing executable "code_4"


- Next, copy the folders "Params" and "Output" into the same directory as the executable

- Make sure that the "Params" folder contains at least the following files:

    - "simulation_parameters.txt"
    - "kinetic_parameters.txt"
    - "initial_configuration_parameters.txt"
    - "structural_parameters_1.txt"
    - "hydro_parameters_inner_2_1.txt"
    - "hydro_parameters_outer_2_1.txt"


- Make sure that the "Output" folder contains the following subfolders:

    - "3D"
    - "DP_distrib"
    - "enzyme_activity"
    - "enzyme_fraction"
    - "enzyme_concentration"
    - "Nbr_reactions"
    - "saccharification"


- You can now run the code via the terminal command

    ./code_4


## Customizing the parameters

    
- You can adjust the parameters within "simulation_parameters.txt", "kinetic_parameters.txt" and "initial_configuration_parameters.txt" in
 the "Params" folder to suit your requirements. The parameters are ordered as follows:
 
### simulation_parameters:

1.) and 2.) sum up to yield the maximum number of steps (there are two parameters for this in case of possible future added effects).
    Example: the first parameter is set to 10000, the second is set to 5000. Then, no more than 15000 Gillespie algorithm steps can be performed

3.) Maximal number of snapshots which may be taken for the purpose of creating animation data.

4.) Maximal value for the real time within the simulation, as calculated by the Gillespie algorithm.

5.) Currently unused, but needs to be set to a positive integer value

6.) Number of steps between system snapshots, if "-vid" command is used (see below). If this is set to 1, a snapshot is taken at every step of the Gillespie algorithm.

7.) mode_code: Determines the shape of the microfibril. Currently restricted to a value of 1

8.) mode_hemi: Determines, on which sides of the microfibril hemicellulose can anchor

9.) mode_lign: Determines, on which sides of the microfibril lignin can anchor

10.) mode_inhib: Determines, whether inhibition is active (1) or inactive (-1)

11.) mode_lignin_glue: Determines, whether the lignin gluing effect is active (1) or inactive (-1)

12.) Currently unused, but needs to be set to a positive integer value

13.) mode_enzyme_size : Determines, whether the enzyme size is included (1) or inactive (-1).

14.) Nruns: Determines the number of simulation runs

15.) enzyme radius: Determines the radius of the enzymes



### kinetic_parameters:

1.) EG rate of reaction

2.) CBH rate of reaction

3.) BGL rate of reaction

4.) XYL rate of reaction

5.) Lignin rate of adhesion

6.) binding affinity of cellobiose to EG

7.) binding affinity of cellobiose to CBH

8.) binding affinity of glucose to BGL

9.) ratio between digestibility of "crystalline" and "amorphous" cellulose

10.) ratio between digestibility of "crystalline" and "amorphous" hemicellulose



### initial_configuration_parameters:
 
1.) initial EG concentration

2.) initial CBH concentration

3.) initial BGL concentration

4.) initial XYL concentration

5.) Length of the microfibril in bonds 

6.) Percentage of xylose in hemicellulose; currently restricted to a value of 1

7.) Percentage of hemicellulose within the microfibril

8.) Percentage of lignin within the microfibril

9.) Percentage of acetylated hemicellulose within the hemicellulose fraction; currently restricted to a value of 1

10.) Percentage of "crystalline" cellulose

11.) Percentage of "crystalline" hemicellulose

12.) Currently unused, set to 0.3 per default



## Running the customized code

- As above, you may run the code via the terminal command

    ./code_4


- There are several additional command line arguments which may be added. For command line arguments X1, X2, ... Xn, the code would be run as follows:

    ./code_4 X1 X2 ... Xn
    
    **IMPORTANT: the very first command line argument determines, from which file the initial configuration parameters are read. Specifically this means that for the command**

    ./code_4 example
    
    **the file "initial_parameters_example.txt" is loaded. This also leads to a change in the name of the output files:**
    **the files in the saccharification folder for example will be named "saccharification_example_1.txt", "saccharification_example_2.txt", etc. (the number represents the number of the individual simulation run).**
    **If there are no command line arguments, the default "initial_configuration_parameters.txt" file is loaded instead, and the output files will simply be called "saccharification_1.txt", "saccharification_2.txt", etc.**

    The currently available command line arguments are as follows:
    

    "-verbose": This will cause more detailed information of what is happening during the individual simulation runs to be written to the terminal.

    "-heatmap": this will enable the Output of files documenting the distribution of glucose within cellulose polymers of all lengths within the "DP_distrib" folder.
    As this is somewhat computationally expensive, it was chosen as an optional component

    "-vid": This will cause the first simulation run to document the 3D positioning of all cellulose, hemicellulose and lignin bonds at step intervals, which are determined by a parameter in "simulation_parameters.txt" (see above).
    They are saved in the folder "3D" and can be used to create a visualization of the microfibril digestion

    "-fixed_seed": This will fix the seed of the pseudo-random number generator used within the simulation. The seed which will be set is currently specified within the source code, in the file "code_4.cpp"
    
    "-suppress_output": This will suppress all data output from the simulations. It is mainly used for debugging
    
    "-timer": This leads to the output of files which document the the average computation time taken for the individual enzyme reactions. This is still in development

        


## Dealing with the output

After the code has finished running, it will have generated data in the subfolders of the "Output" folder. This will be a single file per subfolder, if the number of runs in "simulation_parameters.txt" is set to 1.
Otherwise it will be multiple files. There are several types of output data, which are summarized here first, before providing instructions on how to average multiple raw data files.

### Types of output data

#### Saccharification

These files track the saccharification data over one simulation and contain 4 columns of data:

- First column: simulated time
- second column: Percentage of released glucose
- third column: Percentage of glucose contained in cellobiose
- fourth column: Percentage of released xylose

#### DP_distrib

These files track the entire distribution of polymer lengths over the course of one simulation and contain 3 columns of data:

- First column: simulated time
- second column: Polymer length
- third column: Percentage of glucose contained in the respective polymer length

Keep in mind that there are multiple values for the same time point. Therefore the ordering of the rows is also important here:
For each time point, there are

N = L_fibril + 1

rows, each containing the respective values for polymer lengths from 0 to L_fibril respectively.

#### enzyme_activity

These files track the fraction of previously occoured Gillespie reactions associated to EG, CBH, BGL, XYL and lignin adhesion over the course of a simulation and contain 6 columns of data:

- First column: simulated time
- second column: EG fraction
- third column: CBH fraction
- fourth column: BGL fraction
- fifth column: XYL fraction
- sixth column: lignin adhesion fraction

#### enzyme_fraction

These files track the fraction within the propensity table associated to EG, CBH, BGL, XYL and lignin adhesion over the course of a simulation and countain 6 column of data:

- First column: simulated time
- second column: EG fraction
- third column: CBH fraction
- fourth column: BGL fraction
- fifth column: XYL fraction
- sixth column: lignin adhesion fraction

#### enzyme_concentration

These files track the enzyme concentration over the course of a simulation (this may change due to adhesion to lignin) and contain 5 columns:

- First column: simulated time
- second column: EG concentration
- third column: CBH concentration
- fourth column: BGL concentration
- fifth column: XYL concentration

#### 3D

Each of these files contain the coordinates of all bonds of the respective type (cellulose, hemicellulose, lignin) at a single time point during a simulation and contain 3 columns:

- First column: x-coordinates
- second column: y-coordinates
- third column: z-coordinates


### Averaging none-3D data and none-heatmap_data

Within the "Output" folder there exists a python script "calc_mean_interpolate.py", which is used for calculating the average over all files within a subfolder.
It was tested using Python version 3.7 and is executed via the terminal command

    python3 calc_mean_interpolate.py N_files common_name outputfilename location_of_experimental_data

Here, 

- N_files denotes the number of files over which the average is calculated

- common_name denotes the file location up to the name common to all files

- outputfilename denotes the name of the averaged file

- location_of_experimental_data denotes the location of possible experimental data. This determines the simulated time up to which the average is calculated. If you do not want this functionality, simply add any word which is not a file path here

As an example, this would be the command used to calculate the average over 100 saccharification data files:

    python3 calc_mean_interpolate.py 100 saccharification/saccharification_ mean_saccharification.txt expe_data/expe_saccharification_

Alternatively, the raw data can also be used directly, depending on preference.

### Averaging heatmap data

There is another python script within the "Output" folder, which ist called "calc_mean_heatmap.py". Use it via the terminal command

    python3 calc_mean_heatmap.py N_files common_name outputfilename N_points length_fibril
    
Here,

- N_files denotes the number of files over which the average is calculated

- common_name denotes the file location up to the name common to all files

- outputfilename denotes the name of the averaged file

- N_points denotes the number of time points to be considered. This is required because the averaging procedure differs from the previous one

- length_fibril denotes the length of the microfibril within the simulations leading to the raw data

As an example, this would be the command used to calculate the average over 100 DP_distribution data files:

    python3 calc_mean_heatmap.py 100 DP_distrib/DP_distrib_ mean_DP_distrib.txt 2000 100
    
The averaged files are structured in the same way as the raw data (see above).




# Code structure


The simulation code is split into multiple files: 

- code_4.cpp
- functions.cpp
- functions.hpp
- structs.hpp

Each of them is summarized briefly below. See the files themselves for further documentation


## structs.hpp

This file stores two classes and two structs:

- class bList is used to generate polymer objects for cellulose, hemicellulose and lignin.

- class DPList is used to record the distribution of glucose within cellulose polymers

- struct TList is used to generate reaction table objects to be used for the enzymatic reactions within the Gillespie algorithm

- struct params is used to generate an object containing those system parameters which are common to each individual simulation run.
  They are then easier to pass between functions.


## functions.cpp and functions.hpp

These files store all functions within the code except for the functions main() and run(). The task of each function is summarized within the code files.

## code_4.cpp

This file contains the main() function, as well as a function run(). Within the main function, the parameters specified in the "Params" folder are loaded,
and the data generated by the individual simulation runs is transferred to output files.
An individual simulation run is carried out within the function run().
These are further documented within their source code






