# Simulated saccharification of a lignocellulose microfibril

This README provides instructions on how to run the simulation code, and a brief documentation of its structure.

## System requirements

The "code_4" binary was compiled using GCC version 10.2.1 20210110 (Debian 10.2.1-6) on 64-bit linux, so it can be run as is on such systems.

For windows and mac you need to recompile the source code (see below). Please be aware that in this case, file paths in the source code may need to be updated in order to comply with the respective system (e.g. "\" vs "/" for windows vs linux).

## TLDR

You can run the code by executing code_4 via the terminal command

        ./code_4

## Compiling, setup and first run

You can skip the compiling step, if you are using 64-bit linux and already have a file named "code_4"

### Compiling on Linux

- Run compile.sh via the terminal command. 

        ./compile.sh

### Compiling on windows or mac

- compile.sh contains a single command that uses GCC to build the executable from the source files within the "Code" directory. Depending on which compiler you want to use, this command needs to be adjusted accordingly.

- As said above, file paths in the source code may need to be updated for operating systems other than linux.

### Further setup 

- Copy the code executable "code_4" to the directory from which you want to run it.

- Next, copy the folders "Params" and "Output" into the same directory as the executable

- Make sure that the "Params" folder contains at least the following files:

    - "simulation_parameters.txt"
    - "kinetic_parameters.txt"
    - "initial_configuration_parameters.txt"
    - "structural_parameters_1.txt"
    - "structural_parameters_2.txt"
    - "structural_parameters_3.txt"
    - "structural_parameters_4.txt"
    - "structural_parameters_5.txt"
    - "hydro_parameters_2_1_layer_1.txt"
    - "hydro_parameters_2_1_layer_2.txt"
    - "hydro_parameters_2_1_layer_3.txt"
    - "hydro_parameters_2_1_layer_4.txt"
    - "hydro_parameters_2_2_layer_1.txt"
    - "hydro_parameters_2_2_layer_2.txt"
    - "hydro_parameters_2_2_layer_3.txt"
    - "hydro_parameters_2_2_layer_4.txt"
    - "hydro_parameters_2_3_layer_1.txt"
    - "hydro_parameters_2_3_layer_2.txt"
    - "hydro_parameters_2_3_layer_3.txt"
    - "hydro_parameters_2_3_layer_4.txt"
    - "hydro_parameters_2_4_layer_1.txt"
    - "hydro_parameters_2_4_layer_2.txt"
    - "hydro_parameters_2_4_layer_3.txt"
    - "hydro_parameters_2_4_layer_4.txt"
    - "hydro_parameters_2_5_layer_1.txt"
    - "hydro_parameters_2_5_layer_2.txt"
    - "hydro_parameters_2_5_layer_3.txt"
    - "hydro_parameters_2_5_layer_4.txt"

- Make sure that the "Output" folder contains the following subfolders:

    - "3D"
    - "DP_distrib"
    - "enzyme_activity"
    - "enzyme_fraction"
    - "enzyme_concentration"
    - "Nbr_reactions"
    - "saccharification"
    - "expe_data"


- You can now run the code via the terminal command

    ./code_4


## Customizing the parameters

    
- You can adjust the parameters within "simulation_parameters.txt", "kinetic_parameters.txt" and "initial_configuration_parameters.txt" in
 the "Params" folder to suit your requirements. The parameters are ordered as follows:
 
### simulation_parameters:

1.) and 2.) sum up to yield the maximum number of steps (there are two parameters for this in case of possible future added effects).
    Example: the first parameter is set to 10000, the second is set to 5000. Then, no more than 15000 Gillespie algorithm steps can be performed [These parameters are used in the code_4.cpp to set the maximum sumber of Gillespie steps (current value: 200000 + 200000 = 400000) but are not varied and can be fixed in the interface.]

3.) Maximal number of snapshots which may be taken for the purpose of creating animation data. [This is used in the code to limit the maxmimum number of snapshots. This can be fixed in the interface too.]

4.) Maximal value for the real time within the simulation, as calculated by the Gillespie algorithm. [This is the maximum real time for the simulation to run, calculated by Gillespie algorithm. If this is set to a high value, the simulation runs until all digestible substrate has been digested. Else it can be set to lower value by the user to match the duration of the experimental saccharification time-course data available.]

5.) Number of steps between system snapshots, if "-vid" command is used (see below). If this is set to 1, a snapshot is taken at every step of the Gillespie algorithm.

6.) Maximum number of hemicellulose/lignin layers around the cellulose core. Currently cannot be set above 4. (This can be fixed to 4 in the interface. The user doesn't need to change it in the majority of the scenario.) 

7.) mode_hemi: Determines, on which sides of the microfibril hemicellulose can anchor. Currently to be kept at 2. [User doesn't need to vary it. Can be fixed at 2 in the interface.]

8.) mode_lign: Determines, on which sides of the microfibril lignin can anchor. Currently to be kept at 2. [User doesn't need to vary it. Can be fixed at 2 in the interface.]

9.) mode_inhib: Determines, whether inhibition is active (1) or inactive (-1). [For user, value 1 is reccommended as end-product inhibition is real.] 

10.) mode_lignin_glue: Determines, whether the lignin gluing effect is active (1) or inactive (-1). [For user, value 1 is recommended, as Lignin gluing is a known phenomena.]

11.) mode_enzyme_size : Determines, whether the enzyme size is included for all enzymes (1) or only for CBH (-1). [This can be fixed to 1 in the interface for now.]

12.) Nruns: Determines the number of simulation runs. [10 runs are good enough for getting a smooth average saccharification curve.]

13.) mu_lignin_covering: Determines the average fraction of lignin polymers which acts as a structural barrier. Set this to a value between 0 and 1. [Fix at 0.5]

14.) sigma_lignin_covering : Determines the standard deviation around mu_lignin_covering. Set this to a value between 0 and 1. [Fix at 0.2]


### kinetic_parameters:

1.) EG Kcat (diffusive action)      [units: 1/s. Range: 0.001-1000]

2.) EG Km (diffusive action)        [units: mM. Range: 0.001-1000]

3.) CBH processive digestion rate   [units: 1/h. Range: 10-10000]

4.) BGL Kcat (diffusive action)     [units: 1/s. Range: 0.001-1000]

5.) BGL Km (diffusive action)       [units: mM. Range: 0.001-1000]

6.) XYL Kcat (diffusive action)     [units: 1/s. Range: 0.001-1000]

7.) XYL Km (diffusive action)       [units: mM. Range: 0.001-1000]

8.) Lignin_adhesion_rate            [Just a number. Range: 150-350]

9.) CBH Kcat (attachment rate)      [units: 1/s. Range: 0.001-1000]

10.) CBH Km (attachment rate)       [units: mM. Range: 0.001-1000] 

11.) inhibition binding affinity of cellobiose to EG [dimensionless number. Range: 0-1]

12.) inhibition inding affinity of cellobiose to CBH [dimensionless number. Range: 0-1]

13.) inhibition binding affinity of glucose to EG    [dimensionless number. Range: 0-1]

14.) inhibition inding affinity of glucose to CBH    [dimensionless number. Range: 0-1]

15.) inhibition inding affinity of glucose to BGL    [dimensionless number. Range: 0-1]

16. digestibility ratio of "crystalline" and "amorphous" cellulose [dimensionless number. Range: 0.00001-0.1]

17.) digestibility ratio of "crystalline" and "amorphous" hemicellulose [dimensionless number. Range: 0.00001-0.1]

18.) Enzyme size: radius   [Units: nm. Can be fixed at 2.5 nm]



### initial_configuration_parameters:

1.) mode_code: Determines the shape of the microfibril. To use the 36 chain microfibril used by Ding et al., set this to 5 (Spruce wood, 24 polymers, mode_code = 1 or 2; Mung beans, 18 polymers, mode_code = 3 or 4; Maize, 36 polymers, mode_code = 5). [To take value according to plant sample species]

2.) pct_EG [fraction of EG in enzyme coctail. Range: 0-1.]

3.) pct_CBH [fraction of CBH in enzyme coctail. Range: 0-1.]

4.) pct_BGL [fraction of BGL in enzyme coctail. Range: 0-1.] 

5.) pct_XYL [fraction of XYL in enzyme coctail. Range: 0-1.]

6.) total_enz_molecules [Total number of enzyme molecules in the surrounding volume of a subunit. Fix at 50. But can be varied in the range: 50-150, resulting in faster saccharification.]

7.) length_fibril [Fix at 200. Can be increased to higher values to simulate longer microfibril at the expense of increased simulation time.]

8.) boolean_Xyl_or_MLG [Fix it to 1. If set to 1, there is Xylans in Hemicellulose (xylose released). If set to 0, there is Mixed-Linkage Glucans (glucose released)]

9.) pct_xyl: Percentage of xylose in hemicellulose; currently restricted to a value of 1. [Should be fixed at 1 for now. But can be useful later when other hemicellulose sugars are implemented simultaneously.]

10.) pct_cellu: Percentage of cellulose in microfibril composition [Comes from composition data. User should input. Range 0.01-1.]

11.) pct_hemi: Percentage of hemicellulose within the microfibril. [Comes from composition data. User should input. Range 0-1.] 

12.) pct_lign: Percentage of lignin within the microfibril. [Comes from composition data. User should input. Range 0-1.] 

13.) pct_acetyl_hemi: Percentage of acetylated hemicellulose within the hemicellulose fraction; currently restricted to a value of 0. [Fix it at 0 for now. Allows an extra option if acetylation of hemicellulose is considered and modified later.]

14.) pct_crystalline_cellu: Percentage of "crystalline" cellulose. [Fitting parameter. Can also come from user data, if available. Range 0-1.]

15.) pct_crystalline_hemi: Percentage of "crystalline" hemicellulose. [Fitting parameter. Can also come from user data, if available. Range 0-1.]

16.) Mean_defect_size: Mean size of defect/damaged section of Crystalline cellulose. [Implemented for all mode_codes. Takes value 0-0.5. It signfies the fraction of amorphous cellulose (out of len_poly) contained in the defects.]

17.) Nbr_of_defects: Number of defects/damaged sections of Crystalline cellulose. [Implemented for all mode_codes. Takes value 0-1. It signfies the fraction of the number of outer cellulose polymers which have defects.]

18.) r_monomer: Radius around which each bond location is checked for neighbors. To be kept at 0.6 [Fixed for now. No need to vary]




## Running the customized code

- As above, you may run the code via the terminal command

    ./code_4


- There are several additional command line arguments which may be added. For command line arguments X1, X2, ... Xn, the code would be run as follows:

    ./code_4 X1 X2 ... Xn
    
    **IMPORTANT: the very first command line argument determines, from which file the initial configuration parameters are read. Specifically this means that for the command**

    ./code_4 example
    
    **the file "initial_parameters_example.txt" is loaded. This also leads to a change in the name of the output files:**
    **the files in the saccharification folder for example will be named "saccharification_example_1.txt", "saccharification_example_2.txt", etc. (the number represents the number of the individual simulation run).**
    **If there are no command line arguments, the default "initial_configuration_parameters.txt" file is loaded instead, and the output files will simply be called "saccharification_1.txt", "saccharification_2.txt", etc. Alternatively, if you want to use command line arguments, but do not want a file suffix, use the command line argument "no_suffix" first**

    The currently available command line arguments are as follows:

    "-no_suffix": If used, this should be the first command line argument. It will suppress the use of a suffix on the output files, and the file "initial_configuration_parameters.txt" will be used 

    "-verbose": This will cause more detailed information of what is happening during the individual simulation runs to be written to the terminal.

    "-heatmap": this enables the Output of files documenting the distribution of glucose within cellulose polymers of all lengths within the "DP_distrib" folder.
    As this is somewhat computationally expensive, it was chosen as an optional component.

    "-vid": This causes the first simulation run to document the 3D positioning of all cellulose, hemicellulose and lignin bonds at step intervals, which are determined by a parameter in "simulation_parameters.txt" (see above).
    They are saved in the folder "3D" and can be used to create a visualization of the microfibril digestion.

    "-fixed_seed": This fixes the seed of the pseudo-random number generator used within the simulation. The seed which will be set is currently specified within the source code, in the file "code_4.cpp"
    
    "-parallel": This specifies, how many cores should be used by the program, and needs to be followed by an integer greater than or equal to 1. The default number of cores is 1
    
    "-print_polys": This enables the output of remaining bonds after each simulation run (if there are any)

    "-print_ends_blocked": This enables the output of the average number of ends blocked by CBH enzymes over the simulation"

    "-print_CBH_positions": This enables tracking the positions of each attached CBH enzyme over the course of the simulation

    "-toy_model": This disables structural features of the model: the polymers will be located far from each other and therefore be treated as being fully in solution

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

The files named 'visualisation_total_*.txt' contain 5 columns:

- First column: x-coordinates
- second column: y-coordinates
- third column: z-coordinates
- fourth column: polymer type (0: lignin, 1: hemicellulose, 2: cellulose)
- fifth column: timestamp (real time from Gillespie algorithm)


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

There is a C++ program within the "Output" folder, which is called "calc_mean_heatmap.cpp", as well as its compiled version "calc_mean_heatmap". Use it via the terminal command

    ./calc_mean_heatmap N_files common_name outputfilename N_points length_fibril half_point_divider
    
Here,

- N_files denotes the number of files over which the average is calculated

- common_name denotes the file location up to the name common to all files

- outputfilename denotes the name of the averaged file

- N_points denotes the number of time points to be considered. This is required because the averaging procedure differs from the previous one

- length_fibril denotes the length of the microfibril within the simulations leading to the raw data

- half_point_divider sets the resolution of the averaged file over the simulated time. For the maximum simulated time t_end, 0.5*N_points will be used below half_point_divider*t_end, and the other half will be used above. This means that only if half_point_divider is set to 0.5, the time within the averaged file will be an even grid


As an example, this would be the command used to calculate the average over 100 DP_distribution data files:

    ./calc_mean_heatmap 100 DP_distrib/DP_distrib_ mean_DP_distrib.txt 2000 100 0.5
    
The averaged files are structured in the same way as the raw data (see above).


### Creating movie of the microfibril digestion

Creating a movie from the simulation is now done simply by running the shell script **Moviemaker.sh** (to use simply make script executable and run ```$chmod +x Moviemaker.sh && ./Moviemaker.sh```). 

The `Moviemaker.sh` script has a choice of 2 modes (Mode 1, Mode 2). The user is prompted to enter the choice (1,2) on running the script.

Mode 1: Make animation for a theoretical simulation NOT comparing to any experimental data. Simulation runs till all substrate is digested.

Mode 1: Make animation for a simulation WITH comparison to experimental glucose data. Experimental data to be compared, must be saved in the path `Output/expe_data/expe_saccharification_movie_glc.txt`, file must be named `expe_saccharification_movie_glc.txt`. Simulation runs till the maximim time in the data file.

It first runs a single Gillespie simulation run using the parameters from `Params/initial_configuration_parameters_movie.txt` with the command line argument "-vid".  It then uses python matplotlib via the pythons scripts in the path  `scripts_for_movie/frame_*.py`, to generate the frames with rotating camera angle. It finally uses ffmpeg (must be installed on the system) to create a mp4 video from the generated frames. Then it cleans up all the data generated from the simulation to prevent unwanted file storage space shortage.

Requirements: ffmpeg, python libraries (matplotlib, pyfiglet, math, numpy, glob, os, re, tqdm, PIL)


# Code structure


The simulation code is split into multiple files: 

- code_4.cpp
- functions.cpp
- functions.hpp
- params.hpp
- bList.hpp
- TList.hpp
- DPList.hpp
- CBH_enzyme.hpp
- neighborList.hpp
- tuple_hash.hpp

Each of them is summarized briefly below. See the files themselves for further documentation


## code_4.cpp

This file contains the main() function, as well as a function run(). Within the main function, the parameters specified in the "Params" folder are loaded,
and the data generated by the individual simulation runs is transferred to output files.
An individual simulation run is carried out within the function run().
These are further documented within their source code

## functions.cpp and functions.hpp

These files store all functions within the code except for the functions main() and run(). The task of each function is briefly summarized within the code files.

## params.hpp

The struct params is used to generate an object containing those system parameters which are common to each individual simulation run. They are then easier to pass between functions.

## bList.hpp

The class bList is used to generate polymer objects for cellulose, hemicellulose and lignin.

## TList.hpp

The struct TList is used to generate reaction table objects to be used for the enzymatic reactions within the Gillespie algorithm

## DPList.hpp

The class DPList is used to record the distribution of glucose within cellulose polymers

## CBH_enzyme.hpp

The class CBH_enzyme defines the properties of CBH enzymes within the simulation

## neighborList.hpp

The class neighborList is used to record the neighbors of each polymer bond

## tuple_hash

This defines a hash function for tuples, which are used as keys for unordered maps of neighborList objects in the code.









