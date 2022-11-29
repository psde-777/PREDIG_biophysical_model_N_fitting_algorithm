#ifndef PARAMS_HXX_
#define PARAMS_HXX_
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <string.h>
#include <list>
#include <algorithm>
#include <functional>
#include <math.h>
#include <sys/time.h>
#include <tuple>
#include <unordered_map>
#include "neighborList.hpp"
#include "CBH_enzyme.hpp"


struct params {//Struct containing some parameters which are read-in from file



    //====================== Define directory and names of input and output files ==========================
    std::string paramFileDir = "Params/";//Directory in which the input files are stored.
    std::string simuParams = "simulation_parameters";//This file defines the simulation time and parameters related to the main structure of the fibril
    std::string initConfigParams = "initial_configuration_parameters.txt";//This file defines the enzyme concentrations and the structure of the microfibril:
                                                                // its length, the percentage of hemicellulose, the percentage of lignin and the percentage of acetylated hemicellulose
    std::string structParams = "structural_parameters_";//This file defines the configuration of the fibril irrespective of the availability of the bonds (conditioned by the positioning and amounts of lignin and hemicellulose)
    std::string hydroParams = "hydro_parameters";//This file defines the positions on a section where hemicellulose and lignin can anchor respectively.
    std::string kinParams = "kinetic_parameters";//This file defines the kinetic parameters of the enzymes involved in saccharification

    std::string outputFileDir = "Output/";//Directory in which the output files are stored.

    std::string glcFileName = "amount_Glc";//This file stores the amount of glucose released from the microfibril during each snapshot.
    //====================== These parameters are assigned via input files ============================


    int Nruns;//The number of times the whole simulation is supposed to be repeated
    int N_layers;//The number of hemi/lignin layers around the cellulose core
    int mode_code;//This selects the configuration of the substrate: 24 Glc rectangle (1), 24 Glc diamond (2), 18 Glc rectangle (3), 18 Gls diamond (4), or test, which is 42
    int mode_hemi;//This value is 1 if hemicellulose prefers hydrophilic surfaces -1 if it prefers hydrophobic surfaces
    int mode_lign;//This value is 1 if lignin prefers hydrophilic surfaces -1 if it prefers hydrophobic surfaces
    int mode_inhib;// This value is 1 if inhibition of enzymes by glucose is included, and -1 if not
    int DP_print_Freq;//Defines, how often a snapshot of the DP distribution in the system is made
    int mode_lignin_glue;// Determines, whether the lignin gluing effect is active (1) or inactive (-1)
    //int mode_hemi_structure;//Determines, whether the initial hemicellulose distribution can have holes (1) or not (-1)
    int mode_enzyme_size;//Determines, whether the enzyme size is included as an effect
    int min_x_outer;//Boundaries of fibril
    int max_x_outer;
    int min_y_outer;
    int max_y_outer;
    int pict_3D_Freq; //Defines, how often a 3D structure picture is taken during gillespie loop
    int N_free_ends; //Number of cellulose polymer ends in the system, which are attackable by CBH

    int nbr_hemi_monomer; //The number of hemicellulose monomers
    int nbr_monolignol; //The number of lignin monomers
    int nbr_lignin_blocked;//Number of lignin monomers blocked by enzymes adhered to it


    double mid_x;//Coordinates of the middle of the microfibril
    double mid_y;
    double mid_z;





    unsigned int length_fibril;//Length of the whole fibril at time 0, counted in BONDS from 1

    double init_EG;//Initial number of EG = endoglucanases
    double init_CBH;//Initial number of CBH = cellobiohydrolases = exoglucanases
    double N_free_CBH;//Number of CBH enzymes which are not currently attached to a polymer
    double ends_blocked_per_CBH;//Average number of polymer ends blocked by an attached CBH. Currently,this is set equal to the initial number of cellulose polymers, as one CBH enzyme covers the entire area at the end of the fibril
    double init_BGL;//Initial number of BGL = beta-glucosidase
    double init_XYL;//Initial number of XYL = xylanase
    double init_EG_save;//Initial number of EG = endoglucanases
    double init_CBH_save;//Initial number of CBH = cellobiohydrolases = exoglucanases
    double init_BGL_save;//Initial number of BGL = beta-glucosidase
    double init_XYL_save;//Initial number of XYL = xylanase

    unsigned long int T;//Number of iterations for the stationnary phase of the simulation
    unsigned long int Transient;//Number of iterations for transient
    unsigned long int Nbr_picts;//Number of pictures of the fibril to be taken
    double real_time;//The time as calculated via the gillespie algorithm
    double max_time;//The maximum gillespie time (simulation stops, if this is exceeded)


    //CBH speed and association time. Data from Brady et al. 2015, nature communications
    double average_CBH_step_number = 50;//average number of steps taken by an attached CBH
    double average_CBH_speed = 0.25;//nm/s
    double std_dev_CBH_speed = 0.16;//nm/s
    double average_CBH_association_time = 90./3600.;//average time IN HOURS a processive CBH stays attached. In seconds it is 90
    double CBH_reaction_rate = average_CBH_step_number / average_CBH_association_time;
    //v_CBH = (0.25 +- 0.16) nm/s
    double pct_xyl;//percentage of xylose in hemicellulose polymers
    double pct_lign;//Percentage of lignin in microfibril
    double pct_lign_subset;//Percentage of lignin in outer layer
    double pct_hemi;//Percentage of hemicellulose in microfibril
    double pct_hemi_subset;//Percentage of hemicellulose in outer layer
    double pct_acetyl_hemi;//Percentage of acetylated hemicellulose bonds in hemicellulose polymers
    double pct_glc;//Percentage of cellulose. Calculated from given hemi/lignin values
    double pct_crystalline_cellu;//Percentage of crystalline cellulose
    double pct_crystalline_hemi;//Percentage of crystalline hemicellulose
    double dfct_size; // Percentage of total amorphous cellulose surrounded by crystalline
    int N_amor_core; // Number of amorphous cores hidden in the crystalline


    double k1;//Rate of reaction for EG
    double k2;//Rate of reaction for CBH
    double k3;//Rate of reaction for by BGL
    double k4;//Rate of reaction for XYL
    double k5;//Lignin adhesion rate
    double k6;//CBH attachment rate
    

    double inhib_cellobiose_EG;//Inhibition weight for EG :cellobiose
    double inhib_cellobiose_CBH;//Inhibition weight for CBH :cellobiose
    double inhib_glucose_EG;//Inhibition weight for EG :glucose
    double inhib_glucose_CBH;//Inhibition weith for CBH :glucose
    double inhib_glucose_BGL;//Inhibition weight for BGL
//    double inhib_XYL;//Inhibition weight for XYL
    double crystal_modifier_cellu;//Modification factor for proṕensity of crystalline cellulose
    double crystal_modifier_hemi;//Modification factor for proṕensity of crystalline hemicellulose



    double enzyme_radius;//Radius of all enzymes in terms of bonds
    double V_enzyme;//Volume of a single enzyme
    double r_monomer;//Radius of a glucose unit
    double mu_lignin_covering;//Average percentage of covering lignin units per polymer
    double sigma_lignin_covering;//Standard deviation of gaussian distribution of percentage of covering lignin units per polymer
    double lignols_blocked_per_enzyme;//Defines how many monolignols it takes to bind an enzyme


    std::string input_file;
    std::string output_file;
//    std::string output_file;

    bool verbose = false;//If this is set to 1, an increased number of progress messages will be printed
    bool enzyme_timer = false;//If this is set to 1, the time each digestion take is tracked and printed at the end of the run
    bool heatmap_bool = false;//If this is set to 1, a (potentially large) file will be generated, which contains the degree-of-polymerization(DP) distribution during each simulation step
    bool print_polys = false;//If this is set to true, the positions of all remaining bonds (if there are any) will be printed to file after a finished simulation
    bool print_ends_blocked = false;//If this is set to true, the average number of blocked polymer ends by CBH enzymes will be tracked and printed to file
    bool print_CBH_positions = false;//If this is set to true, the positions of all attached CBH enzymes will be tracked
    bool toy_model = false;//If this is set to true, every bond will be accessible from the beginning. This resembles the situation of all polymers freely floating in space
    bool print_time = false;//If this is set to true, the computation time of each run will be printed to a file
    bool xyl_or_mlg = true; //If this is set to true then Xylose in outer shell, if false then MLGs in outer shell

    std::vector<double> propensities; //The vector containing the propensities for each reaction type. ALL REACTIONS IN THE REACTION TABLES CONTAIN POINTERS TO THESE VALUES
    std::vector<double> crystal_propensities; //The vector containing the crystalline propensities for each reaction type. ALL REACTIONS IN THE REACTION TABLES CONTAIN POINTERS TO THESE VALUES

    //Vectors for tracking glucose/cellobiose/xylose release
    std::vector<double> amount_glc_mean;
    std::vector<double> amount_cellobiose_mean;
    std::vector<double> amount_xylose_mean;
    std::vector<double> time_mean;
    std::vector<double> Nbr_reactions;

    //Vectors for tracking enzyme concentrations:
    std::vector<double> EG_conc_mean;
    std::vector<double> CBH_conc_mean;
    std::vector<double> BGL_conc_mean;
    std::vector<double> XYL_conc_mean;
    std::vector<double> enzymes_glued;//enzymes glued by lignin

    int N_enzymes_glued;

    //Vectors for counting the times, each enzyme is used as a fraction of overall step number
    std::vector<double> EG_activity;
    std::vector<double> CBH_activity;
    std::vector<double> BGL_activity;
    std::vector<double> XYL_activity;
    std::vector<double> lign_activity;

    //Vectors for tracking the fraction of reactions of different enzyme types among the reaction tables  
    std::vector<double> EG_fraction;
    std::vector<double> CBH_fraction;
    std::vector<double> BGL_fraction;
    std::vector<double> XYL_fraction;   
    std::vector<double> lign_fraction; 

    std::vector<double> timestamp;//Times at which the fractions of reactions are stored

    std::vector<CBH_enzyme> CBH_enzymes;//The vector containing all CBH_enzyme objects
    std::vector<double> average_ends_occupied_by_CBH; //The vector containing the average number of ends occupied by a CBH enzyme
    std::vector<double> number_of_attached_CBH; //The vector containing the number of attached CBH enzymes at each step
    std::unordered_map<int,std::tuple<int,int>> free_poly_ends;//Stores the polymer and position of polymer ends which are not yet blocked by CBH
//    std::vector<double> CBH_clocks;//Tracks the attachment times of the enzymes
//    std::vector<bool> clock_used;



    void initialize(){
        N_layers = 1;
        mode_code = 0;
        mode_hemi = 0;
        mode_lign = 0;
        mode_inhib = 0;
        length_fibril = 0;
        init_EG = 0;
        init_CBH = 0;
        init_BGL = 0;
        pict_3D_Freq = 1;
        T = 0;
        Transient = 0;
        Nbr_picts = 0;
        r_monomer = 1;
        N_free_ends = 0;
        real_time = 0;
        nbr_monolignol = 0;
        nbr_lignin_blocked = 0;
        nbr_hemi_monomer = 0;
        pct_lign = 0;
        pct_hemi = 0;
        pct_acetyl_hemi = 0;
        DP_print_Freq = 1;
        enzyme_radius = 0;
        V_enzyme = 0;
        input_file = "default";
        output_file = "default";
        mid_x = 0.;
        mid_y = 0.;
        mid_z = 0;
//        output_file = "test";
    }
};




#endif
