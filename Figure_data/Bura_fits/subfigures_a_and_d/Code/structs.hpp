#ifndef STRUCTS_HXX_
#define STRUCTS_HXX_


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



/****************************Definition of Structures*******************/


class bList//Structure of polymers: cellulose, hemicellulose and lignin
{
public:
    int index;//Index of the polymer
    double x;//Abscissa
    double y;//Ordinate
    int starting_z;//z coordinate of bond 0
    std::vector<int> z;//Applicate of the bonds
    std::vector<int> status;//Indicates -2 if the bond is acetylated or lignin, -1 if the bond is not yet accessible, 1 if the bond is accessible, 0 if the polymer is fully digested
    std::vector<int> bond_type;//Indicates 1 for glc-glc, 2 for xyl-glc, 3 for glc-xyl, 4 for xyl-xyl, and -1 for lign-lign. For now this is only important for hemicellulose
    std::vector<bool> crystalline;//If set to true, the corresponding bond is "crystalline", otherwise it is "amorphous"
    std::vector<int> N_blocked_positions;//Indicates the number of occupied positions within a sphere of radius enzyme_size around the bond 
    int len_poly;//Length of the polyssacharide it belongs to

//    bool exists;//If this is set to false, the polymer will "not exist", i.e. not interact or be digestable
    int reactionTable; //Reaction table to which this polymer belongs. This should be the same as index


    void set_z(int z_set){
        starting_z = z_set;
    }

    bList(){
        this->starting_z=0;
        this->index=0;
        this->x=0;
        this->y=0;
        this->z.push_back(0);
        this->z.clear();
        this->status.push_back(0);
        this->status.clear();
        this->bond_type.push_back(0);
        this->bond_type.clear();
        this->crystalline.push_back(false);
        this->crystalline.clear();
        this->N_blocked_positions.push_back(0);
        this->N_blocked_positions.clear();
        this->len_poly=0;
//        std::cout << "z.size() = " << z.size() << std::endl;
    }


    ~bList(){

    }

};

class DPList {//Structure of the DP distribution at a time t

public:
    double real_time;//real_time at which this is taken
    int t0;//step number
    std::vector<double> DP;//this was an array

    DPList(){
        this->real_time = 0;
        this->t0 = 0;
        //DP = new int[1];
        DP.push_back(0);
    }

    DPList(double real_time, int t0, int length_fibril){
        this->real_time = real_time;
        this->t0 = t0;
//        DP = new int[length_fibril];
        for(int i=0; i<length_fibril;i++)
            DP.push_back(0);
        //    DP[i] = 0;

    }

    DPList(int length_fibril){
        this->real_time = 0;
        this->t0 = 0;
//        DP = new int[length_fibril];
        for(int i = 0;i<length_fibril;i++)
            DP.push_back(0);
  //          DP[i] = 0;
    }


    ~DPList(){
     //   delete[] DP;
    }

};


struct TList//Structure of the reaction table
{
//    int originalTable; //Original table from which this polymer originates 
    int nbr_element;//Number of elements of the structure
    double prop_sum;//Sum of the propensities of all reactions inside this structure
    int index_poly;//Index of the polymer targetted
    std::vector<int> num_bond;//Index of the bond inside the "*z" list of the polymer
    std::vector<int> material;//Indicates 1 if cellulose is the substrate 2 if hemicellulose is the substrate
    std::vector<double> indic_action;//Action label: digestion by EG =1, CBH =2, BGL =3;
    std::vector<double> liste_prop;//Propensity associated to this action
    std::vector<double> prop_uninhib; //Propensity without inhibition



    void calcTableProp(){//Calculates the starting propensity sum of the table. Only used during initialization
        prop_sum = 0;
        for(int i=0; i<nbr_element;i++){
            prop_sum += liste_prop[i];
        }
    }
    
    void addProp(double prop){//Increase prop_sum by prop
        prop_sum += prop;
    }

    void subProp(double prop){//Reduce prop_sum by prop
        double prop_sum_before = prop_sum;
        double epsilon = 1e-10;//values below this are considered to be 0
//        std::cout << "Before subtraction: prop_sum = " << prop_sum << std::endl;        
        if(prop_sum > 0)
            prop_sum -= prop;
        if(prop_sum < -epsilon)
            std::cout << "NEGATIVE PROPENSITIES OCCURING!!! prop = " << prop << ", prop_sum before: " << prop_sum_before << ", prop_sum now: " << prop_sum << std::endl;

    }

};



struct params {//Struct containing some parameters which are read-in from file



    //====================== Define directory and names of input and output files ==========================
    std::string paramFileDir = "Params/";//Directory in which the input files are stored.
    std::string simuParams = "simulation_parameters";//This file defines the simulation time and parameters related to the main structure of the fibril
    std::string initConfigParams = "initial_configuration_parameters.txt";//This file defines the enzyme concentrations and the structure of the microfibril:
                                                                // its length, the percentage of hemicellulose, the percentage of lignin and the percentage of acetylated hemicellulose
    std::string structParams = "structural_parameters_";//This file defines the configuration of the fibril irrespective of the availability of the bonds (conditioned by the positioning and amounts of lignin and hemicellulose)
    std::string hydroParams_inner = "hydro_parameters_inner_";//This file defines the positions on a section where hemicellulose and lignin can anchor respectively (inner layer)
    std::string hydroParams_outer = "hydro_parameters_outer_";//This file defines the positions on a section where hemicellulose and lignin can anchor respectively (outer layer)
    std::string kinParams = "kinetic_parameters";//This file defines the kinetic parameters of the enzymes involved in saccharification

    std::string outputFileDir = "Output/";//Directory in which the output files are stored.

    std::string glcFileName = "amount_Glc";//This file stores the amount of glucose released from the microfibril during each snapshot.
    //====================== These parameters are assigned via input files ============================

    int Nruns;
    int free_bonds_req;//Number of free neighbors required for digestability(depends on enzyme_radius)
    int Run_number;
    int mode_code;//This selects the configuration of the substrate: 24 Glc rectangle (1), 24 Glc diamond (2), 18 Glc rectangle (3), 18 Gls diamond (4), or test, which is 42
    int mode_hemi;//This value is 1 if hemicellulose prefers hydrophilic surfaces -1 if it prefers hydrophobic surfaces
    int mode_lign;//This value is 1 if lignin prefers hydrophilic surfaces -1 if it prefers hydrophobic surfaces
    int mode_inhib;// This value is 1 if inhibition of enzymes by glucose is included, and -1 if not
    int DP_print_Freq;//Defines, how often a snapshot of the DP distribution in the system is made
    int mode_lignin_glue;// Determines, whether the lignin gluing effect is active (1) or inactive (-1)
    int mode_hemi_structure;//Determines, whether the initial hemicellulose distribution can have holes (1) or not (-1)
    int mode_enzyme_size;//Determines, whether the enzyme size is included as an effect
    int min_x_outer;//Boundaries of fibril
    int max_x_outer;
    int min_y_outer;
    int max_y_outer;
    int pict_3D_Freq; //Defines, how often a 3D structure picture is taken during gillespie loop

    int nbr_hemi_monomer;
    int nbr_monolignol;
    int nbr_lignin_blocked;//Number of lignin monomers blocked by enzymes adhered to it

    int min_x_cellu;
    int max_x_cellu;
    int min_y_cellu;
    int max_y_cellu;

    int min_x_hemi;
    int max_x_hemi;
    int min_y_hemi; 
    int max_y_hemi;



    int min_x_lign;
    int max_x_lign;
    int min_y_lign;
    int max_y_lign;


    int min_x_overall;
    int max_x_overall;
    int min_y_overall;
    int max_y_overall;

    float mid_x;//Coordinates of the middle of the microfibril
    float mid_y;
    float mid_z;





    unsigned int length_fibril;//Length of the whole fibril at time 0, counted in BONDS from 1
    double init_EG;//Initial number of EG = endoglucanases
    double init_CBH;//Initial number of CBH = cellobiohydrolases = exoglucanases
    double init_BGL;//Initial number of BGL = beta-glucosidase
    double init_XYL;//Initial number of XYL = xylanase
    double init_EG_save;//Initial number of EG = endoglucanases
    double init_CBH_save;//Initial number of CBH = cellobiohydrolases = exoglucanases
    double init_BGL_save;//Initial number of BGL = beta-glucosidase
    double init_XYL_save;//Initial number of XYL = xylanase

    unsigned long int T;//Number of iterations for the stationnary phase of the simulation
    unsigned long int Transient;//Number of iterations for transient
    unsigned long int Nbr_picts;//Number of pictures of the fibril to be taken
    double max_time; 
    float pct_xyl;//percentage of xylose in hemicellulose polymers
    float pct_lign;//Percentage of lignin in microfibril
    float pct_lign_subset;//Percentage of lignin in outer layer
    float pct_hemi;//Percentage of hemicellulose in microfibril
    float pct_hemi_subset;//Percentage of hemicellulose in outer layer
    float pct_acetyl_hemi;//Percentage of acetylated hemicellulose bonds in hemicellulose polymers
    float pct_glc;//Percentage of cellulose. Calculated from given hemi/lignin values
    float pct_crystalline_cellu;//Percentage of crystalline cellulose
    float pct_crystalline_hemi;//Percentage of crystalline hemicellulose
/*
    double save_k1;//Safety values to return to after each run
    double save_k2;
    double save_k3;
    double save_k4;
    double save_k5;*/

    double k1;//Rate of reaction for EG
    double k2;//Rate of reaction for CBH
    double k3;//Rate of reaction for by BGL
    double k4;//Rate of reaction for XYL
    double k5;//Lignin adhesion rate
    

    double inhib_cellobiose_EG;//Inhibition weight for EG
    double inhib_cellobiose_CBH;//Inhibition weight for CBH
    double inhib_glucose_BGL;//Inhibition weight for BGL
    double crystal_modifier_cellu;//Modification factor for proṕensity of crystalline cellulose
    double crystal_modifier_hemi;//Modification factor for proṕensity of crystalline hemicellulose

    double enzyme_radius;//Radius of all enzymes in terms of bonds
    double V_enzyme;//Volume of a single enzyme
    double r_monomer;//Radius of a glucose unit

    char str[50];
    char str2[50];
    char str5[50];
    char str3[50];
    char str4[50];
    char str6[50];
    char str7[50];
    char str8[50];
    char str9[50];
    char input_file[100];
    char output_file[100];

    bool enzyme_timer;
    bool verbose;

    //Vectors for tracking glucose/cellobiose/xylose release
    std::vector<float> amount_glc_mean;
    std::vector<float> amount_cellobiose_mean;
    std::vector<float> amount_xylose_mean;
    std::vector<float> time_mean;
    std::vector<float> Nbr_reactions;

    //Vectors for tracking enzyme concentrations:
    std::vector<float> EG_conc_mean;
    std::vector<float> CBH_conc_mean;
    std::vector<float> BGL_conc_mean;
    std::vector<float> XYL_conc_mean;
    std::vector<float> enzymes_glued;//enzymes glued by lignin

    int N_enzymes_glued;

    //Vectors for counting the times, each enzyme is used as a fraction of overall step number
    std::vector<float> EG_activity;
    std::vector<float> CBH_activity;
    std::vector<float> BGL_activity;
    std::vector<float> XYL_activity;
    std::vector<float> lign_activity;

    //Vectors for tracking the fraction of reactions of different enzyme types among the reaction tables  
    std::vector<float> EG_fraction;
    std::vector<float> CBH_fraction;
    std::vector<float> BGL_fraction;
    std::vector<float> XYL_fraction;   
    std::vector<float> lign_fraction; 

    std::vector<float> timestamp;



    void initialize(){
        min_x_outer = -2;
        max_x_outer = 7;
        min_y_outer = -1;
        max_y_outer = 10;
        Run_number = 0;
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

        nbr_monolignol = 0;
        nbr_lignin_blocked = 0;
        nbr_hemi_monomer = 0;
        pct_lign = 0;
        pct_hemi = 0;
        pct_acetyl_hemi = 0;
        DP_print_Freq = 1;
        enzyme_radius = 0;
        V_enzyme = 0;
    }
};




class neighborList {

//    private:
    public:

        double x,y,z;//Coordinates of bond object
        double dx_mid, dy_mid, dz_mid; //Distance vector components between center of microfibril
        double d_mid; //Actual distance from center of microfibril
        int N_neighbors;//Number of neighbors
        double V_enzyme;//Volume of a single enzyme
        double V_free;//Free volume around bond in radius enzyme
        double r_sphere;//Radius around which neighbors are investigated
        double alpha_max;
        bool outer_bond;//Set to true if this is an outer bond

        std::vector<int> x_neighbors;
        std::vector<int> y_neighbors;
        std::vector<int> z_neighbors;
        std::vector<float> scalar_products;//Projections of the vectors between bond and respective neighbors onto the radial vector outward from the center
//        std::vector<double> dist_neighbors
        std::vector<double> volume_contribution;//volumes of neighbors inside the enzyme sphere


        neighborList(){
            this->x = 0;
            this->y = 0;
            this->z = 0;
            this->dx_mid = 0;
            this->dy_mid = 0;
            this->dz_mid = 0;
            this->d_mid = 0;
            this->N_neighbors = 100;            
            this->V_enzyme = 0;
            this->V_free = 0;
            this->outer_bond = false;
            this->alpha_max = 0.;

            for(int i = 0;i<N_neighbors;i++){
                x_neighbors.push_back(0);
                y_neighbors.push_back(0);                
                z_neighbors.push_back(0);
                volume_contribution.push_back(0);
                scalar_products.push_back(0);
            }

        }
        neighborList(double x, double y, double z, double r_sphere, const params& par, const std::vector<bList>& celluList, const std::vector<bList>& hemiList, const std::vector<bList>& lignList){
            this->x = x;
            this->y = y;
            this->z = z;
            this->dx_mid = x-par.mid_x;
            this->dy_mid = y-par.mid_y;
            this->dz_mid = 0;
            this->d_mid = sqrt(dx_mid*dx_mid + dy_mid*dy_mid + dz_mid*dz_mid);
            this->N_neighbors = 0;            
            this->r_sphere = r_sphere;
            this->V_enzyme = 4.*M_PI*r_sphere*r_sphere*r_sphere/3.;
            this->V_free = V_enzyme;
            this->alpha_max = 0.;
            double x_check= 0;
            double y_check = 0;
            double z_check = 0;

            for(int i=0;i<celluList.size();i++){
                x_check = double(celluList[i].x);
                y_check = double(celluList[i].y);
                for(int j = 0;j < celluList[i].len_poly;j++){
                    z_check = double(celluList[i].z[j]);
                    double dx = this->x - x_check;
                    double dy = this->y - y_check;
                    double dz = this->z - z_check;
                    double d = sqrt(dx*dx + dy*dy + dz*dz);
                    if(d < (par.r_monomer + this->r_sphere) and d != 0){
                        //std::cout << "d = " << d << std::endl;
                        add_neighbor(x_check,y_check,z_check,par);
                    }
                }

            }
            for(int i=0;i<hemiList.size();i++){
                x_check = double(hemiList[i].x);
                y_check = double(hemiList[i].y);
                for(int j = 0;j < hemiList[i].len_poly;j++){
                    z_check = double(hemiList[i].z[j]);
                    double dx = this->x - x_check;
                    double dy = this->y - y_check;
                    double dz = this->z - z_check;
                    double d = sqrt(dx*dx + dy*dy + dz*dz);
                    if(d < (par.r_monomer + this->r_sphere) and d != 0)
                        add_neighbor(x_check,y_check,z_check,par);
                }
                
            }
            for(int i=0;i<lignList.size();i++){
                x_check = double(lignList[i].x);
                y_check = double(lignList[i].y);
                for(int j = 0;j < lignList[i].len_poly;j++){
                    z_check = double(lignList[i].z[j]);
                    double dx = this->x - x_check;
                    double dy = this->y - y_check;
                    double dz = this->z - z_check;
                    double d = sqrt(dx*dx + dy*dy + dz*dz);
                    if(d < (par.r_monomer + this->r_sphere) and d != 0)
                        add_neighbor(x_check,y_check,z_check,par);
                }
            }
            calc_free_vol(par);
            if(par.verbose == true)
                std::cout << "new neighbor number = " << this->N_neighbors << std::endl;
            this->outer_bond = is_outer_bond();
//            for(int i=0; i<N_neighbors; i++){
//                std::cout << this->x-par.mid_x << "\t" << this->y - par.mid_y << std::endl;                
//                std::cout << x_neighbors[i]-this->x << "\t" << y_neighbors[i]-this->y << "\t" << z_neighbors[i]-this->z << "\t" << scalar_products[i] << std::endl;
//            }
//            std::cout << "Coordinates of this bond: " << x << "\t" << y << "\t" << z << std::endl;
//            std::cout << "The fraction of volume which is free is p_free = " << V_free/V_enzyme << "; the number of neighbors is N_neighbors = " << N_neighbors << std::endl; 
        }

        ~neighborList(){
            
        }

        void add_neighbor(int x, int y, int z, const params& par){//Adds a neighbor with coordinates x,y,z

            x_neighbors.push_back(x);
            y_neighbors.push_back(y);
            z_neighbors.push_back(z);
            N_neighbors++;
            float dx = float(x) - float(this->x);
            float dy = float(y) - float(this->y);
            float dz = float(z) - float(this->z);
            float d = sqrt(dx*dx + dy*dy + dz*dz);
            float scalar_prod = 0.;
            volume_contribution.push_back(calc_shared_vol(this->r_sphere,par.r_monomer,d));

            dx /= d;//Calculate projection of distance vector onto radial vector of bond
            dy /= d;
            dz /= d;
//            std::cout << this->d_mid << std::endl;
            if(this->d_mid > 0){
                scalar_prod = dx*this->dx_mid/this->d_mid + dy*this->dy_mid/this->d_mid + dz*this->dz_mid/this->d_mid;
//                if(scalar_prod > 0 and par.verbose == true)
//                    std::cout << scalar_prod << std::endl;
                scalar_products.push_back(scalar_prod);
            }
            else
                scalar_products.push_back(0);
            //if(scalar_products[scalar_products.size()-1] > 0){
//            std::cout << scalar_products[scalar_products.size()-1] << std::endl;
            //}



        }


        void remove_neighbor(int x, int y, int z){//Removes the neighbor with coordinates x,y,z
            if(N_neighbors > 0){
                for(int i=0;i<N_neighbors;i++){
                    if(x_neighbors[i] == x and y_neighbors[i] == y and z_neighbors[i] == z){
                        V_free += volume_contribution[i];
                        x_neighbors.erase(x_neighbors.begin()+i);
                        y_neighbors.erase(y_neighbors.begin()+i);
                        z_neighbors.erase(z_neighbors.begin()+i);
                        volume_contribution.erase(volume_contribution.begin()+i);
                        scalar_products.erase(scalar_products.begin()+i);
                        N_neighbors--;
                    }
                }                
            }
            this->outer_bond = is_outer_bond();
        }

        void remove_neighbor(int i){//Removes the neighbor with index i
//            std::cout << "N_neighbors = " << N_neighbors << "; sizes of lists: " << volume_contribution.size() << "\t" << x_neighbors.size() << "\t" << y_neighbors.size() << "\t" << z_neighbors.size() << "; i = " << i << std::endl;
            if(N_neighbors > 0 and i < N_neighbors){
                V_free += volume_contribution[i];
                x_neighbors.erase(x_neighbors.begin()+i);
                y_neighbors.erase(y_neighbors.begin()+i);
                z_neighbors.erase(z_neighbors.begin()+i);
                volume_contribution.erase(volume_contribution.begin()+i);
                scalar_products.erase(scalar_products.begin()+i);
                N_neighbors--;
            }
            this->outer_bond = is_outer_bond();
        }


        void calc_free_vol(const params& par){
            V_free = V_enzyme;
            for(int i=0;i<N_neighbors;i++){
//                std::cout << "volume_contribution = " << volume_contribution[i] << std::endl;
                V_free -= volume_contribution[i];
            }
            if(V_free < 0)
                V_free = 0;


        }

        double calc_shared_vol(double r1, double r2, double d){
            if(d >= r1+r2)
                return 0.;
            else if(d*d < (r1+r2)*(r1+r2)){
                if (r1<r2)
                    return 4.*M_PI*r1*r1*r1/3.;
                else
                    return 4.*M_PI*r2*r2*r2/3.;
            }
            else{
                return ((M_PI*(r1+r2-d)*(r1+r2-d)*(d*d + 2.*d*r2 - 3.*r2*r2 + 2.*d*r1 + 6.*r2*r1 -3.*r1*r1))/(12.*d));
            }
        }

        bool is_outer_bond(){

        /*    if(par.verbose == true){
                cout << "Entered function is_outer_poly()" << endl;
            }
        */

        if(N_neighbors == 0){
            this->outer_bond = true;
            return true;
        }
        else{
            for(int i=0; i<N_neighbors; i++){
                if(scalar_products[i] > this->alpha_max){
                    this->outer_bond = false;
                    return false;
                }
            }
            this->outer_bond = true;
            return true;
        }

/*

            if(x == par.min_x_overall and y == par.min_y_overall){
                return true;
            }
            if(x == par.min_x_overall and y == par.max_y_overall){
                return true;
            }
            if(x == par.max_x_overall and y == par.min_y_overall){
                return true;
            }
            if(x == par.max_x_overall and y == par.max_y_overall){
                return true;
            }


            int N_directions = 8;

            int distance_x = par.max_x_overall - par.min_x_overall;
            int distance_y = par.max_y_overall - par.min_y_overall;
            int distance = findmax(distance_x,distance_y);

            par.mid_x = par.min_x_overall + 0.5*(par.max_x_overall-par.min_x_overall);
            par.mid_y = par.min_y_overall + 0.5*(par.max_y_overall-par.min_y_overall);
            par.mid_z = this->z;

        //    int stop_count = 0;//Counts, how many neighbors lie in the cardinal directions. If less than 3, it is an outer polymer

            vector<int> dir_x;
            vector<int> dir_y;
            vector<int> dir_z;
            vector<int> unit_dir_x;
            vector<int> unit_dir_y;
            vector<int> unit_dir_z;
            vector<bool> dir_stop;


            for(int i = this->x-1; i <= this->x+1; i++){
                for(int j = this->y-1; j <= this->y+1; j++){
                    for(int k = this->z-1; k <= this->z+1; k++){
                        if(i != this->x or j != this->y or k != this->z)
                        unit_dir_x.push_back(i);
                        unit_dir_y.push_back(j);
                        unit_dir_z.push_back(k);
                        dir_stop.push_back(false);
                    }
                }
            }

            N_directions = dir_stop.size();
            for(int i=0;i<N_directions;i++){
                dir_x.push_back(unit_dir_x[i]);
                dir_y.push_back(unit_dir_y[i]);
                dir_z.push_back(unit_dir_z[i]);
            }


//            if(scalar_prod == len(vec_1)*len(vec_2)) --> parallel

            int neighbor_dir_x = 0;
            int neighbor_dir_y = 0;
            int neighbor_dist_square = 0;
            int scalar_prod = 0;
            int scalar_prod2 = 0;
            int scalar_prod3 = 0;
            float angle = 0;
            vector<float> angles;

            for(int i=0; i<N_neighbors; i++){
                neighbor_dir_x = x_neighbors[i] - this->x;
                neighbor_dir_y = y_neighbors[i] - this->y;
                neighbor_dist_square = neighbor_dir_x * neighbor_dir_x + neighbor_dir_y * neighbor_dir_y;
                for(int j; j<dir_x.size(); j++){
                    scalar_prod = (neighbor_dir_x * (dir_x[j]-this->x) + neighbor_dir_y * (dir_y[j]-this->y));
                    if( scalar_prod*scalar_prod == neighbor_dist_square){//the vectors are parallel
                        dir_stop[j] = true;
                    }
                }
            }



            for(int i=0; i<N_neighbors; i++){
                scalar_prod = (x_neighbors[i]-par.mid_x)*(this->x - par.mid_x) + (y_neighbors[i]-par.mid_y)*(this->y - par.mid_y) + (z_neighbors[i]-par.mid_z)*(this->z - par.mid_z);
                if(scalar_prod > 0){
                    for(int j=0; j<N_neighbors; j++){
                        if(i != j){
                            scalar_prod2 = (x_neighbors[j]-par.mid_x)*(this->x - par.mid_x) + (y_neighbors[j]-par.mid_y)*(this->y - par.mid_y) + (z_neighbors[j]-par.mid_z)*(this->z - par.mid_z);
                            if(scalar_prod2 > 0){
                                scalar_prod3 = (x_neighbors[i]-par.mid_x)*(x_neighbors[j]-par.mid_x) + (y_neighbors[i]-par.mid_y)*(y_neighbors[j]-par.mid_y) + (z_neighbors[i]-par.mid_z)*(z_neighbors[j]-par.mid_z);
                                angle = acos(float(scalar_prod3)/(sqrt(float((x_neighbors[i]-par.mid_x)*(x_neighbors[i]-par.mid_x) + (y_neighbors[i]-par.mid_y)*(y_neighbors[i]-par.mid_y) + (z_neighbors[i]-par.mid_z)*(z_neighbors[i]-par.mid_z)))*sqrt(float((x_neighbors[j]-par.mid_x)*(x_neighbors[j]-par.mid_x) + (y_neighbors[j]-par.mid_y)*(y_neighbors[j]-par.mid_y) + (z_neighbors[j]-par.mid_z)*(z_neighbors[j]-par.mid_z)))));
                                angles.push_back(angle);
                            }
                        }
                    }
                }
            }
            if(angles.size() == 0){
                outer_bond = true;
                return true;
            }
*/

/*
=====================================================================================================================
================================== INSERT CHECK WHETHER THIS IS ALREADY SUFFICIENTLY INSIDE =========================
=====================================================================================================================*/


        //End of function    
        }

};


#endif