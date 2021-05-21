#ifndef TLIST_HXX_
#define TLIST_HXX_
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


struct TList//Structure of the reaction table
{
//    int originalTable; //Original table from which this polymer originates 
    int nbr_element;//Number of elements of the structure
    double prop_sum;//Sum of the propensities of all reactions inside this structure
    int index_poly;//Index of the polymer targetted
    std::vector<int> num_bond;//Index of the bond inside the "*z" list of the polymer
    std::vector<int> material;//Indicates 1 if cellulose is the substrate 2 if hemicellulose is the substrate
    std::vector<int> indic_action;//Action label: digestion by EG =1, CBH =2, BGL =3;
    std::vector<double * > liste_prop;//Propensity associated to this action
    std::vector<double * > prop_uninhib; //Propensity without inhibition
    std::vector<bool> crystalline;//Equal to 1 if the reaction acts on a crystalline bond, otherwise equal to 0
    std::vector<bool> covered;//Equal to 1 if the corresponding bond is covered by a CBH enzyme. Only imortant for action_mu 6 reactions 



    void calcTableProp(){//Calculates the starting propensity sum of the table. Only used during initialization
        prop_sum = 0;
        for(int i=0; i<nbr_element;i++){
            if(covered[i] == false){
                prop_sum += *liste_prop[i];
            }
        }
    }
    
    void addProp(double prop){//Increase prop_sum by prop
        prop_sum += prop;
    }

    void subProp(double prop){//Reduce prop_sum by prop
        double prop_sum_before = prop_sum;
        double epsilon = 1e-10*prop_sum_before;//values below this are considered to be 0
//        std::cout << "Before subtraction: prop_sum = " << prop_sum << std::endl;        
        if(prop_sum > 0)
            prop_sum -= prop;
        if(prop_sum < -epsilon){
            std::cout << "NEGATIVE PROPENSITIES OCCURING!!! prop = " << prop << ", prop_sum before: " << prop_sum_before << ", prop_sum now: " << prop_sum << std::endl;
            std::cout << "Reaction dump: " << std::endl;
            print_all_reactions();
            exit(1);
        }

    }

    void print_reaction(int i){
        std::cout << "===================================" << std::endl;
        std::cout << "Reaction " << i << ": " << std::endl;
        std::cout << "poly: " << index_poly << "; bond: " << num_bond[i] << "; indic_action: " << indic_action[i] << std::endl;
        std::cout << "===================================" << std::endl;
    }

    void print_all_reactions(){
        std::cout << "Table number: " << index_poly << "; number of reactions contained in it: " << nbr_element << "; reactions: " << std::endl;
        std::cout << "===================================" << std::endl;
        for(int i=0; i<nbr_element; i++){
        std::cout << "Reaction " << i << ": " << std::endl;
        std::cout << "poly: " << index_poly << "; bond: " << num_bond[i] << "; indic_action: " << indic_action[i] << std::endl;
        std::cout << "===================================" << std::endl;
        }
    }

};

#endif