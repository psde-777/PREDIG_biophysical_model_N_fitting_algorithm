#ifndef BLIST_HXX_
#define BLIST_HXX_
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
    std::vector<bool> covering;//Indicates whether this is a bond that covers the microfibril(1), or a "transparent" bond (0). This is only important for lignin, as we only want partial blocking for this
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



#endif