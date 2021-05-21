#ifndef CBHENZYME_HXX_
#define CBHENZYME_HXX_
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


class CBH_enzyme{//The processive CBH enzymes are implemented as objects with a defined attachment location and duration. They can attach and detach

    public:
        bool attached;
        double attachment_time;
        double attachment_duration;
        int poly_attached;
        int bond_attached;
        double radius;
        neighborList enzyme_neighbors;


        CBH_enzyme(){
            this->attached = false;
            this->attachment_time = -1;
            this->attachment_duration = -1;
            this->poly_attached = -1;
            this->bond_attached = -1;
            this->enzyme_neighbors = neighborList();
            this->radius = 0.;
        }

        CBH_enzyme(double attachment_time, double attachment_duration, int poly_attached, int bond_attached, double radius, const std::vector<bList>& celluList){
            this->attached = true;
            this->attachment_time = attachment_time;
            this->attachment_duration = attachment_duration;
            this->poly_attached = poly_attached;
            this->bond_attached = bond_attached;
            this->enzyme_neighbors = neighborList();            
            this->radius = radius;
        }
        ~CBH_enzyme(){

        }
        void attach(double attachment_time, double attachment_duration, int poly_attached, int bond_attached, const std::vector<bList>& celluList){
            this->attached = true;
            this->attachment_time = attachment_time;
            this->attachment_duration = attachment_duration;
            this->poly_attached = poly_attached;
            this->bond_attached = bond_attached;
            this->enzyme_neighbors = neighborList();            
        }
        void detach(){
            this->attached = false;
            this->attachment_time = -1;
            this->attachment_duration = -1;
            this->poly_attached = -1;
            this->bond_attached = -1;
            this->enzyme_neighbors = neighborList();            
        }
        void init_neighbors(double x, double y, double z, double r_sphere, const double mid_x, const double mid_y, double r_monomer,
                     const std::vector<bList>& celluList, const std::vector<bList>& hemiList, const std::vector<bList>& lignList)
        {
              this->enzyme_neighbors = neighborList(x,y,z,r_sphere,r_monomer,celluList);
//            this->enzyme_neighbors = neighborList(x,y,z,r_sphere,mid_x,mid_y,r_monomer,celluList,hemiList,lignList);
        }
//box_muller(par.average_CBH_association_time*(par.std_dev_CBH_speed/par.average_CBH_speed),par.average_CBH_association_time)
};



#endif