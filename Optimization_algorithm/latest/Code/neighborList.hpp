#ifndef NEIGHBORLIST_HXX_
#define NEIGHBORLIST_HXX_
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
#include "bList.hpp"


class neighborList {//This stores the positions of each neighboring bond within a distance around the investigated position

    public:

        double x,y;//Coordinates of bond object
        int z;
        double dx_mid, dy_mid, dz_mid; //Distance vector components between center of microfibril
        double d_mid; //Actual distance from center of microfibril
        int N_neighbors;//Number of neighbors
        double V_enzyme;//Volume of a single enzyme CURRENTLY HAS NO INFLUENCE
        double V_free;//Free volume around bond in radius enzyme CURRENTLY HAS NO INFLUENCE
        double r_sphere;//Radius around which neighbors are investigated
        double alpha_max; //Maximal average scalar product
        double r_monomer;//Radius of a monomer
        bool outer_bond;//Set to true if this is an outer bond

        std::vector<double> x_neighbors;//Stores the x coordinates of all neighbors
        std::vector<double> y_neighbors;//Stores the y coordinates of all neighbors
        std::vector<int> z_neighbors;//Stores the z coordinates of all neighbors
        std::vector<int> material_neighbors;//Stores the material type of all neighbors (cellulose = 1, hemicellulose = 2, lignin = 3)
        std::vector<double> scalar_products;//Projections of the vectors between bond and respective neighbors onto the radial vector outward from the center
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
            this->N_neighbors = 10;            
            this->V_enzyme = 0;
            this->V_free = 0;
            this->outer_bond = false;
            this->alpha_max = 0.;

            for(int i = 0;i<N_neighbors;i++){
                x_neighbors.push_back(0);
                y_neighbors.push_back(0);                
                z_neighbors.push_back(0);
                material_neighbors.push_back(1);
                volume_contribution.push_back(0);
                scalar_products.push_back(0);
            }

        }
        neighborList(double x, double y, double z, double r_sphere, const double mid_x, const double mid_y, double r_monomer,
                     const std::vector<bList>& celluList, const std::vector<bList>& hemiList, const std::vector<bList>& lignList)
        {
            this->x = x;
            this->y = y;
            this->z = z;
            this->dx_mid = x-mid_x;
            this->dy_mid = y-mid_y;
            this->dz_mid = 0;
            this->d_mid = sqrt(dx_mid*dx_mid + dy_mid*dy_mid + dz_mid*dz_mid);
            this->N_neighbors = 0;            
            this->r_sphere = r_sphere;
            this->r_monomer = r_monomer;
            this->V_enzyme = 4.*M_PI*r_sphere*r_sphere*r_sphere/3.;
            this->V_free = V_enzyme;
            this->alpha_max = -0.1;
            double x_check= 0; // Placeholder for the x coordinate of polymer bonds being checked for neighbor-status
            double y_check = 0; // Placeholder for the y coordinate of polymer bonds being checked for neighbor-status
            double z_check = 0; // Placeholder for the u coordinate of polymer bonds being checked for neighbor-status
            double epsilon = 0.001*this->r_sphere; //Used to ensure that we don't add a bond as its own neighbor

            int material = 1;//Check cellu neighbors
            for(int i=0;i<celluList.size();i++){
                x_check = double(celluList[i].x);
                y_check = double(celluList[i].y);
                for(int j = 0;j < celluList[i].len_poly;j++){
                    z_check = double(celluList[i].z[j]);
                    double dx = this->x - x_check;
                    double dy = this->y - y_check;
                    double dz = this->z - z_check;
                    double d = sqrt(dx*dx + dy*dy + dz*dz);
                    if(d < (this->r_monomer + this->r_sphere) and d > epsilon){ 
                        //std::cout << "d = " << d << std::endl;
                        add_neighbor(x_check,y_check,z_check,material);
                    }
                }

            }
            material = 2;//Check hemi neighbors
            for(int i=0;i<hemiList.size();i++){
                x_check = double(hemiList[i].x);
                y_check = double(hemiList[i].y);
                for(int j = 0;j < hemiList[i].len_poly;j++){
                    z_check = double(hemiList[i].z[j]);
                    double dx = this->x - x_check;
                    double dy = this->y - y_check;
                    double dz = this->z - z_check;
                    double d = sqrt(dx*dx + dy*dy + dz*dz);
                    if(d < (this->r_monomer + this->r_sphere) and d > epsilon)
                        add_neighbor(x_check,y_check,z_check,material);
                }
                
            }
            material = 3;//Check lignin neighbors
            for(int i=0;i<lignList.size();i++){
                x_check = double(lignList[i].x);
                y_check = double(lignList[i].y);
                for(int j = 0;j < lignList[i].len_poly;j++){
                    if(lignList[i].covering[j] == true){
                        z_check = double(lignList[i].z[j]);
                        double dx = this->x - x_check;
                        double dy = this->y - y_check;
                        double dz = this->z - z_check;
                        double d = sqrt(dx*dx + dy*dy + dz*dz);
                        if(d < (this->r_monomer + this->r_sphere) and d > epsilon)
                            add_neighbor(x_check,y_check,z_check,material);                        
                    }
                }
            }
            calc_free_vol();
            //if(par.verbose == true)
                //std::cout << "new neighbor number = " << this->N_neighbors << std::endl;
            this->outer_bond = is_outer_bond();
        }


        neighborList(double x, double y, double z, double r_sphere, double r_monomer,
                     const std::vector<bList>& celluList)//Only used for checking the number of free ends surrounding a CBH enzyme
        {
            this->x = x;
            this->y = y;
            this->z = z;
            this->dx_mid = 0.;
            this->dy_mid = 0.;
            this->dz_mid = 0;
            this->d_mid = 0.;
            this->N_neighbors = 0;            
            this->r_sphere = r_sphere;
            this->r_monomer = r_monomer;
            this->V_enzyme = 4.*M_PI*r_sphere*r_sphere*r_sphere/3.;
            this->V_free = V_enzyme;
            this->alpha_max = 0.25;//A value of 0.25 means that the angle between each neighbor and the radial vector cannot be smaller than 75 degrees
            double x_check= 0;
            double y_check = 0;
            double z_check = 0;
            double epsilon = 0.001*this->r_sphere;
            double dx = 0.;
            double dy = 0.;
            double dz = 0.;
            double d = 0.;

            int material = 1;//Check cellu neighbors
            for(int i=0;i<celluList.size();i++){
                x_check = double(celluList[i].x);
                y_check = double(celluList[i].y);
                if(celluList[i].len_poly == 2){
                    z_check = double(celluList[i].z[0]);
                    dx = this->x - x_check;
                    dy = this->y - y_check;
                    dz = this->z - z_check;
                    d = sqrt(dx*dx + dy*dy + dz*dz);
                    if(d < (this->r_monomer + this->r_sphere) and d > epsilon){
                        //std::cout << "d = " << d << std::endl;
                        add_neighbor(x_check,y_check,z_check,material);
                    }
                    z_check = double(celluList[i].z[1]);
                    dx = this->x - x_check;
                    dy = this->y - y_check;
                    dz = this->z - z_check;
                    d = sqrt(dx*dx + dy*dy + dz*dz);
                    if(d < (this->r_monomer + this->r_sphere) and d > epsilon){
                        //std::cout << "d = " << d << std::endl;
                        add_neighbor(x_check,y_check,z_check,material);
                    }

                }
                else if(celluList[i].len_poly == 3){
                    z_check = double(celluList[i].z[1]);
                    dx = this->x - x_check;
                    dy = this->y - y_check;
                    dz = this->z - z_check;
                    d = sqrt(dx*dx + dy*dy + dz*dz);
                    if(d < (this->r_monomer + this->r_sphere) and d > epsilon){
                        //std::cout << "d = " << d << std::endl;
                        add_neighbor(x_check,y_check,z_check,material);
                    }
                }
                else if(celluList[i].len_poly > 3){
                    z_check = double(celluList[i].z[1]);
                    dx = this->x - x_check;
                    dy = this->y - y_check;
                    dz = this->z - z_check;
                    d = sqrt(dx*dx + dy*dy + dz*dz);
                    if(d < (this->r_monomer + this->r_sphere) and d > epsilon){
                        //std::cout << "d = " << d << std::endl;
                        add_neighbor(x_check,y_check,z_check,material);
                    }
                    z_check = double(celluList[i].z[celluList[i].len_poly-2]);
                    dx = this->x - x_check;
                    dy = this->y - y_check;
                    dz = this->z - z_check;
                    d = sqrt(dx*dx + dy*dy + dz*dz);
                    if(d < (this->r_monomer + this->r_sphere) and d > epsilon){
                        //std::cout << "d = " << d << std::endl;
                        add_neighbor(x_check,y_check,z_check,material);
                    }
                }

            }
            this->outer_bond = false;
        }

        ~neighborList(){
            
        }

        void add_neighbor(double x, double y, int z, int mat){//Adds a neighbor with coordinates x,y,z

            this->x_neighbors.push_back(x);
            this->y_neighbors.push_back(y);
            this->z_neighbors.push_back(z);
            this->material_neighbors.push_back(mat);
            this->N_neighbors++;
            double dx = double(x) - double(this->x);
            double dy = double(y) - double(this->y);
            double dz = double(z) - double(this->z);
            double d = sqrt(dx*dx + dy*dy + dz*dz);
            double scalar_prod = 0.;
            volume_contribution.push_back(calc_shared_vol(this->r_sphere,this->r_monomer,d));

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


        void remove_neighbor(double x, double y, int z){//Removes the neighbor with coordinates x,y,z
            if(N_neighbors > 0){
                for(int i=0;i<N_neighbors;i++){
                    if(x_neighbors[i] == x and y_neighbors[i] == y and z_neighbors[i] == z){
                        V_free += volume_contribution[i];
                        x_neighbors.erase(x_neighbors.begin()+i);
                        y_neighbors.erase(y_neighbors.begin()+i);
                        z_neighbors.erase(z_neighbors.begin()+i);
                        material_neighbors.erase(material_neighbors.begin()+i);
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
                material_neighbors.erase(material_neighbors.begin()+i);
                volume_contribution.erase(volume_contribution.begin()+i);
                scalar_products.erase(scalar_products.begin()+i);
                N_neighbors--;
            }
            this->outer_bond = is_outer_bond();
        }


        void calc_free_vol(){
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

            if(N_neighbors == 0){
                this->outer_bond = true;
                return true;
            }
            else{
                double average_scalar_product = 0.;
                for(int i=0; i<this->N_neighbors; i++){
                    average_scalar_product += this->scalar_products[i];
//                    if(this->scalar_products[i] > this->alpha_max){
//                        this->outer_bond = false;
//                        return false;
//                    }
                }
                average_scalar_product/=this->N_neighbors;
                if(average_scalar_product > this->alpha_max){
                    this->outer_bond = false;
                    return false;
                }
                else{
                    this->outer_bond = true;
                    return true;                    
                }
            }

/*
=====================================================================================================================
================================== INSERT CHECK WHETHER THIS IS ALREADY SUFFICIENTLY INSIDE =========================
=====================================================================================================================*/
        //End of function    
        }

        void clear(){

            this->x = -1;
            this->y = -1;
            this->z = -1;
            this->dx_mid = -1;
            this->dy_mid = -1;
            this->dz_mid = -1;
            this->d_mid = -1;
            this->N_neighbors = 0;
            this->V_enzyme = 0;
            this->V_free = 0;
            this->r_sphere = 0;
            this->alpha_max = -1;
            this->outer_bond = 0;


            this->x_neighbors.clear();
            this->y_neighbors.clear();
            this->z_neighbors.clear();
            this->material_neighbors.clear();
            this->volume_contribution.clear();

        }

        bool has_neighbor_material(int mat){//Returns true, if at least one neighbor of this bond is of material mat
            for(int i=0; i < this->N_neighbors; i++){
                if(this->material_neighbors[i] == mat){
                    return true;
                }
            }
            return false;
        //End of function
        }
};


#endif