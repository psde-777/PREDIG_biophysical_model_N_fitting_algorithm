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
#include "functions.hpp"
#include "structs.hpp"
#include <unistd.h>
#include <unordered_map>

using namespace std;
using std::ifstream;
using std::ofstream;
using std::ios;



//=========================== Initiation of Structures ====================================
//Initiating the Table list
void initTList(TList& list, int index)
{
    list.nbr_element = 0;
    list.index_poly = index;
    // The vectors do not need to be initialized because the default constructor already provides an empty vector, which we want
//End of function    
}

//======================== Functions ========================================



//Adds a bond to the polymer
void addbond(bList& list,int position,int state, int bond_type, bool crystalline)
{
    list.z.push_back(position);
    list.status.push_back(state);
    list.bond_type.push_back(bond_type);
    list.N_blocked_positions.push_back(0);
    list.crystalline.push_back(crystalline);
    list.len_poly++;
//End of function    
}

//Erase the end of the polymer, until digested bond
void taylorOldPoly(bList& list,int bond_num)
{
    if(list.len_poly > 1){
        if(bond_num >= 0 and bond_num <= list.len_poly-1){
            if(list.len_poly > bond_num){
                list.z.resize(bond_num);
                list.status.resize(bond_num);
                list.bond_type.resize(bond_num);
                list.N_blocked_positions.resize(bond_num);
                list.crystalline.resize(bond_num);
                list.len_poly = bond_num;
            }
            else{
                list.z.clear();
                list.status.clear();
                list.bond_type.clear();
                list.N_blocked_positions.clear();
                list.crystalline.clear();
                list.len_poly = 0;

            }
        }
    }
    else if(list.len_poly == 1){
        list.z.clear();
        list.status.clear();
        list.bond_type.clear();
        list.N_blocked_positions.clear();
        list.crystalline.clear();
        list.len_poly = 0;        
    }
//End of function    
}

//Erase the begining of the polymer, from digested bond
void taylorNewPoly(bList& list,int bond_num)
{
    // erase the first 3 elements:
//    myvector.erase (myvector.begin(),myvector.begin()+3);
//    cout << "Function taylorNewPoly: list.len_poly = " << list.len_poly << "; bond_num = " << bond_num << endl;
    if(list.len_poly > 1){
        if(bond_num == list.len_poly -1){
            list.z.clear();
            list.status.clear();
            list.bond_type.clear();
            list.N_blocked_positions.clear();
            list.crystalline.clear();
            list.len_poly = 0;
        }
        else{
            list.z.erase(list.z.begin(), list.z.begin() + bond_num+1);
            list.status.erase(list.status.begin(), list.status.begin() + bond_num+1);
            list.bond_type.erase(list.bond_type.begin(), list.bond_type.begin() + bond_num+1);
            list.N_blocked_positions.erase(list.N_blocked_positions.begin(), list.N_blocked_positions.begin()+bond_num+1);
            list.crystalline.erase(list.crystalline.begin(), list.crystalline.begin()+bond_num+1);
            list.len_poly = list.len_poly - bond_num - 1;            
        }
    }
    else{
        list.z.erase(list.z.begin());
        list.status.erase(list.status.begin());
        list.bond_type.erase(list.bond_type.begin());
        list.N_blocked_positions.erase(list.N_blocked_positions.begin());
        list.crystalline.erase(list.crystalline.begin());
        list.len_poly = 0;  
    }


    
//End of function    
}

//Adds the new possible reaction to the reactions Table... if confidence == false, the table is checked for the existence of the reaction before it is added
void addreaction(bool confidence, TList& list,int poly,int bond_num,int mat,float act2,double prope)
{
    if(prope > 0){

        if(confidence == false){
            int test = 0;

            for(int i = 0;i<list.nbr_element;i++){
                if(list.num_bond[i] == bond_num and list.material[i] == mat and list.indic_action[i] == act2)
                    test = 1;
                if (test == 1 & list.prop_uninhib[i] != prope)
                    list.liste_prop[i] = prope;
                    list.prop_uninhib[i] = prope;
            }

            if(test == 0){
                list.num_bond.push_back(bond_num);
                list.material.push_back(mat);
                list.indic_action.push_back(act2);
                list.liste_prop.push_back(prope);
                list.prop_uninhib.push_back(prope);
                list.addProp(prope);

                list.nbr_element++;        
            }

        }
        else{
            list.num_bond.push_back(bond_num);
            list.material.push_back(mat);
            list.indic_action.push_back(act2);
            list.liste_prop.push_back(prope);
            list.prop_uninhib.push_back(prope);
            list.addProp(prope);

            list.nbr_element++;                
        }
    }


//End of function    

}


//Fills the reaction table of a newly created polymer
void fill_table(params& par, vector<TList>& allTables, vector<bList>& polyList, int poly_selected, int substrate, bool& error, vector<float>& chem_entities){
    int leftBond = 0;
    int rightBond = 0;
    double propensite = 0;
    double crystal_propensite = 0;

    int reaction_type = 0;
    if(substrate == 1){
        if(polyList[poly_selected].len_poly == 0){
            cout << "Polymer is empty! Stopping" << endl;
            exit(1);
            return;
        }
        else if(polyList[poly_selected].len_poly == 1){
            reaction_type = 3;
            propensite = prop(par, reaction_type, chem_entities);
            crystal_propensite = par.crystal_modifier_cellu*propensite;

            if(propensite > 0 and polyList[poly_selected].status[0] == 1){
                if(polyList[poly_selected].crystalline[0] == false)
                    addreaction(1, allTables[poly_selected],poly_selected,0,substrate,reaction_type,propensite);
                else
                    addreaction(1, allTables[poly_selected],poly_selected,0,substrate,reaction_type,crystal_propensite);
            }
        }
        else if(polyList[poly_selected].len_poly > 1){
            reaction_type = 1;
            propensite = prop(par, reaction_type, chem_entities);
            crystal_propensite = par.crystal_modifier_cellu*propensite;
//            cout << propensite << "\t" << crystal_propensite << endl;

            if(propensite > 0){
                if(polyList[poly_selected].len_poly == 2){//Digestion of this polymer length by EG is currently disabled
/*                    reaction_type = 1;
                    propensite = prop(par, reaction_type, chem_entities);
                    crystal_propensite = par.crystal_modifier_cellu*propensite;
                    if(polyList[poly_selected].status[0] == 1)
                        addreaction(1, allTables[poly_selected],poly_selected,0,substrate,reaction_type,propensite);
                    if(polyList[poly_selected].status[1] == 1)                    
                        addreaction(1, allTables[poly_selected],poly_selected,1,substrate,reaction_type,propensite);*/
                }
                else{
                    for(int i = 2;i<polyList[poly_selected].len_poly-2;i++){
                        if(polyList[poly_selected].status[i] == 1){
                            if(polyList[poly_selected].crystalline[i] == false)
                                addreaction(1, allTables[poly_selected],poly_selected,i,substrate,reaction_type,propensite);
                            else
                                addreaction(1, allTables[poly_selected],poly_selected,i,substrate,reaction_type,crystal_propensite);
                        }
                    }   
                }
            }
            reaction_type = 2;
            propensite = prop(par, reaction_type, chem_entities);
            crystal_propensite = par.crystal_modifier_cellu*propensite;

            if(propensite > 0){
                if(polyList[poly_selected].len_poly == 2){
                    leftBond = 0;
                    rightBond = 1;
                }
                else if(polyList[poly_selected].len_poly > 2){
                    leftBond = 1;
                    rightBond = polyList[poly_selected].len_poly-2;   
                }
                if(polyList[poly_selected].status[leftBond] == 1){
                    if(polyList[poly_selected].crystalline[leftBond] == false)
                        addreaction(1, allTables[poly_selected],poly_selected,leftBond,substrate,reaction_type,propensite);
                    else
                        addreaction(1, allTables[poly_selected],poly_selected,leftBond,substrate,reaction_type,crystal_propensite);
                }
                if(polyList[poly_selected].status[rightBond] == 1){
                    if(polyList[poly_selected].crystalline[rightBond] == false)
                        addreaction(1, allTables[poly_selected],poly_selected,rightBond,substrate,reaction_type,propensite);
                    else
                        addreaction(1, allTables[poly_selected],poly_selected,rightBond,substrate,reaction_type,crystal_propensite);
                }
            }
        }

        for(int i=0;i<allTables[poly_selected].nbr_element;i++){
            if(polyList[poly_selected].len_poly != 2 and (allTables[poly_selected].num_bond[i] == 0 or allTables[poly_selected].num_bond[i] == polyList[poly_selected].len_poly-1) and (allTables[poly_selected].indic_action[i] == 1 or allTables[poly_selected].indic_action[i] == 2)){
                cout << "Filling of cellu tables is making problems: " << endl;
                exit(1);
                return; 
            }
        }

    }
    else if(substrate == 2){
        reaction_type = 4;
        propensite = prop(par, reaction_type, chem_entities);
        crystal_propensite = par.crystal_modifier_hemi*propensite;
        
        for(int i=0;i<polyList[poly_selected].len_poly;i++){
            if(polyList[poly_selected].bond_type[i] == 4 and propensite > 0 and polyList[poly_selected].status[i] == 1){
                if(polyList[poly_selected].crystalline[i] == false)
                    addreaction(1,allTables[poly_selected],poly_selected,i,substrate,reaction_type,propensite);
                else
                    addreaction(1,allTables[poly_selected],poly_selected,i,substrate,reaction_type,crystal_propensite);
            }
        }
    }
//End of function    
}

//Deletes the reaction labelled "action_target" for the polymer "poly" for its bond "bond_target" of the material "mat" in the list of reactions
void deletereaction(TList& list,vector<TList>& Table,int mat,int poly,int bond_target,int action_target)
{

    for(int i=0; i<list.nbr_element; i++) {

        if((list.num_bond[i]==bond_target) and (list.indic_action[i]==action_target) and (list.material[i]==mat)) {// So basically if the reaction fits all the function arguments

            list.num_bond.erase(list.num_bond.begin() + i);
            list.material.erase(list.material.begin() + i);
            list.indic_action.erase(list.indic_action.begin() + i);

            list.subProp(list.liste_prop[i]);

            list.liste_prop.erase(list.liste_prop.begin() + i);
            list.prop_uninhib.erase(list.prop_uninhib.begin() + i);
            list.nbr_element--;
        }
    }
    if(list.liste_prop.size() == 0)
        list.prop_sum = 0;
//End of function    
}
void delete_specific_reaction(TList& list,vector<TList>& Table,int i){
    if(list.nbr_element > 0){
//        cout << "deletin; i = " << i << "; size of table: " << list.nbr_element << endl;
        list.num_bond.erase(list.num_bond.begin() + i);
        list.material.erase(list.material.begin() + i);
        list.indic_action.erase(list.indic_action.begin() + i);

        list.subProp(list.liste_prop[i]);

        list.liste_prop.erase(list.liste_prop.begin() + i);
        list.prop_uninhib.erase(list.prop_uninhib.begin() + i);
        list.nbr_element--;
    }

    if(list.liste_prop.size() == 0)
        list.prop_sum = 0;
//End of function    
}

//Deletes ALL reactions for the bond "bond_target" from the polymer "poly" of the material "mat" in the list of reactions
void deleteAllreaction(TList& list,vector<TList>& Table,int mat,int poly,int bond_target)
{
    for(int i=0; i<list.nbr_element; i++) {

        if((list.num_bond[i]==bond_target) and (list.material[i]==mat)) {// So basically if the reaction fits all the function arguments

            list.num_bond.erase(list.num_bond.begin() + i);
            list.material.erase(list.material.begin() + i);
            list.indic_action.erase(list.indic_action.begin() + i);

            list.subProp(list.liste_prop[i]);
//cout << list.liste_prop.size() << "\t" << list.prop_uninhib.size() << endl;
            list.liste_prop.erase(list.liste_prop.begin() + i);
            list.prop_uninhib.erase(list.prop_uninhib.begin() + i);
//    cout << "searching in functions.cpp" << endl;
            list.nbr_element--;

            i -= 1;//To be sure we scan the right hand side neighbor in the list that just shift from 1 rank
        }
    }    
    if(list.liste_prop.size() == 0)
        list.prop_sum = 0;

//End of function    
}


//Computes the propensity for each reaction to take place
double prop(params& par, float act, vector<float> chem_entities)//For now hemicellulose and cellulose are digested with same rates, but could easily change
{

    double propens;

    if (act==1)//Degradation by EG
        propens=par.k1*chem_entities.at(0);
    else if (act==2)//Degradation by CBH
        propens=par.k2*chem_entities.at(1);
    else if (act==3)//Degradation by BGL
        propens=par.k3*chem_entities.at(2);
    else if (act==4)//Degradation of hemi by XYL
        propens = par.k4*chem_entities.at(3);
    else if (act==5)//Degradation of hemi by XYL
        if (par.mode_lignin_glue == 1)
            propens = par.k5*chem_entities.at(4)*(chem_entities[0]+chem_entities[1]+chem_entities[2]+chem_entities[3]);
        else
            propens = 0;
    else{
        cout << "In function prop(): variable 'act' does not correspond to any of the enzymes included in this simulation; act = " << act << endl;
        return 0;
    }
//    cout << "In function prop(): act = " << act << "; propens = " << propens << endl;
    return propens;
//End of function    
}

//Counts the amount of glucose in the vector list. Depending on the substrate type (cellulose or hemicellulose), xylose and glucose need to be distinguished
double countGlc(const vector<bList>& list, const int substrate){
    int count = 0;
    if(substrate == 1){    
        for(int i = 0; i < list.size(); i++){
            if(list[i].len_poly > 0){
                count+= (list[i].len_poly+1);
            }
        }
    }
    else if(substrate == 2){
        for(int i = 0; i < list.size();i++){
            if(list[i].len_poly > 0){
                if(list[i].bond_type[0] == 1 or list[i].bond_type[0] == 3)
                    count++;
                for(int j=0;j < list[i].bond_type.size();j++){
                    if(list[i].bond_type[j] == 1 or list[i].bond_type[j] == 2){
                        count++;
                    }
                }
            }
        }
    }
    return count;
//End of function    
}


//Counts the amount of xylose in the vector list
double countXyl(const vector<bList>& list){
    int count = 0;
    for(int i = 0; i < list.size();i++){
        if(list[i].len_poly > 0){
            if(list[i].bond_type[0] == 2 or list[i].bond_type[0] == 4)
                count++;
            for(int j=0;j < list[i].bond_type.size();j++){
                if(list[i].bond_type[j] == 3 or list[i].bond_type[j] == 4){
                    count++;
                }
            }
        }
    }
    return count ;
//End of function    
}

double countLign(const std::vector<bList>& list){
    int count = 0;
    for(int i = 0; i<list.size();i++){
        count += list[i].len_poly;
    }
    return count;
//End of function    
}

int countCellobiose(const vector<bList>& list){
    int count = 0;
    for(int i=0; i<list.size();i++){
        if(list[i].len_poly == 1)
            count ++;
    }
    return count;
//End of function    
}

int countReactions(TList& list, int poly_selected, bool& error){
    if(list.nbr_element == list.liste_prop.size())
        return list.nbr_element;
    else{
        cout << "In function countReactions(): list.nbr_element != list.liste_prop.size()" << endl;
        exit(1);
        return 0;
    }
//End of function    
}

void moveReactions(TList& oldList, TList& newList, bList& structOldPoly, bList& structNewPoly, int oldpoly, int bond_selected, int newPoly){

    int count = 0;

    for(int i=0;i<oldList.nbr_element;i++){


        if(oldList.num_bond[i] > bond_selected){
            count++;
            newList.nbr_element++;
            newList.num_bond.push_back(oldList.num_bond[i]-bond_selected-1);
            newList.material.push_back(oldList.material[i]);
            newList.indic_action.push_back(oldList.indic_action[i]);
            newList.liste_prop.push_back(oldList.liste_prop[i]);
            newList.prop_uninhib.push_back(oldList.prop_uninhib[i]);

            oldList.nbr_element--;
            oldList.subProp(oldList.liste_prop[i]);
            oldList.num_bond.erase(oldList.num_bond.begin()+i);
            oldList.material.erase(oldList.material.begin()+i);
            oldList.indic_action.erase(oldList.indic_action.begin()+i);
            oldList.liste_prop.erase(oldList.liste_prop.begin()+i);
            oldList.prop_uninhib.erase(oldList.prop_uninhib.begin()+i);            

            i--; //move i one to the left so as to avoid skipping the next neighbor

        }
        else if(oldList.num_bond[i] == bond_selected){
            oldList.nbr_element--;
            oldList.subProp(oldList.liste_prop[i]);
            oldList.num_bond.erase(oldList.num_bond.begin()+i);
            oldList.material.erase(oldList.material.begin()+i);
            oldList.indic_action.erase(oldList.indic_action.begin()+i);
            oldList.liste_prop.erase(oldList.liste_prop.begin()+i);
            oldList.prop_uninhib.erase(oldList.prop_uninhib.begin()+i);            

            i--; //move i one to the left so as to avoid skipping the next neighbor            
        }

    }


    newList.calcTableProp();
    if(newList.liste_prop.size() != count){
        while(true){
            cout << "In function moveReactions: Too many or too few elements were moved" << endl;
            cout << "(This message plays on a loop; press ctrl-c to stop)" << endl;
            usleep(1000000);                
        }        
    }

    for(int i=0;i<oldList.nbr_element;i++){
        if(oldList.num_bond[i] >= bond_selected){
            while(true){
                cout << "In function moveReactions: not all bonds were successfully moved!" << endl;
                cout << "(This message plays on a loop; press ctrl-c to stop)" << endl;
                usleep(1000000);                

            }
        }
    }
    for(int i=0;i<newList.nbr_element;i++){
        if(newList.num_bond[i] >= structNewPoly.len_poly){
            while(true){
                cout << "In function moveReactions: A bond in a new reaction table is too large for its polymer!" << endl;
                cout << "(This message plays on a loop; press ctrl-c to stop)" << endl;
                usleep(1000000);                

            }
        }
    }    
//End of function    
}

//int countSubst(std::vector<TList> *list, int mat){

//}

void drawPolys(std::vector<bList>& list){
    cout << " ======= Drawing polymers ======== " << endl;

    for(int i=0;i<list.size();i++){
        cout << "POLYMER NUMBER " << i << endl;
        if(list[i].len_poly > 0)
            cout << "glc" ;
        for(int j=0;j<list[i].len_poly;j++){
            cout << " - glc";
        }
        cout << endl;
    }

    cout << " ====== Done drawing polymers ====== " << endl;
//End of function    
}

int findIndex(TList& list, double suma, int test){

    suma -= list.prop_sum;
    int j = 0;

    if(list.liste_prop.size() == 0){
        cout << "In function findIndex(): table contains no reactions!" << endl;
        return -1;
    }
    suma += list.liste_prop[j];
    while(test > suma){
        j++;
//            cout << list.liste_prop[j] << endl;
        suma += list.liste_prop[j];
    }
    return j;
//End of function    
}


void EG_digest(params& par, vector<TList>& allTables, vector<bList>& polyList, int& nbr_poly, int& nbr_xyl_pdt, int& nbr_Glc_pdt,const int bond_selected,const int poly_selected, int len_polyLoopStart, int& nbr_cellobiose, vector<float>& chem_entities, const int substrate,unordered_map<int,neighborList>& bond_neighbors_cellu, unordered_map<int,neighborList>& bond_neighbors_hemi, const bool verbose, bool& error){//EG digestion as single function
    
    if(bond_selected == 0 and polyList[poly_selected].len_poly != 2){
        cout << "EG action on first bond. This should be disabled!" << endl;
        exit(1);
        return;
    }
    if(bond_selected == polyList[poly_selected].len_poly-1 and polyList[poly_selected].len_poly != 2){
        cout << "EG action on last bond. This should be disabled!" << endl;
        exit(1);
        return;
    }
    if(polyList[poly_selected].len_poly == 1){
        cout << "EG is digesting cellobiose! This should be disabled! Stopping..." << endl;
        exit(1);
        return;
    }


    if(polyList[poly_selected].len_poly <= bond_selected){
        cout << "In function EG_digest: the bond of the current reaction is larger than the size of the polymer!" << endl;
        exit(1);
        return;
    }

    int check_radius = 2;
    int dx = 1;//Distance of bonds in coordinates
    int dy = 1;
    int dz = 1;

    if(par.enzyme_radius > 0){
        check_radius = 2*par.enzyme_radius;//Number of bonds which are checked around each position
    }

    int cut_applicate = polyList[poly_selected].z[bond_selected];


    int x_pos = polyList[poly_selected].x;//Exact coordinates of the bond
    int y_pos = polyList[poly_selected].y;
    int z_pos = cut_applicate;
    int neighbor_key = cantor_pair_three(x_pos, y_pos, z_pos);
    int test1; //Used for finding specific reactions in a table
    int indic_cut;//
    int new_poly_index = -1;
    int minSize_BGL = 8;//DP at which BGL starts digesting 
    int x_bond,y_bond,z_bond;// For later freeing of lower bonds

    double propensite; //
    double crystal_propensite;

    int reaction_type = 0;

    bool newPolyFlag = false; //Checks whether a new polymer is generated here

    x_bond = polyList[poly_selected].x;
    y_bond = polyList[poly_selected].y;
    z_bond = polyList[poly_selected].z[bond_selected];

    //Remove all reactions associated to the bond just digested
    deleteAllreaction(allTables[poly_selected],allTables,substrate,poly_selected,bond_selected);


    //======================= Remove the bond and split the polymer into two new polymers =============================

    if(polyList[poly_selected].len_poly==0){
        //Do nothing
    }
    else if(polyList[poly_selected].len_poly == 1){
        //Do nothing
    }
    else if(polyList[poly_selected].len_poly==2)
    {

        if(bond_selected==0)// If the cut was at the beginning of the polymer
        {                        
            taylorNewPoly(polyList[poly_selected],bond_selected);//Shift the indices by one to the left
            allTables[poly_selected].num_bond.clear();
            allTables[poly_selected].material.clear();
            allTables[poly_selected].indic_action.clear();
            allTables[poly_selected].liste_prop.clear();
            allTables[poly_selected].prop_uninhib.clear();
            allTables[poly_selected].nbr_element = 0;            
            

            reaction_type = 3;
            propensite = prop(par,reaction_type,chem_entities);
            crystal_propensite = par.crystal_modifier_cellu*propensite;
            if(propensite > 0 and polyList[poly_selected].status[0] == 1 and polyList[poly_selected].len_poly == 1){
                
                addreaction(0, allTables[poly_selected],poly_selected,0,substrate,reaction_type,propensite);
            }


            nbr_Glc_pdt++;//We release a Glc
        }
        else if(bond_selected == 1){
            taylorOldPoly(polyList[poly_selected],bond_selected);//Erase the end of the polymer, until digested bond
            allTables[poly_selected].num_bond.clear();
            allTables[poly_selected].material.clear();
            allTables[poly_selected].indic_action.clear();
            allTables[poly_selected].liste_prop.clear();
            allTables[poly_selected].prop_uninhib.clear();
            allTables[poly_selected].nbr_element = 0;            

            reaction_type = 3;
            propensite = prop(par,reaction_type,chem_entities);
            crystal_propensite = par.crystal_modifier_cellu*propensite;
            if(propensite > 0 and polyList[poly_selected].status[0] == 1 and polyList[poly_selected].len_poly == 1){
                if(polyList[poly_selected].crystalline[0]==false)
                    addreaction(0, allTables[poly_selected],poly_selected,0,substrate,reaction_type,propensite);
                else
                    addreaction(0, allTables[poly_selected],poly_selected,0,substrate,reaction_type,crystal_propensite);
            }
    
            nbr_Glc_pdt++;//We release a Glc            
        }
    }
    else if(polyList[poly_selected].len_poly>2){
        if((bond_selected !=0) and (bond_selected !=polyList[poly_selected].len_poly-1))
        {

        //====================== CREATE AND INITIALIZE THE NEW POLYMER AND ITS REACTION TABLE =========================0

            polyList.push_back(bList());
            nbr_poly++;            
            new_poly_index = polyList.size()-1;
            newPolyFlag = true;


            polyList[new_poly_index].index= new_poly_index;//Specify the new index of the new polymer
            polyList[new_poly_index].x= polyList[poly_selected].x;
            polyList[new_poly_index].y= polyList[poly_selected].y;
            polyList[new_poly_index].set_z(bond_selected+1);

            allTables.push_back(TList());//Make the new reaction list associated to the new polymer;
            initTList(allTables[new_poly_index], new_poly_index);//Initialize it                        

            for(int i=0; i<polyList[poly_selected].len_poly; i++)//Make the new polymer into a copy of the old one before tayloring
            {
                addbond(polyList[new_poly_index],i,1,1,false);
                polyList[new_poly_index].z[i]=polyList[poly_selected].z[i];
                polyList[new_poly_index].status[i]=polyList[poly_selected].status[i];
                polyList[new_poly_index].bond_type[i]=polyList[poly_selected].bond_type[i];
                polyList[new_poly_index].N_blocked_positions[i] = polyList[poly_selected].N_blocked_positions[i];
                polyList[new_poly_index].crystalline[i] = polyList[poly_selected].crystalline[i];
            }

            taylorOldPoly(polyList[poly_selected],bond_selected);//Erase the end of the old polymer, until digested bond
            taylorNewPoly(polyList[new_poly_index],bond_selected);//Erase the begining of the new polymer, from digested bond



            allTables[poly_selected].num_bond.clear();
            allTables[poly_selected].material.clear();
            allTables[poly_selected].indic_action.clear();
            allTables[poly_selected].liste_prop.clear();
            allTables[poly_selected].prop_uninhib.clear();
            allTables[poly_selected].nbr_element = 0;
            fill_table(par,allTables,polyList,poly_selected,substrate,error,chem_entities);

            //Fill new table
            fill_table(par,allTables, polyList, new_poly_index, substrate, error,chem_entities);



        }

    }

    //====================================== Modify propensities; this is only used when the enzyme size is inactive =====================
    if(par.mode_enzyme_size == -1){
    
        int poly_index = 0;
        int z_index = 0;
        int found_material = 0;
        vector<bList> hemiList,lignList;

        hemiList.push_back(bList());
        hemiList[0].index= 0;//Specify the new index of the new polymer
        hemiList[0].x= 9999;
        hemiList[0].y= 9999;
        hemiList[0].len_poly = 0;
        hemiList[0].set_z(9999);

        lignList.push_back(bList());
        lignList[0].index= 0;//Specify the new index of the new polymer
        lignList[0].x= 9999;
        lignList[0].y= 9999;
        lignList[0].len_poly = 0;
        lignList[0].set_z(9999);


        if(bond_neighbors_cellu.count(neighbor_key) != 1){
            cout << "In function EG_digest: a key is found more than once in bond_neighbors_cellu. This should never happen. Stopping." << endl;
            exit(1);
        }
        else if(bond_neighbors_cellu[neighbor_key].N_neighbors != bond_neighbors_cellu[neighbor_key].x_neighbors.size()){
            cout << "In function EG_digest: a member of bond_neighbors_cellu has an N_neighbors parameter which differs from the size of the corresponding vector. This should never happen. Stopping." << endl;
            exit(1); 
        }
        for(int i=0; i<bond_neighbors_cellu[neighbor_key].N_neighbors; i++){
            if(find_specific_bond(polyList, hemiList, lignList, bond_neighbors_cellu[neighbor_key].x_neighbors[i], bond_neighbors_cellu[neighbor_key].y_neighbors[i], bond_neighbors_cellu[neighbor_key].z_neighbors[i], poly_index, z_index, found_material) == true){
                if(found_material == 1){
    //                cout << "Found" << endl;
                    if(polyList[poly_index].status[z_index] == -1){
                        polyList[poly_index].status[z_index] = 1;
                        if(polyList[poly_index].len_poly == 1){
                            reaction_type = 3;
                            propensite = prop(par,reaction_type,chem_entities);
                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                            if(propensite > 0){
                                if(polyList[poly_index].crystalline[z_index] == false){
                                    addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                }
                                else{
                                    addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                }
                            }
                        }
                        else if(polyList[poly_index].len_poly == 2){
                            reaction_type = 2;
                            propensite = prop(par,reaction_type,chem_entities);
                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                            if(propensite > 0){
                               if(polyList[poly_index].crystalline[z_index] == false){
                                    addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                }
                                else{
                                    addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                }                            
                            }
                        }
                        else if(polyList[poly_index].len_poly > 2){
                            if(z_index > 1 and z_index < polyList[poly_index].len_poly-2){
                                reaction_type = 1;
                                propensite = prop(par,reaction_type,chem_entities);
                                crystal_propensite = par.crystal_modifier_cellu * propensite;
                                if(propensite > 0){
                                   if(polyList[poly_index].crystalline[z_index] == false){
                                        addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                    }
                                    else{
                                        addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                    }                            
                                }

                            }
                            if(z_index == 1 or z_index == polyList[poly_index].len_poly-2){
                                reaction_type = 2;
                                propensite = prop(par,reaction_type,chem_entities);
                                crystal_propensite = par.crystal_modifier_cellu * propensite;
                                if(propensite > 0){
                                   if(polyList[poly_index].crystalline[z_index] == false){
                                        addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                    }
                                    else{
                                        addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                    }                            
                                }
                            }

                        }
                    }

                }
            }
        }
    }


/*
    for(int i=0; i<nbr_poly; i++)//We only scan cellulose because lignin is anyway undigestable and hemicellulose availability is fixed at the begining because it's made of one external layer
    {
        if(polyList[i].len_poly>0)//THIS WAS 0; TESTING!!!
        {
            if(check_position(polyList[i].x, polyList[i].y, polyList[poly_selected].x, polyList[poly_selected].y)==true)
            {
                for(int j=0; j<polyList[i].len_poly; j++)
                {
                    if((polyList[i].z[j]==cut_applicate) or (polyList[i].z[j]==cut_applicate-1) or (polyList[i].z[j]==cut_applicate+1) or ((bond_selected == 1 or bond_selected == polyList[poly_selected].len_poly + polyList[new_poly_index].len_poly -1) and ((polyList[i].z[j]==cut_applicate-2) or (polyList[i].z[j]==cut_applicate+2))) )
                    {
                        polyList[i].status[j]=1;
                        reaction_type = 1;
                        propensite = prop(par,reaction_type,chem_entities);                        
                        crystal_propensite = par.crystal_modifier_cellu*propensite;
                        if(propensite > 0 and j != 0 and j != polyList[i].len_poly-1){
                            if(polyList[i].crystalline[j] == false)
                                addreaction(0, allTables[i],i,j,substrate,reaction_type,propensite);
                            else{
//                                cout << "Adding crystal reaction" << endl;
                                addreaction(0, allTables[i],i,j,substrate,reaction_type,crystal_propensite);                            
                            }
 
                        }
                        reaction_type = 2;
                        propensite = prop(par,reaction_type,chem_entities);                        
                        crystal_propensite = par.crystal_modifier_cellu*propensite;                        

                        if(propensite > 0 and polyList[i].len_poly >= 2){
                            if(polyList[i].len_poly == 2){
                                if(j == 0 or j==1){
                                    if(polyList[i].crystalline[j] == false)
                                        addreaction(0, allTables[i],i,j,substrate,reaction_type,propensite);
                                    else
                                        addreaction(0, allTables[i],i,j,substrate,reaction_type,crystal_propensite);
                                }
                            }
                            else if(polyList[i].len_poly > 2){
                                if(j == 1 or j == polyList[i].len_poly-2){
                                    if(polyList[i].crystalline[j] == false)
                                        addreaction(0, allTables[i],i,j,substrate,reaction_type,propensite);
                                    else
                                        addreaction(0, allTables[i],i,j,substrate,reaction_type,crystal_propensite);
                                }
                            }
                        }
                        reaction_type = 3;
                        propensite = prop(par,reaction_type,chem_entities);                        
                        crystal_propensite = par.crystal_modifier_cellu*propensite;
                        if(propensite > 0 and polyList[i].len_poly == 1){
                            if(polyList[i].crystalline[0] == false)
                                addreaction(0, allTables[i],i,j,substrate,reaction_type,propensite);
                            else
                                addreaction(0, allTables[i],i,j,substrate,reaction_type,crystal_propensite);
                        }
                    }
                }
            }
        }
    }*/

    if(par.mode_enzyme_size == -1){
        remove_neighbor_from_vector(x_pos,y_pos,z_pos,bond_neighbors_cellu,bond_neighbors_hemi);
    }

//    cout << "Before EG " << bond_neighbors_cellu.size() << endl;
//    remove_neighbor_from_vector(x_pos,y_pos,z_pos,bond_neighbors_cellu,bond_neighbors_hemi);
//    cout << "After EG " << bond_neighbors_cellu.size() << endl;

    if(polyList[poly_selected].len_poly == 1){
        nbr_cellobiose++;
    }
    if(newPolyFlag == true){
        if(polyList[new_poly_index].len_poly == 1){
            nbr_cellobiose++;
        }
    }
//End of function    
}

void CBH_digest(params& par, vector<TList>& allTables, vector<bList>& polyList, int& nbr_poly, int& nbr_xyl_pdt, int& nbr_Glc_pdt,const int bond_selected,const int poly_selected, int len_polyLoopStart, int& nbr_cellobiose,vector<float>& chem_entities, const int substrate, unordered_map<int,neighborList>& bond_neighbors_cellu, unordered_map<int,neighborList>& bond_neighbors_hemi,const bool verbose, bool& error){

    if(polyList[poly_selected].len_poly <= bond_selected){
        cout << "In function CBH_digest: the bond of the current reaction is larger than the size of the polymer! bond_selected: " << bond_selected << "; Size of polymer: " << polyList[poly_selected].len_poly << endl;
        exit(1);
        return;
    }
    if(polyList[poly_selected].len_poly < 2){
        cout << "something smaller than len_poly = 2 is being digested. stopping..." << endl;
        exit(1);
        return;
    }

    int check_radius = 2;
    int dx = 1;//Distance of bonds in coordinates
    int dy = 1;
    int dz = 1;
    int x_pos = polyList[poly_selected].x;//Exact coordinates of the bond
    int y_pos = polyList[poly_selected].y;
    int z_pos = polyList[poly_selected].z[bond_selected];
    int neighbor_key = cantor_pair_three(x_pos, y_pos, z_pos);
    if(par.enzyme_radius > 0){
        check_radius = 2*par.enzyme_radius;//Number of bonds which are checked around each position
    }

    int cut_applicate=polyList[poly_selected].z[bond_selected];
    int test1;
    int indic_cut;
    int new_poly_index = -1;//index of new polymer, if cellobiose is created
        

    int reaction_type = 0;
    double propensite;
    double crystal_propensite;

    int newPolyFlag = false;

    if(bond_selected==1){//If the cut takes place at the head of the polymer
        indic_cut=0;
    }
    else if(bond_selected != 0 and bond_selected==polyList[poly_selected].len_poly-2)//If the cut takes place at the tip of the polymer
        indic_cut=1;



    //Remove the bond and split the polymer into two new polymers
    if((bond_selected !=0) and (bond_selected !=polyList[poly_selected].len_poly-1))
    {
        polyList.push_back(bList());//Make the new polymer of cellobiose
        new_poly_index = polyList.size()-1;
        newPolyFlag = true;
        polyList[new_poly_index].index= new_poly_index;//Specify the new index of the new polymer
        polyList[new_poly_index].x= polyList[poly_selected].x;
        polyList[new_poly_index].y= polyList[poly_selected].y;
        polyList[new_poly_index].set_z(bond_selected+1);
        allTables.push_back(TList());//Make the new reaction list associated to the new polymer;
        initTList(allTables[new_poly_index], new_poly_index);//list.originalTable);
        for(int i=0; i<polyList[poly_selected].len_poly; i++)
        {
            addbond(polyList[new_poly_index],i,1,1,false);
            polyList[new_poly_index].z[i]=polyList[poly_selected].z[i];
            polyList[new_poly_index].status[i]=polyList[poly_selected].status[i];
            polyList[new_poly_index].bond_type[i]=polyList[poly_selected].bond_type[i];
            polyList[new_poly_index].N_blocked_positions[i] = polyList[poly_selected].N_blocked_positions[i];
            polyList[new_poly_index].crystalline[i] = polyList[poly_selected].crystalline[i];
        }


        taylorOldPoly(polyList[poly_selected],bond_selected);//Erase the end of the polymer, until digested bond
        taylorNewPoly(polyList[new_poly_index],bond_selected);//Erase the begining of the polymer, from digested bond

        //Adjust old table
        allTables[poly_selected].num_bond.clear();
        allTables[poly_selected].material.clear();
        allTables[poly_selected].indic_action.clear();
        allTables[poly_selected].liste_prop.clear();
        allTables[poly_selected].prop_uninhib.clear();
        allTables[poly_selected].nbr_element = 0;
        fill_table(par,allTables,polyList,poly_selected,substrate,error,chem_entities);

        //Fill new table
        fill_table(par,allTables, polyList, new_poly_index, substrate, error,chem_entities);



//        moveReactions(allTables[poly_selected],allTables[new_poly_index],polyList[poly_selected], polyList[new_poly_index],poly_selected, bond_selected, nbr_poly);

        nbr_poly++;
    }



    //If the cut was at the end of the strand,this can happen when poly_len=2
    else if(bond_selected != 0 and  bond_selected == polyList[poly_selected].len_poly-1)
    {
        if(polyList[poly_selected].len_poly == 1){
            cout << "CBH action on cellobiose. This should currently be turned off! Stopping..." << endl;
            exit(1);
            return;        
        }        
        else if(polyList[poly_selected].len_poly > 2){
            cout << "CBH action on last bond. This should currently be turned off! Stopping..." << endl;
            exit(1);
            return;
        }
        else{

            taylorOldPoly(polyList[poly_selected],bond_selected);//Erase the end of the polymer, until digested bond
            allTables[poly_selected].num_bond.clear();
            allTables[poly_selected].material.clear();
            allTables[poly_selected].indic_action.clear();
            allTables[poly_selected].liste_prop.clear();
            allTables[poly_selected].prop_uninhib.clear();
            allTables[poly_selected].nbr_element = 0;            

            if(prop(par,3,chem_entities) > 0 and polyList[poly_selected].status[0] == 1 and polyList[poly_selected].len_poly == 1){
                reaction_type = 3;
                propensite = prop(par,reaction_type,chem_entities);                        
                crystal_propensite = par.crystal_modifier_cellu*propensite;
                if(polyList[poly_selected].crystalline[0] == false)
                    addreaction(0, allTables[poly_selected],poly_selected,0,substrate,reaction_type,propensite);
                else
                    addreaction(0, allTables[poly_selected],poly_selected,0,substrate,reaction_type,crystal_propensite);
            }
    
            nbr_Glc_pdt++;//We release a Glc
        }

    }

    //If the cut was at the begining of the strand, this can happen when poly_len=2
    else if(bond_selected==0)
    {

        if(polyList[poly_selected].len_poly == 1){
            cout << "CBH action on cellobiose. This should currently be turned off! Stopping..." << endl;
            exit(1);
            return;        
        }
        else if(polyList[poly_selected].len_poly > 2){
            cout << "CBH action on first bond. This should currently be turned off! Stopping..." << endl;
            exit(1);
            return;                    
        }
        else{
            taylorNewPoly(polyList[poly_selected],bond_selected);//Shift the indices by one to the left
            allTables[poly_selected].num_bond.clear();
            allTables[poly_selected].material.clear();
            allTables[poly_selected].indic_action.clear();
            allTables[poly_selected].liste_prop.clear();
            allTables[poly_selected].prop_uninhib.clear();
            allTables[poly_selected].nbr_element = 0;            

            if(prop(par,3,chem_entities) > 0 and polyList[poly_selected].status[0] == 1 and polyList[poly_selected].len_poly == 1){
                reaction_type = 3;
                propensite = prop(par,reaction_type,chem_entities);                        
                crystal_propensite = par.crystal_modifier_cellu*propensite;
                if(polyList[poly_selected].crystalline[0] == false)
                    addreaction(0, allTables[poly_selected],poly_selected,0,substrate,reaction_type,propensite);
                else
                    addreaction(0, allTables[poly_selected],poly_selected,0,substrate,reaction_type,crystal_propensite);
            }

            nbr_Glc_pdt++;//We release a Glc

        }
    }



    //====================================== Modify propensities; this is only used when the enzyme size is inactive =====================
    if(par.mode_enzyme_size == -1){

        int poly_index = 0;
        int z_index = 0;
        int found_material = 0;
        vector<bList> hemiList,lignList;
        if(bond_neighbors_cellu.count(neighbor_key) != 1){
            cout << "In function CBH_digest: a key is found more than once in bond_neighbors_cellu. This should never happen. Stopping." << endl;
            exit(1);
        }
        else if(bond_neighbors_cellu[neighbor_key].N_neighbors != bond_neighbors_cellu[neighbor_key].x_neighbors.size()){
            cout << "In function CBH_digest: a member of bond_neighbors_cellu has an N_neighbors parameter which differs from the size of the corresponding vector. This should never happen. Stopping." << endl;
            exit(1); 
        }
        for(int i=0; i<bond_neighbors_cellu[neighbor_key].N_neighbors; i++){
            if(find_specific_bond(polyList, hemiList, lignList, bond_neighbors_cellu[neighbor_key].x_neighbors[i], bond_neighbors_cellu[neighbor_key].y_neighbors[i], bond_neighbors_cellu[neighbor_key].z_neighbors[i], poly_index, z_index, found_material) == true){
                if(found_material == 1){
                    if(polyList[poly_index].status[z_index] == -1){
                        polyList[poly_index].status[z_index] = 1;
                        if(polyList[poly_index].len_poly == 1){
                            reaction_type = 3;
                            propensite = prop(par,reaction_type,chem_entities);
                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                            if(propensite > 0){
                                if(polyList[poly_index].crystalline[z_index] == false){
                                    addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                }
                                else{
                                    addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                }
                            }
                        }
                        else if(polyList[poly_index].len_poly == 2){
                            reaction_type = 2;
                            propensite = prop(par,reaction_type,chem_entities);
                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                            if(propensite > 0){
                               if(polyList[poly_index].crystalline[z_index] == false){
                                    addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                }
                                else{
                                    addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                }                            
                            }
                        }
                        else if(polyList[poly_index].len_poly > 2){
                            if(z_index > 1 and z_index < polyList[poly_index].len_poly-2){
                                reaction_type = 1;
                                propensite = prop(par,reaction_type,chem_entities);
                                crystal_propensite = par.crystal_modifier_cellu * propensite;
                                if(propensite > 0){
                                   if(polyList[poly_index].crystalline[z_index] == false){
                                        addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                    }
                                    else{
                                        addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                    }                            
                                }

                            }
                            if(z_index == 1 or z_index == polyList[poly_index].len_poly-2){
                                reaction_type = 2;
                                propensite = prop(par,reaction_type,chem_entities);
                                crystal_propensite = par.crystal_modifier_cellu * propensite;
                                if(propensite > 0){
                                   if(polyList[poly_index].crystalline[z_index] == false){
                                        addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                    }
                                    else{
                                        addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                    }                            
                                }
                            }

                        }
                    }

                }
            }
        }
    }



/*
    //Make available the surrounding bonds for digestion
    for(int i=0; i<nbr_poly; i++)//We only scan cellulose because lignin is anyway undigestable and hemicellulose availability is fixed at the begining
    {
        if(polyList[i].len_poly>0)//THIS WAS 0; TESTING!!!
        {
            if(check_position(polyList[i].x, polyList[i].y, polyList[poly_selected].x, polyList[poly_selected].y)==true)//((polyList[i].x==polyList[poly_selected].x) or (polyList[i].x==polyList[poly_selected].x-1) or (polyList[i].x==polyList[poly_selected].x+1)) and ((polyList[i].y==polyList[poly_selected].y) or (polyList[i].y==polyList[poly_selected].y-1) or (polyList[i].y==polyList[poly_selected].y+1)))
            {
                for(int j=0; j<polyList[i].len_poly; j++)
                {
                    if((polyList[i].z[j]==cut_applicate) or (polyList[i].z[j]==cut_applicate-1) or (polyList[i].z[j]==cut_applicate+1) or (polyList[i].z[j]==cut_applicate-2) or (polyList[i].z[j]==cut_applicate+2 ))
                    {
                        polyList[i].status[j]=1;
                        reaction_type = 1;
                        propensite = prop(par,reaction_type,chem_entities);                        
                        crystal_propensite = par.crystal_modifier_cellu*propensite;
                        if(propensite > 0 and j != 0 and j != polyList[i].len_poly-1){
                            if(polyList[i].crystalline[j] == false)
                                addreaction(0, allTables[i],i,j,substrate,reaction_type,propensite);
                            else
                                addreaction(0, allTables[i],i,j,substrate,reaction_type,crystal_propensite);
                        }
                        reaction_type = 2;
                        propensite = prop(par,reaction_type,chem_entities);                        
                        crystal_propensite = par.crystal_modifier_cellu*propensite;
                        if(propensite > 0 and polyList[i].len_poly >= 2){
                            if(polyList[i].len_poly == 2){
                                if(j == 0 or j==1){
                                    if(polyList[i].crystalline[j] == false)
                                        addreaction(0, allTables[i],i,j,substrate,reaction_type,propensite);
                                    else
                                        addreaction(0, allTables[i],i,j,substrate,reaction_type,crystal_propensite);
                                }
                            }
                            else if(polyList[i].len_poly > 2){
                                if(j == 1 or j == polyList[i].len_poly-2){
                                    if(polyList[i].crystalline[j] == false)
                                        addreaction(0, allTables[i],i,j,substrate,reaction_type,propensite);
                                    else
                                        addreaction(0, allTables[i],i,j,substrate,reaction_type,crystal_propensite);
                                }
                            }
                        }
                        reaction_type = 3;
                        propensite = prop(par,reaction_type,chem_entities);                        
                        crystal_propensite = par.crystal_modifier_cellu*propensite;
                        if(propensite > 0 and polyList[i].len_poly == 1){
                            if(polyList[i].crystalline[0] == false)
                                addreaction(0, allTables[i],i,j,substrate,reaction_type,propensite);
                            else
                                addreaction(0, allTables[i],i,j,substrate,reaction_type,crystal_propensite);
                        }
                    }
                }
            }
        }
    }
*/

    if(par.mode_enzyme_size == -1){
        remove_neighbor_from_vector(x_pos,y_pos,z_pos,bond_neighbors_cellu,bond_neighbors_hemi);
    }

    if(polyList[poly_selected].len_poly == 1){
        nbr_cellobiose++;
    }
    if(newPolyFlag == true){
        if(polyList[new_poly_index].len_poly == 1){
            nbr_cellobiose++;
        }
    }
//End of function    
}



void BGL_digest(params& par, vector<TList>& allTables, vector<bList>& polyList, int& nbr_poly, int& nbr_xyl_pdt, int& nbr_Glc_pdt,const int bond_selected,const int poly_selected, int len_polyLoopStart, int& nbr_cellobiose,vector<float>& chem_entities, const int substrate,unordered_map<int,neighborList>& bond_neighbors_cellu, unordered_map<int,neighborList>& bond_neighbors_hemi, const bool verbose, bool& error){

    if(polyList[poly_selected].len_poly > 1){
        cout << "BGL acts on something other than cellobiose! Stopping..." << endl;
        exit(1);
        return;
    }
    if(polyList[poly_selected].len_poly <= bond_selected){
        cout << "In function BGL_digest: the bond of the current reaction is larger than the size of the polymer!" << endl;
        exit(1);
        return;
    }
    int check_radius = 2;
    int dx = 1;//Distance of bonds in coordinates
    int dy = 1;
    int dz = 1;

    int x_pos = polyList[poly_selected].x;//Exact coordinates of the bond
    int y_pos = polyList[poly_selected].y;
    int z_pos = polyList[poly_selected].z[bond_selected];

    int neighbor_key = cantor_pair_three(x_pos, y_pos, z_pos);

    if(par.enzyme_radius > 0){
        check_radius = 2*par.enzyme_radius;//Number of bonds which are checked around each position
    }
    int cut_applicate=polyList[poly_selected].z[bond_selected];
    int test1;
    int indic_cut;
    double propensite;
    double crystal_propensite;
    int reaction_type = 0;
    if(bond_selected==0)//If the cut takes place at the head of the polymer
        indic_cut=0;
    if(bond_selected==polyList[poly_selected].len_poly-1)//If the cut takes place at the tip of the polymer
        indic_cut=1;
    //Remove all reactions associated to the bond just digested
    if(allTables[poly_selected].nbr_element > 1){
        cout << "In function BGL_digest: cellobiose contains more than the one BGL reaction; Reactions: " << endl;
        for(int i=0;i<allTables[poly_selected].nbr_element;i++)
            cout << "Enzyme: " << allTables[poly_selected].indic_action[i] << "; bond: " << allTables[poly_selected].num_bond[i] << "; material: " << allTables[poly_selected].material[i] << endl;
        cout << "Stopping..." << endl;
        exit(1);
        return;
    }
    deleteAllreaction(allTables[poly_selected],allTables,substrate,poly_selected,bond_selected);
    if(polyList[poly_selected].len_poly>1)
    {
        cout << "BGL is digesting something larger than cellobiose!! this should not currently be happening. stopping..." << endl;
        exit(1);
        return;
    }
    else if(polyList[poly_selected].len_poly==1)
    {
        polyList[poly_selected].status[bond_selected]=0;
        polyList[poly_selected].len_poly=0;
        nbr_Glc_pdt=nbr_Glc_pdt+2;
        allTables[poly_selected].prop_sum = 0;
        nbr_cellobiose--;
    }


    //====================================== Modify propensities; this is only used when the enzyme size is inactive =====================
    if(par.mode_enzyme_size == -1){
        int poly_index = 0;
        int z_index = 0;
        int found_material = 0;
        vector<bList> hemiList,lignList;
        if(bond_neighbors_cellu.count(neighbor_key) != 1){
            cout << "In function BGL_digest: a key is found more than once in bond_neighbors_cellu. This should never happen. Stopping." << endl;
            exit(1);
        }
        else if(bond_neighbors_cellu[neighbor_key].N_neighbors != bond_neighbors_cellu[neighbor_key].x_neighbors.size()){
            cout << "In function BGL_digest: a member of bond_neighbors_cellu has an N_neighbors parameter which differs from the size of the corresponding vector. This should never happen. Stopping." << endl;
            exit(1);
        }
        for(int i=0; i<bond_neighbors_cellu[neighbor_key].N_neighbors; i++){
            if(find_specific_bond(polyList, hemiList, lignList, bond_neighbors_cellu[neighbor_key].x_neighbors[i], bond_neighbors_cellu[neighbor_key].y_neighbors[i], bond_neighbors_cellu[neighbor_key].z_neighbors[i], poly_index, z_index, found_material) == true){
                if(found_material == 1){
                    if(polyList[poly_index].status[z_index] == -1){
                        polyList[poly_index].status[z_index] = 1;
                        if(polyList[poly_index].len_poly == 1){
                            reaction_type = 3;
                            propensite = prop(par,reaction_type,chem_entities);
                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                            if(propensite > 0){
                                if(polyList[poly_index].crystalline[z_index] == false){
                                    addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                }
                                else{
                                    addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                }
                            }
                        }
                        else if(polyList[poly_index].len_poly == 2){
                            reaction_type = 2;
                            propensite = prop(par,reaction_type,chem_entities);
                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                            if(propensite > 0){
                               if(polyList[poly_index].crystalline[z_index] == false){
                                    addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                }
                                else{
                                    addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                }                            
                            }
                        }
                        else if(polyList[poly_index].len_poly > 2){
                            if(z_index > 1 and z_index < polyList[poly_index].len_poly-2){
                                reaction_type = 1;
                                propensite = prop(par,reaction_type,chem_entities);
                                crystal_propensite = par.crystal_modifier_cellu * propensite;
                                if(propensite > 0){
                                   if(polyList[poly_index].crystalline[z_index] == false){
                                        addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                    }
                                    else{
                                        addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                    }                            
                                }

                            }
                            if(z_index == 1 or z_index == polyList[poly_index].len_poly-2){
                                reaction_type = 2;
                                propensite = prop(par,reaction_type,chem_entities);
                                crystal_propensite = par.crystal_modifier_cellu * propensite;
                                if(propensite > 0){
                                   if(polyList[poly_index].crystalline[z_index] == false){
                                        addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                    }
                                    else{
                                        addreaction(true, allTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                    }                            
                                }
                            }

                        }
                    }

                }

            }
        }
    }



//Was using the one above
/*

    for(int i=0; i<nbr_poly; i++)//We only scan cellulose because lignin is anyway undigestable and hemicellulose availability is fixed at the begining because it's made of one external layer
    {
        if(polyList[i].len_poly>0)//THIS WAS 0; TESTING!!!
        {
            if(check_position(polyList[i].x, polyList[i].y, polyList[poly_selected].x, polyList[poly_selected].y)==true)//((polyList[i].x==polyList[poly_selected].x) or (polyList[i].x==polyList[poly_selected].x-1) or (polyList[i].x==polyList[poly_selected].x+1)) and ((polyList[i].y==polyList[poly_selected].y) or (polyList[i].y==polyList[poly_selected].y-1) or (polyList[i].y==polyList[poly_selected].y+1)))
            {
                for(int j=0; j<polyList[i].len_poly; j++)
                {
                    if((polyList[i].z[j]==cut_applicate) or (polyList[i].z[j]==cut_applicate-1) or (polyList[i].z[j]==cut_applicate+1))  //or ((bond_selected == 1 or bond_selected == polyList[poly_selected].len_poly + polyList[new_poly_index].len_poly -1) and ((polyList[i].z[j]==cut_applicate-2) or (polyList[i].z[j]==cut_applicate+2))) )
                    {
                        polyList[i].status[j]=1;
                        reaction_type = 1;
                        propensite = prop(par,reaction_type,chem_entities);                        
                        crystal_propensite = par.crystal_modifier_cellu*propensite;
                        if(propensite > 0 and j != 0 and j != polyList[i].len_poly-1){
                            if(polyList[i].crystalline[j] == false)
                                addreaction(0, allTables[i],i,j,substrate,reaction_type,propensite);
                            else
                                addreaction(0, allTables[i],i,j,substrate,reaction_type,crystal_propensite);
                        }
                        reaction_type = 2;
                        propensite = prop(par,reaction_type,chem_entities);                        
                        crystal_propensite = par.crystal_modifier_cellu*propensite;
                        if(propensite > 0 and polyList[i].len_poly >= 2){
                            if(polyList[i].len_poly == 2){
                                if(j == 0 or j==1){
                                    if(polyList[i].crystalline[j] == false)
                                        addreaction(0, allTables[i],i,j,substrate,reaction_type,propensite);
                                    else
                                        addreaction(0, allTables[i],i,j,substrate,reaction_type,crystal_propensite);
                                }
                            }
                            else if(polyList[i].len_poly > 2){
                                if(j == 1 or j == polyList[i].len_poly-2){
                                    if(polyList[i].crystalline[j] == false)
                                        addreaction(0, allTables[i],i,j,substrate,reaction_type,propensite);
                                    else
                                        addreaction(0, allTables[i],i,j,substrate,reaction_type,crystal_propensite);
                                }
                            }
                        }
                        reaction_type = 3;
                        propensite = prop(par,reaction_type,chem_entities);                        
                        crystal_propensite = par.crystal_modifier_cellu*propensite;
                        if(propensite > 0 and polyList[i].len_poly == 1){
                            if(polyList[i].crystalline[0] == false)
                                addreaction(0, allTables[i],i,j,substrate,reaction_type,propensite);
                            else
                                addreaction(0, allTables[i],i,j,substrate,reaction_type,crystal_propensite);
                        }
                    }
                }
            }
        }
    }*/

    if(par.mode_enzyme_size == -1){
        remove_neighbor_from_vector(x_pos,y_pos,z_pos,bond_neighbors_cellu,bond_neighbors_hemi);
    }

//End of function    
}

//Digestion by Xylanase
void XYL_digest(params& par, vector<TList>& allTables, vector<TList>& allCelluTables, vector<bList>& polyList, vector<bList>& celluList, int& nbr_poly, int& nbr_xyl_pdt, int& nbr_Glc_pdt,const int bond_selected,const int poly_selected, int len_polyLoopStart, int& nbr_cellobiose,vector<float>& chem_entities, const int substrate,unordered_map<int,neighborList>& bond_neighbors_cellu, unordered_map<int,neighborList>& bond_neighbors_hemi, const bool verbose, bool& error){
    if(polyList[poly_selected].len_poly <= bond_selected){
        cout << "In function XYL_digest: the bond of the current reaction is larger than the size of the polymer! bond: " << bond_selected << "; len_poly = " << polyList[poly_selected].len_poly << endl;
        cout << " Stopping" << endl;
        exit(1);
        return;
    }


    if(polyList[poly_selected].len_poly <= bond_selected){
        cout << "In function XYL_digest: the bond of the current reaction is larger than the size of the polymer!" << endl;
        exit(1);
        return;
    }


    int check_radius = 2;
    int dx = 1;//Distance of bonds in coordinates
    int dy = 1;
    int dz = 1;

    int x_pos = polyList[poly_selected].x;//Exact coordinates of the bond
    int y_pos = polyList[poly_selected].y;
    int z_pos = polyList[poly_selected].z[bond_selected];

    if(par.enzyme_radius > 0){
        check_radius = 2*par.enzyme_radius;//Number of bonds which are checked around each position
    }    
    int len_poly_before = polyList[poly_selected].len_poly;
    int cut_applicate = polyList[poly_selected].z[bond_selected];
    int test1; //Used for finding specific reactions in a table
    int indic_cut;//
    int new_poly_index = -1;
    int minSize_BGL = 8;//DP at which BGL starts digesting 


    double propensite; //
    double crystal_propensite;
    int reaction_type = 0;


    int neighbor_key = cantor_pair_three(x_pos, y_pos, z_pos);

    bool newPolyFlag = false; //Checks whether a new polymer is generated here

    //Remove all reactions associated to the bond just digested
    deleteAllreaction(allTables[poly_selected],allTables,substrate,poly_selected,bond_selected);


    //======================= Remove the bond and split the polymer into two new polymers =============================

    if(polyList[poly_selected].len_poly==0){
        //Do nothing
    }
    else if(polyList[poly_selected].len_poly == 1){
        polyList[poly_selected].status[bond_selected]=0;
        polyList[poly_selected].len_poly=0;
        nbr_xyl_pdt += +2;
        allTables[poly_selected].prop_sum = 0;
    }
    else if(polyList[poly_selected].len_poly==2)
    {

        if(bond_selected==0)// If the cut was at the beginning of the polymer
        {                        
            taylorNewPoly(polyList[poly_selected],bond_selected);//Shift the indices by one to the left
            allTables[poly_selected].num_bond.clear();
            allTables[poly_selected].material.clear();
            allTables[poly_selected].indic_action.clear();
            allTables[poly_selected].liste_prop.clear();
            allTables[poly_selected].prop_uninhib.clear();
            allTables[poly_selected].nbr_element = 0;            

            reaction_type = 4;
            propensite = prop(par,reaction_type,chem_entities);                        
            crystal_propensite = par.crystal_modifier_hemi*propensite;
            if(propensite > 0 and polyList[poly_selected].status[0] == 1 and polyList[poly_selected].len_poly == 1 and polyList[poly_selected].bond_type[0] == 4){
                if(polyList[poly_selected].crystalline[0] == false)
                    addreaction(1, allTables[poly_selected],poly_selected,0,substrate,reaction_type,propensite);
                else
                    addreaction(1, allTables[poly_selected],poly_selected,0,substrate,reaction_type,crystal_propensite);
            }

            nbr_xyl_pdt++;//We release a Xyl
        }
        else if(bond_selected == 1){
            taylorOldPoly(polyList[poly_selected],bond_selected);//Erase the end of the polymer, until digested bond
            allTables[poly_selected].num_bond.clear();
            allTables[poly_selected].material.clear();
            allTables[poly_selected].indic_action.clear();
            allTables[poly_selected].liste_prop.clear();
            allTables[poly_selected].prop_uninhib.clear();
            allTables[poly_selected].nbr_element = 0;            

            reaction_type = 4;
            propensite = prop(par,reaction_type,chem_entities);                        
            crystal_propensite = par.crystal_modifier_hemi*propensite;
            if(propensite > 0 and polyList[poly_selected].status[0] == 1 and polyList[poly_selected].len_poly == 1 and polyList[poly_selected].bond_type[0] == 4){
                if(polyList[poly_selected].crystalline[0] == false)
                    addreaction(1, allTables[poly_selected],poly_selected,0,substrate,reaction_type,propensite);
                else
                    addreaction(1, allTables[poly_selected],poly_selected,0,substrate,reaction_type,crystal_propensite);
            }

            nbr_xyl_pdt++;//We release a Xyl           
        }
    }
    else if(polyList[poly_selected].len_poly>2){
        if((bond_selected !=0) and (bond_selected !=polyList[poly_selected].len_poly-1))
        {

        //====================== CREATE AND INITIALIZE THE NEW POLYMER AND ITS REACTION TABLE =========================0

            polyList.push_back(bList());
            new_poly_index = polyList.size()-1;
            newPolyFlag = true;


            polyList[new_poly_index].index= new_poly_index;//Specify the new index of the new polymer
            polyList[new_poly_index].x= polyList[poly_selected].x;
            polyList[new_poly_index].y= polyList[poly_selected].y;
            polyList[new_poly_index].set_z(bond_selected+1);

            allTables.push_back(TList());//Make the new reaction list associated to the new polymer;
            initTList(allTables[new_poly_index], new_poly_index);//Initialize it                        

            for(int i=0; i<polyList[poly_selected].len_poly; i++)//Make the new polymer into a copy of the old one before tayloring
            {
                addbond(polyList[new_poly_index],i,1,polyList[poly_selected].bond_type[i],false);
                polyList[new_poly_index].z[i]=polyList[poly_selected].z[i];
                polyList[new_poly_index].N_blocked_positions[i] = polyList[poly_selected].N_blocked_positions[i];
                polyList[new_poly_index].crystalline[i] = polyList[poly_selected].crystalline[i];                

            }

            taylorOldPoly(polyList[poly_selected],bond_selected);//Erase the end of the old polymer, until digested bond
            taylorNewPoly(polyList[new_poly_index],bond_selected);//Erase the begining of the new polymer, from digested bond



            //============================== Transfer the indices of polymer and bond position from the Old to the New polymer in the Reaction Table
            //Adjust old table
            allTables[poly_selected].num_bond.clear();
            allTables[poly_selected].material.clear();
            allTables[poly_selected].indic_action.clear();
            allTables[poly_selected].liste_prop.clear();
            allTables[poly_selected].prop_uninhib.clear();
            allTables[poly_selected].nbr_element = 0;
            fill_table(par,allTables,polyList,poly_selected,substrate,error,chem_entities);

            //Fill new table
            fill_table(par,allTables, polyList, new_poly_index, substrate, error,chem_entities);


            nbr_poly++;
        }
        else if(bond_selected == 0){
            taylorNewPoly(polyList[poly_selected],bond_selected);//Erase the begining of the new polymer, from digested bond            
            //Adjust old table
            allTables[poly_selected].num_bond.clear();
            allTables[poly_selected].material.clear();
            allTables[poly_selected].indic_action.clear();
            allTables[poly_selected].liste_prop.clear();
            allTables[poly_selected].prop_uninhib.clear();
            allTables[poly_selected].nbr_element = 0;
            fill_table(par,allTables,polyList,poly_selected,substrate,error,chem_entities);            
            nbr_xyl_pdt++;
        }
        else if(bond_selected == polyList[poly_selected].len_poly-1){
            taylorOldPoly(polyList[poly_selected],bond_selected);//Erase the end of the polymer, until digested bond
            //Adjust old table
            allTables[poly_selected].num_bond.clear();
            allTables[poly_selected].material.clear();
            allTables[poly_selected].indic_action.clear();
            allTables[poly_selected].liste_prop.clear();
            allTables[poly_selected].prop_uninhib.clear();
            allTables[poly_selected].nbr_element = 0;
            fill_table(par,allTables,polyList,poly_selected,substrate,error,chem_entities);            
            nbr_xyl_pdt++;
        }

    }


    int cellu_substrate = 1;

    //====================================== Modify propensities; this is only used when the enzyme size is inactive =====================
    if(par.mode_enzyme_size == -1){
        int poly_index = 0;
        int z_index = 0;
        int found_material = 0;
        for(int i=0; i<bond_neighbors_hemi[neighbor_key].N_neighbors; i++){
            if(find_specific_bond(celluList, polyList, polyList, bond_neighbors_hemi[neighbor_key].x_neighbors[i], bond_neighbors_hemi[neighbor_key].y_neighbors[i], bond_neighbors_hemi[neighbor_key].z_neighbors[i], poly_index, z_index, found_material) == true){
                if(found_material == 1){
    //                cout << "Found cellu in hemi" << endl;
                    if(celluList[poly_index].status[z_index] == -1){
                        celluList[poly_index].status[z_index] = 1;
                        if(celluList[poly_index].len_poly == 1){
                            reaction_type = 3;
                            propensite = prop(par,reaction_type,chem_entities);
                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                            if(propensite > 0){
                                if(celluList[poly_index].crystalline[z_index] == false){
                                    addreaction(true, allCelluTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                }
                                else{
                                    addreaction(true, allCelluTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                }
                            }
                        }
                        else if(celluList[poly_index].len_poly == 2){
                            reaction_type = 2;
                            propensite = prop(par,reaction_type,chem_entities);
                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                            if(propensite > 0){
                               if(celluList[poly_index].crystalline[z_index] == false){
                                    addreaction(true, allCelluTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                }
                                else{
                                    addreaction(true, allCelluTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                }                            
                            }
                        }
                        else if(celluList[poly_index].len_poly > 2){
                            if(z_index > 1 and z_index < celluList[poly_index].len_poly-2){
                                reaction_type = 1;
                                propensite = prop(par,reaction_type,chem_entities);
                                crystal_propensite = par.crystal_modifier_cellu * propensite;
                                if(propensite > 0){
                                   if(celluList[poly_index].crystalline[z_index] == false){
                                        addreaction(true, allCelluTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                    }
                                    else{
                                        addreaction(true, allCelluTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                    }                            
                                }

                            }
                            if(z_index == 1 or z_index == celluList[poly_index].len_poly-2){
                                reaction_type = 2;
                                propensite = prop(par,reaction_type,chem_entities);
                                crystal_propensite = par.crystal_modifier_cellu * propensite;
                                if(propensite > 0){
                                   if(celluList[poly_index].crystalline[z_index] == false){
                                        addreaction(true, allCelluTables[poly_index],poly_index, z_index,1, reaction_type,propensite);
                                    }
                                    else{
                                        addreaction(true, allCelluTables[poly_index],poly_index, z_index,1, reaction_type,crystal_propensite);
                                    }                            
                                }
                            }

                        }
                    }

                }
                else if(found_material == 2){
                    if(polyList[poly_index].status[z_index] == -1){
                        polyList[poly_index].status[z_index] = 1;
                        reaction_type = 4;
                        propensite = prop(par,reaction_type,chem_entities);
                        crystal_propensite = par.crystal_modifier_hemi * propensite;
                        if(propensite > 0 and z_index < polyList[poly_index].len_poly){
                           if(polyList[poly_index].crystalline[z_index] == false){
                                addreaction(true, allTables[poly_index],poly_index, z_index,2, reaction_type,propensite);
                            }
                            else{
                                addreaction(true, allTables[poly_index],poly_index, z_index,2, reaction_type,crystal_propensite);
                            }                            
                        }
                    }
                }
            }

        }
    }







    if(polyList[poly_selected].len_poly == 1){
//        nbr_xylobiose++; //xylobiose is not included yet

    }
    if(newPolyFlag == true){
        if(polyList[new_poly_index].len_poly == 1){
//            nbr_xylobiose++; //xylobiose is not included yet

        }
    }
    if(newPolyFlag == true){
        if(polyList[poly_selected].len_poly+1 + polyList[new_poly_index].len_poly+1 != len_poly_before+1){
            cout << "In function XYL_digest(): Xylose number is not conserved. Stopping" << endl;
            exit(1);
        }
    }
    else if(polyList[poly_selected].len_poly+1 != len_poly_before){
            cout << "In function XYL_digest(): Xylose number is not conserved. Stopping" << endl;
            exit(1);        
    }

    if(verbose == true){
        for(int i=0; i<allTables.size();i++){
            for(int j=0;j<allTables[i].nbr_element;j++){
                if(polyList[i].len_poly <= allTables[i].num_bond[j]){
                    cout << "In function XYL_digest(), specifically at the end: There are some bonds here which are at higher positions than the length of the polymer." << endl;
                    cout << "Poly: " << i << "; len_poly: " << polyList[i].len_poly << "; num_bond: " << allTables[i].num_bond[j] << endl;
                    exit(1);
                    return;
                }
            }
        }
    }

    if(par.mode_enzyme_size == -1){
        remove_neighbor_from_vector(x_pos,y_pos,z_pos,bond_neighbors_cellu,bond_neighbors_hemi);
    }

//End of function    
}


//Updates reaction tables after a digestion reaction. Called at each step that is not a lignin adhesion step
void update_reactiontables(std::vector<bList>& cellu, std::vector<bList>& hemi, std::vector<bList>& lign, std::vector<TList>& Table_cellu, std::vector<TList>& Table_hemi, params& par, std::vector<float>& chem_entities, std::unordered_map<int, neighborList>& bond_neighbors_cellu, std::unordered_map<int, neighborList>& bond_neighbors_hemi, int x, int y, int z, int substrate, int poly_selected, int bond_selected){
    int bond_key = cantor_pair_three(x,y,z);//Used to find the neighborList object associated to the previously digested bond at coordinates x,y,z
    int neighbor_key = 0;
    int current_poly = 0;
    int current_bond = 0;
    int current_substrate = 0;
    int reaction_type = 0;
    float propensite = 0;
    float crystal_propensite = 0;
    if(substrate == 1){
        if(bond_neighbors_cellu.count(bond_key) != 0){
            for(int i=0;i<bond_neighbors_cellu[bond_key].N_neighbors; i++){
//                if(par.verbose == true)
 //                   cout << "neighbor coordinates: " << bond_neighbors_cellu[bond_key].x_neighbors[i] << ", " << bond_neighbors_cellu[bond_key].y_neighbors[i] << ", " << bond_neighbors_cellu[bond_key].z_neighbors[i] << endl;
                neighbor_key = cantor_pair_three(bond_neighbors_cellu[bond_key].x_neighbors[i],bond_neighbors_cellu[bond_key].y_neighbors[i], bond_neighbors_cellu[bond_key].z_neighbors[i]);
                if(bond_neighbors_cellu.count(neighbor_key) != 0){
//                    cout << "Found neighbor" << endl;
                    bond_neighbors_cellu[neighbor_key].remove_neighbor(x,y,z);
                    if(bond_neighbors_cellu[neighbor_key].outer_bond == true){
                        if(find_specific_bond(cellu, hemi, lign, bond_neighbors_cellu[neighbor_key].x, bond_neighbors_cellu[neighbor_key].y, bond_neighbors_cellu[neighbor_key].z, current_poly, current_bond, current_substrate) == true){
                            if(current_substrate == 1){
                                if(cellu[current_poly].status[current_bond] == -1){
                                    cellu[current_poly].status[current_bond] = 1;
                                    if(cellu[current_poly].len_poly == 1){
                                        reaction_type = 3;
                                        propensite = prop(par,reaction_type,chem_entities);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;
                                        if(propensite > 0){
                                            if(cellu[current_poly].crystalline[current_bond] == false){
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                            }
                                            else{
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                            }
                                        }
                                    }
                                    else if(cellu[current_poly].len_poly == 2){
                                        reaction_type = 2;
                                        propensite = prop(par,reaction_type,chem_entities);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;
                                        if(propensite > 0){
                                           if(cellu[current_poly].crystalline[current_bond] == false){
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                            }
                                            else{
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                            }                            
                                        }
                                    }
                                    else if(cellu[current_poly].len_poly > 2){
                                        if(current_bond > 1 and current_bond < cellu[current_poly].len_poly-2){
                                            reaction_type = 1;
                                            propensite = prop(par,reaction_type,chem_entities);
                                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                                            if(propensite > 0){
                                               if(cellu[current_poly].crystalline[current_bond] == false){
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                                }
                                                else{
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                                }                            
                                            }

                                        }
                                        if(current_bond == 1 or current_bond == cellu[current_poly].len_poly-2){
                                            reaction_type = 2;
                                            propensite = prop(par,reaction_type,chem_entities);
                                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                                            if(propensite > 0){
                                               if(cellu[current_poly].crystalline[current_bond] == false){
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                                }
                                                else{
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                                }                            
                                            }
                                        }

                                    }

                                }
                            }
                            else if(current_substrate == 2){
                                if(hemi[current_poly].status[current_bond] == -1){
                                    hemi[current_poly].status[current_bond] = 1;
//                                    hemi[current_poly].status[current_bond] = 1;
                                    reaction_type = 4;
                                    propensite = prop(par,reaction_type,chem_entities);
                                    crystal_propensite = par.crystal_modifier_hemi * propensite;

                                    if(propensite > 0 and current_bond < hemi[current_poly].len_poly){
                                       if(hemi[current_poly].crystalline[current_bond] == false){
                                            addreaction(true, Table_hemi[current_poly],current_poly, current_bond,2, reaction_type,propensite);
                                        }
                                        else{
                                            addreaction(true, Table_hemi[current_poly],current_poly, current_bond,2, reaction_type,crystal_propensite);
                                        }                            
                                    }
                                }
                            }
                        }
                    }
                }
                else if(bond_neighbors_hemi.count(neighbor_key) != 0){
                    bond_neighbors_hemi[neighbor_key].remove_neighbor(x,y,z);
                    if(bond_neighbors_hemi[neighbor_key].outer_bond == true){
                        if(find_specific_bond(cellu, hemi, lign, bond_neighbors_hemi[neighbor_key].x, bond_neighbors_hemi[neighbor_key].y, bond_neighbors_hemi[neighbor_key].z, current_poly, current_bond, current_substrate) == true){
         //                   cout << "found cellu" << endl;
                            if(current_substrate == 1){
                                if(cellu[current_poly].status[current_bond] == -1){
                                    cellu[current_poly].status[current_bond] = 1;
                                    if(cellu[current_poly].len_poly == 1){
                                        reaction_type = 3;
                                        propensite = prop(par,reaction_type,chem_entities);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;
                                        if(propensite > 0){
                                            if(cellu[current_poly].crystalline[current_bond] == false){
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                            }
                                            else{
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                            }
                                        }
                                    }
                                    else if(cellu[current_poly].len_poly == 2){
                                        reaction_type = 2;
                                        propensite = prop(par,reaction_type,chem_entities);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;
                                        if(propensite > 0){
                                           if(cellu[current_poly].crystalline[current_bond] == false){
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                            }
                                            else{
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                            }                            
                                        }
                                    }
                                    else if(cellu[current_poly].len_poly > 2){
                                        if(current_bond > 1 and current_bond < cellu[current_poly].len_poly-2){
                                            reaction_type = 1;
                                            propensite = prop(par,reaction_type,chem_entities);
                                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                                            if(propensite > 0){
                                               if(cellu[current_poly].crystalline[current_bond] == false){
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                                }
                                                else{
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                                }                            
                                            }

                                        }
                                        if(current_bond == 1 or current_bond == cellu[current_poly].len_poly-2){
                                            reaction_type = 2;
                                            propensite = prop(par,reaction_type,chem_entities);
                                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                                            if(propensite > 0){
                                               if(cellu[current_poly].crystalline[current_bond] == false){
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                                }
                                                else{
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                                }                            
                                            }
                                        }

                                    }

                                }
                            }
                            else if(current_substrate == 2){
                                if(hemi[current_poly].status[current_bond] == -1){
                                    hemi[current_poly].status[current_bond] = 1;
                                    hemi[current_poly].status[current_bond] = 1;
                                    reaction_type = 4;
                                    propensite = prop(par,reaction_type,chem_entities);
                                    crystal_propensite = par.crystal_modifier_hemi * propensite;
                                    
                                    if(propensite > 0 and current_bond < hemi[current_poly].len_poly){
                                       if(hemi[current_poly].crystalline[current_bond] == false){
                                            addreaction(true, Table_hemi[current_poly],current_poly, current_bond,2, reaction_type,propensite);
                                        }
                                        else{
                                            addreaction(true, Table_hemi[current_poly],current_poly, current_bond,2, reaction_type,crystal_propensite);
                                        }                            
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        bond_neighbors_cellu.erase(bond_key);
    }
    else if(substrate == 2){
        if(bond_neighbors_hemi.count(bond_key) != 0){
            for(int i=0;i<bond_neighbors_hemi[bond_key].N_neighbors; i++){
                if(par.verbose == true)
                    cout << "neighbor coordinates: " << bond_neighbors_hemi[bond_key].x_neighbors[i] << ", " << bond_neighbors_hemi[bond_key].y_neighbors[i] << ", " << bond_neighbors_hemi[bond_key].z_neighbors[i] << endl;
                neighbor_key = cantor_pair_three(bond_neighbors_hemi[bond_key].x_neighbors[i],bond_neighbors_hemi[bond_key].y_neighbors[i], bond_neighbors_hemi[bond_key].z_neighbors[i]);
                if(bond_neighbors_cellu.count(neighbor_key) != 0){
    //                cout << "Found neighbor" << endl;
                    if(bond_neighbors_cellu[neighbor_key].outer_bond == false){
    //                    cout << "No outer bond" << endl;
                    }
                    bond_neighbors_cellu[neighbor_key].remove_neighbor(x,y,z);
                    if(bond_neighbors_cellu[neighbor_key].outer_bond == true){
       //                 cout << "outer bond!" << endl;
                        if(find_specific_bond(cellu, hemi, lign, bond_neighbors_cellu[neighbor_key].x, bond_neighbors_cellu[neighbor_key].y, bond_neighbors_cellu[neighbor_key].z, current_poly, current_bond, current_substrate) == true){
         //                   cout << "found hemi" << endl;
                            if(current_substrate == 1){
                                if(cellu[current_poly].status[current_bond] == -1){
                                    cellu[current_poly].status[current_bond] = 1;
                                    if(cellu[current_poly].len_poly == 1){
                                        reaction_type = 3;
                                        propensite = prop(par,reaction_type,chem_entities);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;
                                        if(propensite > 0){
                                            if(cellu[current_poly].crystalline[current_bond] == false){
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                            }
                                            else{
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                            }
                                        }
                                    }
                                    else if(cellu[current_poly].len_poly == 2){
                                        reaction_type = 2;
                                        propensite = prop(par,reaction_type,chem_entities);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;
                                        if(propensite > 0){
                                           if(cellu[current_poly].crystalline[current_bond] == false){
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                            }
                                            else{
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                            }                            
                                        }
                                    }
                                    else if(cellu[current_poly].len_poly > 2){
                                        if(current_bond > 1 and current_bond < cellu[current_poly].len_poly-2){
                                            reaction_type = 1;
                                            propensite = prop(par,reaction_type,chem_entities);
                                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                                            if(propensite > 0){
                                               if(cellu[current_poly].crystalline[current_bond] == false){
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                                }
                                                else{
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                                }                            
                                            }

                                        }
                                        if(current_bond == 1 or current_bond == cellu[current_poly].len_poly-2){
                                            reaction_type = 2;
                                            propensite = prop(par,reaction_type,chem_entities);
                                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                                            if(propensite > 0){
                                               if(cellu[current_poly].crystalline[current_bond] == false){
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                                }
                                                else{
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                                }                            
                                            }
                                        }

                                    }

                                }
                            }
                            else if(current_substrate == 2){
                                if(hemi[current_poly].status[current_bond] == -1){
                                    hemi[current_poly].status[current_bond] = 1;
                                    hemi[current_poly].status[current_bond] = 1;
                                    reaction_type = 4;
                                    propensite = prop(par,reaction_type,chem_entities);
                                    crystal_propensite = par.crystal_modifier_hemi * propensite;

                                    if(propensite > 0 and current_bond < hemi[current_poly].len_poly){
                                       if(hemi[current_poly].crystalline[current_bond] == false){
                                            addreaction(true, Table_hemi[current_poly],current_poly, current_bond,2, reaction_type,propensite);
                                        }
                                        else{
                                            addreaction(true, Table_hemi[current_poly],current_poly, current_bond,2, reaction_type,crystal_propensite);
                                        }                            
                                    }
                                }
                            }
                        }
                    }
                }
                else if(bond_neighbors_hemi.count(neighbor_key) != 0){
                    bond_neighbors_hemi[neighbor_key].remove_neighbor(x,y,z);
                    if(bond_neighbors_hemi[neighbor_key].outer_bond == true){
                        if(find_specific_bond(cellu, hemi, lign, bond_neighbors_hemi[neighbor_key].x, bond_neighbors_hemi[neighbor_key].y, bond_neighbors_hemi[neighbor_key].z, current_poly, current_bond, current_substrate) == true){
                            if(current_substrate == 1){
                                if(cellu[current_poly].status[current_bond] == -1){
                                    cellu[current_poly].status[current_bond] = 1;
                                    if(cellu[current_poly].len_poly == 1){
                                        reaction_type = 3;
                                        propensite = prop(par,reaction_type,chem_entities);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;
                                        if(propensite > 0){
                                            if(cellu[current_poly].crystalline[current_bond] == false){
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                            }
                                            else{
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                            }
                                        }
                                    }
                                    else if(cellu[current_poly].len_poly == 2){
                                        reaction_type = 2;
                                        propensite = prop(par,reaction_type,chem_entities);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;
                                        if(propensite > 0){
                                           if(cellu[current_poly].crystalline[current_bond] == false){
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                            }
                                            else{
                                                addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                            }                            
                                        }
                                    }
                                    else if(cellu[current_poly].len_poly > 2){
                                        if(current_bond > 1 and current_bond < cellu[current_poly].len_poly-2){
                                            reaction_type = 1;
                                            propensite = prop(par,reaction_type,chem_entities);
                                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                                            if(propensite > 0){
                                               if(cellu[current_poly].crystalline[current_bond] == false){
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                                }
                                                else{
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                                }                            
                                            }

                                        }
                                        if(current_bond == 1 or current_bond == cellu[current_poly].len_poly-2){
                                            reaction_type = 2;
                                            propensite = prop(par,reaction_type,chem_entities);
                                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                                            if(propensite > 0){
                                               if(cellu[current_poly].crystalline[current_bond] == false){
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,propensite);
                                                }
                                                else{
                                                    addreaction(true, Table_cellu[current_poly],current_poly, current_bond,1, reaction_type,crystal_propensite);
                                                }                            
                                            }
                                        }

                                    }

                                }
                            }
                            else if(current_substrate == 2){
                                if(hemi[current_poly].status[current_bond] == -1){
                                    hemi[current_poly].status[current_bond] = 1;
                                    hemi[current_poly].status[current_bond] = 1;
                                    reaction_type = 4;
                                    propensite = prop(par,reaction_type,chem_entities);
                                    crystal_propensite = par.crystal_modifier_hemi * propensite;
                                    
                                    if(propensite > 0 and current_bond < hemi[current_poly].len_poly){
                                       if(hemi[current_poly].crystalline[current_bond] == false){
                                            addreaction(true, Table_hemi[current_poly],current_poly, current_bond,2, reaction_type,propensite);
                                        }
                                        else{
                                            addreaction(true, Table_hemi[current_poly],current_poly, current_bond,2, reaction_type,crystal_propensite);
                                        }                            
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        bond_neighbors_hemi.erase(bond_key);
    }
//    remove_neighbor_from_vector(x,y,z,bond_neighbors_cellu,bond_neighbors_hemi);    
}



//Updates the reaction tables regarding inhibition through glucose or cellobiose
void update_reactiontables_inhib(vector<TList>& Table_cellu, vector<bList>& celluList, vector<TList>& Table_hemi, params& par, const int nbr_Glc_pdt, const int nbr_cellobiose, vector<float>& chem_entities){

//    double inhib_nr_EG = par.inhib_cellobiose_EG * nbr_cellobiose;
//    double inhib_nr_CBH = par.inhib_cellobiose_CBH * nbr_cellobiose;
//    double inhib_nr_BGL =  par.inhib_glucose_BGL * nbr_Glc_pdt;

    double inhib_nr_EG = par.inhib_cellobiose_EG * chem_entities[0]* nbr_cellobiose / (chem_entities[0] + chem_entities[1] + nbr_cellobiose);
    double inhib_nr_CBH = par.inhib_cellobiose_CBH * chem_entities[1]* nbr_cellobiose / (chem_entities[0] + chem_entities[1] + nbr_cellobiose);
    double inhib_nr_BGL =  par.inhib_glucose_BGL * chem_entities[2] * nbr_Glc_pdt / (chem_entities[2] + nbr_Glc_pdt);
//  double inhib_nr_XYL = par.inhib_xylose_XYL * chem_entities[3] * nbr_xyl_pdt / (chem_entities[3]+nbr_xyl_pdt);
    int index = 0;
//    cout << "updating reactionTables" << endl;
    for(int i = 0; i < Table_cellu.size(); i++){
        for(int j = 0; j < Table_cellu[i].indic_action.size(); j++){
            if(Table_cellu[i].indic_action[j] == 1){
                index = Table_cellu[i].num_bond[j];
                if(celluList[i].crystalline[index] == false)
                    Table_cellu[i].liste_prop[j] = par.k1*(chem_entities[0]-inhib_nr_EG);//Table_cellu[i].prop_uninhib[j] / double(1+nbr_cellobiose*inhib_cellobiose_EG);
                else
                    Table_cellu[i].liste_prop[j] = par.crystal_modifier_cellu * par.k1*(chem_entities[0]-inhib_nr_EG);//Table_cellu[i].prop_uninhib[j] / double(1+nbr_cellobiose*inhib_cellobiose_EG);
                if(Table_cellu[i].liste_prop[j] < 0){
                    Table_cellu[i].liste_prop[j] = 0;
                }
            }
            else if(Table_cellu[i].indic_action[j] == 2){
                index = Table_cellu[i].num_bond[j];
                if(celluList[i].crystalline[index] == false)
                    Table_cellu[i].liste_prop[j] = par.k2*(chem_entities[1]-inhib_nr_CBH);//Table_cellu[i].prop_uninhib[j] / double(1+nbr_cellobiose*inhib_cellobiose_CBH);
                else
                    Table_cellu[i].liste_prop[j] = par.crystal_modifier_cellu * par.k2*(chem_entities[1]-inhib_nr_CBH);//Table_cellu[i].prop_uninhib[j] / double(1+nbr_cellobiose*inhib_cellobiose_CBH);
                if(Table_cellu[i].liste_prop[j] < 0){
                    Table_cellu[i].liste_prop[j] = 0;
                }
            }
            else{
                index = Table_cellu[i].num_bond[j];
                if(celluList[i].crystalline[index] == false)
                    Table_cellu[i].liste_prop[j] = par.k3*(chem_entities[2]-inhib_nr_BGL);//Table_cellu[i].prop_uninhib[j] / double(1+nbr_Glc_pdt*inhib_glucose_BGL);
                else
                    Table_cellu[i].liste_prop[j] = par.crystal_modifier_cellu * par.k3*(chem_entities[2]-inhib_nr_BGL);//Table_cellu[i].prop_uninhib[j] / double(1+nbr_Glc_pdt*inhib_glucose_BGL);
                if(Table_cellu[i].liste_prop[j] < 0){
                    Table_cellu[i].liste_prop[j] = 0;
                }
            }
//        cout << "without inhib: " << Table_cellu[i].prop_uninhib[j] << "; with inhib: " <<  Table_cellu[i].liste_prop[j] << endl;
        }
        Table_cellu[i].calcTableProp();
    }

//End of function    
}


void chain_length_distribution(vector<bList>& polyList, vector<DPList>& chain_distrib, const int length_fibril, const double real_time, const int t0, const int index){

    chain_distrib[index].t0 = t0;
    chain_distrib[index].real_time = real_time;
    for(int i=0; i < polyList.size();i++){
        chain_distrib[index].DP[polyList[i].len_poly-1] += 1;
    }



//HEAT MAP!!!
//End of function    
}

//Checks whether all bonds with status = 1 have reactions associated to them 
bool check_reactions(vector<TList>& allTables, vector<bList>& polyList, params& par){
    bool checker = true;
    int test = 0;
    for(int i=0;i<polyList.size();i++){
        for(int j=0;j<polyList[i].len_poly;j++){
            if(polyList[i].status[j] == 1){
                test = 0;
                for(int k=0;k<allTables[i].nbr_element;k++){
                    if(allTables[i].num_bond[k] == j)
                        test = 1;
                }
                if(test == 0){
                    if(j > 1 and j < polyList[i].len_poly-1)
                        cout << " In function check_reactions(): Poly: " << i << "; Bond: " << j << "; its status: " << polyList[i].status[j] << "; it has no reaction associated to it" << endl;
                    checker =false;
                }
            }
        }
    }

return checker;
//End of function    
}

//Returns the bond type according to the type of left bond l and right bond r. Used for bond_type distribution on hemicellulose
int find_bond_type(int l, int r){
//1 for glc-glc, 2 for xyl-glc, 3 for glc-xyl, 4 for xyl-xyl, and -1 for lign-lign. For now this is only important for hemicellulose
    if(l == 1){
        if(r == 0){//Only the cases r = 0 and r = 1 need to be checked, because the bonds are distributed from the left and only bond_types 0 or 1 are in the system before
            return 3;
        }
        else{
            return 1;
        }
    }
    else if(l == 2){
        if(r == 0){
            return 3;
        }
        else{
            return 1;
        }
    }
    else if(l == 3){
        if(r == 0){
            return 4;
        }
        else{
            return 2;
        }
    }
    else if(l == 4){
        if(r == 0){
            return 4;
        }
        else{
            return 2;
        }
    }
    else{
        cout << "In function find_bond_type(): l is not 1,2,3 or 4; l = " << l << endl;
        return 42;
    }
//End of function    
}


//Returns true if the bond exists within the polyLists, and fills the reference arguments with the coordinates and the material. If it returns false, the reference arguments should not be used further and are set to -1
bool find_specific_bond(const vector<bList>& celluList, const vector<bList>& hemiList, const vector<bList>& lignList, const int x, const int y, const int z, int& index_poly, int& z_index, int& material){
    for(int i=0;i<hemiList.size();i++){
        if(hemiList[i].x == x and hemiList[i].y == y){
            if(hemiList[i].len_poly > 0){
                if(hemiList[i].z[0] <= z and hemiList[i].z[hemiList[i].len_poly-1] >= z){
                    material = 2;
                    index_poly = i;
                    z_index = z - hemiList[i].z[0];
                    return true;
                }

            }
/*            for(int j=0;j<hemiList[i].len_poly;j++){
                if(hemiList[i].z[j] == z){
                    material = 2;
                    index_poly = i;
                    z_index = j;
                    return true;
                }
            }*/
        }
    }

    for(int i=0;i<lignList.size();i++){
        if(lignList[i].x == x and lignList[i].y == y){
            if(lignList[i].len_poly > 0){
                if(lignList[i].z[0] <= z and lignList[i].z[lignList[i].len_poly-1] >= z){
                    material = 3;
                    index_poly = i;
                    z_index = z - lignList[i].z[0];
                    return true;
                }

            }
/*            for(int j=0;j<lignList[i].len_poly;j++){
                if(lignList[i].z[j] == z){
                    material = 3;
                    index_poly = i;
                    z_index = j;
                    return true;
                }
            }*/
        }
    }

    for(int i=0;i<celluList.size();i++){
        if(celluList[i].x == x and celluList[i].y == y){
            if(celluList[i].len_poly > 0){
                if(celluList[i].z[0] <= z and celluList[i].z[celluList[i].len_poly-1] >= z){
                    material = 1;
                    index_poly = i;
                    z_index = z - celluList[i].z[0];
                    return true;
                }

            }

/*            for(int j=0;j<celluList[i].len_poly;j++){
                if(celluList[i].z[j] == z){
                    material = 1;
                    index_poly = i;
                    z_index = j;
                    return true;
                }
            }*/
        }
    }
    index_poly = -1;
    material = -1;
    z_index = -1;
    return false;
//End of function
}

//Returns the largest number of uninterrupted "true" elements of free_positions 
int largest_gap(bool *free_positions, int size, int& gap_start){
    int gap_size = 0;
    int max_gap_size = 0;
    int i=0;
    while(i<size){
        if(gap_size == 0)
            gap_start = i;
        if(free_positions[i] == true){
            gap_size++;

            max_gap_size = findmax(gap_size,max_gap_size);
        }
        else{
            gap_size = 0;
        }
        i++;
    }
    return gap_size;
}


bool check_position(int check_x, int check_y, int ref_x, int ref_y){

    bool check_bool = false;

    if(check_x == ref_x and (check_y == ref_y + 1 or check_y == ref_y - 1 ))
        check_bool = true;
    else if(check_x == ref_x + 1 and (check_y == ref_y or check_y == ref_y + 1 or check_y == ref_y - 1))
        check_bool = true;
    else if(check_x == ref_x - 1 and (check_y == ref_y or check_y == ref_y + 1 or check_y == ref_y - 1))
        check_bool = true;
/*    else if(check_x == ref_x and check_y == ref_y)
        check_bool = true;*/
    else if(check_x == ref_x and (check_y == ref_y + 2 or check_y == ref_y - 2 ))
        check_bool = true;
    else if(check_x == ref_x + 2 and (check_y == ref_y or check_y == ref_y + 2 or check_y == ref_y - 2))
        check_bool = true;
    else if(check_x == ref_x - 2 and (check_y == ref_y or check_y == ref_y + 2 or check_y == ref_y - 2))
        check_bool = true;
/*    else if(check_x == ref_x and check_y == ref_y)
        check_bool = true;*/

    return check_bool;
//End of function    
}


//Looks through all reaction tables and looks at the relative fraction of reactions associated to each enzyme 
void check_reaction_table_distribution(const std::vector<TList>& Table_cellu, const std::vector<TList>& Table_hemi, TList& Table_lign, params& par){
    double a0 = 0;
    double EG_fraction = 0;
    double CBH_fraction = 0;
    double BGL_fraction = 0;
    double XYL_fraction = 0;
    double lign_fraction = 0;

    for(int i=0;i<Table_cellu.size();i++){
        a0 += Table_cellu[i].prop_sum;
        for(int j=0;j<Table_cellu[i].nbr_element;j++){
            if(Table_cellu[i].indic_action[j] == 1)
                EG_fraction+=Table_cellu[i].liste_prop[j];
            else if(Table_cellu[i].indic_action[j] == 2)
                CBH_fraction+=Table_cellu[i].liste_prop[j];
            else if(Table_cellu[i].indic_action[j] == 3)
                BGL_fraction+=Table_cellu[i].liste_prop[j];
            else{
                cout << "In function check_reaction_table_distribution: Table_cellu contains reactions not associated to EG, CBH or BGL" << endl;
            }            
        }
    }
    for(int i=0;i<Table_hemi.size();i++){
        a0 += Table_hemi[i].prop_sum;
        for(int j=0;j<Table_hemi[i].nbr_element;j++){
            if(Table_hemi[i].indic_action[j] == 4)
                XYL_fraction += Table_hemi[i].liste_prop[j];
            else{
                cout << "In function check_reaction_table_distribution: Table_hemi contains reactions not associated to xylanase" << endl;
            }
        }
    }
    lign_fraction = Table_lign.prop_sum;
    a0 += Table_lign.prop_sum;
    if(a0 > 0){
        EG_fraction /= a0;
        CBH_fraction /= a0;
        BGL_fraction /= a0;
        XYL_fraction /= a0;
        lign_fraction /= a0;
    }
    par.EG_fraction.push_back(EG_fraction);
    par.CBH_fraction.push_back(CBH_fraction);
    par.BGL_fraction.push_back(BGL_fraction);
    par.XYL_fraction.push_back(XYL_fraction);
    par.lign_fraction.push_back(lign_fraction);
    if(par.time_mean.size() > 0)
        par.timestamp.push_back(par.time_mean[par.time_mean.size()-1]);
    else 
        par.timestamp.push_back(0);
//End of function    
}

//Returns the smaller of a and b
int findmin(int a, int b){
    if(a < b)
        return a;
    else if(b < a)
        return b;
    else if(b == a)
        return a;
    else{
        cout << "In function findmin: neither a<b, a>b, or a== b are true" << endl;
        return 0;
    }
//End of function    
}

//Returns the larger of a and b
int findmax(int a, int b){
    if(a < b)
        return b;
    else if(b < a)
        return a;
    else if(b == a)
        return a;
    else{
        cout << "In function findmin: neither a<b, a>b, or a== b are true" << endl;
        return 0;
    }
//End of function    
}


void lignin_glue(params& par, vector<TList>& Table_cellu, vector<TList>& Table_hemi, TList& Table_lign, vector<float>& chem_entities, int nbr_poly_lign){
    int test1 = 0;
    double propensite = 0;
    double test = 0;
    int enzyme_number = 0;
    for(int i=0;i< chem_entities.size()-1;i++){
        if(chem_entities[i] > 0)
            test1 = 1;
    }
    if(test1 == 1){
        test = drand48()*(chem_entities[0]+chem_entities[1]+chem_entities[2]+chem_entities[3]);
        if(par.verbose == true){
            cout << "Lignin glue: " << chem_entities[0] << "\t" << chem_entities[1] << "\t" << chem_entities[2] << "\t" << chem_entities[3] << "; sum = " << (chem_entities[0]+chem_entities[1]+chem_entities[2]+chem_entities[3]) << "; test = " << test << endl;
        }

        if(test <= chem_entities[0])
            enzyme_number = 1;
        else if(test <= chem_entities[0] + chem_entities[1])
            enzyme_number = 2;
        else if(test <= chem_entities[0] + chem_entities[1] + chem_entities[2])
            enzyme_number = 3;
        else //if(test <= chem_entities[3])
            enzyme_number = 4;


        if(enzyme_number == 1){
            chem_entities[enzyme_number-1] -=1;//* chem_entities[enzyme_number]/double(par.init_EG);
            par.N_enzymes_glued++;
        }
        if(enzyme_number == 2){
            chem_entities[enzyme_number-1] -=1;// * chem_entities[enzyme_number]/double(par.init_CBH);
            par.N_enzymes_glued++;
        }
        if(enzyme_number == 3){
            chem_entities[enzyme_number-1] -=1;// * chem_entities[enzyme_number]/double(par.init_BGL);
            par.N_enzymes_glued++;
        }
        if(enzyme_number == 4){
            chem_entities[enzyme_number-1] -=1;// * chem_entities[enzyme_number]/double(par.init_XYL);
            par.N_enzymes_glued++;
        }

    

        if(chem_entities[enzyme_number-1] < 0)
            chem_entities[enzyme_number-1] = 0;
        propensite = prop(par,enzyme_number,chem_entities);
        if(propensite > 0){
            if(Table_cellu.size()>0){
                for(int i=0;i<Table_cellu.size();i++){
                    for(int j=0;j<Table_cellu[i].nbr_element;j++){
                        if(Table_cellu[i].indic_action[j] == enzyme_number){
                            Table_cellu[i].liste_prop[j] = propensite;
                            Table_cellu[i].prop_uninhib[j] = propensite;
                        }
                    }
                    Table_cellu[i].calcTableProp();
                }
            }
            if(Table_hemi.size()>0){
                for(int i=0;i<Table_hemi.size();i++){
                    for(int j=0;j<Table_hemi[i].nbr_element;j++){
                        if(Table_hemi[i].indic_action[j]==enzyme_number){
                            Table_hemi[i].liste_prop[j]= propensite;
                            Table_hemi[i].prop_uninhib[j]= propensite;
                        }
                    }
                    Table_hemi[i].calcTableProp();
                }
            }
        }
        else{
            if(Table_cellu.size()>0){
                for(int i=0;i<Table_cellu.size();i++){
                    for(int j=0;j<Table_cellu[i].liste_prop.size();j++){
                        if(Table_cellu[i].indic_action[j]==enzyme_number){
                            delete_specific_reaction(Table_cellu[i],Table_cellu,j);
                        }
                        Table_cellu[i].calcTableProp();
                    }
                }
            }
            if(Table_hemi.size()>0){
                for(int i=0;i<Table_hemi.size();i++){
                    for(int j=0;j<Table_hemi[i].liste_prop.size();j++){
                        if(Table_hemi[i].indic_action[j]==enzyme_number){
                            delete_specific_reaction(Table_hemi[i],Table_hemi,j);
                        }
                        Table_hemi[i].calcTableProp();
                    }
                }            

            }    
        }
    }
    else{
        cout << "In function lignin_glue: all propensities are 0, except for the glue function" << endl;
        for(int i=0;i<Table_cellu.size();i++){
            Table_cellu[i].nbr_element = 0;
        }
        for(int i=0;i<Table_hemi.size();i++){
            Table_hemi[i].nbr_element = 0;
        }

    }

    int lignins_glued = (2.*float(par.enzyme_radius) + 3.) + int(drand48() * (2.*float(par.enzyme_radius) + 3.));//the enzyme blocks up to 2*r_enzyme + 1 lignin bonds
//    lignins_glued = par.nbr_monolignol;
    if(lignins_glued == 0)
        lignins_glued = 1;
    if(par.verbose == true){
        cout << "lignins glued: " << lignins_glued << endl;
    }
    if(lignins_glued > 0)
        par.nbr_lignin_blocked += lignins_glued;
    else
        par.nbr_lignin_blocked += 1;
    if (par.nbr_lignin_blocked >= par.nbr_monolignol){
        chem_entities.at(4) = 0;
    }
    else{
        chem_entities.at(4)=par.nbr_monolignol - par.nbr_lignin_blocked;//+ chem_entities[0] + chem_entities[1] + chem_entities[2] + chem_entities[3];
    }
    propensite = prop(par,5,chem_entities);
    Table_lign.liste_prop[0] = propensite;
    Table_lign.prop_uninhib[0] = propensite;
    Table_lign.prop_sum = propensite; 
//    cout << Table_lign.liste_prop[0] << endl;
//End of function    
}


//Returns a normally distributed random number between 0 and 1 with mean value at 0.5, according to the box-muller algorithm
double box_muller(){

    double z0 = 0;
    double rnum1 = drand48();
    double rnum2 = drand48();
    double sigma = 0.5;
    double mu = 0.5;


    z0=sqrt(-2*log(rnum1))*cos(2*M_PI*rnum2);
    z0=(sigma*z0+mu);
    while(z0 < 0 or z0 >= (2*mu)){
        rnum1 = drand48();//Generates random numbers between 0 and 1
        rnum2 = drand48();
        z0=sqrt(-2*log(rnum1))*cos(2*M_PI*rnum2);
        z0=(sigma*z0+mu);           
    }
//    cout << z0 << endl;
    return z0;
//End of function    
}

int find_poly(const vector<bList>& polyList, const int x, const int y, const int z){
//        cout << "Finding poly; x,y,z = " << x << "," << y << "," << z <<  endl;
    for(int i=0;i<polyList.size();i++){
        if(polyList[i].x == x and polyList[i].y == y and polyList[i].len_poly > 0){
            if(polyList[i].z[0] < z){
                if(polyList[i].len_poly > 1){
                    if(polyList[i].z[polyList[i].len_poly-1] >= z){
  //                      cout << "Found poly " << endl;
                        return i;
                    }
                }
            }
            else if(polyList[i].z[0] == z)
                return i;
        }
    }
    return -1;//This is the case if there is a lignin polymer at the searched position layer or the searched position is outside the range of possible positions (for example when a hole is created at a corner of the hull)
//End of function    
}



int find_bond(const bList& poly, int z_selected){
    if(poly.len_poly == 0)
        return -1;
    else if(poly.len_poly == 1 and poly.z[0] == z_selected)
        return 0;
 //   cout << "Finding bond; len_poly:  " << poly.len_poly << "; z[0] = " << poly.z[0] << "; z[len-1] = " << poly.z[poly.len_poly-1] << endl;
    for(int i=0;i<poly.len_poly;i++){
        if(poly.z[i] == z_selected){
   //         cout << "Found bond" << endl;
            return i;
        }
    }

    cout << "In function find_bond(): z_selected is not contained in the argument polymer; z_selected = " << z_selected << "; z[0] = " << poly.z[0] << "; z[len_poly-1] = " << poly.z[poly.len_poly-1] << endl;
    return -1;
//End of function    
}

void hole_cutter(vector<bList>& celluList, vector<bList>& hemiList, vector<bList>& lignList, vector<TList>& Table_hemi, params& par, int& difference_hemi, int& difference_lign, int& nbr_poly_hemi, int& nbr_poly_lign, bool& error_bool, int material){



    if(difference_hemi == 0 and difference_lign == 0){
        return;
    }

    unordered_map<int,vector<int>> possible_coordinates;
    vector<int> keys;
    vector<int> coordinates_to_add;
    coordinates_to_add.push_back(0);
    coordinates_to_add.push_back(0);
    int key_to_add = 0;

    for(int i=0; i < nbr_poly_hemi; i++){
        key_to_add = cantor_pair_two(hemiList[i].x, hemiList[i].y);
        if(possible_coordinates.count(key_to_add) == 0){
            keys.push_back(key_to_add);
            coordinates_to_add[0] = hemiList[i].x;
            coordinates_to_add[1] = hemiList[i].y;
            possible_coordinates.insert(make_pair(key_to_add, coordinates_to_add));
        }
    }
    for(int i=0; i < nbr_poly_lign; i++){
        key_to_add = cantor_pair_two(lignList[i].x, lignList[i].y);
        if(possible_coordinates.count(key_to_add) == 0){
            keys.push_back(key_to_add);
            coordinates_to_add[0] = lignList[i].x;
            coordinates_to_add[1] = lignList[i].y;
            possible_coordinates.insert(make_pair(key_to_add, coordinates_to_add));
        }
    }


    int max_hole_size = int(0.1*par.length_fibril);
    if(max_hole_size < 2*par.enzyme_radius)
        max_hole_size = 2*par.enzyme_radius;
    int z_hole_bottom = 0;
    int x_hole, y_hole;
    int z_min = 0;
    int z_max = par.length_fibril - 1;
    int hole_center_diameter = 0;
    int coordinate_index = 0;
    int diameter = 0;
    int left_right_dist = 0;
    vector<vector<int>> hole_coordinates;    


    while(difference_hemi > 0 or difference_lign > 0){
        if(par.verbose == true){
            cout << "hole cutting" << endl;
        }
        //Specify location of hole
        hole_coordinates.clear();
        z_hole_bottom = int(drand48() * (par.length_fibril-1));
        if(z_max - z_hole_bottom < max_hole_size){
            hole_center_diameter = z_max - z_hole_bottom;
        }
        else{
            hole_center_diameter = int(drand48() * max_hole_size);
        }

        coordinate_index = int(drand48() * keys.size());
        x_hole = possible_coordinates[keys[coordinate_index]][0];
        y_hole = possible_coordinates[keys[coordinate_index]][1];        



        //Build the vector containing the possible positions of bonds to be erased
        for(int i = 0; i < hole_center_diameter; i++){
            hole_coordinates.push_back(vector<int>());
            hole_coordinates[i].push_back(x_hole);
            hole_coordinates[i].push_back(y_hole);
            hole_coordinates[i].push_back(z_hole_bottom+i);
        }
        diameter = hole_center_diameter - 2;
        left_right_dist = 0;
         while(diameter > 0){
            left_right_dist++;
            for(int i=left_right_dist; i < diameter; i++){
            hole_coordinates.push_back(vector<int>());
            hole_coordinates[hole_coordinates.size()-1].push_back(x_hole+left_right_dist);
            hole_coordinates[hole_coordinates.size()-1].push_back(y_hole);
            hole_coordinates[hole_coordinates.size()-1].push_back(z_hole_bottom+i);
            hole_coordinates.push_back(vector<int>());
            hole_coordinates[hole_coordinates.size()-1].push_back(x_hole-left_right_dist);
            hole_coordinates[hole_coordinates.size()-1].push_back(y_hole);
            hole_coordinates[hole_coordinates.size()-1].push_back(z_hole_bottom+i);
            hole_coordinates.push_back(vector<int>());
            hole_coordinates[hole_coordinates.size()-1].push_back(x_hole);
            hole_coordinates[hole_coordinates.size()-1].push_back(y_hole+left_right_dist);
            hole_coordinates[hole_coordinates.size()-1].push_back(z_hole_bottom+i);
            hole_coordinates.push_back(vector<int>());
            hole_coordinates[hole_coordinates.size()-1].push_back(x_hole);
            hole_coordinates[hole_coordinates.size()-1].push_back(y_hole-left_right_dist);
            hole_coordinates[hole_coordinates.size()-1].push_back(z_hole_bottom+i);
            }
            diameter -= 2;
        }


        //Cut the hole

        int material = 0;
        int poly_index = 0;
        int z_index = 0;
        int new_poly_index = 0;
        for(int i=0;i<hole_coordinates.size();i++){
            if(find_specific_bond(celluList, hemiList, lignList,hole_coordinates[i][0], hole_coordinates[i][1], hole_coordinates[i][2], poly_index, z_index, material) == true){
                if(material == 2 and difference_hemi > 0){
                    if(z_index == hemiList[poly_index].len_poly-1){
                        taylorOldPoly(hemiList[poly_index],z_index);
                        difference_hemi--;
                    }
                    else if (z_index == 0)
                    {
                        taylorNewPoly(hemiList[poly_index],z_index);
                        difference_hemi--;
                    }
                    else{
                        hemiList.push_back(bList());
                        nbr_poly_hemi++;
                        new_poly_index = hemiList.size()-1;
                        hemiList[new_poly_index].x = hemiList[poly_index].x;
                        hemiList[new_poly_index].y = hemiList[poly_index].y;
                        hemiList[new_poly_index].set_z(hemiList[poly_index].z[z_index]);
                        Table_hemi.push_back(TList());
                        initTList(Table_hemi[Table_hemi.size()-1],Table_hemi.size()-1);        
                        //cout << Table_hemi.size() << endl;
                        hemiList[new_poly_index].reactionTable = Table_hemi.size() - 1;


                        for(int j = 0;j<hemiList[poly_index].len_poly;j++){
                            addbond(hemiList[new_poly_index],j,1,1,0);
                            hemiList[new_poly_index].z[j] = hemiList[poly_index].z[j];
                            hemiList[new_poly_index].status[j] = hemiList[poly_index].status[j];
                            hemiList[new_poly_index].bond_type[j] = hemiList[poly_index].bond_type[j];
                            hemiList[new_poly_index].crystalline[j] = hemiList[poly_index].crystalline[j];
                            hemiList[new_poly_index].N_blocked_positions[j] = hemiList[poly_index].N_blocked_positions[j];
                        }
                        taylorOldPoly(hemiList[poly_index],z_index);
                        taylorNewPoly(hemiList[new_poly_index],z_index);
                    }

                }
                else if(material == 3 and difference_lign > 0){
                    if(z_index == lignList[poly_index].len_poly-1){
                        taylorOldPoly(lignList[poly_index],z_index);
                        difference_lign--;
                    }
                    else if (z_index == 0)
                    {
                        taylorNewPoly(lignList[poly_index],z_index);
                        difference_lign--;
                    }
                    else{
                        lignList.push_back(bList());
                        nbr_poly_lign++;
                        new_poly_index = lignList.size()-1;
                        lignList[new_poly_index].x = lignList[poly_index].x;
                        lignList[new_poly_index].y = lignList[poly_index].y;
                        lignList[new_poly_index].set_z(lignList[poly_index].z[z_index]);


                        for(int j = 0;j<lignList[poly_index].len_poly;j++){
                            addbond(lignList[new_poly_index],j,1,1,0);
                            lignList[new_poly_index].z[j] = lignList[poly_index].z[j];
                            lignList[new_poly_index].status[j] = lignList[poly_index].status[j];
                            lignList[new_poly_index].bond_type[j] = lignList[poly_index].bond_type[j];
                            lignList[new_poly_index].crystalline[j] = lignList[poly_index].crystalline[j];
                            lignList[new_poly_index].N_blocked_positions[j] = lignList[poly_index].N_blocked_positions[j];
                        }
                        taylorOldPoly(lignList[poly_index],z_index);
                        taylorNewPoly(lignList[new_poly_index],z_index);
                    }
                }
            }

        }


        if(difference_hemi < 0){
            difference_hemi = 0;
        }
        if(difference_lign < 0){
            difference_lign = 0;
        }
    }

cout << "End of hole-cutter function; difference_hemi = " << difference_hemi << "; difference_lign = " << difference_lign << endl;
//END OF FUNCTION

}



int free_bonds(int enzyme_radius, int mode_code){
    if(enzyme_radius == 0)
        return 1;
    int square_diameter = 2*enzyme_radius + 1;
    int nbr_holes_per_plane = square_diameter * enzyme_radius;//We only need half of bonds free
    return(0.5*nbr_holes_per_plane * square_diameter);

//End of function    
}


void count_free_neighbors(std::vector<bList>& celluList, std::vector<bList>& hemiList, std::vector<bList>& lignList, params& par, bool& errorbool){
    int check_radius = 1;

    if(par.enzyme_radius > 0){
        check_radius = par.enzyme_radius;//Number of bonds which are checked around each position
    }
    int free_req = par.free_bonds_req;
    int dx = 0.8;//dx,dy ca 0.8 nm, according to Gibson et al "The hierarchical structure and mechanics of plant materials", 2012
    int dy = 0.8;
    int dz = 1;//roughly the diameter of a glucose molecule. This is the distance between two bonds
    int min_x_cellu = 0;
    int max_x_cellu = 0;
    int min_y_cellu = 0;
    int max_y_cellu = 0;
    int nbr_poly_cellu = celluList.size();
    int nbr_poly_hemi = hemiList.size();
    int nbr_poly_lign = lignList.size();
    int x_pos_plus = 0;
    int y_pos_plus = 0;
    int z_pos_plus = 0;
    int x_pos_minus = 0;
    int y_pos_minus = 0;
    int z_pos_minus = 0;

    bool x_plus = false;//This is true if the last bond checked was free
    bool y_plus = false;
    bool z_plus = false;        
    bool x_minus = false;
    bool y_minus = false;
    bool z_minus = false;


    cout << free_req << endl;



    int min_x_hemi = 0;
    int max_x_hemi = 0;
    int min_y_hemi = 0;
    int max_y_hemi = 0;

    min_x_hemi = hemiList[0].x;
    max_x_hemi = hemiList[0].x;
    min_y_hemi = hemiList[0].y;
    max_y_hemi = hemiList[0].y;



    int min_x_lign = 0;
    int max_x_lign = 0;
    int min_y_lign = 0;
    int max_y_lign = 0;



    min_x_lign = lignList[0].x;
    max_x_lign = lignList[0].x;
    min_y_lign = lignList[0].y;
    max_y_lign = lignList[0].y;    




    for(int i=1;i < nbr_poly_lign;i++){
        min_x_lign = findmin(min_x_lign,lignList[i].x);
        max_x_lign = findmax(max_x_lign,lignList[i].x);            
        min_y_lign = findmin(min_y_lign,lignList[i].y);
        max_y_lign = findmax(max_y_lign,lignList[i].y);                        
    }

    min_x_cellu = celluList[0].x;
    max_x_cellu = celluList[0].x;
    min_y_cellu = celluList[0].y;
    max_y_cellu = celluList[0].y;
    for(int i=1;i < nbr_poly_cellu;i++){
        min_x_cellu = findmin(min_x_cellu,celluList[i].x);
        max_x_cellu = findmax(max_x_cellu,celluList[i].x);            
        min_y_cellu = findmin(min_y_cellu,celluList[i].y);
        max_y_cellu = findmax(max_y_cellu,celluList[i].y);                        
    }

    for(int i=1;i < nbr_poly_hemi;i++){
        min_x_hemi = findmin(min_x_hemi,hemiList[i].x);
        max_x_hemi = findmax(max_x_hemi,hemiList[i].x);            
        min_y_hemi = findmin(min_y_hemi,hemiList[i].y);
        max_y_hemi = findmax(max_y_hemi,hemiList[i].y);                        
    }


    int min_x_overall = findmin(min_x_hemi,min_x_lign);
    int max_x_overall = findmax(max_x_hemi,max_x_lign);
    int min_y_overall = findmin(min_y_hemi, min_y_lign);
    int max_y_overall = findmax(max_y_hemi,max_y_lign);

    int min_x_inner = par.min_x_outer + 1;
    int max_x_inner = par.max_x_outer - 1;
    int min_y_inner = par.min_y_outer + 1;
    int max_y_inner = par.max_y_outer -1;




//End of function
}

/*
bool is_outer_poly(const std::unordered_map<int, neighborList>& bond_neighbors_cellu, const std::unordered_map<int, neighborList>& bond_neighbors_hemi, const vector<bList>& hemiList, const vector<bList>& lignList, const vector<bList>& celluList, const int x, const int y, const int z, const int substrate, const params& par){

    if(par.verbose == true){
        cout << "Entered function is_outer_poly()" << endl;
    }


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


    int distance_x = par.max_x_overall - par.min_x_overall;
    int distance_y = par.max_y_overall - par.min_y_overall;
    int distance = findmax(distance_x,distance_y);


//    int stop_count = 0;//Counts, how many neighbors lie in the cardinal directions. If less than 3, it is an outer polymer

    int dir1_x,dir1_y,dir2_x,dir2_y,dir3_x,dir3_y,dir4_x,dir4_y;
    int dir5_x,dir5_y,dir6_x,dir6_y,dir7_x,dir7_y,dir8_x,dir8_y;
    bool dir1_stop = false;
    bool dir2_stop = false;    
    bool dir3_stop = false;
    bool dir4_stop = false;
    bool dir5_stop = false;
    bool dir6_stop = false;
    bool dir7_stop = false;
    bool dir8_stop = false;


    dir1_stop = true;//Currently, only the diagonals are considered
    dir2_stop = true;
    dir3_stop = true;
    dir4_stop = true;



    for(int i=1;i<distance;i++){

        dir1_x = x+i;
        dir1_y = y;
        dir2_x = x-i;
        dir2_y = y;
        dir3_x = x;
        dir3_y = y+i;
        dir4_x = x;
        dir4_y = y-i;

        dir5_x = x+i;
        dir5_y = y+i;
        dir6_x = x-i;
        dir6_y = y-i;
        dir7_x = x+i;
        dir7_y = y-i;
        dir8_x = x-i;
        dir8_y = y+i;

//INCLUDE Z-DISTANCE TO CHECK
        for(int j=0;j<celluList.size();j++){
            if(celluList[j].len_poly > 0){
                if(dir1_stop == false){
                    if(celluList[j].x == dir1_x and celluList[j].y == dir1_y){
                        if(celluList[j].z[0] <= z and celluList[j].z[celluList[j].z.size()-1] >= z){
                            dir1_stop = true;
                        }
                    }
                }
                if(dir2_stop == false){
                    if(celluList[j].x == dir2_x and celluList[j].y == dir2_y){
                        if(celluList[j].z[0] <= z and celluList[j].z[celluList[j].z.size()-1] >= z){
                            dir2_stop = true;
                        }                    
                    }                    
                }
                if(dir3_stop == false){
                    if(celluList[j].x == dir3_x and celluList[j].y == dir3_y){
                        if(celluList[j].z[0] <= z and celluList[j].z[celluList[j].z.size()-1] >= z){
                            dir3_stop = true;
                        }                    
                    }                    
                }
                if(dir4_stop == false){
                    if(celluList[j].x == dir4_x and celluList[j].y == dir4_y){
                        if(celluList[j].z[0] <= z and celluList[j].z[celluList[j].z.size()-1] >= z){
                            dir4_stop = true;
                        }                    
                    }                    
                }
                if(dir5_stop == false){
                    if(celluList[j].x == dir5_x and celluList[j].y == dir5_y){
                        if(celluList[j].z[0] <= z and celluList[j].z[celluList[j].z.size()-1] >= z){
                            dir5_stop = true;
                        }
                    }
                }
                if(dir6_stop == false){
                    if(celluList[j].x == dir6_x and celluList[j].y == dir6_y){
                        if(celluList[j].z[0] <= z and celluList[j].z[celluList[j].z.size()-1] >= z){
                            dir6_stop = true;
                        }                    
                    }                    
                }
                if(dir7_stop == false){
                    if(celluList[j].x == dir7_x and celluList[j].y == dir7_y){
                        if(celluList[j].z[0] <= z and celluList[j].z[celluList[j].z.size()-1] >= z){
                            dir7_stop = true;
                        }                    
                    }                    
                }
                if(dir8_stop == false){
                    if(celluList[j].x == dir8_x and celluList[j].y == dir8_y){
                        if(celluList[j].z[0] <= z and celluList[j].z[celluList[j].z.size()-1] >= z){
                            dir8_stop = true;
                        }                    
                    }                    
                }

            }
        }

        for(int j=0;j<hemiList.size();j++){
            if(hemiList[j].len_poly > 0){
                if(dir1_stop == false){
                    if(hemiList[j].x == dir1_x and hemiList[j].y == dir1_y){
                        if(hemiList[j].z[0] <= z and hemiList[j].z[hemiList[j].z.size()-1] >= z){
                            dir1_stop = true;
                        }
                    }
                }
                if(dir2_stop == false){
                    if(hemiList[j].x == dir2_x and hemiList[j].y == dir2_y){
                        if(hemiList[j].z[0] <= z and hemiList[j].z[hemiList[j].z.size()-1] >= z){
                            dir2_stop = true;
                        }                    
                    }                    
                }
                if(dir3_stop == false){
                    if(hemiList[j].x == dir3_x and hemiList[j].y == dir3_y){
                        if(hemiList[j].z[0] <= z and hemiList[j].z[hemiList[j].z.size()-1] >= z){
                            dir3_stop = true;
                        }                    
                    }                    
                }
                if(dir4_stop == false){
                    if(hemiList[j].x == dir4_x and hemiList[j].y == dir4_y){
                        if(hemiList[j].z[0] <= z and hemiList[j].z[hemiList[j].z.size()-1] >= z){
                            dir4_stop = true;
                        }                    
                    }                    
                }
                if(dir5_stop == false){
                    if(hemiList[j].x == dir5_x and hemiList[j].y == dir5_y){
                        if(hemiList[j].z[0] <= z and hemiList[j].z[hemiList[j].z.size()-1] >= z){
                            dir5_stop = true;
                        }
                    }
                }
                if(dir6_stop == false){
                    if(hemiList[j].x == dir6_x and hemiList[j].y == dir6_y){
                        if(hemiList[j].z[0] <= z and hemiList[j].z[hemiList[j].z.size()-1] >= z){
                            dir6_stop = true;
                        }                    
                    }                    
                }
                if(dir7_stop == false){
                    if(hemiList[j].x == dir7_x and hemiList[j].y == dir7_y){
                        if(hemiList[j].z[0] <= z and hemiList[j].z[hemiList[j].z.size()-1] >= z){
                            dir7_stop = true;
                        }                    
                    }                    
                }
                if(dir8_stop == false){
                    if(hemiList[j].x == dir8_x and hemiList[j].y == dir8_y){
                        if(hemiList[j].z[0] <= z and hemiList[j].z[hemiList[j].z.size()-1] >= z){
                            dir8_stop = true;
                        }                    
                    }                    
                }
            }
        }

        for(int j=0;j<lignList.size();j++){
            if(lignList[j].len_poly > 0){
                if(dir1_stop == false){
                    if(lignList[j].x == dir1_x and lignList[j].y == dir1_y){
                        if(lignList[j].z[0] <= z and lignList[j].z[lignList[j].z.size()-1] >= z){
                            dir1_stop = true;
                        }
                    }
                }
                if(dir2_stop == false){
                    if(lignList[j].x == dir2_x and lignList[j].y == dir2_y){
                        if(lignList[j].z[0] <= z and lignList[j].z[lignList[j].z.size()-1] >= z){
                            dir2_stop = true;
                        }                    
                    }                    
                }
                if(dir3_stop == false){
                    if(lignList[j].x == dir3_x and lignList[j].y == dir3_y){
                        if(lignList[j].z[0] <= z and lignList[j].z[lignList[j].z.size()-1] >= z){
                            dir3_stop = true;
                        }                    
                    }                    
                }
                if(dir4_stop == false){
                    if(lignList[j].x == dir4_x and lignList[j].y == dir4_y){
                        if(lignList[j].z[0] <= z and lignList[j].z[lignList[j].z.size()-1] >= z){
                            dir4_stop = true;
                        }                    
                    }                    
                }
                if(dir5_stop == false){
                    if(lignList[j].x == dir5_x and lignList[j].y == dir5_y){
                        if(lignList[j].z[0] <= z and lignList[j].z[lignList[j].z.size()-1] >= z){
                            dir5_stop = true;
                        }
                    }
                }
                if(dir6_stop == false){
                    if(lignList[j].x == dir6_x and lignList[j].y == dir6_y){
                        if(lignList[j].z[0] <= z and lignList[j].z[lignList[j].z.size()-1] >= z){
                            dir6_stop = true;
                        }                    
                    }                    
                }
                if(dir7_stop == false){
                    if(lignList[j].x == dir7_x and lignList[j].y == dir7_y){
                        if(lignList[j].z[0] <= z and lignList[j].z[lignList[j].z.size()-1] >= z){
                            dir7_stop = true;
                        }                    
                    }                    
                }
                if(dir8_stop == false){
                    if(lignList[j].x == dir8_x and lignList[j].y == dir8_y){
                        if(lignList[j].z[0] <= z and lignList[j].z[lignList[j].z.size()-1] >= z){
                            dir8_stop = true;
                        }                    
                    }                    
                }
            }
        }


        if(dir1_stop == true and dir2_stop == true and dir3_stop == true and dir4_stop == true and dir5_stop == true and dir6_stop == true and dir7_stop == true and dir8_stop == true){
            return false;
        }
    }

    if(dir1_stop == true and dir2_stop == true and dir3_stop == true and dir4_stop == true and dir5_stop == true and dir6_stop == true and dir7_stop == true and dir8_stop == true){
        return false;
    }
    else{
        return true;
    }

//End of function    
}*/

//Returns true if x_pos and y_pos are found as a valid combination inside celluList, hemiList or lignList
bool is_valid_pos(const vector<bList>& celluList, const vector<bList>& hemiList, const vector<bList>& lignList,int x_pos, int y_pos){

    for(int i=0;i<celluList.size();i++){
        if(celluList[i].x == x_pos and celluList[i].y == y_pos){
            return true;
        }
    }
    for(int i=0;i<hemiList.size();i++){
        if(hemiList[i].x == x_pos and hemiList[i].y == y_pos){
            return true;
        }  
    }
    for(int i=0;i<lignList.size();i++){
        if(lignList[i].x == x_pos and lignList[i].y == y_pos){
            return true;
        }
    }


    return false;
//End of function    
}

//Returns true if a bond exists at the specified x-, y- and z coordinates
bool bond_exists(const vector<bList>& celluList, const vector<bList>& hemiList, const vector<bList>& lignList,int x_pos, int y_pos, int z_pos){

    for(int i=0;i<celluList.size();i++){
        if(celluList[i].x == x_pos and celluList[i].y == y_pos){
            if(celluList[i].len_poly > 0){
                if(celluList[i].z[0] <= z_pos and celluList[i].z[celluList[i].z.size()-1] >= z_pos){
                    return true;
                }
            }
        }
    }
    for(int i=0;i<hemiList.size();i++){
        if(hemiList[i].x == x_pos and hemiList[i].y == y_pos){
            if(hemiList[i].len_poly > 0){
                if(hemiList[i].z[0] <= z_pos and hemiList[i].z[hemiList[i].z.size()-1] >= z_pos){
                    return true;
                }
            }
        }  
    }
    for(int i=0;i<lignList.size();i++){
        if(lignList[i].x == x_pos and lignList[i].y == y_pos){
            if(lignList[i].len_poly > 0){
                if(lignList[i].z[0] <= z_pos and lignList[i].z[lignList[i].z.size()-1] >= z_pos){
                    return true;
                }
            }
        }
    }


    return false;
//End of function    

}

//Returns true if there is no outer bond blocking the one at position (x,y,z). This function assumes integer coordinates for the bonds
/*bool is_outer_bond(const std::vector<bList>& celluList, const std::vector<bList>& hemiList, const std::vector<bList>& lignList, const int x, const int y, const int z, const int substrate, const params& par){

    if(z == 0 or z==par.length_fibril-1){
        return true;
    }

    float gradient_x = x-par.mid_x;//x_coordinate of vector pointing outwards from the microfibril center and intersecting with (x,y)
    float gradient_y = x-par.mid_y;//y_coordinate of vector pointing outwards from the microfibril center and intersecting with (x,y)
    float r = sqrt((x-par.mid_x)*(x-par.mid_x) + (y-par.mid_y)*(y-par.mid_y));//Length of gradient vector

    gradient_x/=r;
    gradient_y/=r;


    float x_check = x + gradient_x*(par.enzyme_radius + par.r_monomer);
    float y_check = y + gradient_y*(par.enzyme_radius + par.r_monomer);

    neighborList radial_neighbors = neighborList(x_check, y_check, z, par.enzyme_radius, par, celluList, hemiList, lignList);

    if (radial_neighbors.N_neighbors == 0){
        return true;
    }
    else{
        return false;
    }
//End of function
}*/


void min_max_cellu_hemi_lign(const vector<bList>& celluList, const vector<bList>& hemiList, const vector<bList>& lignList, params& par){

    par.min_x_cellu = 0;
    par.max_x_cellu = 0;
    par.min_y_cellu = 0;
    par.max_y_cellu = 0;
    par.min_x_lign = 0;
    par.max_x_lign = 0;
    par.min_y_lign = 0;
    par.max_y_lign = 0;
    par.min_x_hemi = 0;
    par.max_x_hemi = 0;
    par.min_y_hemi = 0;
    par.max_y_hemi = 0;

    par.mid_x = 0;
    par.mid_y = 0;


    if(celluList.size()>0){

        par.min_x_cellu = celluList[0].x;
        par.max_x_cellu = celluList[0].x;
        par.min_y_cellu = celluList[0].y;
        par.max_y_cellu = celluList[0].y;

        par.mid_x += celluList[0].x;
        par.mid_y += celluList[0].y;

        for(int i=1;i < celluList.size();i++){
            par.min_x_cellu = findmin(par.min_x_cellu,celluList[i].x);
            par.max_x_cellu = findmax(par.max_x_cellu,celluList[i].x);            
            par.min_y_cellu = findmin(par.min_y_cellu,celluList[i].y);
            par.max_y_cellu = findmax(par.max_y_cellu,celluList[i].y);                        
            par.mid_x += celluList[i].x;
            par.mid_y += celluList[i].y;
        }
    }

    if(hemiList.size()>0){

        par.min_x_hemi = hemiList[0].x;
        par.max_x_hemi = hemiList[0].x;
        par.min_y_hemi = hemiList[0].y;
        par.max_y_hemi = hemiList[0].y;
        par.mid_x += hemiList[0].x;
        par.mid_y += hemiList[0].y;

        for(int i=1;i < hemiList.size();i++){
            par.min_x_hemi = findmin(par.min_x_hemi,hemiList[i].x);
            par.max_x_hemi = findmax(par.max_x_hemi,hemiList[i].x);            
            par.min_y_hemi = findmin(par.min_y_hemi,hemiList[i].y);
            par.max_y_hemi = findmax(par.max_y_hemi,hemiList[i].y);                        
            par.mid_x += hemiList[i].x;
            par.mid_y += hemiList[i].y;
        }
    }

    if(lignList.size()>0){

        par.min_x_lign = lignList[0].x;
        par.max_x_lign = lignList[0].x;
        par.min_y_lign = lignList[0].y;
        par.max_y_lign = lignList[0].y;    
        par.mid_x += lignList[0].x;
        par.mid_y += lignList[0].y;


        for(int i=1;i < lignList.size();i++){
            par.min_x_lign = findmin(par.min_x_lign,lignList[i].x);
            par.max_x_lign = findmax(par.max_x_lign,lignList[i].x);            
            par.min_y_lign = findmin(par.min_y_lign,lignList[i].y);
            par.max_y_lign = findmax(par.max_y_lign,lignList[i].y);                        
            par.mid_x += lignList[i].x;
            par.mid_y += lignList[i].y;
        }
    }

    par.min_x_overall = findmin(par.min_x_hemi,par.min_x_lign);
    par.max_x_overall = findmin(par.max_x_hemi,par.max_x_lign);    
    par.min_y_overall = findmin(par.min_y_hemi,par.min_y_lign);
    par.max_y_overall = findmin(par.max_y_hemi,par.max_y_lign);

    if(lignList.size() > 0 or hemiList.size() > 0 or celluList.size() > 0){
        par.mid_x /= (celluList.size() + hemiList.size() + lignList.size());
        par.mid_y /= (celluList.size() + hemiList.size() + lignList.size());        
    }

    par.mid_z = 0.5*par.length_fibril;
//End of function
}

//Returns an integer which is unique for every set (a,b,c)
int cantor_pair_three(int a, int b, int c){
    int offset = 10;
    int two_pair = .5*(a+offset + b+offset)*(a+offset + b+offset + 1) + b+offset;
    return .5*(two_pair + c+offset)*(two_pair + c+offset + 1) + c+offset;    
//End of function
}

int cantor_pair_two(int a, int b){
    int offset = 10;
    int two_pair = 0.5*(a+offset + b+offset)*(a+offset + b+offset + 1) + b+offset;
    return two_pair;
//End of function
}

//Removes a bond from the neighbor vector
void remove_neighbor_from_vector(const int x_pos, const int y_pos, const int z_pos, unordered_map<int,neighborList>& bond_neighbors_cellu, unordered_map<int,neighborList>& bond_neighbors_hemi){


    int index_bond = cantor_pair_three(x_pos, y_pos, z_pos);
    int index_neighbor = 0;
//    char index_neighbor[3];


    if(bond_neighbors_cellu.find(index_bond) != bond_neighbors_cellu.end()){
        if(bond_neighbors_cellu[index_bond].N_neighbors == 0){
//            cout << "In function remove_neighbor_from_vector: a bond has no neighbors" << endl;
//            exit(1);
        }
        for(int j=0;j<bond_neighbors_cellu[index_bond].N_neighbors;j++){
            index_neighbor = cantor_pair_three(bond_neighbors_cellu[index_bond].x_neighbors[j],bond_neighbors_cellu[index_bond].y_neighbors[j],bond_neighbors_cellu[index_bond].z_neighbors[j]);
            if(bond_neighbors_cellu.find(index_neighbor) != bond_neighbors_cellu.end()){
                bond_neighbors_cellu[index_neighbor].remove_neighbor(x_pos,y_pos,z_pos);
            }
            else if(bond_neighbors_hemi.find(index_neighbor) != bond_neighbors_hemi.end()){
                bond_neighbors_hemi[index_neighbor].remove_neighbor(x_pos,y_pos,z_pos);
            }
        }
    //    cout << "Before " << bond_neighbors_cellu.size() << endl;
        bond_neighbors_cellu.erase(index_bond);
   //     cout << "After " << bond_neighbors_cellu.size() << endl;
    }
    else if(bond_neighbors_hemi.find(index_bond) != bond_neighbors_hemi.end()){
        for(int j=0;j<bond_neighbors_hemi[index_bond].N_neighbors;j++){
            index_neighbor = cantor_pair_three(bond_neighbors_hemi[index_bond].x_neighbors[j],bond_neighbors_hemi[index_bond].y_neighbors[j],bond_neighbors_hemi[index_bond].z_neighbors[j]);


            if(bond_neighbors_cellu.find(index_neighbor) != bond_neighbors_cellu.end()){
                bond_neighbors_cellu[index_neighbor].remove_neighbor(x_pos,y_pos,z_pos);
            }
            else if(bond_neighbors_hemi.find(index_neighbor) != bond_neighbors_hemi.end()){
                bond_neighbors_hemi[index_neighbor].remove_neighbor(x_pos,y_pos,z_pos);
            }
        }

        bond_neighbors_hemi.erase(index_bond);

    }






/*    for(int i=0;i<bond_neighbors_cellu.size();i++){
        if(bond_neighbors_cellu[i].x == x_pos and bond_neighbors_cellu[i].y == y_pos and bond_neighbors_cellu[i].z == z_pos){
            bond_neighbors_cellu.erase(bond_neighbors_hemi.begin()+i);
            i--;*/
/*            bond_neighbors_cellu[i].x_neighbors.clear();
            bond_neighbors_cellu[i].y_neighbors.clear();
            bond_neighbors_cellu[i].z_neighbors.clear();
            bond_neighbors_cellu[i].volume_contribution.clear();
            bond_neighbors_cellu[i].V_free = bond_neighbors_cellu[i].V_enzyme;
            bond_neighbors_cellu[i].N_neighbors = 0;*/

/*
        }
        else{
            for(int j=0;j<bond_neighbors_cellu[i].N_neighbors;j++){
                if (bond_neighbors_cellu[i].x_neighbors[j] == x_pos and bond_neighbors_cellu[i].y_neighbors[j] == y_pos and bond_neighbors_cellu[i].z_neighbors[j] == z_pos){
//                    cout << "Removing neighbor j = " << j << endl; 
                    bond_neighbors_cellu[i].remove_neighbor(j);
//                    cout << "Done"<< endl; 
                }
            }
        }
    }
    for(int i=0;i<bond_neighbors_hemi.size();i++){
        if(bond_neighbors_hemi[i].x == x_pos and bond_neighbors_hemi[i].y == y_pos and bond_neighbors_hemi[i].z == z_pos){
            bond_neighbors_hemi.erase(bond_neighbors_hemi.begin()+i);
            i--;
        }
        else{
            for(int j=0;j<bond_neighbors_hemi[i].N_neighbors;j++){
                if (bond_neighbors_hemi[i].x_neighbors[j] == x_pos and bond_neighbors_hemi[i].y_neighbors[j] == y_pos and bond_neighbors_hemi[i].z_neighbors[j] == z_pos){
                    bond_neighbors_hemi[i].remove_neighbor(j);
                }
            }
        }
    }
*/
//End of function
}



