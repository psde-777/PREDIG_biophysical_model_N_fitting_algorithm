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
#include "functions.hpp"
#include "bList.hpp"
#include "TList.hpp"
#include "DPList.hpp"
#include "neighborList.hpp"
#include "CBH_enzyme.hpp"
#include "params.hpp"
#include "tuple_hash.hpp"
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
void addbond(
    bList& list,
    int position,
    int state,
    int bond_type,
    bool crystalline)
{
    list.z.push_back(position);
    list.status.push_back(state);
    list.bond_type.push_back(bond_type);
    list.N_blocked_positions.push_back(0);
    list.crystalline.push_back(crystalline);
    list.covering.push_back(1);
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
void addreaction(params& par, bool confidence, std::vector<TList>& list,int poly,int bond_num,int mat,int act2,double prope, bool crystalline_bool){

    if(mat != 1 and mat != 2){
        cout << "In function addreaction: mat != 1 and mat != 2. This should not occur. Stopping" << endl;
        exit(1);
    }
//    if(act2 == 6){
//        cout << "Prope = " << prope << "; N_free_CBH = " << par.N_free_CBH <<  endl;
//    }
    if(prope >= 0.0){
        if(confidence == false){
            int test = 0;

            for(int i=0; i<list[poly].nbr_element; i++){
                if(list[poly].num_bond[i] == bond_num and list[poly].material[i] == mat and list[poly].indic_action[i] == act2){
                    test = 1;
                    break;
                }
/*                if (test == 1 and  *list[poly].prop_uninhib[i] != prope){
                    if(crystalline_bool == true){
                        list[poly].liste_prop[i] = &par.crystal_propensities[act2-1];
                        list[poly].prop_uninhib[i] = &par.crystal_propensities[act2-1];
                    }
                    else{
                        list[poly].liste_prop[i] = &par.propensities[act2-1];
                        list[poly].prop_uninhib[i] = &par.propensities[act2-1];
                    }

                }*/
            }
            if(test == 0){
                list[poly].num_bond.push_back(bond_num);
                list[poly].material.push_back(mat);
                list[poly].indic_action.push_back(act2);
                if(crystalline_bool == true){
                    list[poly].liste_prop.push_back(&par.crystal_propensities[act2-1]);
                    list[poly].prop_uninhib.push_back(&par.crystal_propensities[act2-1]);
                }
                else{
                    list[poly].liste_prop.push_back(&par.propensities[act2-1]);
                    list[poly].prop_uninhib.push_back(&par.propensities[act2-1]);
                }
                list[poly].crystalline.push_back(crystalline_bool);
                list[poly].covered.push_back(false);
                list[poly].addProp(prope);
                list[poly].nbr_element++;
//                cout << "added; liste[poly].nbr_element = " << list[poly].nbr_element << endl;

            }
            else{
                return;
            }

        }
        else{
            list[poly].num_bond.push_back(bond_num);
            list[poly].material.push_back(mat);
            list[poly].indic_action.push_back(act2);
            if(crystalline_bool == true){
                list[poly].liste_prop.push_back(&par.crystal_propensities[act2-1]);
                list[poly].prop_uninhib.push_back(&par.crystal_propensities[act2-1]);
            }
            else{
                list[poly].liste_prop.push_back(&par.propensities[act2-1]);
                list[poly].prop_uninhib.push_back(&par.propensities[act2-1]);
            }
            list[poly].crystalline.push_back(crystalline_bool);
            list[poly].covered.push_back(false);
            list[poly].addProp(prope);
            list[poly].nbr_element++;
//            cout << "added; liste[poly].nbr_element = " << list[poly].nbr_element << endl;
        }
    }


//End of function
}


//Clears the table with index table from all reactions
void clear_table(params& par, std::vector<TList>& allTables, int table){
    allTables[table].num_bond.clear();
    allTables[table].material.clear();
    allTables[table].indic_action.clear();
    allTables[table].liste_prop.clear();
    allTables[table].prop_uninhib.clear();
    allTables[table].crystalline.clear();
    allTables[table].covered.clear();
    allTables[table].nbr_element = 0;
    allTables[table].prop_sum = 0;
}


//Fills the reaction table of a newly created polymer
void fill_table(params& par, vector<TList>& allTables, vector<bList>& polyList, int poly_selected, int substrate, bool& error, vector<double>& chem_entities, int& nbr_Glc_pdt, int& nbr_cellobiose){
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
            propensite = prop(par, reaction_type, chem_entities,nbr_Glc_pdt,nbr_cellobiose);
            crystal_propensite = par.crystal_modifier_cellu*propensite;

            if(propensite > 0 and polyList[poly_selected].status[0] == 1){
                if(polyList[poly_selected].crystalline[0] == false)
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,propensite,0);
                else
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,crystal_propensite,1);
            }
        }
        else if(polyList[poly_selected].len_poly > 1){
            reaction_type = 1;
            propensite = prop(par, reaction_type, chem_entities,nbr_Glc_pdt,nbr_cellobiose);
            crystal_propensite = par.crystal_modifier_cellu*propensite;
//            cout << propensite << "\t" << crystal_propensite << endl;

            if(propensite > 0){
                if(polyList[poly_selected].len_poly == 2){//Digestion of this polymer length by EG is currently disabled

                }
                else if(polyList[poly_selected].len_poly == 3){

                }
                else if(polyList[poly_selected].len_poly == 4){
                    /*
                    if(polyList[poly_selected].status[1] == 1){
                            if(polyList[poly_selected].crystalline[1] == false)
                                addreaction(par, false, allTables,poly_selected,1,substrate,reaction_type,propensite,0);
                            else
                                addreaction(par, false, allTables,poly_selected,1,substrate,reaction_type,crystal_propensite,1);
                    }
                    if(polyList[poly_selected].status[2] == 1){
                            if(polyList[poly_selected].crystalline[2] == false)
                                addreaction(par, false, allTables,poly_selected,2,substrate,reaction_type,propensite,0);
                            else
                                addreaction(par, false, allTables,poly_selected,2,substrate,reaction_type,crystal_propensite,1);
                    }
                    */
                }
                else{
                    for(int i = 2;i<polyList[poly_selected].len_poly-2;i++){
                        if(polyList[poly_selected].status[i] == 1){
                            if(polyList[poly_selected].crystalline[i] == false)
                                addreaction(par, false, allTables,poly_selected,i,substrate,reaction_type,propensite,0);
                            else
                                addreaction(par, false, allTables,poly_selected,i,substrate,reaction_type,crystal_propensite,1);
                        }
                    }
                }
            }
            //CBH REACTION ADDING IS DONE IN FUNCTION CBH_ATTACHMENT()
            /*reaction_type = 2;
            propensite = prop(par, reaction_type, chem_entities,nbr_Glc_pdt,nbr_cellobiose);
            crystal_propensite = par.crystal_modifier_cellu*propensite;
            propensite = 0;
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
                    if(polyList[poly_selected].crystalline[leftBond] == false){
                        addreaction(par, false, allTables,poly_selected,leftBond,substrate,reaction_type,propensite,enzyme_index);

                    }
                    else{
                        addreaction(par, false, allTables,poly_selected,leftBond,substrate,reaction_type,crystal_propensite,enzyme_index);

                    }
                }
                if(polyList[poly_selected].status[rightBond] == 1 and par.N_free_CBH > 0){
                    if(polyList[poly_selected].crystalline[rightBond] == false){
                        addreaction(par, false, allTables,poly_selected,rightBond,substrate,reaction_type,propensite,enzyme_index);

                    }
                    else{
                        addreaction(par, false, allTables,poly_selected,rightBond,substrate,reaction_type,crystal_propensite,enzyme_index);

                    }
                }
            }*/
        }

        for(int i=0;i<allTables[poly_selected].nbr_element;i++){
            if(polyList[poly_selected].len_poly != 2 and (allTables[poly_selected].num_bond[i] == 0 or allTables[poly_selected].num_bond[i] == polyList[poly_selected].len_poly-1) and (allTables[poly_selected].indic_action[i] == 1 or allTables[poly_selected].indic_action[i] == 2)){
                cout << "Filling of cellu tables is making problems: " << endl;
                exit(1);
                return;
            }
        }
    }
    else if(substrate == 2){    // changed bond_type in hemi to 1: partho
        reaction_type = 4;
        propensite = prop(par, reaction_type, chem_entities,nbr_Glc_pdt,nbr_cellobiose);
        crystal_propensite = par.crystal_modifier_hemi*propensite;

        for(int i=0;i<polyList[poly_selected].len_poly;i++){

            if (par.xyl_or_mlg == true){ // partho added check for MLG or XYL
                if(polyList[poly_selected].bond_type[i] == 4 and propensite > 0 and polyList[poly_selected].status[i] == 1){    //bond type = 1 for mlg, 4 for xyl partho change
                    if(polyList[poly_selected].crystalline[i] == false)
                        addreaction(par, false,allTables,poly_selected,i,substrate,reaction_type,propensite,0);
                    else
                        addreaction(par, false,allTables,poly_selected,i,substrate,reaction_type,crystal_propensite,1);
                }                
            }
            else if (par.xyl_or_mlg == false){  // partho added check for MLG or XYL
                if(polyList[poly_selected].bond_type[i] == 1 and propensite > 0 and polyList[poly_selected].status[i] == 1){    //bond type = 1 for mlg, 4 for xyl partho change
                    if(polyList[poly_selected].crystalline[i] == false)
                        addreaction(par, false,allTables,poly_selected,i,substrate,reaction_type,propensite,0);
                    else
                        addreaction(par, false,allTables,poly_selected,i,substrate,reaction_type,crystal_propensite,1);
                }
            }
        }
    }
//End of function
}

//Deletes the reaction labelled "action_target" for the polymer "poly" for its bond "bond_target" of the material "mat" in the list of reactions
void deletereaction(TList& list,vector<TList>& Table,int mat,int poly,int bond_target,int action_target)
{

    for(int i=0; i<Table[poly].nbr_element; i++) {

        if((Table[poly].num_bond[i]==bond_target) and (Table[poly].indic_action[i]==action_target) and (Table[poly].material[i]==mat)) {// So basically if the reaction fits all the function arguments

            Table[poly].num_bond.erase(Table[poly].num_bond.begin() + i);
            Table[poly].material.erase(Table[poly].material.begin() + i);
            Table[poly].indic_action.erase(Table[poly].indic_action.begin() + i);

            Table[poly].prop_sum -= *Table[poly].liste_prop[i];

            Table[poly].liste_prop.erase(Table[poly].liste_prop.begin() + i);
            Table[poly].prop_uninhib.erase(Table[poly].prop_uninhib.begin() + i);
            Table[poly].crystalline.erase(Table[poly].crystalline.begin() + i);
            Table[poly].covered.erase(Table[poly].covered.begin() + i);

            Table[poly].nbr_element--;
        }
    }
    if(Table[poly].liste_prop.size() == 0)
        Table[poly].prop_sum = 0;
//End of function
}


void delete_specific_reaction(vector<TList>& Table, int table_selected, int i){
    if(Table[table_selected].nbr_element > 0){
//        cout << "deletin; i = " << i << "; size of table: " << Table[table_selected].nbr_element << endl;
        Table[table_selected].num_bond.erase(Table[table_selected].num_bond.begin() + i);
        Table[table_selected].material.erase(Table[table_selected].material.begin() + i);
        Table[table_selected].indic_action.erase(Table[table_selected].indic_action.begin() + i);

        Table[table_selected].prop_sum -= *Table[table_selected].liste_prop[i];


        Table[table_selected].liste_prop.erase(Table[table_selected].liste_prop.begin() + i);
        Table[table_selected].prop_uninhib.erase(Table[table_selected].prop_uninhib.begin() + i);
        Table[table_selected].crystalline.erase(Table[table_selected].crystalline.begin() + i);
        Table[table_selected].covered.erase(Table[table_selected].covered.begin() + i);
        Table[table_selected].nbr_element--;
    }

    if(Table[table_selected].liste_prop.size() == 0)
        Table[table_selected].prop_sum = 0;
//End of function
}

//Deletes ALL reactions for the bond "bond_target" from the polymer "poly" of the material "mat" in the list of reactions
void deleteAllreaction(TList& list,vector<TList>& Table,int mat,int poly,int bond_target)
{
    for(int i=0; i<Table[poly].nbr_element; i++) {

        if((Table[poly].num_bond[i]==bond_target) and (Table[poly].material[i]==mat)) {// So basically if the reaction fits all the function arguments

            Table[poly].num_bond.erase(Table[poly].num_bond.begin() + i);
            Table[poly].material.erase(Table[poly].material.begin() + i);
            Table[poly].indic_action.erase(Table[poly].indic_action.begin() + i);

            Table[poly].prop_sum -= *Table[poly].liste_prop[i];

//cout << list.liste_prop.size() << "\t" << list.prop_uninhib.size() << endl;
            Table[poly].liste_prop.erase(Table[poly].liste_prop.begin() + i);
            Table[poly].prop_uninhib.erase(Table[poly].prop_uninhib.begin() + i);
            Table[poly].crystalline.erase(Table[poly].crystalline.begin() + i);
            Table[poly].covered.erase(Table[poly].covered.begin() + i);
//    cout << "searching in functions.cpp" << endl;
            Table[poly].nbr_element--;

            i -= 1;//To be sure we scan the right hand side neighbor in the list that just shift from 1 rank
        }
    }
    if(Table[poly].liste_prop.size() == 0)
        Table[poly].prop_sum = 0;

//End of function
}


//Computes the propensity for each reaction to take place
double prop(params& par, double act, vector<double> chem_entities, int& nbr_Glc_pdt, int& nbr_cellobiose)//For now hemicellulose and cellulose are digested with same rates, but could easily change
{	// partho changes in propensity for inhibition
    //calculating inhib values for EG, CBH and BGL
//    double inhib_nr_EG = par.inhib_cellobiose_EG * chem_entities[0]* nbr_cellobiose / (chem_entities[0] + chem_entities[1] + nbr_cellobiose); //old inhibition	//EG inhibit by cellobiose
//    double inhib_nr_CBH = par.inhib_cellobiose_CBH * chem_entities[1]* nbr_cellobiose / (chem_entities[0] + chem_entities[1] + nbr_cellobiose); //old inhibition	//CBH inhibit by cellobiose
//    double inhib_nr_BGL =  par.inhib_glucose_BGL * chem_entities[2] * nbr_Glc_pdt / (chem_entities[2] + nbr_Glc_pdt);   //old inhibition				//BGL inhibit by glucose
  //No inhibition of XYL
    double inhib_nr_EG = par.inhib_cellobiose_EG * chem_entities[0] * nbr_cellobiose / (chem_entities[0] + chem_entities[1] + nbr_cellobiose)    +   par.inhib_glucose_EG * chem_entities[0]* nbr_Glc_pdt / (chem_entities[0] + chem_entities[1] + chem_entities[2] + nbr_Glc_pdt);      //partho: EG inhbit by glucose & cellobiose
    double inhib_nr_CBH = par.inhib_cellobiose_CBH * chem_entities[1] * nbr_cellobiose / (chem_entities[0] + chem_entities[1] + nbr_cellobiose)  +   par.inhib_glucose_CBH * chem_entities[1]*nbr_Glc_pdt / (chem_entities[0] + chem_entities[1] + chem_entities[2] + nbr_Glc_pdt);      //partho: CBH inhbit by glucose & cellobiose
    double inhib_nr_BGL =  par.inhib_glucose_BGL * chem_entities[2] * nbr_Glc_pdt / (chem_entities[0] + chem_entities[1] + chem_entities[2] + nbr_Glc_pdt);       //partho: BGL only inhibited by glucose
  //No Inhibiton of XYL

/*
//  Inhibtion of XYL by Glc
    double inhib_nr_EG = par.inhib_cellobiose_EG * chem_entities[0] * nbr_cellobiose / (chem_entities[0] + chem_entities[1] + nbr_cellobiose)    +   par.inhib_glucose_EG * chem_entities[0]* nbr_Glc_pdt / (chem_entities[0] + chem_entities[1] + chem_entities[2] + chem_entities[3] + nbr_Glc_pdt);      //partho: EG inhbit by glucose & cellobiose
    double inhib_nr_CBH = par.inhib_cellobiose_CBH * chem_entities[1] * nbr_cellobiose / (chem_entities[0] + chem_entities[1] + nbr_cellobiose)  +   par.inhib_glucose_CBH * chem_entities[1]*nbr_Glc_pdt / (chem_entities[0] + chem_entities[1] + chem_entities[2] + chem_entities[3] + nbr_Glc_pdt);      //partho: CBH inhbit by glucose & cellobiose
    double inhib_nr_BGL =  par.inhib_glucose_BGL * chem_entities[2] * nbr_Glc_pdt / (chem_entities[0] + chem_entities[1] + chem_entities[2] + chem_entities[3] + nbr_Glc_pdt);       //partho: BGL only inhibited by glucose
    double inhib_nr_XYL = par.inhib_XYL * chem_entities[3] * nbr_Glc_pdt / (chem_entities[0] + chem_entities[1] + chem_entities[2] + chem_entities[3] + nbr_Glc_pdt);       //partho: XYL inhbit by glucose for MLG
*/


//    cout <<"------------------------------------------------------------------------"<<endl;
//    cout << "partho--->"<< "inhib_nr_EG="<< inhib_nr_EG << "--inhib_nr_CBH=" << inhib_nr_CBH <<"--inhib_nr_BGL="<< inhib_nr_BGL<< "---Number of CelloBiose="<< nbr_cellobiose<< "---Number of Glucose="<< nbr_Glc_pdt<<endl;
//    cout <<"========================================================================"<<endl;


    double propens;

    if (par.mode_inhib == -1){
           if(act == 1)//Degradation by EG
               propens=par.k1*chem_entities.at(0);
           else if(act == 2){//Degradation by CBH
               propens = par.k2;//Processive enzyme --> does not care about concentration once it is attached
               //propens=par.k2*chem_entities.at(1);
           }
           else if(act == 3)//Degradation by BGL
               propens = par.k3*chem_entities.at(2);
           else if(act==4)//Degradation of hemi by XYL
               propens = par.k4*chem_entities.at(3);
           else if(act == 5)//Degradation of hemi by XYL
           if (par.mode_lignin_glue == 1){
               propens = par.k5*chem_entities.at(4)*(chem_entities[0]+chem_entities[1]+chem_entities[2]+chem_entities[3]);
           }
           else{
               propens = 0;
           }
           else if(act == 6){
           //        propens = par.k6*par.N_free_CBH;
               if(par.N_free_ends == 0){
                   propens = 0.;
               }
               else{
                   propens = par.k6*par.N_free_CBH*(1.- double((par.CBH_enzymes.size()-par.N_free_CBH)*par.ends_blocked_per_CBH)/double(par.N_free_ends));//*(1.-tanh(par.ends_blocked_per_CBH*(par.CBH_enzymes.size()-par.N_free_CBH)));
                   if (propens < 0){
                       propens = 0.;
                   }
               }
           }
           else{
               cout << "In function prop(): variable 'act' does not correspond to any of the enzymes included in this simulation; act = " << act << endl;
           return 0;
           }
           if (par.verbose==true){
            //    cout <<"partho:--     action:-->     "<< act<< "      ----propens:-->  "<<propens<< endl;
           }

    }
    else if (par.mode_inhib == 1 and nbr_Glc_pdt >=0 and nbr_cellobiose>=0){
            if(act == 1){//Degradation by EG
            propens=par.k1*(chem_entities.at(0) - inhib_nr_EG); // reduced EG propensity partho
//                cout <<"partho :----  action:-->    "<< act<<"--   propens=  "<<par.k1*chem_entities.at(0) << " ---propens_inhib="<< propens<< endl;
            }
            else if(act == 2){//Degradation by CBH
                propens = par.k2;//Processive enzyme --> does not care about concentration once it is attached
                //propens=par.k2*chem_entities.at(1);
            }
            else if(act == 3){//Degradation by BGL
                propens = par.k3*(chem_entities.at(2) - inhib_nr_BGL);
//                cout <<"partho :----  action:-->    "<< act<<"--   propens=  "<<par.k3*chem_entities.at(2) << " ---propens_inhib="<< propens<< endl;
            } // reduced BGL propensity partho

            else if(act==4){//Degradation of hemi by XYL
 //               propens = par.k4*(chem_entities.at(3) - inhib_nr_XYL);
                propens = par.k4*chem_entities.at(3);
//                cout <<"partho :----  action:-->    "<< act<<"--   propens=  "<<par.k4*chem_entities.at(3) << " ---propens_inhib="<< propens<< endl;
            }
            else if(act == 5)//Lignin Glueing
            if (par.mode_lignin_glue == 1){
                propens = par.k5*chem_entities.at(4)*(chem_entities[0]+chem_entities[1]+chem_entities[2]+chem_entities[3]);
            }
            else{
                propens = 0;
            }
            else if(act == 6){
            //        propens = par.k6*par.N_free_CBH;
                if(par.N_free_ends == 0){
                    propens = 0.;
                }
                else{
                    propens = par.k6*(par.N_free_CBH - inhib_nr_CBH)*(1.- double((par.CBH_enzymes.size()-par.N_free_CBH)*par.ends_blocked_per_CBH)/double(par.N_free_ends));//*(1.-tanh(par.ends_blocked_per_CBH*(par.CBH_enzymes.size()-par.N_free_CBH)));  //reduced CBH attachment propensity partho
                    if (propens < 0){
                        propens = 0.;
                    }
//                    cout <<"partho :----  action:-->    "<< act<<"--   propens=  "<<par.k6*(par.N_free_CBH)*(1.- double((par.CBH_enzymes.size()-par.N_free_CBH)*par.ends_blocked_per_CBH)/double(par.N_free_ends)) << " ---propens_inhib="<< propens<< endl;
                }
            }
            else{
                cout << "In function prop(): variable 'act' does not correspond to any of the enzymes included in this simulation; act = " << act << endl;
            return 0;
            }
            if (par.verbose==true){
            //    cout <<"partho:--     action:-->    "<< act<< "      ----propens:-->  "<<propens<< endl;
            }
    }


//    cout << "In function prop(): act = " << act << "; propens = " << propens << endl;
//    cout<< endl;
//    cout <<"partho:--       Inhibitionmode:-->  "<<par.mode_inhib<<endl;
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
                if(list[i].bond_type[0] == 1 or list[i].bond_type[0] == 3 or list[i].bond_type[0] == 5) // Partho new bond
                    count++;
                for(int j=0;j < list[i].bond_type.size();j++){
                    if(list[i].bond_type[j] == 1 or list[i].bond_type[j] == 2 or list[i].bond_type[j] == 5){   //Partho new bond
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


void pick_reaction(params& par,
    const std::vector<TList>& Table_cellu,
    const std::vector<TList>& Table_hemi,
    const TList& Table_lign,
    const std::vector<bList>& cellu,
    const std::vector<bList>& hemi,
    const std::vector<bList>& lign,
    const double a0,
    double& test,
    double& tau,
    int& mu1,
    int& table_selected,
    int& poly_selected,
    int& bond_selected,
    int& action_mu1,
    int& substrate)
{


    double r1 = drand48();//Generates random numbers between 0 and 1
    double r2 = drand48();

    tau = 1/a0*(log(1/r1));//Chooses time interval of the next reaction

    test = r2*a0; //Chooses next reaction to happen
    int j = 0;
    double suma = Table_lign.prop_sum;
    if(test <= suma){
        substrate = 0;
        action_mu1 = Table_lign.indic_action[0];
        mu1 = 0;
        poly_selected=Table_lign.index_poly;//Corresponding polymer
        bond_selected=Table_lign.num_bond[0];//Returns the index = position of the bond in the Glc list of thepoly it belongs to
        //substrate=Table_cellu[table_selected].material[mu1];//Corresponding material: 1= cellulose; 2= hemicellulose
    }
    else{
        if(Table_cellu.size() > 0){
            substrate = 1;
            suma += Table_cellu[j].prop_sum;

            while(test>suma)
            {
                j++;
                if(j==Table_cellu.size())
                    break;
                if(Table_cellu[j].liste_prop.size()>0)
                    suma += Table_cellu[j].prop_sum;
            }
            if(j==Table_cellu.size()){
                substrate = 2;
                j=0;
            }
            else
                substrate = 1;
            if(substrate == 2 and Table_hemi.size()>0){
                j = 0;
                suma +=Table_hemi[0].prop_sum;
                while(test>suma){
                    j++;
                    if(Table_hemi[j].liste_prop.size()>0)
                        suma += Table_hemi[j].prop_sum;
                }

            }

            table_selected = j;

        }
        else if(Table_hemi.size() > 0){
            substrate = 2;
            if(Table_hemi[j].liste_prop.size()>0)
                suma += Table_hemi[j].prop_sum;
            else
                suma = Table_lign.prop_sum;

            while(test > suma){
                j++;
                if(Table_hemi[j].liste_prop.size()>0)
                    suma += Table_hemi[j].prop_sum;
            }
            table_selected = j;
        }
        else {
            cout << "Reaction table empty. Stopping." << endl;
            exit(1);
        }
        if(par.verbose == true)
            cout << "table_selected: " << table_selected << endl;
        //============ Choose reaction within chosen table ========

        if(substrate == 1){
            if(par.verbose == true){
                cout << "substrate = 1" << endl;
            }
            mu1 = findIndex(Table_cellu[table_selected], suma, test);//Index of the reaction chosen
            if(mu1 == -1){
                suma = -1.;
                mu1 = -1;
                table_selected = -1;
                poly_selected = -1;
                bond_selected = -1;
                action_mu1 = -1;
                substrate = -1;
                tau = -1.;
                return;

            }
/*                else if(mu1 == Table_cellu[table_selected].nbr_element-1){
                N_last++;
            }*/
            poly_selected=Table_cellu[table_selected].index_poly;//Corresponding polymer
            bond_selected=Table_cellu[table_selected].num_bond[mu1];//Returns the index = position of the bond in the Glc list of thepoly it belongs to
            //substrate=Table_cellu[table_selected].material[mu1];//Corresponding material: 1= cellulose; 2= hemicellulose
            action_mu1=Table_cellu[table_selected].indic_action[mu1];//Corresponding reaction
        }
        else if(substrate == 2){
            if(par.verbose == true){
                cout << "substrate = 2" << endl;
                cout << "table_selected = " << table_selected << "; size of Table_hemi = " << Table_hemi.size() << "; size of Table_cellu: " << Table_cellu.size() << endl;
            }
            if(table_selected >= Table_hemi.size()){
                cout << "table_selected = " << table_selected << "; Table_hemi.size() = " << Table_hemi.size() << endl;
                cout << "table_selected >= Table_hemi.size()! Stopping" << endl;
                exit(1);
            }
            if(par.verbose == true){
                for(int i=0;i<Table_hemi[table_selected].nbr_element;i++){
                    cout << i << endl;
                    if(Table_hemi[table_selected].num_bond[i] >= hemi[table_selected].len_poly){
                        cout << "Table_hemi[table_selected].num_bond[i] >= hemi[table_selected].len_poly" << endl;
                        cout << "Table_hemi[table_selected].num_bond[i] = " << Table_hemi[table_selected].num_bond[i] << "; hemi[table_selected].len_poly = " << hemi[table_selected].len_poly << endl;
                        exit(1);
                    }
                }
                cout << "Through for" << endl;
            }
            mu1 = findIndex(Table_hemi[table_selected], suma, test);//Index of the reaction chosen
            if(par.verbose == true){
                cout << "found index" << endl;
            }
            if(mu1 == -1){
                suma = -1.;
                mu1 = -1;
                table_selected = -1;
                poly_selected = -1;
                bond_selected = -1;
                action_mu1 = -1;
                substrate = -1;
                tau = -1.;
                return;

            }
/*                else if(mu1 == Table_hemi[table_selected].nbr_element-1){
                N_last++;
            }*/
            poly_selected=Table_hemi[table_selected].index_poly;//Corresponding polymer
            bond_selected=Table_hemi[table_selected].num_bond[mu1];//Returns the index = position of the bond in the Glc list of thepoly it belongs to
            //substrate=Table_hemi[table_selected].material[mu1];//Corresponding material: 1= cellulose; 2= hemicellulose
            action_mu1=Table_hemi[table_selected].indic_action[mu1];//Corresponding reaction

        }
        else if(substrate != 1 and substrate != 2)
            cout << "Something is going wrong in the selection of the reaction. substrate is neither cellu nor hemi!" << endl;

    }


}



int findIndex(const TList& list, double suma, double test){

    suma -= list.prop_sum;

    int j = 0;

    if(list.liste_prop.size() == 0){
        cout << "In function findIndex(): table contains no reactions!" << endl;
        return -1;
    }
    if(list.covered[j] == false){
        suma += *list.liste_prop[j];
    }
    while(test > suma or *list.liste_prop[j] <= 0){
        j++;
//        cout << *list.liste_prop[j] << endl;
        if(j>= list.liste_prop.size()){
            j--;
            break;
        }
        if(list.covered[j] == false){
            suma += *list.liste_prop[j];
        }
    }
    if(*list.liste_prop[j] > 0 and list.covered[j] == false){
        return j;
    }
    else{
        return -1;
//        cout << "Selected a reaction with 0 propensity... j = " << j << "; nbr_element = " << list.nbr_element << endl;
//        cout << "Stopping" << endl;
//        exit(1);
    }
//End of function
}


void EG_digest(params& par, vector<TList>& allTables, vector<bList>& polyList, int& nbr_poly, int& nbr_xyl_pdt, int& nbr_Glc_pdt,const int bond_selected,const int poly_selected, int len_polyLoopStart, int& nbr_cellobiose, vector<double>& chem_entities, const int substrate,unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu, unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_hemi, const bool verbose, bool& error){//EG digestion as single function

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
    std::tuple<double,double,int> neighbor_key = std::make_tuple(x_pos, y_pos, z_pos);
    int test1; //Used for finding specific reactions in a table
    int indic_cut;//
    int new_poly_index = -1;
    int minSize_BGL = 8;//DP at which BGL starts digesting
    int x_bond,y_bond,z_bond;// For later freeing of lower bonds
    int enzyme_index = -1;

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
        cout << "EG digestion on a polymer of length 1. This should not occur. Stopping." << endl;
        exit(1);
        //Do nothing
    }
    else if(polyList[poly_selected].len_poly==2){
        cout << "EG digestion of a polymer of length 2. This should not currently occur. Stopping" << endl;
        exit(1);
/*
        if(bond_selected==0)// If the cut was at the beginning of the polymer
        {
            if(par.free_poly_ends.find(cantor_pair_two(poly_selected,bond_selected)) != par.free_poly_ends.end()){
                par.free_poly_ends.erase(cantor_pair_two(poly_selected,bond_selected));
            }
            taylorNewPoly(polyList[poly_selected],bond_selected);//Shift the indices by one to the left
            clear_table(par,allTables,poly_selected);


            reaction_type = 3;
            propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
            crystal_propensite = par.crystal_modifier_cellu*propensite;
            if(propensite > 0 and polyList[poly_selected].status[0] == 1 and polyList[poly_selected].len_poly == 1){

                addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,propensite,enzyme_index);
            }


            nbr_Glc_pdt++;//We release a Glc
        }
        else if(bond_selected == 1){
            if(par.free_poly_ends.find(cantor_pair_two(poly_selected,bond_selected)) != par.free_poly_ends.end()){
                par.free_poly_ends.erase(cantor_pair_two(poly_selected,bond_selected));
            }
            taylorOldPoly(polyList[poly_selected],bond_selected);//Erase the end of the polymer, until digested bond
            clear_table(par,allTables,poly_selected);

            reaction_type = 3;
            propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
            crystal_propensite = par.crystal_modifier_cellu*propensite;
            if(propensite > 0 and polyList[poly_selected].status[0] == 1 and polyList[poly_selected].len_poly == 1){
                if(polyList[poly_selected].crystalline[0]==false)
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,propensite,enzyme_index);
                else
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,crystal_propensite,enzyme_index);
            }

            nbr_Glc_pdt++;//We release a Glc
        }
*/
    }
    else if(polyList[poly_selected].len_poly == 3){
        cout << "EG digestion of a polymer of length 3. This should not currently occur. Stopping" << endl;
        exit(1);
    }
    else if(polyList[poly_selected].len_poly>3){
        if((bond_selected !=0) and (bond_selected !=polyList[poly_selected].len_poly-1)){

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


            //Check, if the end of the polymer is occupied by a CBH enzyme
            bool insert_new_end = false;
            if(par.free_poly_ends.find(cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2)) != par.free_poly_ends.end()){//Check if we need to move a free polymer end from the old polymer to the new one
                par.free_poly_ends.erase(cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2));//This needs to be re-added later
                par.N_free_ends--;
                insert_new_end = true;
            }

            taylorOldPoly(polyList[poly_selected],bond_selected);//Erase the end of the old polymer, until digested bond
            taylorNewPoly(polyList[new_poly_index],bond_selected);//Erase the begining of the new polymer, from digested bond
            indic_cut = 0;//In this function, the old poly is always the one on the left, and the new one on the right

            clear_table(par,allTables,poly_selected);
            fill_table(par,allTables,polyList,poly_selected,substrate,error,chem_entities,nbr_Glc_pdt,nbr_cellobiose);

            //Fill new table
            fill_table(par,allTables, polyList, new_poly_index, substrate, error,chem_entities,nbr_Glc_pdt,nbr_cellobiose);

            if(insert_new_end == true and polyList[new_poly_index].len_poly > 1){//We do not need to check whether the status of the bond is equal to 1, because we already know that a free poly end was found there
                if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,polyList[new_poly_index].len_poly-2)) == par.free_poly_ends.end()){
                    if(polyList[new_poly_index].len_poly-2 >= 0 and CBH_enzyme_attached(par,new_poly_index,polyList[new_poly_index].len_poly-2 == false)){
                        par.free_poly_ends.insert({cantor_pair_two(new_poly_index,polyList[new_poly_index].len_poly-2),std::make_tuple(new_poly_index,polyList[new_poly_index].len_poly-2)});
                        if(true){
                            addreaction(par, false,allTables,new_poly_index,polyList[new_poly_index].len_poly-2,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[polyList[new_poly_index].len_poly-2]);
                        }
                        par.N_free_ends++;
                    }
                }
            }

            //insert newly free ends at cut location
            if(polyList[poly_selected].len_poly > 1){
                if(par.free_poly_ends.find(cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2)) == par.free_poly_ends.end()){
                    if(polyList[poly_selected].len_poly-2 >= 0 and CBH_enzyme_attached(par, poly_selected, polyList[poly_selected].len_poly-2) == false){
                        par.free_poly_ends.insert({cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2),make_tuple(poly_selected,polyList[poly_selected].len_poly-2)});
                        if(true){
                            addreaction(par, false,allTables,poly_selected,polyList[poly_selected].len_poly-2,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[polyList[poly_selected].len_poly-2]);
                        }
                        par.N_free_ends++;
                    }
                }
            }
            if(newPolyFlag == true){
                if(polyList[new_poly_index].len_poly > 1){
                    if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par,new_poly_index,1) == false){
                        par.free_poly_ends.insert({cantor_pair_two(new_poly_index,1),make_tuple(new_poly_index,1)});
                        if(true){
                            addreaction(par, false,allTables,new_poly_index,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[1]);
                        }
                        par.N_free_ends++;
                    }
                }

            }


        //UPDATE CBH ENZYMES ATTACHED
        int CBH_found = 0;//max two CBH enzymes can be attached to one polymer --> this is a cut-off condition for the following for loop
        for(int i=0; i<par.CBH_enzymes.size();i++){
            if(par.CBH_enzymes[i].attached == true and par.CBH_enzymes[i].poly_attached == poly_selected and par.CBH_enzymes[i].bond_attached == 1){
                if(polyList[poly_selected].len_poly != 1){
                    CBH_found++;
                    addreaction(par, false, allTables,poly_selected,1,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[1]);
                }
                else{
                    par.CBH_enzymes[i].detach();
                    par.N_free_CBH++;
                }


            }
            else if(par.CBH_enzymes[i].attached == true and par.CBH_enzymes[i].poly_attached == poly_selected and (par.CBH_enzymes[i].bond_attached != 1)){
                if(polyList[new_poly_index].len_poly != 1){
                    CBH_found++;
                    par.CBH_enzymes[i].poly_attached = new_poly_index;
                    par.CBH_enzymes[i].bond_attached = polyList[new_poly_index].len_poly-2;
                    if(par.CBH_enzymes[i].attached == true and par.CBH_enzymes[i].bond_attached == 0 and polyList[new_poly_index].len_poly !=2){
                        cout << "In function EG_digest: bond_attached = 0, but len_poly != 2. Stopping" << endl;
                        exit(1);
                    }
                    else if(par.CBH_enzymes[i].attached == true and polyList[new_poly_index].len_poly == 1){
                        par.CBH_enzymes[i].detach();
                        par.N_free_CBH++;
                    }
                    else{
                        if(par.CBH_enzymes[i].attached == true and par.CBH_enzymes[i].bond_attached == 1 and polyList[new_poly_index].len_poly == 3){
                            int count_attached = 0;
                            for(int j=0;j<par.CBH_enzymes.size();j++){
                                if(par.CBH_enzymes[j].attached == true and par.CBH_enzymes[j].poly_attached == new_poly_index and par.CBH_enzymes[j].bond_attached == 1){
                                    count_attached++;
                                }
                            }
                            if(count_attached > 1){
                                par.CBH_enzymes[i].detach();
                                par.N_free_CBH++;
                            }
                            else{
                                addreaction(par, false, allTables,new_poly_index,polyList[new_poly_index].len_poly-2,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[polyList[new_poly_index].len_poly-2]);
                            }
                            //There can only be one poly attached to bond 1
                        }
                        else{
                            addreaction(par, false, allTables,new_poly_index,polyList[new_poly_index].len_poly-2,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[polyList[new_poly_index].len_poly-2]);
                        }
                    }
                }
                else{
                    par.CBH_enzymes[i].detach();
                    par.N_free_CBH++;
                }
            }
            if(CBH_found > 2){
                cout << "In function EG_digest(): CBH_found > 2" << endl;
                exit(1);
            }
        }


        }

    }

    //REACTION TABLE MODIFICATION IS CARRIED OUT IN SEPARATE FUNCTION, EXCEPT FOR CBH ATTACHMENT AND DIGESTION REACTIONS


    if(polyList[poly_selected].len_poly == 1){
        for(int i=0;i<par.CBH_enzymes.size();i++){
            if(par.CBH_enzymes[i].attached == true and par.CBH_enzymes[i].poly_attached == poly_selected){
                par.CBH_enzymes[i].detach();
                par.N_free_CBH++;
            }

        }
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,0)) != par.free_poly_ends.end()){
            par.free_poly_ends.erase(cantor_pair_two(poly_selected,0));
            deletereaction(allTables[poly_selected], allTables, 1, poly_selected, 0, 6);//Delete the reaction that will now be carried out
            par.N_free_ends--;
        }
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) != par.free_poly_ends.end()){
            par.free_poly_ends.erase(cantor_pair_two(poly_selected,1));
            deletereaction(allTables[poly_selected], allTables, 1, poly_selected, 1, 6);//Delete the reaction that will now be carried out
            par.N_free_ends--;
        }
        addreaction(par, false, allTables,poly_selected,0,1,3,prop(par,3,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[0]);
        nbr_cellobiose++;
    }
    if(newPolyFlag == true){

        if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,0)) != par.free_poly_ends.end()){
            par.free_poly_ends.erase(cantor_pair_two(new_poly_index,0));
            deletereaction(allTables[new_poly_index], allTables, 1, new_poly_index, 0, 6);//Delete the reaction that will now be carried out
            par.N_free_ends--;
        }
        if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,1)) != par.free_poly_ends.end()){
            par.free_poly_ends.erase(cantor_pair_two(new_poly_index,1));
            deletereaction(allTables[new_poly_index], allTables, 1, new_poly_index, 1, 6);//Delete the reaction that will now be carried out
            par.N_free_ends--;
        }
        if(polyList[new_poly_index].len_poly == 1){
            for(int i=0;i<par.CBH_enzymes.size();i++){
                if(par.CBH_enzymes[i].attached == true and par.CBH_enzymes[i].poly_attached == new_poly_index){
                    par.CBH_enzymes[i].detach();
                    par.N_free_CBH++;
                }

            }
            addreaction(par, false, allTables,new_poly_index,0,1,3,prop(par,3,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[0]);
            nbr_cellobiose++;
        }
    }

    //Check if either polymer has length 3 or 2, or is missing reactions
    if(polyList[poly_selected].len_poly > 3){
        //First end
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, poly_selected,1) == false){
            //cout << "In function EG_digest(): adding attachment reaction for length 3 polymer (old poly)" << endl;
            par.free_poly_ends.insert({cantor_pair_two(poly_selected,1),make_tuple(poly_selected,1)});
            if(true){
                addreaction(par, false,allTables,poly_selected,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[1]);
            }
            par.N_free_ends++;
        }
        else if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, poly_selected,1) == true){
            addreaction(par, false,allTables,poly_selected,1,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose), polyList[poly_selected].crystalline[1]);
        }
        else if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) != par.free_poly_ends.end()){
            addreaction(par, false,allTables,poly_selected,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[1]);
        }
        //Second end
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, poly_selected,polyList[poly_selected].len_poly-2) == false){
            //cout << "In function EG_digest(): adding attachment reaction for length 3 polymer (old poly)" << endl;
            par.free_poly_ends.insert({cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2),make_tuple(poly_selected,polyList[poly_selected].len_poly-2)});
            if(true){
                addreaction(par, false,allTables,poly_selected,polyList[poly_selected].len_poly-2,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[polyList[poly_selected].len_poly-2]);
            }
            par.N_free_ends++;
        }
        else if(par.free_poly_ends.find(cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, poly_selected,polyList[poly_selected].len_poly-2) == true){
            addreaction(par, false,allTables,poly_selected,polyList[poly_selected].len_poly-2,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose), polyList[poly_selected].crystalline[polyList[poly_selected].len_poly-2]);
        }
        else if(par.free_poly_ends.find(cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2)) != par.free_poly_ends.end()){
            addreaction(par, false,allTables,poly_selected,polyList[poly_selected].len_poly-2,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[polyList[poly_selected].len_poly-2]);
        }

    }
    else if(polyList[poly_selected].len_poly == 3){
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,0)) != par.free_poly_ends.end()){
            par.free_poly_ends.erase(cantor_pair_two(poly_selected,0));
            deletereaction(allTables[poly_selected], allTables, 1, poly_selected, 0, 6);//Delete the reaction that will now be carried out
            par.N_free_ends--;
        }
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,2)) != par.free_poly_ends.end()){
            par.free_poly_ends.erase(cantor_pair_two(poly_selected,2));
            deletereaction(allTables[poly_selected], allTables, 1, poly_selected, 2, 6);//Delete the reaction that will now be carried out
            par.N_free_ends--;
        }
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, poly_selected,1) == false){
            //cout << "In function EG_digest(): adding attachment reaction for length 3 polymer (old poly)" << endl;
            par.free_poly_ends.insert({cantor_pair_two(poly_selected,1),make_tuple(poly_selected,1)});
            if(true){
                addreaction(par, false,allTables,poly_selected,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[1]);
            }
            par.N_free_ends++;
        }
        else if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, poly_selected,1) == true){
            addreaction(par, false,allTables,poly_selected,1,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose), polyList[poly_selected].crystalline[1]);
        }
        else if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) != par.free_poly_ends.end()){
            addreaction(par, false,allTables,poly_selected,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[1]);
        }
    }
    else if(polyList[poly_selected].len_poly == 2){
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,0)) == par.free_poly_ends.end() and CBH_enzyme_attached(par,poly_selected,0) == false){
            //cout << "In function EG_digest(): adding attachment reaction at position 0 for length 2 polymer (old poly)" << endl;
            par.free_poly_ends.insert({cantor_pair_two(poly_selected,0),make_tuple(poly_selected,0)});
            if(true){

                addreaction(par, false,allTables,poly_selected,0,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[0]);
            }
            par.N_free_ends++;
        }
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par,poly_selected,1) == false){
            //cout << "In function EG_digest(): adding attachment reaction at position 1 for length 2 polymer (old poly)" << endl;
            par.free_poly_ends.insert({cantor_pair_two(poly_selected,1),make_tuple(poly_selected,1)});
            if(true){
                addreaction(par, false,allTables,poly_selected,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[1]);
            }
            par.N_free_ends++;
        }
    }

    if(newPolyFlag == true){
        if(polyList[new_poly_index].len_poly > 3){
            //First end
            if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, new_poly_index,1) == false){
                //cout << "In function EG_digest(): adding attachment reaction for length 3 polymer (old poly)" << endl;
                par.free_poly_ends.insert({cantor_pair_two(new_poly_index,1),make_tuple(new_poly_index,1)});
                if(true){
                    addreaction(par, false,allTables,new_poly_index,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[1]);
                }
                par.N_free_ends++;
            }
            else if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, new_poly_index,1) == true){
                addreaction(par, false,allTables,new_poly_index,1,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose), polyList[new_poly_index].crystalline[1]);
            }
            else if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,1)) != par.free_poly_ends.end()){
                addreaction(par, false,allTables,new_poly_index,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[1]);
            }
            //Second end
            if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,polyList[new_poly_index].len_poly-2)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, new_poly_index,polyList[new_poly_index].len_poly-2) == false){
                //cout << "In function EG_digest(): adding attachment reaction for length 3 polymer (old poly)" << endl;
                par.free_poly_ends.insert({cantor_pair_two(new_poly_index,polyList[new_poly_index].len_poly-2),make_tuple(new_poly_index,polyList[new_poly_index].len_poly-2)});
                if(true){
                    addreaction(par, false,allTables,new_poly_index,polyList[new_poly_index].len_poly-2,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[polyList[new_poly_index].len_poly-2]);
                }
                par.N_free_ends++;
            }
            else if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,polyList[new_poly_index].len_poly-2)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, new_poly_index,polyList[new_poly_index].len_poly-2) == true){
                addreaction(par, false,allTables,new_poly_index,polyList[new_poly_index].len_poly-2,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose), polyList[new_poly_index].crystalline[polyList[new_poly_index].len_poly-2]);
            }
            else if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,polyList[new_poly_index].len_poly-2)) != par.free_poly_ends.end()){
                addreaction(par, false,allTables,new_poly_index,polyList[new_poly_index].len_poly-2,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[polyList[new_poly_index].len_poly-2]);
            }

        }
        else if(polyList[new_poly_index].len_poly == 3){
            if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,0)) != par.free_poly_ends.end()){
                par.free_poly_ends.erase(cantor_pair_two(new_poly_index,0));
                deletereaction(allTables[new_poly_index], allTables, 1, new_poly_index, 0, 6);//Delete the reaction that will now be carried out
                par.N_free_ends--;
            }
            if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,2)) != par.free_poly_ends.end()){
                par.free_poly_ends.erase(cantor_pair_two(new_poly_index,2));
                deletereaction(allTables[new_poly_index], allTables, 1, new_poly_index, 2, 6);//Delete the reaction that will now be carried out
                par.N_free_ends--;
            }
            if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, new_poly_index,1) == false){
                //cout << "In function EG_digest(): adding attachment reaction for length 3 polymer (new poly)" << endl;
                par.free_poly_ends.insert({cantor_pair_two(new_poly_index,1),make_tuple(new_poly_index,1)});
                if(true){
                    addreaction(par, false,allTables,new_poly_index,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[1]);
                }
                par.N_free_ends++;
            }
            else if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, new_poly_index,1) == true){
                addreaction(par, false,allTables,new_poly_index,1,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose), polyList[new_poly_index].crystalline[1]);
            }
            else if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,1)) != par.free_poly_ends.end()){
                addreaction(par, false,allTables,new_poly_index,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[1]);
            }
        }
        else if(polyList[new_poly_index].len_poly == 2){
            if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,0)) == par.free_poly_ends.end() and CBH_enzyme_attached(par,new_poly_index,0) == false){
                //cout << "In function EG_digest(): adding attachment reaction at position 0 for length 2 polymer (new poly)" << endl;
                par.free_poly_ends.insert({cantor_pair_two(new_poly_index,0),make_tuple(new_poly_index,0)});
                if(true){
                    addreaction(par, false,allTables,new_poly_index,0,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[0]);
                }
                par.N_free_ends++;
            }
            if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par,new_poly_index,1) == false){
                //cout << "In function EG_digest(): adding attachment reaction at position 1 for length 2 polymer (new poly)" << endl;
                par.free_poly_ends.insert({cantor_pair_two(new_poly_index,1),make_tuple(new_poly_index,1)});
                if(true){
                    addreaction(par, false,allTables,new_poly_index,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[1]);
                }
                par.N_free_ends++;
            }
        }
    }
    if(par.verbose == true){
        cout << "Done with EG reaction" << endl;
    }
//End of function
}

void CBH_digest(params& par, vector<TList>& allTables, vector<bList>& polyList, int& nbr_poly, int& nbr_xyl_pdt, int& nbr_Glc_pdt,const int bond_selected,const int poly_selected, int len_polyLoopStart, int& nbr_cellobiose,vector<double>& chem_entities, const int substrate, unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu, unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_hemi,const bool verbose, bool& error){


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
    std::tuple<double,double,int> neighbor_key = std::make_tuple(x_pos, y_pos, z_pos);

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
    else if(bond_selected==polyList[poly_selected].len_poly-2)//If the cut takes place at the tip of the polymer
        indic_cut=1;
    else{
        cout << "CBH reaction which is neither at position 1 nor len_poly-2. It is at position " << bond_selected << ", and the length of the poly is " << polyList[poly_selected].len_poly << endl;

        cout << "Stopping" << endl;
        exit(1);
    }

    deleteAllreaction(allTables[poly_selected],allTables,substrate,poly_selected,bond_selected);

    //Remove the bond and split the polymer into two new polymers
    if((bond_selected !=0) and (bond_selected !=polyList[poly_selected].len_poly-1)){//This is the case for all polys with len_poly > 2
        polyList.push_back(bList());//Make the new polymer of cellobiose
        new_poly_index = polyList.size()-1;
        newPolyFlag = true;
        polyList[new_poly_index].index= new_poly_index;//Specify the new index of the new polymer
        polyList[new_poly_index].x= polyList[poly_selected].x;
        polyList[new_poly_index].y= polyList[poly_selected].y;
        polyList[new_poly_index].set_z(bond_selected+1);
        allTables.push_back(TList());//Make the new reaction list associated to the new polymer;
        initTList(allTables[new_poly_index], new_poly_index);

        //Check if we need to update a free end
        bool free_end_bool = false;
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2)) != par.free_poly_ends.end()){
            par.free_poly_ends.erase(cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2));
            par.N_free_ends--;
            free_end_bool = true;
        }
        if(indic_cut == 0){//Cut was at beginning
            addbond(polyList[new_poly_index],0,1,1,false);
            polyList[new_poly_index].z[0]=polyList[poly_selected].z[0];
            polyList[new_poly_index].status[0]=polyList[poly_selected].status[0];
            polyList[new_poly_index].bond_type[0]=polyList[poly_selected].bond_type[0];
            polyList[new_poly_index].N_blocked_positions[0] = polyList[poly_selected].N_blocked_positions[0];
            polyList[new_poly_index].crystalline[0] = polyList[poly_selected].crystalline[0];

            //add one free bond at the beginning and use taylorNewPoly on old poly
            taylorNewPoly(polyList[poly_selected],bond_selected);//Erase the beginning of the polymer, until digested bond

            //WE DO NOT HAVE TO UPDATE THE OTHER END, BECAUSE CBH IS STILL ATTACHED THERE (PROCESSIVITY)
        }
        else if(indic_cut == 1){//Cut was at the end

            addbond(polyList[new_poly_index],0,1,1,false);
            polyList[new_poly_index].z[0]=polyList[poly_selected].z[polyList[poly_selected].len_poly-1];
            polyList[new_poly_index].status[0]=polyList[poly_selected].status[polyList[poly_selected].len_poly-1];
            polyList[new_poly_index].bond_type[0]=polyList[poly_selected].bond_type[polyList[poly_selected].len_poly-1];
            polyList[new_poly_index].N_blocked_positions[0] = polyList[poly_selected].N_blocked_positions[polyList[poly_selected].len_poly-1];
            polyList[new_poly_index].crystalline[0] = polyList[poly_selected].crystalline[polyList[poly_selected].len_poly-1];
            taylorOldPoly(polyList[poly_selected],bond_selected);//Erase the end of the polymer, until digested bond
        }


        //Clear old table
        clear_table(par,allTables,poly_selected);
        //Fill old table
        fill_table(par,allTables,polyList,poly_selected,substrate,error,chem_entities,nbr_Glc_pdt,nbr_cellobiose);


        //Fill new table
        fill_table(par,allTables, polyList, new_poly_index, substrate, error,chem_entities,nbr_Glc_pdt,nbr_cellobiose);

        if(free_end_bool == true){
            if(par.free_poly_ends.find(cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2)) == par.free_poly_ends.end() and CBH_enzyme_attached(par,poly_selected,polyList[poly_selected].len_poly-2) == false){
                if(polyList[poly_selected].len_poly-2 >= 0){
                    par.free_poly_ends.insert({cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2),std::make_tuple(poly_selected,polyList[poly_selected].len_poly-2)});
                    if(true){
                        addreaction(par, false,allTables,poly_selected,polyList[poly_selected].len_poly-2,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[polyList[poly_selected].len_poly-2]);
                    }
                    par.N_free_ends++;
                }
            }
        }
        //UPDATE CBH ENZYMES ATTACHED
        int CBH_count = 0;//max two CBH enzymes can be attached to one polymer --> this is a cut-off condition for the following for loop
        for(int i=0; i<par.CBH_enzymes.size();i++){
            if(par.CBH_enzymes[i].attached == true and par.CBH_enzymes[i].poly_attached == poly_selected and par.CBH_enzymes[i].bond_attached == 1){
                CBH_count++;
                if(polyList[poly_selected].len_poly > 1){
//                    (BOND_ATTACHED AND POLY_ATTACHED DO NOT NEED TO BE CHANGED)
                    if(indic_cut == 0){//Cut was at beginning, meaning that the enzyme needs to move up
                        addreaction(par, false, allTables,poly_selected,1,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[1]);
                    }
                    else if(indic_cut == 1){//Cut was at end. The CBH reaction needs to be re-added, because we cleared the table
                        addreaction(par, false, allTables,poly_selected,1,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[1]);
                    }
                }
                else{
                    par.CBH_enzymes[i].detach();
                    par.N_free_CBH++;
                }

            }
            else if(par.CBH_enzymes[i].attached == true and par.CBH_enzymes[i].poly_attached == poly_selected and (par.CBH_enzymes[i].bond_attached != 1)){//== len_poly because the enzyme thinks that the polymer is 2 bonds longer than it is
                CBH_count++;
                if(polyList[poly_selected].len_poly != 1){
                    if(indic_cut == 0 and polyList[poly_selected].len_poly != 3){//Cut was at beginning
                        par.CBH_enzymes[i].bond_attached = polyList[poly_selected].len_poly-2;
                        if(par.CBH_enzymes[i].bond_attached == 0 and polyList[par.CBH_enzymes[i].poly_attached].len_poly != 2){
                            cout << "In CBH_digest: CBH reaction for a bond 0 is about to be added to a polymer that is not of length 2. This should not happen. Stopping" << endl;
                            exit(1);
                        }
                        else{
                            addreaction(par, false, allTables,poly_selected,par.CBH_enzymes[i].bond_attached,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[par.CBH_enzymes[i].bond_attached]);
                        }
                    }
                    else if(indic_cut == 1 and polyList[poly_selected].len_poly != 3){//Cut was at end
                        par.CBH_enzymes[i].bond_attached = polyList[poly_selected].len_poly-2;
                        if(par.CBH_enzymes[i].bond_attached == 0 and polyList[par.CBH_enzymes[i].poly_attached].len_poly != 2){
                            cout << "In CBH_digest: CBH reaction for a bond 0 is about to be added to a polymer that is not of length 2. This should not happen. Stopping" << endl;
                            exit(1);
                        }
                        else{
                            addreaction(par, false, allTables,poly_selected,par.CBH_enzymes[i].bond_attached,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[par.CBH_enzymes[i].bond_attached]);
                        }
                    }
                    else{
                        par.CBH_enzymes[i].detach();
                        par.N_free_CBH++;
                    }
                }
                else{
                    par.CBH_enzymes[i].detach();
                    par.N_free_CBH++;
                }
            }
/*            if(CBH_count > 2){
                cout << "More than two CBH attached to one poly! Stopping..." << endl;
                exit(1);
            }*/
            if(CBH_count == 2){
                break;
            }
        }


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
            clear_table(par,allTables,poly_selected);

            if(prop(par,3,chem_entities,nbr_Glc_pdt,nbr_cellobiose) > 0 and polyList[poly_selected].len_poly == 1){
                polyList[poly_selected].status[0] = 1;
                reaction_type = 3;
                propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                crystal_propensite = par.crystal_modifier_cellu*propensite;
                if(polyList[poly_selected].crystalline[0] == false)
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,propensite,0);
                else
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,crystal_propensite,1);
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
            clear_table(par,allTables,poly_selected);

            if(prop(par,3,chem_entities,nbr_Glc_pdt,nbr_cellobiose) > 0 and polyList[poly_selected].len_poly == 1){
                polyList[poly_selected].status[0] = 1;
                reaction_type = 3;
                propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                crystal_propensite = par.crystal_modifier_cellu*propensite;
                if(polyList[poly_selected].crystalline[0] == false)
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,propensite,0);
                else
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,crystal_propensite,1);
            }

            nbr_Glc_pdt++;//We release a Glc

        }
    }



    //REACTION TABLE MODIFICATION IS CARRIED OUT IN SEPARATE FUNCTION
    if(newPolyFlag == true){
        if(polyList[new_poly_index].len_poly > 1){
            cout << "In function CBH_digest(): a newly generated polymer is longer than cellobiose. This should not happen. Stopping." << endl;
            exit(1);
        }
        else if(polyList[new_poly_index].len_poly <= 0){
             cout << "In function CBH_digest(): a newly generated polymer has zero length. This should not happen. Stopping." << endl;
             exit(1);
        }
    }

//    cout << "CBH_digest: len_poly afterwards: " << polyList[poly_selected].len_poly << endl;
    //Check if either polymer has length 3 or 2
    if(polyList[poly_selected].len_poly > 3){
        //First end
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, poly_selected,1) == false){
            //cout << "In function EG_digest(): adding attachment reaction for length 3 polymer (old poly)" << endl;
            par.free_poly_ends.insert({cantor_pair_two(poly_selected,1),make_tuple(poly_selected,1)});
            if(true){
                addreaction(par, false,allTables,poly_selected,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[1]);
            }
            par.N_free_ends++;
        }
        else if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, poly_selected,1) == true){
            addreaction(par, false,allTables,poly_selected,1,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose), polyList[poly_selected].crystalline[1]);
        }
        else if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) != par.free_poly_ends.end()){
            addreaction(par, false,allTables,poly_selected,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[1]);
        }
        //Second end
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, poly_selected,polyList[poly_selected].len_poly-2) == false){
            //cout << "In function EG_digest(): adding attachment reaction for length 3 polymer (old poly)" << endl;
            par.free_poly_ends.insert({cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2),make_tuple(poly_selected,polyList[poly_selected].len_poly-2)});
            if(true){
                addreaction(par, false,allTables,poly_selected,polyList[poly_selected].len_poly-2,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[polyList[poly_selected].len_poly-2]);
            }
            par.N_free_ends++;
        }
        else if(par.free_poly_ends.find(cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, poly_selected,polyList[poly_selected].len_poly-2) == true){
            addreaction(par, false,allTables,poly_selected,polyList[poly_selected].len_poly-2,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose), polyList[poly_selected].crystalline[polyList[poly_selected].len_poly-2]);
        }
        else if(par.free_poly_ends.find(cantor_pair_two(poly_selected,polyList[poly_selected].len_poly-2)) != par.free_poly_ends.end()){
            addreaction(par, false,allTables,poly_selected,polyList[poly_selected].len_poly-2,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[polyList[poly_selected].len_poly-2]);
        }

    }
    else if(polyList[poly_selected].len_poly == 3){
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,0)) != par.free_poly_ends.end()){
            par.free_poly_ends.erase(cantor_pair_two(poly_selected,0));
            deletereaction(allTables[poly_selected], allTables, 1, poly_selected, 0, 6);//Delete the reaction
            par.N_free_ends--;
        }
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,2)) != par.free_poly_ends.end()){
            par.free_poly_ends.erase(cantor_pair_two(poly_selected,2));
            deletereaction(allTables[poly_selected], allTables, 1, poly_selected, 2, 6);//Delete the reaction
            par.N_free_ends--;
        }
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, poly_selected,1) == false){
            par.free_poly_ends.insert({cantor_pair_two(poly_selected,1),make_tuple(poly_selected,1)});
            addreaction(par, false,allTables,poly_selected,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose), polyList[poly_selected].crystalline[1]);
            par.N_free_ends++;
        }
        else if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, poly_selected,1) == true){
            addreaction(par, false,allTables,poly_selected,1,1,2,prop(par,2,chem_entities,nbr_Glc_pdt,nbr_cellobiose), polyList[poly_selected].crystalline[1]);
        }
        else if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) != par.free_poly_ends.end()){
            addreaction(par, false,allTables,poly_selected,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[1]);
        }
    }
    else if(polyList[poly_selected].len_poly == 2){
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,0)) == par.free_poly_ends.end() and CBH_enzyme_attached(par,poly_selected,0) == false){
//            cout << "In function CBH_digest(): adding attachment reaction at position 0 for length 2 polymer" << endl;
            par.free_poly_ends.insert({cantor_pair_two(poly_selected,0),make_tuple(poly_selected,0)});
            if(true){
                addreaction(par, false,allTables,poly_selected,0,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[0]);
            }
            par.N_free_ends++;
        }
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par,poly_selected,1) == false){
//            cout << "In function CBH_digest(): adding attachment reaction at position 1 for length 2 polymer" << endl;
            par.free_poly_ends.insert({cantor_pair_two(poly_selected,1),make_tuple(poly_selected,1)});
            if(true){
                addreaction(par, false,allTables,poly_selected,1,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[1]);
            }
            par.N_free_ends++;
        }
    }
    else if(polyList[poly_selected].len_poly == 1){
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,0)) != par.free_poly_ends.end()){
            par.free_poly_ends.erase(cantor_pair_two(poly_selected,0));
            deletereaction(allTables[poly_selected], allTables, 1, poly_selected, 0, 6);//Delete the reaction
            par.N_free_ends--;
        }
        if(par.free_poly_ends.find(cantor_pair_two(poly_selected,1)) != par.free_poly_ends.end()){
            par.free_poly_ends.erase(cantor_pair_two(poly_selected,1));
            deletereaction(allTables[poly_selected], allTables, 1, poly_selected, 1, 6);//Delete the reaction
            par.N_free_ends--;
        }
        addreaction(par, false, allTables,poly_selected,0,1,3,prop(par,3,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[poly_selected].crystalline[0]);
        nbr_cellobiose++;
        int count_CBHs = 0;
        for(int i=0;i<par.CBH_enzymes.size();i++){
            if(par.CBH_enzymes[i].poly_attached == poly_selected){
                count_CBHs++;
                par.CBH_enzymes[i].detach();
                par.N_free_CBH++;
            }
            if(count_CBHs == 2){
                break;
            }
        }
    }
    if(newPolyFlag == true){
        if(polyList[new_poly_index].len_poly != 1){
            cout << "In function CBH_digest(): len_poly of new_poly_index != 1. This should not currently be possible. Stopping" << endl;
            exit(1);
        }
        if(polyList[new_poly_index].len_poly == 1){
            if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,0)) != par.free_poly_ends.end()){
                par.free_poly_ends.erase(cantor_pair_two(new_poly_index,0));
                deletereaction(allTables[new_poly_index], allTables, 1, new_poly_index, 0, 6);//Delete the reaction
                par.N_free_ends--;
            }
            if(par.free_poly_ends.find(cantor_pair_two(new_poly_index,1)) != par.free_poly_ends.end()){
                par.free_poly_ends.erase(cantor_pair_two(new_poly_index,1));
                deletereaction(allTables[new_poly_index], allTables, 1, new_poly_index, 1, 6);//Delete the reaction
                par.N_free_ends--;
            }
            addreaction(par, false, allTables,new_poly_index,0,1,3,prop(par,3,chem_entities,nbr_Glc_pdt,nbr_cellobiose),polyList[new_poly_index].crystalline[0]);
            nbr_cellobiose++;
            int count_CBHs = 0;
            for(int i=0;i<par.CBH_enzymes.size();i++){
                if(par.CBH_enzymes[i].poly_attached == new_poly_index){
                    count_CBHs++;
                    par.CBH_enzymes[i].detach();
                    par.N_free_CBH++;
                }
                if(count_CBHs == 2){
                    break;
                }
            }
        }
    }

    if(par.verbose == true){
        cout << "Done with CBH reaction" << endl;
    }
//End of function
}



void BGL_digest(params& par, vector<TList>& allTables, vector<bList>& polyList, int& nbr_poly, int& nbr_xyl_pdt, int& nbr_Glc_pdt,const int bond_selected,const int poly_selected, int len_polyLoopStart, int& nbr_cellobiose,vector<double>& chem_entities, const int substrate,unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu, unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_hemi, const bool verbose, bool& error){

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

    std::tuple<double,double,int> neighbor_key = std::make_tuple(x_pos, y_pos, z_pos);



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

    //REACTION TABLE MODIFICATION IS CARRIED OUT IN SEPARATE FUNCTION

    if(par.verbose == true){
        cout << "Done with BGL reaction" << endl;
    }
//End of function
}

//Digestion by Xylanase
void XYL_digest(params& par, vector<TList>& allTables, vector<TList>& allCelluTables, vector<bList>& polyList, vector<bList>& celluList, int& nbr_poly, int& nbr_xyl_pdt, int& nbr_Glc_pdt,const int bond_selected,const int poly_selected, int len_polyLoopStart, int& nbr_cellobiose,vector<double>& chem_entities, const int substrate,unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu, unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_hemi, const bool verbose, bool& error){
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


    std::tuple<double,double,int> neighbor_key = std::make_tuple(x_pos, y_pos, z_pos);

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

        if (par.xyl_or_mlg == true){
            nbr_xyl_pdt += +2;  // partho count xylose release since Xylans in Hemi
        }
        else if (par.xyl_or_mlg == false){
            nbr_Glc_pdt += +2;  // partho counts glucose release since MLG present   
        }
        
        allTables[poly_selected].prop_sum = 0;
    }
    else if(polyList[poly_selected].len_poly==2)
    {

        if(bond_selected==0)// If the cut was at the beginning of the polymer
        {
            taylorNewPoly(polyList[poly_selected],bond_selected);//Shift the indices by one to the left
            clear_table(par,allTables,poly_selected);

            reaction_type = 4;
            propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
            crystal_propensite = par.crystal_modifier_hemi*propensite;
            
            if (par.xyl_or_mlg == true){    //partho Xylose in hemi... so bondtype 4
                if(propensite > 0 and polyList[poly_selected].status[0] == 1 and polyList[poly_selected].len_poly == 1 and polyList[poly_selected].bond_type[0] == 4){  //bond type = 1 for mlg, 4 for xyl      partho
                if(polyList[poly_selected].crystalline[0] == false)
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,propensite,0);
                else
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,crystal_propensite,1);
                }
                nbr_xyl_pdt++;  //  partho release a Xylose since xylan present in hemi 
            }
            else if (par.xyl_or_mlg == false){  //  partho MLG in hemi... so bondtype 1
                if(propensite > 0 and polyList[poly_selected].status[0] == 1 and polyList[poly_selected].len_poly == 1 and polyList[poly_selected].bond_type[0] == 1){  //bond type = 1 for mlg, 4 for xyl      partho
                if(polyList[poly_selected].crystalline[0] == false)
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,propensite,0);
                else
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,crystal_propensite,1);
                }
                nbr_Glc_pdt++;  //  partho release a Glucose since MLG present in hemi
            }
        }
        else if(bond_selected == 1){
            taylorOldPoly(polyList[poly_selected],bond_selected);//Erase the end of the polymer, until digested bond
            clear_table(par,allTables,poly_selected);

            reaction_type = 4;
            propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
            crystal_propensite = par.crystal_modifier_hemi*propensite;

            if (par.xyl_or_mlg == true){    //  partho Xylose in hemi... so bondtype 4
                if(propensite > 0 and polyList[poly_selected].status[0] == 1 and polyList[poly_selected].len_poly == 1 and polyList[poly_selected].bond_type[0] == 4){  //bond type = 1 for mlg, 4 for xyl      partho
                if(polyList[poly_selected].crystalline[0] == false)
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,propensite,0);
                else
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,crystal_propensite,1);
                }
                nbr_xyl_pdt++;  // partho release a Xylose
            }
            else if (par.xyl_or_mlg == false){  //partho MLG in hemi... so bondtype 1
                if(propensite > 0 and polyList[poly_selected].status[0] == 1 and polyList[poly_selected].len_poly == 1 and polyList[poly_selected].bond_type[0] == 1){  //bond type = 1 for mlg, 4 for xyl      partho
                if(polyList[poly_selected].crystalline[0] == false)
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,propensite,0);
                else
                    addreaction(par, false, allTables,poly_selected,0,substrate,reaction_type,crystal_propensite,1);
                }
                nbr_Glc_pdt++;  //  partho release a Glucose
            }
		    
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
            clear_table(par,allTables,poly_selected);
            fill_table(par,allTables,polyList,poly_selected,substrate,error,chem_entities,nbr_Glc_pdt,nbr_cellobiose);

            //Fill new table
            fill_table(par,allTables, polyList, new_poly_index, substrate, error,chem_entities,nbr_Glc_pdt,nbr_cellobiose);


            nbr_poly++;
        }
        else if(bond_selected == 0){
            taylorNewPoly(polyList[poly_selected],bond_selected);//Erase the begining of the new polymer, from digested bond
            //Adjust old table
            clear_table(par,allTables,poly_selected);
            fill_table(par,allTables,polyList,poly_selected,substrate,error,chem_entities,nbr_Glc_pdt,nbr_cellobiose);

            if (par.xyl_or_mlg == true){
                nbr_xyl_pdt++;  //partho xylose released since only xylan in hemi
            }
            if (par.xyl_or_mlg == false){
                nbr_Glc_pdt++;  //partho glucose released since only MLG in hemi 
            }
        }
        else if(bond_selected == polyList[poly_selected].len_poly-1){
            taylorOldPoly(polyList[poly_selected],bond_selected);//Erase the end of the polymer, until digested bond
            //Adjust old table
            clear_table(par,allTables,poly_selected);
            fill_table(par,allTables,polyList,poly_selected,substrate,error,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
            if (par.xyl_or_mlg == true){
                nbr_xyl_pdt++;  //partho xylose released since only xylan in hemi
            }
            if (par.xyl_or_mlg == false){
                nbr_Glc_pdt++;  //partho glucose released since only MLG in hemi 
            }
        }
    }


    int cellu_substrate = 1;

    //REACTION TABLE MODIFICATION IS CARRIED OUT IN SEPARATE FUNCTION







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


    if(par.verbose == true){
        cout << "Done with XYL reaction" << endl;
    }
//End of function
}

void CBH_attachment(params& par, std::vector<TList>& allTables,
                    std::vector<bList>& celluList,
                    std::vector<bList>& hemiList,
                    std::vector<bList>& lignList,
                    int& nbr_poly,
                    const int bond_selected,
                    const int poly_selected,
                    std::vector<double>& chem_entities,
                    const int substrate,
                    std::unordered_map<std::tuple<double,double,int>,neighborList>& bond_neighbors_cellu,
                    std::unordered_map<std::tuple<double,double,int>,neighborList>& bond_neighbors_hemi,
                    int& nbr_Glc_pdt, int& nbr_cellobiose
                    )
{
    if(par.k6 <= 0.0){
        cout << "CBH attachment reaction was called, even though k6 = 0. Stopping" << endl;
        exit(1);
    }
    double dx = 0.;
    double dy = 0.;
    double dz = 0.;
    int poly_attached = 0;
    int bond_attached = 0;


/*
    for(int k=0;k<par.CBH_enzymes.size();k++){
        if(par.CBH_enzymes[k].attached == true){
            poly_attached = par.CBH_enzymes[k].poly_attached;
            bond_attached = par.CBH_enzymes[k].bond_attached;
            dx = celluList[poly_selected].x - celluList[poly_attached].x;
            dy = celluList[poly_selected].y - celluList[poly_attached].y;
            dz = celluList[poly_selected].z[bond_selected] - celluList[poly_attached].z[bond_attached];
            if(dist(dx,dy,dz) < 2*par.enzyme_radius){
                cout << "In function CBH_attachment(): bond_selected should be covered by another CBH enzyme, i.e. the attachment reaction should not be pickable for this bond. Stopping" << endl;
                exit(1);
            }
        }
    }*/


    deletereaction(allTables[poly_selected],allTables,1,poly_selected,bond_selected,6);//Delete the reaction that is now taking place
    if(par.N_free_CBH > 0){
        for(int i = 0; i<par.CBH_enzymes.size();i++){
            if(par.CBH_enzymes[i].attached == false){
                double attachment_duration = box_muller(par.average_CBH_association_time*(par.std_dev_CBH_speed/par.average_CBH_speed),par.average_CBH_association_time);
                    while(par.real_time + attachment_duration == par.real_time){
                    attachment_duration = box_muller(par.average_CBH_association_time*(par.std_dev_CBH_speed/par.average_CBH_speed),par.average_CBH_association_time);
                }
                if(CBH_enzyme_attached(par, poly_selected, bond_selected) == false){
                    par.CBH_enzymes[i].attach(par.real_time, attachment_duration, poly_selected,bond_selected,celluList);
                    par.CBH_enzymes[i].init_neighbors(celluList[poly_selected].x,celluList[poly_selected].y,celluList[poly_selected].z[bond_selected],par.enzyme_radius,par.mid_x,par.mid_y,par.r_monomer,celluList,hemiList,lignList);
                    par.N_free_CBH--;
                    par.free_poly_ends.erase(cantor_pair_two(poly_selected,bond_selected));
                    par.N_free_ends--;
                    int reaction_type = 2;
                    double propensite = prop(par, reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                    //Add the new digestion reaction for the attached CBH
                    addreaction(par, false, allTables, poly_selected, bond_selected, 1, reaction_type, propensite,celluList[poly_selected].crystalline[bond_selected]);
                    //Propensity update of attachment reaction is done in the function update_CBH_attachment_reactions()
                    return;
                }
                else{
                    cout << "In function CBH_attachment(): There is already a CBH attached at bond_selected in poly_selected. It should not be possible to call this function in this case. Stopping" << endl;
//                    return;
                    exit(1);
                }
            }
        }

    }
    else{
        cout << "CBH_attachment reaction was called, even though there are no free CBHs. This should not happen. Stopping" << endl;
        exit(1);
    }

//End of function
}


//Updates reaction tables after a digestion reaction. Called at each step that is not a lignin adhesion step
void update_reactiontables(std::vector<bList>& cellu, std::vector<bList>& hemi, std::vector<bList>& lign, std::vector<TList>& Table_cellu, std::vector<TList>& Table_hemi, params& par, std::vector<double>& chem_entities, std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu, std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_hemi, double x, double y, int z, int substrate, int poly_selected, int bond_selected, int& nbr_Glc_pdt, int& nbr_cellobiose){
    std::tuple<double,double,int> bond_key = std::make_tuple(x,y,z);//Used to find the neighborList object associated to the previously digested bond at coordinates x,y,z
    std::tuple<double,double,int> neighbor_key = std::make_tuple(0.,0.,0);
    int current_poly = 0;
    int current_bond = 0;
    int current_substrate = 0;
    int reaction_type = 0;
    int enzyme_index = -1;
    double propensite = 0;
    double crystal_propensite = 0;
    //=============================================================================================================================================================================================
    //======================================== IF THE PREVIOUSLY DIGESTED BOND WAS IN A CELLULOSE POLYMER =========================================================================================
    //=============================================================================================================================================================================================
    if(substrate == 1){
        if(bond_neighbors_cellu.count(bond_key) != 0){
            for(int i=0;i<bond_neighbors_cellu[bond_key].N_neighbors; i++){
//                if(par.verbose == true)
 //                   cout << "neighbor coordinates: " << bond_neighbors_cellu[bond_key].x_neighbors[i] << ", " << bond_neighbors_cellu[bond_key].y_neighbors[i] << ", " << bond_neighbors_cellu[bond_key].z_neighbors[i] << endl;
                neighbor_key = std::make_tuple(bond_neighbors_cellu[bond_key].x_neighbors[i],bond_neighbors_cellu[bond_key].y_neighbors[i], bond_neighbors_cellu[bond_key].z_neighbors[i]);
                //=============================================================================================================================================================================================
                //======================================== First we check the cellu neighbors =================================================================================================================
                //=============================================================================================================================================================================================
                if(bond_neighbors_cellu.count(neighbor_key) != 0){

                    bond_neighbors_cellu[neighbor_key].remove_neighbor(x,y,z);//Delete the digested bond from the list of neighbors of the neighbor
                    bond_neighbors_cellu[neighbor_key].outer_bond = bond_neighbors_cellu[neighbor_key].is_outer_bond();
                    if(bond_neighbors_cellu[neighbor_key].outer_bond == true){ //or par.mode_enzyme_size == -1){//See if the neighbor is now digestible
                        //The function find_specific_bond() takes current_poly, current_bond and current_substrate as references and returns the required indices within the cellu vector, if the bond exists
                        if(find_specific_bond(cellu, hemi, lign, bond_neighbors_cellu[neighbor_key].x, bond_neighbors_cellu[neighbor_key].y, bond_neighbors_cellu[neighbor_key].z, current_poly, current_bond, current_substrate) == true){
                            if(current_substrate == 1){
                                if(cellu[current_poly].status[current_bond] == -1){
                                    cellu[current_poly].status[current_bond] = 1;
                                    if(current_bond == 1){//CHECK IF WE HAVE TO ADD FREE ENDS
                                        if(cellu[current_poly].len_poly > 1){
                                            if(CBH_enzyme_attached(par,current_poly,1) == false){
                                                if(par.free_poly_ends.find(cantor_pair_two(current_poly,1)) == par.free_poly_ends.end() and CBH_enzyme_attached(par, current_poly, 1) == false){
                                                    par.free_poly_ends.insert({cantor_pair_two(current_poly,1),std::make_tuple(current_poly,1)});
                                                    par.N_free_ends++;
                                                    if(true){
                                                        addreaction(par, false,Table_cellu,current_poly,current_bond,current_substrate,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),cellu[current_poly].crystalline[current_bond]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else if(current_bond == cellu[current_poly].len_poly-2){
                                        if(cellu[current_poly].len_poly > 1){
                                            if(CBH_enzyme_attached(par,current_poly,cellu[current_poly].len_poly-2) == false){
                                                if(par.free_poly_ends.find(cantor_pair_two(current_poly,cellu[current_poly].len_poly-2)) == par.free_poly_ends.end() and CBH_enzyme_attached(par,current_poly,cellu[current_poly].len_poly-2) == false){
                                                    par.free_poly_ends.insert({cantor_pair_two(current_poly,cellu[current_poly].len_poly-2),std::make_tuple(current_poly,cellu[current_poly].len_poly-2)});
                                                    if(true){
                                                        addreaction(par, false,Table_cellu,current_poly,current_bond,current_substrate,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),cellu[current_poly].crystalline[current_bond]);
                                                    }
                                                    par.N_free_ends++;
                                                }
                                            }
                                        }
                                    }
                                    if(cellu[current_poly].len_poly == 1){
                                        reaction_type = 3;
                                        propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;
                                        if(propensite > 0){
                                            if(cellu[current_poly].crystalline[current_bond] == false){
                                                addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,propensite,0);
                                            }
                                            else{
                                                addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,crystal_propensite,1);
                                            }
                                        }
                                    }
                                    //CBH REACTION ADDING IS DONE IN FUNCTION CBH_ATTACHMENT()
/*                                    else if(cellu[current_poly].len_poly == 2){
                                        reaction_type = 2;
                                        propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;

                                        propensite = 0;
                                        if(propensite > 0){
                                           if(cellu[current_poly].crystalline[current_bond] == false){
                                                addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,propensite,enzyme_index);

                                            }
                                            else{
                                                addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,crystal_propensite,enzyme_index);

                                            }
                                        }
                                    }*/
                                    else if(cellu[current_poly].len_poly > 2){

                                        reaction_type = 1;
                                        propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;

                                        if(current_bond > 1 and current_bond < cellu[current_poly].len_poly-2){
                                            if(propensite > 0){
                                               if(cellu[current_poly].crystalline[current_bond] == false){
                                                    addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,propensite,0);
                                                }
                                                else{
                                                    addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,crystal_propensite,1);
                                                }
                                            }

                                        }
                                        else if(cellu[current_poly].len_poly == 4){
                                        /*
                                            if(current_bond == 1 or current_bond == 2){
                                                if(propensite > 0){
                                                   if(cellu[current_poly].crystalline[current_bond] == false){
                                                        addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,propensite,0);
                                                    }
                                                    else{
                                                        addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,crystal_propensite,1);
                                                    }
                                                }
                                            }
                                        */
                                        }
                                    //CBH REACTION ADDING IS DONE IN FUNCTION CBH_ATTACHMENT()
                                        /*if(current_bond == 1 or current_bond == cellu[current_poly].len_poly-2){
                                            reaction_type = 2;
                                            propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                                            propensite = 0;
                                            if(propensite > 0){
                                               if(cellu[current_poly].crystalline[current_bond] == false){
                                                    addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,propensite,enzyme_index);

                                                }
                                                else{
                                                    addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,crystal_propensite,enzyme_index);

                                                }
                                            }
                                        }*/

                                    }
                                }
                            }
                            else if(current_substrate == 2){
                                cout << "In function update_reaction_tables(): current_substrate = 2 for a neighbor that should be contained in cellu. Stopping..." << endl;
                                exit(1);
                            }
                        }
                    }
                }
                //=============================================================================================================================================================================================
                //======================================== Next we check the hemi neighbors ===================================================================================================================
                //=============================================================================================================================================================================================
                else if(bond_neighbors_hemi.count(neighbor_key) != 0){
                    bond_neighbors_hemi[neighbor_key].remove_neighbor(x,y,z);
                    bond_neighbors_hemi[neighbor_key].outer_bond = bond_neighbors_hemi[neighbor_key].is_outer_bond();
                    if(bond_neighbors_hemi[neighbor_key].outer_bond == true ){//or par.mode_enzyme_size == -1){
                        if(find_specific_bond(cellu, hemi, lign, bond_neighbors_hemi[neighbor_key].x, bond_neighbors_hemi[neighbor_key].y, bond_neighbors_hemi[neighbor_key].z, current_poly, current_bond, current_substrate) == true){
         //                   cout << "found cellu" << endl;
                            if(current_substrate == 1){/**/
                                cout << "In function update_reaction_tables(): current_substrate = 1 for a neighbor that should be contained in hemi. Stopping..." << endl;
                                exit(1);
                            }
                            else if(current_substrate == 2){
                                if(hemi[current_poly].status[current_bond] == -1){
                                    hemi[current_poly].status[current_bond] = 1;
                                    hemi[current_poly].status[current_bond] = 1;
                                    reaction_type = 4;
                                    propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                                    crystal_propensite = par.crystal_modifier_hemi * propensite;

                                    if(propensite > 0 and current_bond < hemi[current_poly].len_poly){
                                       if(hemi[current_poly].crystalline[current_bond] == false){
                                            addreaction(par, false, Table_hemi,current_poly, current_bond,2, reaction_type,propensite,0);
                                        }
                                        else{
                                            addreaction(par, false, Table_hemi,current_poly, current_bond,2, reaction_type,crystal_propensite,1);
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
    //=============================================================================================================================================================================================
    //======================================== IF THE PREVIOUSLY DIGESTED BOND WAS IN A HEMICELLULOSE POLYMER =========================================================================================
    //=============================================================================================================================================================================================
    else if(substrate == 2){
        if(bond_neighbors_hemi.count(bond_key) != 0){
            for(int i=0;i<bond_neighbors_hemi[bond_key].N_neighbors; i++){
                if(par.verbose == true)
                    cout << "neighbor coordinates: " << bond_neighbors_hemi[bond_key].x_neighbors[i] << ", " << bond_neighbors_hemi[bond_key].y_neighbors[i] << ", " << bond_neighbors_hemi[bond_key].z_neighbors[i] << endl;
                neighbor_key = std::make_tuple(bond_neighbors_hemi[bond_key].x_neighbors[i],bond_neighbors_hemi[bond_key].y_neighbors[i], bond_neighbors_hemi[bond_key].z_neighbors[i]);

                //=============================================================================================================================================================================================
                //======================================== First we check the cellu neighbors =================================================================================================================
                //=============================================================================================================================================================================================
                if(bond_neighbors_cellu.count(neighbor_key) != 0){
    //                cout << "Found neighbor" << endl;
                    if(bond_neighbors_cellu[neighbor_key].outer_bond == false or par.mode_enzyme_size == -1){
    //                    cout << "No outer bond" << endl;
                    }
                    bond_neighbors_cellu[neighbor_key].remove_neighbor(x,y,z);
                    bond_neighbors_cellu[neighbor_key].outer_bond = bond_neighbors_cellu[neighbor_key].is_outer_bond();
                    if(bond_neighbors_cellu[neighbor_key].outer_bond == true){//g or par.mode_enzyme_size == -1){
       //                 cout << "outer bond!" << endl;
                        if(find_specific_bond(cellu, hemi, lign, bond_neighbors_cellu[neighbor_key].x, bond_neighbors_cellu[neighbor_key].y, bond_neighbors_cellu[neighbor_key].z, current_poly, current_bond, current_substrate) == true){
         //                   cout << "found hemi" << endl;
                            if(current_substrate == 1){
                                if(cellu[current_poly].status[current_bond] == -1){
                                    cellu[current_poly].status[current_bond] = 1;
                                    if(current_bond == 1){
                                        if(cellu[current_poly].len_poly > 1){
                                            if(CBH_enzyme_attached(par,current_poly,1) == false){
                                                if(par.free_poly_ends.find(cantor_pair_two(current_poly,1)) == par.free_poly_ends.end()){
                                                    par.free_poly_ends.insert({cantor_pair_two(current_poly,1),std::make_tuple(current_poly,1)});
                                                    par.N_free_ends++;
                                                    if(true){
                                                        addreaction(par, false,Table_cellu,current_poly,1,current_substrate,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),cellu[current_poly].crystalline[1]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else if(current_bond == cellu[current_poly].len_poly-2){
                                        if(cellu[current_poly].len_poly > 1){
                                            if(CBH_enzyme_attached(par,current_poly,cellu[current_poly].len_poly-2) == false){
                                                if(par.free_poly_ends.find(cantor_pair_two(current_poly,cellu[current_poly].len_poly-2)) == par.free_poly_ends.end()){
                                                    par.free_poly_ends.insert({cantor_pair_two(current_poly,cellu[current_poly].len_poly-2),std::make_tuple(current_poly,cellu[current_poly].len_poly-2)});
                                                    par.N_free_ends++;
                                                    if(true){
                                                        addreaction(par, false,Table_cellu,current_poly,cellu[current_poly].len_poly-2,current_substrate,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),cellu[current_poly].crystalline[cellu[current_poly].len_poly-2]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if(cellu[current_poly].len_poly == 1){
                                        reaction_type = 3;
                                        propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;
                                        if(propensite > 0){
                                            if(cellu[current_poly].crystalline[current_bond] == false){
                                                addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,propensite,0);
                                            }
                                            else{
                                                addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,crystal_propensite,1);
                                            }
                                        }
                                    }
                                    //CBH REACTION ADDING IS DONE IN FUNCTION CBH_ATTACHMENT()
/*                                    else if(cellu[current_poly].len_poly == 2){
                                        reaction_type = 2;
                                        propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;
                                        propensite = 0;
                                        if(propensite > 0){
                                           if(cellu[current_poly].crystalline[current_bond] == false){
                                                addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,propensite,enzyme_index);

                                            }
                                            else{
                                                addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,crystal_propensite,enzyme_index);

                                            }
                                        }
                                    }*/
                                    else if(cellu[current_poly].len_poly > 2){

                                        reaction_type = 1;
                                        propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                                        crystal_propensite = par.crystal_modifier_cellu * propensite;

                                        if(current_bond > 1 and current_bond < cellu[current_poly].len_poly-2){
                                            if(propensite > 0){
                                               if(cellu[current_poly].crystalline[current_bond] == false){
                                                    addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,propensite,0);
                                                }
                                                else{
                                                    addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,crystal_propensite,1);
                                                }
                                            }

                                        }
                                        else if(cellu[current_poly].len_poly == 4){
                                            /*
                                            if(current_bond == 1 or current_bond == 2){
                                                if(propensite > 0){
                                                   if(cellu[current_poly].crystalline[current_bond] == false){
                                                        addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,propensite,0);
                                                    }
                                                    else{
                                                        addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,crystal_propensite,1);
                                                    }
                                                }
                                            }
                                            */
                                        }
                                    //CBH REACTION ADDING IS DONE IN FUNCTION CBH_ATTACHMENT()
                                        /*if(current_bond == 1 or current_bond == cellu[current_poly].len_poly-2){
                                            reaction_type = 2;
                                            propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                                            crystal_propensite = par.crystal_modifier_cellu * propensite;
                                            propensite = 0;
                                            if(propensite > 0){
                                               if(cellu[current_poly].crystalline[current_bond] == false){
                                                    addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,propensite,enzyme_index);

                                                }
                                                else{
                                                    addreaction(par, false, Table_cellu,current_poly, current_bond,1, reaction_type,crystal_propensite,enzyme_index);

                                                }
                                            }
                                        }*/

                                    }

                                }
                            }
                            else if(current_substrate == 2){
                                cout << "In function update_reaction_tables(): current_substrate = 2 for a neighbor that should be contained in cellu. Stopping..." << endl;
                                exit(1);
                            }
                        }
                    }
                }
                //=============================================================================================================================================================================================
                //======================================== Next we check the hemi neighbors ===================================================================================================================
                //=============================================================================================================================================================================================
                else if(bond_neighbors_hemi.count(neighbor_key) != 0){
                    bond_neighbors_hemi[neighbor_key].remove_neighbor(x,y,z);
                    bond_neighbors_hemi[neighbor_key].outer_bond = bond_neighbors_hemi[neighbor_key].is_outer_bond();
                    if(bond_neighbors_hemi[neighbor_key].outer_bond == true){// or par.mode_enzyme_size == -1){
                        if(find_specific_bond(cellu, hemi, lign, bond_neighbors_hemi[neighbor_key].x, bond_neighbors_hemi[neighbor_key].y, bond_neighbors_hemi[neighbor_key].z, current_poly, current_bond, current_substrate) == true){
                            if(current_substrate == 1){
                                cout << "In function update_reaction_tables(): current_substrate = 1 for a neighbor that should be contained in hemi. Stopping..." << endl;
                                exit(1);
                            }
                            else if(current_substrate == 2){
                                if(hemi[current_poly].status[current_bond] == -1){
                                    hemi[current_poly].status[current_bond] = 1;
                                    reaction_type = 4;
                                    propensite = prop(par,reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                                    crystal_propensite = par.crystal_modifier_hemi * propensite;

                                    if(propensite > 0 and current_bond < hemi[current_poly].len_poly and current_bond >= 0){
                                       if(hemi[current_poly].crystalline[current_bond] == false){
                                            addreaction(par, false, Table_hemi,current_poly, current_bond,2, reaction_type,propensite,0);
                                        }
                                        else{
                                            addreaction(par, false, Table_hemi,current_poly, current_bond,2, reaction_type,crystal_propensite,1);
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
void update_reactiontables_inhib(vector<TList>& Table_cellu, vector<bList>& celluList, vector<TList>& Table_hemi, params& par, const int nbr_Glc_pdt, const int nbr_cellobiose, vector<double>& chem_entities){
/*
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
*/
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
bool find_specific_bond(const vector<bList>& celluList, const vector<bList>& hemiList, const vector<bList>& lignList, const double x, const double y, const int z, int& index_poly, int& z_index, int& material){
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
                EG_fraction+=*Table_cellu[i].liste_prop[j];
            else if(Table_cellu[i].indic_action[j] == 2)
                CBH_fraction+=*Table_cellu[i].liste_prop[j];
            else if(Table_cellu[i].indic_action[j] == 3)
                BGL_fraction+=*Table_cellu[i].liste_prop[j];
            else if(Table_cellu[i].indic_action[j] != 6){
                cout << "In function check_reaction_table_distribution: Table_cellu contains reactions not associated to EG, CBH or BGL" << endl;
            }
        }
    }
    for(int i=0;i<Table_hemi.size();i++){
        a0 += Table_hemi[i].prop_sum;
        for(int j=0;j<Table_hemi[i].nbr_element;j++){
            if(Table_hemi[i].indic_action[j] == 4)
                XYL_fraction += *Table_hemi[i].liste_prop[j];
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


//Returns the absolute value of the vector (dx,dy,dz)
double dist(double dx, double dy, double dz){
    return sqrt(dx*dx + dy*dy + dz*dz);
}


void lignin_glue(params& par, vector<bList>& cellu, vector<TList>& Table_cellu, vector<TList>& Table_hemi, TList& Table_lign, vector<double>& chem_entities, int nbr_poly_lign, int& nbr_Glc_pdt, int& nbr_cellobiose){
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
        else if(enzyme_number == 2){
            chem_entities[enzyme_number-1] -=1;// * chem_entities[enzyme_number]/double(par.init_CBH);
            par.propensities[enzyme_number-1] = prop(par,enzyme_number,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
            par.crystal_propensities[enzyme_number-1] = par.crystal_modifier_cellu * prop(par,enzyme_number,chem_entities,nbr_Glc_pdt,nbr_cellobiose);//Since we use pointers for the propensities, this is all we need to change
            par.N_enzymes_glued++;
            int glued_enzyme = int(drand48()*par.CBH_enzymes.size());
            if(par.CBH_enzymes[glued_enzyme].attached == true){
                deletereaction(Table_cellu[par.CBH_enzymes[glued_enzyme].poly_attached], Table_cellu, 1, par.CBH_enzymes[glued_enzyme].poly_attached, par.CBH_enzymes[glued_enzyme].bond_attached,2);
                if(par.free_poly_ends.find(cantor_pair_two(par.CBH_enzymes[glued_enzyme].poly_attached,par.CBH_enzymes[glued_enzyme].bond_attached)) == par.free_poly_ends.end() ){
                    par.free_poly_ends.insert({cantor_pair_two(par.CBH_enzymes[glued_enzyme].poly_attached,par.CBH_enzymes[glued_enzyme].bond_attached),make_tuple(par.CBH_enzymes[glued_enzyme].poly_attached,par.CBH_enzymes[glued_enzyme].bond_attached)});
                    par.N_free_ends++;
                    if(true){
                        addreaction(par, false,Table_cellu,par.CBH_enzymes[glued_enzyme].poly_attached,par.CBH_enzymes[glued_enzyme].bond_attached,1,6,prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose),-1);
                    }
                }
                par.CBH_enzymes.erase(par.CBH_enzymes.begin() + glued_enzyme);
            }
            else{
                par.CBH_enzymes.erase(par.CBH_enzymes.begin() + glued_enzyme);
                par.N_free_CBH--;
                par.propensities[5] = prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
                par.crystal_propensities[5] = par.crystal_modifier_cellu * prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose);//Since we use pointers for the propensities, this is all we need to change
/*                if(propensite > 0){
                    for(int i = 0; i<Table_cellu.size();i++){
                        for(int j=0; j<Table_cellu[i].nbr_element; j++){
                            if(Table_cellu[i].indic_action[j] == 6){
                                Table_cellu[i].liste_prop[j] = propensite;
                                Table_cellu[i].prop_uninhib[j] = propensite;
                            }
                        }
                    }
                }
                else{
                    for(int i = 0; i<Table_cellu.size();i++){
                        for(int j=0; j<Table_cellu[i].nbr_element; j++){
                            if(Table_cellu[i].indic_action[j] == 6){
                                delete_specific_reaction(Table_cellu,i,j);
                                j--;
                            }
                        }
                    }
                }*/
            }
        }
        else if(enzyme_number == 3){
            chem_entities[enzyme_number-1] -=1;// * chem_entities[enzyme_number]/double(par.init_BGL);
            par.propensities[enzyme_number-1] = prop(par,enzyme_number,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
            par.crystal_propensities[enzyme_number-1] = par.crystal_modifier_cellu * prop(par,enzyme_number,chem_entities,nbr_Glc_pdt,nbr_cellobiose);//Since we use pointers for the propensities, this is all we need to change
            par.N_enzymes_glued++;
        }
        else if(enzyme_number == 4){
            chem_entities[enzyme_number-1] -=1;// * chem_entities[enzyme_number]/double(par.init_XYL);
            par.propensities[enzyme_number-1] = prop(par,enzyme_number,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
            par.crystal_propensities[enzyme_number-1] = par.crystal_modifier_cellu * prop(par,enzyme_number,chem_entities,nbr_Glc_pdt,nbr_cellobiose);//Since we use pointers for the propensities, this is all we need to change
            par.N_enzymes_glued++;
        }



        if(chem_entities[enzyme_number-1] < 0)
            chem_entities[enzyme_number-1] = 0;

        if(enzyme_number != 4){
            for(int i=0;i<Table_cellu.size();i++){
                Table_cellu[i].calcTableProp();
            }
        }
        else{
            for(int i=0;i<Table_hemi.size();i++){
                Table_hemi[i].calcTableProp();
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

    int lignins_glued = int(par.lignols_blocked_per_enzyme * box_muller(0.25,1.));//int(0.5*4.*par.enzyme_radius) + int(drand48() * 0.5* 3. * double((4.*par.enzyme_radius)));// 2*enzyme_diameter /2nm <= N_glued < 8*enzyme_diameter / 2nm. 2nm is the assumed diameter of one monolignol
    while(lignins_glued < 1){
        lignins_glued = int(par.lignols_blocked_per_enzyme * box_muller(0.25,1.));
    }
    if(par.verbose == true){
        cout << "lignins glued: " << lignins_glued << endl;
    }
    if(lignins_glued > 0)
        par.nbr_lignin_blocked += lignins_glued;
    else{
        par.nbr_lignin_blocked += 1;
    }
    if (par.nbr_lignin_blocked >= par.nbr_monolignol){
        chem_entities.at(4) = 0;
    }
    else{
        chem_entities.at(4)=par.nbr_monolignol - par.nbr_lignin_blocked;//+ chem_entities[0] + chem_entities[1] + chem_entities[2] + chem_entities[3];
    }
    par.propensities[4] = prop(par,5,chem_entities,nbr_Glc_pdt,nbr_cellobiose);

    Table_lign.prop_sum = par.propensities[4];
//    cout << Table_lign.liste_prop[0] << endl;
//End of function
}


//Returns a normally distributed random number between 0 and 1 with mean value at 0.5, according to the box-muller algorithm
double box_muller(double sigma, double mu ){



    double z0 = 0;
    double rnum1 = drand48();
    double rnum2 = drand48();
//    double sigma = 0.5;
//    double mu = 0.5;


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

void hole_cutter(std::vector<bList>& polyList, std::vector<TList>& Table_hemi, params& par, int& difference_monomers, int& nbr_poly, int material){
//void hole_cutter(vector<bList>& celluList, vector<bList>& hemiList, vector<bList>& lignList, vector<TList>& Table_hemi, params& par, int& difference_hemi, int& difference_lign, int& nbr_poly_hemi, int& nbr_poly_lign, bool& error_bool, int material){

    if(difference_monomers == 0){
        return;
    }

    int max_hole_size = int(0.1*par.length_fibril);
    if(max_hole_size < 2*par.enzyme_radius)
        max_hole_size = 2*par.enzyme_radius;

    int N_mono_before = 0;//Stores the number of xylose molecules contained in the polymer before it was shortened
    int hole_size = 0;//Stores the hole sizes during each loop
    int hole_bottom = 0;//Stores the position of the hole during each loop
    int poly_selected = 0;//Stores the polymer index that is cut during each loop
    int new_poly_index = 0;//Index of the new poly

    int loop_count = 0;
    while(difference_monomers > 0){
        loop_count++;
        if(loop_count > 10000){
            cout << "In function hole_cutter(): too many loops. Stopping" << endl;
            exit(1);
        }

        hole_size = findmin(int(drand48()*max_hole_size), difference_monomers);

        do{
            poly_selected = int(drand48() * nbr_poly);
        }
        while(polyList[poly_selected].len_poly <= hole_size+2);

        do{
            hole_bottom = int(drand48() * (polyList[poly_selected].len_poly-1 - hole_size));
        }
        while(hole_bottom <= 0);


        N_mono_before = polyList[poly_selected].len_poly+1;


        //============================= Make a copy of the polymer that is shortened ====================================

        polyList.push_back(bList());
        nbr_poly++;
        new_poly_index = polyList.size()-1;


        polyList[new_poly_index].index= new_poly_index;//Specify the new index of the new polymer
        polyList[new_poly_index].x= polyList[poly_selected].x;
        polyList[new_poly_index].y= polyList[poly_selected].y;
        polyList[new_poly_index].set_z(polyList[poly_selected].z[0]);

        if(material == 2){
            Table_hemi.push_back(TList());//Make the new reaction list associated to the new polymer, if it is a hemicellulose polymer
            initTList(Table_hemi[new_poly_index], new_poly_index);//Initialize it
        }

        for(int i=0; i<polyList[poly_selected].len_poly; i++)//Make the new polymer into a copy of the old one before tayloring
        {
            addbond(polyList[new_poly_index],i,1,1,false);
            polyList[new_poly_index].z[i]=polyList[poly_selected].z[i];
            polyList[new_poly_index].status[i]=polyList[poly_selected].status[i];
            polyList[new_poly_index].bond_type[i]=polyList[poly_selected].bond_type[i];
            polyList[new_poly_index].N_blocked_positions[i] = polyList[poly_selected].N_blocked_positions[i];
            polyList[new_poly_index].crystalline[i] = polyList[poly_selected].crystalline[i];
        }
        //============================= Shorten the polymers ============================================================

        taylorOldPoly(polyList[poly_selected],hole_bottom);//Erase the end of the polymer, until beginning of hole

        taylorNewPoly(polyList[new_poly_index],hole_bottom+hole_size);//Erase the beginning of the polymer, until end of hole

        difference_monomers -= (N_mono_before - (polyList[poly_selected].len_poly + 1) - (polyList[new_poly_index].len_poly + 1));



    }
    if(difference_monomers < 0){
        difference_monomers = 0;
    }

}




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




void min_max_cellu_hemi_lign(const vector<bList>& celluList, params& par){
    par.mid_z = 0.5*par.length_fibril;
    par.mid_x = 0.;
    par.mid_y = 0.;
    for(int i=0;i<celluList.size();i++){
        par.mid_x += celluList[i].x;
        par.mid_y += celluList[i].y;
    }
    par.mid_x /= celluList.size();
    par.mid_y /= celluList.size();
    if(par.verbose == true){
        cout << "mid_x = " << par.mid_x << "; mid_y = " << par.mid_y << endl;
    }
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
void remove_neighbor_from_vector(const int x_pos, const int y_pos, const int z_pos, unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu, unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_hemi){


    std::tuple<double,double,int> index_bond = std::make_tuple(x_pos, y_pos, z_pos);
    std::tuple<double,double,int> index_neighbor = std::make_tuple(0.,0.,0);
//    char index_neighbor[3];


    if(bond_neighbors_cellu.find(index_bond) != bond_neighbors_cellu.end()){
        if(bond_neighbors_cellu[index_bond].N_neighbors == 0){
//            cout << "In function remove_neighbor_from_vector: a bond has no neighbors" << endl;
//            exit(1);
        }
        for(int j=0;j<bond_neighbors_cellu[index_bond].N_neighbors;j++){
            index_neighbor = std::make_tuple(bond_neighbors_cellu[index_bond].x_neighbors[j],bond_neighbors_cellu[index_bond].y_neighbors[j],bond_neighbors_cellu[index_bond].z_neighbors[j]);
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
            index_neighbor = std::make_tuple(bond_neighbors_hemi[index_bond].x_neighbors[j],bond_neighbors_hemi[index_bond].y_neighbors[j],bond_neighbors_hemi[index_bond].z_neighbors[j]);


            if(bond_neighbors_cellu.find(index_neighbor) != bond_neighbors_cellu.end()){
                bond_neighbors_cellu[index_neighbor].remove_neighbor(x_pos,y_pos,z_pos);
            }
            else if(bond_neighbors_hemi.find(index_neighbor) != bond_neighbors_hemi.end()){
                bond_neighbors_hemi[index_neighbor].remove_neighbor(x_pos,y_pos,z_pos);
            }
        }

        bond_neighbors_hemi.erase(index_bond);

    }


}

//Returns true, if a CBH enzyme is attached to the bond bond at the poly poly
bool CBH_enzyme_attached(params& par, int poly, int bond){
    for(int i=0; i<par.CBH_enzymes.size(); i++){
        if(par.CBH_enzymes[i].attached == true and par.CBH_enzymes[i].poly_attached == poly and par.CBH_enzymes[i].bond_attached == bond){
            return true;
        }
    }
    return false;
}


void update_attachment_reactions(params& par, std::vector<bList>& cellu, std::vector<TList>& Table_cellu, std::vector<double>& chem_entities, int& nbr_Glc_pdt, int& nbr_cellobiose){
    par.propensities[5] = prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
    par.crystal_propensities[5] = par.crystal_modifier_cellu * prop(par,6,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
    for(int i=0;i<Table_cellu.size();i++){
        Table_cellu[i].calcTableProp();
    }
    int reaction_type = 2;
    double propensite = prop(par, reaction_type,chem_entities,nbr_Glc_pdt,nbr_cellobiose);
    for(int i=0;i<par.CBH_enzymes.size();i++){
        if(par.CBH_enzymes[i].attached == true){
            deletereaction(Table_cellu[par.CBH_enzymes[i].poly_attached],Table_cellu, 1, par.CBH_enzymes[i].poly_attached, par.CBH_enzymes[i].bond_attached,6);
            //Check, if there is a missing digestion reaction for CBH
            addreaction(par, false, Table_cellu, par.CBH_enzymes[i].poly_attached, par.CBH_enzymes[i].bond_attached, 1, reaction_type, propensite,cellu[par.CBH_enzymes[i].poly_attached].crystalline[par.CBH_enzymes[i].bond_attached]);

        }
    }


//End of function
}


void distribute_lignin_covering(params& par,std::vector<bList>& lign){
    int N_covered_bonds = 0;//Determines the number of covered bonds
    int cover_position = 0;//Determines the positioning of the middle of the covering area along the polymer
    int covering_lower = 0;//Lower bound of covered area
    int covering_upper = 0;//Upper bound of covered area
    for(int i=0;i<lign.size();i++){

        N_covered_bonds = int(lign[i].len_poly * box_muller(par.sigma_lignin_covering, par.mu_lignin_covering));
        //int(lign[i].len_poly - lign[i].len_poly*cos(drand48()));
        if(N_covered_bonds <= 0){
            N_covered_bonds = 1;
        }
        else if(N_covered_bonds >= lign[i].len_poly){
            N_covered_bonds = lign[i].len_poly-1;
        }
        if(N_covered_bonds == lign[i].len_poly-1){
            for(int j=0;j<lign[i].len_poly;j++){
                lign[i].covering[j] = true;
            }
        }
        else{

            cover_position = int(drand48()*lign[i].len_poly);

            int count_covered = lign[i].len_poly;


            covering_lower = cover_position - int(0.5*N_covered_bonds-1);
            covering_upper = cover_position + int(0.5*N_covered_bonds);

            while(covering_lower < 0){
                covering_lower++;
                covering_upper++;
            }
            while(covering_upper >= lign[i].len_poly){
                covering_upper--;
            }


            int upper_index = lign[i].len_poly-1;
            int lower_index = 0;
            while(count_covered > N_covered_bonds){
                if(upper_index > covering_upper){
                    lign[i].covering[upper_index] = false;
                    count_covered--;
                    upper_index--;
                }
                if(count_covered > N_covered_bonds){
                    if(lower_index < covering_lower and lower_index != upper_index){
                        lign[i].covering[lower_index] = false;
                        count_covered--;
                        lower_index++;
                    }
                }
                if(upper_index < lower_index or (upper_index <= covering_upper and lower_index >= covering_lower))
                    break;
            }
        }

    }
    //End of function
}


//Checks whether action_mu1 = 6 reactions need to be blocked from occuring in this step
void check_for_blocked_ends(std::vector<bList>& cellu, std::vector<TList>& Table_cellu, params& par){
    bool test = false;//This checks whether the covered bool is true within a reaction
    int bond_selected = 0;
    int poly_attached = 0;
    int bond_attached = 0;
    double dx = 0.;
    double dy = 0.;
    double dz = 0.;
    for(int i=0; i<Table_cellu.size();i++){
        for(int j=0;j<Table_cellu[i].nbr_element;j++){
            if(Table_cellu[i].indic_action[j] !=2){
                test = false;
                bond_selected = Table_cellu[i].num_bond[j];
                for(int k=0;k<par.CBH_enzymes.size();k++){
                    if(par.CBH_enzymes[k].attached == true){
                        poly_attached = par.CBH_enzymes[k].poly_attached;
                        bond_attached = par.CBH_enzymes[k].bond_attached;
                        dx = cellu[i].x - cellu[poly_attached].x;
                        dy = cellu[i].y - cellu[poly_attached].y;
                        dz = cellu[i].z[bond_selected] - cellu[poly_attached].z[bond_attached];
                        if(dist(dx,dy,dz) < 2*par.enzyme_radius){
                            test = true;
                            break;
                        }
                    }
                }
                if(test == true){
                    Table_cellu[i].covered[j] = true;
                }
                else{
                    Table_cellu[i].covered[j] = false;
                }
            }
        }
        Table_cellu[i].calcTableProp();
    }
    //End of function
}
