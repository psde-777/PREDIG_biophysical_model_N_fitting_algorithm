#ifndef FUNCTIONS_HPP_
#define FUNCTIONS_HPP_

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
#include "bList.hpp"
#include "TList.hpp"
#include "DPList.hpp"
#include "neighborList.hpp"
#include "CBH_enzyme.hpp"
#include "params.hpp"
#include "tuple_hash.hpp"
#include <unordered_map>


/****************************Initiation of Structures*******************/

/**
 * Initiates the Table list TODO: move this into the constructor of TList
 *
 * @param list The reaction table
 * @param origin The index of the polymer to which the reaction table belongs
 *
 * @return void
 * **/
void initTList(
	TList& list,
 	int origin
 );


/****************************Functions*******************/



/**
 * Adds a bond to the polymer "list"
 *
 * @param list, The bList polymer object to which the bond will be added
 * @param position The z-coordinate of the new bond
 * @param state: The digestion state of the new bond Indicates -2 if the bond is acetylated or lignin,
 *  -1 if the bond is not yet accessible,
 *   1 if the bond is accessible,
 *   0 if the polymer is fully digested
 * @param bond_type, The type of the new bond
 * @param crystalline); The crystallinity state of the bond
 *
 * @return void
 *
**/
void addbond(
	bList& list,
	int position,
	int state,
	int bond_type,
	bool crystalline);


/**
 * Erases the end of the polymer "list", until digested bond
 *
 * @param list, The bList polymer object whose end will be erased
 * @param bond_num The bond up to which the polymer is erased
 *
 * @return void
 *
**/
void taylorOldPoly(
	bList& list,
    int bond_num);


/**
 * Erases the begining of the polymer "list", from digested bond
 *
 * @param list, The bList polymer object whose end will be erased
 * @param bond_num The bond up to which the polymer is erased
 *
 * @return void
 **/
void taylorNewPoly(
	bList& list,
    int bond_num);


/**
 * Adds a new reaction to the reaction table "list".
 *
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param confidence If you are sure that this reaction is not already contained in the table
 *  set confidence to true
 *  otherwise set it to false
 * @param list The vector containing all reaction tables of the current material
 * @param poly the index of the reaction table in list
 * @param bond_num The bond number in the polymer to which the reaction belongs
 * @param mat The material of the current polymer
 * @param act2 The type of reaction. 1 for EG, 2 for CBH, 3 for BGL, 4 for XYL, 5 for lignin adhesion, 6 for CBH attachment
 * @param prope The propensity of the reaction
 * @param crystalline_bool The crystallinity status of the bond in question
 *
 * @return void
 **/
void addreaction(
	params& par,
    bool confidence,
    std::vector<TList>& list,
    int poly,
    int bond_num,
    int mat,
    int act2,
    double prope,
    bool crystalline_bool);


/**
 * Clears the table with index table from all reactions
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param allTables The vector containing all reaction table objects
 * @param table The index of the reaction table in allTables to be cleared
 *
 * @return void
 **/
void clear_table(
	params& par,
    std::vector<TList>& allTables,
    int table);


/**
 * Fills the reaction table of a newly created polymer
 *
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param allTables The vector containing all reaction table objects of the same type (cellu, hemi or lign)
 * @param polyList The vector containing all polymer objects
 * @param poly_selected The index of the selected polymer both in allTables and polyList
 * @param substrate The type of the polymer
 * @param error for error checking outside of this function. This is set to true if an error is detected within this function
 * @param chem_entities The vector containing the enzyme concentrations
 * @param nbr_Glc_pdt Number of glucose produced
 * @param nbr_cellobiose Number of cellobiose in the system
 *
 * @return void
 **/
void fill_table(
	params& par,
    std::vector<TList>& allTables,
    std::vector<bList>& polyList,
    int poly_selected,
    int substrate,
    bool& error,
    std::vector<double>& chem_entities, int& nbr_Glc_pdt, int& nbr_cellobiose);



/**
 * Deletes the reaction labeled "action_target" for the polymer "poly" for its bond "bond_target" of the material "mat" in the list of reactions
 * @param list The reaction table from which the reaction will be deleted
 * @param Table The vector of all reaction tables
 * @param mat The material of the polymer to which the reaction table belongs
 * @param poly The index of the polymer from which the reactions will be deleted
 * @param bond_target The bond from which the reaction will be deleted
 * @param action_target The type of the reaction; 1 for EG, 2 for CBH, 3 for BGL, 4 for XYL, 5 for lignin adhesion, 6 for CBH attachment
 *
 * @return void
 **/
void deletereaction(
	TList& list,std::vector<TList>& Table,int mat,int poly,int bond_target,int action_target);



/**
 * Deletes a reaction of known index in a reaction table
 * @param Table The vector containing all reaction tables
 * @param table_selected The index of the reaction table in Table
 * @param i The reaction to be deleted
 *
 * @return void
 **/
void delete_specific_reaction(
	std::vector<TList>& Table,
    int table_selected,
    int i);



/**
 * Deletes ALL reactions for the bond "bond_target" from the polymer "poly" of the material "mat" in the list of reactions
 * @param Table The vector containing all reaction tables
 * @param mat The material of the polymer to which the reaction table belongs
 * @param poly The index of the polymer from which the reactions will be deleted
 * @param bond_target The bond from which the reaction(s) will be deleted
 *
 * @return void
 **/
void deleteAllreaction(
	TList& list,std::vector<TList>& Table,int mat,int poly,int bond_target);



/**
 * Computes the propensity for a reaction to take place, depending on its type "act", and the corresponding input parameters
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param act The type of the reaction
 * @param chem_entities The vector containing all enzyme concentrations
 * @param nbr_Glc_pdt Number of glucose produced
 * @param nbr_cellobiose Number of cellobiose in the system
 *
 * @return The propensity of the reaction
 **/
double prop(
	params& par,
    double act,
    std::vector<double> chem_entities,int& nbr_Glc_pdt, int& nbr_cellobiose); //partho changed added inhib in prop function


/**
 * Counts the amount of glucose in the bList vector list
 * @param list The vector containing all polymers
 * @param substrate The type of polymer (cellu or hemi)
 *
 * @return Overall number of glucose monomers contained in list
 **/
double countGlc(
	const std::vector<bList>& list,
    const int substrate);


/**
 * Counts the amount of xylose in the bList vector list
 * @param list The vector containing all polymers
 * @param substrate The type of polymer (cellu or hemi)
 *
 * @return Overall number of xylose monomers contained in list
 **/
double countXyl(
	const std::vector<bList>& list);


/**
 * Counts the amount of monolignols in the bList vector list
 * @param list The vector containing all lignin polymers
 *
 * @return Overall number of monolignols contained in list
 **/
double countLign(
	const std::vector<bList>& list);


/**
 * TODO: Update this function to include checking of hemicellulose, e.g. if we have MLGs
 * Counts the amount of cellobiose, i.e. polymers of length 1 (we count in bonds!)
 * @param list The vector containing all polymers
 *
 * @return Overall number of cellobiose, contained in list
 **/
int countCellobiose(
	const std::vector<bList>& list);


/**
 * Counts reactions inside a reaction table, which are associated to the polymer poly_selected
 * @param list The reaction table to be counted (NOT the vector of all reaction tables)
 * @param poly_selected The index of the polymer TODO: this is possibly redundant
 * @param error Set to true if an error is detected in this function. Used for errorhandling outside
 *
 * @return The number of reactions in the reaction table
 **/
int countReactions(
	TList& list,
    int poly_selected,
    bool& error);


/**
 * Prints a visual representation of all polymers contained in list to the terminal
 * @param list The reaction table to be counted (NOT the vector of all reaction tables)
 *
 * @return void
 **/
void drawPolys(
	std::vector<bList>& list);



/**
 * Finds the reaction chosen for each gillespie step
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param Table_cellu The vector containing all reaction tables for cellulose
 * @param Table_hemi The vector containing all reaction tables for hemicellulose
 * @param Table_lign The vector containing all reaction tables for lignin
 * @param cellu The vector containing all cellu polymers
 * @param hemi The vector containing all hemi polymers
 * @param lign The vector containing all lign polymers
 * @param a0 The propensity value picked between 0 and the sum of all propensities
 * @param test
 * @param tau The time calculated for the gillespie algorithm
 * @param mu1 The index of the picked reaction in the reaction table
 * @param table_selected The selected reaction table
 * @param poly_selected The selected poly TODO: This and table_selected is redundant
 * @param bond_selected The bond on which the picked reaction acts
 * @param action_mu1 The type of the picked reaction
 * @param substrate The type of substrate on which the reaction acts
 *
 * @return void
 **/
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
	int& substrate);



/**
 * Returns the index of the reaction chosen for each Gillespie step from the previously selected reaction table list (this function is only called in the function pick_reaction())
 * @param list The reaction table to be searched
 * @param suma The sum of the propensities of previously searched tables
 * @param test The propensity sum within the reaction table (increases as we search through the table)
 *
 * @return index of the reaction
 **/
int findIndex(
	const TList& list,
    double suma,
    double test);






/**
 * Carries out digestions via endoglucanase
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param allTables The vector containing all reaction tables
 * @param polyList The vector containing all polymers of the substrate which is digested
 * @param nbr_poly The number of polymers of the given substrate
 * @param nbr_xyl_pdt Overall number of released xylose molecules
 * @param nbr_Glc_pdt Overall number of released glucose molecules
 * @param bond_selected Bond to be digested
 * @param poly_selected Polymer in which the bond to be digested is located
 * @param len_polyLoopStart The length of the polymer at the beginning of the function (for error checking)
 * @param nbr_cellobiose Overall nummber of cellobiose "poly"mers in the system
 * @param chem_entities The vector containing the enzyme concentrations
 * @param substrate The substrate
 * @param bond_neighbors_cellu The unordered_map containing all neighbors of each bond in the cellu polymers
 * @param bond_neighbors_hemi The unordered_map containing all neighbors of each bond in the hemi polymers
 * @param verbose Parameter for switching on verbose output
 * @param error Used for error checking. Set to true if a tracked error is detected
 *
 * @return void
 **/
void EG_digest(
	params& par,
    std::vector<TList>& allTables,
    std::vector<bList>& polyList,
    int& nbr_poly,
    int& nbr_xyl_pdt,
    int& nbr_Glc_pdt,
    const int bond_selected,
    const int poly_selected,
    int len_polyLoopStart,
    int& nbr_cellobiose,
    std::vector<double>& chem_entities,
    const int substrate,
    std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>,neighborList>& bond_neighbors_hemi,
    const bool verbose,
    bool& error);


/**
 * Carries out digestions via cellobiohydrolase
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param allTables The vector containing all reaction tables
 * @param polyList The vector containing all polymers of the substrate which is digested
 * @param nbr_poly The number of polymers of the given substrate
 * @param nbr_xyl_pdt Overall number of released xylose molecules
 * @param nbr_Glc_pdt Overall number of released glucose molecules
 * @param bond_selected Bond to be digested
 * @param poly_selected Polymer in which the bond to be digested is located
 * @param len_polyLoopStart The length of the polymer at the beginning of the function (for error checking)
 * @param nbr_cellobiose Overall nummber of cellobiose "poly"mers in the system
 * @param chem_entities The vector containing the enzyme concentrations
 * @param substrate The substrate
 * @param bond_neighbors_cellu The unordered_map containing all neighbors of each bond in the cellu polymers
 * @param bond_neighbors_hemi The unordered_map containing all neighbors of each bond in the hemi polymers
 * @param verbose Parameter for switching on verbose output
 * @param error Used for error checking. Set to true if a tracked error is detected
 *
 * @return void
 **/
void CBH_digest(
	params& par,
    std::vector<TList>& allTables,
    std::vector<bList>& polyList,
    int& nbr_poly,
    int& nbr_xyl_pdt,
    int& nbr_Glc_pdt,
    const int bond_selected,
    const int poly_selected,
    int len_polyLoopStart,
    int& nbr_cellobiose,
    std::vector<double>& chem_entities,
    const int substrate,
    std::unordered_map<std::tuple<double,double,int>,neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>,neighborList>& bond_neighbors_hemi,
    const bool verbose,
    bool& error);


/**
 * Carries out digestions via betaglucosidase
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param allTables The vector containing all reaction tables
 * @param polyList The vector containing all polymers of the substrate which is digested
 * @param nbr_poly The number of polymers of the given substrate
 * @param nbr_xyl_pdt Overall number of released xylose molecules
 * @param nbr_Glc_pdt Overall number of released glucose molecules
 * @param bond_selected Bond to be digested
 * @param poly_selected Polymer in which the bond to be digested is located
 * @param len_polyLoopStart The length of the polymer at the beginning of the function (for error checking)
 * @param nbr_cellobiose Overall nummber of cellobiose "poly"mers in the system
 * @param chem_entities The vector containing the enzyme concentrations
 * @param substrate The substrate
 * @param bond_neighbors_cellu The unordered_map containing all neighbors of each bond in the cellu polymers
 * @param bond_neighbors_hemi The unordered_map containing all neighbors of each bond in the hemi polymers
 * @param verbose Parameter for switching on verbose output
 * @param error Used for error checking. Set to true if a tracked error is detected
 *
 * @return void
 **/
void BGL_digest(
	params& par,
    std::vector<TList>& allTables,
    std::vector<bList>& polyList,
    int& nbr_poly,
    int& nbr_xyl_pdt,
    int& nbr_Glc_pdt,
    const int bond_selected,
    const int poly_selected,
    int len_polyLoopStart,
    int& nbr_cellobiose,std::vector<double>& chem_entities,
    const int substrate,
    std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_hemi,
    const bool verbose,
    bool& error);


/**
 * Carries out digestions via xylanase
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param allTables The vector containing all reaction tables of hemi
 * @param allCelluTables The vector containing all reaction tables of cellu
 * @param polyList The vector containing all hemi polymers
 * @param celluList The vector containing all cellu polymers
 * @param nbr_poly The number of hemi polymers
 * @param nbr_xyl_pdt Overall number of released xylose molecules
 * @param nbr_Glc_pdt Overall number of released glucose molecules
 * @param bond_selected Bond to be digested
 * @param poly_selected Polymer in which the bond to be digested is located
 * @param len_polyLoopStart The length of the polymer at the beginning of the function (for error checking)
 * @param nbr_cellobiose Overall nummber of cellobiose "poly"mers in the system
 * @param chem_entities The vector containing the enzyme concentrations
 * @param substrate The substrate
 * @param bond_neighbors_cellu The unordered_map containing all neighbors of each bond in the cellu polymers
 * @param bond_neighbors_hemi The unordered_map containing all neighbors of each bond in the hemi polymers
 * @param verbose Parameter for switching on verbose output
 * @param error Used for error checking. Set to true if a tracked error is detected
 *
 * @return void
 **/
void XYL_digest(
	params& par,
    std::vector<TList>& allTables,
    std::vector<TList>& allCelluTables,
    std::vector<bList>& polyList,
    std::vector<bList>& celluList,
    int& nbr_poly,
    int& nbr_xyl_pdt,
    int& nbr_Glc_pdt,
    const int bond_selected,
    const int poly_selected,
    int len_polyLoopStart,
    int& nbr_cellobiose,std::vector<double>& chem_entities,
    const int substrate,
    std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_hemi,
    const bool verbose,
    bool& error);



/**
 * Attaches a CBH to the given bond of the given poly and adds the digestion reaction
  * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param allTables The vector containing all reaction tables of hemi
 * @param celluList The vector containing all cellu polymers
 * @param hemiList The vector containing all hemi polymers
 * @param lignList The vector containing all lign polymers
 * @param nbr_poly The number of polymers of the given substrate
 * @param bond_selected Bond to be digested
 * @param poly_selected Polymer in which the bond to be digested is located
 * @param chem_entities The vector containing the enzyme concentrations
 * @param substrate The substrate
 * @param bond_neighbors_cellu The unordered_map containing all neighbors of each bond in the cellu polymers
 * @param bond_neighbors_hemi The unordered_map containing all neighbors of each bond in the hemi polymers
 * @param nbr_Glc_pdt Number of glucose produced
 * @param nbr_cellobiose Number of cellobiose in the system
 *
 * @return void
 **/
void CBH_attachment(
	params& par,
    std::vector<TList>& allTables,
    std::vector<bList>& celluList,
    std::vector<bList>& hemiList,
    std::vector<bList>& lignList,
    int& nbr_poly,
    const int bond_selected,
    const int poly_selected,
    std::vector<double>& chem_entities,
    const int substrate,
    std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_hemi, int& nbr_Glc_pdt, int& nbr_cellobiose);


/**
 * Updates the reaction tables after a digestion has been carried out
 * @param cellu The vector containing all cellu polymers
 * @param hemi The vector containing all hemi polymers
 * @param lign The vector containing all lign polymers
 * @param Table_cellu The vector containing all cellu reaction tables
 * @param Table_hemi The vector containing all hemi reaction tables
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param chem_entities The vector containing all enzyme concentrations
 * @param bond_neighbors_cellu The vector containing all neighbors of each cellu bond
 * @param bon_neighbors_hemi The vector containing all neighbors of each hemi bond
 * @param x The x coordinate of the previously carried out reaction
 * @param y The y coordinate of the previously carried out reaction
 * @param z The z coordinate of the previously carried out reaction
 * @param substrate The current substrate
 * @param poly_selected The polymer on which the previous reaction took place
 * @param bond_selected The bond on which the previous reaction took place
 * @param nbr_Glc_pdt Number of glucose produced
 * @param nbr_cellobiose Number of cellobiose in the system
*
 * @return void
 **/
void update_reactiontables(
	std::vector<bList>& cellu,
    std::vector<bList>& hemi,
    std::vector<bList>& lign,
    std::vector<TList>& Table_cellu,
    std::vector<TList>& Table_hemi,
    params& par,
    std::vector<double>& chem_entities,
    std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_hemi,
    double x,
    double y,
    int z,
    int substrate,
    int poly_selected,
    int bond_selected, int& nbr_Glc_pdt, int& nbr_cellobiose);



/**
 * Updates the reaction tables regarding inhibition through glucose or cellobiose
 * @param Table_cellu The vector containing all cellu reaction tables
 * @param celluList The vector containing all cellu polymers
 * @param Table_hemi The vector containing all hemi reaction tables
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param nbr_Glc_pdt The overall number of released glucose molecules in the system
 * @param nbr_cellobiose The number of cellobiose "poly"mers in the system
 * @param chem_entities The vector containing the enzyme concentrations
 *
 * @return void
 **/
void update_reactiontables_inhib(
	std::vector<TList>& Table_cellu,
    std::vector<bList>& celluList,
    std::vector<TList>& Table_hemi,
    params& par,
    const int nbr_Glc_pdt,
    const int nbr_cellobiose,
    std::vector<double>& chem_entities);


/**
 * Counts the number of chains for different degrees of polymerization
 * @param polyList The vector containing all cellu polymers
 * @param chain_distrib The vector containing the DP distribution for all timesteps
 * @param length_fibril The lenght of the microfibril
 * @param real_time The real time as calculated by the gillespie algorithm
 * @param t0 The step number of the gillespie algorithm
 * @param index The index of the chain_distrib vector to be filled
 *
 * @return void
 **/
void chain_length_distribution(
	std::vector<bList>& polyList,
    std::vector<DPList>& chain_distrib,
    const int length_fibril,
    const double real_time,
    const int t0,
    const int index);


/**
 * Checks whether all bonds with status = 1 have reactions associated to them
 * @param allTables The vector containing all reaction tables of the given substrate
 * @param polyList The vector containing all polys of the given substrate
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 *
 * @return true if all bonds with status = 1 have reactions, and false otherwise
 **/
bool check_reactions(
	std::vector<TList>& allTables,
    std::vector<bList>& polyList,
    params& par);



/**
 * Returns the bond type according to the type of left bond l and right bond r. Used for bond_type distribution on hemicellulose
 * @param l The type of the bond to the left
 * @param r The type of the bond to the right
 *
 * @return The bond type
 **/
int find_bond_type(
	int l,
    int r);


//
/**
 * Returns true if the bond exists within the polyLists, and fills the reference arguments with the polymer index and the material. If it returns false, the reference arguments should not be used further and are set to -1
 * @param celluList The vector containing all cellu polymers
 * @param hemiList The vector containing all hemi polymers
 * @param lignList The vector containing all lign polymers
 * @param x The x coordinate of the bond to be found
 * @param y The y coordinate of the bond to be found
 * @param z The z coordinate of the bond to be found
 * @param index_poly The index of the poly. Is set if the bond is found
 * @param z_index The index of the bond within the polymer
 * @param material The type of material in which the bond is located
 *
 * @return true of the bond exists within the polylists, and false otherwise
 **/
bool find_specific_bond(
	const std::vector<bList>& celluList,
    const std::vector<bList>& hemiList,
    const std::vector<bList>& lignList,
    const double x,
    const double y,
    const int z,
    int& index_poly,
    int& z_index,
    int& material);


/**
 * Looks through all reaction tables and looks at the relative fraction of reactions associated to each enzyme
 * @param Table_cellu The vector containing all cellu reaction tables
 * @param Table_hemi The vector containing all hemi reaction tables
 * @param Table_lign The vector containing all lign reaction tables
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 *
 * @return void
 **/
void check_reaction_table_distribution(
	const std::vector<TList>& Table_cellu,
    const std::vector<TList>& Table_hemi,
    TList& Table_lign,
    params& par);


/**
 * Returns the smaller of a and b
 * @param a The first number
 * @param b The second number
 *
 * @return minimum
 **/
int findmin(
	int a,
    int b);


/**
 * Returns the larger of a and b
 * @param a The first number
 * @param b The second number
 *
 * @return maximum
 **/
int findmax(
	int a,
    int b);



/**
 * Returns the absolute value of the vector (dx,dy,dz)
 * @param dx The x coordinate of the vector
 * @param dy The y coordinate of the vector
 * @param dz The z coordinate of the vector
 *
 * @returns the absolute value
 **/
double dist(double dx, double dy, double dz);


/**
 * Adjusts propensities of all tables by "glueing enzymes"
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param cellu The vector containing all cellu polymers
 * @param Table_cellu The vector containing all cellu reaction tables
 * @param Table_hemi The vector containing all hemi reaction tables
 * @param Table_lign The lign reaction table
 * @param chem_entities The vector containing all enzyme concentrations
 * @param nbr_poly_lign The number of lignin polymers
 * @param nbr_Glc_pdt Number of glucose produced
 * @param nbr_cellobiose Number of cellobiose in the system
 *
 * @return void
 **/
void lignin_glue(
	params& par,
    std::vector<bList>& cellu,
    std::vector<TList>& Table_cellu,
    std::vector<TList>& Table_hemi,
    TList& Table_lign,
     std::vector<double>& chem_entities,
    int nbr_poly_lign, int& nbr_Glc_pdt, int& nbr_cellobiose);


//Returns a normally distributed random number between 0 and 2*mu, according to the box-muller algorithm
/**
 * Clears the table with index table from all reactions
 * @param sigma The standard deviation of the normal distribution from which we pick
 * @param mu The mean of the normal distribution from which we pick
 *
 * @return random number between 0 and 2*mu
 **/
double box_muller(
	double sigma,
    double mu);



/**
 * Returns the index of the polymer of coordinates x,y, which also contains a bond at position z
 * @param Table The vector of all polymers to be searched
 * @param x The x-coordinate of the bond to be found
 * @param y The y-coordinate of the bond to be found
 * @param z The z-coordinate of the bond to be found
 *
 * @return the index of the found polymer
 **/
int find_poly(
	const std::vector<bList>& Table,
    const int x,
    const int y,
    const int z);


/**
 * Returns the index of the bond within a poly whose z coordinate is z_selected
 * @param poly The polymer to be searched
 * @param z_selected The z coordinate of the bond
 *
 * @return the index of the bond
 **/
int find_bond(
	const bList& poly,
    int z_selected);




/**
 * Simple hole cutter function. Cuts 1D holes
 * @param polyList The vector containing all polymers to be cut
 * @param Table_hemi The vector containing all reaction tables fo the polymers to be cut
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param difference_monomers The overall number of monomers to be removed
 * @param nbr_poly The number of polymers
 * @param material The material that is being cut
 *
 * @return void
 **/
void hole_cutter(
                 std::vector<bList>& polyList,
                 std::vector<TList>& Table_hemi,
                 params& par,
                 int& difference_monomers,
                 int& nbr_poly,
                 int material);



/**
 * Returns true if a bond exists at the specified x-, y- and z coordinates
 * @param celluList The vector containing all cellu polymers
 * @param hemiList The vector containing all hemi polymers
 * @param lignList The vector containing all lign polymers
 * @param x_pos The x-coordinate of the bond
 * @param y_pos The y-coordinate of the bond
 * @param z_pos The z-coordinate of the bond
 *
 * @return true if bond exists, false otherwise
 **/
bool bond_exists(
	const std::vector<bList>& celluList,
    const std::vector<bList>& hemiList,
    const std::vector<bList>& lignList,
    int x_pos,
    int y_pos,
    int z_pos);


/**
 * Returns true if x_pos and y_pos are found as a valid combination inside celluList, hemiList or lignList
 * @param celluList The vector containing all cellu polymers
 * @param hemiList The vector containing all hemi polymers
 * @param lignList The vector containing all lign polymers
 * @param x_pos The x-coordinate of the bond
 * @param y_pos The y-coordinate of the bond
 *
 * @return true if combination is valid, false otherwise
 **/
bool is_valid_pos(
	const std::vector<bList>& celluList,
    const std::vector<bList>& hemiList,
    const std::vector<bList>& lignList,
    int x_pos,
    int y_pos);






/**
 * Sets min_x,max_x,min_y,max_y for hemi and lign
 * @param celluList The vector containing all cellu polymers
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 **/
void min_max_cellu_hemi_lign(
	const std::vector<bList>& celluList,
    params& par);


/**
 * Returns an integer which is unique for every set (a,b,c), with a,b,c > 0
 * @param a The first number
 * @param b The second number
 * @param c The third number
 *
 * @return unique integer
 **/
int cantor_pair_three(
	int a,
    int b,
    int c);


/**
 * Returns an integer which is unique for every set (a,b), with a,b > 0
 * @param a The first number
 * @param b The second number
 *
 * @return unique integer
 **/
int cantor_pair_two(
	int a,
    int b);



/**
 * Removes a bond from the neighbor vector
 * @param x_pos The x-coordinate of the bond
 * @param y_pos The y-coordinate of the bond
 * @param z_pos The z-coordinate of the bond
 * @param bond_neighbors_cellu The unordered map containing all neighbors of each cellu bond
 * @param bond_neighbors_hemi The unordered map containing all neighbors of each hemi bond
 *
 * @return void
 **/
void remove_neighbor_from_vector(
	const int x_pos,
    const int y_pos,
    const int z_pos,
    std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_hemi);



/**
 * Returns true, if a CBH enzyme is attached to the bond bond at the poly poly, and false otherwise
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param poly The polymer to be checked
 * @param bond The bond to be checked
 *
 * @return true if enzyme is attached, false otherwise
 **/
bool CBH_enzyme_attached(
	params& par,
    int poly,
    int bond);


//
/**
 * Updates CBH attachment reactions if there is a change in the number of free CBH enzymes
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param cellu The vector containing all cellu tables
 * @param Table_cellu The vector containing all cellu reaction tables
 * @param chem_entities The vector containing all enzyme concentrations
 * @param nbr_Glc_pdt Number of glucose produced
 * @param nbr_cellobiose Number of cellobiose in the system
 *
 * TODO: Update this function to also check hemi, if we have mixed-linkage glucans
 *
 * @return void
 **/
void update_attachment_reactions(
	params& par,
    std::vector<bList>& cellu,
    std::vector<TList>& Table_cellu,
    std::vector<double>& chem_entities, int& nbr_Glc_pdt, int& nbr_cellobiose);



/**
 * Sets the covering parameters of each lignin polymer. Only called once in the beginning
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 * @param lign The vector containing all lignin polymers
 *
 * @return void
 **/
void distribute_lignin_covering(
	params& par,
	std::vector<bList>& lign);



/**
 * Checks whether action_mu1 = 6 reactions need to be blocked from occuring in this step
 * @param cellu The vector containing all cellu polymers
 * @param Table_cellu The vector containing all cellu reaction tables
 * @param par The Param struct in which most parameters are stored. Should ALWAYS be passed as a reference!
 *
 * @return void
 **/
void check_for_blocked_ends(std::vector<bList>& cellu,
                            std::vector<TList>& Table_cellu,
                            params& par);

#endif
