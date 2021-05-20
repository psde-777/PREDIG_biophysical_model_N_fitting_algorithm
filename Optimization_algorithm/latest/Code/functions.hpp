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
//Initiates the Table list
void initTList(
	TList& list,
 	int origin
 );


/****************************Functions*******************/


//Adds a bond to the polymer "list"
void addbond(
	bList& list,
	int position,
	int state,
	int bond_type,
	bool crystalline);

//Erases the end of the polymer "list", until digested bond
void taylorOldPoly(
	bList& list,
    int bond_num);


//Erases the begining of the polymer "list", from digested bond
void taylorNewPoly(
	bList& list,
    int bond_num);


//Adds a new reaction to the reaction table "list". If you are sure that this reaction is not already contained in the table, set confidence to true, otherwise set it to false
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


//Clears the table with index table from all reactions
void clear_table(
	params& par,
    std::vector<TList>& allTables,
    int table);


//Fills the reaction table of a newly created polymer
void fill_table(
	params& par,
    std::vector<TList>& allTables,
    std::vector<bList>& polyList,
    int poly_selected,
    int substrate,
    bool& error,
    std::vector<double>& chem_entities);


//Deletes the reaction labeled "action_target" for the polymer "poly" for its bond "bond_target" of the material "mat" in the list of reactions
void deletereaction(
	TList& list,std::vector<TList>& Table,int mat,int poly,int bond_target,int action_target);


//Deletes a reaction of known index in a reaction table
void delete_specific_reaction(
	std::vector<TList>& Table,
    int table_selected,
    int i);


//Deletes ALL reactions for the bond "bond_target" from the polymer "poly" of the material "mat" in the list of reactions
void deleteAllreaction(
	TList& list,std::vector<TList>& Table,int mat,int poly,int bond_target);


//Computes the propensity for a reaction to take place, depending on its type "act", and the corresponding input parameters
double prop(
	params& par,
    double act,
    std::vector<double> chem_entities);


//Counts the amount of glucose in the bList vector list
double countGlc(
	const std::vector<bList>& list,
    const int substrate);


//Counts the amount of xylose in the bList vector list
double countXyl(
	const std::vector<bList>& list);


//Counts the amount of monolignols in the bList vector list
double countLign(
	const std::vector<bList>& list);


//Counts the amount of cellobiose in the bList vector list
int countCellobiose(
	const std::vector<bList>& list);


//Counts reactions inside a reaction table, which are associated to the polymer poly_selected
int countReactions(
	TList& list,
    int poly_selected,
    bool& error);


//Prints a visual representation of all polymers contained in list to the terminal
void drawPolys(
	std::vector<bList>& list);



//Finds the reaction chosen for each gillespie step
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


//Returns the index of the reaction chosen for each Gillespie step from the previously selected reaction table list (this function is only called in the function pick_reaction())
int findIndex(
	const TList& list,
    double suma,
    double test);





//Carries out digestions via endoglucanase
void EG_digest(
	params& par,
    std::vector<TList>& allTables,
    std::vector<bList>& polyList,
    int& nbr_poly,
    int& nbr_xyl_pdt,
    int& nbr_Glc_pdt,const int bond_selected,const int poly_selected,
    int len_polyLoopStart,
    int& nbr_cellobiose,std::vector<double>& chem_entities,
    const int substrate,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_hemi,
    const bool verbose,
    bool& error);


//Carries out digestions via cellobiohydrolase
void CBH_digest(
	params& par,
    std::vector<TList>& allTables,
    std::vector<bList>& polyList,
    int& nbr_poly,
    int& nbr_xyl_pdt,
    int& nbr_Glc_pdt,const int bond_selected,const int poly_selected,
    int len_polyLoopStart,
    int& nbr_cellobiose,std::vector<double>& chem_entities,
    const int substrate,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_hemi,
    const bool verbose,
    bool& error);


//Carries out digestions via betaglucosidase
void BGL_digest(
	params& par,
    std::vector<TList>& allTables,
    std::vector<bList>& polyList,
    int& nbr_poly,
    int& nbr_xyl_pdt,
    int& nbr_Glc_pdt,const int bond_selected,const int poly_selected,
    int len_polyLoopStart,
    int& nbr_cellobiose,std::vector<double>& chem_entities,
    const int substrate,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_hemi,
    const bool verbose,
    bool& error);


//Carries out digestions via xylanase
void XYL_digest(
	params& par,
    std::vector<TList>& allTables,
    std::vector<TList>& allCelluTables,
    std::vector<bList>& polyList,
    std::vector<bList>& celluList,
    int& nbr_poly,
    int& nbr_xyl_pdt,
    int& nbr_Glc_pdt,const int bond_selected,const int poly_selected,
    int len_polyLoopStart,
    int& nbr_cellobiose,std::vector<double>& chem_entities,
    const int substrate,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_hemi,
    const bool verbose,
    bool& error);


//Attaches a CBH to the given bond of the given poly and adds the digestion reaction
void CBH_attachment(
	params& par,
    std::vector<TList>& allTables,
    std::vector<bList>& celluList,
    std::vector<bList>& hemiList,
    std::vector<bList>& lignList,
    int& nbr_poly,
    const int bond_selected,const int poly_selected,
    std::vector<double>& chem_entities,
    const int substrate,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_hemi);

//Updates the reaction tables after a digestion has been carried out
void update_reactiontables(
	std::vector<bList>& cellu,
    std::vector<bList>& hemi,
    std::vector<bList>& lign,
    std::vector<TList>& Table_cellu,
    std::vector<TList>& Table_hemi,
    params& par,
    std::vector<double>& chem_entities,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_hemi,
    double x,
    double y,
    int z,
    int substrate,
    int poly_selected,
    int bond_selected);


//Updates the reaction tables regarding inhibition through glucose or cellobiose
void update_reactiontables_inhib(
	std::vector<TList>& Table_cellu,
    std::vector<bList>& celluList,
    std::vector<TList>& Table_hemi,
    params& par,
    const int nbr_Glc_pdt,
    const int nbr_cellobiose,
    std::vector<double>& chem_entities);


//Counts the number of chains for different degrees of polymerization
void chain_length_distribution(
	std::vector<bList>& polyList,
    std::vector<DPList>& chain_distrib,
    const int length_fibril,
    const double real_time,
    const int t0,
    const int index);


//Checks whether all bonds with status = 1 have reactions associated to them 
bool check_reactions(
	std::vector<TList>& allTables,
    std::vector<bList>& polyList,
    params& par);

//Returns the bond type according to the type of left bond l and right bond r. Used for bond_type distribution on hemicellulose
int find_bond_type(
	int l,
    int r);


//Returns true if the bond exists within the polyLists, and fills the reference arguments with the polymer index and the material. If it returns false, the reference arguments should not be used further and are set to -1
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


//Looks through all reaction tables and looks at the relative fraction of reactions associated to each enzyme 
void check_reaction_table_distribution(
	const std::vector<TList>& Table_cellu,
    const std::vector<TList>& Table_hemi,
    TList& Table_lign,
    params& par);


//Returns the smaller of a and b
int findmin(
	int a,
    int b);


//Returns the larger of a and b
int findmax(
	int a,
    int b);


//Returns the absolute value of the vector (dx,dy,dz)
double dist(double dx, double dy, double dz);


//Adjusts propensities of all tables via "glueing enzymes"
void lignin_glue(
	params& par,
    std::vector<bList>& cellu,
    std::vector<TList>& Table_cellu,
    std::vector<TList>& Table_hemi,
    TList& Table_lign,
     std::vector<double>& chem_entities,
    int nbr_poly_lign);


//Returns a normally distributed random number between 0 and 1 with mean value at 0.5, according to the box-muller algorithm
double box_muller(
	double sigma,
    double mu);


//Returns the index of the polymer of coordinates x,y, which also contains a bond at position z
int find_poly(
	const std::vector<bList>& Table,
    const int x,
    const int y,
    const int z);


//Returns the index of the bond whose z coordinate is z_selected
int find_bond(
	const bList& poly,
    int z_selected);


//Returns the largest number of uninterrupted "true" elements of free_positions 
int largest_gap(
	bool *free_positions,
    int size,
    int& gap_start);


//Simple hole cutter function. Cuts 1D holes
void hole_cutter(std::vector<bList>& polyList, std::vector<TList>& Table_hemi, params& par, int& difference_monomers, int& nbr_poly, int material);

/*
//Cuts holes into outer wall. The holes are oval shaped. Currently disabled
void hole_cutter(
	std::vector<bList>& celluList,
    std::vector<bList>& hemiList,
    std::vector<bList>& lignList,
    std::vector<TList>& Table_hemi,
    params& par,
    int& difference_hemi,
    int& difference_lign,
    int& nbr_poly_hemi,
    int& nbr_poly_lign,
    bool& error_bool,
    int material);
*/


//Returns true if a bond exists at the specified x-, y- and z coordinates
bool bond_exists(
	const std::vector<bList>& celluList,
    const std::vector<bList>& hemiList,
    const std::vector<bList>& lignList,
    int x_pos,
    int y_pos,
    int z_pos);


//Returns true if there is no shell polymer blocking the one at position (x,y)
//bool is_outer_poly(const std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_cellu, const std::unordered_map<std::tuple<double,double,int>, neighborList>& bond_neighbors_hemi, const std::vector<bList>& hemiList, const std::vector<bList>& lignList, const std::vector<bList>& celluList, const int x, const int y, const int z, const int substrate, const params& par);


//Returns true if x_pos and y_pos are found as a valid combination inside celluList, hemiList or lignList
bool is_valid_pos(
	const std::vector<bList>& celluList,
    const std::vector<bList>& hemiList,
    const std::vector<bList>& lignList,int x_pos,
    int y_pos);


//Returns true if there is no outer bond blocking the one at position (x,y,z)
//bool is_outer_bond(const std::vector<bList>& celluList, const std::vector<bList>& hemiList, const std::vector<bList>& lignList, const int x, const int y, const int z, const int substrate, const params& par);


//Sets min_x,max_x,min_y,max_y for hemi and lign
void min_max_cellu_hemi_lign(
	const std::vector<bList>& celluList,
    params& par);//Set coordinates of boundaries for hemi and lign


//Returns an integer which is unique for every set (a,b,c)
int cantor_pair_three(
	int a,
    int b,
    int c);


//Returns an integer which is unique for every set (a,b)
int cantor_pair_two(
	int a,
    int b);


//Removes a bond from the neighbor vector
void remove_neighbor_from_vector(
	const int x_pos,
    const int y_pos,
    const int z_pos,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_cellu,
    std::unordered_map<std::tuple<double,double,int>,
    neighborList>& bond_neighbors_hemi);


//Returns true, if a CBH enzyme is attached to the bond bond at the poly poly
bool CBH_enzyme_attached(
	params& par,
    int poly,
    int bond);


//Updates CBH attachment reactions if there is a change in the number of free CBH enzymes
void update_attachment_reactions(
	params& par,
    std::vector<bList>& cellu,
    std::vector<TList>& Table_cellu,
    std::vector<double>& chem_entities);


//Sets the covering parameters of each lignin polymer. Only called once in the beginning
void distribute_lignin_covering(
	params& par,
	std::vector<bList>& lign);


//Checks whether action_mu1 = 6 reactions need to be blocked from occuring in this step
void check_for_blocked_ends(std::vector<bList>& cellu, std::vector<TList>& Table_cellu, params& par);

#endif
