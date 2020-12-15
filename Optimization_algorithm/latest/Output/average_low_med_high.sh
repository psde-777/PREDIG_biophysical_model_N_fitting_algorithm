#$1: number of runs

python calc_mean_interpolate.py $1 enzyme_fraction/enzyme_fraction_low_ mean_enzyme_fraction_low.txt &
python calc_mean_interpolate.py $1 enzyme_fraction/enzyme_fraction_medium_ mean_enzyme_fraction_medium.txt &
python calc_mean_interpolate.py $1 enzyme_fraction/enzyme_fraction_high_ mean_enzyme_fraction_high.txt &

python calc_mean_interpolate.py $1 saccharification/saccharification_low_ mean_saccharification_low.txt &
python calc_mean_interpolate.py $1 saccharification/saccharification_medium_ mean_saccharification_medium.txt &
python calc_mean_interpolate.py $1 saccharification/saccharification_high_ mean_saccharification_high.txt &

python calc_mean_interpolate.py $1 Nbr_reactions/Nbr_reactions_low_ mean_Nbr_substrate_low.txt &
python calc_mean_interpolate.py $1 Nbr_reactions/Nbr_reactions_medium_ mean_Nbr_substrate_medium.txt &
python calc_mean_interpolate.py $1 Nbr_reactions/Nbr_reactions_high_ mean_Nbr_substrate_high.txt &