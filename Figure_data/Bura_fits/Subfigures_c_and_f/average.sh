name="low"
python3 calc_mean_interpolate.py $1 Output/saccharification/saccharification_${name}_ mean_saccharification_${name}.txt Output/expe_data/expe_saccharification_${name}_ &
name="medium"
python3 calc_mean_interpolate.py $1 Output/saccharification/saccharification_${name}_ mean_saccharification_${name}.txt Output/expe_data/expe_saccharification_${name}_ &
name="high"
python3 calc_mean_interpolate.py $1 Output/saccharification/saccharification_${name}_ mean_saccharification_${name}.txt Output/expe_data/expe_saccharification_${name}_ &
