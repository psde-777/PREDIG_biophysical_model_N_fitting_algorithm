name="low"
cp best_init_specs_${name}.txt Params/initial_configuration_parameters_${name}.txt
name="medium"
cp best_init_specs_${name}.txt Params/initial_configuration_parameters_${name}.txt
name="high"
cp best_init_specs_${name}.txt Params/initial_configuration_parameters_${name}.txt

cp best_kin_specs.txt Params/kinetic_parameters.txt

./code_4 low &
./code_4 medium &
./code_4 high &
