# PREDIG: Lignocellulose Enzymatic Digestion Model (forked from Plant Cell Wall Saccharification Model (PCWSM))

***IMPORTANT: Please read the NOTICE and LICENSE files in this repository for important copyright and licensing information regarding the contents of this project.***

Within this repository, you can find the code of our model for the saccharification of a lignocellulosic plant cell wall sub-unit, a microfibril composed of cellulose surrounded by hemicellulose and lignin.

## Content

The repository is split into two folders, each of which are summarised below.

### Folder "Simulation_code"

This folder contains the newest release version of our biophysical model's code, written in C++. It also contains scripts for averaging the output data & generating animations of the saccharification process, written in Python (version used 3.9.2). A further README for running all of these codes is provided within the folder.

### Folder "Optimization_algorithm"

This folder contains our algorithm for fitting the parameters of the biophysical model, as well as a README for running it.
