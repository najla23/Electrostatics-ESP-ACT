# SI-Accurate-Electrostatics-for-Physics-Based-Force-Fields
Supporting information for reproducing results in the scientific article
[_Accurate Electrostatics for Physics-Based Force Fields through Machine Learning_].


## Directory Layout

- `Analytical_Fitting/` - Analytical fitting for ESP
- `Antechamber/` - Force field files 
- `Charge_Models/` - Reproduction of Table 2
- `Database/` - Database for training 
- `ESP_Alkali_Halides/` - ESP fitting and electrostatic calculations for different charge models, e.g. a (positive) point charge with either one Gaussian or 1S Slater distributed charge, or a 
   point charge with two Gaussian charges or a point charge with a 1S and a 2S Slater charge.
- `Fig_Script/` - Reproduction of the figure 2 presented in the article
- `AlexandriaFF/` - Trained force field files on SAPT using Alexandria
- `ForceFields/` - Available force field files such as TIP3P, TIP4P, GAFF,... 
- `SAPT_Alkali_Halides/` - SAPT calculations  
- `Selection/` - Data sets and a script to generate random test and train compounds 
- `Tab_Script/` - Reproduction of tables presented in the article

## Requirements

To run the scripts, you'll need Python installed on your computer, if you want to run on your own computer,
install python using e.g. [Miniconda](https://conda.io/miniconda.html) or [Anaconda](https://docs.conda.io))
