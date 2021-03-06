# Example simulations

These are examples of setting up simulations in Python.  

* *no_interventions.py* - Run the simulations without testing, tracing, or social distancing.
* *classrooms.py* - Classrooms are split into "platoons" that meet physically in one class per week, masks are worn, and classrooms are less dense.
* *households.py* - Housing density is reduced in addition to the interventions in classrooms.py. 
* *testing.py* - Add weekly testing to the population in addition to symptomatic people reporting themselves to get tested.  
* *tracing.py* - Add contact tracing interventions are added for roommates, housemates, floormates, and classmates.

The simulations read the following environment variables.  The TEST variable is used in testing.py and tracing.py.  The TRACE variables are only
read by the tracing.py simulation:

* PLOT_DIR - The output directory for simulation files. The simulation's name (stored as the variable "sim_name") is appended to this path. The default is "../../Data/results".
* N_SIM_RUNS - The number of copies of the simulation to run. The default is 1.
* NSLOTS - The number of CPU cores to use for parallel computations. The default is 4 or the number of physical cores, whichever is lower.
* N_IMPORTS - The Poisson average number of daily exogeneous infections.  The default is 1.
* BETA_val - Base transmission rate, i.e. beta.  The default is 0.02.
* BETA_roommate - Beta adjustment for roommates.  Effective beta is BETA_roommate/BETA_val.  The default is 0.7.
* BETA_household - Beta adjustment for households.  Effective beta is BETA_household/BETA_val.  The default is 0.08.
* BETA_floor - Beta adjustment for roommates.  Effective beta is BETA_floor/BETA_val.  The default is 0.03.
* TEST_sensitivity - Covid-19 test sensitivity.  Default is 0.9.
* TRACE_sensitivity - Floormate and classmate tracing sensitivity. Default is 0.3.
* TRACE_specificity - Floormate and classmate tracing specificity.  The default is 0.9.
 
# Running the simulations

After bu_covid is installed, run the simulations by setting any required environment variables and then by calling the Python interpreter.

Here is an example using a Linux bash shell which overrides default values for 4 parameters:
```
export PLOT_DIR=~/covid_sims
export N_SIM_RUNS=100
export NSLOTS=8
export N_IMPORTS=3
python tracing.py
# output will be written to:  ~/covid_sims/tracing
```
