# Simulation results

This directory contains results of running five simulations:

* *no_intervention* - baseline simulation;
* *classrooms* - social distancing and classroom masks intervention;
* *household* - in addition to the *classroom* intervention, reducing the size of a household;
* *testing* - in addition to the previous interventions, implementing regular testing;
* *tracing* - in addition to the above intervention, tracing intervention is performed. 

Each subdirectory contains a png file from Covasim output and seven csv files. sim_results_{simulation_name}.csv files contain the data used by R scripts to generate images to summarize the results and to compare various interventions:


![Single Run Cummulative Results](./Cummulative.png)

![Single Run Daily Results](./Daily.png)

