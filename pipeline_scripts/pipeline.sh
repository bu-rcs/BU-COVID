#!/bin/bash


# ===== Step 1 =====
# generate networks
cd ../R
Rscript pre-process_house_data.R
Rscript pre-process_classroom_data.R


# ===== Step 2 =====
# run simulations

export N_SIM_RUNS=100
cd ../python/bu_covid
python no_interventions.py
python social_distancing.py
python self_attestation.py
python testing.py

# ===== Step 3 =====
# generate output plots

cd ../../R
Rscript plot_no_intervention_results.R
Rscript plot_compare_interventions.R


