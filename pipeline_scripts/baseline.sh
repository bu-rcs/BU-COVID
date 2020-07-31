#!/bin/bash


# ===== Step 1 =====
# generate networks
cd ../R
Rscript pre-process_house_data.R
Rscript pre-process_classroom_data.R


# ===== Step 2 =====
# install COVASIM and run simulations
cd ../python/bu_covid
python setup.py install --user


