#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## ---------------------------
##
## Example of a baseline (no interventions) simulation.
##
## Authors: Brian Gregor, Wenrui Li 
##          Boston University
##
## Date Created: 2020-07-31
##
## Email: bgregor@bu.edu
##
## ---------------------------


import matplotlib
# With this line the plots will not show up on the screen
# but they will be saved. Comment it out if you want to see plots
# on the screen.  This MUST be called before the "import sciris" line.
#matplotlib.use('Agg')




import covasim as cv
from igraph import Graph
import sciris as sc
import os 
from matplotlib import pyplot as plt


# this ***MUST*** be imported after covasim.
import bu_covid as bu



# Location of classroom networks
net_dir = '../../Data/networks'  
# Location of housing networks
hnet_dir = '../../Data/networks'  
# Location of pop_info.csv
pop_info_path = '../../Data/input/pop_info.csv'

# Name of this simulation for output files:
sim_name = 'no_interventions'

# =============================================================================
#   baseline
#
#        Per-layer beta multiplicative factors
#        Classroom and housing networks in place.
#        classroom networks 'on' only on class days
#        Platooning turned off
#        Testing turned off
#        Tracing off
# ============================================================================= 

# Set a destination for plots and graphml files
plot_dir = '../../Data/results/no_interventions'
if 'PLOT_DIR' in os.environ:  
    plot_dir = os.environ['PLOT_DIR']

# Set the number of parallel runs using the environment
n_runs = 1
if 'N_SIM_RUNS' in os.environ:
    n_runs=int(os.environ['N_SIM_RUNS'])

n_imports=1
if 'N_IMPORTS' in os.environ:
    n_imports = int(os.environ['N_IMPORTS'])

# Base transmission rate
beta_val = 0.02
if 'BETA_val' in os.environ:
    beta_roommate = os.environ['BETA_val']
    
beta_roommate = 0.7
if 'BETA_roommate' in os.environ:
    beta_roommate = os.environ['BETA_roommate']

beta_household=0.08
if 'BETA_household' in os.environ:
    beta_household = os.environ['BETA_household']

beta_floor=0.03
if 'BETA_floor' in os.environ:
    beta_floor = os.environ['BETA_floor']


# Make sure the plot directory exists
os.makedirs(plot_dir, exist_ok = True) 

#%%
# =============================================================================
#  Build class networks
# =============================================================================
class_files = ['ClassNetworkU.graphml',
               'ClassNetworkM.graphml',
               'ClassNetworkT.graphml',
               'ClassNetworkW.graphml',
               'ClassNetworkR.graphml',
               'ClassNetworkF.graphml',
               'ClassNetworkS.graphml',
               ]
graphs = [Graph.Read_GraphML(os.path.join(net_dir, cf)) for cf in class_files]
   
class_layer_names = ['sun','mon','tue','wed','thu','fri','sat']     
# The classroom networks are modified into 3 groups per day 
#   students per class
#   an extra network for students who share 2 classes per day
#   an extra network for students who share 3+ classes per day.
# These are supersets, i.e. 1 classes contains all 1, 2 and 3+ students,
# 2 contains all 2 and 3+, 3+ contains 3+.
class_contacts = bu.get_all_class_contacts_dicts(graphs, class_layer_names, 3)

#%%
# =============================================================================
#  Build housing networks
# =============================================================================
h_files = ['RoomNet.graphml',
           'HouseholdNet_baseline.graphml',
           'FloorNet.graphml',
           'Building_net.graphml']
hgraphs=[Graph.Read_GraphML(os.path.join(hnet_dir, hf)) for hf in h_files]
 

#if not bu.validate_graphs(hgraphs):
#    raise Exception('Graph vertices do not match!  Exiting.')
    
hlayer_names = ['roommate','household','floor','building']
household_contacts = bu.get_all_housing_contacts_dict(hgraphs, hlayer_names)

#%%
# =============================================================================
#  Build BU population dictionary
# =============================================================================
# Use net_dir to load the all_age.txt and all_sex.txt files
people_lut = bu.load_people_info(pop_info_path)
BU_pop = bu.gen_BU_pop2(people_lut, class_contacts,household_contacts)


#%%

# =============================================================================
#  Build intervention
# =============================================================================
# The semester starts
start_day = '2020-09-02'
end_day = '2020-12-19'

# Total number of simulation days
num_days = (bu.make_dt(end_day) - bu.make_dt(start_day)).days

# Get the pop size from the number of vertices in any of the graphs (they're all the same)
pop_size=BU_pop['uid'].shape[0]

# Set up simulation parameters.  
pars = dict(pop_size = pop_size,
            pop_infected = 0,     # Start with nobody infected. 
            beta = beta_val,      # Base transmission rate
            start_day = start_day,
            end_day = end_day,
            asymp_factor=0.5,
            n_imports=n_imports) # n_imports - average spontaneous infections per day.


#%% =============================================================================
#     INTERVENTIONS    
# =============================================================================

# Our list of interventions.
interventions = [] 


# Set up varied beta per layer intervention
beta_layer = {}
# There is no change in beta for classrooms here, so the
# multiplier is 1.
for key in class_contacts:
    beta_layer[key] = 1
    
# Adjust for housing layers to reflect living together.
beta_layer['roommate'] = beta_roommate/beta_val,
beta_layer['household'] = beta_household/beta_val,
beta_layer['floor'] = beta_floor/beta_val,
beta_layer['building'] = beta_household/beta_val

# Add to the list of interventions so these will be applied.
for key in beta_layer:
    interventions.append(cv.change_beta(days=0, changes=beta_layer[key], layers=key))

# On class days turn on the class networks and off on non-class days.
# Add as covasim clip_edges interventions.
class_interventions = []
for day_num, day in enumerate(class_layer_names):
    # In the layer_contact_lists find a matching day name
    for key in class_contacts:
        if key.find(day) >= 0:
            # [1.0,0.0] - this is a multiplier: 1.0 when there's class, 0 where there isn't.
            # The second number could be altered if students interact on non-class days.
            class_interventions.append(bu.gen_daily_interventions(start_day,end_day,[key],[day_num,],[1.0,0.0],cv.clip_edges))
# and add the class schedule to our main intervention list.
interventions += class_interventions 


verbose = False # keep the simulations quiet
if n_runs == 1:
    verbose = True  #unles there's just 1.

#%% Create the list of daily snapshot objects
# These track numerous statistics. See get_BU_snapshots() for details.
analyzers = bu.get_BU_snapshots(num_days)

#%%
sim=cv.Sim(pars=pars, popfile=BU_pop, verbose = verbose, load_pop = True, 
           interventions = interventions, analyzers=analyzers)
# Update the sim with our extra categories
bu.update_sim_people(sim,BU_pop)
#%%

# Covasim uses the sciris library for parallelization, which provides a wrapper 
# on top of python.multiprocessing.  We were getting memory errors when running
# on 16+ cores simultaneously so we have a custom parallelization implemented
# here. It has been tested without issues on a 64-core system.
# get_n_cores() --> get the number of cores to run on.  Please read the docstring
# for this function to see how to use it.
sims_complete = bu.parallel_run_sims(sim, n_runs = n_runs, n_cores = bu.get_n_cores())

# Add the list of completed sims to a MultiSim object.
msim = cv.MultiSim(sims_complete)
# Calculate statistics.
msim.reduce()


print("Creating dataframes") 
#%% Get the snapshots dataframe
snapshots_df = bu.snapshots_to_df(sims_complete) 
 #%%
infection_df=bu.infection_count_to_df(sims_complete)

#%%
diag2iso_df=bu.diag2iso_count_to_df(sims_complete)

severe_df=bu.severe_count_to_df(sims_complete)

critical_df=bu.critical_count_to_df(sims_complete)

dead_df=bu.dead_count_to_df(sims_complete)

#%% And the simulation results dataframe.  
sim_results_df = bu.sim_results_to_df(sims_complete)

# =============================================================================
#  Output
# =============================================================================

# All the results are now in memory.  Write out the dataframes as CSV
# files.  Generate one plot using Covasim and save it.

print('Generating plots')

os.makedirs(plot_dir,exist_ok = True)

# Write out the snapshots_df
with  open(os.path.join(plot_dir,'snapshots_%s.csv' % sim_name), 'w') as f:
    snapshots_df.to_csv(f)
# And the sim_results_df
with  open(os.path.join(plot_dir,'sim_results_%s.csv' % sim_name), 'w') as f:
    sim_results_df.to_csv(f)

# Write out the infection df
with  open(os.path.join(plot_dir,'infection_counts_%s.csv' % sim_name), 'w') as f:
    infection_df.to_csv(f)

# Write out the diag2iso df
with  open(os.path.join(plot_dir,'diag2iso_counts_%s.csv' % sim_name), 'w') as f:
    diag2iso_df.to_csv(f)

with  open(os.path.join(plot_dir,'severe_counts_%s.csv' % sim_name), 'w') as f:
    severe_df.to_csv(f)

with  open(os.path.join(plot_dir,'critical_counts_%s.csv' % sim_name), 'w') as f:
    critical_df.to_csv(f)
    
with  open(os.path.join(plot_dir,'dead_counts_%s.csv' % sim_name), 'w') as f:
    dead_df.to_csv(f)   
    
plots1 = sc.odict({
                'Total counts': [
                    'cum_infections',
                    #'cum_diagnoses',
                    'cum_recoveries',
                    #'cum_quarantined'
                ],
                'Daily counts': [
                    'new_infections',
                   # 'new_diagnoses',
                    'new_recoveries',
                    'new_deaths',
                   # 'n_quarantined',
                #    'n_diagnosed',
                   # 'new_diagnoses'
                ],
                'Health outcomes': [
                    'cum_severe',
                    'cum_critical',
                    'cum_deaths',
                ],
        })

msim.plot(to_plot=plots1, show_args={'interventions':False},do_save=True,fig_path=os.path.join(plot_dir,'Sim_%s.png' % sim_name))


    
