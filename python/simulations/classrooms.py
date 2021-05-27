#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## ---------------------------
##
## Example of adding social distancing in the classrooms to the simulation.
## This is done by splitting the classrooms into sub-groups ("platoons") 
## who attend their class once per week.
##
## Authors: Brian Gregor, Wenrui Li 
##          Boston University
##
## Date Created: 2020-07-31
##
## Email: bgregor@bu.edu
##
## ---------------------------


#%%



import matplotlib
# With this line the plots will not show up on the screen
# but they will be saved. Comment it out if you want to see plots
# on the screen.  This MUST be called before the "import sciris" line.
matplotlib.use('Agg')




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
sim_name = 'classrooms'

# =============================================================================
#   classrooms
#
#        Per-layer beta multiplicative factors
#        Classroom and housing networks in place.
#        classroom networks 'on' only on class days
#        Platooning turned on
#        Housing density unchanged
#        self-attestation for testing is off
#        Testing is off
#        Tracing is off
# ============================================================================= 

# Set a destination for plots and graphml files
plot_dir = '../../Data/results'

if 'PLOT_DIR' in os.environ:  
    plot_dir = os.environ['PLOT_DIR']  

# Set the number of parallel runs using the environment
n_runs = 1
if 'N_SIM_RUNS' in os.environ:
    n_runs=int(os.environ['N_SIM_RUNS'])

n_imports=1
if 'N_IMPORTS' in os.environ:
    n_imports = float(os.environ['N_IMPORTS'])
    
beta_roommate = 0.7
if 'BETA_roommate' in os.environ:
    beta_roommate = float(os.environ['BETA_roommate'])

beta_household=0.08
if 'BETA_household' in os.environ:
    beta_household = float(os.environ['BETA_household'])

beta_floor=0.03
if 'BETA_floor' in os.environ:
    beta_floor = float(os.environ['BETA_floor'])
# Make sure the plot directory exists

test_sensitivity = 0.9
if 'TEST_sensitivity' in os.environ:
    test_sensitivity = float(os.environ['TEST_sensitivity'])

plot_dir = os.path.join(plot_dir,sim_name)
os.makedirs(plot_dir, exist_ok = True) 

#%%
# =============================================================================
#  Build class networks
# =============================================================================
class_files = ['ClassNetPlatoonU.graphml',
               'ClassNetPlatoonM.graphml',
               'ClassNetPlatoonT.graphml',
               'ClassNetPlatoonW.graphml',
               'ClassNetPlatoonR.graphml',
               'ClassNetPlatoonF.graphml',
               'ClassNetPlatoonS.graphml',              
               ]
graphs = [Graph.Read_GraphML(os.path.join(net_dir, cf)) for cf in class_files]

class_layer_names = ['sun','mon','tue','wed','thu','fri','sat']     
class_contacts = bu.get_all_class_contacts_dicts(graphs, class_layer_names, 3)

#%%
# =============================================================================
#  Build housing networks.  
# =============================================================================
h_files = ['RoomNet.graphml',
           'HouseholdNet_baseline.graphml',
           'FloorNet.graphml',
           'Building_net.graphml']
hgraphs=[Graph.Read_GraphML(os.path.join(hnet_dir, hf)) for hf in h_files]
    
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
 
# The semester starts
start_day = '2020-09-02'
end_day = '2020-12-19'

# Total number of simulation days
num_days = (bu.make_dt(end_day) - bu.make_dt(start_day)).days

# The population size
pop_size=BU_pop['uid'].shape[0]

beta_val = 0.02


# Set the quarantine and isolation beta multipliers
quar_factor = {}
iso_factor = {} 
for key in {**household_contacts, **class_contacts}:
    quar_factor[key] = 0.0
    iso_factor[key] = 0.0
    

# Set up simulation parameters.  Get the pop size from the number of 
# vertices in any of the graphs (they're all the same)
pars = dict(pop_size = pop_size,
            #pop_infected = int(round(pop_size*0.0025)), # Initial number of infected people
            pop_infected = 0, 
            beta = beta_val,      # Base transmission rate
            start_day = start_day,
            end_day = end_day,
            quar_factor = quar_factor,
            iso_factor = iso_factor,
            asymp_factor=0.5,
            n_imports=n_imports) # n_imports - average spontaneous infections per day.





# =============================================================================
#     Set up INTERVENTIONS
# =============================================================================

interventions = []

# Set beta for different layers
beta_layer = {}
# The classrooms are not as dense and masks are being worn.
# Reduce the classroom beta by half.  beta_layer is a multiplier.
for key in class_contacts:
    beta_layer[key] = 0.5
beta_layer['roommate'] = beta_roommate/beta_val,
beta_layer['household'] = beta_household/beta_val,
beta_layer['floor'] = beta_floor/beta_val
for key in beta_layer:
    interventions.append(cv.change_beta(days=0, changes=beta_layer[key], layers=key))

# Schedule classes by day of the week.    
class_interventions = []
for day_num, day in enumerate(class_layer_names):
    # In the layer_contact_lists find a matching day name
    for key in class_contacts:
        if key.find(day) >= 0:
            class_interventions.append(bu.gen_daily_interventions(start_day,end_day,[key],[day_num,],[1.0,0.0],cv.clip_edges))

interventions += class_interventions





#%%
# =============================================================================
#  Run the simulations
# =============================================================================


verbose = False # keep the simulations quiet
if n_runs == 1:
    verbose = True  #unles there's just 1.

#%% Create the list of daily snapshot objects
analyzers = bu.get_BU_snapshots(num_days)

#%%
sim=cv.Sim(pars=pars, popfile=BU_pop, verbose = verbose, load_pop = True, 
           interventions = interventions, analyzers=analyzers)
# Update the sim with our extra categories
bu.update_sim_people(sim,BU_pop)
#%%

sims_complete = bu.parallel_run_sims(sim, n_runs = n_runs, n_cores = bu.get_n_cores())

msim = cv.MultiSim(sims_complete)
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

print('Generating plots')

os.makedirs(plot_dir,exist_ok = True)

# Write out the snapshots_df
with  open(os.path.join(plot_dir,'snapshots_%s.csv' % sim_name), 'w') as f:
    snapshots_df.to_csv(f)
# And the sim_results_df
with  open(os.path.join(plot_dir,'sim_results_%s.csv' % sim_name), 'w') as f:
    sim_results_df.to_csv(f)

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
                    'cum_diagnoses',
                    'cum_recoveries',
                    #'cum_quarantined'
                ],
                'Daily counts': [
                    'new_infections',
                    'new_diagnoses',
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
 
