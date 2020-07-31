#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 12:40:29 2020

@author: bgregor
"""

''' To run on the SCC...

 source /rprojectnb/bucovid/Covasim/sourceme.sh
 
 python /rprojectnb/bucovid/brian/code/graphml/run_multisim_baseline.py


ADDS The intervention for people self-reporting themselves and getting tested. Line 218.

'''


#%%



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
net_dir = '/rprojectnb/bucovid/Networks'
# Location of housing networks
hnet_dir = '/rprojectnb/bucovid/Networks'
# Location of pop_info.csv
pop_info_path = '/rprojectnb/bucovid/Data/pop_info.csv'

# Name of this simulation for output files:
sim_name = 'baseline'

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
plot_dir = '/rprojectnb/bucovid/wenrui/Covasim/result'
if 'PLOT_DIR' in os.environ:  
    plot_dir = os.environ['PLOT_DIR']

# Set the number of parallel runs using the environment
n_runs = 1
if 'N_SIM_RUNS' in os.environ:
    n_runs=int(os.environ['N_SIM_RUNS'])

n_imports=1
if 'N_IMPORTS' in os.environ:
    n_imports = int(os.environ['N_IMPORTS'])
    
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
#if not bu.validate_graphs(graphs):
#    raise Exception('Graph vertices do not match!  Exiting.')
    
    
class_layer_names = ['sun','mon','tue','wed','thu','fri','sat']     
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

# The population size
pop_size=BU_pop['uid'].shape[0]

beta_val = 0.02


    
# Set up simulation parameters.  Get the pop size from the number of 
# vertices in any of the graphs (they're all the same)
pars = dict(pop_size = pop_size,
            #pop_infected = int(round(pop_size*0.0025)), # Initial number of infected people
            pop_infected = 0, 
            beta = beta_val,      # Base transmission rate
            start_day = start_day,
            end_day = end_day,
            asymp_factor=0.5,
            n_imports=n_imports) # n_imports - average spontaneous infections per day.




#%% =============================================================================
#     INTERVENTIONS    
# =============================================================================


base_interventions = [] 


# Set up varied beta per layer intervention
#  This is for platoons.
beta_layer = {}
for key in class_contacts:
    beta_layer[key] = 1
beta_layer['roommate'] = beta_roommate/beta_val,
beta_layer['household'] = beta_household/beta_val,
beta_layer['floor'] = beta_floor/beta_val,
beta_layer['building'] = beta_household/beta_val
for key in beta_layer:
    base_interventions.append(cv.change_beta(days=0, changes=beta_layer[key], layers=key))

# On class days turn on the class networks and reduce on non-class days.
class_interventions = []
for day_num, day in enumerate(class_layer_names):
    # In the layer_contact_lists find a matching day name
    for key in class_contacts:
        if key.find(day) >= 0:
            class_interventions.append(bu.gen_daily_interventions(start_day,end_day,[key],[day_num,],[1.0,0.0],cv.clip_edges))

base_interventions += class_interventions 


interventions = base_interventions.copy()

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
#snapshots_df.plot(x="days", y=["n_res_quar", "n_res_diag","n_nonres_quar","n_nonres_diag"])
#plt.show() 
#%% And the simulation results dataframe.  
sim_results_df = bu.sim_results_to_df(sims_complete)

#%%  Example of making an age histogram
#agehist = cv.age_histogram(sim=sims_complete[0])
#agehist.plot()



#%%


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
plots2 = sc.odict({ 'R_eff':['r_eff'] })
msim.plot(to_plot=plots2, show_args={'interventions':False},do_save=True,fig_path=os.path.join(plot_dir,'Sim_%s_Reff.png' % sim_name))



#%%
#print('Writing the graphml ZIP file')

#bu.write_graphml(sims_complete, plot_dir, sim_name, BU_pop)



    
