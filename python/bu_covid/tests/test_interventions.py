# placeholder. 

# How to test that interventions worked...create a simulation object with a 
# function call.  Add an intervention.  Run the simulation for a day (or 
# however many is needed for a test). Look for results.  Caution...these
# interventions might involve statistical quantities!

from bu_covid import interventions as bv
import bu_covid as bu
import covasim as cv
import pytest

@pytest.fixture
def make_test_sim():
    import covasim as cv
    from igraph import Graph
    import sciris as sc
    import os 
    import bu_covid as bu
    
    # Location of classroom networks
    net_dir = '../../../Data/networks'
    # Location of housing networks
    hnet_dir = '../../../Data/networks'
    # Location of pop_info.csv
    pop_info_path = '../../../Data/input/pop_info.csv'

    # Name of this simulation for output files:
    sim_name = 'testing'

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
    #  Build housing networks
    # =============================================================================
    h_files = ['RoomNet.graphml',
               'HouseholdNet_10.graphml',
               'FloorNet.graphml']
    hgraphs=[Graph.Read_GraphML(os.path.join(hnet_dir, hf)) for hf in h_files]
     
    hlayer_names = ['roommate','household','floor']
    household_contacts = bu.get_all_housing_contacts_dict(hgraphs, hlayer_names)

    #%%
    # =============================================================================
    #  Build BU population dictionary
    # =============================================================================
    # Use net_dir to load the all_age.txt and all_sex.txt files
    people_lut = bu.load_people_info(pop_info_path)
    BU_pop = bu.gen_BU_pop2(people_lut, class_contacts,household_contacts)

  
    # Run 1 day interventions by default
    start_day = '2020-01-01'
    end_day = '2020-01-02'

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
                pop_infected = 0, # Initial number of infected people
                beta = beta_val,      # Base transmission rate
                start_day = start_day,
                end_day = end_day,
                quar_factor = quar_factor,
                iso_factor = iso_factor,
                asymp_factor=0.5,
                n_imports=0) # n_imports - average spontaneous infections per day.
    
    return pars, BU_pop


def set_test_sim_date_range(pars,start_date,end_date):
    ''' Update the pars dict with a new date range.'''
    pars['start_day'] = start_date
    pars['end_day'] = end_date
    return pars


def test_pulsed_infections_network_const(make_test_sim):
    ''' Test a fixed number of infections'''
    pars, BU_pop = make_test_sim
    # Turn off infections from person-to-person
    pars['beta'] = 0 
    # Infect 10 people.
    num_days = (bu.make_dt(pars['end_day']) - bu.make_dt(pars['start_day'])).days
    days = (num_days + 1) * [0]
    # All infections on the first day.
    days[0] = 1
    NUM_INFECTIONS = 10
    interv = bv.pulsed_infections_network(BU_pop['contacts'],days,  layer_name='household', pulse_sampling='const', pulse_count=NUM_INFECTIONS)
    sim=cv.Sim(pars=pars, popfile=BU_pop, verbose = False, load_pop = True, 
               interventions = [interv], analyzers = [])
    # Update the sim with our extra categories
    bu.update_sim_people(sim,BU_pop)
    # run it
    sim.run()
    # Check that the correct number of infections on day 0 were created
    assert sim.results['new_infections'].values[0] == NUM_INFECTIONS
    # and none on day 1
    assert sim.results['new_infections'].values[1] == 0
  
def test_pulsed_infections_network_poisson(make_test_sim):
    ''' Test a poisson sampled number of infections'''
    pars, BU_pop = make_test_sim
    # Turn off infections from person-to-person
    pars['beta'] = 0 
    # Infect 10 people.
    num_days = (bu.make_dt(pars['end_day']) - bu.make_dt(pars['start_day'])).days
    days = (num_days + 1) * [0]
    # All infections on the first day.
    days[0] = 1
    NUM_INFECTIONS = 10  
    interv = bv.pulsed_infections_network(BU_pop['contacts'],days, layer_name='household', pulse_sampling='poisson', pulse_count=NUM_INFECTIONS)
    sim=cv.Sim(pars=pars, popfile=BU_pop, verbose = False, load_pop = True, 
               interventions = [interv], analyzers = [])
    # Update the sim with our extra categories
    bu.update_sim_people(sim,BU_pop)
    # Fix the random seed
    cv.utils.set_seed(987654321)
    # run it
    sim.run()
    # Check that the correct number of infections on day 0 were created
    assert sim.results['new_infections'].values[0] == 9
    # and none on day 1
    assert sim.results['new_infections'].values[1] == 0



def test_gen_periodic_testing_interventions_real(make_test_sim):
    # Generate tests for every 3.5 days and every 7 days.  Run the sim
    # for a week.  Check that the sum of sim.results['new_tests'] is
    # equal to the number of people in the sim.  beta is set to zero
    # to prevent any infections/symptomatic attestation going on and
    # triggering more tests.
    pars, BU_pop = make_test_sim
    # Turn off infections from person-to-person
    pars['beta'] = 0 
    # Set the dates. One week plus a day so we get a week of testing.
    pars['start_day'] = '2020-01-01'
    pars['end_day'] = '2020-01-08'
    num_days = (bu.make_dt(pars['end_day']) - bu.make_dt(pars['start_day'])).days
    intervs = bu.gen_periodic_testing_interventions_real(BU_pop, num_days, start_date = pars['start_day'])
    sim=cv.Sim(pars=pars, popfile=BU_pop, verbose = False, load_pop = True, 
               interventions = intervs, analyzers = [])
    # Update the sim with our extra categories
    bu.update_sim_people(sim,BU_pop)
    # Fix the random seed
    cv.utils.set_seed(987654321)
    # run it
    sim.run()    
    # Total number of tests
    assert int(sum(sim.results['new_tests'])) == 10818
    # each day
    assert list(sim.results['new_tests'])[0:-1] == [1683.0, 1697.0, 1747.0, 1119.0, 1128.0, 1697.0, 1747.0]
 
