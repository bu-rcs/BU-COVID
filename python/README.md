# Python bu_covid library

The bu_covid library is a set of functions that process the input GraphML and CSV files in the data structures required for Covasim simulations. Functions have extensive comments and docstrings that document how to call them. 

### __init__.py
Aside from the usual work of importing multiple files for a Python library, __init__.py applies changes to the Covasim PeopleMeta class to attach 5 extra attributes.  These are used to handle unique person ID's, undergraduate status, type of resident (off- or on-campus), and so on.  

### parallel.py
Parallel evaulations of simulations are implemented here.  Parallel computations are implemented using the Python multiprocessing library. As each computation is completed it is pickled and written to a temporary directory.  After the simulations are complete the pickled simulations are loaded back into memory. 

Two functions are exported:
* get_n_cores - gets the number of cores that can be used for parallel simulations.
* parallel_run_sims - the function that is called externally to run simulations in parallel.

### interventions.py
This implements numerous Covasim simulations including periodic testing, specifically testing undergraduates or residents of larger dorms, contact tracing, and pulsed infection events.  Note that not all of these are demonstrated in the example simulations.

The exported functions are:
* test_post_quar - Do an extra test of someone when they are exiting the quarantine period.
* gen_periodic_testing_interventions - Generate a set of interventions that test the population at regular intervals.
* gen_large_dorm_testing_interventions - Generate residential housing testing that allows for a higher frequency of testing for large dormitories.
* test_household_when_pos - Creates an intervention that tests and entire household when someone in the household has tested positive.
* gen_undergrad_testing_interventions - Create undergraduate-only periodic testing.
* contact_tracing_sens_spec - Do contact tracing for classroom and floor layers with specified tracing sensitivity and specificity.
* pulsed_infections_network - Add a "pulse" infection event on a day or set of days.  This is used to simulate say 10 people from a set of households who attend a party together on a weekend and get infected.  

### housing_and_classrooms.py
This provides a set of functions that are used in contact tracing and the setup of housing and classroom networks.

* gen_daily_interventions - Used to turn classroom networks on and off to match the class schedule. Typically used with the Covasim clip_edges intervention.
* get_contact_list - Return a list of contact lists from a graph.
* get_shared_class_contact_lists - Return a list of contact lists from a graph where classes are shared between students.
* get_all_class_contacts_dicts - Return a complete set of dictionaries that define Covasim networks for all of the classes.
* get_housing_contacts - Generate the contact lists for a housing network.
* get_all_housing_contacts_dict - Generate contact lists for all housing networks.

### misc.py
A grab bag of various useful functions related to file I/O. The exported functions are:

* load_people_info - Reads the required pop_info.csv file.
* make_dt - Make a Python datetime object from a specific date string.
* get_daynum - Returns the number of the day of the week from a datetime object.
* gen_BU_pop2 - Generate the population dictionary that's fed into Covasim.  Why the 2 suffix?  This was the 2nd version of this function and we never changed the name.
* validate_graphs - Double checks that two graphs have identical vertices.  Used for double-checking earlier in our development process.
* update_graph - Update a networkx graph with info from the BU_pop.  Useful if using the Covasim Networkx graph output.
* write_graphml - Create and save the transmission trees for all sims in a set of completed simulations.
* update_sim_people - Add our extra parameters to the sim.people object. 
* sim_results_to_df - Converts a list of completed sims to a set of Pandas dataframes.

### snapshots.py
Covasim allows for the gathering of statistics on a regular basis using what it calls snapshots.  This is a collectio of custom ones we found to be useful.  The snapshots are defined as Python classes.  Rather than interact with them directly the get_BU_snapshots function returns a set of snapshot objects that are attached to simulations.

* snapshots_to_df - Take a list of completed simulations with our snapshot analyzers are convert them to Pandas dataframes.
* get_BU_snapshots - Return a list of snapshots to be used with the simulations.  The order here is specific and must match that in snapshots_to_df. 
The remaining exported functions can convert different snapshot classes to dataframes: infection_count_to_df, diag2iso_count_to_df, severe_count_to_df, critical_count_to_df, and dead_count_to_df.


