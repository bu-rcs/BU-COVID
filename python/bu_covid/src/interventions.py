#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 15:21:02 2020

@author: bgregor
"""

import covasim as cv
import numpy as np
import sciris as sc
import covasim.utils as cvu
import datetime
from .misc import make_dt, get_daynum

__all__=['test_post_quar', 'gen_periodic_testing_interventions','gen_large_dorm_testing_interventions',
         'test_household_when_pos','pulsed_infections_network','gen_undegrad_testing_interventions','contact_tracing_sens_spec']

class test_post_quar(cv.Intervention):
    '''
    Test people daily for a fixed number of days AFTER they are finished with quarantine.  

    Args:
        n_days_post_quar (int)   : number of days to test people post-quarantine
        sensitivity (float) : test sensitivity
        loss_prob (float): Probability of loss to follow-up
        test_delay  (int)   : days for test result to be known
        start_day   (int)   : day the intervention starts
        end_day     (int)   : day the intervention ends
        kwargs      (dict)  : passed to Intervention                                                                                   ( )

    **Examples**::

        interv = cv.test_post_quar(n_days_post_quar=5)
    '''

    def __init__(self, n_days_post_quar = 5, subtarget = None, sensitivity=1.0, loss_prob=0, test_delay=0,
                 start_day=0, end_day=None,**kwargs):
        super().__init__(**kwargs) # Initialize the Intervention object
        self._store_args() # Store the input arguments so the intervention can be recreated
        self.n_days_post_quar = n_days_post_quar # Should be a list of length matching time
        self.sensitivity = sensitivity
        self.loss_prob  = loss_prob
        self.test_delay  = test_delay
        self.start_day   = start_day
        self.end_day     = end_day
        return


    def initialize(self, sim):
        ''' Fix the dates and number of tests '''

        # Handle days
        self.start_day   = sim.day(self.start_day)
        self.end_day     = sim.day(self.end_day)
        self.days        = [self.start_day, self.end_day]

        self.initialized = True
        
        return


    def apply(self, sim):

        t = sim.t
        if t < self.start_day:
            return
        elif self.end_day is not None and t > self.end_day:
            return

        test_probs = np.zeros(sim.n) # Begin by assigning zero probability to everyone
        
        # Get the indices of everyone who left quarantine less than n_days_post_quar ago.
        ended_quar_inds = np.where(np.isfinite(sim.people.date_end_quarantine))
        inds = []
        # this could be vectorized
        for eqi in ended_quar_inds[0]:
            tmp = int(sim.t - sim.people.date_end_quarantine[eqi])
            # Between 1 and n_days_post_quar days after quarantining
            if tmp >= 1 and tmp <= self.n_days_post_quar:
                inds.append( eqi )
        #if len(inds) > 0:
        #    import pdb ; pdb.set_trace()
        # Let's test those people.
        test_probs[inds] = 1.0

        # Don't re-diagnose people, i.e. sometimes people leave quarantine because 
        # they went into isolation.
        diag_inds  = cvu.true(sim.people.diagnosed)
        test_probs[diag_inds] = 0.

        # Now choose who gets tested and test them
        test_inds = np.where(test_probs > 0)
        sim.people.test(test_inds, self.sensitivity, loss_prob=self.loss_prob, test_delay=self.test_delay)

        return



#%%
# =============================================================================
#   Create the interventions to test everyone every 7 days OR
#      big-dorm residents every 3 days and everyone else every 7 days.
# =============================================================================

def build_subtarget_list(i, idx_list, to_ndarray=True):
    ''' From the list of index arrays in idx_list, generate
        a list of indices not in list i.  Sort and
        return it as an ndarray'''
    subtarget = []
    for j,idx in enumerate(idx_list):
        if j != i:
            subtarget += list(idx)
    subtarget.sort()
    if to_ndarray:
        return np.asarray(subtarget, dtype=np.int32)
    return subtarget



def gen_periodic_testing_interventions(BU_pop, num_days,  test_period=7, **kwargs):
    ''' Test everyone in the BU population every test_period days.
        
        BU_pop:  BU population dictionary
        num_days: length of simulation
        test_period:  Number of days between tests
        ...and also all the regular arguments to cv.test_num
        
        Returns: A list of cv.test_num interventions
    '''
    # Split up the population into test_period pieces.
    uids = BU_pop['uid'].copy()
    # shuffle it randomly
    np.random.shuffle(uids)
    # Now split into test_period random parts.
    rand_idx = np.array_split(uids,test_period)
    # ran_idx is now a list of test_period ndarrays of indices
    intervs = []  # The interventions to return.

    for i,idx in enumerate(rand_idx):
        # Set up a multiplier for the tests per day argument.
        test_days = ((i + np.arange(num_days)) % test_period < 1).astype(np.int32)
        num_tests = idx.shape[0]
        daily_tests = num_tests * test_days
        # Get the subtarger to exclude
        exclude_idx = build_subtarget_list(i, rand_idx)
        # And the matching zeros to zero-out their probability of being tested
        exclude_vals = np.zeros(exclude_idx.shape, dtype=np.float32)
        intervs.append(cv.test_num(daily_tests=daily_tests,
                                   subtarget={'inds':exclude_idx, 'vals':exclude_vals},
                                   **kwargs))
    return intervs
 
def gen_large_dorm_testing_interventions(BU_pop, num_days, test_period_large = 3,  test_period=7, **kwargs):
    ''' Test everyone in the large dorm BU population every test_period_large days.
        Everyone else is tested every test_period days.
        
        BU_pop:  BU population dictionary
        num_days: length of simulation
        test_period_large: 
        test_period:  Number of days between tests
        ...and also all the regular arguments to cv.test_num
        
        Returns: A list of cv.test_num interventions
    '''
    # Split the uids into 2 populations - in larger dorms and out of them
    uids_large = BU_pop['uid'][np.where(BU_pop['campResident'] > 1)].copy()
    uids_other = BU_pop['uid'][np.where(BU_pop['campResident'] < 2)].copy()
    # shuffle 'em up
    np.random.shuffle(uids_large)
    np.random.shuffle(uids_other)

    # Split into two sets of indices
    rand_large = np.array_split(uids_large,test_period_large)
    rand_other = np.array_split(uids_other,test_period)
    
    intervs = []
    
    # we'll need all of the rand_other indices for exclusions and vice versa
    rand_others_all = []
    rand_large_all = []
    for ro in rand_other:
        rand_others_all += list(ro)
    for rl in rand_large:
        rand_large_all += list(rl)
        
    # Do the large dorm ones first.
    for i,idx in enumerate(rand_large):
        # Set up a multiplier for the tests per day argument.
        test_days = ((i + np.arange(num_days)) % test_period_large < 1).astype(np.int32)
        num_tests = idx.shape[0]
        daily_tests = num_tests * test_days
        # Get the subtarger to exclude
        exclude_idx = build_subtarget_list(i, rand_large, to_ndarray = False)
        exclude_idx += rand_others_all
        exclude_idx.sort()
        exclude_idx = np.asarray(exclude_idx, dtype=np.int32)
        # And the matching zeros to zero-out their probability of being tested
        exclude_vals = np.zeros(exclude_idx.shape, dtype=np.float32)
        intervs.append(cv.test_num(daily_tests=daily_tests,
                                   subtarget={'inds':exclude_idx, 'vals':exclude_vals},
                                   **kwargs))

    for i,idx in enumerate(rand_other):
        # Set up a multiplier for the tests per day argument.
        test_days = ((i + np.arange(num_days)) % test_period < 1).astype(np.int32)
        num_tests = idx.shape[0]
        daily_tests = num_tests * test_days
        # Get the subtarger to exclude
        exclude_idx = build_subtarget_list(i, rand_other, to_ndarray = False)
        exclude_idx += rand_large_all
        exclude_idx.sort()
        exclude_idx = np.asarray(exclude_idx, dtype=np.int32)
        # And the matching zeros to zero-out their probability of being tested
        exclude_vals = np.zeros(exclude_idx.shape, dtype=np.float32)
        intervs.append(cv.test_num(daily_tests=daily_tests,
                                   subtarget={'inds':exclude_idx, 'vals':exclude_vals},
                                   **kwargs))
        
    # and that's all, folks
    return intervs
    

#%%   7/14/2020   Intervention designed to: 
#     1.	
class test_household_when_pos(cv.Intervention):
    '''
    If a person tests positive, quarantine the roommate + test the household.

    Args:
        contacts (dict)     : contacts dictionary from BU_pop
        sensitivity (float) : test sensitivity
        loss_prob (float): Probability of loss to follow-up
        test_delay  (int)   : days for test result to be known
        start_day   (int)   : day the intervention starts
        end_day     (int)   : day the intervention ends
        household_test_delay (int) : delay in days to get the household tested after the 
                                     household member has tested positive. Default of 1.
        household_layer_name (str) : name of the household layer.  Default is 'household'
        roommate_layer_name (str) : name of the roommate layer.  Default is 'roommate' 
        kwargs      (dict)  : passed to Intervention                                                                                   ( )

    **Examples**::

        interv = cv.test_household_when_pos(BU_pop['contacts'], household_test_delay=2)
    '''

    def __init__(self, contacts, household_test_delay = 1, subtarget = None, sensitivity=1.0, loss_prob=0, test_delay=0,
                 start_day=0, household_layer_name = 'household', roommate_layer_name = 'roommate', end_day=None,**kwargs):
        super().__init__(**kwargs) # Initialize the Intervention object
        self._store_args() # Store the input arguments so the intervention can be recreated
        self.household_test_delay = household_test_delay  
        self.rmln = roommate_layer_name
        self.hhln = household_layer_name
        self.sensitivity = sensitivity
        self.loss_prob  = loss_prob
        self.test_delay  = test_delay
        self.start_day   = start_day
        self.end_day     = end_day
        self.roommates, self.households = self.get_rm_hh(contacts)
        return


    def initialize(self, sim):
        ''' Fix the dates and number of tests '''

        # Handle days
        self.start_day   = sim.day(self.start_day)
        self.end_day     = sim.day(self.end_day)
        self.days        = [self.start_day, self.end_day]

        self.initialized = True
        return


    def apply(self, sim):
        # Boilerplate - apply this intervention in the specified interval
        t = sim.t
        if t < self.start_day:
            return
        elif self.end_day is not None and t > self.end_day:
            return

        test_probs = np.zeros(sim.n) # Begin by assigning zero probability to everyone
        
        # Loop through the household networks.  For each person, check to see if their
        # date_diagnosed is household_test_delay days previous.  If so, assign the 
        # entire household to the test list and test them.
        diag_inds = np.where(sim.t - sim.people.date_diagnosed == self.household_test_delay)[0]
        # find anyone?
        if diag_inds.shape[0] > 0:
            # Fetch their household network.  It's in a list so just loop.
            test_inds = set() # Make it a unique set of inds
            for diag in diag_inds:
                for hh in self.households[diag]:
                    test_inds.add(hh)
                # Remove the roommate from the test_inds as they'll get tested
                # automatically when they get quarantined.
                for rm in self.roommates[diag]:
                    test_inds.remove(rm)
            
            if len(test_inds) > 0:
                # Filter the test_inds set.  Remove anyone who is not susceptible
                tmp = test_inds.copy()
                for ti in tmp:
                    if not sim.people.susceptible[ti]:
                        test_inds.remove(ti)
                           
                # Now test those people.
                if len(test_inds) > 0:
                        sim.people.test(np.array(list(test_inds),dtype=np.int32), self.sensitivity, loss_prob=self.loss_prob, test_delay=self.test_delay)

        return
    
    def get_rm_hh(self, contacts):
        ''' contacts is a list, per-person, containing their contact networks.
            Extract out the roommate and household contact lists into a tuple.'''
        rm_lst = []
        hh_lst = []
        for person in contacts:
            rm_lst.append(person[self.rmln].copy())        
            hh_lst.append(person[self.hhln].copy())
        return rm_lst, hh_lst
  
    
  
#%%   7/15/2020   Pulsed infections
#  	
class pulsed_infections_network(cv.Intervention):
    '''
    Infects a set of people together in a contact network.  
    
    To use, set the number of people to infect in a pulse and the number of 
    networks in a layer to choose them from.  The days is a list or array of 0/1 or 
    Booleans.  When 0 the pulse is not applied, when 1 it is.  The network(s) is selected
    at random.  The people in the network(s) willbe filtered to exclude those in quarantine,
    hospital, etc. Then the remaining number is infected. 

    Args:
        contacts (dict)     : contacts dictionary from BU_pop
        days (list) : list or array of 0 or 1 (or Booleans), one per day.  When 1 the pulse is applied.
        layer_name (str) : The layer name to use for infections.  Default is 'household'
        pulse_count (int) : The number of people to infect.  Default is 10.
        pulse_sampling (str) : 'const' or 'poisson'.  Default is poisson.
        kwargs      (dict)  : passed to Intervention                                                                                   ( )

    **Examples**::

        # days_arr is 1 for every Friday. Infect 10 people.
        interv = bu.pulsed_infections(BU_pop['contacts'],days_arr, layer_name='household', pulse_count=10, network_count=1)
    '''
    def __init__(self, contacts, days, layer_name = 'household', pulse_count = 10, pulse_sampling = 'poisson', **kwargs):
        super().__init__(**kwargs) # Initialize the Intervention object
        self._store_args() # Store the input arguments so the intervention can be recreated
        if pulse_count >= 1:
            self.pulse_count = pulse_count
        else:
            raise Exception('pulse_count must be >= 1')
        self.days = np.array(days, dtype = np.bool)
        self.networks = self.get_network(contacts, layer_name)
        self.pulse_sampling = pulse_sampling.strip().lower()
        self.start_day = None
        self.end_day = None
        return


    def initialize(self, sim):
        ''' Fix the dates and number of tests '''
        # Handle days
        self.start_day   = sim.day(self.start_day)
        self.end_day     = sim.day(self.end_day)
        self.initialized = True
        return

    def get_network(self, contacts, layer_name):
        ''' Retrieve the selected network from the contacts dictionary'''
        nx = set()
        # For each person, add them to their own household network of contacts,
        # sort the list, and add it to the set. 
        for idnum, person in enumerate(contacts):
            if len(person['household']) > 0:
                tmp = person['household'].copy()
                tmp.append(idnum)
                tmp.sort()
                nx.add(tuple(tmp))
        # Convert to a list and return
        return list(nx)

    def apply(self, sim):
        # If today is not a pulse day then return.
        if not self.days[sim.t]:
            return
        # list of those who will be infected
        targets = []
        # limit the number of loop iterations for safety
        max_loops = 100
        iters = 0
        # Index list for the networks
        net_inds = list(range(len(self.networks)))
        
        # Use either constant infection count or a Poisson average
        pulse_count = self.pulse_count
        if self.pulse_sampling == 'poisson':
            pulse_count = np.random.poisson(self.pulse_count)  
        
        while len(targets) < pulse_count and iters < max_loops:
            # Pick a household. Find everyone in it who is susceptible.
            net = np.random.choice(net_inds)
            # Remove this from the possible choices
            del net_inds[net_inds.index(net)]
            # e-z: loop and add to targets if susceptible
            for n in self.networks[net]:
                if sim.people.susceptible[n]:
                    targets.append(n)
            iters += 1
        # Trim any extras
        if len(targets) > pulse_count:
            del targets[pulse_count:]
        # Now infect those targets
        sim.people.infect(inds=np.array(targets, dtype = np.int32))
        return
    
    @staticmethod
    def get_days_arr(start_day, end_day, pulse_day):
        ''' Returns a numpy array for the days to apply a pulse.
        start_day, end_day are date strings.
        pulse_day is th enumber of the day of the week to apply
        the pulse.  Sunday=0.'''
        num_days = (make_dt(end_day) - make_dt(start_day)).days + 1

        days_arr = np.zeros(num_days, dtype=np.int32)
        # datetime for the start of the sim
        sdt = make_dt(start_day)
        # And a 1 day time delta
        td = datetime.timedelta(days=1)
        # What day of the week (number) do you want pulses on?
        # Friday = 5 (sunday is 0)
        for i in range(num_days):
            today = sdt + i * td
            if get_daynum(today) == pulse_day:
                days_arr[i] = 1
        return days_arr
        

# Test undergrad every test_period_large days. Everyone else is tested every test_period days.     
def gen_undegrad_testing_interventions(BU_pop, num_days, test_period_large = 3,  test_period=7, **kwargs):
    ''' Test undergrad every test_period_large days.
        Everyone else is tested every test_period days.
        
        BU_pop:  BU population dictionary
        num_days: length of simulation
        test_period_large: 
        test_period:  Number of days between tests
        ...and also all the regular arguments to cv.test_num
        
        Returns: A list of cv.test_num interventions
    '''
    # Split the uids into 2 populations - in larger dorms and out of them
    uids_large = BU_pop['uid'][np.where(BU_pop['undegrad'] > 0)].copy()
    uids_other = BU_pop['uid'][np.where(BU_pop['undegrad'] < 1)].copy()
    # shuffle 'em up
    np.random.shuffle(uids_large)
    np.random.shuffle(uids_other)

    # Split into two sets of indices
    rand_large = np.array_split(uids_large,test_period_large)
    rand_other = np.array_split(uids_other,test_period)
    
    intervs = []
    
    # we'll need all of the rand_other indices for exclusions and vice versa
    rand_others_all = []
    rand_large_all = []
    for ro in rand_other:
        rand_others_all += list(ro)
    for rl in rand_large:
        rand_large_all += list(rl)
        
    # Do the large dorm ones first.
    for i,idx in enumerate(rand_large):
        # Set up a multiplier for the tests per day argument.
        test_days = ((i + np.arange(num_days)) % test_period_large < 1).astype(np.int32)
        num_tests = idx.shape[0]
        daily_tests = num_tests * test_days
        # Get the subtarger to exclude
        exclude_idx = build_subtarget_list(i, rand_large, to_ndarray = False)
        exclude_idx += rand_others_all
        exclude_idx.sort()
        exclude_idx = np.asarray(exclude_idx, dtype=np.int32)
        # And the matching zeros to zero-out their probability of being tested
        exclude_vals = np.zeros(exclude_idx.shape, dtype=np.float32)
        intervs.append(cv.test_num(daily_tests=daily_tests,
                                   subtarget={'inds':exclude_idx, 'vals':exclude_vals},
                                   **kwargs))

    for i,idx in enumerate(rand_other):
        # Set up a multiplier for the tests per day argument.
        test_days = ((i + np.arange(num_days)) % test_period < 1).astype(np.int32)
        num_tests = idx.shape[0]
        daily_tests = num_tests * test_days
        # Get the subtarger to exclude
        exclude_idx = build_subtarget_list(i, rand_other, to_ndarray = False)
        exclude_idx += rand_large_all
        exclude_idx.sort()
        exclude_idx = np.asarray(exclude_idx, dtype=np.int32)
        # And the matching zeros to zero-out their probability of being tested
        exclude_vals = np.zeros(exclude_idx.shape, dtype=np.float32)
        intervs.append(cv.test_num(daily_tests=daily_tests,
                                   subtarget={'inds':exclude_idx, 'vals':exclude_vals},
                                   **kwargs))
        
    # and that's all, folks
    return intervs 
        


# contact tracing with sensitivity and specificity
class contact_tracing_sens_spec(cv.Intervention):
    '''
    Contact tracing of positive people.

    Args:
        trace_sensitivity (dict): probability of tracing for infeced contacts, per layer
        trace_specificity (dict): probability of not tracing for not infeced contacts, per layer
        trace_time  (dict): days required to trace, per layer
        start_day   (int):  intervention start day
        end_day     (int):  intervention end day
        test_delay  (int):  number of days a test result takes
        presumptive (bool): whether or not to begin isolation and contact tracing on the presumption of a positive diagnosis
        kwargs      (dict): passed to Intervention()
    '''
    def __init__(self, trace_sensitivity=None, trace_specificity=None, trace_time=None, start_day=0, end_day=None, presumptive=False, **kwargs):
        super().__init__(**kwargs) # Initialize the Intervention object
        self._store_args() # Store the input arguments so the intervention can be recreated
        self.trace_sensitivity = trace_sensitivity
        self.trace_specificity = trace_specificity
        self.trace_time  = trace_time
        self.start_day   = start_day
        self.end_day     = end_day
        self.presumptive = presumptive
        return


    def initialize(self, sim):
        ''' Fix the dates and dictionaries '''
        self.start_day = sim.day(self.start_day)
        self.end_day   = sim.day(self.end_day)
        self.days      = [self.start_day, self.end_day]
        if self.trace_sensitivity is None:
            self.trace_sensitivity = 1.0
        if self.trace_specificity is None:
            self.trace_specificity = 1.0    
        if self.trace_time is None:
            self.trace_time = 0.0
        if sc.isnumber(self.trace_sensitivity):
            val = self.trace_sensitivity
            self.trace_sensitivity = {k:val for k in sim.people.layer_keys()}
        if sc.isnumber(self.trace_specificity):
            val = self.trace_specificity
            self.trace_specificity = {k:val for k in sim.people.layer_keys()}
        if sc.isnumber(self.trace_time):
            val = self.trace_time
            self.trace_time = {k:val for k in sim.people.layer_keys()}
        self.initialized = True
        return


    def apply(self, sim):
        t = sim.t
        if t < self.start_day:
            return
        elif self.end_day is not None and t > self.end_day:
            return

        # Figure out whom to test and trace
        if not self.presumptive:
            trace_from_inds = cvu.true(sim.people.date_diagnosed == t) # Diagnosed this time step, time to trace
        else:
            just_tested = cvu.true(sim.people.date_tested == t) # Tested this time step, time to trace
            trace_from_inds = cvu.itruei(sim.people.exposed, just_tested) # This is necessary to avoid infinite chains of asymptomatic testing

        if len(trace_from_inds): # If there are any just-diagnosed people, go trace their contacts
           traceable_layers = {k:v for k,v in self.trace_sensitivity.items() if v != 0.} # Only trace if there's a non-zero tracing probability
           for lkey,this_trace_sensitivity in traceable_layers.items():
               if sim.people.pars['beta_layer'][lkey]: # Skip if beta is 0 for this layer
                   this_trace_time = self.trace_time[lkey]

                # Find all the contacts of these people
                   inds_list = []
                   for k1,k2 in [['p1','p2'],['p2','p1']]: # Loop over the contact network in both directions -- k1,k2 are the keys
                       in_k1 = np.isin(sim.people.contacts[lkey][k1], trace_from_inds).nonzero()[0] # Get all the indices of the pairs that each person is in
                       inds_list.append(sim.people.contacts[lkey][k2][in_k1]) # Find their pairing partner
                   edge_inds = np.unique(np.concatenate(inds_list)) # Find all edges
                
                # Find all infected contacts
                   edge_inds_infected = np.array([],dtype=np.int32)
                   for x,y in enumerate(edge_inds):
                       if not sim.people.susceptible[x]:
                           edge_inds_infected = np.append(edge_inds_infected,y)
                

                # Check contacts
                   contact_inds_infected = cvu.binomial_filter(this_trace_sensitivity, edge_inds_infected) # Filter the indices according to the probability of being able to trace this layer
                   if len(contact_inds_infected):
                       sim.people.known_contact[contact_inds_infected] = True
                       sim.people.date_known_contact[contact_inds_infected]  = np.fmin(sim.people.date_known_contact[contact_inds_infected], sim.people.t+this_trace_time)

           traceable_layers = {k:v for k,v in self.trace_specificity.items() if v != 1.} # Only trace if there's a non-zero tracing probability
           for lkey,this_trace_specificity in traceable_layers.items():
               if sim.people.pars['beta_layer'][lkey]: # Skip if beta is 0 for this layer
                   this_trace_time = self.trace_time[lkey]

                # Find all the contacts of these people
                   inds_list = []
                   for k1,k2 in [['p1','p2'],['p2','p1']]: # Loop over the contact network in both directions -- k1,k2 are the keys
                       in_k1 = np.isin(sim.people.contacts[lkey][k1], trace_from_inds).nonzero()[0] # Get all the indices of the pairs that each person is in
                       inds_list.append(sim.people.contacts[lkey][k2][in_k1]) # Find their pairing partner
                   edge_inds = np.unique(np.concatenate(inds_list)) # Find all edges
                

                # Find all not infected contacts
                   edge_inds_noninfected = np.array([],dtype=np.int32)
                   for x,y in enumerate(edge_inds):
                       if sim.people.susceptible[x]:
                           edge_inds_noninfected = np.append(edge_inds_noninfected,y)

                # Check contacts
                   contact_inds_noninfected = cvu.binomial_filter(1-this_trace_specificity, edge_inds_noninfected) # Filter the indices according to the probability of being able to trace this layer                                               
     
                   if len(contact_inds_noninfected):
                       sim.people.known_contact[contact_inds_noninfected] = True
                       sim.people.date_known_contact[contact_inds_noninfected]  = np.fmin(sim.people.date_known_contact[contact_inds_noninfected], sim.people.t+this_trace_time)

        return
    