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
from .exception import TerrierException
import collections
 
IMPORT_EXOGENOUS = np.array([0.00000000, 0.00000000, 0.28571429, 0.13222111, 0.00000000, 0.00000000 ,0.00000000 ,0.07142857, 0.07114942,
                             0.21406296, 0.40228266, 0.60476382, 0.00000000 ,0.19232560])

__all__=['test_post_quar', 'gen_periodic_testing_interventions','gen_large_dorm_testing_interventions',
         'test_household_when_pos','pulsed_infections_network','gen_undergrad_testing_interventions','contact_tracing_sens_spec',
         'import_infections_network','import_infections_percent_network','pulsed_infections_diffuse','gen_periodic_testing_interventions_real',
         'IMPORT_EXOGENOUS','import_infections_exogenous_network','infect_specific_people','contact_tracing_sens_spec_log']

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
        interv = bu.pulsed_infections(BU_pop['contacts'],days_arr, layer_name='household', pulse_count=10)
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
        hosp_max = sim.people.count('severe')   > sim['n_beds_hosp'] if sim['n_beds_hosp'] else False # Check for acute bed constraint
        icu_max  = sim.people.count('critical') > sim['n_beds_icu']  if sim['n_beds_icu']  else False # Check for ICU bed constraint
        sim.people.infect(inds=np.array(targets, dtype = np.int32),hosp_max = hosp_max,icu_max  = icu_max,layer='importation')
        #sim.people.exogenous[np.array(targets, dtype = np.int32)] = np.full(len(targets),1,np.int32)
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
def gen_undergrad_testing_interventions(BU_pop, num_days, test_period_large = 3,  test_period=7, **kwargs):
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
    uids_large = BU_pop['uid'][np.where(BU_pop['undergrad'] > 0)].copy()
    uids_other = BU_pop['uid'][np.where(BU_pop['undergrad'] < 1)].copy()
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
        
#%%
def gen_periodic_testing_interventions_real(BU_pop, num_days,
                                             test_dist={0: 0.12, 1: 0.152, 2: 0.152, 3: 0.152, 4: 0.152, 5: 0.152, 6: 0.12},
                                             start_date='2020-09-02',
                                             **kwargs):
    ''' Test undergrads every 3-4 days
        Everyone else is tested 1x/week
        
        This is a variation on gen_undergrad_testing_interventions.  For that one the testing is evenly
        spread out around the week.  This one sets a testing target by the day of the week, tries to 
        meet that target, and then spreads out any remaining testing evenly.  The idea is a testing
        distribution that matches the one seen for real where there are more tests M-F than on the weekends.
        
        BU_pop:  BU population dictionary
        num_days: length of simulation
        
        test_dist: Dictionary for days of the week that specifies the target
            distribution. Keys are day numbers (0 is Sunday, etc). Sum of test_dist
            should be 1. 
            
        start_date : date strings of the intervention range. Req'd to get day of week calculations done.
        ...and also all the regular arguments to cv.test_num
        
        Returns: A list of cv.test_num interventions
    '''
    # Verify that test_dist sums close to 1
    if sum(test_dist.values())-1 > 1e-6:
        raise TerrierException('test_dist failed to sum to a value close to 1.')
    
    # Split the uids into 2 populations - frequent testes and everyone else who's test 1/week.
    # As this was originally for large dorms wherever it says "large" 
    # think "frequent"
    # 96 hours == 4 days
    # 168 hours == 7 days
    # Ignore anyone else.
    uids_large = BU_pop['uid'][np.where(BU_pop['labTestEveryNHours'] == 96)].copy()
    uids_other = BU_pop['uid'][np.where(BU_pop['labTestEveryNHours'] == 168)].copy()
    base_excludes = BU_pop['uid'][np.isin(BU_pop['labTestEveryNHours'],[96,168],invert=True)].copy()
    
    # shuffle 'em up
    np.random.shuffle(uids_large)
    np.random.shuffle(uids_other)

    # This is fixed to make the programming quicker.  The undergrads get 
    # tested every 3-4 days.  Everyone else is 1x/week. Most testing must
    # happen on weekedays per the test_dist distribution.
    # undergrad testing groups:
    # M-R, T-F, W-Sa, Su-W.
    thresh = 0.33 # Less on the weekend but wednesday is the same as the other days of the week.
    # A threshold of 0.33 balances things out without creating a jump or dip on Wednesday.
    # loop through the large group, and assign them to the testing groups.
    # In hindsight this is probably not the most straightforward way to handle 
    # this but it's too much effort to fix at the moment.
    large_test_groups = {'m_r':[],
                    't_f':[],
                    'w_s':[],
                    'u_w':[],}
    group_lut={'m_r':[0,1,0,0,1,0,0],
               't_f':[0,0,1,0,0,1,0],
               'w_s':[0,0,0,1,0,0,1],
               'u_w':[1,0,0,1,0,0,0],}
    for uid in uids_large:
        if np.random.uniform() < thresh:
            # weekend group
            if np.random.uniform() < 0.5:
                large_test_groups['w_s'].append(uid)
            else:
                large_test_groups['u_w'].append(uid)
        else:
            if np.random.uniform() < 0.5:
                large_test_groups['m_r'].append(uid)
            else:
                large_test_groups['t_f'].append(uid)            
    # Now for everyone who is NOT an undergraduate test 1x/week.
    # Again bias the testing towards weekdays. Sort everyone into
    # 7 testing bins.
    weekly_test_groups = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: []} 
    # This can be sorted using numpy binning...
    # retrieve daily fractions:
    daily = np.array([0] + list(test_dist.values()))
    # Wednesday is over-sampled due to the two weekend undergrad groups
    # using it.  Reduce it by 25%
    #daily[4] *= 0.75
    # Now renormalize the daily distribution and create sorting bins
    bins = np.cumsum(daily / np.sum(daily))
    # make the last bin a wee higher than 1 just to make sure the sorting
    # works if a value 1.0 is pulled from the RNG
    bins[-1] += 0.005
    dice = np.random.uniform(size = len(uids_other))
    # now find indices for all uids_other for the group they belong to.
    inds = np.digitize(dice,bins) - 1
    # For every uid, sort into a group using the inds
    for uid,idx in zip(uids_other,inds):
        weekly_test_groups[idx].append(uid)
    
    # A list of Interventions    
    intervs = []
    
    # Do the frequent testers first.  We need to figure out the day of the 
    # week from the start date and the number of days.

    # Python returns this as Monday=0. We used Sunday=0. Fix that.
    ring = collections.deque([0,1,2,3,4,5,6],maxlen=7)
    ring.rotate(1)
    ring = tuple(ring) # Convert to a ring    
    
    for group in large_test_groups:
        # for each day of the simulation fill into an ndarray 
        # whether or not we're testing today.
        test_days = np.zeros(num_days,dtype=np.int32)
        cur_day = make_dt(start_date) 
        for i in range(num_days):
            weekday = ring.index(cur_day.weekday())
            if group_lut[group][weekday]:
                test_days[i] = 1
            cur_day += datetime.timedelta(days=1)
        # Multiple by the size of this group
        daily_tests = test_days * len(large_test_groups[group])
        # Build a list of every other person who is not in this subgroup
        # and exclude them.
        excludes = list(base_excludes)
        for subgroup in large_test_groups:
            if subgroup == group:
                continue
            excludes += large_test_groups[subgroup]
        # and now the other test groups
        for w in weekly_test_groups:
            excludes += weekly_test_groups[w]
        excludes.sort()
        excludes = np.array(excludes)
        exclude_vals = np.zeros(excludes.shape,dtype=np.float32)
        # And now create the intervention.
        intervs.append(cv.test_num(daily_tests=daily_tests,
                                   subtarget={'inds':excludes, 'vals':exclude_vals},
                                   **kwargs))
    # And now do something similar for the 1/week people
    for day in weekly_test_groups:
        test_days = np.zeros(num_days,dtype=np.int32)
        cur_day = make_dt(start_date) 
        # this could be vectorized but this works
        for i in range(num_days):
            weekday = ring.index(cur_day.weekday())
            if weekday == day:
                test_days[i] = 1
            cur_day += datetime.timedelta(days=1)        
        daily_tests = test_days * len(large_test_groups[group])
        excludes = list(base_excludes)
        # exclude anyone in the high frequency group
        for group in large_test_groups:
            excludes += large_test_groups[subgroup]
        # and any weekly testers who test on a different day.
        for w in weekly_test_groups:
            if w == day:
                continue
            excludes += weekly_test_groups[w]
        excludes.sort()
        excludes = np.array(excludes)
        exclude_vals = np.zeros(excludes.shape,dtype=np.float32)
        # And now create the intervention.
        intervs.append(cv.test_num(daily_tests=daily_tests,
                                   subtarget={'inds':excludes, 'vals':exclude_vals},
                                   **kwargs))
    return intervs 

#%%

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
                   edge_inds_infected = [] #np.array([],dtype=np.int32)
                   for x,y in enumerate(edge_inds):
                       if not sim.people.susceptible[y]:
                           edge_inds_infected.append(y) # = np.append(edge_inds_infected,y)
                   edge_inds_infected = np.array(edge_inds_infected, dtype=np.int32)

                # Check contacts
                   contact_inds_infected = cvu.binomial_filter(this_trace_sensitivity, edge_inds_infected) # Filter the indices according to the probability of being able to trace this layer
                   if len(contact_inds_infected):
                       sim.people.known_contact[contact_inds_infected] = True
                       sim.people.date_known_contact[contact_inds_infected]  = np.fmin(sim.people.date_known_contact[contact_inds_infected], sim.people.t+this_trace_time)

           traceable_layers = {k:v for k,v in self.trace_specificity.items() if v != 0.} # Only trace if there's a non-zero tracing probability
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
                   edge_inds_noninfected = [] #np.array([],dtype=np.int32)
                   for x,y in enumerate(edge_inds):
                       if sim.people.susceptible[y]:
                           edge_inds_noninfected.append(y) # = np.append(edge_inds_noninfected,y)
                   edge_inds_noninfected = np.array(edge_inds_noninfected)
                # Check contacts
                   contact_inds_noninfected = cvu.binomial_filter(1-this_trace_specificity, edge_inds_noninfected) # Filter the indices according to the probability of being able to trace this layer                                               
     
                   if len(contact_inds_noninfected):
                       sim.people.known_contact[contact_inds_noninfected] = True
                       sim.people.date_known_contact[contact_inds_noninfected]  = np.fmin(sim.people.date_known_contact[contact_inds_noninfected], sim.people.t+this_trace_time)

        return
  
  

class contact_tracing_sens_spec_log(contact_tracing_sens_spec):
    ''' modified, non-vectorized version of contact_tracing_sens_spec that
        creates a per-person contact tracing log per day.  Slower and uses
        a wee bit more memory, but required for the BU_quarantine_network_count snapshot.
        
        Adds the quarantining of those who are contact traced. This is a complete replacement for 
        the regular covasim contact_trace() intervention.
        '''
    def __init__(self, trace_sensitivity=None, trace_specificity=None, trace_time=None, start_day=0, end_day=None, presumptive=False, **kwargs):
        super().__init__(**kwargs) # Initialize the parent object
        # Every day the intervention runs record all contact tracing.  This is
        # a dictionary indexed by day number.
        self.contact_log = {}
        
        
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

        contact_dt = {}
        if trace_from_inds.size == 0: 
            return # nothing to do
        # If there are any just-diagnosed people, go trace their contacts
        traceable_layers = {k:v for k,v in self.trace_sensitivity.items() if v != 0.} # Only trace if there's a non-zero tracing probability

        for lkey,this_trace_sensitivity in traceable_layers.items():
           if sim.people.pars['beta_layer'][lkey] == 0: 
                # Skip if beta is 0 for this layer
                continue
           this_trace_time = self.trace_time[lkey]

           # Find all the contacts of these people. Non-vectorized search.
           inds_list = []
           for k1,k2 in [['p1','p2'],['p2','p1']]: # Loop over the contact network in both directions -- k1,k2 are the keys
               for ti in trace_from_inds:
                   if ti not in contact_dt:
                       contact_dt[ti]={'all':[],'contacts':[]}
                   in_k1 = np.isin(sim.people.contacts[lkey][k1], [ti]).nonzero()[0] # Get all the indices of the pairs that each person is 
                   tmp = sim.people.contacts[lkey][k2][in_k1] # Find their pairing partner
                   if tmp.size > 0:
                       contact_dt[ti]['all'] += list(np.unique(tmp))
                   inds_list.append(tmp) 
           edge_inds = np.unique(np.concatenate(inds_list)) # Find all edges
        
           
           # Find all infected contacts
           edge_inds_infected = [] #np.array([],dtype=np.int32)
           for x,y in enumerate(edge_inds):
               if not sim.people.susceptible[y]:
                   edge_inds_infected.append(y) # = np.append(edge_inds_infected,y)
           edge_inds_infected = np.array(edge_inds_infected, dtype=np.int32)

           # Check contacts
           contact_inds_infected = cvu.binomial_filter(this_trace_sensitivity, edge_inds_infected) # Filter the indices according to the probability of being able to trace this layer
           if len(contact_inds_infected):
               sim.people.known_contact[contact_inds_infected] = True
               sim.people.date_known_contact[contact_inds_infected]  = np.fmin(sim.people.date_known_contact[contact_inds_infected], sim.people.t+this_trace_time)
               # Schedule quarantine for the notified people to start on the date they will be notified. Note that the quarantine duration is based on the time since last contact, rather than time since notified
               sim.people.schedule_quarantine(contact_inds_infected, sim.t+this_trace_time, sim.people.pars['quar_period']-this_trace_time) 
               
               # Take the contact_inds_infected. For each person in contact_dt check their contact list to see if someone
               # was contacted. If so add that person to their contact_inf list.
               for ind in contact_dt:
                  tmp = np.unique(np.array(contact_dt[ind]['all']))
                  contacted = tmp[np.where(np.isin(tmp,contact_inds_infected))] 
                  if contacted.size > 0:
                      contact_dt[ind]['contacts'].append(contacted)

        traceable_layers = {k:v for k,v in self.trace_specificity.items() if v != 0.} # Only trace if there's a non-zero tracing probability
        
        # Reset the "all" list on contact_dt for any existing keys
        for ind in contact_dt:
            contact_dt[ind]['all'] = []
        
        for lkey,this_trace_specificity in traceable_layers.items():
           if sim.people.pars['beta_layer'][lkey] ==0: 
                # Skip if beta is 0 for this layer
                continue
           this_trace_time = self.trace_time[lkey]

           # Find all the contacts of these people
           inds_list = []
           for k1,k2 in [['p1','p2'],['p2','p1']]: # Loop over the contact network in both directions -- k1,k2 are the keys
               for ti in trace_from_inds:
                    if ti not in contact_dt:
                        contact_dt[ti]={'all':[],'contacts':[]}
                    in_k1 = np.isin(sim.people.contacts[lkey][k1], [ti]).nonzero()[0] # Get all the indices of the pairs that each person is 
                    tmp = sim.people.contacts[lkey][k2][in_k1] # Find their pairing partner
                    if tmp.size > 0:
                        contact_dt[ti]['all'] += list(np.unique(tmp))
                    inds_list.append(tmp) 
           edge_inds = np.unique(np.concatenate(inds_list)) # Find all edges
        

           # Find all not infected contacts
           edge_inds_noninfected = [] #np.array([],dtype=np.int32)
           for x,y in enumerate(edge_inds):
               if sim.people.susceptible[y]:
                   edge_inds_noninfected.append(y) # = np.append(edge_inds_noninfected,y)
           edge_inds_noninfected = np.array(edge_inds_noninfected)
           # Check contacts
           contact_inds_noninfected = cvu.binomial_filter(1-this_trace_specificity, edge_inds_noninfected) # Filter the indices according to the probability of being able to trace this layer                                               
 
           if len(contact_inds_noninfected):
               sim.people.known_contact[contact_inds_noninfected] = True
               sim.people.date_known_contact[contact_inds_noninfected]  = np.fmin(sim.people.date_known_contact[contact_inds_noninfected], sim.people.t+this_trace_time)
               sim.people.schedule_quarantine(contact_inds_noninfected, sim.t+this_trace_time, sim.people.pars['quar_period']-this_trace_time) 

               # Take the contact_inds_noninfected. For each person in contact_dt check their contact list to see if someone
               # was contacted. If so add that person to their contact_inf list.
               for ind in contact_dt:
                  tmp = np.unique(np.array(contact_dt[ind]['all']))
                  contacted = tmp[np.where(np.isin(tmp,contact_inds_noninfected))] 
                  if contacted.size > 0:
                      contact_dt[ind]['contacts'].append(contacted)
                      
        # Strip the "all" list from the contact_dt as it's no longer needed.
        for ind in contact_dt:
            del contact_dt[ind]['all']

        # Append today's contact_dt to the contact_log
        self.contact_log[sim.t] = contact_dt
        return    

  
class infect_specific_people(cv.Intervention):
    '''
    Infects a specified set of people on specified days.
    '''
    
    def __init__(self, people_to_infect,  **kwargs):
        ''' people_to_infect: A dictionary of { days: [people_id,...]}
            days should be integers. '''
        super().__init__(**kwargs) # Initialize the Intervention object
        self._store_args() # Store the input arguments so the intervention can be recreated
        self.people_to_infect = people_to_infect
        return


    def initialize(self, sim):
        self.initialized = True
        print('INFET INIT')
        return

    def apply(self, sim):
        # Check to see if today is a day we infect some unfortunate souls
        if not sim.t in self.people_to_infect:
            return

        # Get the list of people to infect. Call the infect() method.
        inds = np.array(self.people_to_infect[sim.t])        
        # Check for acute bed constraint     
        hosp_max = sim.people.count('severe')   > sim['n_beds_hosp'] if sim['n_beds_hosp'] else False 
        # Check for ICU bed constraint
        icu_max  = sim.people.count('critical') > sim['n_beds_icu']  if sim['n_beds_icu']  else False 
        # And....that's all!  The simulation will roll the dice on the
        # people infected today as to their outcomes etc. 
        sim.people.infect(inds=inds,hosp_max=hosp_max, icu_max=icu_max, layer='importation')
        return
    


class import_infections_network(cv.Intervention):
    '''
    Infects a set of people across networks.  
    '''
    
    def __init__(self, days, import_count = 1, import_sampling = 'poisson', **kwargs):
        super().__init__(**kwargs) # Initialize the Intervention object
        self._store_args() # Store the input arguments so the intervention can be recreated
        if import_count >= 1:
            self.import_count = import_count
        else:
            raise Exception('import_count must be >= 1')
        self.days = np.array(days, dtype = np.bool)
        self.import_sampling = import_sampling.strip().lower()
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

    def apply(self, sim):
        # If today is not a pulse day then return.
        if not self.days[sim.t]:
            return
        
        # Use either constant infection count or a Poisson average
        import_count = self.import_count
        if self.import_sampling == 'poisson':
            import_count = cvu.poisson(self.import_count)  
        hosp_max = sim.people.count('severe')   > sim['n_beds_hosp'] if sim['n_beds_hosp'] else False # Check for acute bed constraint
        icu_max  = sim.people.count('critical') > sim['n_beds_icu']  if sim['n_beds_icu']  else False # Check for ICU bed constraint
        #susceptible_inds = cvu.true(sim.people.susceptible)   
        #susceptible_importation_inds = cvu.choose(max_n=len(susceptible_inds), n=import_count)
        #importation_inds = susceptible_inds[susceptible_importation_inds]
        importation_inds = cvu.choose(max_n=len(sim.people), n=import_count)
        sim.people.infect(inds=importation_inds,hosp_max=hosp_max, icu_max=icu_max, layer='importation')
        return
    
    @staticmethod
    def get_days_arr(start_day, end_day, import_day):
        ''' Returns a numpy array for the days to apply a pulse.
        start_day, end_day are date strings.
        import_day is th enumber of the day of the week to apply
        the import.  Sunday=0.'''
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
            if get_daynum(today) == import_day:
                days_arr[i] = 1
        return days_arr




class import_infections_percent_network(cv.Intervention):
    '''
    Infects a set of people across networks.  
    '''
    
    def __init__(self, days, import_count = 20, import_percent = .5,import_sampling = 'poisson', **kwargs):
        super().__init__(**kwargs) # Initialize the Intervention object
        self._store_args() # Store the input arguments so the intervention can be recreated
        if import_count >= 1:
            self.import_count = import_count
        else:
            raise Exception('import_count must be >= 1')
        self.import_percent = import_percent
        self.days = np.array(days, dtype = np.bool)
        self.import_sampling = import_sampling.strip().lower()
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

    def apply(self, sim):
        # If today is not a pulse day then return.
        if not self.days[sim.t]:
            return
        
        # Use either constant infection count or a Poisson average
        import_count_oncampus = self.import_count*self.import_percent
        import_count_offcampus = self.import_count-import_count_oncampus
        if self.import_sampling == 'poisson':
            import_count_oncampus = cvu.poisson(import_count_oncampus)  
            import_count_offcampus = cvu.poisson(import_count_offcampus)
        hosp_max = sim.people.count('severe')   > sim['n_beds_hosp'] if sim['n_beds_hosp'] else False # Check for acute bed constraint
        icu_max  = sim.people.count('critical') > sim['n_beds_icu']  if sim['n_beds_icu']  else False # Check for ICU bed constraint

        susceptible_check = sim.people.susceptible
        # for on-campus
        at_bu = sim.people.campResident > 0
        at_bu_inds = cvu.true(at_bu & susceptible_check)   
        importation_inds = cvu.choose(max_n=len(at_bu_inds), n=import_count_oncampus )
        sim.people.infect(inds=at_bu_inds[importation_inds],hosp_max=hosp_max, icu_max=icu_max, layer='importation')
        # for off-campus student
        not_at_bu = sim.people.campResident < 1
        student = sim.people.category < 2
        not_at_bu_student = not_at_bu & student
        not_at_bu_inds = cvu.true(not_at_bu_student & susceptible_check)
        importation_inds = cvu.choose(max_n=len(not_at_bu_inds), n=import_count_offcampus )
        sim.people.infect(inds=not_at_bu_inds[importation_inds],hosp_max=hosp_max, icu_max=icu_max, layer='importation')
        
        return
    
    @staticmethod
    def get_days_arr(start_day, end_day, import_day):
        ''' Returns a numpy array for the days to apply a pulse.
        start_day, end_day are date strings.
        import_day is th enumber of the day of the week to apply
        the import.  Sunday=0.'''
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
            if get_daynum(today) == import_day:
                days_arr[i] = 1
        return days_arr
        
        
class pulsed_infections_diffuse(pulsed_infections_network):
    ''' A modification to pulsed_infections_network.  Instead of trying to infect whole 
        households or roommates together, infect the layer in a diffuse manner. '''
    
    def apply(self, sim):
        # If today is not a pulse day then return.
        if not self.days[sim.t]:
            return
        # set of those who will be infected
        targets = set()

        iters = 0
        # Index list for the networks
        net_inds = list(range(len(self.networks)))
        
        # Use either constant infection count or a Poisson average
        pulse_count = self.pulse_count
        if self.pulse_sampling == 'poisson':
            pulse_count = np.random.poisson(self.pulse_count)  

        # limit the number of loop iterations for safety
        max_loops = 500 * pulse_count
        
        while len(targets) < pulse_count and iters < max_loops:
            # Pick a household/roommates/etc.
            net = np.random.choice(net_inds)
            # Now pick a random person from that network.
            person = np.random.choice(self.networks[net])
            # If they're not in the set and susceptible add them
            # to the target infection list.
            if person not in targets and sim.people.susceptible[person]:
                targets.add(person)

            iters += 1
        targets = list(targets)
        # Now infect those targets
        sim.people.infect(inds=np.array(targets, dtype = np.int32))
        return    


class import_infections_exogenous_network(cv.Intervention):
    '''
    Infects a set of people across networks.  
    '''
    
    def __init__(self, days, import_count = 1, import_exogenous = IMPORT_EXOGENOUS,import_sampling = 'poisson', **kwargs):
        super().__init__(**kwargs) # Initialize the Intervention object
        self._store_args() # Store the input arguments so the intervention can be recreated
        if import_count >= 1:
            self.import_count = import_count
        else:
            raise Exception('import_count must be >= 1')
        self.import_exogenous = import_exogenous
        self.days = np.array(days, dtype = np.bool)
        self.import_sampling = import_sampling.strip().lower()
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

    def apply(self, sim):
        # If today is not a pulse day then return.
        if not self.days[sim.t]:
            return
        
        # Use either constant infection count or a Poisson average
        import_count = self.import_count
        if self.import_sampling == 'poisson':
            import_count = cvu.poisson(import_count)  
        hosp_max = sim.people.count('severe')   > sim['n_beds_hosp'] if sim['n_beds_hosp'] else False # Check for acute bed constraint
        icu_max  = sim.people.count('critical') > sim['n_beds_icu']  if sim['n_beds_icu']  else False # Check for ICU bed constraint

        susceptible_check = sim.people.susceptible
        # rescale import_exogenous to to be sum = 1
        import_exogenous = self.import_exogenous
        multinomial_p = import_exogenous/import_exogenous.sum()
        # multinomial
        import_subgroup = np.random.multinomial(import_count, multinomial_p, size=1)
        
        # define subgroup 
        at_BUMC = sim.people.campus > 1
        at_CRC = sim.people.campus < 2
        Affiliate = sim.people.category > 3 
        Faculty = (sim.people.category < 3) & (sim.people.category >1)
        Staff = ( sim.people.category < 4) & (sim.people.category > 2)
        Student  = sim.people.category < 2 
        Oncampus = sim.people.campResident > 0
        Offcampus = sim.people.campResident < 1
        Grad = sim.people.undergrad < 1
        Undergrad = sim.people.undergrad > 0
 
        # define exclusive subgroup 
        Affiliate_BUMC = Affiliate & at_BUMC & susceptible_check
        Faculty_BUMC = Faculty & at_BUMC & susceptible_check
        Staff_BUMC = Staff & at_BUMC & susceptible_check
        Grad_offcampus_BUMC = Student & Grad & Offcampus & at_BUMC & susceptible_check
        Undergrad_offcampus_BUMC = Student & Undergrad & Offcampus & at_BUMC & susceptible_check
        Grad_oncampus_BUMC = Student & Grad & Oncampus & at_BUMC & susceptible_check
        Undergrad_oncampus_BUMC = Student & Undergrad & Oncampus & at_BUMC & susceptible_check
        
        Affiliate_CRC = Affiliate & at_CRC & susceptible_check
        Faculty_CRC = Faculty & at_CRC & susceptible_check
        Staff_CRC = Staff & at_CRC & susceptible_check
        Grad_offcampus_CRC = Student & Grad & Offcampus & at_CRC & susceptible_check
        Undergrad_offcampus_CRC = Student & Undergrad & Offcampus & at_CRC & susceptible_check
        Grad_oncampus_CRC = Student & Grad & Oncampus & at_CRC & susceptible_check
        Undergrad_oncampus_CRC = Student & Undergrad & Oncampus & at_CRC & susceptible_check
        
        subgroup_list =[Affiliate_BUMC,Faculty_BUMC,Staff_BUMC,Grad_offcampus_BUMC,Undergrad_offcampus_BUMC,Grad_oncampus_BUMC,Undergrad_oncampus_BUMC,
                        Affiliate_CRC,Faculty_CRC,Staff_CRC,Grad_offcampus_CRC,Undergrad_offcampus_CRC,Grad_oncampus_CRC,Undergrad_oncampus_CRC]
        
        # infect people
        for i in range(import_subgroup.shape[1]):
            tmp_inds = cvu.true(subgroup_list[i])   
            if len(tmp_inds) > 0: 
                importation_inds = cvu.choose(max_n=len(tmp_inds), n=min(import_subgroup[0,i],len(tmp_inds)) )
                sim.people.infect(inds=tmp_inds[importation_inds],hosp_max=hosp_max, icu_max=icu_max, layer='importation')
    
        
        return
    
    @staticmethod
    def get_days_arr(start_day, end_day, import_day):
        ''' Returns a numpy array for the days to apply a pulse.
        start_day, end_day are date strings.
        import_day is th enumber of the day of the week to apply
        the import.  Sunday=0.'''
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
            if get_daynum(today) == import_day:
                days_arr[i] = 1
        return days_arr
