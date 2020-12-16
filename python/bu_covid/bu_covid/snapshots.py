## ---------------------------
##
## Snapshot analyzers for BU specific info
##
## Authors: Brian Gregor, Wenrui Li 
##          Boston University
##
## Date Created: 2020-07-31
##
## Email: bgregor@bu.edu
##
## ---------------------------

 
# Default number of days that quarantine rooms are cleaned for after
# someone leaves quarantine:
QUAR_CLEANING_DAYS = 1   
ISO_CLEANING_DAYS = 1 
ISO_DAYS = 10

__all__=['snapshots_to_df','get_BU_snapshots','infection_count_to_df',
         'diag2iso_count_to_df','severe_count_to_df','critical_count_to_df',
         'dead_count_to_df','diagnosed_count_to_df','recovered_count_to_df',
         'quarantined_count_to_df','quarantined_end_count_to_df','QUAR_CLEANING_DAYS',   
         'ISO_CLEANING_DAYS','ISO_DAYS','safe_nan_compare','BU_pos_test_date_count',
         'BU_quarantine_network_count']

import covasim as cv
import sciris as sc
import numpy as np
import pandas as pd
from collections import deque

import itertools as it 
import operator as op

from .interventions import contact_tracing_sens_spec_log

def safe_nan_compare(x,y,op):
    ''' do an operation "op" between x & y. np.nan values
        always return a False.'''
    result = np.full_like(x,False,dtype=np.bool)
    if isinstance(x,np.ndarray) and isinstance(y,np.ndarray):
        good_ind= np.isfinite(x) & np.isfinite(y)
        result[good_ind] = op(x[good_ind],y[good_ind])
    elif isinstance(x,np.ndarray):
        good_ind= np.isfinite(x)
        result[good_ind] = op(x[good_ind],y)
    else:
        good_ind= np.isfinite(y)
        result[good_ind] = op(x,y[good_ind])
    return result 

class BU_res_quarantine_count(cv.Analyzer):
    '''
    Analyzer that takes a "snapshot" of n_quarantined people who are in quarantine

        sim = cv.Sim(analyzers=BU_res_quarantine_count('2020-04-04', '2020-04-14'))
        sim.run()
        snapshot = sim['analyzers'][0]
        people = snapshot.snapshots[0]            # Option 1
        people = snapshot.snapshots['2020-04-04'] # Option 2
        people = snapshot.get('2020-04-14')       # Option 3
        people = snapshot.get(34)                 # Option 4
        people = snapshot.get()                   # Option 5
        
        The results are stored in a dictionary.
    '''

    def __init__(self, days, *args, **kwargs):
        super().__init__(**kwargs) # Initialize the Analyzer object
        days = sc.promotetolist(days) # Combine multiple days
        days.extend(args) # Include additional arguments, if present
        self.days      = days # Converted to integer representations
        self.dates     = None # String representations
        self.start_day = None # Store the start date of the simulation
        self.snapshots = None # Store the actual snapshots
        return


    def initialize(self, sim):
        self.start_day = sim['start_day'] # Store the simulation start day
        self.days = cv.interventions.process_days(sim, self.days) # Ensure days are in the right format
        self.dates = [sim.date(day) for day in self.days] # Store as date strings
        self.initialized = True
        self.snapshots = {} # Store the snapshots in a dictionary.
        self.yesterday_quarantine = sim.people.quarantined.copy()
        return


    def apply(self, sim):
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            # if sim didn't start on day 0 back-fill in missing days with 0.
            if sim.t > 0 and len(self.snapshots) == 0:
                for i in range(sim.t):
                    self.snapshots[self.dates[i]] = 0

            # campResident is: 0 off campus, 1 on campus, 2 large dorm on campus
            #sim.people.quarantined is a Boolean array of everyone who is quarantined
            self.snapshots[date] = np.sum(np.logical_and(sim.people.quarantined,sim.people.campResident > 0)) 
        return


    def get(self, key=None):
        ''' Retrieve a snapshot from the given key (int, str, or date) '''
        if key is None:
            key = self.days[0]
        day  = cv.misc.day(key, start_day=self.start_day)
        date = cv.misc.date(day, start_date=self.start_day, as_date=False)
        if date in self.snapshots:
            snapshot = self.snapshots[date]
        else:
            dates = ', '.join(list(self.snapshots.keys()))
            errormsg = f'Could not find snapshot date {date} (day {day}): choices are {dates}'
            raise sc.KeyNotFoundError(errormsg)
        return snapshot


class BU_iso_count(BU_res_quarantine_count):
    def __init__(self, days, *args, iso_days=ISO_DAYS, **kwargs):
        super().__init__(days, **kwargs) # Initialize the Analyzer object
        self.iso_days = iso_days
    
    
    ''' Add a method to calculate the number of people
        in isolation.'''
    def iso_count(self, ppl, sim_t):
        ''' ppl: sim.people
            sim_t: simulation day'''
        # Avoid numpy errors comparing against nan values
        diagnosed = safe_nan_compare(ppl.date_diagnosed, sim_t, op.le)

        # not diagnosed 
        not_diagnosed = ~diagnosed
        
        recov = ppl.recovered
        not_recov = ~recov
        alive = ~ppl.dead
        not_at_bu = ppl.campResident < 1
        at_bu = ppl.campResident > 0
        not_severe = ~ppl.severe
        not_critical = ~ppl.critical

        # In isolation today based on health
        iso_today_health = not_recov & alive & not_severe & not_critical & diagnosed 
        # In isolation due to iso_days waiting period
        iso_today_waiting = diagnosed & recov & not_severe & not_critical & safe_nan_compare(ppl.date_diagnosed, (sim_t - self.iso_days), np.greater)
        return iso_today_health + iso_today_waiting
 
class BU_res_iso_count(BU_iso_count):
    ''' Snapshot the number of people isolated (diagnosed and not recovered),
        who are resident at BU.  Just override the apply() function '''
    def apply(self, sim):
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
         #   if sim.t == 0:
         #       self.snapshots[date] = 0
          #      return  # 1st day no one can be isolated.
            # if sim didn't start on day 0 back-fill in missing days with 0.
            if sim.t > 0 and len(self.snapshots) == 0:
                for i in range(sim.t):
                    self.snapshots[self.dates[i]] = 0
            # campResident is: 0 off campus, 1 on campus, 2 large dorm on campus
            at_bu = sim.people.campResident > 0
            self.snapshots[date] = np.sum(self.iso_count(sim.people, sim.t) & at_bu)
        return    
 

class BU_nonres_quarantine_count(BU_res_quarantine_count):
    ''' Snapshot the number of people in quarantine who are not 
        BU residents. '''
    def apply(self, sim):
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            # if sim didn't start on day 0 back-fill in missing days with 0.
            if sim.t > 0 and len(self.snapshots) == 0:
                for i in range(sim.t):
                    self.snapshots[self.dates[i]] = 0

            # campResident is: 0 off campus, 1 on campus, 2 large dorm on campus
            #sim.people.quarantined is a Boolean array of everyone who is quarantined
            self.snapshots[date] = np.sum(np.logical_and(sim.people.quarantined,sim.people.campResident < 1)) 
        return    
    
class BU_nonres_iso_count(BU_iso_count):
    ''' Snapshot the number of people isolated (diagnosed and not recovered),
        who are not resident at BU.  '''
    def apply(self, sim):
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
       #     import pdb ; pdb.set_trace()
        #    if sim.t == 0:
      #          self.snapshots[date] = 0
      #          return  # 1st day no one can be isolated.            
            # if sim didn't start on day 0 back-fill in missing days with 0.
            if sim.t > 0 and len(self.snapshots) == 0:
                for i in range(sim.t):
                    self.snapshots[self.dates[i]] = 0     
            # campResident is: 0 off campus, 1 on campus, 2 large dorm on campus
            not_at_bu = sim.people.campResident < 1
            self.snapshots[date] = np.sum(self.iso_count(sim.people, sim.t) & not_at_bu)
        return    


class BU_infection_count(BU_res_quarantine_count):
    ''' Snapshot the demographics of anyone who is 
        infected on any given day '''
    def apply(self,sim):
        ppl = sim.people
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            
            # They're infected today if their exposure date equals today's date.
            #today_diag = np.where(ppl.date_diagnosed.astype(np.int32) == sim.t)
            today_diag = np.where(ppl.date_exposed.astype(np.int32) == sim.t)
            # This stores several quantities on each date:
            # list of ages of infected
            # list of group (student/faculty/etc)
            # etc. Only store them if there is something found.
            if len(today_diag[0]) > 0:
                self.snapshots[date] = {}
                self.snapshots[date]['age'] = ppl.age[today_diag]
                self.snapshots[date]['test_cat'] = ppl.test_cat[today_diag]
                self.snapshots[date]['campResident'] = ppl.campResident[today_diag]
                self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_diag]
                self.snapshots[date]['category'] = ppl.category[today_diag]
                self.snapshots[date]['campus'] = ppl.campus[today_diag]
                self.snapshots[date]['undergrad'] = ppl.undergrad[today_diag]
                self.snapshots[date]['exogenous'] = np.ones(len(today_diag[0]),dtype=np.uint8)
                self.snapshots[date]['GI'] = np.zeros(len(today_diag[0]),dtype=np.uint8)
                self.snapshots[date]['source'] = np.zeros(len(today_diag[0]),dtype=np.uint32)
                source = [item['source'] for item in ppl.infection_log if item['target'] in today_diag[0]]
                for ind, val in enumerate(source):
                    if val is not None:
                        self.snapshots[date]['GI'][ind] = sim.t - [item['date'] for item in ppl.infection_log if item['target'] == val][0]
                        self.snapshots[date]['exogenous'][ind] = 0 
                        self.snapshots[date]['source'][ind] = ppl.full_info_id[val]
                     
                     
class BU_diag2iso_count(BU_res_quarantine_count):
    ''' Count every time someone in quarantine is diagnosed '''
    def apply(self,sim):
        ppl = sim.people
        if sim.t > 0:
            for ind in cv.interventions.find_day(self.days, sim.t):
                date = self.dates[ind]
                # Find everyone who is diagnosed today who was in quarantine yesteday
                today_diag = np.where(ppl.date_diagnosed.astype(np.int32) == sim.t)
                if len(today_diag[0]) > 0:
                    # temporary so that the boolean logic works.
                    tmp = np.full(self.yesterday_quarantine.shape,False)
                    tmp[today_diag] = True
                    delta_diag = tmp & self.yesterday_quarantine
                    
                    # Consider only those in BU housing
                    res_delta_diag = delta_diag #& (sim.people.campResident > 0)
                    # Get their indices, we'll store demographic info about them
                    res_diag_inds = np.where(res_delta_diag)
                    if len(res_diag_inds[0]) > 0:
                        self.snapshots[date] = {}
                        self.snapshots[date]['age'] = ppl.age[res_diag_inds]
                        self.snapshots[date]['test_cat'] = ppl.test_cat[res_diag_inds]
                        self.snapshots[date]['campResident'] = ppl.campResident[res_diag_inds]
                        self.snapshots[date]['full_info_id'] = ppl.full_info_id[res_diag_inds]
                        self.snapshots[date]['category'] = ppl.category[res_diag_inds]
        # Store today's quarantine numbers for next time
        self.yesterday_quarantine = sim.people.quarantined.copy()
            
class BU_severe_count(BU_res_quarantine_count):
    ''' Snapshot the demographics of anyone who is 
        infected on any given day '''
    def apply(self,sim):
        ppl = sim.people
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            # They're diagnosed today if their diagnosis date equals today's date.
            #today_diag = np.where(ppl.date_diagnosed.astype(np.int32) == sim.t)
            today_diag = np.where(ppl.date_severe.astype(np.int32) == sim.t)
            # This stores several quantities on each date:
            # list of ages of infected
            # list of group (student/faculty/etc)
            # etc. Only store them if there is something found.
            if len(today_diag[0]) > 0:
                self.snapshots[date] = {}
                self.snapshots[date]['age'] = ppl.age[today_diag]
                self.snapshots[date]['test_cat'] = ppl.test_cat[today_diag]
                self.snapshots[date]['campResident'] = ppl.campResident[today_diag]
                self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_diag]
                self.snapshots[date]['category'] = ppl.category[today_diag]

class BU_critical_count(BU_res_quarantine_count):
    ''' Snapshot the demographics of anyone who is 
        infected on any given day '''
    def apply(self,sim):
        ppl = sim.people
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            # They're diagnosed today if their diagnosis date equals today's date.
            #today_diag = np.where(ppl.date_diagnosed.astype(np.int32) == sim.t)
            today_diag = np.where(ppl.date_critical.astype(np.int32) == sim.t)
            # This stores several quantities on each date:
            # list of ages of infected
            # list of group (student/faculty/etc)
            # etc. Only store them if there is something found.
            if len(today_diag[0]) > 0:
                self.snapshots[date] = {}
                self.snapshots[date]['age'] = ppl.age[today_diag]
                self.snapshots[date]['test_cat'] = ppl.test_cat[today_diag]
                self.snapshots[date]['campResident'] = ppl.campResident[today_diag]
                self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_diag]
                self.snapshots[date]['category'] = ppl.category[today_diag]

class BU_dead_count(BU_res_quarantine_count):
    ''' Snapshot the demographics of anyone who is 
        infected on any given day '''
    def apply(self,sim):
        ppl = sim.people
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            # They're diagnosed today if their diagnosis date equals today's date.
            #today_diag = np.where(ppl.date_diagnosed.astype(np.int32) == sim.t)
            today_diag = np.where(ppl.date_dead.astype(np.int32) == sim.t)
            # This stores several quantities on each date:
            # list of ages of infected
            # list of group (student/faculty/etc)
            # etc. Only store them if there is something found.
            if len(today_diag[0]) > 0:
                self.snapshots[date] = {}
                self.snapshots[date]['age'] = ppl.age[today_diag]
                self.snapshots[date]['test_cat'] = ppl.test_cat[today_diag]
                self.snapshots[date]['campResident'] = ppl.campResident[today_diag]
                self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_diag]
                self.snapshots[date]['category'] = ppl.category[today_diag]

class BU_diagnosed_count(BU_res_quarantine_count):
    ''' Snapshot the demographics of anyone who is 
        infected on any given day '''
    def apply(self,sim):
        ppl = sim.people
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            # They're diagnosed today if their diagnosis date equals today's date.
            #today_diag = np.where(ppl.date_diagnosed.astype(np.int32) == sim.t)
            today_diag = cv.true(ppl.date_diagnosed.astype(np.int32) == sim.t)
            # This stores several quantities on each date:
            # list of ages of infected
            # list of group (student/faculty/etc)
            # etc. Only store them if there is something found.
            if  today_diag.size > 0:
                self.snapshots[date] = {}
                self.snapshots[date]['age'] = ppl.age[today_diag]
                self.snapshots[date]['test_cat'] = ppl.test_cat[today_diag]
                self.snapshots[date]['campResident'] = ppl.campResident[today_diag]
                self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_diag]
                self.snapshots[date]['category'] = ppl.category[today_diag]
                self.snapshots[date]['campus'] = ppl.campus[today_diag]
                self.snapshots[date]['undergrad'] = ppl.undergrad[today_diag]
                self.snapshots[date]['exogenous'] = np.ones(today_diag.size,dtype=np.uint8)
                self.snapshots[date]['source'] = np.zeros(today_diag.size,dtype=np.uint32)
                source = [item['source'] for item in ppl.infection_log if item['target'] in today_diag]
                for ind, val in enumerate(source):
                    if val is not None:
                        self.snapshots[date]['exogenous'][ind] = 0 
                        self.snapshots[date]['source'][ind] = ppl.full_info_id[val]

class BU_recovered_count(BU_res_quarantine_count):
    ''' Snapshot the demographics of anyone who is 
        recovered on any given day '''
    def apply(self,sim):
        ppl = sim.people
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            # They're diagnosed today if their diagnosis date equals today's date.
            #today_diag = np.where(ppl.date_diagnosed.astype(np.int32) == sim.t)
            today_diag = np.where(ppl.date_recovered.astype(np.int32) == sim.t)
            # This stores several quantities on each date:
            # list of ages of infected
            # list of group (student/faculty/etc)
            # etc. Only store them if there is something found.
            if len(today_diag[0]) > 0:
                self.snapshots[date] = {}
                self.snapshots[date]['age'] = ppl.age[today_diag]
                self.snapshots[date]['test_cat'] = ppl.test_cat[today_diag]
                self.snapshots[date]['campResident'] = ppl.campResident[today_diag]
                self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_diag]
                self.snapshots[date]['category'] = ppl.category[today_diag]
                self.snapshots[date]['asymptomatic'] = np.zeros(len(today_diag[0]),dtype=np.uint8)
                self.snapshots[date]['isolation_period'] = np.zeros(len(today_diag[0]),dtype=np.uint8)
                date_symptomatic = ppl.date_symptomatic[today_diag]
                date_diagnosed = ppl.date_diagnosed[today_diag]
                for ind, val in enumerate(date_symptomatic):
                    if np.isnan(val):
                        self.snapshots[date]['asymptomatic'][ind] = 1
                for ind, val in enumerate(date_diagnosed):
                    if not np.isnan(val):
                        self.snapshots[date]['isolation_period'][ind] = sim.t-val                        
                        
class BU_quarantined_count(BU_res_quarantine_count):
    ''' Snapshot the demographics of anyone who is 
        quarantined on any given day '''
    def apply(self,sim):
        ppl = sim.people
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            # They're diagnosed today if their diagnosis date equals today's date.
            #today_diag = np.where(ppl.date_diagnosed.astype(np.int32) == sim.t)
            today_quar = np.where(ppl.date_quarantined.astype(np.int32) == sim.t)
            # This stores several quantities on each date:
            # list of ages of infected
            # list of group (student/faculty/etc)
            # etc. Only store them if there is something found.
            if len(today_quar[0]) > 0:
                self.snapshots[date] = {}
                self.snapshots[date]['age'] = ppl.age[today_quar]
                self.snapshots[date]['test_cat'] = ppl.test_cat[today_quar]
                self.snapshots[date]['campResident'] = ppl.campResident[today_quar]
                self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_quar]
                self.snapshots[date]['category'] = ppl.category[today_quar]
 

class BU_quarantined_end_count(BU_res_quarantine_count):
    ''' Count every time someone in quarantine is diagnosed '''
    def apply(self,sim):
        ppl = sim.people
        if sim.t > 0:
            for ind in cv.interventions.find_day(self.days, sim.t):
                date = self.dates[ind]
            # Find everyone who is diagnosed today who was in quarantine yesteday
                not_quarantined = ~ ppl.quarantined
                res_delta_diag = not_quarantined & self.yesterday_quarantine
            # Get their indices, we'll store demographic info about them
                res_diag_inds = np.where(res_delta_diag)
                if len(res_diag_inds[0]) > 0:
                    self.snapshots[date] = {}
                    self.snapshots[date]['age'] = ppl.age[res_diag_inds]
                    self.snapshots[date]['test_cat'] = ppl.test_cat[res_diag_inds]
                    self.snapshots[date]['campResident'] = ppl.campResident[res_diag_inds]
                    self.snapshots[date]['full_info_id'] = ppl.full_info_id[res_diag_inds]
                    self.snapshots[date]['category'] = ppl.category[res_diag_inds]
    # Store today's quarantine numbers for next time
        self.yesterday_quarantine = sim.people.quarantined.copy()

 
class BU_cleaning_rooms_count(BU_iso_count):
    ''' Count the number of isolations rooms in use, including ones
    undergoing cleaning for a specified number of days. This one will record a 4-element
    dictionary for each day - watch for that in the dataframe conversion.
    
    The daily snapshot dictionary keys are: 
        n_oncampus_quarantine_cleaning, n_offcampus_quarantine_cleaning,
        n_oncampus_isolation_cleaning, n_offcampus_isolation_cleaning
    '''
    def __init__(self, days, *args, quar_cleaning_days=QUAR_CLEANING_DAYS, 
                 iso_cleaning_days=ISO_CLEANING_DAYS, iso_days=ISO_DAYS,**kwargs):
        super().__init__(days,**kwargs) # Initialize the BU_res_quarantine_count object
        # Number of days it takes to clean a quarantine room
        self.quar_cleaning_days = quar_cleaning_days 
        self.iso_cleaning_days = iso_cleaning_days
        self.iso_days = ISO_DAYS # Req number of days in BU isolation rooms.
        
        # Number of quar rooms in the cleaning state plus those leaving today.
        self.quar_cleaning = deque(maxlen=quar_cleaning_days) 
        # Number of iso rooms in the cleaning state plus those leaving today
        self.iso_cleaning = deque(maxlen=iso_cleaning_days) 

        # Initialize the deque.
        self.quar_cleaning_keys = ['n_oncampus_quarantine_cleaning',
                                   'n_offcampus_quarantine_cleaning']
        self.iso_cleaning_keys = ['n_oncampus_isolation_cleaning', 
                                 'n_offcampus_isolation_cleaning']
        self.all_keys =  self.quar_cleaning_keys + self.iso_cleaning_keys + ['n_leaving_on_iso_today']
        
        for i in range(len(self.quar_cleaning_keys)):
            self.quar_cleaning.append(dict(zip(self.quar_cleaning_keys ,it.repeat(0))))
        for i in range(len(self.iso_cleaning_keys)):
            self.iso_cleaning.append(dict(zip(self.iso_cleaning_keys ,it.repeat(0))))
        return

    def cleaning_sum(self, deque, keys):
        ''' Sum up the  a deque by dictionary keys'''
        count_sum = dict(zip(keys, it.repeat(0)))
        for c in deque: # get each dict
            for k in keys: # for each key add it
                count_sum[k] += c[k]
        return count_sum

    def apply(self,sim):
        ind = cv.interventions.find_day(self.days, sim.t)
        if len(ind) < 1:
            return
        ind = ind[0]
        date = self.dates[ind]
        # On day 1 do nothing.
        if not hasattr(self, 'yesterday_date_end_quarantine') or not hasattr(self, 'yesterday_quarantine'):
            quar_sum = self.cleaning_sum(self.quar_cleaning, self.quar_cleaning_keys)
            iso_sum = self.cleaning_sum(self.iso_cleaning, self.iso_cleaning_keys)
            self.snapshots[date] = {**quar_sum, **iso_sum}
            self.snapshots[date]['n_leaving_on_iso_today'] = 0     
            self.yesterday_date_end_quarantine = sim.people.date_end_quarantine.copy()
            self.yesterday_quarantine = sim.people.quarantined.copy()       
            # Now, in the case where the initial timestep starts at a time > 0, backfill
            # the dates into the snapshot so the output works  correctly.
            if sim.t > 0:
                for i in range(0,sim.t):
                    tmp = sim.date(i)
                    self.snapshots[tmp] = {**quar_sum, **iso_sum}
                    self.snapshots[tmp]['n_leaving_on_iso_today'] = 0   
            return
        # Now the simulation is running...carry on.
        ppl = sim.people
        diagnosed = safe_nan_compare(ppl.date_diagnosed, sim.t, op.le)
        not_diagnosed = ~diagnosed
        # Is this the day they are released?
        quar_freedom = self.yesterday_date_end_quarantine == sim.t
        # not diagnosed 
        not_diagnosed = ~diagnosed
        # Quarantined yesterday?
        was_quar = self.yesterday_quarantine
        alive = ~sim.people.dead
        not_at_bu = sim.people.campResident < 1
        at_bu = sim.people.campResident > 0
         
        # Find all people who are leaving quarantine today who are not 
        # recovered (i.e. never got sick - quarantine only) and count them.
        # Do this twice for on and off campus residents.  Don't AND against
        # ppl.quarantined because that will be false by the time this is called!
        on_quar_leaving = np.sum(quar_freedom & was_quar & not_diagnosed & at_bu & alive)
        off_quar_leaving = np.sum(quar_freedom & was_quar & not_diagnosed & not_at_bu & alive )

        # Now for leaving isolation:  isolation means that they are diagnosed, NOT IN QUARANTINE, recovered, 
        # were not critical, not severe
        # remember: critical --> ICU, severe --> hospital
        # So if they're leaving isolation use recovered==True as the date_end_quarantine gets incremented daily 
        # until they've recovered in which case it jumps forward by 14.
        
        # Those leaving isolation are those who were in isolation yesterday but today have
        # finished 10 days post-diagnosis, gone to the hospital or have died. 
        # However if they have gone to severe or critical symptoms they go immediately to the hospital
        # and leave isolation that day. 
        # 3 cases. OR them all together
        iso_leaving_1 = safe_nan_compare(ppl.date_diagnosed, sim.t - self.iso_days, op.eq) & (safe_nan_compare(ppl.date_recovered, sim.t, op.le)) 
        iso_leaving_2 = safe_nan_compare(ppl.date_diagnosed, sim.t - self.iso_days, op.lt) & (safe_nan_compare(ppl.date_recovered, sim.t, op.eq)) 
        iso_leaving_3 = safe_nan_compare(ppl.date_severe, sim.t, op.eq) | \
                        safe_nan_compare(ppl.date_critical, sim.t, op.eq) |  \
                        safe_nan_compare(ppl.date_dead, sim.t, op.eq)
        iso_leaving = iso_leaving_1 | iso_leaving_2 | iso_leaving_3
        
        on_iso_leaving = np.sum(iso_leaving & at_bu)
        off_iso_leaving = np.sum(iso_leaving & not_at_bu)

        # Finally store today's leaving count to use tomorrow.
        self.quar_cleaning.append(dict(zip(self.quar_cleaning_keys,[on_quar_leaving,off_quar_leaving])))
        self.iso_cleaning.append(dict(zip(self.iso_cleaning_keys,[on_iso_leaving,off_iso_leaving])))
                                      
        # Now from the deques get the sums
        quar_sum = self.cleaning_sum(self.quar_cleaning, self.quar_cleaning_keys)
        iso_sum = self.cleaning_sum(self.iso_cleaning, self.iso_cleaning_keys)

        self.snapshots[date] = {**quar_sum, **iso_sum}
        # Also store the number of people leaving isolation today.
        self.snapshots[date]['n_leaving_on_iso_today'] = on_iso_leaving
        # Copy the quarantine  list  to use tomorrow.
        self.yesterday_quarantine = sim.people.quarantined.copy()
        self.yesterday_date_end_quarantine = sim.people.date_end_quarantine.copy()


class BU_pos_test_date_count(BU_res_quarantine_count):
    ''' Snapshot to get the daily count of people who tested positive and their
        on or off campus residency '''
    # specify the keys that get used in many places
    key_on_campus='on-campus'
    key_off_campus='off-campus'
    df_keys = ['sim_num','date_pos_test','on-campus-count','off-campus-count']
    # overriede the apply method like the others
    def apply(self, sim):
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            # if sim didn't start on day 0 back-fill in missing days with 0.
            if sim.t > 0 and len(self.snapshots) == 0:
                for i in range(sim.t):
                    self.snapshots[self.dates[i]] = {self.key_on_campus:0, self.key_off_campus:0}
            # Get the indices of everyone who is diagnosed today.
            diag_ind = np.where(np.logical_and(np.isfinite(sim.people.date_diagnosed),sim.people.date_diagnosed==sim.t))[0]
            # just loop thru those and count on & off campus as we go.  This could be vectorized
            # but there's no real benefit as this list will always be small.
            on_campus = off_campus = 0
            for ind in diag_ind:
                # campResident is: 0 off campus, 1 on campus, 2 large dorm on campus
                if sim.people.campResident[ind] > 0:
                    on_campus += 1 
                else:
                    off_campus += 1
            self.snapshots[date] = {self.key_on_campus:on_campus, self.key_off_campus:off_campus}
        return
    
    def to_dict(self, sim_num=0):
        ''' make a dictionary of lists from this set of data.
            sim_num - simulation number for parallel runs
        '''
        data = {}
        for key in self.df_keys:
            data[key] = list() 
        dates = sorted(list(self.snapshots.keys()))
        for date in dates:
            data['sim_num'].append(sim_num)
            data['date_pos_test'].append(date)
            data['on-campus-count'].append(self.snapshots[date][self.key_on_campus])
            data['off-campus-count'].append(self.snapshots[date][self.key_off_campus])
        return data

    @staticmethod
    def to_df(sims_complete):
        ''' From a list of completed sims create a dataframe for this object'''
        data = {}
        for key in BU_pos_test_date_count.df_keys:
            data[key] = list() 
        # For each sim, find the snapshot that is this class type.  Get
        # its data dictionary and add it to the local one.  Finally return
        # a pandas dataframe.
        for num,sim in enumerate(sims_complete):
            for snap in sim['analyzers']:
                if isinstance(snap,BU_pos_test_date_count):
                    sim_data = snap.to_dict(num)
                    try:
                        for key in data:
                            data[key] += sim_data[key]
                    except:
                        zzzz=1
        return pd.DataFrame(data=data)




class BU_quarantine_network_count(BU_res_quarantine_count):
    '''  Wenrui's description of what this is capturing:
        For each trial, each day, we will have 6 numbers. 
             (1) # who enter quarantine on campus who only had on campus index cases  
             (2) # who enter quarantine on campus who only had off campus index cases  
             (3) # who enter quarantine on campus who had both off campus and on campus index cases 
             (4) # who enter quarantine off campus who only had on campus index cases  
             (5) # who enter quarantine off campus who only had off campus index cases  
             (6) # who enter quarantine off campus who had both off campus and on campus index cases (edited) 
        
        This is based on contact tracing info as found in the contact_tracing_sens_spec intervention class, 
        and not on infections which is tracked in the sim.people.infection_log.
        
        The contact_tracing_sens_spec_log class has been modified to create a contact_log which records all daily
        contact tracing. This snapshot will examine that on each day and lookup the residency info in the sim.people
        object to gather its data. 
        
        This is sort of a combination of the BU_pos_test_date_count and BU_quarantined_count snapshots.  The apply 
        method here will look a lot like those two mushed together. '''

    # columns in the output dataframe
    keys = ['sim_num','date_quarantined','q_on_ind_on','q_on_ind_off','q_on_ind_on_off','q_off_ind_on','q_off_ind_off','q_off_ind_on_off']

    def __init__(self, interv_contact_log, days,  **kwargs):
        ''' interv_contact_log is an object of type contact_tracing_sens_spec_log that's part of the
            set of interventions in a simulation.'''
        super().__init__(days, **kwargs) # Initialize the Analyzer object
        self.disable = False
        if not isinstance(interv_contact_log, contact_tracing_sens_spec_log):
            print('class BU_quarantine_network_count: interv_contact_log is not an instance of class contact_tracing_sens_spec_log')
            self.disable = True
        self.interv_contact_log = interv_contact_log
        
        
        
        
    def invert_contact_log(self,contact_log):
        ''' invert the contact log so we can see the contact from each person 
        who got them landed in quarantine. Do as a dict by day then by person
        contacted'''
        lut = {}
        for day in contact_log:
            if len(contact_log[day]) > 0:
                day_dt = {}
                # if there's anything here...
                for source in contact_log[day]:
                    t1 = contact_log[day][source]['contacts']  
                    targets=[]
                    for t in t1:
                        targets += list(t)
                    for target in targets:
                        if target not in day_dt:
                            day_dt[target] = []
                        day_dt[target].append(source)
                if len(day_dt) > 0:
                    lut[day] = day_dt
        return lut
        
    def apply(self, sim):
        if self.disable: # wasn't initialized corretly. 
            return
        ppl = sim.people
        # initialize storage for counters.
        data = {}
        for k in self.keys:
            data[k] = 0
         # if sim didn't start on day 0 back-fill in missing days with 0.
        if sim.t > 0 and len(self.snapshots) == 0:
            for i in range(sim.t):
                self.snapshots[self.dates[i]] = data.copy()           
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            # Get the indices of everyone who is going into quarantine today.
            today_quar = np.where(ppl.date_quarantined.astype(np.int32) == sim.t)[0]
            
            if len(today_quar) > 0:
                # Transform the intervention contact log.  Right now on this day it's
                # the list of people contacted in response to someone's positive test. That
                # needs to be inverted.
                contact_lut = self.invert_contact_log(self.interv_contact_log.contact_log)
                contact_days = sorted(list(contact_lut.keys()), reverse=True)
                # For each person quarantined today, search backwards in the contact_lut
                # to find the sources that put them into quarantine.
                targ_src = {}
                for targ in today_quar:
                    if targ in contact_lut[sim.t]:
                        targ_src[targ] = contact_lut[sim.t][targ] # sources that made this person get quarantined
                # Now build up the data for today.
                for targ in targ_src:
                    if ppl.campResident[targ_src[targ]].size > 0:
                        # on campus target and there exist sources
                        if ppl.campResident[targ] > 0: 
                            # Sources all on campus?
                            if np.all(ppl.campResident[targ_src[targ]] > 0):
                                data['q_on_ind_on'] += 1
                            elif np.all(ppl.campResident[targ_src[targ]] < 1):
                                # Sources all off campus
                                data['q_on_ind_off'] += 1
                            else:
                                # it's a mix!
                                data['q_on_ind_on_off'] += 1
                        else: # off campus targer
                            # Sources all on campus?
                            if np.all(ppl.campResident[targ_src[targ]] > 0):
                                data['q_off_ind_on'] += 1
                            elif np.all(ppl.campResident[targ_src[targ]] < 1):
                                # Sources all off campus
                                data['q_off_ind_off'] += 1
                            else:
                                # it's a mix!
                                data['q_off_ind_on_off'] += 1                        
    
            # Finally...store this day's snapshot
            self.snapshots[date] = data
        return
    
    def to_dict(self, sim_num=0):
        ''' make a dictionary of lists from this set of data.
            sim_num - simulation number for parallel runs
        '''
        data = {}
        for key in self.keys:
            data[key] = list() 
        dates = sorted(list(self.snapshots.keys()))
        for date in dates:
            data['sim_num'].append(sim_num)
            data['date_quarantined'].append(date)
            data['q_on_ind_on'].append(self.snapshots[date]['q_on_ind_on'])
            data['q_on_ind_off'].append(self.snapshots[date]['q_on_ind_off'])
            data['q_on_ind_on_off'].append(self.snapshots[date]['q_on_ind_on_off'])
            data['q_off_ind_on'].append(self.snapshots[date]['q_off_ind_on'])
            data['q_off_ind_off'].append(self.snapshots[date]['q_off_ind_off'])
            data['q_off_ind_on_off'].append(self.snapshots[date]['q_off_ind_on_off'])
        return data

    @staticmethod
    def to_df(sims_complete):
        ''' From a list of completed sims create a dataframe for this object'''
        data = {}
        for key in BU_quarantine_network_count.keys:
            data[key] = list() 
        # For each sim, find the snapshot that is this class type.  Get
        # its data dictionary and add it to the local one.  Finally return
        # a pandas dataframe.
        for num,sim in enumerate(sims_complete):
            for snap in sim['analyzers']:
                if isinstance(snap,BU_quarantine_network_count):
                    sim_data = snap.to_dict(num)
                    for key in data:
                        data[key] += sim_data[key]
        return pd.DataFrame(data=data)
    


def snapshots_to_df(sims_complete):
    ''' Take a list of completed simulations with analyzers in the 
        order:  BU_res_quarantine_count, BU_res_diag_count, BU_nonres_quarantine_count, BU_nonres_diag_count.
    
        Return a pandas dataframe:
           sim_num  date  day n_quar n_iso
           
        The sim_num column is the index of the simulations.'''
    data={'sim_num':[], 'dates':[], 'days':[], 'n_res_quar':[], 'n_res_iso':[],
          'n_nonres_quar':[], 'n_nonres_iso':[], 
          'n_oncampus_quarantine_cleaning':[],
          'n_offcampus_quarantine_cleaning':[],
          'n_oncampus_isolation_cleaning':[], 
          'n_offcampus_isolation_cleaning':[],
          'n_leaving_on_iso_today':[]}
    for i, sim in enumerate(sims_complete):    
        # Get the snapshots
        BU_quar = sim['analyzers'][0]
        BU_iso = sim['analyzers'][1]
        non_BU_quar = sim['analyzers'][2]
        non_BU_iso = sim['analyzers'][3]
        BU_quar_rooms_count = sim['analyzers'][13]
        # Get the dates from bu_quar
        data['dates'] += BU_quar.dates
        sim_days = list( BU_quar.days )
        data['days'] += sim_days
        # Extract the data from all snapshots
        data['n_res_quar'] += [BU_quar.snapshots[x] for x in BU_quar.dates]
        data['n_res_iso'] += [BU_iso.snapshots[x] for x in BU_quar.dates]
        data['n_nonres_quar'] += [non_BU_quar.snapshots[x] for x in BU_quar.dates]
        data['n_nonres_iso'] += [non_BU_iso.snapshots[x] for x in BU_quar.dates]        
        # Now..the count of rooms in cleaning is a little different because it snapshots a 4
        # element dictionary.  Break out the dictionaries into columns here.
        cleaning_keys = BU_quar_rooms_count.all_keys
        for key in cleaning_keys:
            data[key] += [BU_quar_rooms_count.snapshots[x][key] for x in BU_quar.dates]        
        # Now fill in the sim_num
        data['sim_num'] += len(sim_days) * [i]

    return pd.DataFrame(data=data)

#%%
def infection_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'test_cat':[],
          'campResident':[], 'full_info_id':[], 'category':[],'campus':[],'undergrad':[],'exogenous':[],'GI':[],'source':[]}
    for i, sim in enumerate(sims_complete):    
        # Get the snapshot
        BU_infect = sim['analyzers'][4]
        # Loop through the dates and add everything found to data.
        # Not the most elegant code, but it gets the job done.
        count = 0 
        for j,date in enumerate(BU_infect.dates):
            days = BU_infect.days[j]
            # Everything that happened on this day gets saved.
            if date in BU_infect.snapshots:
                for k in range(BU_infect.snapshots[date]['age'].shape[0]):
                    data['dates'].append(date)
                    data['days'].append(days)
                    data['age'].append(BU_infect.snapshots[date]['age'][k])
                    data['test_cat'].append(BU_infect.snapshots[date]['test_cat'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    data['campus'].append(BU_infect.snapshots[date]['campus'][k])
                    data['undergrad'].append(BU_infect.snapshots[date]['undergrad'][k])
                    data['exogenous'].append(BU_infect.snapshots[date]['exogenous'][k])
                    data['GI'].append(BU_infect.snapshots[date]['GI'][k])
                    data['source'].append(BU_infect.snapshots[date]['source'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)


#%%
def diag2iso_count_to_df(sims_complete):
    ''' The diag2iso count is the index 5 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'test_cat':[],
          'campResident':[], 'full_info_id':[], 'category':[]}
    for i, sim in enumerate(sims_complete):    
        # Get the snapshot
        BU_infect = sim['analyzers'][5]
        # Loop through the dates and add everything found to data.
        # Not the most elegant code, but it gets the job done.
        count = 0 
        for j,date in enumerate(BU_infect.dates):
            days = BU_infect.days[j]
            # Everything that happened on this day gets saved.
            if date in BU_infect.snapshots:
                for k in range(BU_infect.snapshots[date]['age'].shape[0]):
                    data['dates'].append(date)
                    data['days'].append(days)
                    data['age'].append(BU_infect.snapshots[date]['age'][k])
                    data['test_cat'].append(BU_infect.snapshots[date]['test_cat'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)

def severe_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'test_cat':[],
          'campResident':[], 'full_info_id':[], 'category':[]}
    for i, sim in enumerate(sims_complete):    
        # Get the snapshot
        BU_infect = sim['analyzers'][6]
        # Loop through the dates and add everything found to data.
        # Not the most elegant code, but it gets the job done.
        count = 0 
        for j,date in enumerate(BU_infect.dates):
            days = BU_infect.days[j]
            # Everything that happened on this day gets saved.
            if date in BU_infect.snapshots:
                for k in range(BU_infect.snapshots[date]['age'].shape[0]):
                    data['dates'].append(date)
                    data['days'].append(days)
                    data['age'].append(BU_infect.snapshots[date]['age'][k])
                    data['test_cat'].append(BU_infect.snapshots[date]['test_cat'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)

def critical_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'test_cat':[],
          'campResident':[], 'full_info_id':[], 'category':[]}
    for i, sim in enumerate(sims_complete):    
        # Get the snapshot
        BU_infect = sim['analyzers'][7]
        # Loop through the dates and add everything found to data.
        # Not the most elegant code, but it gets the job done.
        count = 0 
        for j,date in enumerate(BU_infect.dates):
            days = BU_infect.days[j]
            # Everything that happened on this day gets saved.
            if date in BU_infect.snapshots:
                for k in range(BU_infect.snapshots[date]['age'].shape[0]):
                    data['dates'].append(date)
                    data['days'].append(days)
                    data['age'].append(BU_infect.snapshots[date]['age'][k])
                    data['test_cat'].append(BU_infect.snapshots[date]['test_cat'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)

def dead_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'test_cat':[],
          'campResident':[], 'full_info_id':[], 'category':[]}
    for i, sim in enumerate(sims_complete):    
        # Get the snapshot
        BU_infect = sim['analyzers'][8]
        # Loop through the dates and add everything found to data.
        # Not the most elegant code, but it gets the job done.
        count = 0 
        for j,date in enumerate(BU_infect.dates):
            days = BU_infect.days[j]
            # Everything that happened on this day gets saved.
            if date in BU_infect.snapshots:
                for k in range(BU_infect.snapshots[date]['age'].shape[0]):
                    data['dates'].append(date)
                    data['days'].append(days)
                    data['age'].append(BU_infect.snapshots[date]['age'][k])
                    data['test_cat'].append(BU_infect.snapshots[date]['test_cat'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)

def diagnosed_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'test_cat':[],
          'campResident':[], 'full_info_id':[], 'category':[],'campus':[],'undergrad':[],'exogenous':[],'source':[]}    
    for i, sim in enumerate(sims_complete):    
        # Get the snapshot
        BU_infect = sim['analyzers'][9]
        # Loop through the dates and add everything found to data.
        # Not the most elegant code, but it gets the job done.
        count = 0 
        for j,date in enumerate(BU_infect.dates):
            days = BU_infect.days[j]
            # Everything that happened on this day gets saved.
            if date in BU_infect.snapshots:
                for k in range(BU_infect.snapshots[date]['age'].shape[0]):
                    data['dates'].append(date)
                    data['days'].append(days)
                    data['age'].append(BU_infect.snapshots[date]['age'][k])
                    data['test_cat'].append(BU_infect.snapshots[date]['test_cat'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    data['campus'].append(BU_infect.snapshots[date]['campus'][k])
                    data['undergrad'].append(BU_infect.snapshots[date]['undergrad'][k])
                    data['exogenous'].append(BU_infect.snapshots[date]['exogenous'][k])
                    data['source'].append(BU_infect.snapshots[date]['source'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)
    
#%%             
def recovered_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'test_cat':[],
          'campResident':[], 'full_info_id':[], 'category':[],'asymptomatic':[],'isolation_period':[]}
    for i, sim in enumerate(sims_complete):    
        # Get the snapshot
        BU_infect = sim['analyzers'][10]
        # Loop through the dates and add everything found to data.
        # Not the most elegant code, but it gets the job done.
        count = 0 
        for j,date in enumerate(BU_infect.dates):
            days = BU_infect.days[j]
            # Everything that happened on this day gets saved.
            if date in BU_infect.snapshots:
                for k in range(BU_infect.snapshots[date]['age'].shape[0]):
                    data['dates'].append(date)
                    data['days'].append(days)
                    data['age'].append(BU_infect.snapshots[date]['age'][k])
                    data['test_cat'].append(BU_infect.snapshots[date]['test_cat'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    data['asymptomatic'].append(BU_infect.snapshots[date]['asymptomatic'][k])
                    data['isolation_period'].append(BU_infect.snapshots[date]['isolation_period'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)
    
def quarantined_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'test_cat':[],
          'campResident':[], 'full_info_id':[], 'category':[]}
    for i, sim in enumerate(sims_complete):    
        # Get the snapshot
        BU_infect = sim['analyzers'][11]
        # Loop through the dates and add everything found to data.
        # Not the most elegant code, but it gets the job done.
        count = 0 
        for j,date in enumerate(BU_infect.dates):
            days = BU_infect.days[j]
            # Everything that happened on this day gets saved.
            if date in BU_infect.snapshots:
                for k in range(BU_infect.snapshots[date]['age'].shape[0]):
                    data['dates'].append(date)
                    data['days'].append(days)
                    data['age'].append(BU_infect.snapshots[date]['age'][k])
                    data['test_cat'].append(BU_infect.snapshots[date]['test_cat'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)
    
def quarantined_end_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'test_cat':[],
          'campResident':[], 'full_info_id':[], 'category':[]}
    for i, sim in enumerate(sims_complete):    
        # Get the snapshot
        BU_infect = sim['analyzers'][12]
        # Loop through the dates and add everything found to data.
        # Not the most elegant code, but it gets the job done.
        count = 0 
        for j,date in enumerate(BU_infect.dates):
            days = BU_infect.days[j]
            # Everything that happened on this day gets saved.
            if date in BU_infect.snapshots:
                for k in range(BU_infect.snapshots[date]['age'].shape[0]):
                    data['dates'].append(date)
                    data['days'].append(days)
                    data['age'].append(BU_infect.snapshots[date]['age'][k])
                    data['test_cat'].append(BU_infect.snapshots[date]['test_cat'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)    

 
    
def get_BU_snapshots(num_days, quar_cleaning_days=QUAR_CLEANING_DAYS,
                     iso_cleaning_days=ISO_CLEANING_DAYS,
                     iso_days=ISO_DAYS, interv_contact_log = None):
    ''' Return a list of snapshots to be used with the simulations.  The order here
        is specific and must match that in snapshots_to_df 
        
        interv_contact_log is the contact_tracing_sens_spec_log object used in
        the simulation interventions.'''
    day_lst = list(range(num_days))
    return [BU_res_quarantine_count(day_lst),
            BU_res_iso_count(day_lst, iso_days),
            BU_nonres_quarantine_count(day_lst),
            BU_nonres_iso_count(day_lst, iso_days),
            BU_infection_count(day_lst),
            BU_diag2iso_count(day_lst),
            BU_severe_count(day_lst),
            BU_critical_count(day_lst),
            BU_dead_count(day_lst),
            BU_diagnosed_count(day_lst),
            BU_recovered_count(day_lst),
            BU_quarantined_count(day_lst),
            BU_quarantined_end_count(day_lst),
            BU_cleaning_rooms_count(day_lst, quar_cleaning_days,iso_cleaning_days,iso_days),
            BU_pos_test_date_count(day_lst),
            BU_quarantine_network_count(interv_contact_log, day_lst)]
    
