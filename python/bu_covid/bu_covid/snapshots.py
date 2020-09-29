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
         'ISO_CLEANING_DAYS','ISO_DAYS','safe_nan_compare']

import covasim as cv
import sciris as sc
import numpy as np
import pandas as pd
from collections import deque
import copy
import covasim.utils as cvu
import itertools as it 
import numba
import operator 

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
        diagnosed = safe_nan_compare(ppl.date_diagnosed, sim_t, operator.le)

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
        iso_today_waiting = diagnosed & recov & not_severe & not_critical & safe_nan_compare(ppl.date_diagnosed, (sim_t - self.iso_days), operator.gt)
        return iso_today_health + iso_today_waiting
 
class BU_res_iso_count(BU_iso_count):
    ''' Snapshot the number of people isolated (diagnosed and not recovered),
        who are resident at BU.  Just override the apply() function '''
    def apply(self, sim):
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            if sim.t == 0:
                self.snapshots[date] = 0
                return  # 1st day no one can be isolated.
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
           # import pdb  ; pdb.set_trace()
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
            if sim.t == 0:
                self.snapshots[date] = 0
                return  # 1st day no one can be isolated.            
                
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
            # They're diagnosed today if their diagnosis date equals today's date.
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
        infected on any given day '''
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
        infected on any given day '''
    def apply(self,sim):
        ppl = sim.people
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            # They're diagnosed today if their diagnosis date equals today's date.
            #today_diag = np.where(ppl.date_diagnosed.astype(np.int32) == sim.t)
            today_diag = np.where(ppl.date_quarantined.astype(np.int32) == sim.t)
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

#class BU_quarantined_end_count(BU_res_quarantine_count):
#    ''' Snapshot the demographics of anyone who is 
#        infected on any given day '''
#    def apply(self,sim):
#        ppl = sim.people
#        for ind in cv.interventions.find_day(self.days, sim.t):
#            date = self.dates[ind]
#            # They're diagnosed today if their diagnosis date equals today's date.
#            today_diag = cvu.true(ppl.date_end_quarantine==sim.t)
#            #today_diag = np.where(ppl.date_end_quarantine.astype(np.int32) == sim.t)
#            # This stores several quantities on each date:
#            # list of ages of infected
#            # list of group (student/faculty/etc)
#            # etc. Only store them if there is something found.
#            if len(today_diag) > 0:
#                self.snapshots[date] = {}
#                self.snapshots[date]['age'] = ppl.age[today_diag]
#                self.snapshots[date]['test_cat'] = ppl.test_cat[today_diag]
#                self.snapshots[date]['campResident'] = ppl.campResident[today_diag]
#                self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_diag]
#                self.snapshots[date]['category'] = ppl.category[today_diag]

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


# For the output, I think we could add 4 extra columns in the snapshots_int*.csv: 
#     number of quarantine room unavailable due to cleaning (on-campus), number of 
#     quarantine room unavailable due to cleaning (off-campus), number of isolation
#     room unavailable due to cleaning (on-campus), number of isolation room 
#     unavailable due to cleaning (off-campus).
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
        if sim.t == 0:
            quar_sum = self.cleaning_sum(self.quar_cleaning, self.quar_cleaning_keys)
            iso_sum = self.cleaning_sum(self.iso_cleaning, self.iso_cleaning_keys)
            self.snapshots[date] = {**quar_sum, **iso_sum}
            self.snapshots[date]['n_leaving_on_iso_today'] = 0     
            self.yesterday_date_end_quarantine = sim.people.date_end_quarantine.copy()
            self.yesterday_quarantine = sim.people.quarantined.copy()       
            return
        # Now the simulation is running...carry on.
        ppl = sim.people
        diagnosed = safe_nan_compare(ppl.date_diagnosed, sim.t, operator.le)
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
        iso_leaving_1 = safe_nan_compare(ppl.date_diagnosed, sim.t - self.iso_days, operator.eq) & (safe_nan_compare(ppl.date_recovered, sim.t, operator.le)) 
        iso_leaving_2 = safe_nan_compare(ppl.date_diagnosed, sim.t - self.iso_days, operator.lt) & (safe_nan_compare(ppl.date_recovered, sim.t, operator.eq)) 
        iso_leaving_3 = safe_nan_compare(ppl.date_severe, sim.t, operator.eq) | \
                        safe_nan_compare(ppl.date_critical, sim.t, operator.eq) |  \
                        safe_nan_compare(ppl.date_dead, sim.t, operator.eq)
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
                     iso_days=ISO_DAYS):
    ''' Return a list of snapshots to be used with the simulations.  The order here
        is specific and must match that in snapshots_to_df '''
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
            BU_cleaning_rooms_count(day_lst, quar_cleaning_days,iso_cleaning_days,iso_days)]
    
