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


__all__=['snapshots_to_df','get_BU_snapshots','infection_count_to_df','diag2iso_count_to_df','severe_count_to_df','critical_count_to_df','dead_count_to_df','diagnosed_count_to_df','recovered_count_to_df','quarantined_count_to_df']

import covasim as cv
import sciris as sc
import numpy as np
import pandas as pd
import copy

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
        super().__init__(**kwargs) # Initialize the Intervention object
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


class BU_res_iso_count(BU_res_quarantine_count):
    ''' Snapshot the number of people isolated (diagnosed and not recovered),
        who are resident at BU.  Just override the apply() function '''
    def apply(self, sim):
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            if sim.t == 0:
                self.snapshots[date] = 0
                return  # 1st day no one can be isolated.
            # campResident is: 0 off campus, 1 on campus, 2 large dorm on campus
            #sim.people.quarantined is a Boolean array of everyone who is quarantined
            # Some shorthand to make the numpy statements easier to read
            not_recov = ~sim.people.recovered
            alive = ~sim.people.dead
            at_bu = sim.people.campResident > 0
            #sick =  sim.people.diagnosed
            # Is anyone diagnosed with returned test results?
            diagnosed_inds = np.where(np.isfinite(sim.people.date_diagnosed))
            diag_today = np.full(sim.people.diagnosed.shape, np.iinfo(np.int32).max)
            diag_today[diagnosed_inds] = sim.people.date_diagnosed[diagnosed_inds] 
            diag_today = diag_today <= sim.t
            not_severe = ~sim.people.severe
            not_critical = ~sim.people.critical
            self.snapshots[date] = np.sum(not_recov & alive & at_bu  & not_severe & not_critical & diag_today)
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
    
class BU_nonres_iso_count(BU_res_quarantine_count):
    ''' Snapshot the number of people isolated (diagnosed and not recovered),
        who are not resident at BU.  '''
    def apply(self, sim):
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            if sim.t == 0:
                self.snapshots[date] = 0
                return  # 1st day no one can be isolated.            # campResident is: 0 off campus, 1 on campus, 2 large dorm on campus
            #sim.people.quarantined is a Boolean array of everyone who is quarantined
            not_recov = ~sim.people.recovered
            alive = ~sim.people.dead
            not_at_bu = sim.people.campResident < 1
            # Is anyone diagnosed with returned test results?
            diagnosed_inds = np.where(np.isfinite(sim.people.date_diagnosed))
            diag_today = np.full(sim.people.diagnosed.shape, np.iinfo(np.int32).max)
            diag_today[diagnosed_inds] = sim.people.date_diagnosed[diagnosed_inds] 
            diag_today = diag_today <= sim.t
            not_severe = ~sim.people.severe
            not_critical = ~sim.people.critical
            self.snapshots[date] = np.sum(not_recov & alive & not_at_bu & not_severe & not_critical & diag_today)
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
                self.snapshots[date]['group'] = ppl.group[today_diag]
                self.snapshots[date]['campResident'] = ppl.campResident[today_diag]
                self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_diag]
                self.snapshots[date]['category'] = ppl.category[today_diag]
                self.snapshots[date]['exogenous'] = np.ones(len(today_diag[0]),dtype=np.uint8)
                self.snapshots[date]['GI'] = np.zeros(len(today_diag[0]),dtype=np.uint8)
                self.snapshots[date]['source'] = np.zeros(len(today_diag[0]),dtype=np.uint8)
                source = [item['source'] for item in ppl.infection_log if item['target'] in today_diag[0]]
                for ind, val in enumerate(source):
                    if val is not None:
                        self.snapshots[date]['GI'][ind] = sim.t - [item['date'] for item in ppl.infection_log if item['target'] == val][0]
                        self.snapshots[date]['exogenous'][ind] = 0 
                        self.snapshots[date]['source'][ind] = val
                        
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
                    res_delta_diag = delta_diag & (sim.people.campResident > 0)
                    # Get their indices, we'll store demographic info about them
                    res_diag_inds = np.where(res_delta_diag)
                    if len(res_diag_inds[0]) > 0:
                        self.snapshots[date] = {}
                        self.snapshots[date]['age'] = ppl.age[today_diag]
                        self.snapshots[date]['group'] = ppl.group[today_diag]
                        self.snapshots[date]['campResident'] = ppl.campResident[today_diag]
                        self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_diag]
                        self.snapshots[date]['category'] = ppl.category[today_diag]
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
                self.snapshots[date]['group'] = ppl.group[today_diag]
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
                self.snapshots[date]['group'] = ppl.group[today_diag]
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
                self.snapshots[date]['group'] = ppl.group[today_diag]
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
            today_diag = np.where(ppl.date_diagnosed.astype(np.int32) == sim.t)
            # This stores several quantities on each date:
            # list of ages of infected
            # list of group (student/faculty/etc)
            # etc. Only store them if there is something found.
            if len(today_diag[0]) > 0:
                self.snapshots[date] = {}
                self.snapshots[date]['age'] = ppl.age[today_diag]
                self.snapshots[date]['group'] = ppl.group[today_diag]
                self.snapshots[date]['campResident'] = ppl.campResident[today_diag]
                self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_diag]
                self.snapshots[date]['category'] = ppl.category[today_diag]
                self.snapshots[date]['exogenous'] = np.ones(len(today_diag[0]),dtype=np.uint8)
                source = [item['source'] for item in ppl.infection_log if item['target'] in today_diag[0]]
                for ind, val in enumerate(source):
                    if val is not None:
                        self.snapshots[date]['exogenous'][ind] = 0 

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
                self.snapshots[date]['group'] = ppl.group[today_diag]
                self.snapshots[date]['campResident'] = ppl.campResident[today_diag]
                self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_diag]
                self.snapshots[date]['category'] = ppl.category[today_diag]
                self.snapshots[date]['asymptomatic'] = np.zeros(len(today_diag[0]),dtype=np.uint8)
                date_symptomatic = ppl.date_symptomatic[today_diag]
                for ind, val in enumerate(date_symptomatic):
                    if np.isnan(val):
                        self.snapshots[date]['asymptomatic'][ind] = 1
                        
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
                self.snapshots[date]['group'] = ppl.group[today_diag]
                self.snapshots[date]['campResident'] = ppl.campResident[today_diag]
                self.snapshots[date]['full_info_id'] = ppl.full_info_id[today_diag]
                self.snapshots[date]['category'] = ppl.category[today_diag]



                        
def snapshots_to_df(sims_complete):
    ''' Take a list of completed simulations with analyzers in the 
        order:  BU_res_quarantine_count, BU_res_diag_count, BU_nonres_quarantine_count, BU_nonres_diag_count.
    
        Return a pandas dataframe:
           sim_num  date  day n_quar n_iso
           
        The sim_num column is the index of the simulations.'''
    data={'sim_num':[], 'dates':[], 'days':[], 'n_res_quar':[], 'n_res_iso':[],
          'n_nonres_quar':[], 'n_nonres_iso':[]}
    for i, sim in enumerate(sims_complete):    
        # Get the snapshots
        BU_quar = sim['analyzers'][0]
        BU_iso = sim['analyzers'][1]
        non_BU_quar = sim['analyzers'][2]
        non_BU_iso = sim['analyzers'][3]
        # Get the dates from bu_quar
        data['dates'] += BU_quar.dates
        sim_days = list( BU_quar.days )
        data['days'] += sim_days
        # Extract the data from both snapshots
        data['n_res_quar'] += [BU_quar.snapshots[x] for x in BU_quar.dates]
        data['n_res_iso'] += [BU_iso.snapshots[x] for x in BU_quar.dates]
        data['n_nonres_quar'] += [non_BU_quar.snapshots[x] for x in BU_quar.dates]
        data['n_nonres_iso'] += [non_BU_iso.snapshots[x] for x in BU_quar.dates]        
        # Now fill in the sim_num
        data['sim_num'] += len(sim_days) * [i]

    return pd.DataFrame(data=data)

#%%
def infection_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'group':[],
          'campResident':[], 'full_info_id':[], 'category':[],'exogenous':[],'GI':[],'source':[]}
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
                    data['group'].append(BU_infect.snapshots[date]['group'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
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
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'group':[],
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
                    data['group'].append(BU_infect.snapshots[date]['group'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)

def severe_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'group':[],
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
                    data['group'].append(BU_infect.snapshots[date]['group'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)

def critical_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'group':[],
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
                    data['group'].append(BU_infect.snapshots[date]['group'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)

def dead_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'group':[],
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
                    data['group'].append(BU_infect.snapshots[date]['group'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)

def diagnosed_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'group':[],
          'campResident':[], 'full_info_id':[], 'category':[],'exogenous':[]}
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
                    data['group'].append(BU_infect.snapshots[date]['group'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    data['exogenous'].append(BU_infect.snapshots[date]['exogenous'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)
#%%             
def recovered_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'group':[],
          'campResident':[], 'full_info_id':[], 'category':[],'asymptomatic':[]}
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
                    data['group'].append(BU_infect.snapshots[date]['group'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    data['asymptomatic'].append(BU_infect.snapshots[date]['asymptomatic'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)
    
def quarantined_count_to_df(sims_complete):
    ''' The infection count is the index 4 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'age':[], 'group':[],
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
                    data['group'].append(BU_infect.snapshots[date]['group'][k])
                    data['campResident'].append(BU_infect.snapshots[date]['campResident'][k])
                    data['full_info_id'].append(BU_infect.snapshots[date]['full_info_id'][k])
                    data['category'].append(BU_infect.snapshots[date]['category'][k])
                    count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)
    
def get_BU_snapshots(num_days):
    ''' Return a list of snapshots to be used with the simulations.  The order here
        is specific and must match that in snapshots_to_df '''
    return [BU_res_quarantine_count(list(range(num_days))),
            BU_res_iso_count(list(range(num_days))),
            BU_nonres_quarantine_count(list(range(num_days))),
            BU_nonres_iso_count(list(range(num_days))),
            BU_infection_count(list(range(num_days))),
            BU_diag2iso_count(list(range(num_days))),
            BU_severe_count(list(range(num_days))),
            BU_critical_count(list(range(num_days))),
            BU_dead_count(list(range(num_days))),
            BU_diagnosed_count(list(range(num_days))),
            BU_recovered_count(list(range(num_days))),
            BU_quarantined_count(list(range(num_days))),]
    
