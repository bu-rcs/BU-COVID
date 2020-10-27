#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 11:16:58 2020

@author: bgregor

Add code to support modifying a Sim object before it's run to reflect the
current state of the community.
"""

import numpy as np
from .misc import make_dt
import covasim.utils as cvu

__all__ = ['import_community_test_info']


def parse_test_info_file(info_file, start_day):
    ''' Read the CSV pop_test_info.csv into a dict of lists.
        state_file: CSV file of the current state of the community
        start_day: start date string for the entire simulation. 
    '''
    # header: id,TestCategory,LabTestEveryNHours,PositiveInd,PositiveCollectDate,LastTestCollectDate,LastTestResultDate,DiagnoseDate,RecoverDate,EarliestSymptomaticDate,LatestSymptomaticDate,EndQuarantineDate,DateKnownContact,DateOfDeath
    start_dt = make_dt(start_day)

    data = {}
    with open(info_file) as f:
        # Get the headers
        headers = [x.strip() for x in f.readline().split(',')]
        # Initialize the storage
        for h in headers:
            data[h] = []
        # read each line, 
        for i,line in enumerate(f):
            line = line.split(',')
            try:
                col = 0
                for h in headers:
                    data[h].append(line[col].strip())
                    col += 1
            except:
                print('Error reading %s at line %s. Incorrect number of columns.' % (info_file,i))
    # For the date string columns, convert to date numbers.
    # for non-date string columns convert to int.
    # The PositiveInd column is blank, N, or Y. Convert to -1,0,1

    for i,val in enumerate(data['PositiveInd']):
        if val.upper() == 'N':
            data['PositiveInd'][i] = 0
        elif val.upper() == 'Y':
            data['PositiveInd'][i] = 1
        else:
            data['PositiveInd'][i] = -1
    keys = list(data.keys())
    keys.remove('PositiveInd') 
    # If a key contains the word Date it's a date, otherwise an int.
    for key in keys:
        if key.lower().find('date') >= 0:
            for i,val in enumerate(data[key]):
                if val != '':
                    data[key][i] = (make_dt(data[key][i]) - start_dt).days  
                else:
                    data[key][i] = None
        else:
            for i,val in enumerate(data[key]):
                if val.isdigit():
                    data[key][i] = int(data[key][i])
                else:
                    data[key][i] = None
    # Refactor into a dictionary of id's with an inner dictionary of values.
    headers.remove('id')
    id_data={}
    for i,key in enumerate(data['id']):
        id_data[key] = {}
        for h in headers:
            id_data[key][h] = data[h][i]
    return id_data

def import_community_test_info(sim, info_file):
    ''' sim: a Simulation object that has not been run.
        info_file: a CSV file that describes the current state of
                    each individual in the Sim.People object. 
    '''
    id_data = parse_test_info_file(info_file, sim.date(0))
    # id_data is a dictionary of people id's.  The matching field is
    # sim.people.full_info_id
    # Loop through sim.people.full_info_id. For each id, look it up
    # in id_data.  Fill in any needed fields in a big set of 'if' statements.
    
    # Here's the association of pop_test_info.csv columns to 
    # sim.people columns.  Fill in the date and test status.
    # TestCategory              category
    # LabTestEveryNHours	        labTestEveryNHours
    # PositiveCollectDate	    date_pos_test
    # LastTestCollectDate	    date_tested
    # LastTestResultDate      	none
    # DiagnoseDate	            date_diagnosed
    # RecoverDate	           date_recovered
    # EarliestSymptomaticDate	date_symptomatic
    # LatestSymptomaticDate	date_symptomatic
    # EndQuarantineDate	date_end_quarantine
    # DateKnownContact	date_known_contact
    # DateOfDeath	date_dead
    
    # These columns in sim.people will need to be calculated 
    # to make sure they're consistent with the dates
    # quarantined, recov, date_exposed
    count = 0
    count2 = 0 
    for i,info_id in enumerate(sim.people.full_info_id):
        if info_id not in id_data:
            print('ERROR: full_info_id %s not found in %s' % (info_id,info_file))
            continue
        person = id_data[info_id]
        
        # Guard each with an if statement - missing dates are None values in 
        # the id_data dictionaries
        if person['PositiveCollectDate'] != None:
            # If there's a positive collection date before the start of the simulation 
            # set it to nan so it will be ignored.
           # if person['PositiveCollectDate'] < 0:
           #     sim.people.date_pos_test[i] = np.nan 
            # Check when they were positive. If more than 14 days
            # before the start date then mark as recovered.
            # This is copied from people.infect(). This person needs to be offically
            # infected. We'll use asymptomatic for the import.
            dur_asym2rec = cvu.sample(**sim.pars['dur']['asym2rec'], size=1)
            date_recovered = person['PositiveCollectDate'] + dur_asym2rec 
            # If the positive collection date is too far past OR their recovery 
            # date ends up being before the simulation has started, set them as
            # recovered.
            if person['PositiveCollectDate'] <= -sim.pars['quar_period'] or date_recovered < 0:
                sim.people.recovered[i]   = True
                sim.people.susceptible[i] = False
                sim.people.exposed[i]     = False
                sim.people.infectious[i]  = False
                sim.people.symptomatic[i] = False
                sim.people.severe[i]      = False
                sim.people.critical[i]    = False
                sim.people.quarantined[i]   = False
                sim.people.symptomatic[i]   = False
                sim.people.date_recovered[i] = np.nan  
                sim.people.date_pos_test[i] = np.nan
#                sim.people.date_recovered[i] = -sim.pars['quar_period'] #np.nan   
            elif person['PositiveCollectDate'] > -sim.pars['quar_period'] \
                 and person['PositiveCollectDate'] <= 0:
                # INFECTED PERSON
                sim.people.date_pos_test[i] = person['PositiveCollectDate']
                sim.people.exposed[i]     = True
                sim.people.diagnosed[i]   = True   
                sim.people.severe[i]      = False  # assumption 
                sim.people.critical[i]    = False  # assumption
                sim.people.infectious[i]  = True   # assumption
                sim.people.recovered[i]   = False
                sim.people.quarantined[i]   = False                     
                # Estimated date they recover
                sim.people.date_recovered[i] = date_recovered
                sim.people.dur_disease[i] =  dur_asym2rec
                if  person['DiagnoseDate'] != None:
                    sim.people.date_diagnosed[i] = person['DiagnoseDate']
                else:
                    # No diagnose date is an error. Assign one.
                    sim.people.date_diagnosed[i] = person['PositiveCollectDate'] + 1
                # Maybe the exposure date should be half the test interval plus the days it takes
                # for a test to be positive.
                date_exposed = person['PositiveCollectDate'] - int(sim.people.labTestEveryNHours[i] / 2 / 24)
                sim.people.date_exposed[i] = date_exposed
                count += 1
            elif person['PositiveCollectDate'] > 0: 
                print('future')
                # They're going to get a positive collection date when the sim gets there.
                # *******************************
                # This would be best handled via an intervention that infects people on the assigned
                # day.
                # *******************************

                #sim.people.date_pos_test[i] = person['PositiveCollectDate']
                ### EDIT TO CONVERT TO DAYS NOT HOURS
                #date_exposed = person['PositiveCollectDate'] - int(sim.people.labTestEveryNHours[i] / 2 / 24)
                # Was this exposure date at or after the start of the simulation?
                #if date_exposed < 0:
                #    sim.people.date_exposed[i] = np.nan 
                #else:
                #    sim.people.date_exposed[i] = date_exposed
                pass
        if person['RecoverDate'] != None:
            # Set to 0.
            sim.people.date_recovered[i] = 0 #person['RecoverDate']
            # And make sure that covasim knows they recovered
            sim.people.recovered[i]   = True
            sim.people.susceptible[i]    = False
            sim.people.exposed[i]     = False
            sim.people.infectious[i]  = False
            sim.people.symptomatic[i] = False
            sim.people.severe[i]      = False
            sim.people.critical[i]    = False
            sim.people.diagnosed[i]   = False
            sim.people.quarantined[i] = False
        if person['LastTestCollectDate'] != None:
                sim.people.date_tested[i] = person['LastTestCollectDate']
                
    #    if person['DiagnoseDate'] != None:
      #      if person['PositiveCollectDate']  != None:
     #           if person['PositiveCollectDate'] > -sim.pars['quar_period']:
     #               count2 += 1
     #               sim.people.date_diagnosed[i] = person['DiagnoseDate']
     #           else:
                    
        if person['LatestSymptomaticDate'] != None:
         #   if person['LatestSymptomaticDate'] < 0:
         #       sim.people.date_symptomatic[i] = person['LatestSymptomaticDate']
         #   else:
                sim.people.date_symptomatic[i] = person['LatestSymptomaticDate']
        if person['EndQuarantineDate'] != None:
            if person['EndQuarantineDate'] < 0:
                sim.people.date_end_quarantine[i] = person['EndQuarantineDate']
                # NO date_end_quarantine set here - this person left before
                # the start of the simulation so leave that value be.
                sim.people.quarantined[i] = False
            else:
                sim.people.date_end_quarantine[i] = person['EndQuarantineDate']
                sim.people.quarantined[i] = True
        if person['DateKnownContact'] != None:
            sim.people.known_contact[i] = True
            # If they have not tested positive
            if person['PositiveCollectDate'] is None:
                # They're not positive
                sim.people.date_known_contact[i] = person['DateKnownContact']
                # Quarantine them IF their quarantine started < 14 days ago
                if person['DateKnownContact'] > -sim.pars['quar_period']:
                    if person['DateKnownContact'] == 0:
                        # On the start day use the built-in mechanism to quarantine someone.
                        sim.people.schedule_quarantine([i], start_date=person['DateKnownContact'])
                    else: 
                        # Use schedule_quarantine but reduce the period from the default because they started
                        # quarantine before the simuation started.
                        quar_days = person['DateKnownContact'] + sim.pars['quar_period']
                        sim.people.schedule_quarantine([i], start_date=0,period = quar_days)
                        #  self.date_end_quarantine[ind] = max(self.date_end_quarantine[ind], end_day) # Extend quarantine if required
                        # # quarantine manually
                        # sim.people.quarantined[i] = True 
                        # # And set their end quarantine day 
                        # sim.people.date_end_quarantine[i] = 0 #person['DateKnownContact'] + sim.pars['quar_period']
        if person['DateOfDeath'] != None:
                sim.people.date_dead[i] = person['DateOfDeath']
                sim.people.exposed[i]     = False
                sim.people.infectious[i]  = False
                sim.people.symptomatic[i] = False
                sim.people.severe[i]      = False
                sim.people.critical[i]    = False
                sim.people.recovered[i]   = False
                sim.people.dead[i]        = True      
         #   else:
         #       sim.people.date_dead[i] = person['DateOfDeath']  

             