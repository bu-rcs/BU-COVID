#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 11:16:58 2020

@author: bgregor

Add code to support modifying a Sim object before it's run to reflect the
current state of the community.
"""

import covasim as cv
import sciris as sc
import numpy as np
from .misc import make_dt

__all__ = ['import_community_test_info']


def parse_test_info_file(info_file, start_day):
    ''' Read the CSV pop_test_info.csv into a dict of lists.
        state_file: CSV file of the current state of the community
        start_day: start date string for the entire simulation. 
    '''
    # header: id,TestCategory,LabTestEveryNHours,PositiveInd,PositiveCollectDate,LastTestCollectDate,LastTestResultDate,DiagnoseDate,RecoverDate,EarliestSymptomaticDate,LatestSymptomaticDate,EndQuarantineDate,DateKnownContact,DateOfDeath
    start_ind = 0
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
    start_dt = make_dt(sim.date(0))
    for i,info_id in enumerate(sim.people.full_info_id):
        try:
            if info_id not in id_data:
                print('ERROR: full_info_id %s not found in %' % (info_id,info_file))
                continue
            person = id_data[info_id]
            # Guard each with an if statement - missing dates are None values in 
            # the id_data dictionaries
            if person['PositiveCollectDate']:
                # Mark as not susceptible
                sim.people.susceptible[i] = False
                # If there's a positive collection date before the start of the simulation 
                # set it to zero.
                if person['PositiveCollectDate'] < 0:
                    sim.people.date_pos_test[i] = 0
                # Check when they were positive. If more than 14 days
                # before the start date then mark as recovered.
                if person['PositiveCollectDate'] <= -sim.pars['quar_period']:
                    sim.people.recovered[i]   = True
                    sim.people.exposed[i]     = False
                    sim.people.infectious[i]  = False
                    sim.people.symptomatic[i] = False
                    sim.people.severe[i]      = False
                    sim.people.critical[i]    = False
                    # Set the recovery date to the start of the simulation
                    sim.people.date_recovered[i] = 0   
                    
                else: # They're going to get a positive collection date when the sim gets there. Let's
                    sim.people.date_pos_test[i] = person['PositiveCollectDate']
                    date_exposed = person['PositiveCollectDate'] - int(sim.people.labTestEveryNHours[i] / 2)
                    # Was this exposure date at or after the start of the simulation?
                    if date_exposed < 0:
                        sim.people.date_exposed[i] = 0
                    else:
                        sim.people.date_exposed[i] = date_exposed
            if person['RecoverDate']:
                # They have a recovery date.  If before the start of the simulation,
                # set it to 0. Otherwise, set it.
                if person['RecoverDate'] < 0:
                    sim.people.date_recovered[i] = 0
                else:
                    sim.people.date_recovered[i] = person['RecoverDate']
                # And make sure that covasim knows they recovered
                sim.people.recovered[i]   = True
                sim.people.exposed[i]     = False
                sim.people.infectious[i]  = False
                sim.people.symptomatic[i] = False
                sim.people.severe[i]      = False
                sim.people.critical[i]    = False
    
            if person['LastTestCollectDate']:
                if person['LastTestCollectDate'] < 0:
                    sim.people.date_tested[i] = 0
                else:
                    sim.people.date_tested[i] = person['LastTestCollectDate']
            if person['DiagnoseDate']:
                if person['DiagnoseDate'] < 0:
                    sim.people.date_diagnosed[i] = 0 
                else:
                    sim.people.date_diagnosed[i] = person['DiagnoseDate']
    
            if person['LatestSymptomaticDate']:
                if person['LatestSymptomaticDate'] < 0:
                    sim.people.date_symptomatic[i] = 0
                else:
                    sim.people.date_symptomatic[i] = person['LatestSymptomaticDate']
            if person['EndQuarantineDate']:
                if person['EndQuarantineDate'] < 0:
                    sim.people.quarantined[i] = False
                    # NO date_end_quarantine set here - this person left before
                    # the start of the simulation so leave that value be.
                else:
                    sim.people.date_end_quarantine[i] = person['EndQuarantineDate']
                    sim.people.quarantined[i] = True
            if person['DateKnownContact']:
                sim.people.date_known_contact[i] = person['DateKnownContact']
            if person['DateOfDeath']:
                if person['DateOfDeath'] < 0:
                    sim.people.date_dead[i] = 0
                    sim.people.exposed[i]     = False
                    sim.people.infectious[i]  = False
                    sim.people.symptomatic[i] = False
                    sim.people.severe[i]      = False
                    sim.people.critical[i]    = False
                    sim.people.recovered[i]   = False
                    sim.people.dead[i]        = True      
                else:
                    sim.people.date_dead[i] = person['DateOfDeath']  
        except Exception as e:
            print("Error for person %s: %s" % (info_id,e))