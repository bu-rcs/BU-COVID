from bu_covid import misc
from bu_covid import housing_and_classrooms as hc
import bu_covid.exception as bu_e

import pytest 

import datetime 
import pandas as pd
import tempfile
import os 

def test_make_dt():
    ''' Make a Python datetime object. day is YYYY-MM-DD'''
    now = datetime.datetime.now() 
    # Was a datetime object returned?
    assert type(misc.make_dt('2020-01-01')) == type(now)

def test_make_dt_error():
    ''' Test exception throwing '''
    with pytest.raises(bu_e.TerrierException):
        misc.make_dt('2020') 


def test_get_daynum():
    ''' dt: datetime object. Return the number of the day of the week'''
    today = datetime.datetime(2020, 10, 7, 8, 50, 36, 898767)
    assert misc.get_daynum(today) == 3
    

def test_load_people_info():
    ''' Read the pop_info.csv file. Write a temp file and use that'''
    test_info='''id,age,sex,TestCategory,LabTestEveryNHours,Residence,Affiliation,Campus,Undergrad,Attendance
1,20,0,3,2600,0,1,1,0,1
2,45,0,4,0,0,4,1,0,1
3,20,0,1,96,0,1,1,1,1
4,35,0,4,0,0,2,1,0,1
5,20,0,4,0,0,1,1,0,0
6,20,0,1,96,1,1,1,1,0
7,45,0,2,168,0,3,2,0,1
8,20,0,1,96,1,1,1,1,1
9,35,0,4,0,0,3,1,0,1
10,45,0,4,0,0,1,1,0,1
11,45,0,3,2600,0,3,2,0,1
12,35,0,2,168,0,1,1,0,1
13,45,0,1,168,0,2,1,0,1
14,20,0,1,96,1,1,1,1,1
15,20,0,1,96,1,1,1,1,1
16,20,0,1,96,1,1,1,1,1
17,20,0,1,96,1,1,1,1,1
18,25,0,2,168,0,1,1,0,1
19,25,0,2,168,0,1,1,0,1
20,45,0,4,0,0,4,1,0,1
21,20,0,1,96,0,1,1,1,1
22,20,0,1,96,1,1,1,1,1
23,45,0,4,0,0,4,1,0,1
24,55,0,4,0,0,3,2,0,1
25,20,0,1,96,1,1,1,1,1
26,20,0,1,96,1,1,1,1,1
27,45,0,4,0,0,3,2,0,1
28,35,0,4,0,0,3,2,0,1
29,55,0,4,0,0,2,1,0,1
'''
    tfile = tempfile.NamedTemporaryFile(mode='w',delete=False)
    tfile.write(test_info)
    tfile.close()
    pop = misc.load_people_info(tfile.name)
    os.unlink(tfile.name)
    keys = ['age','sex','TestCategory','LabTestEveryNHours','Residence',
            'Affiliation','Campus','Undergrad','Attendance']
    for person in pop:
        for key in keys:
            assert key in set(pop[person].keys())
            
    # Right number of lines?
    assert len(pop.keys()) == 29


def test_gen_BU_pop2():
    ''' Test the creation of a people dictionary.'''
    pop_info_file=''
    # Load the pop info file
    #pop = misc.load_people_info(pop_info_file)
    # Load networks
    # ...
    
    # Now create a pop dictionary
    #BU_pop = misc.gen_BU_pop2(pop,class_contacts,housing_contacts)
    
    # tests: BU_pop has the correct number of keys
    # spot check some networks to see if the right ones were assigned
    
  
def test_update_sim_people():
    # create a sim object and a BU_pop object.
    # Call update_sim_people() and double-check that
    # the new categories have been added. Also check
    # that the new categories exist properly in the 
    # sim meta attributes.
    pass
    
def test_sim_results_to_df():
    # create a pair of sim objects. Run both, add to a
    # multisim.  Send the multisim to sim_results_to_df()
    # Check the returned pandas dataframe for correctness
    pass

        
