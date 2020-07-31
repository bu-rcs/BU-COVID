#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 17:46:15 2020

Functions used for the BU Covasim modeling.



@author: bgregor
"""

import numpy as np
import datetime
import os 
import pandas as pd
import networkx as nx
import zipfile
import tqdm


import covasim as cv
 

__all__=['load_people_info','make_dt','get_daynum','gen_BU_pop2',
         'validate_graphs',
         'update_graph', 'write_graphml', 
         'update_sim_people','sim_results_to_df', ]


### Use covasim default types
def_int = cv.defaults.default_int
def_flt = cv.defaults.default_float

#%%
def load_people_info(pop_info_path):
    ''' Read the pop info file.  Header:
        
        id,age,sex,group,category,position,campus
        
        The last 2 fields are strings, the rest are ints.
        Turn the into a dictionary with id:{rest of data}
    
    '''
    lut = {} # The lookup table to return.
    with open(pop_info_path,encoding='latin1') as f:
        # read the header row for use as keys.  Strip whitespace if any.
        keys = [x.strip() for x in f.readline().split(',')]
        keys.pop(0) # remove the id header
        for num,line in enumerate(f):
            try:
                # Split the line & convert to integers
                tmp = [x.strip() for x in line.split(',')]
                for i,t in enumerate(tmp):
                    if t.isdigit():
                        tmp[i] = int(t)
                # get the person's uid ("unique id")
                uid = tmp.pop(0)
                # Now fill in the lut
                if uid not in lut:
                    lut[uid] = dict(zip(keys,tmp)) 
            except Exception as e:
                print('ERROR at line: %s' % num)
                print(e)
    return lut

#%%

def make_dt(day):
    ''' Make a Python datetime object. day is YYYY-MM-DD'''
    return datetime.datetime.strptime(day,'%Y-%m-%d') 


def get_daynum(dt):
    ''' dt: datetime object. Return the number of the day of the week'''
    return int(dt.strftime('%w'))




#%%
def gen_BU_pop2(people_lut,class_contacts={}, household_contacts={}):
    ''' Generate the population dictionary that's fed into covasim.
        people_lut: Lookup table:  {uid:{'age':X,'sex':Y...}
        class_contacts: classroom contact networks.
        housing_contacts: housing contact networks.
        return population dictionary:
    '''
    # Step 1. Find out how many people are in the contact lists.
    # Merge the dictionaries into a 3rd using a Python trick
    contact_dicts = {**class_contacts, **household_contacts}
    id_set = set()
    # loop through all keys in all  dictionaries
    for key in contact_dicts:
        # Add all the keys from the inner dictionaries
        id_set |= set(contact_dicts[key].keys())

    # id_set is now the unique list of BU person id's
    # used in the simulation.
    # Turn that into a lookup table for list index values
    # starting from 0. lut converts student id's into indices.
    lut = dict(zip(id_set,range(len(id_set))))
    
    # The length of lut is the population size.
    pop_size = len(lut)
    
    # Initialize the population object. Add our extra fields.
    # Also cut down on the size of our fields where possible.
    BU_pop = {'age':np.full(pop_size,0,dtype=def_flt), 
              'contacts':pop_size * [None], 
              'layer_keys':list(contact_dicts.keys()),
              'sex':np.full(pop_size,0,def_flt), 
              'uid':np.zeros(pop_size,dtype=def_int),
              'group':np.zeros(pop_size,dtype=np.uint8),
              'category':np.zeros(pop_size,dtype=np.uint8),
              'campResident':np.zeros(pop_size,dtype=np.uint8),
              'undergrad':np.zeros(pop_size,dtype=np.uint8),
             # 'position':pop_size * [''],
              #'campus':pop_size * [''],
              'full_info_id':np.zeros(pop_size,dtype=np.uint32)}
    
    # Loop over the people_lut
    # fill in the age & sex of everyone found there.
    for person in people_lut:
        if person in lut:
            # index found.  Fill in values.
            idx = lut[person]
            BU_pop['age'][idx] = people_lut[person]['age']  
            BU_pop['sex'][idx] = people_lut[person]['sex']  
            BU_pop['group'][idx] = people_lut[person]['group']  
            BU_pop['category'][idx] = people_lut[person]['category']
            BU_pop['undergrad'][idx] = people_lut[person]['undergrad']
           # BU_pop['position'][idx] = people_lut[person]['position']  
           # BU_pop['campus'][idx] = people_lut[person]['campus']  
            BU_pop['campResident'][idx] = people_lut[person]['campResident']  
            BU_pop['full_info_id'][idx] = person
            

    # Now loop through the lut.  For each student id fill their
    # contact lists. For each contact list convert from student ID
    # to index value.
    for stu in lut:
        idx = lut[stu]
        BU_pop['uid'][idx]=idx
        contacts = {}
        for i,key in enumerate(BU_pop['layer_keys']):
            if stu in contact_dicts[key]:
                # convert the student id contact list to indices
                tmp = [lut[x] for x in contact_dicts[key][stu]]
                contacts[key] = tmp
            else:
                contacts[key] = []
        BU_pop['contacts'][idx] = contacts                
    
    return BU_pop
 


   
#%%
def validate_graphs(graphs):
    ''' As a sanity check, make sure that each graph vertex
        matches all the others.  Compares the vertices of
        the first graph with the vertices of the same vertex
        number in all the other graphs.
        
        graphs: a list of graphs.
        returns: True or False.
        
    '''
    for g_ind in range(1,len(graphs)):
        for i,vtx in enumerate(graphs[0].vs):
            vtx_g = graphs[g_ind].vs[i]
            if vtx != vtx_g:
                print(vtx)
                print(vtx_g)
                return False
    return True


#%%
def update_graph(g, BU_pop):
    ''' Update a networkx graph with info from the BU_pop'''
    # Delete the pointless None node from the end if it's there
    try:
        g.remove_node(None)
    except:
        pass
    # This uses the networkx call: nx.set_node_attributes. This takes
    # the graph, a dictionary of {node_idx:value} and a label for the value.
    num =  BU_pop['uid'].shape[0]
    # Each person in people_lut has the same keys. Store those
    attrs =  list(BU_pop.keys())
    # Pull out some keys we don't need to add
    attrs.pop(attrs.index('age'))
    attrs.pop(attrs.index('uid'))
    attrs.pop(attrs.index('contacts'))
    attrs.pop(attrs.index('layer_keys'))

    # For each attritubte: generate a dictionary of node indices
    # vs value, and add with the attribute name
    for attr in attrs:
        data = {}
        for i, uid in enumerate(BU_pop['uid']):
            data[i] = BU_pop[attr][i]
        # Add this to the nodes
        nx.set_node_attributes(g, data, attr)
    # Now go through all nodes...for all attributes
    # find any NaN values and set them to -1 for friendlier 
    # processing in R.
    attrs = set([k for n in g.nodes for k in g.nodes[n].keys()])
    for attr in attrs:
        data = {}
        for n in g.nodes:
            if n:
                try:
                    if np.isnan(g.nodes[n][attr]):
                            data[n] = -1
                except:
                    print(n,attr)
                
        # Update the nodes
        nx.set_node_attributes(g, data, attr)
    # ANd that's it.

 
def write_graphml(sims_complete, plot_dir, sim_name, BU_pop):
    # Create and save the transmission trees for all sims.
    graph_names = []
    for i,sim in enumerate(sims_complete):
        ttree = cv.TransTree(sim, to_networkx=True)
        fname = os.path.join(plot_dir,'%s_%s.graphml' % (sim_name,i))
        graph_names.append(fname)
        # Update the ttree graph to replace NaN with -1
        update_graph(ttree.graph, BU_pop)
        nx.readwrite.graphml.write_graphml(ttree.graph,fname)
    
    # Now take all of the graphml  files and pack them into a ZIP
    zipname = os.path.join(plot_dir,'%s_graphmls.zip' % sim_name)
    with zipfile.ZipFile(zipname,'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as myzip:
        for graph in tqdm.tqdm(graph_names):
            myzip.write(graph)
            
    # Now remove all of the graphml files
    for graph in graph_names:
        os.unlink(graph)

 

#%%
def update_sim_people(sim,BU_pop):
    ''' Add our extra parameters to the sim.people object.
        sim should be an unintialized Sim object. '''
    # Initialize to create sim.people
    sim.initialize()
    # Add our parameters
    sim.people['undergrad'] = BU_pop['undergrad']
    sim.people['group'] = BU_pop['group']
    sim.people['category'] = BU_pop['category']
    sim.people['campResident'] = BU_pop['campResident']
    sim.people['full_info_id'] = BU_pop['full_info_id']
    # And that's it...

 
#%%
def sim_results_to_df(sims_complete):
    ''' Takes a list of completed sims.  Returns
        a Pandas dataframe of all of their results. The sim_num
        column indexes the simulations. '''
    data={'sim_num':[], 'dates':[], 'days':[]}
    # Take the 1st simulation and add all of its results keys to the dictionary
    sim0 = sims_complete[0]
    keys = list(sim0.results.keys())
    # Remove day num and date keys
    keys.pop(keys.index('t'))
    keys.pop(keys.index('date'))
    # Add the rest.
    for k in keys:
        data[k] = []
    
    for i, sim in enumerate(sims_complete):
        # Convert all quantities to python lists to avoid the
        # overhead of ndarray concatentation
        for k in keys:
            data[k] += list(sim.results[k])
        data['days'] += list(sim.results['t'])
        data['dates'] += list(sim.results['date'])
        data['sim_num'] += len(sim.results['date']) * [i]
        
    return pd.DataFrame(data=data)





