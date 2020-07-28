# Functions related to processing classroom and housing contact lists

import itertools as it
import datetime

from .misc import *

__all__=['gen_daily_interventions','get_contact_list','get_shared_class_contact_lists','get_all_class_contacts_dicts','get_housing_contacts','get_all_housing_contacts_dict']

#%%
def gen_daily_interventions(start_date,end_date,layers,layer_days,changes,intervention,holidays=None):
    ''' start_date:   starting date of time period in which interventions are generated. YYYY-MM-DD
                This should be before the 1st class day.
        end_date: end date the of the time period.  Inclusive, i.e 12/01/2020 includes 12/01  YYYY-MM-DD
        layers: list of layers that the intervention will apply to
        layer_days: List of numbers of the days when the layer is active. 0: Sunday, 1: Monday,...,6: Saturday.
        changes: A list of 2 numbers, the first is the change value when the layer is on, the 2nd is when the layer is off.
        intervention: an Intervention class
        holidays: list of dates that are holidays where the layer is off.  ['YYYY-MM-DD','YYYY-MM-DD',....] 
        
        return: the intervention
    '''        
    start_dt = make_dt(start_date) 
    end_dt = make_dt(end_date)
    holidays_dt = []
    if holidays:
        holidays_dt = [make_dt(h) for h in holidays]
    # Get the number of days between the start_date and the end_date
    num_days = (end_dt - start_dt).days
    # Make a 1 day time delta
    day_delta = datetime.timedelta(days = 1)
    
    # Initialize a list to hold each day of the intervention
    interv_dates = []
    # Initialize a list to hold the changes.
    interv_changes = []
    for t in range(num_days):
        today = start_dt + t * day_delta
        interv_dates.append(t)
        if get_daynum(today) in layer_days and today not in holidays_dt:
            # class is on
            interv_changes.append(changes[0])
        else:
            # no class
            interv_changes.append(changes[1])
    
    return intervention(interv_dates, interv_changes, layers)
        
#%%
def get_contact_list(g):
    ''' Return a list of contact lists from a graph.
        g is an igraph graph
        layer_name is the layer being constructed, 'sun','mon', etc '''
    edge_list = g.get_edgelist()
    edge_list += [(x[1],x[0]) for x in edge_list] #  Add the reverse connections
    edge_list.sort()
    pop_size = len(g.vs)
    # Initialize the list storage.  Include empty lists for days when the
    # contact list is empty.
    contact_list = pop_size * [[]]
    # The contact list has to
    for key, group in it.groupby(edge_list, lambda x: x[0]): 
        tmp = [x[1] for x in group]
        # This is a list of vertices. The vertex indices are the same for each
        # person across the graphs, i.e. g.vs[0] is always the same person.
        contact_list[key]=tmp
    return contact_list
         
#%%
def get_shared_class_contact_lists(g,layer,max_n):
    ''' Return a list of contact lists from a graph.
        g is an igraph graph 
        layer is the name of the layer
        max_n is the max number of people in a layer
        
        Returns a dictionary of max_n layers '''
    # Get the max weights in the edge list
    max_found = max((x['weight'] for x in g.es))
    if max_n > max_found:
        max_n = int(max_found)
    
    # create a dictionary of layer names
    contacts = {}
    layer_names = []
    pop_size = len(g.vs)
    for i in range(max_n):
        layer_names.append('%s_%s' % (layer,i+1))
        contacts[layer_names[i]] = {} # Build as an inner dictionary
        
    # Now build the edge lists.        
    edge_list = g.get_edgelist()
    els = {}
    for i in range(1,max_n+1):
        els[i] = []
        for j,edge in enumerate(edge_list):
            if int(g.es[j]['weight']) == i:
                els[i].append(edge)

    #Take the union of the lists. 1 shared contact includes the 2 & 3 contacts
    # 2 shared contacts includes 3, etc.

    els_keys = sorted(list(els.keys()))

    for i in range(len(els_keys)-1):
        key1 = els_keys[i]
        for j in range(i+1,len(els_keys)):
            key2 = els_keys[j]
            els[key1] += els[key2]

    # Then connect them both ways.
    for key in els:
        els[key] += [(x[1],x[0]) for x in els[key]]
        els[key] = sorted(els[key])

    # Fill in the contact list dictionary with the correct entries.
    for i in range(1,max_n+1):
        el = els[i]
        contact_list = contacts[layer_names[i-1]]
        for key, group in it.groupby(el, lambda x: x[0]): 
            tmp = [x[1] for x in group]
            # This is a list of vertices. Extract the student id and 
            # convert to integers to save RAM
            stu = int(g.vs[key]['full_info_id'])
            contact_list[stu]=[int(g.vs[x]['full_info_id']) for x in tmp]
    return contacts

 

#%%
def get_all_class_contacts_dicts(graphs, layer_names, max_n=1):
    '''  Return a complete set of dictionaries for all of the classes.
        graphs: a list of graphs.
        layer_names: list of names, one per graph
        max_n: max edge weighting for binning contact lists. '''
    class_contacts = []
    for i,g in enumerate(graphs):
        class_contacts.append(get_shared_class_contact_lists(g,layer_names[i],max_n))        
    all_dict = {}
    for d in class_contacts:
        all_dict.update(d)
    return all_dict


#%%
def get_housing_contacts(g):
    ''' Generate the contact lists for a housing network.
        g: a graph
        
        return: a dictionary of students and their housing contacts '''
    # get the edge lists.        
    edge_list = g.get_edgelist()
    edge_list += [(x[1],x[0]) for x in edge_list] #  Add the reverse connections
    edge_list.sort()

    # Initialize the list storage.  Include empty lists for days when the
    # contact list is empty.
    contacts = {}
    # TEMPORARY: Use one of 2 possible keys for the student id field
    stu_id = 'full_info_id'
    if stu_id not in g.vs.attribute_names():
        stu_id = 'studentID'
    for key, group in it.groupby(edge_list, lambda x: x[0]): 
        tmp = [x[1] for x in group]
        # This is a list of vertices. The vertex indices are the same for each
        # person across the graphs, i.e. g.vs[0] is always the same person.
        stu = int(g.vs[key][stu_id])
        contacts[stu]=[int(g.vs[x][stu_id]) for x in tmp]
    return contacts
         

#%%
def get_all_housing_contacts_dict(graphs, layer_names):
    ''' graphs is a list of housing graphs.
        layer_names: list of names, one per graph.
        
        return: dictionary per layer of all contact lists'''
    housing_contacts={}
    for i,g in enumerate(graphs):
        housing_contacts[layer_names[i]] = get_housing_contacts(g)
    return housing_contacts 


