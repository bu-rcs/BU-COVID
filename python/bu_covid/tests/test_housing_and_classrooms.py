from bu_covid import misc
from bu_covid import housing_and_classrooms as hc
import covasim as cv

import pytest

# Read in sample data structures, they'll be needed in the tests.
# 

# Location of housing networks
hnet_dir = '../../../Data/networks'

# some classroom graphs
@pytest.fixture
def get_classroom_graphs():
    from igraph import Graph
    import os   
    net_dir = '../../../Data/networks'
    class_files = ['ClassNetPlatoonM.graphml','ClassNetPlatoonT.graphml',]
    graphs = [Graph.Read_GraphML(os.path.join(net_dir, cf)) for cf in class_files]
    class_layer_names = ['mon','tue'] 
    return graphs, class_layer_names

@pytest.fixture
def get_housing_graphs():
    from igraph import Graph
    import os   
    # some housing graphs
    h_files = ['RoomNet.graphml','FloorNet.graphml']
    hgraphs=[Graph.Read_GraphML(os.path.join(hnet_dir, hf)) for hf in h_files]
    hlayer_names = ['roommate','floor']
    return hgraphs, hlayer_names


def test_gen_daily_interventions():
    start_date = '2020-01-01'
    end_date = '2020-01-31' 
    layers = ['a','b']
    # Mon, Wed, Fri
    layer_days = [1,3,5]
    # On on layer_days, otherwise off.
    changes = [1.0,0.0] 
    intervention = cv.clip_edges
    holidays = ['2020-01-01','2020-01-20'] # this is a Wednesday and a Monday.
    tmp = hc.gen_daily_interventions(start_date,end_date,layers,layer_days,changes,intervention,holidays = holidays)
    # Check some things...
    # correct layers
    assert layers == tmp.layers
    # correct number of days.
    assert len(tmp.days) == 30
    # Starting days set correctly? New Year's --> days[0] should be 0.
    assert tmp.changes[0:6] == [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
    # MKL day is zero? that's the 20th --> days[19]
    assert tmp.changes[18:25] == [0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0]
 
 
def test_get_contact_list(get_classroom_graphs):
    graphs, class_layer_names = get_classroom_graphs 
    # use the Monday classroom graph
    g = graphs[0]
    tmp = hc.get_contact_list(g)
    assert len(tmp) == 3681
    assert sum([len(y) for y in tmp]) == 27650
    assert tmp[0] == [3, 7, 42, 62, 68, 83, 94, 100, 121, 148, 219, 233, 246, 265, 331, 385, 418, 446, 471, 481, 818, 1004, 1018, 1146, 1153, 1281, 1302, 1356, 1492, 3540]
 

def test_get_shared_class_contact_lists(get_classroom_graphs):
    graphs, class_layer_names = get_classroom_graphs 
    # use the Monday classroom graph
    g = graphs[0]
    tmp = hc.get_shared_class_contact_lists(g,'mon',3)
    # check that everything is of the correct length
    assert list(tmp.keys()) == ['mon_1', 'mon_2']
    assert len(tmp['mon_1'].keys()) == 1844
    # spot check
    assert tmp['mon_1'][3693] == [589, 1085, 1095, 1141, 1142, 1245, 1268, 1287, 1423, 1463]

    
def test_get_all_class_contacts_dicts(get_classroom_graphs):
    graphs, class_layer_names = get_classroom_graphs 
    tmp = hc.get_all_class_contacts_dicts(graphs, class_layer_names, max_n = 3)
    assert list(tmp.keys()) == ['mon_1', 'mon_2', 'tue_1', 'tue_2']
    assert len(tmp['tue_1'].keys()) == 1720
    assert  tmp['tue_2'][1387] == [2223]

def test_get_housing_contacts(get_housing_graphs):
    hgraphs, hlayer_names = get_housing_graphs 
    # use roommates
    tmp = hc.get_housing_contacts(hgraphs[0])
    assert len(tmp.keys()) == 2100
    assert tmp[716] == [216]

def test_get_all_housing_contacts_dict(get_housing_graphs):
    hgraphs, hlayer_names = get_housing_graphs 
    tmp = hc.get_all_housing_contacts_dict(hgraphs, hlayer_names)
    assert list(tmp.keys()) == ['roommate', 'floor'] 
    assert len(tmp['floor']) == 2100
    assert tmp['floor'][3476] == [3012, 3032, 3033, 3034, 3045, 3050, 3054, 3070, 3074, 3094, 3098, 3113, 3142, 3151, 3159, 3163, 3167, 3174, 3183, 3193, 3200, 3215, 3227, 3228, 3238, 3254, 3279, 3291, 3301, 3313, 3316, 3321, 3325, 3338, 3342, 3344, 3352, 3355, 3357, 3383, 3396, 3406, 3413, 3447, 3448, 3459, 3469, 3473, 3500]
    
    

 





