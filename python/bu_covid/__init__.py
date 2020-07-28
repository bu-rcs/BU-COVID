# Load the custom BU code for covasim

# This MUST COME FIRST in order to monkeypatch the covasim People class.
import covasim as cv

# =============================================================================
# Tweak the People class to fit our needs
# =============================================================================
# Add some people characteristics to the PeopleMeta class
BU_attrs = ['undegrad','group','category', 'campResident', 'full_info_id']
cv.defaults.PeopleMeta.person += BU_attrs
cv.defaults.PeopleMeta.all_states += BU_attrs



from .misc  import *
from .housing_and_classrooms import * 
from .parallel  import *
from .snapshots import *
from .interventions import *
 
