'''

To be loaded at the beginning of a python script:

from Modules import *

'''
from __future__ import print_function

# some ROMS tools written in python
import R_tools_gula as tools_g

# some ROMS tools written in .F
import R_tools_fort_gula as toolsF_g

# some ROMS tools written in .F
#import R_tools_fort_cuc as toolsF_cuc

# Compute variables
from R_vars_gula import var


import simulations_old as oldsim
import romstools_old as oldroms



from scipy.ndimage.measurements import label

try:
    from mpl_toolkits.basemap import Basemap
except:
    print('no basemap module installed')
    
    
'''
# Want to modify and reload a module?
del sys.modules["R_vars_gula"]
from R_vars_gula import var
'''

