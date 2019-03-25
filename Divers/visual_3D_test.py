'''
visualisation some partilces data in 3D
'''

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../Modules/')

import pyticles_3d_sig_sa as partF
import pyticles_sig_sa as part
from R_files import load

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

##############################################################################

config = 'Write_uvw'

folder_root = '/home/jeremy/Bureau/Data/Pyticles/'
ncfile = folder_root + config + '/Case_1_' + config +'_12_1550.nc'

start_file = 1550
parameters = 'Case_1 [0,10000,0,10000,[1,100,1]] '+ format(start_file)

simul = load(simul = parameters, floattype=np.float64)


def get_var(var, ncfile, it=0):
    '''
    Returns the netcdf variable 'var' from ncfile
    var is a string : name of pyticiles  variable
    return np.array py_var
    use module netCDF4
    use module numpy as np
    it is used for time slicing 
    '''
    nc = Dataset(ncfile, 'r')
    if var in nc.variables:
        if var in ['u', 'v']:
            py_var = nc.variables[var][it,:,:,:]
        else:
            py_var = nc.variables[var][:]
    else:
        py_var = []
        print(f'Error {var} is not found in file {ncfile}')
    nc.close()
    return py_var
#################################################################

dxm = 1./np.mean(np.mean(simul.pm))
px = get_var('px', ncfile)*dxm/1000
py = get_var('py', ncfile)*dxm/1000
pz = get_var('pz', ncfile)

pdepth = get_var('pdepth', ncfile)

fig = plt.figure()
ax =fig.add_subplot(111, projection='3d')

for ip in range(1,100):
   # ax.plot(px[:,ip], py[:,ip], pdepth[:,ip])
    ax.contour(px,py,pdepth,cmap=cm.coolwarm)
plt.show()



