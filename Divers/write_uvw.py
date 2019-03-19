"""
19-03-2019 JCollin
Script to ensure that vertical velocities written are good
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt



config = 'Write_uvw'

folder_root = '/home/jeremy/Bureau/Data/Pyticles/'
ncfile = folder_root + config + '/Case_1_' + config +'_1_1550.nc'

def get_var(var, ncfile):
    '''
    Returns the netcdf variable 'var' from ncfilei
    var is a string : name of pyticiles  variable
    return np.array py_var
    use module netCDF4
    use module numpy as np
    '''
    nc = Dataset(ncfile, 'r')
    if var in nc.variables:
         py_var = nc.variables[var][:]
    else:
        py_var = []
        print(f'Error {var} is not found in file {ncfile}')
    nc.close()
    return py_var


pdepth = get_var('pdepth', ncfile)
pw = get_var('pw', ncfile)


xpart = np.arange(0,5, 1)

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('time')
ax1.set_ylabel('pdepth', color=color)
ax1.plot(xpart, pdepth[:,0], color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_title('')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('pw', color=color)  # we already handled the x-label with ax1
ax2.plot(xpart, pw[:,0], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()


