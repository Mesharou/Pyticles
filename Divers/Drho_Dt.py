'''
check on DRho/Dt to check if there is not too much variation within the ocean's interior

'''

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as col

from netCDF4 import Dataset

import numpy as np
import sys
from scipy.interpolate import interp1d



sys.path.append("../Modules/")
import seeding_part
from R_files import load
import pyticles_sig_sa as part
import pyticles_3d_sig_sa as partF

roms_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_his.1550.nc'
py_file = '/home/jeremy/Bureau/Data/Pyticles/Iso_surf/Case_1_Iso_surf_12_1550.nc'

#############################################################################
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
#############################################################################

def get_attr(attr, ncfile, it=0):
    '''
    Returns the netcdf attribute 'attr' from ncfile
    attr is a string : name of pyticilesgolabal attribute
    '''
    nc = Dataset(ncfile, 'r')
    if attr in nc.ncattrs():
            py_attr = getattr(nc, attr)
    else:
        py_attr= []
        print(f'Error {attr} is not found in file {ncfile}')
    nc.close()
    return py_attr

##############################################################################

ptemp = get_var('pt', py_file)
psalt = get_var('ps', py_file)
pdepth = get_var('pdepth', py_file)
pz = get_var('pz', py_file)
w_sed0 = get_attr('w_sed0', py_file)
rho = seeding_part.prho(ptemp=ptemp, psalt=psalt, pdepth=pdepth)


############ 
plt.subplot(221)
plt.hist((rho[-1, :] - rho[0, :])/2, bins=20)
plt.title(f'D/Dt rho kg.m^-3.day^-1')
plt.ylabel('rho0 = 1028 kg.m^-3')

plt.subplot(222)
plt.plot(psalt, ptemp)
plt.title('diagram TS along particles trajectories')
plt.ylabel('Temperature')
plt.xlabel('Salinty')
plt.grid(True)

plt.subplot(223)
plt.hist(rho[0,:])
plt.xlabel(f'rho at t=0 kg.m^-3.day^-1')

plt.subplot(224)
plt.hist(pdepth[-1, :])
plt.xlabel('pdepth after 2 days')
plt.show()

###################








