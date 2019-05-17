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
from R_tools import rho1_eos, rho_eos

##############################################################################

roms_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_his.1550.nc'
py_file = '/home/jeremy/Bureau/Data/Pyticles/Rho1_-1.5/Case_1_Rho1_-1.5_6_1510.nc'
start_file = 1550
parameters = 'Case_1 [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

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
px = get_var('px', py_file)
py = get_var('py', py_file)
pz = get_var('pz', py_file)
w_sed0 = get_attr('w_sed0', py_file)
rho = seeding_part.prho(ptemp=ptemp, psalt=psalt, pdepth=pdepth)
######
x_periodic = False
y_periodic = False
ng = 1
[temp, salt] = part.get_ts_io(simul, x_periodic = x_periodic,
                              y_periodic = y_periodic, ng=ng)

[z_r, z_w] = part.get_depths(simul)
roms_rho0 = simul.rho0
rho1 = rho1_eos(temp, salt, z_r, z_w, roms_rho0)

prho1 = part.map_var(simul, rho1, px.reshape(-1), py.reshape(-1),
                pz.reshape(-1), ng=ng).reshape(px.shape)

###########
bins = 20

fig = plt.figure()
ax1 = plt.subplot(221)
ax1.hist(prho1[-1,:]-prho1[0,:], bins=bins)
ax1.set_ylabel('prho1 tend - t0')

ax2 = plt.subplot(222)
ax2.hist(pdepth[0,:], bins=bins)
ax2.set_ylabel('pdepth at = 0')

ax3 = plt.subplot(223)
ax3.plot(psalt, ptemp)
ax3.set_ylabel('TS diagram')

plt.tight_layout()
plt.show()



############ OLD DEPRECIATED
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








