"""
19-03-2019 JCollin
Script to ensure that vertical velocities written are good
To compute Omega we need to take a slice at given time and reshape it Python
!!! With W_sed !=0 you can't relate pdepth = int(pw.dt) ...
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../Modules/')


config = 'Write_uvw'

folder_root = '/home/jeremy/Bureau/Data/Pyticles/'
ncfile = folder_root + config + '/Case_1_' + config +'_12_1550.nc'


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
# ------ TESTING -------

w_sed0 = get_attr('w_sed0',ncfile)
print(w_sed0)

#############################################################################
pdepth = get_var('pdepth', ncfile)
pw = get_var('pw', ncfile)
ip = 15

def show_pw_pdeth(ip):

    xpart = np.arange(0,10, 1)

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('time')
    ax1.set_ylabel('pdepth', color=color)
    ax1.plot(xpart, pdepth[:,ip], color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_title('')

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('pw', color=color)  # we already handled the x-label with ax1
    ax2.plot(xpart, pw[:,ip], color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()


for ip in range(1,20,2):
    print(ip)
    show_pw_pdeth(ip)
# -----------------------------------------------
# Recomputin w using pyticles

import pyticles_3d_sig_sa as partF
import pyticles_sig_sa as part
from R_files import load

start_file = 1550
parameters = 'Case_1 [0,10000,0,10000,[1,100,1]] '+ format(start_file)

simul = load(simul = parameters, floattype=np.float64)
coord = simul.coord
[z_r,z_w] = part.get_depths(simul,coord=coord)
z_w2 = part.get_depths_w(simul,coord=coord)
np.shape(z_w2)

pm, pn = simul.pm, simul.pn

romsfile = folder_root + 'chaba_his.1550.nc' 

u = get_var('u', romsfile)
v = get_var('v', romsfile)

x_periodic , y_periodic = 0,0 
ng = get_attr('ng', ncfile)

nxp = 510; nxm = 500; nyp = 700; nym = 710
tight_coord = [nxp, nxm, nyp, nym]

[upart,vpart,wpart] = part.get_vel_io(simul,x_periodic=x_periodic,
        y_periodic=y_periodic,ng=ng,coord=tight_coord)

w = partF.get_omega(u,v,z_r,z_w,pm,pn)







