# file seeding_check.py
# J.Collin 05-03-2019
# 
# Script to ensure the seeding particles is made correctly
#
# Check whether initial position px, py, pz of seeding particules is written
# In netCDF file using seeding_part module routine 
# 
# Also check wether pdeth0 is the real pdeth0
# Numerical error of interpolation scheme may look like particles were advected
# 2 silation with different w_sed (25 m/s and 250 m/s) are compared
# They have the exact same pdeth0 proving depth0 is a numerical error
# For instance when seeding depth = -500 m  interpolation on sigma level and
# back to z(m) raise a +/- 20 m numerical error using linear interpolation 

#=================== LOADING MODULES #########################################
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.append("../Modules/")

#Specific modules needed for pyticles
import pyticles_sig_sa as part
import pyticles_3d_sig_sa as partF

import seeding_part

from R_files import load

from copy import *

#=========== USER PARAM =================
config = 'Seed_Test'

folder_root = '/home/jeremy/Bureau/Data/Pyticles/'
folder_save = folder_root + config

ncfile_p3 = folder_root + config + '/Case_1_' + config  + '_12_1550.nc'
ncfile_bis =  folder_root + config + '/Case_1_Seed_Test_WSed_25m_12_1550.nc'

# ========== Func =======================
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
# 

ng = 1 
adv3d = True
x_periodic = False
y_periodic = False


# reconstructing particles seeding from pyticles
#####################
# Particles initial location in the horizontal
# Need: to load pm form roms

nqmx = 100
start_file = 1550
parameters = 'Case_1 [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

# Mask
mask = simul.mask
maskrho = copy(mask)
maskrho[np.isnan(maskrho)] = 0.

depths = simul.coord[4]
nz = len(depths)
print(f'depths = {depths}')
ng = 1 


#
[ic,jc] = [600,800] #= part.find_points(simul.x,simul.y,-32.28,37.30)

# distance between 2 particles [in m]
dx_m = 1000.

dx0 = dx_m * simul.pm[ic,jc] # conversion in grid points

iwd  = 50.* dx0 # half width of seeding patch [in grid points]
jwd  = 50.* dx0 # half width of seeding patch [in grid points]

#########
# density of pyticles (1 particle every n grid points)
nnx=dx0*5
nny=dx0
nnlev=1.

(nx,ny) = simul.pm.shape

##########################
# first and last vertical level to fill with pyticlesi
intial_depth = True
depths = -500
lev0= 1; lev1= 1 
#if not adv3d: lev0 = -1; lev1 = lev0
#########
# define initial vertical position using depth
initial_depth = True
depths0 = [-500] # [-50, -100, -200]
if initial_depth:
    lev1 = lev0 + len(depths0) - 1
    nnlev = 1

###################################################################################
###################################################################################
# INITIALIZATION
###################################################################################
# Note that px,py,pz are ZERO-BASED !!!
# px,py,pz directly correspond to the indices of a python array (var[px,py,pz])
###################################################################################
restart = False
if not restart:
    ###################################################################################
    # Define initial px,py,pz pyticles position (Fast .py version_ fill in order x,y,z)
    ###################################################################################

    z,y,x = np.mgrid[lev0:lev1+1:nnlev,np.max([jc-jwd,1]):np.min([jc+jwd+np.min([1.,jwd]),ny]):nny, np.max([ic-iwd,1]):np.min([ic+iwd+np.min([1.,iwd]),nx]):nnx]

    if initial_depth: #initial vertical position = depths0
        from scipy.interpolate import interp1d
        z_w = part.get_depths_w(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)
        z = seeding_part.ini_depth(maskrho,simul,depths0,x,y,z,z_w)

    nq = np.min([len(x.reshape(-1)),nqmx])


        ###################################################################################
    ''' no need for topocheck anymore as we are using sigma levels'''
    ''' but we have to remove pyticles which are in a masked area'''

    ipmx = 0; px0,py0,pz0 = [],[],[]

    topolim=0

    if not adv3d: topolim = np.nanmax([topolim,-advdepth])

    # if you want to add a condition based on temp and/or salt:
    #[temp,salt] = part.get_ts_io(simul)

    ptopo = part.map_topo(simul,x.reshape(-1),y.reshape(-1))
    pmask = part.map_var2d(simul,maskrho,x.reshape(-1),y.reshape(-1))
    #ptemp = part.map_var(simul,temp,x.reshape(-1),y.reshape(-1),z.reshape(-1))
    #psalt = part.map_var(simul,salt,x.reshape(-1),y.reshape(-1),z.reshape(-1))

    for ip in range(len(x.reshape(-1))):
        if (ptopo[ip]>topolim) and (pmask[ip]>=1.) and (ipmx<nq):
            px0.append(x.reshape(-1)[ip])
            py0.append(y.reshape(-1)[ip])
            pz0.append(z.reshape(-1)[ip])
            ipmx +=1

    #del temp,salt
    nq = ipmx

    del x,y,z

########### Comparison between t[0] form output and from model

px_file = get_var('px', ncfile_p3)
py_file = get_var('py', ncfile_p3)
pz_file = get_var('pz', ncfile_p3)
pdepth_file = get_var('pdepth', ncfile_p3)

px_file.shape
err = np.empty(px_file.shape[1])

max_err = max(abs(px_file[0,:] - px0[:] ))
print(f'max error px = {max_err}')

max_err = max(abs(py_file[0,:] - py0[:] ))
print(f'max error py = {max_err}')

max_err = max(abs(pz_file[0,:] - pz0[:] ))
print(f'max error pz = {max_err}')


########### ONLY TO CHECK WETHER PDEPTH 0 is ok ##########
pdepth_test = get_var('pdepth', ncfile_bis)
diff_pdepth = pdepth_test - pdepth_file
np.max(abs(diff_pdepth[0,:]))
ptopo = get_var('ptopo', ncfile_bis)

xpart = np.arange(0, 100, 1)

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('particles')
ax1.set_ylabel('pdepth(t=0)', color=color)
ax1.plot(xpart, pdepth_test[0,:], color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_title('Numerical Error at seeding time for depths0 = -500 m')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('ptopo', color=color)  # we already handled the x-label with ax1
ax2.plot(xpart, ptopo[0,:], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()





