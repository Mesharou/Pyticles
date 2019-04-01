# file seeding_check.py
# J.Collin 28-03-2019
# 
# Script to ensure the seeding particles is made correctly
#
# Also check wether pdeth0 is the real pdeth0
# Numerical error of interpolation scheme may look like particles were advected
# 2 silation with different w_sed (25 m/s and 250 m/s) are compared
# They have the exact same pdeth0 proving depth0 is a numerical error
# For instance when seeding depth = -500 m  interpolation on sigma level and
# back to z(m) raise a +/- 20 m numerical error using linear interpolation 
#
#
# Methdology 
# Firstly create a seeding patch that is exactly at psi_w points 
# Then we interpolate z_w (rho_points) onto psi_w 
# Then we compute px0, py0, pz0 using good values z_w_psi
# In the end we use partf.inerpolate_wpsi in order to get pdepth
#
# Then we check that the Error is really small even around topographic gaps
# 
# Secondly generalize the method in order to support any px,py...
#
# Update
# First of all z_w was initialized at rho points not psi_rho
# Then it seems that using ng = 0 works fine...
# Another issue to solve : it only works for x,y,z on psi_rho points
# Therefore we need to use some 3D interpolation to raise this issue
# 
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

nqmx = 10000
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
#dx0 = 1

iwd  = 20 # half width of seeding patch [in grid points]
jwd  = 20  # half width of seeding patch [in grid points]

#########
# density of pyticles (1 particle every n grid points)
nnx=dx0
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
    ###################################################################################
    # Define initial px,py,pz pyticles position (Fast .py version_ fill in order x,y,z)
    ###################################################################################

z,y,x = np.mgrid[lev0:lev1+1:nnlev,
        np.max([jc-jwd,1]):np.min([jc+jwd+np.min([1.,jwd]),ny]):nny,
        np.max([ic-iwd,1]):np.min([ic+iwd+np.min([1.,iwd]),nx]):nnx]

from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator

z_w = part.get_depths_w(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)
z_tmpX = (z_w[:-1,:] + z_w[1:,:]) / 2
z_wpsi = (z_tmpX[:,:-1] + z_tmpX[:,1:]) / 2    
 
#xx = np.arange(simul.coord[2], simul.coord[3] - 1)
#yy = np.arange(simul.coord[0], simul.coord[1] - 1)
#zz = np.arange(0,nz+1)
#my_interpolating_function = RegularGridInterpolator((xx, yy, zz), z_wpsi)
#z_part = my_interpolating_function(x, y, z)
z = seeding_part.ini_depth(maskrho,simul,depths0,x,y,z,z_wpsi=z_wpsi)

nq = np.min([len(x.reshape(-1)),nqmx])


        ###################################################################################
''' no need for topocheck anymore as we are using sigma levels'''
''' but we have to remove pyticles which are in a masked area'''

ipmx = 0; px0,py0,pz0,ptopo0 = [],[],[],[]

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
       ptopo0.append(ptopo[ip])
       ipmx +=1
    
    #del temp,salt
nq = ipmx
def plot_some():
    plt.figure
    plt.subplot(221)
    plt.plot(px0)
    plt.subplot(222)
    plt.plot(py0)
    plt.subplot(223)
    plt.plot(pz0)
    plt.subplot(224)
    plt.plot(ptopo0)
    plt.show()
#    del x,y,z

#---------- Interpolating z_rho on p


#JC check numerical error
coord = simul.coord
i0=coord[2]; j0=coord[0];
k0 = 0
#z_wpsi[:,:,:] = z_wpsi[5,5]
ng = 0
pdepth_test = partF.interp_3d_psiw(px0,py0,pz0,z_wpsi,ng,nq,i0,j0,k0)
ptopo_test = partF.interp_2d(px0,py0,simul.topo,0,nq,i0,j0)

print(f'JC DEBUG ============')
print(f' pdepth_test = {pdepth_test}')
    
plt.figure

plt.subplot(121)
n, bins, patches = plt.hist(pdepth_test,20)
plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title('Histogram of pdepth à t=0')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
# plt.axis([.03])
plt.grid(True)

plt.subplot(122)
n, bins, patches = plt.hist(ptopo_test,20)
plt.xlabel('Smarts')
plt.title('Histogram of ptopo at t = 0')
plt.grid(True)

plt.show()

# Error

def some_plot_bis():
    plt.subplot(121)
    plot_some()
    plt.subplot(122)
    plt.plot(pdepth_test)
    return


# -- Another test
depth_again = []
for ip in np.arange(nq):
    g = interp1d(list(range(nz+1)), z_wpsi[int(px0[ip]),int( py0[ip]) ], kind = 'cubic')
    depth_again.append((np.float(g(pz0[ip]))))
#    print(ip)

plt.plot(depth_again)
plt.show()



