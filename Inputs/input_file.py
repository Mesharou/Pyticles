'''
###############################################################################
###############################################################################
THE FOLLOWING CONTAINS ROMS AND PARTICLES SETTINGS TO BE EDITED BY USER
###############################################################################
###############################################################################

'''
##############################################################################
# Import some libraries
##############################################################################
import sys
import os
import numpy as np
from copy import *
import time as tm

sys.path.append("../Modules/")
from R_files import load

##############################################################################

debug = True # Increase verbosity to help debug

################################################################################
# ROMS INPUTS
################################################################################
# if meanflow = True Roms data are not updated (used for climatology)
meanflow = False
# in case of periodic channel
x_periodic = False
y_periodic = False
ng = 1 #number of Ghostpoints _ 1 is enough for linear interp _ 2 for other interp

# dfile is frequency for the use of the ROMS outputs
# (default is 1 = using all outputs files)
dfile = 1
start_file = 1510
end_file = 1512

######
# only if part_trap=True, time index in trap_file to start backward simulation
# itime_trap = -1 : last time index in forward simulation
itime_trap = -11 
trap_file = '/home/jeremy/Bureau/Data/Pyticles/Trap_fwd/' \
                    + 'Case_1_Trap_fwd_adv200.0m_6_1510.nc'

###### Restart from a Pyticles output file
# user should not change start_file
# restart_time : number of time step since start_file
restart = False
restart_time = 7 #nb of time steps in the restart_file
restart_file = '/home/jeremy/Bureau/Data/Pyticles/' \
               +'/Cubic_adv/Case_1_Cubic_adv_4_1510.nc'

if not restart:
    restart_time = 0
else:
    start_file += restart_time

# Load simulation
# parameters = my_simul + [0,nx,0,ny,[1,nz,1]] ; nx, ny, nz Roms domain's shape 
my_simul = 'Case_1'
parameters = my_simul + ' [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

##############################################################################
# Pyticles numerical schemes (TO BE EDITED)
#
#time-stepping Default is RK4
timestep = 'RK4' # Choices are
               # FE (forward-Euler)
               # RK2, RK4 (Runge-Kutta 2nd and 4th order)
               # AB2, AB3, AB4 (Adams-Bashforth 2,3,4th order)
               # ABM4 (Adams-Bashforth 4th order + Adams-Moulton corrector).

nsub_steps = 360 # Number of time steps between 2 roms time steps

# Spatial interpolation
# Default is linear
# Available : #define CUBIC_INTERPOLATION
#             #define CRSPL_INTERPOLATION
#             #define WENO_INTERPOLATION
# Beware these higher order schemes have not been rigorously tested
# To define them, in Modules/interp_3d_for_pyticles.F
# Activate ccp keys : NEW_VERSION and chosen numerical scheme
# Compile cpp keys use make command

nadv = 1 # deprecated


##############################################################################
# Particles Dynamics
##############################################################################
# 3D advection
adv3d = True
advzavg = False

if advzavg:
    z_thick = 100. # water column thickness to average 2D velocity field around
                   # Around advdepth
# Else 2D advection using (u,v) interpolated at advdepth 
if not adv3d:
    advdepth = -200.

'''
        NOTE that advdepths is used as follows:

        advdepth <= 0 means depths in meters
        advdepth = 0 means surface
        advdepth > 0 means sigma-level (1 = bottom [0 in netcdf file],...
                             ..., Nz = surface [Nz-1 in netcdf file])
'''
# sedimentation of denser particles (not supported in 2D case)
sedimentation = False
w_sed0 = -40 # vertical velocity for particles sedimentation (m/s)

if not adv3d:
    sedimentation = False
    w_sed0 = 0. # JC no sedimentation for 2D advection

##############################################################################
# Pyticles Outputs
##############################################################################

#Write lon,lat,topo,depth
write_lonlat = True
write_depth = True
write_topo = True
if advzavg: 
    write_topo = True # Needed to keep track when water column intersects with
                      # bathymetry (topo > |advdepth| - z_thick/2)
write_uv = True
write_ts = True
write_uvw = True
if write_uvw:
    write_uv = False

#Write only Temperature (for simulations with no S)
write_t = False
if write_t: write_ts = False

# name of your configuration (used to name output files)
config = 'tmp_debug_test'
folderout = '/home/jeremy/Bureau/Data/Pyticles/' + config + '/'
# create folder if does not exist
if not os.path.exists(folderout):
    os.makedirs(folderout)

#################################################################
# This section should not be edited by users
#################################################################
#name of the simulation (used for naming plots and output files)
simulname = '_' + config
if (not adv3d) and (advdepth > 0):
    simulname = simulname + '_adv' + '{0:04}'.format(advdepth) + 'sig'
elif (not adv3d) and (advdepth <= 0):
    simulname = simulname + '_adv' + '{0:04}'.format(-advdepth) + 'm'
    write_depth = False

####
simulname = simul.simul +  simulname
simul.dtime = np.sign(dfile) * np.ceil(np.abs(dfile))
# Resolution (dx,dy)
dx, dy   = 1./np.mean(simul.pm), 1./np.mean(simul.pn)
# Total size of the domain (nx,ny,nz)
(nx, ny) = simul.pm.shape
nx += 2*ng; ny += 2*ng # add ghost points to array size
depths = simul.coord[4]
nz = len(depths)
k0 = 0
mask = simul.mask
maskrho = np.copy(mask)
maskrho[np.isnan(maskrho)] = 0.
nsub_x, nsub_y = 1,1 #subtiling, will be updated later automatically

if not adv3d: maskrho[simul.topo<-advdepth] = 0.

topo = simul.topo
filetime = simul.filetime
timerange = np.round(np.arange(start_file,end_file,dfile),3)
#for timing purpose
tstart = tm.time()
#Time all subparts of the code 
timing = True
subtstep = np.int(nsub_steps * np.abs(dfile))

################################################################################
# Define Particle seeding (to be edited)
################################################################################

#Initial Particle release
nqmx = 100000  # maximum number of particles
maxvel0 = 5    # Expected maximum velocity (will be updated after the first time step)

###########
# Patch's center in grid points 
# (if continuous injection: user may vary its center Directly in Pyticles.py) 
[ic, jc] = [361, 168] #= part.find_points(simul.x,simul.y,-32.28,37.30)
barycentric = False  # Automatically modifies patch's center to previsously seeded
                    # Particles After being advected over one time step 

dx_m = 2000. # distance between 2 particles [in m]
dx0 = dx_m * simul.pm[ic,jc] # conversion in grid points
iwd  = 10* dx0 # half width of seeding patch [in grid points
jwd  = 10* dx0 # half width of seeding patch [in grid points]

#########
# density of pyticles (n*dx0: particle every n grid points)
# 
nnx = 2 * dx0
nny = 2 * dx0
nnlev = 1

#########
# define initial vertical position using:
# - depth if initial_depth = True
# - if initial_cond you have to define your condition directly in Pyticles.py
# - if initial_surf (seeding particles on isosurface of a variable) 
#   Typically on isopycnal surface, define isopycnal initial value rho0
#   surface is retrieved using potential density anomy with rho1_eos from ROMS 

# Initialize seeding particles using a boolean condition ini_cond
# Inside a box_grid as usual
# Box grid supports: an horizontal (ic, jc), denisty nnx, nny
#                    minimum maximum sigma levels 
# Ex : temp < 5°C  
# Does not support vertical condition
# i.e can't state pcond = True and depth = z0
# Therefore if ini_cond = True: initial_depth = False

initial_cond = False # 1036 in Pyticles.py
initial_depth = True
initial_surf = False

# start from a Pyticles netcdf file, with forward advection
# runs backward 3D advection from here
# Option to choose are left to user
part_trap = False 

if initial_cond:
   initial_depth = False

depths0 = [-100]
rho0 = [-1.5]

# if True release particles continuously
# if False only one release at initial time-step
continuous_injection = True
if continuous_injection:
    dt_injection = 1 #(1 = injection every time step,
                     # 10 = injection every 10 time steps)
    N_injection = 1 + np.int(timerange.shape[0] / dt_injection)

#########################################
# NOT TO BE EDITED
#########################################
# bottom to top vertical levels in sigma coordinate
lev0= 0
lev1= len(depths)

##########
# 2D advection at advdepth 
if not adv3d:
    initial_depth = True
    lev0 = -1
    lev1 = lev0
    depths0 = [advdepth]
    write_uvw = False
    write_uv = True

#########
# define initial vertical position using depth
if initial_depth:
    lev1 = lev0 + len(depths0) - 1
    nnlev = 1

#########
# boolean matrix condition to define seeding patch
if initial_cond:
    lev1 = len(depths)
    nnlev = 1




