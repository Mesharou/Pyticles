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
import pyticles_sig_sa as part

##############################################################################

debug = True # Increase verbosity to help debug

################################################################################
# ROMS INPUTS
################################################################################
# if analytical = True, we use analytical velocities, not CROCO data
analytical = True

# if meanflow = True Roms data are not updated (used for climatology)
meanflow = True
# in case of periodic channel
x_periodic = True
y_periodic = True
ng = 1 #number of Ghostpoints _ 1 is enough for linear interp _ 2 for other interp



##############################################################################
# Particles Dynamics
##############################################################################
# 3D advection
adv3d = False # 3d or 2d

'''
for 2d advection, there are several options:
any arbitrary 2d surface, vertically integrated velocity fields, etc.
- adv3d==False and advsurf==True, means that the velocity files contain only 
  one vertical level
- adv3d==False and advsurf==False, means that the velocity files are 3-d;
  the choice of the velocity field
  will depend on the value of advdepth (depths or sigma-level)
'''

if not adv3d:
    advsurf = True
    advzavg = False
else:
    advsurf = False
    advzavg = False

if advzavg:
    z_thick = 100. # water column thickness to average 2D velocity field around
                   # Around advdepth

# Else 2D advection using (u,v) interpolated at advdepth 
if not adv3d:
    advdepth = -4000.
    if advsurf:  advdepth = 0

'''
        NOTE that advdepths is used as follows:

        advdepth <= 0 means depths in meters
        advdepth = 0 means surface
        advdepth > 0 means sigma-level (1 = bottom [0 in netcdf file],...
                             ..., Nz = surface [Nz-1 in netcdf file])
'''

################################################################################
"""
- dfile: Pyticles output time step in units of CROCO outputs time step. 
  - dfile = 1, Pyticles output frequency = CROCO output frequency
  - dfile = 1/2, Pyticles output frequency = 2 x CROCO output frequency
  - dfile = -1, backward experiment

-  integers, start_file and end_file: CROCO output time step for initial and
   last Pyticles time step.

- (optional) boolean, restart: If True, restart a Pyticles simulation (typically to re-run
  a simulation that crashed)

- (optional) integer > 0, restart_time: number of CROCO time steps since initialization
  used to restart. Note that if simulation backward the sign of dfile is 
  automatically applied.

- (optional) itime_trap and trap_file: use to initialize a Pyticles simulation from a 
  former Pyticles simulation. Very specific

- str, my_simul: name of CROCO simulation specified in R_files.py to retrieve
 CROCO files

"""

dfile = 1
start_file = 0
end_file = 1000

######
# only if part_trap=True, time index in trap_file to start backward simulation
# itime_trap = -1 : last time index in forward simulation
itime_trap = -11 
trap_file = '/home/wang/Bureau/Data/Pyticles/Trap_fwd/' \
                    + 'Case_1_Trap_fwd_adv200.0m_6_1510.nc'

###### Restart from a Pyticles output file
# user should not change start_file
# restart_time : number of time step since start_file
restart = False
restart_time = 28 #nb of time steps in the restart_file
restart_file = '/home2/datawork/lwang/IDYPOP/Data/Pyticles/debug_high_freq/' \
               + 'apero_hfo3h_bk3d_06winter_trap1000m_sed50_28_3740.nc' 
              

if not restart:
    restart_time = 0
else:
    start_file += restart_time * int(np.sign(dfile))



if not analytical:
    
    # Load simulation
    # parameters = my_simul + [0,nx,0,ny,[1,nz,1]] ; nx, ny, nz Roms domain's shape 
    my_simul = 'polgyr'

    ##########
    if 'surf' in my_simul or advsurf: 
        advsurf = True
        light = True # do not load unnecessary files for a pure surface advection
    else:
        light = False
    #########

    # user may add my_simul in Module/R_files.py to indicate roms output path and
    # parameters
    parameters = my_simul + ' [0,15000,0,15000,[1,300,1]] '+ format(start_file)
    simul = load(simul = parameters, light = light, floattype=np.float64)
    


else:
    
    #domain
    dx, dy  = 1., 1.
    nxi,nyi = 100,100 # interior points, not counting ghost points
    nz = 1; dz = 1
    
    parameters, simul = part.ana_simul(nxi,nyi,nz,dx,dy,dz)
    
    #Velocity field:
    flow = [0,1,0,0] # [div,rot,S1,S2]
    
    light = True


# Names of horizontal velocity fields in the simulation files
u_name = 'u'; v_name = 'v'


##############################################################################
# Pyticles numerical schemes (TO BE EDITED)
#
#time-stepping Default is RK4
timestep = 'RK4' # Choices are
               # FE (forward-Euler)
               # RK2, RK4 (Runge-Kutta 2nd and 4th order)
               # AB2, AB3, AB4 (Adams-Bashforth 2,3,4th order)
               # ABM4 (Adams-Bashforth 4th order + Adams-Moulton corrector).

####
# -- Number of time steps between 2 roms time steps determined by CFL condition 
# sum(u_i/dx_i) < cmax

inline_cfl = False
dzmin = 1
cmax = 1

# substeps computed at each pyticles time-step
if inline_cfl:
    umax = None
    vmax = None
    wmax = None
# substeps computed at the beginning
else:
    umax = 2
    vmax = 2
    wmax = 2*1e-3
    nsub_steps = part.get_nsub_steps(simul=simul, cmax=cmax, umax=umax,
                                     vmax=vmax, wmax=wmax, dzmin=dzmin) 

####
"""
# Spatial interpolation
Default is linear
Available : #define CUBIC_INTERPOLATION
             #define CRSPL_INTERPOLATION
             #define WENO_INTERPOLATION

Beware these higher order schemes have not been rigorously tested
To define them, in Modules/interp_3d_for_pyticles.F
Activate ccp keys : NEW_VERSION and chosen numerical scheme
Compile cpp keys use make command
"""

nadv = 1 # deprecated


##############################################################################

# sedimentation of denser particles (not supported in 2D case)
sedimentation = False
sedimentation_only = False
w_sed0 = -10 # vertical velocity for particles sedimentation (m/d)

if sedimentation_only:
    sedimentation = False

if not adv3d:
    sedimentation = False
    w_sed0 = 0. # JC no sedimentation for 2D advection

# Remove particles below/above (below=True/False) a certain sigma level (klim):
remove = False
if  remove:
    below = True
    nz = len(simul.coord[4])
    klim = 75
    print(f"klim is {klim}")
else:
    klim = -1 # default value
    below = False # default value
    
##############################################################################
# Pyticles Outputs
##############################################################################
plot_part = True # Plot figures showing particles positions during execution

#Write lon,lat,topo,depth
write_lonlat = False

if adv3d:
    write_depth = True
else:
    write_depth = False
    

write_topo = False
if advzavg: 
    write_topo = True # Needed to keep track when water column intersects with
                      # bathymetry (topo > |advdepth| - z_thick/2)
elif light: 
    write_topo = False
    

write_uv = False
write_ts = False
write_uvw = False
# True : pw is w-vertical velocity in z-coordinates
# False : pw is omega velocity in sigma-coordinates
cartesian = True

if write_uvw:
    write_uv = False

if adv3d and write_uv:
    write_uv = False
    write_uvw = True    
    
#Write only Temperature (for simulations with no S)
write_t = False

if write_t: write_ts = False

# name of your configuration (used to name output files)
config = 'birotexp'

folderout = './out/'
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
(nxi, nyi) = simul.pm.shape
nx = nxi + 2*ng; ny = nyi + 2*ng # add ghost points to array size
depths = simul.coord[4]
nz = len(depths)
k0 = 0
maskrho = simul.mask
maskrho[np.isnan(maskrho)] = 0.
nsub_x, nsub_y = 1,1 #subtiling, will be updated later automatically

if not adv3d and not advsurf: maskrho[simul.topo<-advdepth] = 0.

if adv3d: advtopo = simul.topo
filetime = simul.filetime
timerange = np.round(np.arange(start_file, end_file + dfile, dfile),3)
#for timing purpose
tstart = tm.time()
#Time all subparts of the code 
timing = True

# number of sub-time steps for advection
if not inline_cfl:
    subtstep = int(nsub_steps * np.abs(dfile))
else:
    subtstep = 0

################################################################################
# Define Particle seeding (to be edited)
################################################################################

#Initial Particle release
nqmx = 100000000 # maximum number of particles
maxvel0 = 5    # Expected maximum velocity (will be updated after the first time step)

###########
# Patch's center in grid points 
# (if continuous injection: user may vary its center Directly in Pyticles.py) 
#[ic, jc] = [1450, 930] #= part.find_points(simul.x,simul.y,-32.28,37.30)

[ic, jc] = [50, 50]
print('ic,jc is ',ic,jc)

barycentric = False  # Automatically modifies patch's center to previsously seeded
                     # Particles After being advected over one time step 
    
# Size of the patch and distance between particles in meters are conserved
# even when box's center moves during simulation
preserved_meter = False

# --> injection on lon lat
spheric_injection = False

if preserved_meter:
    dx_box = 1000  # horizontal particles spacing meters
    nx_box = 100 # number of intervals in x-dir
    ny_box = 100      
    nnlev = 1  
else:
    #dx_m = 1000. # distance between 2 particles [in m]
    #dx0 = dx_m * simul.pm[ic,jc] # conversion in grid points
    dx0 = 1
    iwd  = 40 * dx0 # half width of seeding patch [in grid points
    jwd  = 40 * dx0 # half width of seeding patch [in grid points]
    # density of pyticles (n*dx0: particle every n grid points)
    nnx = 0.1 * dx0
    nny = 0.1 * dx0
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
initial_iso = False

# 2D advection
if not adv3d:
    initial_iso = False
    initial_depth = False

# start from a Pyticles netcdf file, with forward advection
# runs backward 3D advection from here
# Option to choose are left to user
part_trap = False 

if initial_cond:
   initial_depth = False


# depths0 < 0 : depth in meter
# depths0 > 0 : sigma layer
# depths0 = 0 : release at surface
depths0 = [-100]
rho0 = [-1.5]

# if True release particles continuously
# if False only one release at initial time-step
continuous_injection = False
if continuous_injection:
    dt_injection = 1 #(1 = injection every time step,
                     # 10 = injection every 10 time steps)
    N_injection = 1 + int(timerange.shape[0] / dt_injection)

    # --> compute only once initial particle position to save time
    static_injection = False

#########################################
# NOT TO BE EDITED
#########################################
# bottom to top vertical levels in sigma coordinate
lev0 = len(depths)
lev1 = len(depths)

##########
# 2D advection at advdepth 
if not adv3d:
    lev0 = -1
    lev1 = lev0
    depths0 = [advdepth]
    write_uvw = False
    write_uv = False

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




