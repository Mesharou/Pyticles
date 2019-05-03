'''
################################################################################
################################################################################
THE FOLLOWING CONTAINS DEFINITION OF THE PARTICLES SETTINGS TO BE EDITED BY USER
################################################################################
################################################################################

'''
##############################################################################
# Import some libraries
##############################################################################
import sys
import numpy as np

sys.path.append("../Modules/")
from R_files import load


##############################################################################

debug = True # Increase verbosity to help debug


################################################################################
# ROMS outputs
################################################################################

# if meanflow = True Roms data are not updated (used for climatology)
meanflow = False
# in case of periodic channel
x_periodic = False
y_periodic = False

# dfile is frequency for the use of the ROMS outputs
# (default is 1 = using all outputs files)
dfile = 1
start_file = 1510
end_file = 1535

###### Restart from a Pyticles output file
# user should not change start_file
# restart_time : number of time step since start_file
restart = False
restart_time = 25 #nb of time steps in the restart_file
restart_file = '/home/jeremy/Bureau/Data/Pyticles/' \
               +'/Visual_test/Case_1_Visual_test_12_1510.nc'

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
# Particles Dynamcis
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
sedimentation = True
w_sed0 = -0 # vertical velocity for particles sedimentation (m/s)

##############################################################################
# Pyticles Outputs
##############################################################################

#Write lon,lat,topo,depth
write_lonlat = True
write_depth = True
write_topo = True
write_uv = True
write_ts = True
write_uvw = True
if write_uvw:
    write_uv = False

#Write only Temperature (for simulations with no S)
write_t = False
if write_t: write_ts = False

# name of your configuration (used to name output files)
config = 'Visual_2_depths'
folderout = '/home/jeremy/Bureau/Data/Pyticles/' + config + '/'


################################################################################
# Define Particle seeding
################################################################################

#Initial Particle release
nqmx = 25000   # maximum number of particles
maxvel0 = 5    # Expected maximum velocity (will be updated after the first time step)

###########
# Patch's center in grid points 
# (if continuous injection: user may vary its center Directly in Pyticles.py) 
[ic,jc] = [120,400] #= part.find_points(simul.x,simul.y,-32.28,37.30)
barycentric = False  # Automatically modifies patch's center to previsously seeded
                    # Particles After being advected over one time step 

dx_m = 1000. # distance between 2 particles [in m]
dx0 = dx_m * simul.pm[ic,jc] # conversion in grid points
iwd  = 100.* dx0 # half width of seeding patch [in grid points
jwd  = 100.* dx0 # half width of seeding patch [in grid points]

#########
# density of pyticles (n*dx0: particle every n grid points)
# 
nnx = 20 * dx0
nny = 50 * dx0
nnlev = 1

#########
# define initial vertical position using:
# - depth if initial_depth = True
# - if initial_cond you have to define your condition manually in Pyticles.py
# - if initial_surf (seeding particles on isosurface of a variable) 


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
eddy_center = False

eddy_file = ''

if initial_cond:
   initial_depth = False

depths0 = [-50, -500]
surf0 = [1028]

# if True release particles continuously
# if False only one release at initial time-step
continuous_injection = False

##############################################################################
# Pyticles numerical schemes
# 
#time-stepping Default is RK4
timestep = 'RK4' # Choices are 
               # FE (forward-Euler)
               # RK2, RK4 (Runge-Kutta 2nd and 4th order)
               # AB2, AB3, AB4 (Adams-Bashforth 2,3,4th order)
               # ABM4 (Adams-Bashforth 4th order + Adams-Moulton corrector).

nsub_steps = 360 # Number of time steps between 2 roms time steps

##############################################################################



