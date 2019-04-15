'''
Pyticles parameters
'''

debug = True # Increase verbosity to help debug

##############################################################################
# Particles Dynamcis
##############################################################################
#3D advection
adv3d = True
# if adv3d = False then the particles are advected in 2d using horizontal velocities at advdepth

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
#Write only Temperature (for simulations with no S)
write_t = False
if write_t: write_ts = False

# name of your configuration (used to name output files)
config='Iso_surf'
folderout= '/home/jeremy/Bureau/Data/Pyticles/' + config + '/'

############# Restart from a Pyticles output file
restart = False
restart_time = 4 #nb of time steps in the restart_file
restart_file = '/home/jeremy/Bureau/Data/Pyticles/home/jeremy/Bureau/Data/Pyticles/Port_Test_P3/Case_1_Port_Test_P3_12_1550.nc'

################################################################################
# Define Particle seeding
################################################################################

#Initial Particle release
nqmx = 25000   # maximum number of particles
maxvel0 = 5    # Expected maximum velocity (will be updated after the first time step)

# patch's center (in continuous injection case : user may vary its center 
# Directly in Pyticles.py) 
[ic,jc] = [600,800] #= part.find_points(simul.x,simul.y,-32.28,37.30)

dx_m = 1000. # distance between 2 particles [in m]
dx0 = dx_m * simul.pm[ic,jc] # conversion in grid points
iwd  = 5.* dx0 # half width of seeding patch [in grid points
jwd  = 5.* dx0 # half width of seeding patch [in grid points]

#########
# density of pyticles (n*dx0: particle every n grid points)
# 
nnx = 1 * dx0
nny = 1 * dx0
nnlev = 1

#########
# define initial vertical position using:
# - depth if initial_depth = True
# - if initial_cond you have to define your condition manually in Pyticles.py
# - if initial_surf (seeding particles on isosurface of a variable) 

initial_cond = False
initial_depth = False
initial_surf = True

depths0 = [-50, -500]
surf0 = [0]

# if True release particles continuously
# if False only one release at initial time-step
continuous_injection = False

##############################################################################
# Pyticles numerical schemes
# 
#time-stepping Default is RK4
timestep='RK4' # Choices are 
               # FE (forward-Euler)
               # RK2, RK4 (Runge-Kutta 2nd and 4th order)
               # AB2, AB3, AB4 (Adams-Bashforth 2,3,4th order)
               # ABM4 (Adams-Bashforth 4th order + Adams-Moulton corrector).

nsub_steps = 360 # Number of time steps between 2 roms time steps

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
start_file = 1550
end_file = 1559

# Number of time steps between frames
# This should not be here anymore...
subtstep = np.int(nsub_steps * np.abs(dfile))  

# Load simulation
# [mysimul is the name of the simul as defined in Modules/R_files.py]
# parameters = mysimul + [0,nx,0,ny,[1,nz,1]] ; nx, ny, nz Roms domain's shape 
mysimul = 'Case_1'
parameters = mysimul + ' [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

'''
simul attributes are:
['Cs_r', 'Cs_w', 'Forder',  'coord', 'coordmax', 'cst', 'date', 'day', 'domain', 'dt', 'dtime', 'f', 'filetime', 'floattype', 'g', 'get_domain', 'hc', 'hour', 'infiletime', 'load_file', 'mask', 'min', 'month', 'ncfile', 'ncname', 'oceandate', 'oceantime', 'pm', 'simul.pn', 'rdrg', 'rho0', 'simul', 'time', 'time0', 'topo', 'update', 'variables_grd', 'varnames', 'x', 'y', 'year']

simul.ncname attributes are:
['Z0',  'frc', 'grd', 'his', 'tend', 'tfile', 'tstart', 'wind']

'''




