#!/usr/bin/env python
'''

Version 2.0 of pyticles script

- Run on sigma levels
- New Horiz Interpolation options 
- New Time Interpolation options 

designed for intensive memory usage

run the script using:
python Pyticles.py -1
or
./Pyticles.py -1


!---------------------------------------------------------------------------------------------
! 18/05/18:
!     - add online computation of lon,lat,depth,topo
! 27/04/17:
!     - Change condition for update
! 08/02/17:
!     - shift continuous release by -1 time-step
! 18/12/16:
!     - change time_int to time in output file
!     - write outputs at the beginning of the time loop and start time at time0 
!       [time_int was shifted by 1 in earlier versions]
!     - add capability of releasing particles and writing outputs at subtime-steps period
! 14/09/16:
!     - Compute maxvel using full 3d velocity
! 26/08/16:
!     - add 2d/3d advection capability
!     - add write_uv option
! 24/08/16:
!     - add restart capability
! 05/10/16:
!     - add mask as a shared array and use it in the time_step routines
!     - add the # define CHECK_MASK in pyticles_3d_sig_sa.F
! 04/19/16:
!     - add routines for old S-Coordinate system 
! 03/25/16:
!     - add capabilities of dt_injection > 1 (continuous injection of particles every dt_injection )
!     - Correct bug in postprocessing scripts for ng>0
! 02/03/16:
!     - Modify mask handling (0 instead of nan*s for u,v,w)
!     - add tests for nan in time_step for RK4 (particles escaping an open-boundary domain)
! 01/26/16:
!     - Add periodic capabilities by adding ghost points (size=ng) at horizontal boundaries
!       [if ng>0 at an open boundary, we put nan so that particles going out are removed]
!     - subtiling of the domain deactivated while using a periodic domain 
!     - add option write_t to compute and interpolate temp (ptemp) and write it in the output
!----------------------------------------------------------------------------------------------

'''
###################################################################################
# Load all useful modules
##################################################################################

#Plotting modules
import matplotlib
matplotlib.use('Agg') #Choose the backend (needed for plotting inside subprocess)
import matplotlib.pyplot as plt
import matplotlib.colors as col

#Some standard modules
import sys, os
import numpy as np
import time as tm
from netCDF4 import Dataset
import multiprocessing as mp
import ctypes 
import queue
import resource

#add the Modules folder in your python PATH
#sys.path.remove("/home2/datahome/jgula/Python_Modules") #just for JG
sys.path.append("./Modules/") 

#Specific modules needed for pyticles
import pyticles_sig_sa as part
import pyticles_3d_sig_sa as partF

#Simulations (path, data...)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy.ma as ma

#from scipy.interpolate import interp1d

#Specific modules needed for pyticles
import pyticles_sig_sa as part
import pyticles_3d_sig_sa as partF
#import seeding_part
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from R_files import load

#For nested loop
from itertools import product
from copy import *

from scipy.interpolate import interp1d
import seeding_part


##################################################################################
# Parameters for multiprocessing
##################################################################################

if len(sys.argv)<=1:
    #no proc number has been specified as argv... choose 0
    nproc=0
else:
    #define nproc as specified
    nproc = int(sys.argv[1])

print('-----------------------------------')



#Use multiprocess version of the code   
print('Parallel version')

#If nproc<0: use all available processors
if nproc<=0: 
    nproc=mp.cpu_count();
else:
    nproc = np.min([mp.cpu_count(),nproc])
print(nproc, ' processors will be used')

print('-----------------------------------')

###################################################################################
###################################################################################
# THE FOLLOWING CONTAINS DEFINITION OF THE PARTICLES SETTINGS TO BE EDITED BY USER
###################################################################################
###################################################################################


###################################################################################
# Define location of outputs
##################################################################################

# name of your configuration (used to name output files)
config='Rho_Seed'

folderout= '/home/jeremy/Bureau/Data/Pyticles/' + config + '/'

if not os.path.exists(folderout):
    os.makedirs(folderout)

###################################################################################
# Load simulations parameters (date, subsection, static fields, files path ...etc..)
# [you can check "dir(simul)" to see what has been loaded]
###################################################################################

#time-stepping
timestep='RK4' # Choices are FE (forward-Euler)
               #             RK2, RK4 (Runge-Kutta 2nd and 4th order)
               #             AB2, AB3, AB4 (Adams-Bashforth 2,3,4th order)
               #             ABM4 (Adams-Bashforth 4th order + Adams-Moulton corrector).

nsub_x, nsub_y = 1,1 #subtiling, will be updated later automatically
nadv = 1 # number of extra-points needed for interpolation, 0 for linear, 1 for higher order, etc.
debug = True

#############
x_periodic = False
y_periodic = False
ng = 1 #number of Ghostpoints _ 1 is enough for linear interp _ 2 for other interp
#############
#3D advection
adv3d = True 
# if adv3d = False then the particles are advected in 2d using horizontal velocities at advdepth
if len(sys.argv)>=3:
    advdepth = np.int(sys.argv[2]) 
else:
    advdepth = 0 
'''
        NOTE that advdepths is used as follows:

        advdepth <= 0 means depths in meters
        advdepth = 0 means surface
        advdepth > 0 means sigma-level (1 = bottom [0 in netcdf file],...
                             ..., Nz = surface [Nz-1 in netcdf file])
'''
#############
meanflow=False # if True the velocity field is not updated in time
#############    python Pyticles.py 14 $depth > output_case1
# JC modif

# Initialize seeding particles using a boolean condition ini_cond
# Inside a box_grid as usual
# Box grid supports: an horizontal (ic, jc), denisty nnx, nny
#                    minimum maximum sigma levels 
# Ex : temp < 5°C  
# Does not support vertical condition
# i.e can't state pcond = True and depth = z0
# Therefore if ini_cond = True: initial_depth = False

continuous_injection = False 
# if True release particles continuously, if False only one release at initial time-step

initial_cond = True
initial_depth = True

if initial_cond:
    initial_depth = False
   # continuous_injection = False

sedimentation=True
w_sed0 = -0 # vertical velocity for particles sedimentation (m/s)

#name of the simulation (used for naming plots and output files)
simulname = '_' + config
if (not adv3d) and (advdepth > 0):
    simulname = simulname + '_adv' + '{0:04}'.format(advdepth) + 'sig'
elif (not adv3d) and (advdepth <= 0):     
    simulname = simulname + '_adv' + '{0:04}'.format(-advdepth) + 'm'
    sedimentaion=False
    w_sed0 = 0. # JC no sedimentation for 2D advection


#Write T,S at each particle position directly in output file
write_ts=True

#Write only Temperature (for simulations with no S)
write_t=False
if write_t: write_ts = False

#Write lon,lat,topo,depth
write_lonlat=True
write_depth=True
write_topo=True

#Write u,v at each particle position directly in output file
write_uv=True
if adv3d: write_uv=False #not implemented yet for 3d


###################################################################################
# ROMS outputs
###################################################################################

# dfile is frequency for the use of the ROMS outputs (default is 1 = using all outputs files)
dfile = 1
start_file = 1550
end_file = 1560

#############

restart = False
restart_time = 4 #nb of time steps in the restart_file
restart_file = '/home/jeremy/Bureau/Data/Pyticles/home/jeremy/Bureau/Data/Pyticles/Port_Test_P3/Case_1_Port_Test_P3_12_1550.nc'

if not restart: 
    restart_time = 0
else:
    start_file += restart_time

#############

# Load simulation [mysimul is the name of the simul as defined in Modules/R_files.py]
parameters = 'Case_1 [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

'''
simul attributes are:
['Cs_r', 'Cs_w', 'Forder',  'coord', 'coordmax', 'cst', 'date', 'day', 'domain', 'dt', 'dtime', 'f', 'filetime', 'floattype', 'g', 'get_domain', 'hc', 'hour', 'infiletime', 'load_file', 'mask', 'min', 'month', 'ncfile', 'ncname', 'oceandate', 'oceantime', 'pm', 'simul.pn', 'rdrg', 'rho0', 'simul', 'time', 'time0', 'topo', 'update', 'variables_grd', 'varnames', 'x', 'y', 'year']

simul.ncname attributes are:
['Z0',  'frc', 'grd', 'his', 'tend', 'tfile', 'tstart', 'wind']

'''

simulname = simul.simul +  simulname

simul.dtime = np.sign(dfile) * np.ceil(np.abs(dfile))

# Resolution (dx,dy)
dx, dy   = 1./np.mean(simul.pm), 1./np.mean(simul.pn)

# Total size of the domain (nx,ny,nz)
(nx,ny) = simul.pm.shape

nx += 2*ng; ny += 2*ng # add ghost points to array size
depths = simul.coord[4]
nz = len(depths)
k0 = 0

mask = simul.mask
maskrho = copy(mask)
maskrho[np.isnan(maskrho)] = 0.

if not adv3d: maskrho[simul.topo<-advdepth] = 0.

topo = simul.topo
filetime = simul.filetime

timerange = np.round(np.arange(start_file,end_file,dfile),3)


#for timing purpose
tstart = tm.time()

#Time all subparts of the code 
timing=True


###################################################################################
# Define Particle seeding
###################################################################################

if True:

    #Initial Particle release
    subtstep = np.int(360 * np.abs(dfile))    # Number of time steps between frames
    nqmx = 25000   # maximum number of particles
    maxvel0 = 5    # Expected maximum velocity (will be updated after the first time step)
    
    ##########################
    # Particles initial location in the horizontal

    #if config=='TEST_Py3':
    [ic,jc] = [600,800] #= part.find_points(simul.x,simul.y,-32.28,37.30)

    # distance between 2 particles [in m]
    dx_m = 1000.
    
    dx0 = dx_m * simul.pm[ic,jc] # conversion in grid points

    iwd  = 5.* dx0 # half width of seeding patch [in grid points
    jwd  = 5.* dx0 # half width of seeding patch [in grid points]

    #########
    # density of pyticles (1 particle every n grid points)
    nnx=dx0
    nny=dx0
    nnlev=1

    ##########################
    # first and last vertical level to fill with pyticles
    lev0= 0
    lev1= len(depths)
    
    if not adv3d: lev0 = -1; lev1 = lev0

    #########
    # define initial vertical position using depth
    depths0 = [-50, -500] # [-50, -100, -200]
    if initial_depth:
        lev1 = lev0 + len(depths0) - 1
        nnlev = 1


    # OLD VERSION
    # define sigma vertical position using a condition on roms output
    # condition pcond is defined later
    # possibility to pass to coord as an argument to go faster
    if initial_cond:

        lev1 = len(depths)
        nnlev = 1

        #temp = part.get_t_io(simul, x_periodic=x_periodic,
        #        y_periodic=y_periodic, ng=ng)
        #ini_cond = (temp > 14.) & (temp < 16.)
     #   print('------------------------')
     #   print(f'ini_cond.shape {ini_cond.shape}')

###########

if continuous_injection:
    dt_injection = 1 #(1 = injection every time step, 10 = injection every 10 time steps)
    N_injection = 1 + np.int(timerange.shape[0]/dt_injection)


###################################################################################
###################################################################################
# THE FOLLOWING SHOULD NOT BE EDITED BY USER
###################################################################################
###################################################################################


def shared_array(nx,ny=1,nz=1,nt=1,prec='double',value=np.nan):
    '''
    Function used to create shared variables compatible with numpy and fortran ordered
    '''
    if prec=='float':
        shared_array_base = mp.Array(ctypes.c_float, nx*ny*nz*nt )
    elif prec=='double':
        shared_array_base = mp.Array(ctypes.c_double, int(nx*ny*nz*nt), lock=True)
    elif prec=='int':
        shared_array_base = mp.Array(ctypes.c_int32, nx*ny*nz*nt )        
    var = np.ctypeslib.as_array(shared_array_base.get_obj())
    
    if nt>1:   
        var = var.reshape(nx,ny,nz,nt, order='F');
    elif nz>1:
        var = var.reshape(nx,ny,nz, order='F');
    elif ny>1:
        var = var.reshape(nx,ny, order='F');
    else:
        var = var.reshape(-1, order='F');
        
    var[:] = value
    
    return var

###################################################################################
###################################################################################
# INITIALIZATION
###################################################################################
# Note that px,py,pz are ZERO-BASED !!!
# px,py,pz directly correspond to the indices of a python array (var[px,py,pz])
###################################################################################

if not restart:
    ###################################################################################
    # Define initial px,py,pz pyticles position (Fast .py version_ fill in order x,y,z)
    ###################################################################################

    #z,y,x = np.mgrid[lev0:lev1+1:nnlev,
    #        np.max([jc-jwd,1]):np.min([jc+jwd+np.min([1.,jwd]),ny]):nny,
    #        np.max([ic-iwd,1]):np.min([ic+iwd+np.min([1.,iwd]),nx]):nnx]

    z, y, x = seeding_part.seed_box(ic=ic, jc=jc, lev0=lev0,
            lev1=lev1, iwd=iwd, jwd=jwd, nx=nx, ny=ny, nnx=nnx, nny=nny,
            nnlev=nnlev)

   # if (z == z_mod).all():
   #     print(f'Z is ok')
   # if (y == y_mod).all():
   #     print(f'Y is ok')
   # if (x == x_mod).all():
   #     print(f'X is ok')


    if initial_depth: #initial vertical position = depths0
        z_w = part.get_depths_w(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)
        z = seeding_part.ini_depth(maskrho,simul,depths0,x,y,z,z_w,ng=ng)
    
    nq = np.min([len(x.reshape(-1)),nqmx])
    if (nq / nproc) < 10:
        print('----------------------------------------------------')
        print(f'WARNING : Multithreading Issue')
        print('number of particles too small relatively to nprocs')
        print(f'in subranges')
        print('use less procs')
        print('----------------------------------------------------')
    ###################################################################################
    ''' no need for topocheck anymore as we are using sigma levels'''
    ''' but we have to remove pyticles which are in a masked area'''
    # ipmx : count seeding partciles
    # px0, py0, pz0 : initial position for seeding
    # topolim : used in ADV_2D to prevent particles from being seeded below seafloor
    
    ipmx = 0; px0,py0,pz0 = [],[],[]
    topolim=0

    if not adv3d: topolim = np.nanmax([topolim,-advdepth])

    # initializing px0, py0, pz0
    if initial_cond:
        # only true for rho variables
        # maybe not safe nor useful to use i0, j0 ,k0 especially if only load
        # the subdomain for cond (see above)
        i0, j0, k0 = 0, 0, 0
        '''
        A matter of functionality... We take temp, salt, z_w from simulation 
        Then interpolate each varialbe on px0,py0,pz0
        Finally we compute potential density

        Maybe we could compute rho from temp, salt, z_w and finally interpolate it...
        '''
        [temp, salt] = part.get_ts_io(simul, x_periodic = x_periodic,
                y_periodic = y_periodic, ng=ng)
        
        ptemp = partF.interp_3d(x.reshape(-1), y.reshape(-1), z.reshape(-1),
                temp, ng, nq, i0, j0, k0)
        
        psalt = partF.interp_3d(x.reshape(-1), y.reshape(-1), z.reshape(-1),
                salt, ng, nq, i0, j0, k0)
        
        z_w = part.get_depths_w(simul, x_periodic = x_periodic,
                y_periodic = y_periodic, ng=ng)
        
        pdepth = partF.interp_3d_w(x.reshape(-1), y.reshape(-1), z.reshape(-1),
                z_w, ng, nq, i0, j0, k0)
        
        prho = seeding_part.prho1(ptemp = ptemp, psalt = psalt, pdepth = pdepth)
        
        rho_min = 1026
        rho_max = 1027
        pcond = (prho > rho_min) & (prho < rho_max)
        
        # WE will keep someting like this for future
        #
        #pcond = partF.interp_3d(x.reshape(-1), y.reshape(-1), z.reshape(-1),
        #        ini_cond, ng, nq, i0, j0, k0)
        
        ipmx = seeding_part.remove_mask(simul, topolim, x, y, z, px0, py0, pz0, nq,
                ng=ng, pcond=pcond)
    else:
        ipmx = seeding_part.remove_mask(simul, topolim, x, y, z, px0, py0, pz0,
                nq, ng=ng)

    if (not initial_cond) and (not continuous_injection): del x,y,z
    # Else we need the grid box to compute px0, py0, pz0 at each injection time

    nq = ipmx
    print('nq = {nq}')
    ###################################################################################

    if continuous_injection:
        '''
        we'll multiply nq by the number of injections (still in the limit of nqmx)
        and leave nan for future pyticles which will be released at the same positions every time step
        '''
        nq_injection = nq
        nq = np.nanmin([nq_injection*(N_injection+1),nqmx])
        print('it would take', nq_injection*N_injection - nqmx, ' more pyticles')
        print('to be able to release through all the simulation')
        nq_1=nq_injection
        if debug:
            print('---------------------------------------')
            print(f'nq_injection = {nq_injection}')
            print(f'nq = {nq}')


    else:
        nq_1=-1  

    ###################################################################################
    # nq: total number of particles
    
    px = shared_array(nq,prec='double')
    py = shared_array(nq,prec='double')
    pz = shared_array(nq,prec='double')

    print('px.shape = {px.shape}')
    if continuous_injection:
        px[:nq_1]=px0; py[:nq_1]=py0; pz[:nq_1]=pz0

    else: 
        px[:]=px0; py[:]=py0; pz[:]=pz0
        del px0,py0,pz0


    ###################################################################################
# restart = True
else: 

    if not continuous_injection:
        # just load px,py,pz from restart_file
        nc = Dataset(restart_file, 'r')
        px0 = nc.variables['px'][restart_time,:]
        py0 = nc.variables['py'][restart_time,:]
        pz0 = nc.variables['pz'][restart_time,:]
        nc.close()

        nq = len(px0)
        
        px = shared_array(nq,prec='double')
        py = shared_array(nq,prec='double')
        pz = shared_array(nq,prec='double')

        px[:]=px0; py[:]=py0; pz[:]=pz0
        del px0,py0,pz0
        
    else:
        # load px,py,pz from restart_file
        nc = Dataset(restart_file, 'r')
        px0 = nc.variables['px'][restart_time,:]
        py0 = nc.variables['py'][restart_time,:]
        pz0 = nc.variables['pz'][restart_time,:]

        nq = len(px0)
        
        px = shared_array(nq,prec='double')
        py = shared_array(nq,prec='double')
        pz = shared_array(nq,prec='double')

        px[:]=px0; py[:]=py0; pz[:]=pz0
        del px0,py0,pz0

        ##################################
        # determine px0,py0,pz0,nq_injection
        # JC not True in case of initital_cond and continuous injection
        # because px0, py0, pz0 may vary upon time
        nq_injection = np.argmax(np.isnan(nc.variables['px'][0,:]))

        px0 = nc.variables['px'][0,:nq_injection]
        py0 = nc.variables['py'][0,:nq_injection]
        pz0 = nc.variables['pz'][0,:nq_injection]
        nc.close()
        
        nq_1=np.nanmin([nq_injection*((restart_time)//dt_injection+1),nqmx])
        
        '''
        if (nq_1<nqmx) and (restart_time%dt_injection)==0:

            nq_0 = nq_1
            nq_1 = np.nanmin([nq_injection*(restart_time/dt_injection+1),nqmx])
            px[nq_0:nq_1]=px0[:nq_1-nq_0]; py[nq_0:nq_1]=py0[:nq_1-nq_0];
            if adv3d: pz[nq_0:nq_1]=pz0[:nq_1-nq_0]
        '''

###################################################################################

# Time between 2 frames (in seconds)
delt   = shared_array(2,value=simul.dt*np.abs(dfile)) # Expected time between frames

maxvel = shared_array(2,prec='double',value=maxvel0)

# Total number of time steps:
istep = shared_array(1,prec='int',value=-1)

# Index of the  previous (itim[0]) and next(itim[1]) time-step for u,v,w,dz
itim = shared_array(2,prec='int')
itim[:]=[0,1]

###################################################################################
#If using a Adams-Bashforth method we need to have access to previous vel. values
# So we create some shared arrays contaning dpx, dpy, dpz for the previous (ab_order-1) time-steps

if timestep[:2]=='AB':
    ab_order = np.int(timestep[-1])
    dpx = shared_array(nq,ab_order,prec='double')
    dpy = shared_array(nq,ab_order,prec='double')
    dpz = shared_array(nq,ab_order,prec='double')
    iab = shared_array(ab_order,prec='int',value=0)
    iab[:]=range(ab_order) # init values of iab



###################################################################################
# Call a function as a subprocess
###################################################################################


def run_process(my_func):
    ''' 
    We do not pass arguments nor declare any global/local variable so that 
    all variables will be accessible through the subroutine, 
    but only variables declared as global shared arrays (px,py,pz,ptemp,psalt)
    can be modified by the subroutine...
    '''
    results = mp.Queue(); i=0
    proc=mp.Process(target=my_func, args=())
    proc.start(); proc.join()
    
    
    return results


###################################################################################
# Update px,py,px
###################################################################################

def update_xyz():
     exec(compile(open('Pyticles_subroutines/update_xyz_largemem.py').read(),\
             'Pyticles_subroutines/update_xyz_largemem.py', 'exec'))

###################################################################################
# Compte T,S at each pyticles positions -> ptemp,psalt
###################################################################################

def update_ts():   
    exec(compile(open('Pyticles_subroutines/update_ts.py').read(),\
            'Pyticles_subroutines/update_ts.py', 'exec'))
    
###################################################################################
# Compte T at each pyticles positions -> ptemp
###################################################################################

def update_t():   

    exec(compile(open('Pyticles_subroutines/update_t.py').read(), \
            'Pyticles_subroutines/update_t.py', 'exec'))

###################################################################################
# Compte T at each pyticles positions -> ptemp
###################################################################################

def update_depth():   

    exec(compile(open('Pyticles_subroutines/update_depth.py').read(), \
            'Pyticles_subroutines/update_depth.py', 'exec'))

###################################################################################
# Compte T at each pyticles positions -> ptemp
###################################################################################

def update_lonlat():   
    exec(compile(open('Pyticles_subroutines/update_lonlat.py').read(), \
            'Pyticles_subroutines/update_lonlat.py', 'exec'))

###################################################################################
# Compte T at each pyticles positions -> ptemp
###################################################################################

def update_topo():   
    exec(compile(open('Pyticles_subroutines/update_topo.py').read(), \
            'Pyticles_subroutines/update_topo.py', 'exec'))

###################################################################################
# Compte u,v at each pyticles positions -> pu,pv
###################################################################################

def update_uv_2d():   
    exec(compile(open('Pyticles_subroutines/update_uv_2d.py').read(), \
            'Pyticles_subroutines/update_uv_2d.py', 'exec'))

###################################################################################
# Create output file and write time, px,py,pz ( and ptemp,psalt)
###################################################################################

def write_output():   
    exec(compile(open('Pyticles_subroutines/write_output.py').read(), \
            'Pyticles_subroutines/write_output.py', 'exec'))
    
###################################################################################

def linear(var1,var2,alpha):
    return alpha * var2 + (1.-alpha) * var1
    
###################################################################################
# Plot pyticles on a SST map (is done at each time-step)
###################################################################################

subtightcoord_saves=[]; subcoord_saves=[];  subsubrange_saves=[]


def plot_rect(rect,line='k'):
    plt.plot([rect[2],rect[2],rect[3],rect[3],rect[2]],\
             [rect[0],rect[1],rect[1],rect[0],rect[0]],line, linewidth=2)
    
###################################################################################


def plot_selection(alldomain=True):

    plt.figure(figsize=(6.0,4.0))  
    ax1 = plt.subplot(1,1,1)

    nc = Dataset(simul.ncfile, 'r')
    if alldomain:
        if adv3d or advdepth==0:
            sst = np.squeeze(simul.Forder(nc.variables['temp'][simul.infiletime,-1,:,:]))
        else:
            [z_r,z_w] = part.get_depths(simul)
            temp =  simul.Forder(np.squeeze(nc.variables['temp'][simul.infiletime,:,:,:]))
            sst = part.vinterp(temp, [advdepth], z_r,z_w)[:,:,0]
            sst[sst==0.] = np.nan
        sst *= simul.mask
        topo = simul.topo
    else:
        [ny1,ny2,nx1,nx2] = np.array(coord)-ng

        if adv3d or advdepth==0:
            sst = np.squeeze(simul.Forder(nc.variables['temp'][simul.infiletime,-1,ny1:ny2,nx1:nx2]))
        else:
            [z_r,z_w] = part.get_depths(simul,coord=[ny1,ny2,nx1,nx2])
            temp =  simul.Forder(np.squeeze(nc.variables['temp'][simul.infiletime,:,ny1:ny2,nx1:nx2]))
            sst = part.vinterp(temp, [advdepth], z_r,z_w)[:,:,0]
            sst[sst==0.] = np.nan
        sst *= simul.mask[nx1:nx2,ny1:ny2]
        topo = simul.topo[nx1:nx2,ny1:ny2]
    nc.close()

    sst[sst<0] = 0.
    #plt.imshow(sst[:,::1].T); plt.colorbar(shrink=0.25);
    plt.pcolormesh(ma.masked_invalid(sst[:,:].T),cmap='jet',rasterized=True);
    plt.colorbar(shrink=0.25);

    if not adv3d and advdepth<-topo.min():
        plt.contourf(topo.T,[0,-advdepth],cmap = col.LinearSegmentedColormap.from_list('my_colormap',['white','lightgray'],256))
        plt.contourf(topo.T,[0,topo.min()],cmap = col.LinearSegmentedColormap.from_list('my_colormap',['Gainsboro','gray'],256))
        plt.contour(topo.T,[topo.min()],colors=('k',),linewidths=(0.5,));
        plt.contour(topo.T,[-advdepth],colors=('k',),linewidths=(0.2,));

    if alldomain:
        plt.plot(px[::1]+0.5,(py[::1]+0.5),'o', markersize=2, markeredgecolor='k', markerfacecolor='k');
        plt.axis([0,sst.shape[0]-1,0,sst.shape[1]-1]);
    else:
        plt.plot(px[::1]-coord[2]+0.5,(py[::1]-coord[0]+0.5),'o', markersize=1, markerfacecolor='white');
        #plt.axis('scaled');
        #plt.axis([nx1-coord[2]+0.5,nx2-coord[2]+0.5,ny1-coord[0]+0.5,ny2-coord[0]+0.5]);
    
    '''
    for isub,jsub in product(range(nsub_x),range(nsub_y)):
        try:
            plot_rect(subtightcoord_saves[jsub+isub*nsub_y],'r')
            plot_rect(subcoord_saves[jsub+isub*nsub_y],'k--')
        except:
            print 'no subdomains to plot'
    '''

    color = 'w'; box = 'round,pad=0.1'; props = dict(boxstyle=box, fc=color, ec='k', lw=1, alpha=1.)
  #JC  ax1.text(0.95,0.05,simul.date[:-8], horizontalalignment='right', verticalalignment='bottom', bbox=props, transform=ax1.transAxes)

    plt.title(format(np.sum(px>0)) + ' pyticles ' )
    plt.savefig(folderout + simulname + '_' + format(nproc) + '_' + '{0:04}'.format(time+dfile) +'.png', size=None, figure=None, magnification='auto', dpi=150,bbox_inches='tight'); plt.clf()
    


###################################################################################

def plot_selection_sub(alldomain=True):
    
    nc = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
    if alldomain:
        sst = np.squeeze(simul.Forder(nc.variables['temp'][simul.infiletime,-1,:,:]))
        sst *= simul.mask
    else:
        sst = np.squeeze(simul.Forder(nc.variables['temp'][simul.infiletime,-1,ny1:ny2,nx1:nx2]))
        sst *= simul.mask[nx1:nx2,ny1:ny2]

    nc.close()

    for sub in range(len(subtightcoord_saves)):
        plt.imshow(sst[:,::1].T);
        #plt.contourf(sst.T,np.arange(np.floor(np.nanmin(sst)),np.ceil(np.nanmax(sst))+0.1,0.01*10));    
        plt.colorbar(shrink=0.25);
        
        plt.plot(px[subsubrange_saves[sub]]+0.5,(py[subsubrange_saves[sub]]+0.5),'o', markersize=1, markerfacecolor='white');
        plt.axis([0,sst.shape[0],0,sst.shape[1]]); 

        plot_rect(subtightcoord_saves[sub],'r')
        plot_rect(subcoord_saves[sub],'k--')
        
        plt.title(format(np.sum(px>0)) + ' pyticles ' )
        plt.savefig(folderout + simulname + '_sub_ '+ format(sub) +'_' + '{0:04}'.format(time+dfile) +'.png', size=None, figure=None, magnification='auto', dpi=150,bbox_inches='tight'); plt.clf()

###################################################################################
# Create output netcdf file for pyticles
###################################################################################

#File name
if not restart:
    newfile = folderout + simulname + '_' + format(nproc) + '_' + '{0:04}'.format(filetime) + '.nc'
else:
    newfile = restart_file

print('newfile', newfile)

###################################################################################
#START OF THE TIME LOOP
###################################################################################

print (' ')
print (' ')
tstart = tm.time()

###############################
time = timerange[0]-dfile;
coord= part.subsection(px,py,nx=nx,ny=ny,offset=50)
run_process(plot_selection)
###############################




#Initialization
pm_s=np.array([]); 

itime=restart_time

for time in timerange:

    print('time is ', time)

    alpha_time = time - np.floor(time)

    ###################################################################################
    # Define domainstimerange
    # (Find index range (zero based) in which the particles are expected to stay until the next frame)
    ###################################################################################
    if debug: print('max. vel. is ',  maxvel)

    tightcoord= part.subsection(px,py,dx,dy,maxvel*0.,delt[0],nx,ny,ng)
    coord= part.subsection(px,py,dx,dy,maxvel,delt[0],nx,ny,ng, nadv= nadv)

    nx_s,ny_s = coord[3]-coord[2], coord[1]-coord[0]
    i0=coord[2]; j0=coord[0]

    print('coord is ', coord)
    
    print('Compute coord............................', tm.time()-tstart)
    tstart = tm.time()   



    ###################################################################################
    # Compute ptemp,psalt
    
    # (we are calling a subprocess in order to keep a clean memory after computation)
    ###################################################################################


    if write_ts:
        
        ptemp = shared_array(nq,prec='double'); psalt = shared_array(nq,prec='double')
        r = run_process(update_ts)

        print('get T,S..................................', tm.time()-tstart)
        tstart = tm.time()   

    elif write_t:

        ptemp = shared_array(nq,prec='double')
        r = run_process(update_t)

        print('get T....................................', tm.time()-tstart)
        tstart = tm.time()   


    if write_uv:
        
        pu = shared_array(nq,prec='double'); pv = shared_array(nq,prec='double')
        r = run_process(update_uv_2d)

        print('get u,v..................................', tm.time()-tstart)
        tstart = tm.time()   


    if write_lonlat:

        plon = shared_array(nq,prec='double'); plat = shared_array(nq,prec='double')
        r = run_process(update_lonlat)

        print('get lon,lat..............................', tm.time()-tstart)
        tstart = tm.time()   


    if write_depth:

        pdepth = shared_array(nq,prec='double'); 
        r = run_process(update_depth)

        print('get depth................................', tm.time()-tstart)
        tstart = tm.time()   


    if write_topo:

        ptopo = shared_array(nq,prec='double');
        r = run_process(update_topo)

        print('get topo................................', tm.time()-tstart)
        tstart = tm.time()   


    ###################################################################################
    # Write in file
    ###################################################################################   

    r = run_process(write_output)

    if write_ts: del ptemp,psalt
    elif write_t: del ptemp

    if write_lonlat: del plon,plat
    if write_depth: del pdepth
    if write_topo: del ptopo

    if write_uv: del pu,pv

    print('Write in file............................', tm.time()-tstart)
    tstart = tm.time()

    if debug: print('memory usage', \
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)
    


    ###################################################################################
    # Update px,py,pz 
    # (we are calling it as a subprocess in order to keep a clean memory after computation)
    ###################################################################################

    if debug: print('memory usage',\
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)
    
    ###################################################################################
    # Devide domain into subdomains to diminish use of ram
    ###################################################################################
    #Automatic division:
    
    nsub_x = 1
    nsub_y = 1
    nsub_x = 1 + (coord[3] - coord[2]) // 1000
    nsub_y = 1 + (coord[1] - coord[0]) // 1000

    # if domain is periodic, don't divide into subdomains because code cannot handle it yet!
    if x_periodic or y_periodic: nsub_x,nsub_y = 1,1 

    #no need for 2d advection
    if not adv3d: nsub_x,nsub_y = 1,1 
    ###################################################################################
    # the subtightcoord def needs to be checked to see if ng is needed)
    subtightcoord_saves=[]; subcoord_saves=[];  subsubrange_saves=[]

    if nsub_x*nsub_y==1:
        subcoord = coord
        subsubrange = list(range(nq))
        r = run_process(update_xyz);
        
    else:
        subtightcoord = [0,0,0,0]; subcoord = [0,0,0,0]; 
        for isub,jsub in product(list(range(nsub_x)),list(range(nsub_y))):
            
            #subtightcoord will be used for test if particle should be advected or not
            subtightcoord[0] = max(0, tightcoord[0]
                    + jsub*(tightcoord[1]+1-tightcoord[0])/nsub_y)
            subtightcoord[1] = min(ny , tightcoord[0] 
                    + (jsub+1)*(tightcoord[1]+1-tightcoord[0])/nsub_y)
            subtightcoord[2] = max(0, tightcoord[2] 
                    + isub*(tightcoord[3]+1-tightcoord[2])/nsub_x)
            subtightcoord[3] = min(nx , tightcoord[2] 
                    + (isub+1)*(tightcoord[3]+1-tightcoord[2])/nsub_x)
            subtightcoord_saves.append(copy(subtightcoord))
              
            #subcoord will be used to compute u,v,w (same than coord)
            subcoord[0] = max(0, subtightcoord[0] - int(np.abs(maxvel[1]*delt[0]/dy)) )
            subcoord[1] = min(ny ,subtightcoord[1] + int(np.abs(maxvel[1]*delt[0]/dy)) )
            subcoord[2] = max(0, subtightcoord[2] - int(np.abs(maxvel[0]*delt[0]/dx)) )
            subcoord[3] = min(nx ,subtightcoord[3] + int(np.abs(maxvel[0]*delt[0]/dx)) )
            subcoord_saves.append(copy(subcoord))
            
            subsubrange=[]
            #select pyticles inside subtightcoord only
            for i in range(px.shape[0]):
                if (subtightcoord[0]<=py[i]+0.5+ng<subtightcoord[1]) and\
                (subtightcoord[2]<=px[i]+0.5+ng<subtightcoord[3]):
                    subsubrange.append(i)
            
            subsubrange_saves.append(copy(subsubrange))
                    
        # run process for each subcoord and corresponding subsubrange of points
        for isub,jsub in product(list(range(nsub_x)),list(range(nsub_y))):
            subtightcoord = subtightcoord_saves[jsub+isub*nsub_y]
            subcoord = subcoord_saves[jsub+isub*nsub_y]
            subsubrange = subsubrange_saves[jsub+isub*nsub_y]
            r = run_process(update_xyz);

    ###################################################################################

    #if not meanflow and (time+dfile)%1<np.abs(dfile)*1e-2: simul.update(np.int(np.floor(time)+simul.dtime));
    if not meanflow and ( np.round(time+dfile)-(time+dfile)<=np.abs(dfile)*1e-2) : simul.update(np.int(np.floor(time)+simul.dtime));

    print('Total computation of px,py,pz............', tm.time()-tstart)
    tstart = tm.time()
        
    if debug: print('memory usage', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)

    ###################################################################################
    # Realease new pyticles (if continuous injection)
    ###################################################################################

    if (continuous_injection) and (nq_1<nqmx) and ((itime+1)%dt_injection)==0:
        nq_0 = nq_1
        nq_1 = np.nanmin([nq_injection*((itime + 1)//dt_injection + 1), nqmx])
       # JC ======================================================================= 
       # ic = ic + 10
       # jc = jc - 50
        z, y, x = seeding_part.seed_box(ic=ic, jc=jc, lev0=lev0,
            lev1=lev1, iwd=iwd, jwd=jwd, nx=nx, ny=ny, nnx=nnx, nny=nny,
            nnlev=nnlev)
        print(f'z : {z.shape}')

        if initial_depth: #initial vertical position = depths0
            z_w = part.get_depths_w(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)
            z = seeding_part.ini_depth(maskrho,simul,depths0,x,y,z,z_w,ng=ng)
        
        print(f'z : {z.shape}')
        ipmx = 0; px0,py0,pz0 = [],[],[]

        # initializing px0, py0, pz0
        if initial_cond:
            # only true for rho variables
            # maybe not safe nor useful to use i0, j0 ,k0 especially if only load
            # the subdomain for cond (see above)
            #i0, j0, k0 = 0, 0, 0
            temp = part.get_t_io(simul, x_periodic=x_periodic,
                   y_periodic=y_periodic, ng=ng, coord=coord)
            ini_cond = (temp > 14.) & (temp < 16.)

            pcond = partF.interp_3d(x.reshape(-1), y.reshape(-1), z.reshape(-1),
                ini_cond, ng, nq, i0, j0, k0)
            ipmx = seeding_part.remove_mask(simul, topolim, x, y, z, px0, py0, pz0, nq,
                    ng=ng, pcond=pcond)
            nq_1 = np.nanmin([nq_0 + ipmx, nqmx]) 
        else:
            ipmx = seeding_part.remove_mask(simul, topolim, x, y, z, px0, py0, pz0,
                    nq, ng=ng)

        del x,y,z
        if debug:
            print('-----------------------')
            print(f' i0, j0, k0 = {i0}, {j0}, {k0}')
            print(f'temp.shape = {temp.shape}')
            print(f'ini_cond.shape : {ini_cond.shape}')
            print(f'pcond.shape : {pcond.shape}')
            print(f'ipmx = {ipmx}')
            print(f'nq_0 = {nq_0}; nq_1 = {nq_1}')
            print(f'px0.shape = {np.shape(px0)}')
            print('-----------------------')
        #JC ==========================================================================
        px[nq_0:nq_1] = px0[:nq_1-nq_0]
        py[nq_0:nq_1] = py0[:nq_1-nq_0]
 
        if adv3d: pz[nq_0:nq_1] = pz0[:nq_1-nq_0]
        


    ###################################################################################
    # Plot particles position (+ SST)
        
    if (time+dfile)%1<np.abs(dfile)*1e-2: run_process(plot_selection)
    
    ###################################################################################
    itime += 1
    ###################################################################################

    #wait = raw_input("PRESS ENTER TO CONTINUE.")

    print('Plot selection...........................', tm.time()-tstart)
    tstart = tm.time()    
        
    print(' ')
    print(' ')
        
    if debug: print('memory usage', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)
    if debug: print('error',np.sqrt((px[0]-31.)**2+(py[0]-60.)**2))

        
        
        
        
    
    
    
        
        
