#!/usr/bin/env python3
'''

Version 3.0 of pyticles script

- Run on sigma levels
- New Horiz Interpolation options 
- New Time Interpolation options 

designed for intensive memory usage

run the script using:
python Pyticles.py -1
or
./Pyticles.py -1

or in interactive 
python -i Pyticles.py nprocs
nprocs : numbers of processors 

! Bugs to be fixed:
    Seems to be an issue regarding the number of released particles at the last
    time step of backward simulation in case of continuous seeding of particles
    Bug was detetcted at 1020 time in Lu simulation

!---------------------------------------------------------------------------------------------
! 2019-11-29 [jeremy collin]:
!     - add key preserved_meter (input_file)
!       horizontal spacing between particles is defined with dx (m)
!       dx is preserved as the patch center moves (barycentric,...)
!       dx[ic,jc] using linear interpolation allows float coordinates
! 2019-06-12 [jeremy collin]:
!     - add capacity to seed particles from pyticles netcdf output file trap_file
!       starting from time index itime_trap
!       Purpose is to run a forward 2D simulation with output: trap_file
!       Then a backward 3D simulation to identify origin of sediment found in
!       trap_file
! 2019-05-27 [jeremy collin]:
!     - add correction for 3D xy_periodic (2D still unstable)
! 2019-05-23 [jeremy collin]:
!     - add possibility for vertically averaged zonal velocity for advection
!       (advzavg) pu being <u>(z=advdepth +/- z_thick)
!     - add possibility for bathymetric intersection in advzavg
! 2019-04-24 [jeremy collin]:
!     - add possibility for seeded particles to be released at barycenter of
!       previously released particiles after being advected for one time_step
!       barycentric = True in input_file
! 2019-04-24 [jeremy collin]:
!     - add possibilty to release particles on isopycanl surface
!       More generally any iso-surface since varaible is defined at-rho points
!       initial_surf = True in input_file
!       USER have to define rho-matrix for surface in Pyticles.py (for now...)
! 2019-04-15 [jeremy collin]:
!     - add input_file: contains all Pyticles parameters
! 2019-04-12 [jeremy collin]:
!     - converted R_tools in python 3
! 2019-04-09 [jeremy collin]:
!     - Writes more output for reusability
! 2019-04-01 [jeremy collin]:
!     - fixed initial depth numerical depth interpolation issue
!       seeding_part.ini_depth
! 2019-03-14 [jeremy collin]:
!     - New module for particle seeding
!       seeding_part.py
! 2019-03-25 [jeremy collin]:
!     - writes pw (write_uvw = True)
! 2019-03-06 [jeremy collin]:
!     - converted Pyticles project to Python 3
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
import sys, os
import numpy as np
import numpy.ma as ma
import time as tm
import multiprocessing as mp
import ctypes 
import queue
import resource

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#For nested loop
from itertools import product
from copy import *
from netCDF4 import Dataset
#from scipy.interpolate import interp1d

# Specific modules needed for pyticles
# add the Modules folder in your python PATH
# sys.path.remove("/home2/datahome/jgula/Python_Modules") #just for JG
sys.path.append("./Modules/") 
sys.path.append("./Inputs/")

# Specific modules needed for pyticles
import pyticles_sig_sa as part
import pyticles_3d_sig_sa as partF

from R_files import load
#For nested loop
from itertools import product

import seeding_part
#from R_tools_fort import rho1_eos # JC 
# Input file with Pyticles parameters
from input_file import *

################################################################################
# Parameters for multiprocessing
################################################################################

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

################################################################################
################################################################################
# THE FOLLOWING CONTAINS DEFINITION OF THE PARTICLES SETTINGS TO BE EDITED BY USER
# MOST OF PARAMTERS ARE DEFINED IN INPUTS/INPUT_FILE
################################################################################
################################################################################


################################################################################
# Load simulations parameters (date, subsection, static fields, files path ...etc..)
# [you can check "dir(simul)" to see what has been loaded]
################################################################################
#'nsub_x, nsub_y = 1,1 #subtiling, will be updated later automatically'

################################################################################
################################################################################
# THE FOLLOWING SHOULD NOT BE EDITED BY USER
################################################################################
################################################################################
# PRINTING SIMULATION PARAMETERS
print('=' * 60)
print(f'MEANFLOW : {meanflow} ')
print(f'x_periodic : {x_periodic}')
print(f'y_periodic : {y_periodic}')
print(f'dfile : {dfile}')
print(f'adv3d : {adv3d}')
print(f'advzavg : {advzavg}')
if not adv3d: print(f' advdepth: {advdepth}')
if advzavg: print(f'z_thick : {z_thick}')
if sedimentation: print(f'w_sed0 : {w_sed0} (m/s)')
print(f'initial_cond : {initial_cond}')
if initial_surf: print(f'initial surface :  {rho0}')
if initial_depth: print(f'initial depth :  {depths0}')


##############################################################################

def shared_array(nx, ny=1, nz=1, nt=1, prec='double', value=np.nan):
    '''
    Function used to create shared variables compatible with numpy and fortran
    ordered
    '''
    if prec=='float':
        shared_array_base = mp.Array(ctypes.c_float, nx*ny*nz*nt )
    elif prec=='double':
        shared_array_base = mp.Array(ctypes.c_double, int(nx*ny*nz*nt), lock=True)
    elif prec=='int':
        shared_array_base = mp.Array(ctypes.c_int32, nx*ny*nz*nt )        
    var = np.ctypeslib.as_array(shared_array_base.get_obj())
    
    if nt>1:
        try:
            var = var.reshape(nx,ny,nz,nt, order='F');
        except TypeError:
            print('------------------------------------')
            print('Type Error')
            print(f' nx = {nx} {type(nx)}')
            print(f' ny = {ny} {type(ny)}')
            print(f' nz = {nz} {type(nz)}')
            print(f' nt = {nt} {type(nt)}')

    elif nz>1:
        var = var.reshape(nx,ny,nz, order='F');
    elif ny>1:
        var = var.reshape(nx,ny, order='F');
    else:
        var = var.reshape(-1, order='F');
        
    var[:] = value
    
    return var

################################################################################
################################################################################
# INITIALIZATION
################################################################################
# Note that px,py,pz are ZERO-BASED !!!
# px,py,pz directly correspond to the indices of a python array (var[px,py,pz])
################################################################################

if not restart:
    ###########################################################################
    # Define initial px,py,pz pyticles position 
    # (Fast .py version_ fill in order x,y,z)
    ###########################################################################
    if preserved_meter:
        z, y, x = seeding_part.seed_meter(ic=ic, jc=jc, lev0=lev0, lev1=lev1,
                    nnlev=nnlev, nx_box=nx_box, ny_box=ny_box, dx_box=dx_box,
                    simul=simul, ng=ng)
    else:
        z, y, x = seeding_part.seed_box(ic=ic, jc=jc, lev0=lev0,
                lev1=lev1, iwd=iwd, jwd=jwd, nx=nx, ny=ny, nnx=nnx, nny=nny,
                nnlev=nnlev)

    ##############################
    # isosurface
    if initial_surf:
        lev1 = len(rho0) - 1

    ##############################
    # initial vertical position = depths0
    if initial_depth: 
        z_w = part.get_depths_w(simul, x_periodic=x_periodic,
                                y_periodic=y_periodic, ng=ng)
        z = seeding_part.ini_depth(maskrho, simul, depths0, x, y, z, z_w,
                                   ng=ng)
   
    ##############################
    # Release particles on iso-surfaces of a variable
    # Typically isopycnals
    if initial_surf:
        [temp, salt] = part.get_ts_io(simul, x_periodic = x_periodic,
                                      y_periodic = y_periodic, ng=ng)
        [z_r, z_w] = part.get_depths(simul, x_periodic=x_periodic,
                                     y_periodic=y_periodic, ng=ng)
        roms_rho0 = simul.rho0
        rho = rho1_eos(temp, salt, z_r, z_w, roms_rho0)
        ## temporary box used for sigma-interpolation onto surf0 
        lev1 = rho.shape[2] #  Needed to get all levels
        if preserved_meter:
            z, y, x = seeding_part.seed_meter(ic=ic, jc=jc, lev0=lev0, lev1=lev1,
                        nnlev=nnlev, nx_box=nx_box, ny_box=ny_box, dx_box=dx_box,
                        simul=simul, ng=ng)
        else:
            z_box, y_box, x_box = seeding_part.seed_box(ic=ic, jc=jc, lev0=lev0,
                    lev1=lev1, nnx=nnx, nny=nny, iwd=iwd, jwd=jwd, nx=nx, ny=ny)
        
        map_rho = part.map_var(simul, rho, x_box.reshape(-1), y_box.reshape(-1),
                z_box.reshape(-1), ng=ng).reshape(x_box.shape)
       
        del z
        z = np.ndarray(x.shape)
        z = seeding_part.ini_surf(simul, rho0, x, y, z, map_rho, ng=ng)


    ##############################
    # Total number of particles                            
    nq = np.min([len(x.reshape(-1)),nqmx])
    
    ############################################################################
    ''' no need for topocheck anymore as we are using sigma levels'''
    ''' but we have to remove pyticles which are in a masked area'''
    # ipmx : count seeding particles
    # px0, py0, pz0 : initial position for seeding
    # topolim : used in ADV_2D to prevent particles from being seeded below
    #           seafloor
    
    ipmx = 0;
    px0, py0, pz0 = [], [], []
    topolim=0

    if (not adv3d) and (not advzavg): topolim = np.nanmax([topolim, -advdepth])
    elif advzavg: topolim = np.nanmax([topolim, -advdepth - z_thick/2])

    # initializing px0, py0, pz0
    if initial_cond:
        # only true for variables at rho-points
        # maybe not safe nor useful to use i0, j0 ,k0 especially if only load
        # the subdomain for cond (see above)
        i0, j0, k0 = 0, 0, 0
        '''
        A matter of functionality... We take temp, salt, z_w from simulation 
        Then interpolate each variable on px0,py0,pz0
        Finally we compute potential density

        '''
        #########################
        # Data to compute condition
        [temp, salt] = part.get_ts_io(simul, x_periodic=x_periodic,
                                      y_periodic=y_periodic, ng=ng)
        ptemp = partF.interp_3d(x.reshape(-1), y.reshape(-1), z.reshape(-1),
                                temp, ng, nq, i0, j0, k0)
        psalt = partF.interp_3d(x.reshape(-1), y.reshape(-1), z.reshape(-1),
                                salt, ng, nq, i0, j0, k0)
        
        z_w = part.get_depths_w(simul, x_periodic = x_periodic,
                                y_periodic = y_periodic, ng=ng)
        pdepth = partF.interp_3d_w(x.reshape(-1), y.reshape(-1), z.reshape(-1),
                                   z_w, ng, nq, i0, j0, k0)
        #########################
        # boolean condition
        rho_min = -2.5
        rho_max = -0.5
        prho1 = rho1_eos(ptemp, psalt, pdepth, pdepth, simul.rho0)
        pcond = (prho1 > rho_min) & (prho1 < rho_max)
        #########################
        # Remove particles that do not match condition
        ipmx = seeding_part.remove_mask(simul, topolim, x, y, z, px0, py0, pz0, nq, ng=ng, pcond=pcond)
    
    elif part_trap:
        '''
        Particles are released using position of a forward simulation
        Then advected backwards
        
        Here ipmx is used using all particles at it, including nan

        Need to define a time index corresponding to start of backward
        simulation

        '''
        nq_1save, ipmx, px0, py0, pz0 = seeding_part.ini_trap(trap_file,
                                simul, maskrho, itime_fwd=itime_trap,
                               x_periodic=x_periodic, y_periodic=y_periodic, ng=ng)

    else:
        '''
        Returns px0.... not in args
        '''
        
        ipmx = seeding_part.remove_mask(simul, topolim, x, y, z, px0, py0, pz0, nq, ng=ng)

    del x, y, z
    # Else we need the grid box to compute px0, py0, pz0 at each injection time
    nq = ipmx
    print(f'nq = {nq}')

    ############################################################################

    if continuous_injection:
        '''
        we'll multiply nq by the number of injections (still in the limit of
        nqmx) and leave nan for future pyticles which will be released at
        the same positions every time step
        '''
            
        nq_injection = nq
        nq = np.nanmin([nq_injection * N_injection, nqmx])
        print('it would take', nq_injection * N_injection - nqmx, ' more pyticles')
        print('to be able to release through all the simulation')
        nq_1 = nq_injection
        nq_0 = 0
        if debug:
            print('---------------------------------------')
            print(f'nq_injection = {nq_injection}')
            print(f'nq = {nq}')
    else:
        nq_1 = -1  

    ############################################################################
    # nq: total number of particles
    
    px = shared_array(nq,prec='double')
    py = shared_array(nq,prec='double')
    pz = shared_array(nq,prec='double')
    if continuous_injection:
        if part_trap: nq_1 = nq_1save
        px[:nq_1] = px0
        py[:nq_1] = py0
        pz[:nq_1] = pz0
    else: 
        px[:] = px0
        py[:] = py0
        pz[:] = pz0
        del px0,py0,pz0

    ############################################################################
else: # restart = True

    if not continuous_injection:
        # load px,py,pz from restart_file
        print('#################################################')
        print(f'restart_file is {restart_file}')
        nc = Dataset(restart_file, 'r')
        px0 = nc.variables['px'][restart_time, :]
        py0 = nc.variables['py'][restart_time, :]
        pz0 = nc.variables['pz'][restart_time, :]
        nc.close()

        nq = len(px0)
        
        px = shared_array(nq,prec='double')
        py = shared_array(nq,prec='double')
        pz = shared_array(nq,prec='double')

        px[:]=px0
        py[:]=py0
        pz[:]=pz0
        del px0, py0, pz0
        
    else:
        # load px,py,pz from restart_file
        nc = Dataset(restart_file, 'r')
        px0 = nc.variables['px'][restart_time, :]
        py0 = nc.variables['py'][restart_time, :]
        pz0 = nc.variables['pz'][restart_time, :]

        nq = len(px0)
        
        px = shared_array(nq,prec='double')
        py = shared_array(nq,prec='double')
        pz = shared_array(nq,prec='double')

        px[:]=px0
        py[:]=py0
        pz[:]=pz0
        del px0,py0,pz0

        ##################################
        # determine px0,py0,pz0,nq_injection
        # JC not True in case of initital_cond and continuous injection
        # because px0, py0, pz0 may vary upon time
        nq_injection = np.argmax(np.isnan(nc.variables['px'][0, :]))

        px0 = nc.variables['px'][0, :nq_injection]
        py0 = nc.variables['py'][0, :nq_injection]
        pz0 = nc.variables['pz'][0, :nq_injection]
        nc.close()
        
        nq_1 = np.nanmin([nq_injection*((restart_time)//dt_injection+1),nqmx])
        nq_0 = np.nanmin([nq_injection*((restart_time-1)//dt_injection+1),nqmx])
        '''
        if (nq_1<nqmx) and (restart_time%dt_injection)==0:

            nq_0 = nq_1
            nq_1 = np.nanmin([nq_injection*(restart_time/dt_injection+1),nqmx])
            px[nq_0:nq_1]=px0[:nq_1-nq_0]; py[nq_0:nq_1]=py0[:nq_1-nq_0];
            if adv3d: pz[nq_0:nq_1]=pz0[:nq_1-nq_0]
        '''

################################################################################

# Time between 2 frames (in seconds)
delt   = shared_array(2,value=simul.dt*np.abs(dfile)) 
maxvel = shared_array(2,prec='double',value=maxvel0)

# Total number of time steps:
istep = shared_array(1,prec='int',value=-1)

# Index of the  previous (itim[0]) and next(itim[1]) time-step for u,v,w,dz
itim = shared_array(2, prec='int')
itim[:] = [0, 1]

################################################################################
#If using a Adams-Bashforth method we need to have access to previous velocity
# values
# So we create some shared arrays contaning dpx, dpy, dpz for the previous (ab_order-1) time-steps

if timestep[:2]=='AB':
    ab_order = np.int(timestep[-1])
    dpx = shared_array(nq,ab_order,prec='double')
    dpy = shared_array(nq,ab_order,prec='double')
    dpz = shared_array(nq,ab_order,prec='double')
    iab = shared_array(ab_order,prec='int',value=0)
    iab[:]=range(ab_order) # init values of iab


################################################################################
# Call a function as a subprocess
################################################################################


def run_process(my_func):
    ''' 
    We do not pass arguments nor declare any global/local variable so that 
    all variables will be accessible through the subroutine, 
    but only variables declared as global shared arrays (px,py,pz,ptemp,psalt)
    can be modified by the subroutine...
    '''
    results = mp.Queue(); i=0
    proc = mp.Process(target=my_func, args=())
    proc.start(); proc.join()
    
    return results

################################################################################
# Update px,py,px
################################################################################

def update_xyz():
     exec(compile(open('Pyticles_subroutines/update_xyz_largemem.py').read(),\
             'Pyticles_subroutines/update_xyz_largemem.py', 'exec'))

################################################################################
# Compute T,S at each pyticles positions -> ptemp,psalt
################################################################################

def update_ts():   
    exec(compile(open('Pyticles_subroutines/update_ts.py').read(),\
            'Pyticles_subroutines/update_ts.py', 'exec'))
    
################################################################################
# Compute T at each pyticles positions -> ptemp
################################################################################

def update_t():   
    exec(compile(open('Pyticles_subroutines/update_t.py').read(), \
            'Pyticles_subroutines/update_t.py', 'exec'))

################################################################################
# Compute depth at each pyticles positions -> pdepth
################################################################################

def update_depth():   
    exec(compile(open('Pyticles_subroutines/update_depth.py').read(), \
            'Pyticles_subroutines/update_depth.py', 'exec'))

################################################################################
# Compute lon,lat at each pyticles positions -> plon,plat
################################################################################

def update_lonlat():   
    exec(compile(open('Pyticles_subroutines/update_lonlat.py').read(), \
            'Pyticles_subroutines/update_lonlat.py', 'exec'))

################################################################################
# Compute topo at each pyticles positions -> ptopo
################################################################################

def update_topo():   
    exec(compile(open('Pyticles_subroutines/update_topo.py').read(), \
            'Pyticles_subroutines/update_topo.py', 'exec'))

################################################################################
# Compute u,v at each pyticles positions -> pu,pv
################################################################################

def update_uv_2d():   
    exec(compile(open('Pyticles_subroutines/update_uv_2d.py').read(), \
            'Pyticles_subroutines/update_uv_2d.py', 'exec'))

################################################################################
# Compute u,v,w at each pyticles positions -> pu,pv,pw
################################################################################

def update_uvw_3d():
    exec(compile(open('Pyticles_subroutines/update_uvw_3d.py').read(), \
            'Pyticles_subroutines/update_uvw_3d.py', 'exec'))

################################################################################
# Create output file and write time, px,py,pz ( and ptemp,psalt)
################################################################################

def write_output():   
    exec(compile(open('Pyticles_subroutines/write_output.py').read(), \
            'Pyticles_subroutines/write_output.py', 'exec'))
    
################################################################################

def linear(var1,var2,alpha):
    return alpha * var2 + (1.-alpha) * var1
    
################################################################################
# Plot pyticles on a SST map (is done at each time-step)
################################################################################

subtightcoord_saves=[]; subcoord_saves=[];  subsubrange_saves=[]


def plot_rect(rect,line='k'):
    plt.plot([rect[2],rect[2],rect[3],rect[3],rect[2]],\
             [rect[0],rect[1],rect[1],rect[0],rect[0]],line, linewidth=2)
    
################################################################################


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
    #plt.savefig(folderout + simulname + '_' + format(nproc) + '_' + '{0:04}'.format(time+dfile) +'.png', size=None, figure=None, magnification='auto', dpi=150,bbox_inches='tight'); plt.clf()
    


################################################################################

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

################################################################################
# Create output netcdf file for pyticles
################################################################################

#File name
if not restart:
    newfile = folderout + simulname + '_' + format(nproc) + '_' + '{0:04}'.format(filetime) + '.nc'
else:
    newfile = restart_file

print('newfile', newfile)

################################################################################
#START OF THE TIME LOOP
################################################################################
if (nq / nproc) < 10:
    print('----------------------------------------------------')
    print(f'WARNING : Multiprocessing Issue')
    print('number of particles too small relatively to nprocs')
    print(f'in subranges')
    print('use less procs')
    print('----------------------------------------------------')
    #  [nq/nprocs > nprocs-1]
print (' ')
print (' ')
tstart = tm.time()

###############################
time = timerange[0]-dfile;
coord= part.subsection(px, py, nx=nx, ny=ny, offset=50)
run_process(plot_selection)
###############################

#Initialization
pm_s = np.array([]); 

itime = restart_time

for time in timerange:
    print('--------------------------------------------------------------------')
    print(' time is ', time)
    print('--------------------------------------------------------------------')

    alpha_time = time - np.floor(time)

    ############################################################################
    # Define domainstimerange
    # (Find index range (zero based) in which the particles are expected
    # to stay until the next frame)
    ############################################################################
    if debug: print('max. vel. is ',  maxvel)

    tightcoord= part.subsection(px, py, dx, dy, maxvel*0., delt[0], nx, ny, ng)
    coord= part.subsection(px, py, dx, dy, maxvel, delt[0], nx, ny, ng,
                           nadv= nadv)

    nx_s, ny_s = coord[3]-coord[2], coord[1]-coord[0]
    i0 = coord[2]; j0 = coord[0]; 

    print('coord is ', coord)
    
    print('Compute coord............................', tm.time()-tstart)
    tstart = tm.time()   


    ############################################################################
    # Compute ptemp,psalt
    
    # (we are calling a subprocess in order to keep a clean memory after
    #  computation)
    ############################################################################


    if write_ts:
        
        ptemp = shared_array(nq, prec='double')
        psalt = shared_array(nq, prec='double')
        r = run_process(update_ts)

        print('get T,S..................................', tm.time()-tstart)
        tstart = tm.time()   

    elif write_t:

        ptemp = shared_array(nq, prec='double')
        r = run_process(update_t)

        print('get T....................................', tm.time()-tstart)
        tstart = tm.time()   


    if write_uv:
        
        pu = shared_array(nq, prec='double')
        pv = shared_array(nq, prec='double')
        r = run_process(update_uv_2d)

        print('get u,v..................................', tm.time()-tstart)
        tstart = tm.time()   
   
   
    if write_uvw:
        pu = shared_array(nq, prec='double')
        pv = shared_array(nq, prec='double')
        pw = shared_array(nq, prec='double')
        r = run_process(update_uvw_3d)

        print('get u,v,w..................................', tm.time()-tstart)
        tstart = tm.time()


    if write_lonlat:

        plon = shared_array(nq, prec='double')
        plat = shared_array(nq, prec='double')
        r = run_process(update_lonlat)

        print('get lon,lat..............................', tm.time()-tstart)
        tstart = tm.time()   


    if write_depth:

        pdepth = shared_array(nq, prec='double'); 
        r = run_process(update_depth)

        print('get depth................................', tm.time()-tstart)
        tstart = tm.time()   


    if write_topo:

        ptopo = shared_array(nq, prec='double');
        r = run_process(update_topo)

        print('get topo................................', tm.time()-tstart)
        tstart = tm.time()   


    ###########################################################################
    # Write in file
    ###########################################################################   

    r = run_process(write_output)

    if write_ts: del ptemp, psalt
    elif write_t: del ptemp

    if write_lonlat: del plon, plat
    if write_depth: del pdepth
    if write_topo: del ptopo

    if write_uv: del pu, pv
    if write_uvw: del pu, pv, pw

    print('Write in file............................', tm.time()-tstart)
    tstart = tm.time()

    if debug: print('memory usage', \
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)


    ###########################################################################
    # Update px,py,pz 
    # (we are calling it as a subprocess in order to keep a clean memory after
    # computation)
    ###########################################################################

    if debug: print('memory usage',\
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)
    
    ############################################################################
    # Devide domain into subdomains to diminish use of ram
    ############################################################################
    #Automatic division:
    
    nsub_x = 1 + (coord[3] - coord[2]) // 1000
    nsub_y = 1 + (coord[1] - coord[0]) // 1000

    # if domain is periodic, don't divide into subdomains because
    # code cannot handle it yet!
    if x_periodic or y_periodic: nsub_x, nsub_y = 1, 1 

    #no need for 2d advection
    if not adv3d: nsub_x,nsub_y = 1, 1 
    ############################################################################
    # the subtightcoord def needs to be checked to see if ng is needed)
    subtightcoord_saves=[]; subcoord_saves=[];  subsubrange_saves=[]

    if nsub_x*nsub_y == 1:
        subcoord = coord
        subsubrange = list(range(nq))
        r = run_process(update_xyz);
        
    else:
        subtightcoord = [0, 0, 0, 0]
        subcoord = [0 ,0, 0, 0] 
        for isub, jsub in product(list(range(nsub_x)), list(range(nsub_y))):
            
            # subtightcoord will be used for test if particle should be advected
            # or not
            subtightcoord[0] = max(0, tightcoord[0]
                             + jsub*(tightcoord[1] + 1 - tightcoord[0])//nsub_y)
            subtightcoord[1] = min(ny , tightcoord[0] 
                       + (jsub+1) * (tightcoord[1] + 1 - tightcoord[0])//nsub_y)
            subtightcoord[2] = max(0, tightcoord[2] 
                          + isub * (tightcoord[3] + 1 - tightcoord[2])//nsub_x)
            subtightcoord[3] = min(nx , tightcoord[2] 
                       + (isub+1) * (tightcoord[3] + 1 - tightcoord[2])//nsub_x)
            subtightcoord_saves.append(copy(subtightcoord))
              
            #subcoord will be used to compute u,v,w (same than coord)
            subcoord[0] = max(0,
                        subtightcoord[0] - int(np.abs(maxvel[1]*delt[0]/dy)) )
            subcoord[1] = min(ny, 
                        subtightcoord[1] + int(np.abs(maxvel[1]*delt[0]/dy)) )
            subcoord[2] = max(0,
                          subtightcoord[2] - int(np.abs(maxvel[0]*delt[0]/dx)) )
            subcoord[3] = min(nx,
                          subtightcoord[3] + int(np.abs(maxvel[0]*delt[0]/dx)) )
            subcoord_saves.append(copy(subcoord))
            
            subsubrange=[]
            #select pyticles inside subtightcoord only
            for i in range(px.shape[0]):
                if (subtightcoord[0] <= py[i]+0.5+ng < subtightcoord[1]) \
                and (subtightcoord[2] <= px[i]+0.5+ng < subtightcoord[3]):
                    subsubrange.append(i)
            
            subsubrange_saves.append(copy(subsubrange))
                    
        # run process for each subcoord and corresponding subsubrange of points
        for isub, jsub in product(list(range(nsub_x)), list(range(nsub_y))):
            subtightcoord = subtightcoord_saves[jsub + isub*nsub_y]
            subcoord = subcoord_saves[jsub + isub*nsub_y]
            subsubrange = subsubrange_saves[jsub + isub*nsub_y]
            r = run_process(update_xyz);

    ############################################################################

    #if not meanflow and (time+dfile)%1<np.abs(dfile)*1e-2: simul.update(np.int(np.floor(time)+simul.dtime));
    if not meanflow and ( np.round(time+dfile)-(time+dfile)<=np.abs(dfile)*1e-2):
        simul.update(np.int(np.floor(time)+simul.dtime));

    print('Total computation of px,py,pz............', tm.time()-tstart)
    tstart = tm.time()
        
    if debug: print('memory usage',
                    resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)

    ###########################################################################
    # Realease new pyticles (if continuous injection)
    ###########################################################################

    if (continuous_injection) and (nq_1<nqmx) and ((itime+1)%dt_injection)==0:
        
        ###############################################
        # modify seeding patch center to barycenter of
        # previously released particules after a time_step of advection

        if barycentric:
            [ic, jc] = [np.nanmean(px[nq_0:nq_1]), np.nanmean(py[nq_0:nq_1])]
        ###############################
        if preserved_meter:
            z, y, x = seeding_part.seed_meter(ic=ic, jc=jc, lev0=lev0, lev1=lev1,
                        nnlev=nnlev, nx_box=nx_box, ny_box=ny_box, dx_box=dx_box,
                        simul=simul, ng=ng, debug=debug)

        else:
            z, y, x = seeding_part.seed_box(ic=ic, jc=jc, lev0=lev0, lev1=lev1,
                                            iwd=iwd, jwd=jwd, nx=nx, ny=ny,
                                             nnx=nnx, nny=nny, nnlev=nnlev)
        ###############################
        # Release particles at depths0
        if initial_depth: 
            z_w = part.get_depths_w(simul, x_periodic=x_periodic,
                                    y_periodic=y_periodic, ng=ng)
            z = seeding_part.ini_depth(maskrho, simul, depths0, x, y, z,
                                       z_w, ng=ng)
        ipmx = 0
        px0 = []
        py0 = []
        pz0 = []
        ##############################
        # Release particles on iso-surfaces of a variable
        # Typically isopycnals
        if initial_surf:
            [temp, salt] = part.get_ts_io(simul, x_periodic = x_periodic,
                                    y_periodic = y_periodic, ng=ng, coord=coord)
            [z_r, z_w] = part.get_depths(simul, x_periodic=x_periodic,
                                     y_periodic=y_periodic, ng=ng, coord=coord)
            rho = seeding_part.prho(ptemp=temp, psalt=salt, pdepth=z_r)

            ## temporary box used for sigma-interpolation onto surf0 vector
            lev1 = rho.shape[2] #  Needed to get all levels
            if preserved_meter:
                z_box, y_box, x_box = seeding_part.seed_meter(ic=ic, jc=jc, lev0=lev0,
                                      lev1=lev1, nnlev=nnlev, nx_box=nx_box,
                                      ny_box=ny_box, dx_box=dx_box,simul=simul, ng=ng)
            else:
                z_box, y_box, x_box = seeding_part.seed_box(ic=ic, jc=jc, lev0=lev0,
                        lev1=lev1, nnx=nnx, nny=nny, iwd=iwd, jwd=jwd, nx=nx, ny=ny)

            map_rho = part.map_var(simul, rho, x_box.reshape(-1),
                                   y_box.reshape(-1),
                    z_box.reshape(-1), ng=ng, coord=coord).reshape(x_box.shape)

            del z
            z = np.ndarray(x.shape)
            z = seeding_part.ini_surf(simul, rho0, x, y, z, map_rho, ng=ng)

        ################################
        # Add a boolean condition at rho points 
        # Only particles statisfying condition are released
        # Should be consistent with initial release
        if initial_cond:
            temp = part.get_t_io(simul, x_periodic=x_periodic,
                                 y_periodic=y_periodic, ng=ng, coord=coord)
            ini_cond = (temp > 14.) & (temp < 16.)

            pcond = partF.interp_3d(x.reshape(-1), y.reshape(-1), z.reshape(-1),
                                    ini_cond, ng, nq, i0, j0, k0)
            ipmx = seeding_part.remove_mask(simul, topolim, x, y, z, px0, py0,
                                            pz0, nq, ng=ng, pcond=pcond)
        if part_trap:
            '''
            Particles are released using position of a forward simulation
            Then advected backwards
         
            Here ipmx is used using all particles at it, including nan
            ipmx has to be changed to nq_1save

            Need to define a time index corresponding to start of backward
            simulation

            '''

            nq_1save, ipmx, px0, py0, pz0 = seeding_part.ini_trap(trap_file,
                        simul, maskrho, itime_fwd=itime_trap-itime,
                        x_periodic=False, y_periodic=False, ng=ng)
            ipmx = nq_1save

        else:    
            ipmx = seeding_part.remove_mask(simul, topolim, x, y, z, px0, py0,
                                            pz0, nq, ng=ng)

        del x, y, z
        
        ####################################################################
        #nq_0 = nq_1 + 1 
        nq_0 = nq_1
        nq_1 = np.nanmin([nq_0 + ipmx, nqmx]) 
        px[nq_0: nq_1] = px0
        py[nq_0: nq_1] = py0
        if adv3d: pz[nq_0: nq_1] = pz0
        
    ############################################################################
    # Plot particles position (+ SST)
        
    if (time+dfile)%1<np.abs(dfile)*1e-2: run_process(plot_selection)
    
    ############################################################################
    itime += 1
    ############################################################################

    #wait = raw_input("PRESS ENTER TO CONTINUE.")

    print('Plot selection...........................', tm.time()-tstart)
    tstart = tm.time()    
        
    print(' ')
    print(' ')
        
    if debug: print('memory usage', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)
    if debug: print('error',np.sqrt((px[0]-31.)**2+(py[0]-60.)**2))

        
