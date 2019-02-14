#!/usr/bin/env python
'''

Interpolate depth in meters at particles positions

and write it in the netcdf file (pdepth)

!---------------------------------------------------------------------------------------------
! 18/12/16:
!     - add capability of computing outputs at subtime-steps period
!---------------------------------------------------------------------------------------------

'''
###################################################################################
# Load all useful modules
##################################################################################

#Plotting modules 
import matplotlib
matplotlib.use('Agg') #Choose the backend (needed for plotting inside subprocess)
import matplotlib.pyplot as plt

#Some standard modules
import sys, os
import numpy as np
import time as tm
from netCDF4 import Dataset
import multiprocessing as mp; import ctypes   

#add the Modules folder in your python PATH
sys.path.remove("/home2/datahome/jgula/Python_Modules")
sys.path.append("../Modules/")
sys.path.append("/home2/datahome/jgula/Python_Modules")
#try:
#    sys.path.remove("/home/gula/Desktop/Work_capella/Python/Python_Modules")
    #sys.path.remove("/home2/datahome/jgula/Python_Modules")
#except:
#    sys.path.remove("/home/gula/Python_Modules")
    
#Specific modules needed for pyticles
import pyticles_sig_sa as part
import pyticles_3d_sig_sa as partF

#Simulations (path, data...)
from R_files import load
from R_netcdf import ionetcdf



###################################################################################

ncfile = sys.argv[-1]; #'/home2/datawork/jgula/Pyticles/Pyticles_nesea/test/nesed_avg_test_adv0000m_14_0060.nc'

meanflow = False

###################################################################################

print 'Loading file'
nc = Dataset(ncfile, 'r')
parameters = nc.simulation
base = nc.base

#check if periodic
if nc.x_periodic==1: x_periodic=True
elif nc.x_periodic==0: x_periodic=False

if nc.y_periodic==1: y_periodic=True
elif nc.y_periodic==0: y_periodic=False

ng = 0 # ghost points in the horizontal

nq = len(nc.dimensions['nq'])
ntime = len(nc.dimensions['time'])
nc.close()

###################################################################################
# load simulation parameters
###################################################################################

print 'Loading simul'
simul = load(simul = parameters)
depths = simul.coord[4]


###################################################################################
# get time
###################################################################################

ocean_time = ionetcdf.get(ncfile,'ocean_time',simul)

# old version (output generated before 16/12/20)
#time = np.round(ionetcdf.get(ncfile,'time_int',simul),2)
#time += time[2] - time[1]
# new version
time = np.round(ionetcdf.get(ncfile,'time',simul),3)

dtime = time[2] - time[1]
nqmx  = nq

###################################################################################
# Create new variable in pyticle file
###################################################################################

print 'Creating variable'
nc = Dataset(ncfile, 'a')

if 'pdepth' not in nc.variables.keys(): 
        nc.createVariable('pdepth','d',('time','nq',))

nc.close()

tstart = tm.time()

###################################################################################

def linear(var1,var2,alpha):
    return alpha * var2 + (1.-alpha) * var1

###################################################################################
for filetime in time[:]:
###################################################################################

    itime = np.int(np.abs(time[0] - filetime) / np.abs(dtime))
    print itime, time[itime]
    if not meanflow: simul.update(np.floor(filetime));
    
    alpha_time = filetime - np.floor(filetime)
    
    
    ###################################################################################
    # Get positions (px,py,pz) from pyticle file
    ###################################################################################
    #nc = Dataset(ncfile, 'r')
    
    px = ionetcdf.get(ncfile,'px',simul,time=itime)
    py = ionetcdf.get(ncfile,'py',simul,time=itime)
    pz = ionetcdf.get(ncfile,'pz',simul,time=itime)

    #nc.close()
    print 'load px,py,pz.................', tm.time()-tstart
    tstart = tm.time()
    
    coord = part.subsection(px,py,ny=simul.coordmax[1],nx=simul.coordmax[3],offset=10)

    
    ###################################################################################
    # Get depth (z_w) from simulation file
    ###################################################################################
    [_,z_w] = part.get_depths(simul,coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)
    del _
    
    if not meanflow and alpha_time != 0:
        simul.update(np.ceil(filetime))
        [_,z_w2] = part.get_depths(simul,coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)
        del _
        simul.update(np.floor(filetime))

        z_w = linear(z_w,z_w2,alpha_time)
        del z_w2
    
    
    print 'get var.......................', tm.time()-tstart
    tstart = tm.time() 

    ###################################################################################
    # Interpolate depth to particles positions (px,py,pz)
    ###################################################################################
    pz_w = part.map_var(simul,z_w,px,py,pz,ng,coord=coord)
    
    print 'interpolate to pyticles.....', tm.time()-tstart
    tstart = tm.time()
    
    ###################################################################################
    # write in file
    ###################################################################################
    nc = Dataset(ncfile, 'a')
    nc.variables['pdepth'][itime,:]=pz_w
    nc.close()

    print 'write in file...............', tm.time()-tstart
    tstart = tm.time()
















