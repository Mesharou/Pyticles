#!/usr/bin/env python
# # Interpolate CROCO variables onto particles position (R_tools/Netcdf4)
#
#
# - Interpolate 3d CROCO variables at particles positions
# and save it to Pyticles netcdf file.
#
# - If Pyticles outputs have a higher frequency than CROCO,
# a time-linear interpolation is applied.
#
# - Both 2D (iso-depth) and 3D Pyticles experiments are supported
#
#
# To run in production mode it is advized to test your code with this notebook,
# then convert it to Python script (Jupytext/jupyter-nbcvonert...)
#  
# ___
#      
#  - 18/12/16:
#      - add capability of computing outputs at subtime-steps period
# ___

# %matplotlib inline

# Load all useful modules

# +
#Plotting modules 
import matplotlib
matplotlib.use('Agg') #Choose the backend (needed for plotting inside subprocess)
import matplotlib.pyplot as plt

#Some standard modules
import sys, os
import numpy as np
import time as tm
from netCDF4 import Dataset
import multiprocessing as mp
import ctypes   

#add the Modules folder in your python PATH
sys.path.append("../Modules/")

#Specific modules needed for pyticles
import pyticles_sig_sa as part
import pyticles_3d_sig_sa as partF

#sys.path.append("/home/gula/Desktop/Work_capella/Python/Python_Modules")
#sys.path.append("/home/gula/Python_Modules")
#sys.path.append("/home2/datawork/jgula/Pyticles/Pyticles_nesea/Modules") 

#Simulations (path, data...)
from R_files import load
from R_netcdf import ionetcdf
from R_vars import var

# to avoid OSError: [Errno -101] NetCDF: HDF error:
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'


# -

# local functions
def linear(var1, var2, alpha):
    "Linear interpolation"
    return alpha * var2 + (1 - alpha) * var1


# get Pyticles file
ncfile = '/home2/datawork/jcollin/Pyticles/tutorials/polgyr_dynamic_injection_1_1000.nc'
meanflow = False

print('Loading file')
with Dataset(ncfile, 'r') as nc:
    parameters = nc.simulation
    base = nc.base
    ng = 0. #nc.ng # ghost points in the horizontal
    nq = len(nc.dimensions['nq'])
    ntime = len(nc.dimensions['time'])

# load simulation parameters

print('Loading simul')
simul = load(simul=parameters)
depths = simul.coord[4]

# get time

# +
ocean_time = ionetcdf.get(ncfile, 'ocean_time', simul)

# old version (output generated before 16/12/20)
#time = np.round(ionetcdf.get(ncfile,'time_int',simul),2)
#time += time[2] - time[1]
# new version
ptime = np.round(ionetcdf.get(ncfile,'time',simul),3)
dtime = ptime[1] - ptime[0]
nqmx  = nq
# -

# Loop on variables

#for varname in ['u','v','rho1','temp','salt']:
for varname in ['rho1']:
    pvar = np.zeros(nq)
    pvarname = 'p' + varname
    
    nc = Dataset(ncfile, 'a')
    if pvarname not in list(nc.variables.keys()): 
        nc.createVariable(pvarname, 'd', ('time','nq',))
    nc.close

    ###################################################################################
    # loop on Pyticles time indices
    # alpha time manage Pyticles experiment with higher frequency outputs than CROCO
    ###################################################################################
    for filetime in ptime[:]:
        itime = int(np.abs(ptime[0] - filetime) / np.abs(dtime))
        print(itime, ptime[itime])
        if not meanflow: simul.update(np.floor(filetime));
        
        alpha_time = filetime - np.floor(filetime)
        
        ###############################################################################
        # Get positions (px,py,pz) from pyticle file
        ###############################################################################
        px = ionetcdf.get(ncfile, 'px', simul, time=itime)
        py = ionetcdf.get(ncfile, 'py', simul, time=itime)

        try:
            adv3d = False
            nc = Dataset(ncfile, 'r')
            advdepth = nc.depth
            nc.close()
        except:
            adv3d = True
            pz = ionetcdf.get(ncfile, 'pz', simul, time=itime)

        print('filetime, itime', filetime, itime)

        coord = part.subsection(px, py, ny=simul.coord[1],
                                nx=simul.coord[3], offset=10)
        tstart = tm.time()

        ######################################################################
        # Get var from simulation file
        ######################################################################
        #compute your variable here:
        if adv3d:
            myvar = var(varname, simul, coord=coord).data
        else:
            myvar = var(varname, simul, depths=[advdepth], coord=coord).data
        
        if not meanflow and alpha_time != 0:
            simul.update(np.ceil(filetime))
            if adv3d:
                myvar2 = var(varname, simul, coord=coord).data
            else:
                myvar2 = var(varname, simul, depth=[advdepth], coord=coord).data
            
            simul.update(np.floor(filetime))
            myvar = linear(myvar, myvar2, alpha_time)
            del myvar2
        
        ######################################################################
        print('get var.......................', tm.time()-tstart)
        tstart = tm.time()
        
        #if ng>0:
        #    myvar = part.periodize3d_fromvar(simul,myvar,coord,
        #               x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)

        ######################################################################
        # Interpolate var to particles positions (px,py,pz)
        ######################################################################
        if adv3d:
            pvar = part.map_var(simul, myvar, px, py, pz, ng, coord=coord)
        else:
            pvar = part.map_var2d(simul, myvar, px, py, ng, coord=coord) 

        print('interpolate to particules.....', tm.time()-tstart)
        tstart = tm.time() 

        ######################################################################
        # write in file
        # for older netcdf version use 
        #  with Dataset(ncfile, 'a') as nc:
        ######################################################################
        with Dataset(ncfile, 'a', format='NETCDF4') as nc:
            nc.variables[pvarname][itime, :] = pvar
            
