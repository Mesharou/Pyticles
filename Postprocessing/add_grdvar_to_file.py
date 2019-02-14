#!/usr/bin/env python
'''

Interpolate grd variables (lon,lat,topo,etc.) at particles positions

and write it in the pyticles netcdf file


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

#sys.path.append("/home/gula/Desktop/Work_capella/Python/Python_Modules")
#sys.path.append("/home/gula/Python_Modules")
#sys.path.append("/home2/datawork/jgula/Pyticles/Pyticles_nesea/Modules") 

#Simulations (path, data...)
from R_files import load
from R_netcdf import ionetcdf
from R_vars import var

###################################################################################

ncfile = sys.argv[-1]; #'/home2/datawork/jgula/Pyticles/Pyticles_nesea/test/nesed_avg_test_adv0000m_14_0060.nc'

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

ng = 0. #nc.ng # ghost points in the horizontal

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

dtime = time[1] - time[0]
nqmx  = nq

###################################################################################

def linear(var1,var2,alpha):
    return alpha * var2 + (1.-alpha) * var1

###################################################################################
# Loop on variables
###################################################################################

for varname in ['lon','lat','topo']:

    pvar = np.zeros(nq)
    
    nc = Dataset(ncfile, 'a')

    if 'p'+varname not in nc.variables.keys(): 
        nc.createVariable('p'+varname,'d',('time','nq',))

    nc.close()
    #nc = Dataset(ncfile, 'a'); nc.variables['pt'][:]; nc.close()

    ###################################################################################
    # loop on time indices
    ###################################################################################
    for filetime in time[:]:
    ###################################################################################

        itime = np.int(np.abs(time[0] - filetime) / np.abs(dtime))
        print itime, time[itime]
        
        alpha_time = filetime - np.floor(filetime)
        
        ###############################################################################
        # Get positions (px,py,pz) from pyticle file
        ###############################################################################
        px = ionetcdf.get(ncfile,'px',simul,time=itime)
        py = ionetcdf.get(ncfile,'py',simul,time=itime)

        print 'filetime, itime', filetime, itime

        coord = part.subsection(px,py,ny=simul.coordmax[1],nx=simul.coordmax[3],offset=10)
        tstart = tm.time()

        ###################################################################################
        # Get var from simulation file
        ###################################################################################
        #compute your variable here:
        if varname=='lon':
            myvar = simul.x
        elif varname=='lat':
            myvar = simul.y
        elif varname=='topo':
            myvar = simul.topo

        ###################################################################################

        print 'get var.......................', tm.time()-tstart
        tstart = tm.time()

        ###################################################################################
        # Interpolate var to particles positions (px,py,pz)
        ###################################################################################
        pvar = part.map_var2d(simul,myvar,px,py,ng) 
        # for test
        # ptopo = part.map_topo(simul,px,py,ng=0)

        print 'interpolate to particules.....', tm.time()-tstart
        tstart = tm.time() 

        ###################################################################################
        # write in file
        ###################################################################################
        nc = Dataset(ncfile, 'a')
        nc.variables['p'+varname][itime,:]=pvar    
        nc.close()















