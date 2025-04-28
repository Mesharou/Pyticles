#!/usr/bin/env python
'''

Interpolate a 3d variable at particles positions

and write it in the pyticles netcdf file

u and v are implemented as examples

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
#sys.path.remove("/home2/datahome/jgula/Python_Modules_p3")
sys.path.append("../Modules/")
    
#Specific modules needed for pyticles
import pyticles_sig_sa as part
import pyticles_3d_sig_sa as partF

#Simulations (path, data...)
from R_files import load
from R_netcdf import ionetcdf
from R_vars import var

###################################################################################

ncfile = sys.argv[-1]; #'/home2/datawork/jgula/Pyticles/Pyticles_nesea/test/nesed_avg_test_adv0000m_14_0060.nc'
meanflow = False

###################################################################################

print('Loading file')
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

print('Loading simul')
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

for varname in ['u','v','rho1','temp','salt']:
#for varname in ['rho','rho1']:

    pvar = np.zeros(nq)
    
    nc = Dataset(ncfile, 'a')

    if 'p'+varname not in list(nc.variables.keys()): 
        nc.createVariable('p'+varname,'d',('time','nq',))

    nc.close()
    #nc = Dataset(ncfile, 'a'); nc.variables['pt'][:]; nc.close()

    ###################################################################################
    # loop on time indices
    ###################################################################################
    for filetime in time[:]:
    ###################################################################################

        itime = np.int(np.abs(time[0] - filetime) / np.abs(dtime))
        print(itime, time[itime])
        if not meanflow: simul.update(np.floor(filetime));
        
        alpha_time = filetime - np.floor(filetime)
        
        ###############################################################################
        # Get positions (px,py,pz) from pyticle file
        ###############################################################################
        px = ionetcdf.get(ncfile,'px',simul,time=itime)
        py = ionetcdf.get(ncfile,'py',simul,time=itime)


        nc = Dataset(ncfile, 'r')
        adv3d = nc.adv3d

        if adv3d:
            pz = ionetcdf.get(ncfile,'pz',simul,time=itime)
        else:
            advdepth = nc.depth
        nc.close()

        print('filetime, itime', filetime, itime)

        coord = part.subsection(px,py,ny=simul.coordmax[1],nx=simul.coordmax[3],offset=10)
        tstart = tm.time()

        ###################################################################################
        # Get var from simulation file
        ###################################################################################
        #compute your variable here:
        if adv3d:
            myvar = var(varname,simul,coord=coord).data
        else:
            myvar = var(varname,simul,depths=[advdepth],coord=coord).data
        
        if not meanflow and alpha_time != 0:
            simul.update(np.ceil(filetime))
            if adv3d:
                myvar2 = var(varname,simul,coord=coord).data
            else:
                myvar2 = var(varname,simul,depth=[advdepth],coord=coord).data
            simul.update(np.floor(filetime))

            myvar = linear(myvar,myvar2,alpha_time)
            del myvar2
        
        ###################################################################################

        print('get var.......................', tm.time()-tstart)
        tstart = tm.time()
        
        #if ng>0:
        #    myvar = part.periodize3d_fromvar(simul,myvar,coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)

        ###################################################################################
        # Interpolate var to particles positions (px,py,pz)
        ###################################################################################
        if adv3d:
            pvar = part.map_var(simul,myvar,px,py,pz,ng,coord=coord)
        else:
            pvar = part.map_var2d(simul,myvar,px,py,ng,coord=coord) 

        print('interpolate to particules.....', tm.time()-tstart)
        tstart = tm.time() 

        ###################################################################################
        # write in file
        ###################################################################################
        nc = Dataset(ncfile, 'a')
        nc.variables['p'+varname][itime,:]=pvar    
        nc.close()















