#!/usr/bin/env python
# python -i test.py filam [0,1000,0,1600,[-100,0,10]] temp,w,buoy 190 1

###################################################################################
# Load all useful modules (see /home/gula/Desktop/python_v2/Modules/)
###############################################
####################################

#add the Modules folder in your python PATH
import sys
sys.path.append("../Modules/") 
sys.path.remove("/home/gula/Desktop/Work_capella/Python/Python_Modules")


#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

#Some standard modules
import sys, os
import numpy as np
import time as tm
from netCDF4 import Dataset
import numpy.ma as ma

#Specific modules needed for pyticles
import pyticles_sig_sa as part
import pyticles_3d_sig_sa as partF

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#Simulations (path, data...)
from R_files import load
from R_netcdf import ionetcdf

###################################################################################
ncfile = '../WOES_newcode2_sig_WOES_RK4_8_0143.nc' #the netcdf file containing particles data
###################################################################################
itime = 17

nc = Dataset(ncfile, 'r')
parameters = nc.simulation
base = nc.base

px = nc.variables['px'][itime,:]
py = nc.variables['py'][itime,:]
pz = nc.variables['pz'][itime,:]
'''
px = nc.variables['px'][:,117]
py = nc.variables['py'][:,117]
pz = nc.variables['pz'][:,117]
'''
nc.close()

#############
x_periodic = False
y_periodic = False
ng = 0 #number of Ghostpoints _ 1 is enough for linear interp _ 2 for other interp
#############

fifig = './case1'
config = 'case1'
###################################################################################
# load simulation parameters
###################################################################################

print('Loading simul')
simul = load(simul = parameters)
depths = simul.coord[4]

###################################################################################
# get depth, lon, lat
###################################################################################
coord =simul.coord[:4] #= part.subsection(px,py,ny=simul.coordmax[1],nx=simul.coordmax[3],offset=10)

#simul.VertCoordType='OLD'
#[z_r,z_w] = part.get_depths(simul,coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)

simul.VertCoordType='NEW'
[z_r_new,z_w_new] = part.get_depths(simul,coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)

#simul.VertCoordType='OLD'

'''
#--------------------------del sys.modules["-
plt.subplot(1,2,1)
plt.plot(z_w_old[480,213,:].T); 
plt.subplot(1,2,2)
plt.plot(z_w_new[480,213,:].T); 

plt.savefig('test.png'); plt.clf()

#---------------------------
'''



###################################################################################




pdepth = part.map_var(simul,z_w,px,py,pz,ng,coord=coord)
#del z_r,z_w

[plon,plat] = part.map_lonlat(simul,px,py,ng,coord=coord)

pdepth[plon==0] = np.nan
plat[plon==0] = np.nan
plon[plon==0] = np.nan

print('mean depth is', np.nanmean(pdepth))

ptopo = part.map_topo(simul,px,py,ng,coord=coord)


#############################################

print(np.nansum(pdepth<-ptopo))
print(np.nanargmin(pdepth+ptopo))

'''
#print px[2703],py[2703],pz[2703],pdepth[2703],ptopo[2703]
i0=0
j0=0
k0=0

partF.interp_3d_w(px[1],py[1],pz[1]*0,z_w,ng,px.shape[0],i0,j0,k0)
partF.interp_2d(px[1],py[1],z_w[:,:,0],ng,px.shape[0],i0,j0)
'''

from copy import copy
mask = copy(simul.mask)
mask[np.isnan(mask)] = 0.

pmask = part.map_var2d(simul,mask,px,py,ng,coord=coord)

np.nanargmin(pmask)
ip=107

[u,v,w]=part.get_vel_io(simul,pm=simul.pm,pn=simul.pn,coord=coord)

#pu = part.map_var(simul,u,px,py,pz,ng,coord=coord)


topo = copy(simul.topo)
topo[simul.topo==50] = 0.

###################################################################################
# plot a vertical section showing topography and particles positions
# centered around px0,py0 along the x-axis
###################################################################################

py0 = np.int(py[ip]) ; px0 =np.int(px[ip]); res=1.
nx1 = px0 - 45; nx2 = px0 + 45
ny1 = py0 - 45; ny2 = py0 + 45

plt.plot((px[np.abs(py-py0)<1.]-px0)*res,pdepth[np.abs(py-py0)<1.],'o')

plt.plot((np.arange(nx1,nx2)-px0)*res,5000*u[nx1:nx2,py0,-1]-25.,'--oy',lw=5.)
plt.plot((np.arange(nx1,nx2)-px0)*res,5000*u[nx1:nx2,py0+1,-1]-25.,'-oy',lw=5.)

plt.plot((np.arange(nx1,nx2)-px0)*res-0.5,-topo[nx1:nx2,py0],'--or',lw=5.)
plt.plot((np.arange(nx1,nx2)-px0)*res-0.5,-topo[nx1:nx2,py0+1],'-or',lw=5.)

#plt.plot((np.arange(nx1,nx2)-px0)*res,-part.rho2psi(simul.topo)[nx1:nx2,py0],'--k',lw=5.)
plt.ylim([-100, 0])
plt.xlim([-2, 2])

###################################################################################

plt.savefig(fifig + 'topovertx_' + config + '_.png', size=None, figure=None, magnification='auto', transparent='true', dpi=400,bbox_inches='tight'); plt.clf()

###################################################################################
# plot a vertical section showing topography and particles positions
# centered around px0,py0 along the y-axis
###################################################################################

plt.plot((py[np.abs(px-px0)<1.]-py0)*res,pdepth[np.abs(px-px0)<1.],'o')

plt.plot((np.arange(ny1,ny2)-py0)*res,5000*v[px0,ny1:ny2,-1]-25.,'--y',lw=5.)
plt.plot((np.arange(ny1,ny2)-py0)*res,5000*v[px0+1,ny1:ny2,-1]-25.,'-y',lw=5.)

plt.plot((np.arange(ny1,ny2)-py0)*res-0.5,-topo[px0,ny1:ny2],'--r',lw=5.)
plt.plot((np.arange(ny1,ny2)-py0)*res-0.5,-topo[px0+1,ny1:ny2],'-r',lw=5.)

plt.ylim([-100, 0.])
plt.xlim([-2,2])

###################################################################################

plt.savefig(fifig + 'topoverty_' + config + '_.png', size=None, figure=None, magnification='auto', transparent='true', dpi=400,bbox_inches='tight'); plt.clf()

###################################################################################

