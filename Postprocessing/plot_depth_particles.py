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
from mpl_toolkits.basemap import Basemap


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
ncfile = '/net/libra/local/tmp/1/gula/particles/ATLBIG/atlbigsig_ATLBIG_dx__RK4_8_0024.nc' #the netcdf file containing particles data
###################################################################################
itime = 0

nc = Dataset(ncfile, 'r')
parameters = nc.simulation
base = nc.base

px = nc.variables['px'][itime,:]
py = nc.variables['py'][itime,:]
pz = nc.variables['pz'][itime,:]

nc.close()

#############
x_periodic = False
y_periodic = False
ng = 0 #number of Ghostpoints _ 1 is enough for linear interp _ 2 for other interp
#############

fifig = './'
config = 'ATLBIG'
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

[z_r,z_w] = part.get_depths(simul,coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)

###################################################################################

pdepth = part.map_var(simul,z_w,px,py,pz,ng,coord=coord)


[plon,plat] = part.map_lonlat(simul,px,py,ng,coord=coord)

pdepth[plon==0] = np.nan
plat[plon==0] = np.nan
plon[plon==0] = np.nan

print('mean depth is', np.nanmean(pdepth))
#ptopo = part.map_topo(simul,px,py,ng,coord=coord)

#############################################

from copy import copy
topo = copy(simul.topo)

###################################################################################
# PLOT
###################################################################################


res             = 'i'
Lx,Ly           = 90e3,90e3 # extent in km
lonLS,latLS     = -32.28,37.30
stride          = 0.2
scaleunit       = 10          #scale unit in km

#lonmin,lonmax   = -41,-24
#latmin,latmax   = 31,44

xmin,xmax = -32.5,-32.1
zmin,zmax   = -2100,-500

fs          =8        # fontsize
symbol      = 'o'       # symbol for LS site
ms          = 4       # marker size for LS site

levels_bck  = np.arange(0,4200,200)   # near fields
cmap_bck    = plt.cm.gray_r

# --- particles ---
ptsize      = 5         # pointsize
ln          = 0         # linewidth

#########################################
# --- select particles along section --- 
sel = []
nq_injection = np.argmax(np.isnan(px))
jsec = np.int(np.nanmean(py[:nq_injection])) ; px0 =np.int(np.nanmean(px[:nq_injection])); res=1.

for pp in np.arange(nq_injection):
    if py[pp] == jsec: sel.append(pp)

#########################################


print(' ... make plot ... ')
fig = plt.figure(figsize=(7,3))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5)

#------------ map --------------------------
print(' ... subplot 1 ... ')
ax = plt.subplot(121,axisbg='gray')
m = Basemap(projection='lcc',resolution='i',\
                    lon_0=lonLS,lat_0=latLS,\
                     width=Lx,height=Ly)

m.drawcoastlines()
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-60,70,stride),labels=[1,0,0,0],linewidth=0,fontsize=fs)
m.drawmeridians(np.arange(-100,100,stride),labels=[0,0,0,1],linewidth=0,fontsize=fs)
m.drawmapscale(m.lonmin+0.5*(m.lonmax-m.lonmin),\
               m.latmin+0.075*(m.latmax-m.latmin),\
               lonLS,latLS,scaleunit,units='km',\
               barstyle='fancy',yoffset=0.015*(m.ymax-m.ymin),fontsize=fs-2)
xr,yr = m(simul.x,simul.y)
m.plot(lonLS,latLS,symbol,color='yellow',markersize=ms,latlon=True)
ctf = m.contourf(xr,yr,simul.topo,levels=levels_bck,cmap=cmap_bck,extend='max')

#for pc in ctf.collections: pc.set_rasterized(True)

scat = m.scatter(plon[:],plat[:],c='r',s=ptsize,linewidths=ln,latlon=True)
ax.tick_params(labelsize=fs)

m.plot(part.rho2psi(simul.x)[:,jsec],part.rho2psi(simul.y)[:,jsec],'y--',lw=0.5,latlon=True)

#------------ section --------------------------
print(' ... subplot 2 ... ')
ax = plt.subplot(122,axisbg='white')
step = 1

z_w = part.rho2psi(z_w)

for k in np.arange(step,50+step,step):
    plt.plot(part.rho2psi(simul.x)[:,jsec],z_w[:,jsec,k],color='gray',lw=0.55)

plt.plot(part.rho2psi(simul.x)[:,jsec],z_w[:,jsec,0],color='k',lw=1)

plt.fill_between(part.rho2psi(simul.x)[:,jsec],-3000*np.ones(simul.x.shape[0]-1),z_w[:,jsec,0],facecolor='k',edgecolor='k',alpha=0.2)
plt.scatter(plon[sel],pdepth[sel],c='r',s=ptsize,linewidths=ln)
plt.xlabel('Longitude [$^{\circ}$E]',fontsize=fs)
plt.xticks((-32.4,-32.2))

plt.xlim((xmin,xmax))
plt.ylim((zmin,zmax))
plt.ylabel('Depth [m]',fontsize=fs)

ax.tick_params(labelsize=fs)
plt.savefig(fifig + 'topo_' + config + '_CV' +'.png',bbox_inches='tight'); plt.clf()






###################################################################################
# plot a vertical section showing topography and particles positions
# centered around px0,py0 along the x-axis
###################################################################################
plt.subplot(1,2,1)

nq_injection = np.argmax(np.isnan(px))

py0 = np.int(np.nanmean(py[:nq_injection])) ; px0 =np.int(np.nanmean(px[:nq_injection])); res=1.
nx1 = px0 - 10; nx2 = px0 + 10
ny1 = py0 - 10; ny2 = py0 + 10

plt.plot((px[np.abs(py-py0)<1.]-px0)*res,pdepth[np.abs(py-py0)<1.],'o')

plt.plot((np.arange(nx1,nx2)-px0)*res-0.5,-topo[nx1:nx2,py0],'-r',lw=5.)
#plt.plot((np.arange(nx1,nx2)-px0)*res-0.5,-topo[nx1:nx2,py0+1],'-or',lw=5.)

###################################################################################
plt.subplot(1,2,2)

plt.plot((py[np.abs(px-px0)<1.]-py0)*res,pdepth[np.abs(px-px0)<1.],'o')

plt.plot((np.arange(ny1,ny2)-py0)*res-0.5,-topo[px0,ny1:ny2],'-r',lw=5.)
#plt.plot((np.arange(ny1,ny2)-py0)*res-0.5,-topo[px0+1,ny1:ny2],'-r',lw=5.)

###################################################################################

plt.savefig(fifig + 'topo_' + config + '_.png', size=None, figure=None, magnification='auto', transparent='true', dpi=400,bbox_inches='tight'); plt.clf()




###################################################################################

