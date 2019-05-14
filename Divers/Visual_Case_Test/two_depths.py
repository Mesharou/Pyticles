'''
JCollin 05-2019

In this scenario particles were seeded at different depths: depths0 = [-50, -500]

using a simple condition upon pdepth 2 sets of particles are plotted using
different marker shapes

In color a variable pvar is shown at particles location

For a matter of readability, particles are shown every dtime

To keep track of particles, their trajectory is plotted using full time resolution

Particles initial position is plotted in a different color 

'''

from importlib import reload
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.collections import LineCollection
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import time as tm


import sys
sys.path.append('../../Modules/')
# only for JC
sys.path.append('home/jeremy/Bureau/Project/Pyticles/')
from R_files import load
import visual_tools as vt

##############################################################################
# INPUT PARAMETERS
##############################################################################
start_file = 1510
end_file = 1535

my_simul = 'Case_1'
parameters = my_simul + ' [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

ncfile = '/home/jeremy/Bureau/Data/Pyticles/Visual_2_depths/' \
         + 'Case_1_Visual_2_depths_1_1510.nc'
roms_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_his.1550.nc'
grd_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_grd.nc'
##############################################################################
# COMPUTING TEST ZONE
##############################################################################
# Case Adv3d
# Given a variable at particle location and a 2D var over wall domain
zeta = vt.get_var('zeta', roms_file, itime=0, ndims=3)
temp = vt.get_var('temp', roms_file, itime=0, ndims=4)
salt = vt.get_var('salt', roms_file, itime=0, ndims=4)
tend = tm.time() - tstart

topo_roms = vt.get_var('h', grd_file)
#zeta.mask = mask
px = vt.get_var('px', ncfile).data
py = vt.get_var('py', ncfile).data
pz = vt.get_var('pz', ncfile).data
pu = vt.get_var('pu', ncfile).data
pw = vt.get_var('pw', ncfile).data
pt = vt.get_var('pt', ncfile).data
psalt = vt.get_var('ps', ncfile).data
pdepth = vt.get_var('pdepth', ncfile).data
ntraj = 20
####################
# Mask version
mask = simul.mask
mask[np.isnan(mask)] = 0.
#if not adv3d: maskrho[simul.topo<-advdepth] = 0.
topo = simul.topo

coord = simul.coord
sst = temp[-1, :, :]
######################
'''
First scenario:
Particles were released at depths [-500, -50] m
plot pw in seismic diverging color, with different marker shape to distinguish
background variable is topographic contours
dtime : int
'''

dtime = 4
save_fig = False
zlevs =20


indx = pdepth[0, :] > - 400

ng = 2

xmin = int(np.floor(np.nanmin(px))) - ng
xmax = int(np.ceil(np.nanmax(px))) + ng
ymin = int(np.floor(np.nanmin(py))) - ng
ymax = int(np.ceil(np.nanmax(py))) + ng

save_path = '/home/jeremy/Bureau/Data/Pyticles/RESU/Visual_tools/Map_traj/'
fname = save_path + 'test_pw_dtime_' + str(dtime) 

fig = plt.figure(figsize=[8.,8.])

# background 2D variable
plt.contourf(topo_roms, cmap='Greys', levels=zlevs, zorder=1, extend='max')
cbar0 = plt.colorbar(orientation='vertical',shrink=0.5)

# trajectories above 2D background variable
vt.plot_traj(px, py, zorder=2)

# particles vertical velocity pw every dtime
# marker Plus for depths0 = -50 m
# marker 'O' for depths0 = -500 m
vt.plot_part(pw[:,indx], px[:,indx], py[:,indx], ng=2, cmap='seismic',
             vmin=-2e-3, vmax=2e-3, marker='P', zorder=3, size=16, dtime=dtime)
vt.plot_part(pw[:,~indx], px[:,~indx], py[:,~indx], ng=2, cmap='seismic',
             vmin=-2e-3, vmax=2e-3, marker='o',size=16, dtime=dtime, zorder=3)
cbar = plt.colorbar(orientation='horizontal', shrink=0.25, pad=0.1)

# initial position 
vt.plot_seeding(px, py, c='g', marker='X', s=36, zorder=30 )

plt.axis([xmin, xmax, ymin, ymax])
if(save_fig):
    plt.savefig(fname, format='png')

plt.show()

##############################################################################
'''
Psalt
'''

save_fig = True
fname = save_path + 'test_psalt_dtime_' + str(dtime)
figsize = [8., 8.]
dtime=6

fig = plt.figure(figsize=figsize)

# background 2D variable
plt.contourf(topo_roms, cmap='Greys', levels=zlevs, zorder=1, extend='max')
cbar0 = plt.colorbar(orientation='vertical',shrink=0.5)

# trajectories above 2D background variable
vt.plot_traj(px, py, zorder=2)

# particles salinty psalt every dtime
# marker Plus for depths0 = -50 m
# marker 'O' for depths0 = -500 m
cmap = 'gist_rainbow'
vt.plot_part(psalt[:,indx], px[:,indx], py[:,indx], ng=2, cmap=cmap,
             marker='P', zorder=3, size=16, dtime=dtime)
vt.plot_part(psalt[:,~indx], px[:,~indx], py[:,~indx], ng=2, cmap=cmap,
             marker='o',size=16, dtime=dtime, zorder=3)
cbar = plt.colorbar(orientation='horizontal', shrink=0.25, pad=0.1)

# initial position 
vt.plot_seeding(px, py, c='g', marker='X', s=36, zorder=30 )

# Editing axes property
plt.axis([xmin, xmax, ymin, ymax])
# Text
plt.title('Salinity in Psu with topography')

if(save_fig):
    plt.savefig(fname)

plt.show()

