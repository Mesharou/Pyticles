from importlib import reload
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.collections import LineCollection
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import time as tm


import sys
sys.path.append('../Modules/')
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

##############################################################################
# Case Adv3d
# Given a variable at particle location and a 2D var over wall domain
tstart = tm.time()
zeta = vt.get_var('zeta', roms_file, itime=0, ndims=3)
temp = vt.get_var('temp', roms_file, itime=0, ndims=4)
salt = vt.get_var('salt', roms_file, itime=0, ndims=4)
tend = tm.time() - tstart
print(f'time using ndims : {tend}')

topo_roms = vt.get_var('h', grd_file)
#zeta.mask = mask
px = vt.get_var('px', ncfile).data
py = vt.get_var('py', ncfile).data
pz = vt.get_var('pz', ncfile).data
pu = vt.get_var('pu', ncfile).data
pw = vt.get_var('pw', ncfile).data
pt = vt.get_var('pt', ncfile).data
pdepth = vt.get_var('pdepth', ncfile).data
ntraj = 20
####################
# Function trial map_traj
#####################
# Mask version
mask = simul.mask
mask[np.isnan(mask)] = 0.
#if not adv3d: maskrho[simul.topo<-advdepth] = 0.
topo = simul.topo

coord = simul.coord
sst = temp[-1, :, :]
##################
# TEST plot_traj : independant function
'''
Need to ensure we can edit figure outside function

Now we want to observe 2 sets of particles 
EX : we find that some particles have a strange pdepth at = 0
we define indx = pdepth[0, :] > - 49 (supposed to be initialized around -50 m)
Using plot_traj kwargs we specify another marker for these particles
'''


zlevs = [5, 10, 25, 50, 75, 100, 200, 300, 400, 500, 700, 800, 900, 1000, 1100,
1200]
indx = pdepth[0, :] > - 400

ng = 2

xmin = int(np.floor(np.nanmin(px))) - ng
xmax = int(np.ceil(np.nanmax(px))) + ng
ymin = int(np.floor(np.nanmin(py))) - ng
ymax = int(np.ceil(np.nanmax(py))) + ng
		
fig = plt.figure()
vt.plot_traj(pw[:,indx], px[:,indx], py[:,indx], ng=2, cmap='seismic',
             vmin=-1e-3, vmax=1e-3, marker='X', zorder=20, size=50)
vt.plot_traj(pw[:,~indx], px[:,~indx], py[:,~indx], ng=2, cmap='seismic', 
             vmin=-2e-3, vmax=2e-3, marker='o')
cbar = plt.colorbar(orientation='horizontal', shrink=0.25, pad=0.1)

levs=20
plt.contourf(topo_roms, cmap='Greys', levels=zlevs,  extend='max')
cbar0 = plt.colorbar(orientation='vertical',shrink=0.5)
plt.axis([xmin, xmax, ymin, ymax])
plt.show()



toto
##################
# pdepth & zeta
vt.map_traj(pdepth, px=px, py=py, map_var=zeta, ng=2, mask=simul.mask.T,
           xlabel='pdepth', ylabel='zeta')
# pu & sst
vt.map_traj(pu, px=px, py=py, map_var=sst, ng=2, mask=simul.mask.T, coord=coord)
# pdepth & topo
vt.map_traj(pdepth, px=px, py=py, map_var=topo, ng=2, cmap1='jet', cmap2='terrain')
# pw and sst 
vt.map_traj(pw, px=px, py=py, map_var=sst, ng=2, cmap1='seismic',
            cmap2='ocean', vmin1=-2e-3, vmax1=2e-3, mask=simul.mask.T)
# pt and topo
vt.map_traj(pt, px=px, py=py, map_var=topo.T, ng=2, cmap1='jet',
            cmap2='terrain', method='contourf', vmax2=2000)
plt.plot(px, py, color='grey', alpha = 1, zorder=2)

plt.show()

plt.close('all')

####################
# Contourf test
'''
Issue: when truncating contour after using xlim ylim
       when using levs = int 
	   contours are delimited using the whole domain
	   resulting in a loss of graphical resolution in subdomaine coord =[ymin, ymax...]
'''
np.min(topo_roms)
np.max(topo_roms)
pvar = pw
cmap1='seismic'
vmin1=-2e-3; vmax1=2e-3

xmin = int(np.floor(np.nanmin(px))) - ng
xmax = int(np.ceil(np.nanmax(px))) + ng
ymin = int(np.floor(np.nanmin(py))) - ng
ymax = int(np.ceil(np.nanmax(py))) + ng


levs= [6000, 5000, 4000, 2000, 1000, 500, 200, 100, 50, 5]
levs=20
plt.contourf(topo_roms[ymin:ymax, xmin:xmax], cmap='Greys', levels=levs, extend='min')
cbar0 = plt.colorbar(orientation='horizontal',shrink=0.5,)

plt.scatter(px-xmin, py-ymin, c=pvar, cmap=cmap1, alpha=0.75, vmin=vmin1,
                vmax=vmax1, zorder=10)
cbar1 = plt.colorbar()
#plt.axis([ymin, ymax, xmin, xmax])

plt.show()

#####################

figsize=(8, 6)

plt.figure()
plt.subplot(221)
vt.map_traj(pdepth, px=px, py=py, map_var=zeta, ng=2, mask=simul.mask.T,
           xlabel='pdepth', ylabel='zeta')
plt.plot(px, py, color='grey', alpha = 1, zorder=2)

plt.subplot(222)
vt.map_traj(pu, px=px, py=py, map_var=sst, ng=2, mask=simul.mask.T)
plt.plot(px, py, color='grey', alpha = 1, zorder=2)

plt.subplot(223)
vt.map_traj(pw, px=px, py=py, map_var=sst, ng=2, cmap1='seismic',
            cmap2='ocean', vmin1=-2e-3, vmax1=2e-3, mask=simul.mask.T)
plt.plot(px, py, color='grey', alpha = 1, zorder=2)

plt.subplot(224)
vt.map_traj(pt, px=px, py=py, map_var=topo, ng=2, cmap1='jet',
            cmap2='terrain')
plt.plot(px, py, color='grey', alpha = 1, zorder=2)

plt.show()

######################
# pt map topo 

plt.figure()
vt.map_traj(pt, px=px, py=py, map_var=topo, ng=2, cmap1='jet',
            cmap2='terrain')
plt.plot(px, py, color='grey', alpha = 1, zorder=2)
plt.show()


######################
# pdepth map zeta 

plt.figure()
vt.map_traj(pdepth, px=px, py=py, map_var=zeta, ng=2, mask=simul.mask.T,
           xlabel='pdepth', ylabel='zeta')
plt.plot(px, py, color='grey', alpha = 1, zorder=2)

plt.show()



######################
# pw map sst 

plt.figure()
plt.plot(px, py, color='grey', alpha = 1, zorder=2)
vt.map_traj(pw, px=px, py=py, map_var=sst, ng=2, cmap1='seismic',
            cmap2='ocean', vmin1=-2e-3, vmax1=2e-3, mask=simul.mask.T)
plt.show()


#####################
# Mask version
mask = simul.mask
mask[np.isnan(mask)] = 0.
#if not adv3d: maskrho[simul.topo<-advdepth] = 0.
topo = simul.topo

it = 0
xmin = int(np.floor(np.nanmin(px)))
xmax = int(np.ceil(np.nanmax(px)))
ymin = int(np.floor(np.nanmin(py)))
ymax = int(np.ceil(np.nanmax(py)))


zeta_plot = zeta[it, :, :]
zeta_plot.mask = simul.mask.T
zeta_plot.mask = ~zeta_plot.mask

plt.scatter(px[:,:], py[:,:], c=pdepth[:,:], cmap='bone',
             alpha=0.75, vmin=-100, zorder=10)
cbar1 = plt.colorbar(orientation='horizontal', shrink=0.25, pad=0.05)
plt.pcolormesh(ma.masked_invalid(zeta_plot[:,:]), cmap='jet', rasterized=True)
cbar2 = plt.colorbar(shrink=0.25);
plt.axis([xmin -2, xmax + 2, ymin -2, ymax +2 ])

plt.show()
plt.close()


######################
# Scatter test
xmin = int(np.floor(np.nanmin(px)))
xmax = int(np.ceil(np.nanmax(px)))
ymin = int(np.floor(np.nanmin(py)))
ymax = int(np.ceil(np.nanmax(py)))

it = 0

zeta_plot = zeta[it, :, :]
zeta_plot *= simul.mask.T
axs = plt.axes
plt.scatter(px[:,:], py[:,:], c=pdepth[:,:], cmap='cool',
             alpha=0.75, vmin=-100, zorder=10)
plt.plot(px[:, :], py[:, :], color='grey', alpha = 1, zorder=2)
plt.contourf(zeta_plot, cmap='autumn', alpha=1, zorder=1)
plt.axis([xmin, xmax, ymin, ymax])
plt.colorbar(shrink=0.25);

#plt.colorbar()
plt.show()

plt.close('all')

######################
# Colored segments
my_rgb = (1, 1, 1)
[nt, npart] = pz.shape
pzmax = np.nanmax(pz)
pzmin = np.nanmin(pz)

## small test
ip = 0
x = px[:, ip]
y = py[:, ip]
z = pz[:, ip]


# Create a set of line segments so that we can color them individually
# This creates the points as a N x 1 x 2 array so that we can stack points
# together easily to get the segments. The segments array for line collection
# needs to be (numlines) x (points per line) x 2 (for x and y)
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)


fig, axs = plt.subplots(2, 1, sharex=True, sharey=True)
# Create a continuous norm to map from data points to colors
norm = plt.Normalize(z.min(), z.max())
lc = LineCollection(segments, cmap='jet', norm=norm)
# Set the values used for colormapping
lc.set_array(z)
lc.set_linewidth(2)
line = axs[0].add_collection(lc)
fig.colorbar(line, ax=axs[0])

# Use a boundary norm instead
#cmap = ListedColormap(['r', 'g', 'b'])
#norm = BoundaryNorm([-1, -0.5, 0.5, 1], cmap.N)
#lc = LineCollection(segments, cmap=cmap, norm=norm)
#lc.set_array(z)
#lc.set_linewidth(2)
#line = axs[1].add_collection(lc)
#fig.colorbar(line, ax=axs[1])

axs[0].set_xlim(x.min(), x.max())
axs[0].set_ylim(y.min(), y.max())
plt.show()

plt.close('all')











