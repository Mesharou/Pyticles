import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.collections import LineCollection
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset


import sys
sys.path.append('../Modules/')
from R_files import load





##############################################################################
def get_var(var, ncfile):
    '''
    Returns the netcdf variable 'var' from ncfilei
    var is a string : name of pyticiles  variable
    return np.array py_var
    use module netCDF4
    use module numpy as np
    '''
    nc = Dataset(ncfile, 'r')
    if var in nc.variables:
         py_var = nc.variables[var][:]
    else:
        py_var = []
        print(f'Error {var} is not found in file {ncfile}')
    nc.close()
    return py_var

#############################################################################

def get_attr(attr, ncfile, it=0):
    '''
    Returns the netcdf attribute 'attr' from ncfile
    attr is a string : name of pyticilesgolabal attribute
    '''
    nc = Dataset(ncfile, 'r')
    if attr in nc.ncattrs():
            py_attr = getattr(nc, attr)
    else:
        py_attr= []
        print(f'Error {attr} is not found in file {ncfile}')
    nc.close()
    return py_attr

##############################################################################
def map_traj(pvar, px=0, py=0, map_var=[], ng=2, mask=[]):
    '''
    Plot particles trajectory with pvar in color 
    Superimposed with another 2D variable 

    parameters : pvar variable at particles points
                 px, py: particles horizontal pozition
                 map_var: 2D variable
         	 mask : land mask : only if map_var is a masked array
    '''
    xmin = int(np.floor(np.nanmin(px))) - ng
    xmax = int(np.ceil(np.nanmax(px))) + ng
    ymin = int(np.floor(np.nanmin(py))) - ng
    ymax = int(np.ceil(np.nanmax(py))) + ng
    if mask.any():
        map_var.mask = mask
        map_var.mask = ~map_var.mask 
    plt.figure()
    vmin = np.nanmin(pvar)
    vmax = np.nanmax(pvar)
    plt.scatter(px, py, c=pdepth, cmap='bone', alpha=0.75, vmin=vmin,
                vmax=vmax, zorder=10)
    cbar1 = plt.colorbar(orientation='horizontal', shrink=0.25, pad=0.05)
    plt.pcolormesh(ma.masked_invalid(map_var), cmap='jet',
    rasterized=True)
    cbar2 = plt.colorbar(shrink=0.25);
    plt.axis([xmin, xmax, ymin, ymax])
    plt.show()
    return 






##############################################################################
# INPUT PARAMETERS
##############################################################################
start_file = 1510
end_file = 1535


my_simul = 'Case_1'
parameters = my_simul + ' [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

ncfile = '/home/jeremy/Bureau/Data/Pyticles/Visual_test/Case_1_Visual_test_12_1510.nc'
roms_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_his.1550.nc'
##############################################################################
# COMPUTING TEST ZONE
##############################################################################
zeta = get_var('zeta', roms_file)
temp = get_var('temp', roms_file)
#zeta.mask = mask
px = get_var('px', ncfile).data
py = get_var('py', ncfile).data
pz = get_var('pz', ncfile).data
pu = get_var('pu', ncfile).data

pdepth = get_var('pdepth', ncfile).data
ntraj = 20
####################
# Function trial map_traj
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
sst = temp[it, -1, :, :]

map_traj(pdepth, px=px, py=py, map_var=zeta_plot, ng=2, mask=simul.mask.T)

map_traj(pu, px=px, py=py, map_var=sst, ng=2, mask=simul.mask.T)









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











