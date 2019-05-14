import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation

from importlib import reload
import time as tm

import sys
sys.path.append('../../Modules/')
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
zeta = vt.get_var('zeta', roms_file, itime=0, ndims=3)
temp = vt.get_var('temp', roms_file, itime=0, ndims=4)
salt = vt.get_var('salt', roms_file, itime=0, ndims=4)

ocean_time = vt.get_var('time', roms_file)
delta_t = ocean_time[1] - ocean_time[0]

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
##############################################################################
'''
TODO:

a) change axes to have xmin, xmax, ymin, ymax
    Start with max px, py at tend

b)  

'''

roms_path = '/home/jeremy/Bureau/Data/Pyticles/'

# Fixed axes limit 
#xmin = np.int(np.floor(np.nanmin(px)))
#xmax = np.ceil(np.nanmax(px))
#ymin = np.floor(np.nanmin(py))
#ymax = np.ceil(np.nanmax(py))

npts = 5

xmin, xmax, ymin, ymax = vt.get_xy_lim(px, py, npts, nx=1202, ny=1404)

zeta = vt.get_var('zeta', roms_file, itime=0, ndims=3).data
map_var = zeta[ymin:ymax, xmin:xmax]

###########
tstart = tm.time()

fig = plt.figure(figsize=(6, 4), facecolor='white')
ax1 = plt.subplot()
quad1 = ax1.pcolormesh(map_var, shading='gouraud', rasterized=True)
ax1.set_xlabel('px')
ax1.set_ylabel('py')
#scat = ax1.scatter(px[0, :]-xmin, py[0, :]-ymin, zorder=2)
scat = vt.plot_seeding(px-xmin, py-ymin, zorder=4)
scat = plt.scatter(px[0,:]-xmin, py[0,:]-ymin, c=pt[0,:], cmap='jet', zorder=4)
#ax1.set_xlim([200, 400])
#ax1.set_ylim([300, 500])
cb1 = fig.colorbar(quad1, ax=ax1)
#ax1.set_xbound(lower = xmin)
#ax1.set_ybound(lower = ymin)
rgb= [1,1,1]
tend = tm.time() - tstart
print(tend)

#plt.show()

###############
def init():
    quad1.set_array([])
    scat.set_offsets([])
    return quad1, scat

def animate(frame):
    print('*' * 60)
    print(frame)
    roms_file = roms_path + 'chaba_his.' + str(start_file + (frame//5)*5) +'.nc'
    print(roms_file)
    zeta = vt.get_var('zeta', roms_file, itime=frame%5, ndims=3)
    map_var = zeta[ymin:ymax, xmin:xmax]
    quad1.set_array(map_var.ravel())
    X = np.c_[px[frame, :]-xmin+npts, py[frame, :]-ymin+npts]
    #scat.set_offsets(X)
    scat = plt.scatter(px[frame,:]-xmin, py[frame,:]-ymin, c=pt[frame,:], cmap='jet',
                       zorder=4)
    ax1.set_title('time ' + str(frame*delta_t/3600/24) + ' (days)')
    #scat = vt.plot_seeding(px[frame, :]-xmin+npts, py[frame, :]-ymin+npts)
    vt.plot_traj(px-xmin+npts, py-ymin+npts, tend=frame, ax=ax1)
    return quad1

anim = animation.FuncAnimation(fig,animate,frames=25,blit=False,
       repeat=False)
plt.show()


