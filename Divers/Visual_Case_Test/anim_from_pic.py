'''
JCollin 05-2019
Generates a series of images using a loop

Example particles advected by vertically averaged velocity around advdepth

'''
##############################################################################
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import matplotlib.colors as col

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
advdepth = -200

my_simul = 'Case_1'
parameters = my_simul + ' [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

save_fig = True
save_dir = '/home/jeremy/Bureau/Data/Pyticles/RESU/Visual_tools/Anim/'
gen_name = 'fig_zavdg_bis_'
fmt = '.png'

ncfile = '/home/jeremy/Bureau/Data/Pyticles/Visual_ZAVG_/' + 'Case_1_Visual_ZAVG__adv200.0m_1_1510.nc'
roms_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_his.1550.nc'
grd_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_grd.nc'
##############################################################################
# COMPUTING TEST ZONE
##############################################################################
zeta = vt.get_var('zeta', roms_file, itime=0, ndims=3)
temp = vt.get_var('temp', roms_file, itime=0, ndims=4)
salt = vt.get_var('salt', roms_file, itime=0, ndims=4)

ocean_time = vt.get_var('time', roms_file)
delta_t = ocean_time[1] - ocean_time[0]

topo_roms = vt.get_var('h', grd_file)
#zeta.mask = mask
px = vt.get_var('px', ncfile).data
py = vt.get_var('py', ncfile).data
pu = vt.get_var('pu', ncfile).data
pt = vt.get_var('pt', ncfile).data
####################
# Function trial map_traj
#####################
# Mask version
mask = simul.mask
mask[np.isnan(mask)] = 0.
#if not adv3d: maskrho[simul.topo<-advdepth] = 0.
topo = simul.topo

fig = plt.figure
#ax = plt.contourf(topo.T, 
#                 [0,-advdepth], cmap = \
#                 col.LinearSegmentedColormap.from_list('my_colormap',['white','lightgray'],256))
#plt.contourf(topo.T,[0,topo.min()],cmap =
#col.LinearSegmentedColormap.from_list('my_colormap',['Gainsboro','gray'],256))
#plt.contour(topo.T,[topo.min()],colors=('k',),linewidths=(0.5,));
#plt.contour(topo.T,[-advdepth],colors=('k',),linewidths=(0.2,));
ax = plt.subplot(111)
vt.mask_contour(ax, topo, advdepth=advdepth)

plt.show()


coord = simul.coord
sst = temp[-1, :, :]
#########################
# Particle age
indx = ~np.isnan(px)
ti = 5
age = np.sum(indx[:ti,:], axis=0)

##############################################################################
'''
Todo: 

vertically averaged (u,v) used for horizontal advection around depths0 = -200 m

background variable : SST

scatter variable: psalt


'''
##################################
# Pdepth
roms_path = '/home/jeremy/Bureau/Data/Pyticles/'
npts = 5
xmin, xmax, ymin, ymax = vt.get_xy_lim(px, py, npts=npts, nx=1202, ny=1404)

scat_min = -700
scat_max = 0 

for itime in range(px.shape[0]):
    fig = plt.figure(figsize=(6, 4), facecolor='white')
    ax1 = plt.subplot()
    roms_file = roms_path + 'chaba_his.' + str(start_file + (itime//5)*5) + '.nc'
    # mask land
    #vt.mask_contour(ax1, topo[xmin:xmax, ymin:ymax], advdepth=advdepth)
    vt.mask_contour(ax1, topo, advdepth=advdepth)

    # map_var
    #map_var = topo[xmin:xmax, ymin:xmax]
    map_var = topo
    quad1 = ax1.contourf(np.ma.masked_invalid(map_var.T), cmap='Greys',
                         zorder=1)
    # text
    ax1.set_xlabel('px')
    ax1.set_ylabel('py')
    ax1.set_title('time ' + str(itime*delta_t/3600/24) + ' (days)')
    # scatter
    #scat = vt.anim_scat(pt, px-xmin+npts, py-ymin+npts, itime=itime, zorder=3,
    #                    size=10)
    scat = vt.anim_scat(pt, px, py, itime=itime, zorder=3,
                        size=10, vmin=10)

    cb0 = plt.colorbar(orientation='vertical', pad=0.1, shrink=0.5)
    cb1 = fig.colorbar(quad1, ax=ax1, shrink=0.5)
    # traj
    #vt.plot_traj(px-xmin+npts, py-ymin+npts, tend=itime, linewidth=0.25, alpha=0.5)
    vt.plot_traj(px, py, tend=itime, linewidth=0.25, alpha=0.5)
    ax1.set_xlim(left=xmin, right=xmax)
    ax1.set_ylim(bottom=ymin, top=ymax)

    if save_fig:
        fname = save_dir + gen_name + str(itime) + fmt 
        plt.savefig(fname)
    plt.close('all')
    
##################################

