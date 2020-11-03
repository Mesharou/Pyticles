'''
Quick series of plot to illustrate outputs from a pyticles run using 
isopycnal release

- scatter depth & topo as background contour


'''
##############################################################################
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
start_file = 1020 
end_file = 1040 

my_simul = 'polygr_apero'
parameters = my_simul + ' [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

ncfile = '/home/wang/Bureau/Data/Pyticles/fisrt_test/' \
         + 'polygr_apero_fisrt_test_1_1020.nc'
roms_file = '/home/wang/Desktop/Pyticles/ROMS/polgyr_his.1020.nc'
grd_file = '/home/wang/Desktop/Pyticles/ROMS/polgyr_grd.nc'

# OUTPUT PARAM
save_fig = True
save_dir = '/home/wang/Desktop/IDYPOP/Plot_2D/'
gen_name = 'fig_zeta_pdepth_topo'
fmt = '.png'

##############################################################################
# COMPUTING TEST ZONE
##############################################################################
# Case Adv3d
# Given a variable at particle location and a 2D var over wall domain
zeta = vt.get_var('zeta', roms_file, ndims=3)
temp = vt.get_var('temp', roms_file,itime=0, ndims=4)
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
ps = vt.get_var('ps', ncfile).data
pdepth = vt.get_var('pdepth', ncfile).data
ntraj = 20
# Mask version
mask = simul.mask
mask[np.isnan(mask)] = 0.
#if not adv3d: maskrho[simul.topo<-advdepth] = 0.
topo = simul.topo

coord = simul.coord
sst = temp[-1, :, :]

##############################################################################
# Pdepth
roms_path = '/home/wang/Desktop/Pyticles/ROMS/'
npts = 5
xmin, xmax, ymin, ymax = vt.get_xy_lim(px, py, npts=npts, nx=1202, ny=1404)
#scat_min = -700
#scat_max = 0
#zlevels = [-1, -0.5, -0.35, -0.32, -0.3, 0, 0.5 ]
print('=' * 40)
advdepth = -1
for itime in range(px.shape[0]):
    fig = plt.figure(figsize=(6, 4), facecolor='white')
    ax1 = plt.subplot()
    roms_file = roms_path + 'polgyr_his.' + str(start_file + (itime//20)*20) +'.nc'
    # mask land
    vt.mask_contour(ax1, topo, advdepth=advdepth)
    # map_var
    #map_var = topo
    map_var = zeta[itime,:,:]
    quad1 = ax1.contourf(np.ma.masked_invalid(map_var.T), cmap='jet',
                         zorder=1, zlevels=10)
    # text
    ax1.set_xlabel('px')
    ax1.set_ylabel('py')
    ax1.set_title('time ' + str(itime*delta_t/3600/24) + ' (days)')
    # scatter
    scat = vt.anim_scat(pdepth, px, py, itime=itime,
                        zorder=3, size=10)
    cb0 = plt.colorbar(orientation='vertical', pad=0.1, shrink=0.5)
    cb1 = fig.colorbar(quad1, ax=ax1, shrink=0.5)
    # traj
    vt.plot_traj(px, py, tend=itime, linewidth=0.25,
    alpha=0.5)
    ax1.set_xlim(left=xmin, right=xmax)
    ax1.set_ylim(bottom=ymin, top=ymax)
    plt.tight_layout()
    if save_fig:
        fname = save_dir + gen_name + str(itime) + fmt
        plt.savefig(fname)
    plt.close(fig)

