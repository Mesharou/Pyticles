'''
JCollin 05-2019
Generates a series of images using a loop

anim_barycentric.py
   example were particles at t(i+1) are released at the 
   barycenter of previous of particles released at t(i)
   after 1 * dt of advection
                  
   topography is displayed as a background map
   particles age in color


'''
##############################################################################
# Importing python library
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation

from importlib import reload
import time as tm

import sys

# import local library
sys.path.append('../../Modules/')
sys.path.append('home/jeremy/Bureau/Project/Pyticles/')
from R_files import load
import visual_tools as vt
##############################################################################
# INPUT PARAMETERS
# start_file, end_file: int, file indexes
# my_simul : same as the one used in R_files (used to retrieve ROMS outputs)
# save_dir : directory to save figures
# gen_name : generic figures name
##############################################################################
start_file = 1510
end_file = 1535

my_simul = 'Case_1'
parameters = my_simul + ' [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

save_fig = True
save_dir = '/home/jeremy/Bureau/Data/Pyticles/RESU/Visual_tools/Anim/'
gen_name = 'fig_barycentric_'
fmt = '.png'

ncfile = '/home/jeremy/Bureau/Data/Pyticles/Barycentric/' \
         + 'Case_1_Barycentric_12_1550.nc'
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
#########################
# Particle age
indx = ~np.isnan(px)
ti = 5
age = np.sum(indx[:ti,:], axis=0)

##############################################################################
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
    # map_var
    map_var = topo[ymin:ymax, xmin:xmax]
    quad1 = ax1.contourf(np.ma.masked_invalid(map_var),cmap='Greys', rasterized=True,
                         zorder=1)
    # text
    ax1.set_xlabel('px')
    ax1.set_ylabel('py')
    ax1.set_title('time ' + str(itime*delta_t/3600/24) + ' (days)')
    # scatter
    scat = vt.anim_scat(pdepth, px-xmin, py-ymin, itime=itime, zorder=3,
                        size=10)
    cb0 = plt.colorbar(orientation='vertical', pad=0.1, shrink=0.5)
    cb1 = fig.colorbar(quad1, ax=ax1, shrink=0.5)
    # traj
    vt.plot_traj(px-xmin, py-ymin, tend=itime, linewidth=0.25, alpha=0.5)
    if save_fig:
        fname = save_dir + gen_name + str(itime) + fmt 
        plt.savefig(fname)
    plt.close('all')
##################################
# Age
gen_name = 'barycentric_age_'
part_age = np.ndarray(px.shape)
save_fig = True
for itime in range(px.shape[0]):
    part_age[itime, :] = vt.age(px, itime)
    fig = plt.figure(figsize=(6, 4), facecolor='white')
    ax1 = plt.subplot()
    roms_file = roms_path + 'chaba_his.' + str(start_file + (itime//5)*5) + '.nc'
    # map_var
    map_var = topo[ymin:ymax, xmin:xmax]
    quad1 = ax1.contourf(np.ma.masked_invalid(map_var),cmap='Greys', rasterized=True,
                         zorder=1)
    # text
    ax1.set_xlabel('px')
    ax1.set_ylabel('py')
    ax1.set_title('time ' + str(itime*delta_t/3600/24) + ' (days)')
    # scatter
    scat = vt.anim_scat(part_age, px-xmin, py-ymin, itime=itime, zorder=3,
                        size=10, vmax=px.shape[0])
    cb0 = plt.colorbar(orientation='vertical', pad=0.1, shrink=0.5)
    cb1 = fig.colorbar(quad1, ax=ax1, shrink=0.5)
    # traj
    vt.plot_traj(px-xmin, py-ymin, tend=itime, linewidth=0.25, alpha=0.5)
    if save_fig:
        fname = save_dir + gen_name + str(itime) + fmt
        plt.savefig(fname)
    plt.close('all')

