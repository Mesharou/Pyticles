'''
JCollin 05-2019

Estimation of numerical Efficiency between linear and cubic spatial
interpolation

'''

##############################################################################
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np

from importlib import reload

import sys

sys.path.append('../../Modules/')
sys.path.append('home/jeremy/Bureau/Project/Pyticles/')
sys.path.append('')
from R_files import load
import visual_tools as vt
from R_tools import rho1_eos
##############################################################################
# INPUT PARAMETERS
##############################################################################
start_file = 1510
end_file = 1535

my_simul = 'Case_1'
parameters = my_simul + ' [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

save_fig = False
save_dir = '/home/jeremy/Bureau/Data/Pyticles/RESU/Visual_tools/' \
         + 'Numerical_Schemes/'
gen_name = 'fig_num_scheme_'
fmt = '.png'

nc_linear = '/home/jeremy/Bureau/Data/Pyticles/Linear_adv_salt/' \
          + 'Case_1_Linear_adv_salt_6_1510.nc'
nc_cubic = '/home/jeremy/Bureau/Data/Pyticles/Cubic_adv_salt/' \
          + 'Case_1_Cubic_adv_salt_6_1510.nc'

#nc_new_linear = '/home/jeremy/Bureau/Data/Pyticles/New_Version_Linear/' \
#        + 'Case_1_New_Version_Linear_6_1510.nc'

roms_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_his.1550.nc'
grd_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_grd.nc'

#############################################################################
itime = 3
# new linear
'''
px_new = vt.get_var('px', nc_new_linear, itime=itime)
py_new = vt.get_var('py', nc_new_linear, itime=itime)
pt_new = vt.get_var('pt', nc_new_linear, itime=itime)
pdepth_new = vt.get_var('pdepth', nc_new_linear, itime=itime)
'''
# old linear and new cubic
px_lin = vt.get_var('px', nc_linear)
px_cub = vt.get_var('px', nc_cubic)
py_lin = vt.get_var('py', nc_linear)
py_cub = vt.get_var('py', nc_cubic)
pt_lin = vt.get_var('pt', nc_linear)
pt_cub = vt.get_var('pt', nc_cubic)
ps_lin = vt.get_var('ps', nc_linear)
ps_cub = vt.get_var('ps', nc_cubic)

pdepth_lin = vt.get_var('pdepth', nc_linear)
pdepth_cub = vt.get_var('pdepth', nc_cubic)
ocean_time = vt.get_var('ocean_time', nc_cubic)
# comparing cubic and old linear
np.std(pdepth_lin - pdepth_cub)
np.std(pt_lin - pt_cub)


num_bins = 20

fig, ax = plt.subplots()

# the histogram of the data
#n, bins, patches = ax.hist(px_cub-px_lin, num_bins, density=1)
n, bins, patches = ax.hist(pdepth_cub[-1,:]-pdepth_lin[-1,:], num_bins, density=1,
histtype='step');

# add a 'best fit' line
#y = mlab.normpdf(bins, mu, sigma)
#ax.plot(bins, y, '--')
ax.set_xlabel('px');
ax.set_ylabel('Probability density');
ax.set_title('Difference between linear and cubic interpolation');

# Tweak spacing to prevent clipping of ylabel
fig.tight_layout();
plt.show()
##################################
# 3D histogram

#####################
# comparing old and new linear
'''
num_bins = 55

fig, ax = plt.subplots()

# the histogram of the data
#n, bins, pates = ax.hist(px_cub-px_lin, num_bins, density=1)
n, bins, patches = ax.hist(pdepth_new-pdepth_lin, num_bins, density=1,
histtype='step');

# add a 'best fit' line
#y = mlab.normpdf(bins, mu, sigma)
#ax.plot(bins, y, '--')
ax.set_xlabel('px');
ax.set_ylabel('Probability density');
ax.set_title('Difference between linear and cubic interpolation');

# Tweak spacing to prevent clipping of ylabel
fig.tight_layout();
plt.show()

###
plt.plot(pdepth_new-pdepth_lin)
plt.show()
# pt 
plt.plot(pt_new-pt_lin)
plt.show()
# pw 
plt.plot(pw_new-pw_lin)
plt.show()
'''
###########################################################
'''
subplot std of numerical schemes difference in function of time
variables : px/py || pt
            ------------
            pdepth|| EMPTY

MISSING psalt and prho1 (Simulations are running though)
'''
var_lin = px_lin 
var_cub = px_cub
time_std = np.std(var_lin-var_cub, axis=1)
time_mean = np.mean(var_lin - var_cub, axis=1)
roms_rho0 = simul.rho0
pm_bar = np.mean(simul.pm)
day_time = (ocean_time- ocean_time[0])/3600/24
prho1 = rho1_eos(pt_cub, ps_cub, pdepth_cub, pdepth_cub, roms_rho0)

fig = plt.figure
ax = plt.subplot(221)
plt.plot(day_time, time_std/pm_bar*1e-3)
plt.plot(day_time, np.std(py_lin - py_cub, axis=1 )/pm_bar*1e-3 )
plt.legend(('px', 'py'))
plt.ylabel('std in km')

ax = plt.subplot(222)
plt.plot(day_time, np.std(pt_lin - pt_cub, axis=1 )  )
#plt.plot( np.std(py_lin - py_cub, axis=1 ))
plt.legend(('std(pt)',))

ax=plt.subplot(223)
plt.plot(day_time, np.std(pdepth_lin - pdepth_cub, axis=1))
plt.ylabel('std(pdepth) (m)')
plt.xlabel('time in days')


ax = plt.subplot(224)
plt.tight_layout()
plt.show()
############################################################
roms_path = '/home/jeremy/Bureau/Data/Pyticles/'
npts = 5
xmin, xmax, ymin, ymax = vt.get_xy_lim(px_lin, py_lin, npts=npts, nx=1202, ny=1404)
scat_min = -700
scat_max = 0
save_dir = '/home/jeremy/Bureau/Data/Pyticles/RESU/Numerical_Schemes/'
delta_t = ocean_time[1] - ocean_time[0]

save_fig = True

for itime in range(px_cub.shape[0]):
    fig = plt.figure(figsize=(12, 8), facecolor='white')
    ax1 = plt.subplot()
    roms_file = roms_path + 'chaba_his.' + str(start_file + (itime//5)*5) + '.nc'
    # map_var
    zeta= vt.get_var('zeta', roms_file, itime=itime%5)
    map_var = zeta[ymin:ymax, xmin:xmax]
    quad1 = ax1.pcolormesh(map_var, shading='ocean', rasterized=True, zorder=1)
    # text
    ax1.set_xlabel('px')
    ax1.set_ylabel('py')
    ax1.set_title('time ' + str(itime*delta_t/3600/24) + ' (days)')
    # scatter
    scat2 = vt.anim_scat(pdepth_lin, px_lin-xmin, py_lin-ymin, itime=itime, zorder=3,
                        size=10)
    scat = vt.anim_scat(pdepth_cub, px_cub-xmin, py_cub-ymin, itime=itime,
                        zorder=4,size=20, marker = '+')
   # cb0axes = fig.add_axes([left, bottom, width, height])
    cb0 = plt.colorbar(orientation='vertical', pad=0.1, shrink=0.5)
    cb1 = fig.colorbar(quad1, ax=ax1, shrink=0.5)
    # traj
    vt.plot_traj(px_lin-xmin, py_lin-ymin, tend=itime, linewidth=0.25, alpha=0.5)
    if save_fig:
        fname = save_dir + gen_name + str(itime) + fmt
        plt.savefig(fname)



