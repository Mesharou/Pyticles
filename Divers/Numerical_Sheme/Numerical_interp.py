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

save_fig = False
save_dir = '/home/jeremy/Bureau/Data/Pyticles/RESU/Visual_tools/' \
         + 'Numerical_Schemes/'
gen_name = 'fig_num_scheme_'
fmt = '.png'

nc_linear = '/home/jeremy/Bureau/Data/Pyticles/Linear_adv/' \
          + 'Case_1_Linear_adv_4_1510.nc'
nc_cubic = '/home/jeremy/Bureau/Data/Pyticles/Cubic_adv/' \
          + 'Case_1_Cubic_adv_4_1510.nc'

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
pdepth_lin = vt.get_var('pdepth', nc_linear)
pdepth_cub = vt.get_var('pdepth', nc_cubic)

# comparing cubic and old linear
np.std(pdepth_lin - pdepth_cub)
np.std(pt_lin - pt_cub)


num_bins = 55

fig, ax = plt.subplots()

# the histogram of the data
#n, bins, patches = ax.hist(px_cub-px_lin, num_bins, density=1)
n, bins, patches = ax.hist(pdepth_cub-pdepth_lin, num_bins, density=1,
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

#####################
# comparing old and new linear
'''
num_bins = 55

fig, ax = plt.subplots()

# the histogram of the data
#n, bins, patches = ax.hist(px_cub-px_lin, num_bins, density=1)
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
############################################################
# plot mean(var_cub - var_std) +-  std(var_cub - var_lin)(t)
'''
First we create 2 large datasets var_cub, var_lin shape (nq, nt)

then we plot mean(var_cub - var_lin)(t) in red
superimposed with std in shaded black

'''
nq = 1000
nt = 100
var_cub = np.random.rand(nq, nt) + 1/1000*np.arange(100)
var_lin = np.random.rand(nq, nt)

fig = plt.figure
plt.plot(np.mean(var_cub - var_lin, axis=0))
plt.plot(np.mean(var_cub - var_lin, axis=0) + np.std(var_cub-var_lin, axis=0),
        'r')
plt.plot(np.mean(var_cub - var_lin, axis=0) - np.std(var_cub-var_lin, axis=0), 
        'r')

plt.show()
