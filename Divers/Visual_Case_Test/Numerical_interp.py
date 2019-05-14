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

nc_linear = '/home/jeremy/Bureau/Data/Pyticles/Linear_interp/' \
          + 'Case_1_Linear_interp_6_1510.nc'
nc_cubic = '/home/jeremy/Bureau/Data/Pyticles/Cubic_interp/' \
          + 'Case_1_Cubic_interp_1_1510.nc'

roms_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_his.1550.nc'
grd_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_grd.nc'

#############################################################################
itime = 14

px_lin = vt.get_var('px', nc_linear, itime=itime)
px_cub = vt.get_var('px', nc_cubic, itime=itime)

num_bins = 50

fig, ax = plt.subplots()

# the histogram of the data
n, bins, patches = ax.hist(px_cub-px_lin, num_bins, normed=1)

# add a 'best fit' line
#y = mlab.normpdf(bins, mu, sigma)
#ax.plot(bins, y, '--')
ax.set_xlabel('Smarts')
ax.set_ylabel('Probability density')
ax.set_title(r'Histogram of IQ: $\mu=100$, $\sigma=15$')

# Tweak spacing to prevent clipping of ylabel
fig.tight_layout()
plt.show()





