'''
JCollin 05-2019
Issue in Modules/pyticles_sig_sa.py
    def get_vel_io
    Index Error occur when computing w

'''
##############################################################################
import numpy as np
import sys
from netCDF4 import Dataset
#copy data
from copy import copy

#ROMSTOOLS
#import R_tools_fort as toolsF

#Simulations (path, data...)
#import R_vars as va

#for plotting
import matplotlib.pyplot as plt

import time as tm
sys.path.append('../Inputs/')
sys.path.append('../Modules/')
from input_file import *
import pyticles_3d_sig_sa as partF

import pyticles_sig_sa as part

##############################################################################
roms_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_his.1550.nc'
py_file = '/home/jeremy/Bureau/Data/Pyticles/Rho1_-1.5/' \
         + 'Case_1_Rho1_-1.5_6_1510.nc'
start_file = 1550
parameters = 'Case_1 [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

##############################################################################

mask = copy(simul.mask)
print(mask.shape)

ng = 1
coord = simul.coord
[z_r, z_w] = part.get_depths(simul, coord=coord, x_periodic=x_periodic,
                               y_periodic=y_periodic, ng=ng)

pm = simul.pm
