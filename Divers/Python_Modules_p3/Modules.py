'''

To be loaded at the beginning of a python script:

from Modules import *

'''
 
 
 
 
###################################################################################
#Load standard python modules
###################################################################################
#del sys.modules["romstools"]


import os,sys

#for netcdf files
#from Scientific.IO.NetCDF import *
from netCDF4 import Dataset

#import dask.array as da

'''
Switching to netcdf4 module 

files needs to be open with:
    ncfile  = Dataset(ncname, 'r', format='NETCDF3_CLASSIC')
instead of
    ncfile = NetCDFFile(ncname,'r')
(http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html)
'''


#for numeric functions
import numpy as np
import numpy.ma as ma # for masked arrays
import scipy.interpolate as interp
import scipy.integrate as integrate

#convert seconds to date
import time as tm

#copy data
from copy import copy

#for plotting
import matplotlib.pyplot as plt
import matplotlib.pyplot as py # backward compatibility
import matplotlib.colors as col
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter

#py.rcParams['text.usetex']=True
#py.rcParams['text.latex.unicode']=True

#For nested loop
from itertools import product


#For multiprocessing
import multiprocessing as mp


###################################################################################
#Load custom python modules
###################################################################################

# some ROMS tools written in python
import R_tools as tools

# some ROMS tools written in .F
import R_tools_fort as toolsF
#del sys.modules["R_tools_fort"]

# Compute variables
from R_vars import var


# some ROMS tools written in python
#import R_tools_libra as tools

# Compute variables
#from R_vars_libra import var

# Simulations (path, data...)
from R_files import load

# netcdf io related functions
from R_netcdf import ionetcdf

# Plotting functions
from R_plot import plot

# Smoothing and interpolation functions
import R_smooth as sm

# Statistics functions
#import R_stats as st

###################################################################################
# Define machine independant path
###################################################################################

user = os.getenv('USER') + '/'

celtic='/mnt/celtic/' + user 
avatar='/mnt/avatar/'  + user       
shirma='/mnt/shirma/' + user 
cherokee='/mnt/cherokee/' + user 
comanche='/mnt/comanche/' + user 
kamiya='/mnt/kamiya/' + user 
goth='/mnt/goth/' + user

if os.getenv('HOSTNAME') =='shirma.atmos.ucla.edu':
    shirma='/shirma/' + user 
elif  os.getenv('HOSTNAME')=='celtic.atmos.ucla.edu':
    celtic='/celtic/' + user 
elif  os.getenv('HOSTNAME')=='avatar.atmos.ucla.edu':
    avatar='/avatar/' + user 
elif  os.getenv('HOSTNAME')=='cherokee.atmos.ucla.edu':
    cherokee='/cherokee/' + user
elif  os.getenv('HOSTNAME')=='goth.atmos.ucla.edu':
    goth='/goth/' + user


libra = '/net/libra/local/tmp/1/' + user

###################################################################################
















