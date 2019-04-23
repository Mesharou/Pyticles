
###################################################################################
#Load standard python modules
###################################################################################

import sys

import os

#for netcdf files
#from Scientific.IO.NetCDF import *
from netCDF4 import Dataset

#for numeric functions
import numpy as np

#for numeric functions
import scipy.interpolate as interp

#copy data
from copy import copy

#for plotting
import matplotlib.pyplot as py
import matplotlib.colors as col
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#import streamlines as st
#import streamplot as sp

#for 3D plotting
#from enthought.mayavi import mlab

#shell scripting
import subprocess as sub

#convert seconds to date
import time as tm


#For nested loop
from itertools import product

###################################################################################
#Load custom python modules
###################################################################################


#ROMSTOOLS
import romstools_old as roms

#Simulations (path, data...)
import simulations_old as sim

#Smoothing
import R_smooth as sm

#from guppy import hpy
#hp = hpy()
#hp.setrelheap()


import R_tools as tools







