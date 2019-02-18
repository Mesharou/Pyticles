#!/usr/bin/env python
'''
Example Script to analyze paticle trajectories
'''
###################################################################################
# Load all useful modules (see /home/gula/Desktop/python_v2/Modules/)
##################################################################################
#import matplotlib
#matplotlib.use('Agg')
import sys

#add the Modules folder in your python PATH
sys.path.append("../Modules/")
sys.path.remove("/home/gula/Desktop/Work_capella/Python/Python_Modules")

import pyticles_3d_sig_sa as partF
import pyticles_sig_sa as part

sys.path.append("/home/gula/Desktop/Work_capella/Python/Python_Modules")

from Modules import *
from Modules_gula import *


#import particules_3d_sa as partF_z
#import particules_sa as part_z

import matplotlib.pyplot as plt

###################################################################################
ncfile = '/net/libra/local/tmp/1/gula/particles/ATLBIG/atlbigsig_Eurofloat_RK4_8_0024.nc' #the netcdf file containing
fifig = './'
###################################################################################



print('Loading file')
nc = Dataset(ncfile, 'r', format='NETCDF3_CLASSIC')
parameters = nc.simulation
base = nc.base

#check if periodic
if nc.x_periodic==1: x_periodic=True
elif nc.x_periodic==0: x_periodic=False

if nc.y_periodic==1: y_periodic=True
elif nc.y_periodic==0: y_periodic=True

ng = nc.ng # ghost points in the horizontal

nq = len(nc.dimensions['nq'])
ntime = len(nc.dimensions['time'])
nc.close()

###################################################################################
# load simulation parameters
###################################################################################

print('Loading simul')
simul = load(simul = parameters)
depths = simul.coord[4]


###################################################################################
# get time
###################################################################################

time = ionetcdf.get(ncfile,'ocean_time',simul)
time[0]=time[1]-(time[2]-time[1])
time_int=list(range(simul.time0,simul.ncname.tend+1,simul.dtime))[:len(time)]

nqmx  = nq

# where
lim1,lim2 = 0,1000

###################################################################################

def load_p(name):  
    px = ionetcdf.get(ncfile,name,simul)
    time_add = 0
    return px[:,lim1:lim2],time_add-lim1

#############
    
def load_time():  
    time = ionetcdf.get(ncfile,'ocean_time',simul)
    time_int=list(range(simul.time0,simul.ncname.tend+1,simul.dtime))[:len(time)]
    #############
    return time[lim1:lim2],time_int[lim1:lim2]

###################################################################################

###########################
# get time
###########################

time,time_int = load_time()
simul.time0=time_int[0]

###########################
# Get px,py
###########################

px,time_add = load_p('px')
py,time_add = load_p('py')
pz,time_add = load_p('pz')
pt,time_add = load_p('pt')
ps,time_add = load_p('ps')


pt[pt==0] = np.nan
ps[np.isnan(pt)] = np.nan

pu,time_add = load_p('pu')
pv,time_add = load_p('pv')

###########################
# Compute pdepth (actual depth) from pz
###########################

print('depth...')
try:
    pdepth,time_add= load_p('pdepth')
    pdepth[np.isnan(pt)] = np.nan
except:   
    pdepth=copy(pz); pdepth[:]=np.nan
    for it in range(pdepth.shape[1]):
        [z_r,z_w] = tools.get_depths(simul)
        pdepth[:,it] = part.map_varw(simul,z_w,px[:,it],py[:,it],pz[:,it])   
    del pz



###########################     
# Compute ptopo (topography at each particule position)
###########################

print('topo...')
ptopo = copy(pdepth); ptopo[:]=np.nan
for it in range(ptopo.shape[1]):
    ptopo[:,it] = part.map_topo(simul,px[:,it],py[:,it])

ptopo[ptopo==0] = np.nan
ptopo[np.isnan(pt)] = np.nan



###########################     
# Compute pmask (mask at each particule position)
###########################

print('mask...')
pmask = copy(pdepth); pmask[:]=np.nan
mymask = copy(simul.mask)
mymask[np.isnan(mymask)] = 0.

for it in range(pmask.shape[1]):
    pmask[:,it] = part.map_var2d(simul,mymask,px[:,it],py[:,it])


###################################################################################

def plot_selection(filetime,ipart):
    simul.update(time_int[-1]); coord = part.subsection(px[:,-1],py[:,-1],offset=10); 
    salt = var('temp',simul,coord=coord,depths=[whichdepth]).data
    plt.contourf(salt.T*simul.mask[coord[2]:coord[3],coord[0]:coord[1]].T,100); plt.colorbar(); #plt.plot(px[:,0]-np.nanmin(px[:,0]),py[:,0]-np.nanmin(py[:,0]),'d');
    plt.plot(px[ipart,-1]-coord[2]-base,py[ipart,-1]-coord[0]-base,'o', markersize=5, markerfacecolor='white')
    del salt
    plt.title('depth = ' + format(whichdepth)); 
    plt.savefig(fifig + 'selection.png')

###################################################################################

def plot_part(filetime, ipart):
    exec(compile(open('subplot_part.py').read(), 'subplot_part.py', 'exec'))

###################################################################################
    
def make_plot(which_plot,filetime=0,ipart=[0]):
    proc=mp.Process(target=which_plot, args=(filetime,ipart))
    proc.daemon=True; proc.start(); 
    print(proc, proc.is_alive())
    proc.join()
    print(proc, proc.is_alive())

###################################################################################

def nanint(a,int_val=0):
    f = np.zeros((a.shape),dtype=int)
    for i in range(len(f)):
        if np.isnan(a[i]):
            f[i]=int_val
        else:
            f[i]=np.int_(a[i])      
    return f
    
###################################################################################
# Define the mask (which particules to plot)
###################################################################################

mask=[];
whichdepth = 0

#for iq in range(nq):
#    if np.nanmin(ptopo[iq,:])<55.:
#        mask.append(iq)

for iq in range(nq):
  if not np.isnan(px[iq,-1]*py[iq,-1]):
    if (np.abs(px[iq,-1]-px[iq,-3])<1e-2) and (np.abs(py[iq,-1]-py[iq,-3])<1e-2):
        mask.append(iq)

print('nb of particules in mask is' , len(mask))

make_plot(plot_selection,0,mask)


###################################################################################
# Plot things
###################################################################################

#name used to save figures
which = 'part_' +  format(len(mask)) + '_allmask'


for iq in mask:
  for filetime in time_int:
    make_plot(plot_part,filetime,[iq])

















    



