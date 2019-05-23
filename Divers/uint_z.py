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
import visual_tools as vt 

##############################################################################

def rho2u_3d(var_rho):

    var_u = 0.5*(var_rho[1:,:,:]+var_rho[:-1,:,:])

    return var_u

##############################################################################

def rho2u_2d(var_rho):
    var_u = 0.5*(var_rho[1:,:]+var_rho[:-1,:])
    return var_u


#######################################################
#Transfert a field at rho points to v points
#######################################################

def rho2v(var_rho):

    if np.ndim(var_rho)==1:
        var_v = 0.5*(var_rho[1:]+var_rho[:-1])
    elif np.ndim(var_rho)==2:
        var_v = rho2v_2d(var_rho)
    else:
        var_v = rho2v_3d(var_rho)

    return var_v


##############################

def periodize2d_fromnc(simul,variable,coord,x_periodic=False,y_periodic=False,ng=0):

    [ny1tot,ny2tot,nx1tot,nx2tot] = simul.coord[0:4]
    [ny1,ny2,nx1,nx2] = coord[0:4]

    nc = Dataset(simul.ncfile, 'r')

    ################################
    mask = copy(simul.mask)
    mask[np.isnan(mask)]=0
    ################################

    nw = min(ng,nx1); ne = min(ng,nx2tot-nx1tot+2*ng-nx2)
    ns = min(ng,ny1); nn = min(ng,ny2tot-ny1tot+2*ng-ny2)

    myvar = np.zeros((nx2-nx1,ny2-ny1))*np.nan

    myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn] = simul.Forder(np.squeeze(nc.variables[variable][simul.infiletime,ny1-ns:ny2-2*ng+nn,nx1-nw:nx2-2*ng+ne]))

    myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn] = (myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn].T * (mask[nx1-nw:nx2-2*ng+ne,ny1-ns:ny2-2*ng+nn]).T).T

    if nw<ng and x_periodic:
        for i in range(ng):
            myvar[ng-nw-1-i,ng-ns:ny2-ny1-ng+nn] = simul.Forder(np.squeeze(nc.variables[variable][simul.infiletime,ny1-ns:ny2-2*ng+nn,nx2tot-1-i]))
        nw=ng

    if ne<ng and x_periodic:
        for i in range(ng):
            myvar[nx2-nx1-ng+ne+i,ng-ns:ny2-ny1-ng+nn] = simul.Forder(np.squeeze(nc.variables[variable][simul.infiletime,ny1-ns:ny2-2*ng+nn,nx1tot+i]))
        ne=ng

    if ns<ng and y_periodic:
        for i in range(1,ng):
            myvar[ng-nw:nx2-nx1-ng+ne,ng-ns-1-i] = simul.Forder(np.squeeze(nc.variables[variable][simul.infiletime,ny2tot-1-i,nx1-nw:nx2-2*ng+ne]))

    if nn<ng and y_periodic:
        for i in range(1,ng):
            myvar[ng-nw:nx2-nx1-ng+ne,ny2-ny1-ng+nn+i] = simul.Forder(np.squeeze(nc.variables[variable][simul.infiletime,ny1tot+i,nx1-nw:nx2-2*ng+ne]))

    return myvar

##############################################################################

def periodize2d_fromvar(simul,var2d,coord,x_periodic=False,y_periodic=False,ng=0):

    [ny1tot,ny2tot,nx1tot,nx2tot] = simul.coord[0:4]
    [ny1,ny2,nx1,nx2] = coord[0:4]

    nw = min(ng,nx1); ne = min(ng,nx2tot-nx1tot+2*ng-nx2)
    ns = min(ng,ny1); nn = min(ng,ny2tot-ny1tot+2*ng-ny2)

    myvar = np.zeros((nx2-nx1,ny2-ny1))*np.nan

    myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn] = var2d[nx1-nw:nx2-2*ng+ne,ny1-ns:ny2-2*ng+nn]

    if nw<ng and x_periodic:
        for i in range(ng):
            myvar[ng-nw-1-i,ng-ns:ny2-ny1-ng+nn] = var2d[nx2tot-1-i,ny1-ns:ny2-2*ng+nn]
        nw=ng

    if ne<ng and x_periodic:
        for i in range(ng):
            myvar[nx2-nx1-ng+ne+i,ng-ns:ny2-ny1-ng+nn] = var2d[nx1tot+i,ny1-ns:ny2-2*ng+nn]
        ne=ng

    if ns<ng and y_periodic:
        for i in range(1,ng):
            myvar[ng-nw:nx2-nx1-ng+ne,ng-ns-1-i] = var2d[nx1-nw:nx2-2*ng+ne,ny2tot-1-i]

    if nn<ng and y_periodic:
        for i in range(1,ng):
            myvar[ng-nw:nx2-nx1-ng+ne,ny2-ny1-ng+nn+i] = var2d[nx1-nw:nx2-2*ng+ne,ny1tot+i]

    return myvar

#################################################
# get_depths (from setdepth.F in romsucla)
#################################################


def get_depths(simul,x_periodic=False,y_periodic=False,ng=0,**kwargs):


    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else:
        coord = simul.coord

    [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
    [ny1,ny2,nx1,nx2] = simul.coord[0:4]
    #topo = np.asfortranarray(simul.topo[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
    topo = periodize2d_fromvar(simul,simul.topo,coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)

    if hasattr(simul, 'zeta'):
        zeta=periodize2d_fromvar(simul,simul.zeta,coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)
    else:
        zeta=periodize2d_fromnc(simul,'zeta',coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)

    if simul.ncname.model=='ucla':
        (z_r,z_w) = partF.zlevs(topo, zeta, simul.hc, simul.Cs_r, simul.Cs_w)
    else:
        if simul.VertCoordType == 'NEW':
            print('using NEW_S_COORD')
            (z_r,z_w) = partF.zlevs_croco_new(topo, zeta, simul.hc, simul.Cs_r, simul.Cs_w, simul.sc_r, simul.sc_w)
        else:
            print('using OLD_S_COORD')
            (z_r,z_w) = partF.zlevs_croco_old(topo, zeta, simul.hc, simul.Cs_r, simul.Cs_w, simul.sc_r, simul.sc_w)


    return [z_r,z_w]
#############################################################################

def interp(simul, varname, coord, advdepth, z_thick):

    coordz = copy(coord)

    if varname=='u':
        imin = 1
        jmin = 0
        coordz[3] += 1
    elif varname=='v':
        imin = 0
        jmin = 1
        coordz[1] += 1

    [z_r, z_w] = get_depths(simul, coord=coordz)
    var0  = simul.Forder(  \
                        np.squeeze(nc.variables[varname][simul.infiletime, :,
                                   coord[0]:coord[1], coord[2]:coord[3]]))
    if varname == 'u':
        h = rho2u_3d(np.max(np.min(z_w, advdepth+z_thick), advdepth-z_thick) \
             - np.min(np.max(z_w, advdepth-z_thick), advdepth+z_thick))
    elif varname == 'v':
        h = rho2v_3d(np.max(np.min(z_w, advdepth+z_thick), advdepth-z_thick) \
            - np.min(np.max(z_w, advdepth-z_thick), advdepth+z_thick))

    varz = np.sum(h*var0, axis = 0)

    return varz
##############################################################################
z_thick = 100.
advdepth = -200

###### MAIN ######
timing = True
ng = 1
coord = simul.coord[0:4]
[ny1tot, ny2tot, nx1tot, nx2tot] = simul.coord[0:4]

if timing: tstart2 = tm.time()

nc = Dataset(simul.ncfile, 'r')
[ny1, ny2, nx1, nx2] = coord

################################

mask = copy(simul.mask)
mask[np.isnan(mask)] = 0

################################

u = np.zeros((nx2-nx1-1, ny2-ny1))*np.nan
v = np.zeros((nx2-nx1, ny2-ny1-1))*np.nan

nw = min(ng,nx1); ne = min(ng,nx2tot-nx1tot+2*ng-nx2)
ns = min(ng,ny1); nn = min(ng,ny2tot-ny1tot+2*ng-ny2)

if x_periodic:
    iper=1
else:
    iper = 0
###########################################################################
# Works but do not take into acccount whether z_wu is above or below
# critical z values : advdepth +/- z_thick/2
varname = 'u'
coordz = copy(coord)
[z_r, z_w] = get_depths(simul, coord=coordz)
z_wu = rho2u_3d(z_w)
z_ru = rho2u_3d(z_r)
h_u = z_wu[:, :, 1:] - z_wu[:, :, :-1]
index = (z_ru <= advdepth + z_thick/2) & (z_ru >= advdepth - z_thick/2)
h_u[~index] = 0
var0  = simul.Forder(  \
                        np.squeeze(nc.variables[varname][simul.infiletime, :,
                                   coord[0]:coord[1], coord[2]:coord[3]]))



ubar = np.sum( h_u * var0, axis=-1) / z_thick
plt.contourf(ubar.T)
plt.clim(-1, 1)
plt.colorbar()
plt.show()
############################################################ 
# Works and takes into account issue related to outrange    
del h_u, index, ubar

indexp = z_wu > advdepth + z_thick/2.
indexm = z_wu < advdepth - z_thick/2
z_wu[indexp] = advdepth + z_thick/2.
z_wu[indexm] = advdepth - z_thick/2.
h_u = z_wu[:, :, 1:] - z_wu[:, :, :-1]
ubar = np.sum( h_u * var0, axis=-1) / z_thick

plt.contourf(ubar.T)
plt.clim(-1, 1)
plt.colorbar()
plt.show()
##############################################################################
# TEST using an analytical velocity field
[nx, ny, nz] = var0.shape
u_ana = np.ndarray(var0.shape) + 1
u_ana = (np.ndarray(var0.shape) + 1) * np.arange(nz)

u_anabar = np.sum( h_u * u_ana, axis=-1) / z_thick

plt.contourf(u_anabar.T)
#plt.clim(-1, 1)
plt.colorbar()
plt.show()



plt.contour(u_ana[1000, :, :])
plt.colorbar()
plt.show()
#################################################################################
'''
Issue found in case where ptopo < -advdepth + z_thick/2
i.e bottom of water column intersects with bathymetry

at the moment ubar/zthick however it should be divided by
-advdepth + zthick/2 - ptopo

We can check usinh analytical velocity field u_ana = 0*np.ndarray(var0.shape) + 1
'''

advdepth = -200
z_thick = 50

u_ana = 0*np.ndarray(var0.shape) + 1
varname = 'u'
coordz = copy(coord)
if varname=='u':
    imin = 1
    jmin = 0
    coordz[3] += 1
elif varname=='v':
    imin = 0
    jmin = 1
    coordz[1] += 1

[z_r, z_w] = get_depths(simul, coord=coordz)
if varname == 'u':
    z_wvar = rho2u_3d(z_w)
elif varname == 'v':
    z_wvar = rho2v(z_w)

plt.figure
plt.subplot(2,1,1)
plt.plot(z_wvar[5,5,:])
plt.subplot(2,1,2)
plt.plot(z_w[5,5,:])
plt.plot(z_w[6,5,:])
plt.show()

# indexp and indexm are here to truncate and therefore ensure the water column
# thickeness does not exceed z_thick when multiplied to var
indexp = z_wvar > advdepth + z_thick/2.
indexm = z_wvar < advdepth - z_thick/2.
z_wvar[indexp] = advdepth + z_thick/2.
z_wvar[indexm] = advdepth - z_thick/2.
h_var = z_wvar[:, :, 1:] - z_wvar[:, :, :-1]
#var_bar = np.sum(h_var * u_ana, axis=-1) / z_thick
var_bar = np.sum(h_var * u_ana, axis=-1) / np.sum(h_var, axis=2)
var_bar[np.isnan(var_bar)] = 0.


# topo at u points
topo = simul.topo
u_topo = rho2u_2d(topo)


#var_bar_ok = np.sum(h_var * u_ana, axis=-1) / (-advdepth + z_thick/2 -u_topo)
'''
# indexes where topo is below, in between and above water column for average
# These are used in case where topography intersects advdepth - zthick/2
# In this region we have to divide
ind_low = u_topo > -advdepth - z_thick/2
ind_mid = (u_topo <= -advdepth - z_thick/2) & (u_topo <= -advdepth +  z_thick/2)  
ind_upp = u_topo < -advdepth + z_thick/2

plt.figure
plt.subplot(311)
plt.pcolormesh(ind_low.T)
plt.colorbar()
plt.subplot(312)
plt.pcolormesh(ind_mid.T)
plt.colorbar()
plt.subplot(313)
plt.pcolormesh(ind_upp.T)
plt.colorbar()
plt.show()


h_thick = np.ndarray(u_topo.shape)
h_thick[ind_low] = z_thick
h_thick[ind_mid] = advdepth + z_thick/2 - u_topo[ind_mid]
h_thick[ind_upp] =  
'''




plt.figure
plt.subplot(211)
plt.pcolormesh(var_bar.T)
plt.colorbar()
plt.subplot(212)
plt.pcolormesh(topo.T)
plt.colorbar()

plt.show()
plt.close('all')


plt.figure
plt.subplot(211)
plt.pcolormesh(var_bar_ok.T)
plt.colorbar()
plt.subplot(212)
plt.pcolormesh(np.sum(h_var, axis=2).T)
plt.colorbar()

plt.show()
plt.close('all')

###########################################################################
'''
Checking whether it works or not...

We used a simulation with advdepth = -200.; z_thick = 100.

particles were seeded close to -200 m bathymetric contour in order to have 
some in middle of water column

'''

################################################################################
# ROMS INPUTS
################################################################################

# if meanflow = True Roms data are not updated (used for climatology)
meanflow = False
# in case of periodic channel
x_periodic = False
y_periodic = False

# dfile is frequency for the use of the ROMS outputs
# (default is 1 = using all outputs files)
dfile = 1
start_file = 1510
end_file = 1535

###### Restart from a Pyticles output file
# user should not change start_file
# restart_time : number of time step since start_file
restart = False
restart_time = 7 #nb of time steps in the restart_file
restart_file = '/home/jeremy/Bureau/Data/Pyticles/' \
               +'/Linear_interp/Case_1_Linear_interp_6_1510.nc'

if not restart:
    restart_time = 0
else:
    start_file += restart_time

# Load simulation
# parameters = my_simul + [0,nx,0,ny,[1,nz,1]] ; nx, ny, nz Roms domain's shape 
my_simul = 'Case_1'
parameters = my_simul + ' [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)
##############################################################################
py_file = '/home/jeremy/Bureau/Data/Pyticles/Zavg_v2.0/Case_1_Zavg_v2.0_adv200.0m_1_1510.nc'
roms_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_his.1510.nc'

pu = vt.get_var('pu', py_file, itime=0).data
pv = vt.get_var('pv', py_file, itime=0).data
px = vt.get_var('px', py_file, itime=0).data
py = vt.get_var('py', py_file, itime=0).data
pz = vt.get_var('pz', py_file, itime=0).data
ptopo = vt.get_var('ptopo', py_file, itime=0)

u = vt.get_var('u', roms_file, itime=0, ndims=4).data 
v = vt.get_var('v', roms_file, itime=0, ndims=4)


indx = ptopo < -advdepth - z_thick/2

[z_r, z_w] = get_depths(simul, coord=simul.coord)

####################################################
ip = -9
px0 = px[ip]
py0 = py[ip]

xu = np.int(np.floor(px0))
yu = np.int(np.floor(py0))

pu_test = pu[ip]

plt.figure
plt.plot(z_r[xu, yu, :], u[:, yu, xu]-pu_test)
plt.plot(z_r[xu, yu+1, :], u[:, yu+1, xu]-pu_test)

plt.show()













