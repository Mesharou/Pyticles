import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as col

from netCDF4 import Dataset

import numpy as np
import sys
from scipy.interpolate import interp1d



sys.path.append("../Modules/")
sys.path.append("Python_Modules/")
import seeding_part
from R_files import load
import pyticles_sig_sa as part
import pyticles_3d_sig_sa as partF

from R_tools import rho1_eos, rho_eos

ncfile = '/home/jeremy/Bureau/Data/Pyticles/Rho1_Seed/Case_1_Rho1_Seed_12_1550.nc'
roms_file = '/home/jeremy/Bureau/Data/Pyticles/chaba_his.1550.nc'

#############################################################################

def get_var(var, ncfile, it=0):
    '''
    Returns the netcdf variable 'var' from ncfile
    var is a string : name of pyticiles  variable
    return np.array py_var
    use module netCDF4
    use module numpy as np
    it is used for time slicing 
    '''
    nc = Dataset(ncfile, 'r')
    if var in nc.variables:
        if var in ['u', 'v']:
            py_var = nc.variables[var][it,:,:,:]
        else:
            py_var = nc.variables[var][:]
    else:
        py_var = []
        print(f'Error {var} is not found in file {ncfile}')
    nc.close()
    return py_var
#############################################################################

def get_attr(attr, ncfile, it=0):
    '''
    Returns the netcdf attribute 'attr' from ncfile
    attr is a string : name of pyticilesgolabal attribute
    '''
    nc = Dataset(ncfile, 'r')
    if attr in nc.ncattrs():
            py_attr = getattr(nc, attr)
    else:
        py_attr= []
        print(f'Error {attr} is not found in file {ncfile}')
    nc.close()
    return py_attr



##############################################################################

##############################################################################
def ini_surf(mask, simul, surf0, x, y, z, rho, ng=0):
    for k in range(len(surf0)):
        for i in range(x.shape[2]):
            for j in range(x.shape[1]):
                if mask[i,j]==1:
                    f = interp1d(rho[:, j, i], list(range(rho.shape[0])), kind='cubic')
                    z[k, j, i] = f(surf0[k])
                else:
                    z[k, j, i] = 0.
    return z
##############################################################################



start_file = 1550
parameters = 'Case_1 [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

x_periodic = get_attr('x_periodic', ncfile)
y_periodic = get_attr('y_periodic', ncfile)
ng = get_attr('ng', ncfile)

#########
# box parameters

ic = 400
jc = 200
lev0 = 0
lev1 = 50
iwd = 1
jwd = 20
nx = simul.coord[3]
ny = simul.coord[1]
nz = len(simul.coord[4])
nny = 1
nnx = 1
nnlev = 1
surf0 = [1025] #- simul.rho0
mask = simul.mask

lev1 = len(surf0) - 1

##### First box in order to interpolate sigma onto rho levels
lev1 = 50 # Needed to get all levels
z, y, x = seeding_part.seed_box(ic=ic, jc=jc, lev0=lev0, lev1=lev1, nnx=nnx,
                                nny=nny, iwd=iwd, jwd=jwd, nx=nx, ny=ny)
#i0 = np.int(np.floor(np.min(x))) - 2
#j0 = np.int(np.floor(np.min(y))) -2 
k0 = 0
nq = len(x.reshape(-1))
i0 = 0
j0 = 0

## Long to process no need to reload this
[temp, salt] = part.get_ts_io(simul, x_periodic = x_periodic,
                y_periodic = y_periodic, ng=ng)
[z_r, z_w] = part.get_depths(simul)
                       
#####
# Good algorithm : takes rho1 from R_tools using rho0 from roms model
# Take temp, salt, z_r, z_w
# interpolate rho1 onto sub_box : x, y, z
# Compute surface vertical level at surf0 values ex 
# 
# rho 1 density perturbation
roms_rho0 = simul.rho0
roms_rho1 = rho1_eos(temp, salt, z_r, z_w, roms_rho0)
map_rho1 = (part.map_var(simul, roms_rho1, x.reshape(-1), y.reshape(-1),
                        z.reshape(-1), ng=ng)).reshape(x.shape)
## rho density anomly
rho0 = surf0
rho = rho_eos(temp, salt , z_r, z_w, rho0)
map_rho = (part.map_var(simul, rho, x.reshape(-1), y.reshape(-1),
                        z.reshape(-1), ng=ng)).reshape(x.shape)


###
fig = plt.figure

ax1 = plt.subplot(211)
ax1.contour(roms_rho1[200:-1, 600, :].T, [-1.5], linewidths=2.)
im1 = ax1.contourf(roms_rho1[200:-1, 600, :].T)

cbar1 = plt.colorbar(im1, ax=ax1)
ax1.set_xlabel('X (grid points)')
plt.title(r'density anomaly \rho')

ax2 = plt.subplot(212)
ax2.contour(roms_rho1[500, :1000, :].T, [-1.5], linewidths=2.)

im2 = ax2.contourf(roms_rho1[500, :1000, :].T)
cbar2 = plt.colorbar(im2, ax=ax2)
ax2.set_xlabel('Y (grid points)')

plt.tight_layout()
plt.show()
###
fig = plt.figure
ip = 2
plt.subplot(121)
plt.contourf(map_rho[:, :, ip])
cbar = plt.colorbar()

plt.subplot(122)
plt.contourf(rho[ic - iwd + ip +1, jc-jwd:jc+jwd, :].T)
cbar = plt.colorbar()

plt.show()

del x, y, z

#####
# Getting iso_surface
lev1 = len(surf0) - 1
z, y, x = seeding_part.seed_box(ic=ic, jc=jc, lev0=lev0, lev1=lev1, nnx=nnx,
                                nny=nny, iwd=iwd, jwd=jwd, nx=nx, ny=ny)
z = ini_surf(mask, simul, surf0, x, y, z, map_rho1, ng=ng)
zrho = ini_surf(mask, simul, 0, x, y, z, map_rho, ng=ng)
########
# diags
#
plt.contour(roms_rho1[ic+200,0:1000,:].T , levels=10)
cbar = plt.colorbar()
plt.scatter(y,z)
plt.show()




plt.subplot(211)
plt.contour(roms_rho1[ic+200,0:1000,:].T, levels=10)
cbar = plt.colorbar()
plt.title(f'rho1 with rho0 = {roms_rho0}')
plt.ylabel('sigma levels')

plt.subplot(212)
plt.contour(rho[ic+200,0:1000,:].T, levels=10)
cbar = plt.colorbar()
plt.ylabel('sigma levels')
plt.xlabel('Y(km)')
plt.title('rho')
#plt.caxis([-1,  1])
plt.show()

##########
plt.contourf(rho[200, 250:1000 , :].T)
plt.colorbar()
plt.show()




############
plt.hist(rho_check - rho0)
plt.show()

plt.contourf(z[0,:,:])
cbar = plt.colorbar()
plt.show()
     
plt.plot(z_r[ic, jc, :])
plt.grid(True)
plt.show()
toto

plt.figure
plt.contourf(new_rho[:, :, 10], 10, cmap=plt.cm.bone)
cbar = plt.colorbar()
plt.show
