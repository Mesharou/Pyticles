import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as col

from netCDF4 import Dataset

import numpy as np
import sys
from scipy.interpolate import interp1d



sys.path.append("../Modules/")
import seeding_part
from R_files import load
import pyticles_sig_sa as part
import pyticles_3d_sig_sa as partF


ptemp = 3
psalt = 35.5
Z = -5000 
rho = 1050.3639165364

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

start_file = 1550
parameters = 'Case_1 [0,10000,0,10000,[1,100,1]] '+ format(start_file)
simul = load(simul = parameters, floattype=np.float64)

x_periodic = get_attr('x_periodic', ncfile)
y_periodic = get_attr('y_periodic', ncfile)
ng = get_attr('ng', ncfile)

[temp, salt] = part.get_ts_io(simul, x_periodic = x_periodic,
                y_periodic = y_periodic, ng=ng)

[z_r, z_w] = part.get_depths(simul)

rho = seeding_part.prho(ptemp=temp, psalt=salt, pdepth=z_r)

# box parameters

ic = 10
jc = 20
lev0 = 0
lev1 = 51
iwd = 5
jwd = 2
nx = rho.shape[0]
ny = rho.shape[1]
nz = rho.shape[2]
nny = 1
nnx = 1
nnlev = 1

z, y, x = seeding_part.seed_box(ic=ic, jc=jc, lev0=lev0, lev1=lev1, nnx=nnx,
                                nny=nny, iwd=iwd, jwd=jwd, nx=nx, ny=ny)
i0 = np.int(np.floor(np.min(x)))
j0 = np.int(np.floor(np.min(y)))
k0 = 0
nq = len(x.reshape(-1))

remap_rho = np.ndarray(x.shape)

ip = 0
print('------ENTERING LOOP-------------------')
#for k in range(x.shape[0]):
#    for j in range(x.shape[1]):
#        for i in range(x.shape[2]):
#            remap_rho[i, j, k] = partF.interp_3d(x[i, j ,k], y[i, j, k],
#                                                   z[i, j ,k], rho, ng, nq,
#                                                   i0, j0, k0)
#            ip += 1
#            print(f'{ip}')



remap_rho  = partF.interp_3d(x.reshape(-1), y.reshape(-1), z.reshape(-1), rho,
                            ng, nq,i0, j0, k0)
new_rho = remap_rho.reshape(x.shape)
mask = simul.mask
rho0 = [1209]

del x, y, z

lev1 = len(rho0) - 1
z, y, x = seeding_part.seed_box(ic=ic, jc=jc, lev0=lev0, lev1=lev1, nnx=nnx,
                                nny=nny, iwd=iwd, jwd=jwd, nx=nx, ny=ny)

z = seeding_part.ini_depth(mask, simul, rho0, x, y, z, new_rho, ng=ng)

##############################################################################
def ini_surf(maskrho, simul, rho0, x, y, z, rho, ng=0):
    for k in range(len(rho0)):
        for i in range(x.shape[2]):
            for j in range(x.shape[1]):
                ix = np.int(np.floor(x[k, j, i])) + ng - 5
                iy = np.int(np.floor(y[k, j, i])) + ng - 18
                ix = i
                iy = j
                print(f'ix = {ix}')
                print(f'iy = {iy}')
                if maskrho[ix,iy]==1:
                    print(f' rho = {rho[:, iy, ix]}')
                    print(f'rho.shape')
                    f = interp1d(rho[:, iy, ix], list(range(rho.shape[0])), kind='cubic')
                    z[k, j, i] = f(rho0[k])
                else:
                    z[k, j, i] = 0.
    return z
##############################################################################

z = ini_surf(mask, simul, rho0, x, y, z, new_rho, ng=ng)


     


toto

plt.figure
plt.contourf(new_rho[:, :, 10], 10, cmap=plt.cm.bone)
cbar = plt.colorbar()
plt.show
