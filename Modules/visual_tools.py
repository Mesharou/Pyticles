import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.collections import LineCollection
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset


import sys
sys.path.append('../Modules/')
from R_files import load





##############################################################################
def get_var(var, ncfile, **kwargs):
    '''
    Returns the netcdf variable 'var' from ncfile
    var is a string : name of variable
    optionnal: itime time index for slicing
    return np.array py_var
    '''
    nc = Dataset(ncfile, 'r')
    if var in nc.variables:
        if 'itime' in kwargs:
            itime = kwargs['itime']
       	    ndims = len(nc.variables[var][:].shape)
            if ndims == 3:
                py_var = nc.variables[var][itime, :, :]
            elif ndims == 4:
                py_var = nc.variables[var][itime, :, : ,:]
            else: 
                print(f'Error dimension {ndim-1} is not supported')
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
def map_traj(pvar, px, py, map_var=[], ng=2, mask=[], cmap1='bone', 
             cmap2='jet', **kwargs):
    '''
    Plot particles trajectory with pvar in color 
    Superimposed with another 2D variable 

    parameters : 
        pvar variable at particles points
        px, py: particles horizontal pozition
        map_var: 2D variable
        mask : land mask : only if map_var is a masked array
        cmap1/cmap2 : colormap for scatter/map

    keyword arguments:
        coord = [ymin, ymax, xmin, max] : plot axes bounds

    '''
    if 'coord' in kwargs:
        coord = kwargs['coord']
        [ymin, ymax, xmin, xmax] = coord[0:4]
    else:
        xmin = int(np.floor(np.nanmin(px))) - ng
        xmax = int(np.ceil(np.nanmax(px))) + ng
        ymin = int(np.floor(np.nanmin(py))) - ng
        ymax = int(np.ceil(np.nanmax(py))) + ng
    if len(mask) > 0:
        map_var.mask = mask
        map_var.mask = ~map_var.mask 
    if ('vmin1' in kwargs) and ('vmax1' in kwargs):
        vmin1 = kwargs['vmin1']
        vmax1 = kwargs['vmax1']
    else:
        vmin1 = np.nanmin(pvar)
        vmax1 = np.nanmax(pvar)
     
    fig = plt.figure()
    vmin = np.nanmin(pvar)
    vmax = np.nanmax(pvar)
    plt.scatter(px, py, c=pvar, cmap=cmap1, alpha=0.75, vmin=vmin1,
                vmax=vmax1, zorder=10)
    cbar1 = plt.colorbar(orientation='horizontal', shrink=0.25, pad=0.05)
    plt.pcolormesh(ma.masked_invalid(map_var), cmap=cmap2,
                   rasterized=True)
    cbar2 = plt.colorbar(shrink=0.25);
    plt.axis([xmin, xmax, ymin, ymax])
    plt.show()
    return fig
##############################################################################















