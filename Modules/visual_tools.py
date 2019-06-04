import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.collections import LineCollection
import matplotlib.colors as col

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

    keyword argument : ndims, number of dimension of variable in netcdf file
                       ex: 4 for w(t,z,y,x)
                       Much faster when itime is defined
    return np.array py_var
    '''
    nc = Dataset(ncfile, 'r')
    if var in nc.variables:
        if 'itime' in kwargs:
            itime = kwargs['itime']
       	    if 'ndims' not in kwargs:
                ndims = len(nc.variables[var][:].shape)
            else:
                ndims = kwargs['ndims']

            if ndims == 2:
                py_var = nc.variables[var][itime, :]
            elif ndims == 3:
                py_var = nc.variables[var][itime, :, :]
            elif ndims == 4:
                py_var = nc.variables[var][itime, :, : ,:]
            else: 
                print(f'Error dimension {ndims-1} is not supported')
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
def map_traj(pvar, px, py, map_var=[], method='pcolormesh', ng=2, mask=[],
             cmap1='bone', cmap2='jet', **kwargs):
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
        vmin1, vmax1 : colorbar extent for pvar
		vmin2, vmax2 : colorbar extent for map_var
		title
		xlabel
		ylabel

    '''
	# 2D plot extents
    if 'coord' in kwargs:
        coord = kwargs['coord']
        [ymin, ymax, xmin, xmax] = coord[0:4]
    else:
        xmin = int(np.floor(np.nanmin(px))) - ng
        xmax = int(np.ceil(np.nanmax(px))) + ng
        ymin = int(np.floor(np.nanmin(py))) - ng
        ymax = int(np.ceil(np.nanmax(py))) + ng
    # Mask lands
    if len(mask) > 0:
	    # in case where map_var is not a masked array raise an exception ?
        map_var.mask = mask
        map_var.mask = ~map_var.mask
	# pvar color scale	
    if ('vmin1' in kwargs) and ('vmax1' in kwargs):
        vmin1 = kwargs['vmin1']
        vmax1 = kwargs['vmax1']
    else:
        vmin1 = np.nanmin(pvar)
        vmax1 = np.nanmax(pvar)
	# map_var color scale
    if ('vmin2' in kwargs) and ('vmax2' in kwargs):
        vmin2 = kwargs['vmin2']
        vmax2 = kwargs['vmax2']
    else:
        vmin2 = np.nanmin(map_var)
        vmax2 = np.nanmax(map_var)
    # Text
    if 'title' in kwargs:
	    title = kwargs['title']
    else:
        title = ''
    if 'xlabel' in kwargs:
        xlabel = kwargs['xlabel']
    else:
        xlabel = ''
    if 'ylabel' in kwargs:
        ylabel = kwargs['ylabel']
    else:
        ylabel = ''

    fig = plt.figure
    #vmin = np.nanmin(pvar)
    #vmax = np.nanmax(pvar)
    plt.scatter(px, py, c=pvar, cmap=cmap1, alpha=0.75, vmin=vmin1,
                vmax=vmax1, zorder=10)
    cbar1 = plt.colorbar(orientation='horizontal', shrink=0.25, pad=0.1)
    if method == 'pcolormesh':
        plt.pcolormesh(ma.masked_invalid(map_var), cmap=cmap2,
                       rasterized=True, alpha=0.5, vmin=vmin2, vmax=vmax2)
    elif method == 'contourf':
	    
        plt.contourf(ma.masked_invalid(map_var), cmap=cmap2, rasterized=True,
    	             alpha=0.5, vmin=vmin2, vmax=vmax2)
    plt.xlabel(xlabel)
    cbar2 = plt.colorbar(shrink=0.25);
    plt.axis([xmin, xmax, ymin, ymax])
    plt.ylabel(ylabel)
    return fig
##############################################################################
def plot_part(pvar, px, py, dtime=1, cmap='jet', marker='o', zorder=1,
              **kwargs):
    '''
    Plot particles trajectory along (px,py) with pvar in color 

    parameters : 
        pvar variable at particles points
        px, py: particles horizontal pozition
        dtime: integer, plot particles every dtime
        cmap : colormap for scatter plot
        marker : string, marker shape of scatter plot 
        zorder: float 
		        Set the zorder for the artist.
				Artists with lower zorder values are drawn first. 
	
	keyword arguments:
        size : scalar or array_like, shape (n, ), optional
		       The marker size in points**2.
		       Default is rcParams['lines.markersize'] ** 2. 
		coord = [ymin, ymax, xmin, max] : plot axes bounds
        vmin, vmax :floats, colorbar extent for pvar

    '''
    # Marker size
    if 'size' in kwargs:
        size = kwargs['size']
    else:
        size = plt.rcParams['lines.markersize'] ** 2
    # pvar color scale  
    if ('vmin' in kwargs) and ('vmax' in kwargs):
        vmin = kwargs['vmin']
        vmax = kwargs['vmax']
    else:
        vmin = np.nanmin(pvar)
        vmax = np.nanmax(pvar)

    plt.scatter(px[0:-1:dtime, :], py[0:-1:dtime, :], c=pvar[0:-1:dtime, :],
                cmap=cmap, alpha=0.75, vmin=vmin, vmax=vmax, zorder=10,
                marker=marker, s=size)
    return 
##############################################################################
def anim_scat(pvar, px, py, itime=0, cmap='jet', marker='o', alpha=0.75,
              zorder=1, **kwargs):
    '''
    Plot particles trajectory along (px,py) with pvar in color 

    parameters : 
        pvar variable at particles points
        px, py: particles horizontal pozition
        dtime: integer, plot particles every dtime
        cmap : colormap for scatter plot
        marker : string, marker shape of scatter plot 
        zorder: float 
                Set the zorder for the artist.
                Artists with lower zorder values are drawn first. 
    
    keyword arguments:
        size : scalar or array_like, shape (n, ), optional
               The marker size in points**2.
               Default is rcParams['lines.markersize'] ** 2. 
        coord = [ymin, ymax, xmin, max] : plot axes bounds
        vmin, vmax :floats, colorbar extent for pvar

    '''
    # Marker size
    if 'size' in kwargs:
        size = kwargs['size']
    else:
        size = plt.rcParams['lines.markersize'] ** 2
    # pvar color scale  
    if ('vmin' in kwargs):
         vmin = kwargs['vmin']
    else:
        vmin = np.nanmin(pvar)
    if ('vmax' in kwargs):
        vmax = kwargs['vmax']
    else:
        vmax = np.nanmax(pvar)
    if 'ax' in kwargs:
        ax = kwargs['ax']
        scat = ax.scatter(px[itime, :], py[itime, :], c=pvar[itime, :],
                cmap=cmap, alpha=aplha, vmin=vmin, vmax=vmax, zorder=zorder,
                marker=marker, s=size)
    else:
        scat = plt.scatter(px[itime, :], py[itime, :], c=pvar[itime, :],
                cmap=cmap, alpha=alpha, vmin=vmin, vmax=vmax, zorder=zorder,
                marker=marker, s=size)
    return

##############################################################################
def plot_seeding(px, py, marker='X', c='g', s=36, zorder=3, **kwargs):
    '''
    plot particles origin px, py at t=0
    
    parameters:
        marker: MarkerStyle optional
                 Marker shape
        
        c: color, sequence, or sequence of color, optional

        s: The marker size in points**2.
           Default is rcParams['lines.markersize'] ** 2.

        zorder: poltting order, high value means on top 

    '''
    if(len(px.shape) == 2):
        px = px[0, :]
        py = py[0, :]

    indx = np.isnan(px)
    if 'ax' in kwargs:
        ax = kwargs['ax']
        scat = ax.scatter(px[~indx], py[~indx], marker=marker, c=c, s=s,
                          zorder=zorder)
    else:
        scat = plt.scatter(px[~indx], py[~indx], marker=marker, c=c, s=s,
                           zorder=zorder)  

    return scat

##############################################################################
def plot_traj(px, py, tend=-1, c='grey', linewidth=None, zorder=2, alpha=0.5,
              **kwargs):
    '''
    Line plot of particles trajectory
    Helpfull along with high dtime value
    
    parameter: tend integer, last time-index for trajectory

    optional ax, axis 
            used to reference axis plot in animation
    '''
    if 'ax' in kwargs:
        ax = kwargs['ax']
        ax.plot(px[0:tend], py[0:tend], c=c, zorder=zorder, alpha=alpha)
    else:
        plt.plot(px[0:tend], py[0:tend], c=c, zorder=zorder, alpha=alpha)

    return

##############################################################################
def get_xy_lim(px, py, nx=2000, ny=2000, npts=5, **kwargs):
    '''
    returns integer maximum partciles displacement xmin, xmax, ymin, ymax
    
    parameters nx, ny should be domain's extent
               npts integer : extra points for a better plot

    '''
    xmin = np.max([0, np.int(np.floor(np.nanmin(px.reshape(-1)))) - npts])
    xmax = np.min([nx, np.int(np.ceil(np.nanmax(px.reshape(-1)))) + npts])
    ymin = np.max([0, np.int(np.floor(np.nanmin(py.reshape(-1)))) - npts])
    ymax = np.min([ny, np.int(np.ceil(np.nanmax(py.reshape(-1)))) + npts])
    return xmin, xmax, ymin, ymax
##############################################################################
def age(px, **kwargs):
    '''
    returns particles age in indexes 

    parameter px particle position
    
    optionnal itime integer, time index

    '''

    indx = ~np.isnan(px)
    if 'itime' in kwargs:
        age = np.sum(indx[:itime, :], axis=0)
    else:
        for itime in range(px.shape[0]):
            age[itime, :] = np.sum(indx[:itime, :], axis=0)
    return age

##############################################################################
def mask_contour(ax, topo, advdepth=0, **kwargs):
    '''
    contourf plot of land mask
    
    parameters: ax, plot axis
                topo, topography
    
    optionnal:  advdepth, advection depth (m < 0) in case of adv2D or zadvg
    '''

    ax = plt.contourf(topo.T, [0,-advdepth], 
                    cmap = col.LinearSegmentedColormap.from_list('my_colormap',
                                                         ['white','lightgray'],
                                                                   256))
    ax = plt.contourf(topo.T,[0,topo.min()], cmap =
                    col.LinearSegmentedColormap.from_list('my_colormap',
                                                          ['Gainsboro', 'gray']
                                                          , 256))
    ax = plt.contour(topo.T,[topo.min()],colors=('k',),linewidths=(0.5,));
    ax = plt.contour(topo.T,[-advdepth],colors=('k',),linewidths=(0.2,));
    return

    
    

