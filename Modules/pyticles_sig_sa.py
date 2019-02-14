###################################################################################
# Pyticles routines
###################################################################################
"""

!----------------------------------------------------------------------------------------------
! Python Routines for Pyticles
!
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
! 05/10/16:
!     - modify map_var to handle automatically u-,v- grids
! 16/01/26:
!     - Add periodic capabilities by adding ghost points (size=ng) at horizontal boundaries
!       [if ng>0 at an open boundary, we put nan so that particles going out are removed]
!     - The periodize... functions are used to add the ghost points for any type of variable
!----------------------------------------------------------------------------------------------

"""

###################################################################################
#Load modules
###################################################################################

#for numeric functions
import numpy as np
np.seterr(invalid='ignore') #remove anoying "Runtimewarning : invalid value encountered in absolute"

#for netcdf files
from netCDF4 import Dataset

#for timing
import time as tm

#pyticles routines
import pyticles_3d_sig_sa as partF

#copy data
from copy import copy
###################################################################################






###################################################################################
#@profile
def subsection(px,py,dx=1,dy=1,maxvel=[1,1],delt=1,nx=2000,ny=2000,ng=0,nadv=0,**kwargs):
    """ Finds index subrange to move particles around"""

    if 'offset' in kwargs:
       offset_x = kwargs['offset']
       offset_y = kwargs['offset']
    else:
       offset_x = np.abs(maxvel[0]*delt/dx)
       offset_y = np.abs(maxvel[1]*delt/dy)
       
    i0 = int(np.floor(np.nanmin(px) + ng  - offset_x) - nadv) -1
    i0 = max(0,i0)
    i1 = int(np.ceil(np.nanmax(px) + ng  + offset_x) + nadv) + 2
    i1 = min(nx,i1)
    j0 = int(np.floor(np.nanmin(py) + ng  - offset_y) - nadv) - 1
    j0 = max(0,j0)
    j1 = int(np.ceil(np.nanmax(py)+ ng  + offset_y) + nadv)  + 2
    j1 = min(ny,j1)

    #########################
    #if subdomain falls inside ghost points, make it include all of the ghost points (ng)
    if i0<ng: i0=0
    if i1>nx-ng: i1=nx
    if j0<ng: j0=0
    if j1>ny-ng: j1=ny
    #########################

    return [j0,j1,i0,i1]

###################################################################################

   
#@profile   
def cull(px,py,pz,nx,ny,nz,x_periodic=False,y_periodic=False,ng=0):
    """ Set particle positions that are out of range to nan"""

    if x_periodic:
        px[px<-1] = px[px<-1] + nx
        px[px>nx-1] = px[px>nx-1] - nx
    else:
        px[px<=1-ng] = np.nan; px[px>=nx-3+ng] = np.nan; 

    if y_periodic:
        py[py<-1] = py[py<-1] + ny
        py[py>ny-1] = py[py>ny-1] - ny
    else:
        py[py<=1-ng] = np.nan; py[py>=ny-3+ng] = np.nan; 
        #pz[pz<0] = np.nan; pz[pz>nz] = np.nan; 


    py[np.isnan(px)] = np.nan
    px[np.isnan(py)] = np.nan
    pz[np.isnan(py)] = np.nan

    return [px,py,pz]
   
    
#@profile   
def cull2d(px,py,nx,ny,x_periodic=False,y_periodic=False,ng=0):
    """ Set particle positions that are out of range to nan"""

    if x_periodic:
        px[px<-1] = px[px<-1] + nx
        px[px>nx-1] = px[px>nx-1] - nx
    else:
        px[px<=1-ng] = np.nan; px[px>=nx-3+ng] = np.nan; 

    if y_periodic:
        py[py<-1] = py[py<-1] + ny
        py[py>ny-1] = py[py>ny-1] - ny
    else:
        py[py<=1-ng] = np.nan; py[py>=ny-3+ng] = np.nan; 
        #pz[pz<0] = np.nan; pz[pz>nz] = np.nan; 


    py[np.isnan(px)] = np.nan
    px[np.isnan(py)] = np.nan

    return [px,py]  


###################################################################################

   
#@profile   
def kick(pz,nz):
    """ Give a kick to particles trapped at surface/bottom  """

    eps = 0.01
    #eps = 0.5
    
    pz[pz<=eps] = eps; pz[pz>=nz-eps] = nz-eps; 


    return [pz]
   
###################################################################################
  
#@profile   
def get_vel_io(simul,pm=None,pn=None,timing=False,x_periodic=False,y_periodic=False,ng=0,**kwargs):  


    if 'coord' in  kwargs:
        coord = kwargs['coord'][0:4]
    else: 
        coord = simul.coord[0:4]
    
    [ny1tot,ny2tot,nx1tot,nx2tot] = simul.coord[0:4]
    
    if timing: tstart2 = tm.time()     
  
    nc = Dataset(simul.ncfile, 'r')
    [ny1,ny2,nx1,nx2] = coord
    nz = len(simul.coord[4])

    ################################
    mask = copy(simul.mask)
    mask[np.isnan(mask)]=0
    ################################

    u = np.zeros((nx2-nx1-1,ny2-ny1,nz))*np.nan
    v = np.zeros((nx2-nx1,ny2-ny1-1,nz))*np.nan

    ################################

    nw = min(ng,nx1); ne = min(ng,nx2tot-nx1tot+2*ng-nx2)
    ns = min(ng,ny1); nn = min(ng,ny2tot-ny1tot+2*ng-ny2)


    # Fill inside points [if x periodic shift one index right in netcdf file]
    if x_periodic: iper=1
    else: iper = 0
    u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1-ns:ny2-2*ng+nn,nx1+iper-nw:nx2-1+iper-2*ng+ne]))
    
    u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn,:] = (u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn,:].T * (mask[nx1+1-nw:nx2-2*ng+ne,ny1-ns:ny2-2*ng+nn]*mask[nx1-nw:nx2-1-2*ng+ne,ny1-ns:ny2-2*ng+nn]).T).T

    # Fill inside points [if y periodic shift one index north in netcdf file]
    if y_periodic: jper=1
    else: jper = 0
    v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny1-ns+jper:ny2-1+jper-2*ng+nn,nx1-nw:nx2-2*ng+ne]))

    v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn,:] = (v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn,:].T * (mask[nx1-nw:nx2-2*ng+ne,ny1+1-ns:ny2-2*ng+nn]*mask[nx1-nw:nx2-2*ng+ne,ny1-ns:ny2-1-2*ng+nn]).T).T

    ################################
    # Filling Ghost points
    ################################

    if nw<ng and x_periodic:
        u[ng-nw-1,ng-ns:ny2-ny1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1-ns:ny2-2*ng+nn,nx1tot]))
        for i in range(1,ng):
            u[ng-nw-1-i,ng-ns:ny2-ny1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1-ns:ny2-2*ng+nn,nx2tot-i]))
        for i in range(ng):
            v[ng-nw-1-i,ng-ns:ny2-ny1-1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny1-ns+jper:ny2-1+jper-2*ng+nn,nx2tot-1-i]))
        nw=ng 

    if ne<ng and x_periodic:
        for i in range(ng):
            u[nx2-nx1-1-ng+ne+i,ng-ns:ny2-ny1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1-ns:ny2-2*ng+nn,nx1tot+i]))
        for i in range(ng):
            v[nx2-nx1-ng+ne+i,ng-ns:ny2-ny1-1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny1-ns+jper:ny2-1+jper-2*ng+nn,nx1tot+i]))
        ne=ng

    if ns<ng and y_periodic:
        v[ng-nw:nx2-nx1-ng+ne,ng-ns-1,:] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny1tot,nx1-nw:nx2-2*ng+ne]))
        for i in range(1,ng):
            v[ng-nw:nx2-nx1-ng+ne,ng-ns-1-i,:] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny2tot-i,nx1-nw:nx2-2*ng+ne]))
        for i in range(1,ng):
            u[ng-nw:nx2-nx1-1-ng+ne,ng-ns-1-i,:] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny2tot-1-i,nx1+iper-nw:nx2-1+iper-2*ng+ne]))

    if nn<ng and y_periodic:
        for i in range(ng):
            v[ng-nw:nx2-nx1-ng+ne,ny2-ny1-1-ng+nn+i,:] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny1tot+i,nx1-nw:nx2-2*ng+ne]))
        for i in range(1,ng):
            u[ng-nw:nx2-nx1-1-ng+ne,ny2-ny1-ng+nn+i,:] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1tot+i,nx1+iper-nw:nx2-1+iper-2*ng+ne]))


    ################################

    if timing: print 'get u,v from file....', tm.time()-tstart2
    if timing: tstart2 = tm.time()    

    ################################

    try:
        w = periodize3d_fromnc(simul,'omega',coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)
    except:
        print 'no omega in file, computing'
        [z_r,z_w] = get_depths(simul,coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)
        w = partF.get_omega(u,v,z_r,z_w,pm,pn)
        w[np.isnan(w)] = 0.
        if x_periodic and nx1<ng and nx2>nx2tot-nx1tot+ng: 
            w[0,:,:] = w[nx2tot,:,:]
            w[-1,:,:] = w[nx1tot+2*ng-1,:,:]
        if y_periodic and ny1<ng and ny2>ny2tot-ny1tot+ng: 
            w[:,0,:] = w[:,ny2tot,:]
            w[:,-1,:] = w[:,ny1tot+2*ng-1,:]

    nc.close()
    
    if timing: print 'get w from file....', tm.time()-tstart2
    if timing: tstart2 = tm.time()    

    return u,v,w

###################################################################################
'''
#@profile   
def get_vel_io_2d(simul,pm=None,pn=None,timing=False,x_periodic=False,y_periodic=False,ng=0, **kwargs):  


    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
    
    [ny1tot,ny2tot,nx1tot,nx2tot] = simul.coord[0:4]
    
    if timing: tstart2 = tm.time()     
  
    nc = Dataset(simul.ncfile, 'r')
    [ny1,ny2,nx1,nx2] = coord

    ################################
    mask = copy(simul.mask)
    mask[np.isnan(mask)]=0
    ################################

    u = np.zeros((nx2-nx1-1,ny2-ny1))*np.nan
    v = np.zeros((nx2-nx1,ny2-ny1-1))*np.nan

    ################################

    nw = min(ng,nx1); ne = min(ng,nx2tot-nx1tot+2*ng-nx2)
    ns = min(ng,ny1); nn = min(ng,ny2tot-ny1tot+2*ng-ny2)

    # Fill inside points [if x periodic shift one index right in netcdf file]
    if x_periodic: iper=1
    else: iper = 0
    u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,-1,ny1-ns:ny2-2*ng+nn,nx1+iper-nw:nx2-1+iper-2*ng+ne]))
    
    u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn] = (u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn].T * (mask[nx1+1-nw:nx2-2*ng+ne,ny1-ns:ny2-2*ng+nn]*mask[nx1-nw:nx2-1-2*ng+ne,ny1-ns:ny2-2*ng+nn]).T).T

    # Fill inside points [if y periodic shift one index north in netcdf file]
    if y_periodic: jper=1
    else: jper = 0
    v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,-1,ny1-ns+jper:ny2-1+jper-2*ng+nn,nx1-nw:nx2-2*ng+ne]))

    v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn] = (v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn].T * (mask[nx1-nw:nx2-2*ng+ne,ny1+1-ns:ny2-2*ng+nn]*mask[nx1-nw:nx2-2*ng+ne,ny1-ns:ny2-1-2*ng+nn]).T).T

    ################################
    # Filling Ghost points
    ################################

    if nw<ng and x_periodic:
        u[ng-nw-1,ng-ns:ny2-ny1-ng+nn] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,-1,ny1-ns:ny2-2*ng+nn,nx1tot]))
        for i in range(1,ng):
            u[ng-nw-1-i,ng-ns:ny2-ny1-ng+nn] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,-1,ny1-ns:ny2-2*ng+nn,nx2tot-i]))
        for i in range(ng):
            v[ng-nw-1-i,ng-ns:ny2-ny1-1-ng+nn] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,-1,ny1-ns+jper:ny2-1+jper-2*ng+nn,nx2tot-1-i]))
        nw=ng 

    if ne<ng and x_periodic:
        for i in range(ng):
            u[nx2-nx1-1-ng+ne+i,ng-ns:ny2-ny1-ng+nn] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,-1,ny1-ns:ny2-2*ng+nn,nx1tot+i]))
        for i in range(ng):
            v[nx2-nx1-ng+ne+i,ng-ns:ny2-ny1-1-ng+nn] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,-1,ny1-ns+jper:ny2-1+jper-2*ng+nn,nx1tot+i]))
        ne=ng

    if ns<ng and y_periodic:
        v[ng-nw:nx2-nx1-ng+ne,ng-ns-1] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,-1,ny1tot,nx1-nw:nx2-2*ng+ne]))
        for i in range(1,ng):
            v[ng-nw:nx2-nx1-ng+ne,ng-ns-1-i] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,-1,ny2tot-i,nx1-nw:nx2-2*ng+ne]))
        for i in range(1,ng):
            u[ng-nw:nx2-nx1-1-ng+ne,ng-ns-1-i] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,-1,ny2tot-1-i,nx1+iper-nw:nx2-1+iper-2*ng+ne]))

    if nn<ng and y_periodic:
        for i in range(ng):
            v[ng-nw:nx2-nx1-ng+ne,ny2-ny1-1-ng+nn+i] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,-1,ny1tot+i,nx1-nw:nx2-2*ng+ne]))
        for i in range(1,ng):
            u[ng-nw:nx2-nx1-1-ng+ne,ny2-ny1-ng+nn+i] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,-1,ny1tot+i,nx1+iper-nw:nx2-1+iper-2*ng+ne]))


    ################################

    if timing: print 'get u,v from file....', tm.time()-tstart2
    if timing: tstart2 = tm.time()    

    ################################

    return u,v
'''
###################################################################################


###################################################################################

#@profile   
def get_vel_io_2d(simul,pm=None,pn=None,timing=False,x_periodic=False,y_periodic=False,ng=0, advdepth = -1,**kwargs):  

    '''
    get horizontal velocities at depth advdepth (default = surface)

    '''
    if 'coord' in  kwargs:
        coord = kwargs['coord'][0:4]
    else: 
        coord = simul.coord[0:4]
    
    [ny1tot,ny2tot,nx1tot,nx2tot] = simul.coord[0:4]
    
    if timing: tstart2 = tm.time()     
  
    nc = Dataset(simul.ncfile, 'r')
    [ny1,ny2,nx1,nx2] = coord

    ################################

    mask = copy(simul.mask)
    mask[np.isnan(mask)]=0

    ################################

    u = np.zeros((nx2-nx1-1,ny2-ny1))*np.nan
    v = np.zeros((nx2-nx1,ny2-ny1-1))*np.nan

    ################################

    def interp(simul,varname,coord=coord,advdepth=advdepth):

        coordz = copy(coord)

        if varname=='u': imin=1;jmin=0; coordz[3]+=1
        elif varname=='v': imin=0;jmin=1; coordz[1]+=1

        if advdepth==0:
            # use surface field
            varz = simul.Forder(np.squeeze(nc.variables[varname][simul.infiletime,-1,coord[0]:coord[1],coord[2]:coord[3]]))
        else:
            [z_r,z_w] = get_depths(simul,coord=coordz)
            var0  = simul.Forder(np.squeeze(nc.variables[varname][simul.infiletime,:,coord[0]:coord[1],coord[2]:coord[3]]))
            varz = vinterp(var0, [advdepth], z_r, z_w, imin=imin,jmin=jmin)[:,:,0]

        return varz

    ################################

    nw = min(ng,nx1); ne = min(ng,nx2tot-nx1tot+2*ng-nx2)
    ns = min(ng,ny1); nn = min(ng,ny2tot-ny1tot+2*ng-ny2)

    # Fill inside points [if x periodic shift one index right in netcdf file]
    if x_periodic: iper=1
    else: iper = 0
    u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn] = interp(simul,'u',coord=[ny1-ns,ny2-2*ng+nn,nx1+iper-nw,nx2-1+iper-2*ng+ne])
    
    u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn] = (u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn].T * (mask[nx1+1-nw:nx2-2*ng+ne,ny1-ns:ny2-2*ng+nn]*mask[nx1-nw:nx2-1-2*ng+ne,ny1-ns:ny2-2*ng+nn]).T).T

    # Fill inside points [if y periodic shift one index north in netcdf file]
    if y_periodic: jper=1
    else: jper = 0
    v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn] = interp(simul,'v',coord=[ny1-ns+jper,ny2-1+jper-2*ng+nn,nx1-nw,nx2-2*ng+ne])

    v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn] = (v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn].T * (mask[nx1-nw:nx2-2*ng+ne,ny1+1-ns:ny2-2*ng+nn]*mask[nx1-nw:nx2-2*ng+ne,ny1-ns:ny2-1-2*ng+nn]).T).T

    ################################
    # Filling Ghost points
    ################################

    if nw<ng and x_periodic:
        u[ng-nw-1,ng-ns:ny2-ny1-ng+nn] = interp(simul,'u',coord=[ny1-ns,ny2-2*ng+nn,nx1tot,nx1tot+1])
        for i in range(1,ng):
            u[ng-nw-1-i,ng-ns:ny2-ny1-ng+nn] = interp(simul,'u',coord=[ny1-ns,ny2-2*ng+nn,nx2tot-i,nx2tot-i+1])
        for i in range(ng):
            v[ng-nw-1-i,ng-ns:ny2-ny1-1-ng+nn] = interp(simul,'v',coord=[ny1-ns+jper,ny2-1+jper-2*ng+nn,nx2tot-1-i,nx2tot-1-i+1])
        nw=ng 

    if ne<ng and x_periodic:
        for i in range(ng):
            u[nx2-nx1-1-ng+ne+i,ng-ns:ny2-ny1-ng+nn] = interp(simul,'u',coord=[ny1-ns,ny2-2*ng+nn,nx1tot+i,nx1tot+i+1])
        for i in range(ng):
            v[nx2-nx1-ng+ne+i,ng-ns:ny2-ny1-1-ng+nn] = interp(simul,'v',coord=[ny1-ns+jper,ny2-1+jper-2*ng+nn,nx1tot+i,nx1tot+i+1])
        ne=ng

    if ns<ng and y_periodic:
        v[ng-nw:nx2-nx1-ng+ne,ng-ns-1] = interp(simul,'v',coord=[ny1tot,ny1tot+1,nx1-nw,nx2-2*ng+ne])
        for i in range(1,ng):
            v[ng-nw:nx2-nx1-ng+ne,ng-ns-1-i] = interp(simul,'v',coord=[ny2tot-i,ny2tot-i+1,nx1-nw,nx2-2*ng+ne])
        for i in range(1,ng):
            u[ng-nw:nx2-nx1-1-ng+ne,ng-ns-1-i] = interp(simul,'u',coord=[ny2tot-1-i,ny2tot-1-i+1,nx1+iper-nw,nx2-1+iper-2*ng+ne])

    if nn<ng and y_periodic:
        for i in range(ng):
            v[ng-nw:nx2-nx1-ng+ne,ny2-ny1-1-ng+nn+i] = interp(simul,'v',coord=[ny1tot+i,ny1tot+i+1,nx1-nw,nx2-2*ng+ne])
        for i in range(1,ng):
            u[ng-nw:nx2-nx1-1-ng+ne,ny2-ny1-ng+nn+i] = interp(simul,'u',coord=[ny1tot+i,ny1tot+i+1,nx1+iper-nw,nx2-1+iper-2*ng+ne])


    ################################

    if timing: print 'get u,v from file....', tm.time()-tstart2
    if timing: tstart2 = tm.time()    

    ################################


    return u,v

###################################################################################

#@profile   
def get_vel_io_surf(simul,pm=None,pn=None,timing=False,x_periodic=False,y_periodic=False,ng=0,**kwargs):  

    '''
    get surface horizontal velocities

    '''
    if 'coord' in  kwargs:
        coord = kwargs['coord'][0:4]
    else: 
        coord = simul.coord[0:4]
    
    [ny1tot,ny2tot,nx1tot,nx2tot] = simul.coord[0:4]
    
    if timing: tstart2 = tm.time()     
  
    nc = Dataset(simul.ncfile, 'r')
    [ny1,ny2,nx1,nx2] = coord

    ################################
    mask = copy(simul.mask)
    mask[np.isnan(mask)]=0
    ################################

    u = np.zeros((nx2-nx1-1,ny2-ny1))*np.nan
    v = np.zeros((nx2-nx1,ny2-ny1-1))*np.nan

    ################################

    nw = min(ng,nx1); ne = min(ng,nx2tot-nx1tot+2*ng-nx2)
    ns = min(ng,ny1); nn = min(ng,ny2tot-ny1tot+2*ng-ny2)

    # Fill inside points [if x periodic shift one index right in netcdf file]
    if x_periodic: iper=1
    else: iper = 0
    u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn] = simul.Forder(np.squeeze(nc.variables['surf_u'][simul.infiletime,ny1-ns:ny2-2*ng+nn,nx1+iper-nw:nx2-1+iper-2*ng+ne]))
    
    u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn] = (u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn].T * (mask[nx1+1-nw:nx2-2*ng+ne,ny1-ns:ny2-2*ng+nn]*mask[nx1-nw:nx2-1-2*ng+ne,ny1-ns:ny2-2*ng+nn]).T).T

    # Fill inside points [if y periodic shift one index north in netcdf file]
    if y_periodic: jper=1
    else: jper = 0
    v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn] = simul.Forder(np.squeeze(nc.variables['surf_v'][simul.infiletime,ny1-ns+jper:ny2-1+jper-2*ng+nn,nx1-nw:nx2-2*ng+ne]))

    v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn] = (v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn].T * (mask[nx1-nw:nx2-2*ng+ne,ny1+1-ns:ny2-2*ng+nn]*mask[nx1-nw:nx2-2*ng+ne,ny1-ns:ny2-1-2*ng+nn]).T).T

    ################################
    # Filling Ghost points
    ################################

    if nw<ng and x_periodic:
        u[ng-nw-1,ng-ns:ny2-ny1-ng+nn] = simul.Forder(np.squeeze(nc.variables['surf_u'][simul.infiletime,ny1-ns:ny2-2*ng+nn,nx1tot]))
        for i in range(1,ng):
            u[ng-nw-1-i,ng-ns:ny2-ny1-ng+nn] = simul.Forder(np.squeeze(nc.variables['surf_u'][simul.infiletime,ny1-ns:ny2-2*ng+nn,nx2tot-i]))
        for i in range(ng):
            v[ng-nw-1-i,ng-ns:ny2-ny1-1-ng+nn] = simul.Forder(np.squeeze(nc.variables['surf_v'][simul.infiletime,ny1-ns+jper:ny2-1+jper-2*ng+nn,nx2tot-1-i]))
        nw=ng 

    if ne<ng and x_periodic:
        for i in range(ng):
            u[nx2-nx1-1-ng+ne+i,ng-ns:ny2-ny1-ng+nn] = simul.Forder(np.squeeze(nc.variables['surf_u'][simul.infiletime,ny1-ns:ny2-2*ng+nn,nx1tot+i]))
        for i in range(ng):
            v[nx2-nx1-ng+ne+i,ng-ns:ny2-ny1-1-ng+nn] = simul.Forder(np.squeeze(nc.variables['surf_v'][simul.infiletime,ny1-ns+jper:ny2-1+jper-2*ng+nn,nx1tot+i]))
        ne=ng

    if ns<ng and y_periodic:
        v[ng-nw:nx2-nx1-ng+ne,ng-ns-1] = simul.Forder(np.squeeze(nc.variables['surf_v'][simul.infiletime,ny1tot,nx1-nw:nx2-2*ng+ne]))
        for i in range(1,ng):
            v[ng-nw:nx2-nx1-ng+ne,ng-ns-1-i] = simul.Forder(np.squeeze(nc.variables['surf_v'][simul.infiletime,ny2tot-i,nx1-nw:nx2-2*ng+ne]))
        for i in range(1,ng):
            u[ng-nw:nx2-nx1-1-ng+ne,ng-ns-1-i] = simul.Forder(np.squeeze(nc.variables['surf_u'][simul.infiletime,ny2tot-1-i,nx1+iper-nw:nx2-1+iper-2*ng+ne]))

    if nn<ng and y_periodic:
        for i in range(ng):
            v[ng-nw:nx2-nx1-ng+ne,ny2-ny1-1-ng+nn+i] = simul.Forder(np.squeeze(nc.variables['surf_v'][simul.infiletime,ny1tot+i,nx1-nw:nx2-2*ng+ne]))
        for i in range(1,ng):
            u[ng-nw:nx2-nx1-1-ng+ne,ny2-ny1-ng+nn+i] = simul.Forder(np.squeeze(nc.variables['surf_u'][simul.infiletime,ny1tot+i,nx1+iper-nw:nx2-1+iper-2*ng+ne]))


    ################################

    if timing: print 'get u,v from file....', tm.time()-tstart2
    if timing: tstart2 = tm.time()    

    ################################

    return u,v


###################################################################################


def periodize3d_fromnc(simul,variable,coord,x_periodic=False,y_periodic=False,ng=0):

    [ny1tot,ny2tot,nx1tot,nx2tot] = simul.coord[0:4]
    [ny1,ny2,nx1,nx2] = coord[0:4]
    nz = len(simul.coord[4])
    nc = Dataset(simul.ncfile, 'r')

    ################################
    mask = copy(simul.mask)
    mask[np.isnan(mask)]=0
    ################################

    nw = min(ng,nx1); ne = min(ng,nx2tot-nx1tot+2*ng-nx2)
    ns = min(ng,ny1); nn = min(ng,ny2tot-ny1tot+2*ng-ny2)
    
    if variable=='omega':
        myvar = np.zeros((nx2-nx1,ny2-ny1,nz+1))*np.nan
    else:
        myvar = np.zeros((nx2-nx1,ny2-ny1,nz))*np.nan

    myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables[variable][simul.infiletime,:,ny1-ns:ny2-2*ng+nn,nx1-nw:nx2-2*ng+ne]))

    myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn,:] = (myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn,:].T * (mask[nx1-nw:nx2-2*ng+ne,ny1-ns:ny2-2*ng+nn]).T).T

    if nw<ng and x_periodic:
        for i in range(ng):
            myvar[ng-nw-1-i,ng-ns:ny2-ny1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables[variable][simul.infiletime,:,ny1-ns:ny2-2*ng+nn,nx2tot-1-i]))
        nw=ng

    if ne<ng and x_periodic:
        for i in range(ng):
            myvar[nx2-nx1-ng+ne+i,ng-ns:ny2-ny1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables[variable][simul.infiletime,:,ny1-ns:ny2-2*ng+nn,nx1tot+i,]))
        ne=ng

    if ns<ng and y_periodic:
        for i in range(1,ng):
            myvar[ng-nw:nx2-nx1-ng+ne,ng-ns-1-i,:] = simul.Forder(np.squeeze(nc.variables[variable][simul.infiletime,:,ny2tot-1-i,nx1-nw:nx2-2*ng+ne,]))

    if nn<ng and y_periodic:
        for i in range(1,ng):
            myvar[ng-nw:nx2-nx1-ng+ne,ny2-ny1-ng+nn+i,:] = simul.Forder(np.squeeze(nc.variables[variable][simul.infiletime,:,ny1tot+i,nx1-nw:nx2-2*ng+ne,]))

    return myvar

###################################################################################


def periodize3d_fromvar(simul,var3d,coord,x_periodic=False,y_periodic=False,ng=0):

    [ny1tot,ny2tot,nx1tot,nx2tot] = simul.coord[0:4]
    [ny1,ny2,nx1,nx2] = coord[0:4]
    nz = var3d.shape[2]
    nc = Dataset(simul.ncfile, 'r')

    ################################
    mask = copy(simul.mask)
    mask[np.isnan(mask)]=0
    ################################

    nw = min(ng,nx1); ne = min(ng,nx2tot-nx1tot+2*ng-nx2)
    ns = min(ng,ny1); nn = min(ng,ny2tot-ny1tot+2*ng-ny2)
    
    myvar = np.zeros((nx2-nx1,ny2-ny1,nz))*np.nan

    myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn,:] = var3d[nx1-nw:nx2-2*ng+ne,ny1-ns:ny2-2*ng+nn,:]

    myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn,:] = (myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn,:].T * (mask[nx1-nw:nx2-2*ng+ne,ny1-ns:ny2-2*ng+nn]).T).T

    if nw<ng and x_periodic:
        for i in range(ng):
            myvar[ng-nw-1-i,ng-ns:ny2-ny1-ng+nn,:] = var3d[nx2tot-1-i,ny1-ns:ny2-2*ng+nn,:]
        nw=ng

    if ne<ng and x_periodic:
        for i in range(ng):
            myvar[nx2-nx1-ng+ne+i,ng-ns:ny2-ny1-ng+nn,:] = var3d[nx1tot+i,ny1-ns:ny2-2*ng+nn,:]
        ne=ng

    if ns<ng and y_periodic:
        for i in range(1,ng):
            myvar[ng-nw:nx2-nx1-ng+ne,ng-ns-1-i,:] = var3d[nx1-nw:nx2-2*ng+ne,ny2tot-1-i,:]

    if nn<ng and y_periodic:
        for i in range(1,ng):
            myvar[ng-nw:nx2-nx1-ng+ne,ny2-ny1-ng+nn+i,:] = var3d[nx1-nw:nx2-2*ng+ne,ny1tot+i,:]

    return myvar

###################################################################################


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

###################################################################################


def periodize2d_fromnc_depth(simul,variable,coord,x_periodic=False,y_periodic=False,ng=0,advdepth=-1):
    '''
    A tester
    '''
    [ny1tot,ny2tot,nx1tot,nx2tot] = simul.coord[0:4]
    [ny1,ny2,nx1,nx2] = coord[0:4]

    nc = Dataset(simul.ncfile, 'r')

    ################################
    mask = copy(simul.mask)
    mask[np.isnan(mask)]=0
    ################################

    def interp(simul,varname,coord=coord,advdepth=advdepth):

        coordz = copy(coord)

        if varname=='u': imin=1;jmin=0; coordz[3]+=1
        elif varname=='v': imin=0;jmin=1; coordz[1]+=1
        else: imin=0; jmin=0; 

        if advdepth==0:
            # use surface field
            varz = simul.Forder(np.squeeze(nc.variables[varname][simul.infiletime,-1,coord[0]:coord[1],coord[2]:coord[3]]))
        else:
            [z_r,z_w] = get_depths(simul,coord=coordz)
            var0  = simul.Forder(np.squeeze(nc.variables[varname][simul.infiletime,:,coord[0]:coord[1],coord[2]:coord[3]]))
            varz = vinterp(var0, [advdepth], z_r, z_w, imin=imin,jmin=jmin)[:,:,0]

        return varz

    ################################

    nw = min(ng,nx1); ne = min(ng,nx2tot-nx1tot+2*ng-nx2)
    ns = min(ng,ny1); nn = min(ng,ny2tot-ny1tot+2*ng-ny2)

    myvar = np.zeros((nx2-nx1,ny2-ny1))*np.nan

    myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn] = interp(simul,variable,coord=[ny1-ns,ny2-2*ng+nn,nx1-nw,nx2-2*ng+ne])

    myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn] = (myvar[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-ng+nn].T * (mask[nx1-nw:nx2-2*ng+ne,ny1-ns:ny2-2*ng+nn]).T).T

    if nw<ng and x_periodic:
        for i in range(ng):
            myvar[ng-nw-1-i,ng-ns:ny2-ny1-ng+nn] = interp(simul,variable,coord=[ny1-ns,ny2-2*ng+nn,nx2tot-1-i,nx2tot-1-i+1])
        nw=ng

    if ne<ng and x_periodic:
        for i in range(ng):
            myvar[nx2-nx1-ng+ne+i,ng-ns:ny2-ny1-ng+nn] =  interp(simul,variable,coord=[ny1-ns,ny2-2*ng+nn,nx1tot+i,nx1tot+i+1])
        ne=ng

    if ns<ng and y_periodic:
        for i in range(1,ng):
            myvar[ng-nw:nx2-nx1-ng+ne,ng-ns-1-i] =  interp(simul,variable,coord=[ny2tot-1-i,ny2tot-1-i+1,nx1-nw,nx2-2*ng+ne])

    if nn<ng and y_periodic:
        for i in range(1,ng):
            myvar[ng-nw:nx2-nx1-ng+ne,ny2-ny1-ng+nn+i] =  interp(simul,variable,coord=[ny1tot+i,ny1tot+i+1,nx1-nw,nx2-2*ng+ne])

    return myvar

###################################################################################

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

###################################################################################

    
#@profile
def get_ts_io(simul,x_periodic=False,y_periodic=False,ng=0,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]

    nc = Dataset(simul.ncfile, 'r')
    [ny1,ny2,nx1,nx2] = coord; 

    temp = periodize3d_fromnc(simul,'temp',coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)

    try:
        salt = periodize3d_fromnc(simul,'salt',coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)
    except:
        salt = np.zeros(temp.shape)

    nc.close()

    return temp,salt
###################################################################################


#@profile
def get_ts_io_surf(simul,x_periodic=False,y_periodic=False,ng=0,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]

    nc = Dataset(simul.ncfile, 'r')
    [ny1,ny2,nx1,nx2] = coord; 

    temp = periodize2d_fromnc(simul,'surf_t',coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)

    try:
        salt = periodize2d_fromnc(simul,'surf_s',coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)
    except:
        salt = np.zeros(temp.shape)

    nc.close()

    return temp,salt
###################################################################################


#@profile
def get_ts_io_2d(simul,x_periodic=False,y_periodic=False,ng=0,advdepth=-1,**kwargs):

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else:
        coord = simul.coord[0:4]

    nc = Dataset(simul.ncfile, 'r')
    [ny1,ny2,nx1,nx2] = coord;

    ################################

    temp = periodize2d_fromnc_depth(simul, 'temp', coord, x_periodic=x_periodic, y_periodic=y_periodic, ng=ng, advdepth=advdepth)

    try:
        salt = periodize2d_fromnc_depth(simul, 'salt', coord, x_periodic=x_periodic, y_periodic=y_periodic, ng=ng, advdepth=advdepth)
    except:
        salt = np.zeros(temp.shape)

    nc.close()

    return temp,salt



###################################################################################
   
    
#@profile
def get_t_io(simul,x_periodic=False,y_periodic=False,ng=0,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]

    nc = Dataset(simul.ncfile, 'r')
    [ny1,ny2,nx1,nx2] = coord; 

    temp = periodize3d_fromnc(simul,'temp',coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)

    nc.close()

    return temp

###################################################################################
# Compute T,S at each particle position
###################################################################################
#@profile


def map_ts(simul,temp,salt,px,py,pz,ng=0,**kwargs):  
    '''
    T,S on horizontal and vertical rho-grid
    '''

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    [j0,j1,i0,i1]=coord; k0=0; nqmx=px.shape[0]
    
    [ptemp,psalt] = partF.interp_3d_ts(px,py,pz,temp,salt,ng,nqmx,i0,j0,k0)

    return [ptemp,psalt]


###################################################################################
# Compute a 3D variable at each particle position
###################################################################################
#@profile


def map_var(simul,var0,px,py,pz,ng=0,**kwargs):  
    '''
    var on horizontal and vertical rho-grid
    '''

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    [j0,j1,i0,i1]=coord; k0=0; nqmx=px.shape[0]
   
    print var0.shape
    print [i1-i0,j1-j0,len(simul.coordmax[4])]
   
    if var0.shape==(i1-i0,j1-j0,len(simul.coordmax[4])):
        pvar0 = partF.interp_3d(px,py,pz,var0,ng,nqmx,i0,j0,k0)
    elif var0.shape==(i1-i0-1,j1-j0-1,len(simul.coordmax[4])):
        pvar0 = partF.interp_3d_psi(px,py,pz,var0,ng,nqmx,i0,j0,k0)
    elif var0.shape==(i1-i0-1,j1-j0,len(simul.coordmax[4])):
        pvar0 = partF.interp_3d_u(px,py,pz,var0,ng,nqmx,i0,j0,k0)
    elif var0.shape==(i1-i0,j1-j0-1,len(simul.coordmax[4])):
        pvar0 = partF.interp_3d_v(px,py,pz,var0,ng,nqmx,i0,j0,k0)
    elif var0.shape==(i1-i0,j1-j0,len(simul.coordmax[4])+1):
        pvar0 = partF.interp_3d_w(px,py,pz,var0,ng,nqmx,i0,j0,k0)
    elif var0.shape==(i1-i0-1,j1-j0-1,len(simul.coordmax[4])+1):
        pvar0 = partF.interp_3d_psiw(px,py,pz,var0,ng,nqmx,i0,j0,k0)
        

    return pvar0



'''
###################################################################################
# Compute a 3D variable at each particle position
###################################################################################
#@profile


def map_varw(simul,var0,px,py,pz,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    [j0,j1,i0,i1]=coord; k0=0; nqmx=px.shape[0]
    
    pvar0 = partF.interp_3dw(px,py,pz,var0,nqmx,i0,j0,k0)


    return pvar0
'''



###################################################################################
# Compute topo at each particle position
###################################################################################
#@profile


def map_topo(simul,px,py,ng=0,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    [j0,j1,i0,i1]=coord;
    
    ptopo = partF.interp_2d(px,py,simul.topo,ng,px.shape[0],i0,j0)

    return ptopo
    
    
    
###################################################################################
# Compute var2d at each particle position
###################################################################################
#@profile


def map_var2d(simul,var2d,px,py,ng=0,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    [j0,j1,i0,i1]=coord;

    if var2d.shape==(i1-i0,j1-j0):
        pvar = partF.interp_2d(px,py,var2d,ng,px.shape[0],i0,j0)
    elif var2d.shape==(i1-i0-1,j1-j0-1):
        pvar = partF.interp_2d_psi(px,py,var2d,ng,px.shape[0],i0,j0)
    elif var2d.shape==(i1-i0-1,j1-j0):
        pvar = partF.interp_2d_u(px,py,var2d,ng,px.shape[0],i0,j0)
    elif var2d.shape==(i1-i0,j1-j0-1):
        pvar = partF.interp_2d_v(px,py,var2d,ng,px.shape[0],i0,j0)

        

    #

    return pvar
    
    
###################################################################################
# Compute lon,lat at each particle position
###################################################################################
#@profile

def map_lonlat(simul,px,py,ng=0,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    [j0,j1,i0,i1]=coord;
    
    plon = partF.interp_2d(px,py,simul.x,ng,px.shape[0],i0,j0)
    plat = partF.interp_2d(px,py,simul.y,ng,px.shape[0],i0,j0)

    return plon,plat

###################################################################################
# Some analytical veloticy field
#######################################################
 
  
#@profile   
def ana_vel(nx,ny,nz,dxyz=[1.,1.,1.],flow=[0,1,0,0],norm=[1.,1.,0.],timing=False,config='rot',**kwargs):  
    

    if timing: tstart2 = tm.time()   
    
    [div,rot,S1,S2] = flow
    [u0,v0,w0] = norm
    [dx,dy,dz] = dxyz
    
    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = [0,ny,0,nx]
        
    [ny1,ny2,nx1,nx2] = coord;
    
    if 'xyz' in  kwargs:
        x,y,z =  kwargs['xyz']
    else:
        x,y,z = np.mgrid[0:nx,0:ny,0:nz+1]
        x,y,z = x*dx,y*dy,z*dz

########################################################


    if config=='rot':
        ########################################################
        ##Solid body rotation
        
        x0,y0 = nx/2.,ny/2.
        #u = u0*(0.5*S1*(x-x0)+0.5*S2*(y-y0)+0.5*div*(x-x0)-rot/2.*(y-y0))*np.pi/20.
        #v = v0*(-0.5*S1*(y-y0)+0.5*S2*(x-x0)+0.5*div*(y-y0)+rot/2.*(x-x0))*np.pi/20.
        
        u = -0.5*(y-y0)*np.pi/4.
        v = 0.5*(x-x0)*np.pi/4.  
        
        w = z*w0
    
    elif config=='rotexp': 
        ########################################################
        ##Solid body rotation with decreasing exp
        
        x0,y0 = nx/2,ny/2
        r = np.sqrt((x-x0)**2 + (y-y0)**2)
        r0=10
        
        u = u0*(0.5*S1*(x-x0)+0.5*S2*(y-y0)+0.5*div*(x-x0)-rot/2.*(y-y0))*np.exp(-r/r0)*np.pi/4.
        v = v0*(-0.5*S1*(y-y0)+0.5*S2*(x-x0)+0.5*div*(y-y0)+rot*(x-x0)-rot/2.*(x-x0))*np.exp(-r/r0)*np.pi/4.
        w = z*w0   
        
    
    elif config=='birotexp': 
        #######################################################
        #double rotation with decreasing exp
        
        x0,y0 = nx/4,ny/2
        r = np.sqrt((x-x0)**2 + (y-y0)**2)
        r0=5
        
        u1 = u0*(0.5*S1*(x-x0)+0.5*S2*(y-y0)+0.5*div*(x-x0)-rot/2.*(y-y0))*np.exp(-r/r0)
        v1 = v0*(-0.5*S1*(y-y0)+0.5*S2*(x-x0)+0.5*div*(y-y0)+rot*(x-x0)-rot/2.*(x-x0))*np.exp(-r/r0)
        w1 = z*w0   
        
        x0,y0 = 3*nx/4,ny/2  
        r = np.sqrt((x-x0)**2 + (y-y0)**2)
        r0=5
        
        #rot=-rot
        u2 = u0*(0.5*S1*(x-x0)+0.5*S2*(y-y0)+0.5*div*(x-x0)-rot/2.*(y-y0))*np.exp(-r/r0)
        v2 = v0*(-0.5*S1*(y-y0)+0.5*S2*(x-x0)+0.5*div*(y-y0)+rot*(x-x0)-rot/2.*(x-x0))*np.exp(-r/r0)
        w2 = z*w0      
        
        u = u1 + u2
        v = v1 + v2
        w = w1 + w2



    elif config=='front':
        ########################################################
        ##Jeroen's front
        u0,v0,w0,z0 = np.zeros((nx-1,ny-1,nz+1)),np.zeros((nx-1,ny-1,nz+1)),np.zeros((nx-1,ny-1,nz+1)),np.zeros((nx-1,ny-1,nz+1))

        
        for i in range(nx-1):
          for j in range(ny-1):                           
            for k in range(nz+1):  
              u0[i,j,k],v0[i,j,k],w0[i,j,k],z0[i,j,k]= partF.ana_front(i,j,k,nx,ny,nz)
              
              

    elif config=='jet':
        ########################################################
     
        x0,y0 = nx/2.,ny/2.

        u = 0.*x + 1.
        v = 2.*(-1)**x
        
        w = z*w0

    #######################################################
    # Put u,v on u,v grids
    if config=='front':
        u,v,w,dz = np.zeros((nx-1,ny,nz)),np.zeros((nx,ny-1,nz)),np.zeros((nx,ny,nz+1)),np.zeros((nx,ny,nz))

        u[:,1:-1,:] = 0.25*(u0[:,1:,1:] + u0[:,:-1,1:] + u0[:,1:,:-1] + u0[:,:-1,:-1])
        v[1:-1,:,:] = 0.25*(v0[1:,:,1:] + v0[:-1,:,1:] + v0[1:,:,:-1] + v0[:-1,:,:-1])  
        w[1:-1,1:-1,:] = 0.25*(w0[1:,1:,:] + w0[:-1,1:,:] + w0[1:,:-1,:] + w0[:-1,:-1,:])  
        
        dz[1:-1,1:-1,:] = 0.25*(z0[1:,1:,1:] + z0[:-1,1:,1:] + z0[1:,:-1,1:] + z0[:-1,:-1,1:]) \
             -  0.25*(z0[1:,1:,:-1] + z0[:-1,1:,:-1] + z0[1:,:-1,:-1] + z0[:-1,:-1,:-1])
        
    else:
        u = 0.25*(u[1:,:,1:] + u[:-1,:,1:] + u[1:,:,:-1] + u[:-1,:,:-1])
        v = 0.25*(v[:,1:,1:] + v[:,:-1,1:] + v[:,1:,:-1] + v[:,:-1,:-1])   

        u = u[nx1:nx2-1,ny1:ny2,:]
        v = v[nx1:nx2,ny1:ny2-1,:]
        w = w[nx1:nx2,ny1:ny2,:]
        
        dz = np.ones((nx,ny,nz))*dz
        
    #######################################################

    if timing: print 'create  u,v analytically....', tm.time()-tstart2
    if timing: tstart2 = tm.time()    

    return u,v,w,dz
    
    
#######################################################
#######################################################
 
#Some routines copied from the R_tools modules:

#######################################################
#######################################################


#######################################################
#Transfert a field at u points to the rho points
#######################################################

def v2rho(var_v):


    if np.ndim(var_v)<3:
        var_rho = v2rho_2d(var_v)
    else:
        var_rho = v2rho_3d(var_v)

    return var_rho

#######################################################

def v2rho_2d(var_v):

    [Mp,L]=var_v.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[:,1:L]=0.5*(var_v[:,0:Lm]+var_v[:,1:L])
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Lp-1]=var_rho[:,L-1]
    return var_rho

#######################################################

def v2rho_3d(var_v):

    [Mp,L,N]=var_v.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((Mp,Lp,N))
    var_rho[:,1:L,:]=0.5*(var_v[:,0:Lm,:]+var_v[:,1:L,:])
    var_rho[:,0,:]=var_rho[:,1,:]
    var_rho[:,Lp-1,:]=var_rho[:,L-1,:]
    return var_rho


#######################################################
#Transfert a 2 or 2-D field at u points to the rho points
#######################################################

def u2rho(var_u):

    if np.ndim(var_u)<3:
        var_rho = u2rho_2d(var_u)
    else:
        var_rho = u2rho_3d(var_u)

    return var_rho

#######################################################

def u2rho_2d(var_u):

    [M,Lp]=var_u.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[1:M,:]=0.5*(var_u[0:Mm,:]+var_u[1:M,:])
    var_rho[0,:]=var_rho[1,:]
    var_rho[Mp-1,:]=var_rho[M-1,:]

    return var_rho

#######################################################

def u2rho_3d(var_u):

    [M,Lp,N]=var_u.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Mp,Lp,N))
    var_rho[1:M,:,:]=0.5*(var_u[0:Mm,:]+var_u[1:M,:,:])
    var_rho[0,:,:]=var_rho[1,:,:]
    var_rho[Mp-1,:,:]=var_rho[M-1,:,:]

    return var_rho
    
    
    
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
            print 'using NEW_S_COORD'
            (z_r,z_w) = partF.zlevs_croco_new(topo, zeta, simul.hc, simul.Cs_r, simul.Cs_w, simul.sc_r, simul.sc_w)
        else:
            print 'using OLD_S_COORD'
            (z_r,z_w) = partF.zlevs_croco_old(topo, zeta, simul.hc, simul.Cs_r, simul.Cs_w, simul.sc_r, simul.sc_w)


    return [z_r,z_w]
    
       
#################################################
# get_depths (from setdepth.F in romsucla)
#################################################


def get_depths_w(simul,x_periodic=False,y_periodic=False,ng=0,**kwargs):


    if 'coord' in  kwargs: 
        coord = kwargs['coord']
    else: 
        coord = simul.coord
        
    [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
    [ny1,ny2,nx1,nx2] = simul.coord[0:4]
    topo = periodize2d_fromvar(simul,simul.topo,coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng) 

    if hasattr(simul, 'zeta'):
        zeta=periodize2d_fromvar(simul,simul.zeta,coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng) 
    else: 
        zeta=periodize2d_fromnc(simul,'zeta',coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng) 

    if simul.ncname.model=='ucla':
        (z_w) = partF.zlevs_w(topo, zeta, simul.hc,  simul.Cs_w)
    else:
        if simul.VertCoordType == 'NEW':
            print 'using NEW_S_COORD'
            (z_w) = partF.zlevs_croco_new_w(topo, zeta, simul.hc, simul.Cs_w, simul.sc_w)
        else:
            print 'using OLD_S_COORD'
            (z_w) = partF.zlevs_croco_old_w(topo, zeta, simul.hc, simul.Cs_w, simul.sc_w)

    return z_w


#######################################################
#interpolate a 3D variable on horizontal levels of constant depths
#######################################################


def vinterp(var, depths, z_r, z_w=None, mask=None,imin=0,jmin=0,kmin=1, floattype=np.float64,interp=0):


    if mask is None:  mask = np.ones((z_r.shape[0],z_r.shape[1]), order='F', dtype=floattype); mask[z_r[:,:,-1]==0] = 0

    if z_w is None: 
        print 'no z_w specified'
        z_w=np.zeros((z_r.shape[0],z_r.shape[1],z_r.shape[2]+1), order='F')

    if np.ndim(depths)==1: newz = np.asfortranarray(np.zeros((z_r.shape[0],z_r.shape[1],len(depths))) + depths, dtype=floattype)
    else: newz = depths

    if interp==1:
        #print "data will be interpolated below ground up to topo-below"
        below=100 #depths[1]-depths[0]
        vnew=partF.sigma_to_z_intr_bot(z_r, z_w,mask,var,newz,below,imin,jmin,kmin,0.)
    else:
        #print "no interpolation below ground"   
        vnew=partF.sigma_to_z_intr_sfc(z_r, z_w,mask,var,newz,imin,jmin,kmin,0.)

    return vnew

      
#######################################################
#Find indices corresponding to lon,lat (nearest neighbours)
#######################################################


#import matplotlib.nxutils as nxutils 
import scipy.spatial as sp

def find_points(x,y,lon,lat):
    '''
    return indices of the nearest points
    '''

    if isinstance(lon,int) or isinstance(lon,float):
        lon=np.array([lon]); lat=np.array([lat])

    # Shoe-horn existing data for entry into KDTree routines
    combined_x_y_arrays = np.dstack([y.ravel(),x.ravel()])[0]
    points_list = np.dstack([lat.ravel(),lon.ravel()])[0]

    def do_kdtree(combined_x_y_arrays,points):
        mytree = sp.cKDTree(combined_x_y_arrays)
        dist, indexes = mytree.query(points)
        return indexes
    results = do_kdtree(combined_x_y_arrays,points_list)

    i,j = results/x.shape[1], results%x.shape[1]

    print i,j

    return i[0],j[0]



#######################################################
#Transfert a field at rho points to psi points
#######################################################

def rho2psi(var_rho):

    if np.ndim(var_rho)<3:
        var_psi = rho2psi_2d(var_rho)
    else:
        var_psi = rho2psi_3d(var_rho)

    return var_psi


##############################

def rho2psi_2d(var_rho):

    var_psi = 0.25*(var_rho[1:,1:]+var_rho[1:,:-1]+var_rho[:-1,:-1]+var_rho[:-1,1:])

    return var_psi

#############################

def rho2psi_3d(var_rho):

    var_psi = 0.25*(var_rho[1:,1:,:]+var_rho[1:,:-1,:]+var_rho[:-1,:-1,:]+var_rho[:-1,1:,:])

    return var_psi


    
