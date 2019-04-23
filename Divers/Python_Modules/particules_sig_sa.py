###################################################################################
# Particules routines
###################################################################################
"""





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

#particules routines
import particules_3d_sig_sa as partF

#copy data
from copy import copy
###################################################################################






###################################################################################
#@profile
def subsection(px,py,dx=1,maxvel=[1,1],delt=1,nx=2000,ny=2000,**kwargs):
    """ Finds index subrange to move particles around"""

    if 'offset' in kwargs:
       offset_x = kwargs['offset']
       offset_y = kwargs['offset']
    else:
       offset_x = np.abs(maxvel[0]*delt/dx)
       offset_y = np.abs(maxvel[1]*delt/dx)     
       
    i0 = int(np.floor(np.nanmin(px) - offset_x))
    i0 = max(0,i0)
    i1 = int(np.ceil(np.nanmax(px)  + offset_x))
    i1 = min(nx,i1)
    j0 = int(np.floor(np.nanmin(py) - offset_y))
    j0 = max(0,j0)
    j1 = int(np.ceil(np.nanmax(py)  + offset_y))
    j1 = min(ny,j1)

    return [j0,j1,i0,i1]

###################################################################################

   
#@profile   
def cull(px,py,pz,nx,ny,nz,periodic=False):
    """ Set particle positions that are out of range to nan"""

    if periodic:
        px[px<-0.5] = px[px<-0.5] + nx
        px[px>=nx-0.5] = px[px>=nx-0.5] - nx
    else:
        px[px<=1] = np.nan; px[px>=nx-3] = np.nan; 
        py[py<=1] = np.nan; py[py>=ny-3] = np.nan; 
        #pz[pz<0] = np.nan; pz[pz>nz] = np.nan; 

    py[np.isnan(px)] = np.nan
    px[np.isnan(py)] = np.nan
    pz[np.isnan(px)] = np.nan
    
    #pz[np.isnan(py)] = np.nan
    #px[np.isnan(pz)] = np.nan
    #py[np.isnan(pz)] = np.nan

    return [px,py,pz]
   
   
###################################################################################

   
#@profile   
def kick(pz,nz):
    """ Give a kick to particules trapped at surface/bottom  """

    eps = 0.01

    pz[pz<=eps] = eps; pz[pz>=nz-eps] = nz-eps; 


    return [pz]
   
###################################################################################
  
#@profile   
def get_vel_io(simul,timing=False,**kwargs):  


    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    if timing: tstart2 = tm.time()     
  
    nc = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
    [ny1,ny2,nx1,nx2] = coord

    if simul.ncname.model == 'agrif_jc':
        u = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1:ny2,nx1+1:nx2]))
    else:
        u = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1:ny2,nx1:nx2-1]))


    u = (u.T * 0.5*(simul.mask[nx1+1:nx2,ny1:ny2]+simul.mask[nx1:nx2-1,ny1:ny2]).T).T
    if timing: print('get u from file....', tm.time()-tstart2)
    if timing: tstart2 = tm.time()    

    v = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny1:ny2-1,nx1:nx2]))
    v = (v.T * 0.5*(simul.mask[nx1:nx2,ny1+1:ny2]+simul.mask[nx1:nx2,ny1:ny2-1]).T).T
    if timing: print('get v from file....', tm.time()-tstart2)
    if timing: tstart2 = tm.time()    

    try:
        w = simul.Forder(np.squeeze(nc.variables['omega'][simul.infiletime,:,ny1:ny2,nx1:nx2]))
    except:
        print('no omega in file, computing')
        [z_r,z_w] = get_depths(simul,coord=coord)
        w = partF.get_omega(u,v,z_r,z_w,simul.pm[nx1:nx2,ny1:ny2],simul.pn[nx1:nx2,ny1:ny2])

    w = (w.T * simul.mask[nx1:nx2,ny1:ny2].T).T
    nc.close()
    
    if timing: print('get w from file....', tm.time()-tstart2)
    if timing: tstart2 = tm.time()    

    return u,v,w



###################################################################################
   
    
#@profile
def get_ts_io(simul,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]

    nc = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
    [ny1,ny2,nx1,nx2] = coord; 

    temp = np.squeeze(simul.Forder(nc.variables['temp'][simul.infiletime,:,ny1:ny2,nx1:nx2]))
    temp = (temp.T * simul.mask.T[ny1:ny2,nx1:nx2]).T

    try:
        salt = np.squeeze(simul.Forder(nc.variables['salt'][simul.infiletime,:,ny1:ny2,nx1:nx2]))
        salt = (salt.T * simul.mask.T[ny1:ny2,nx1:nx2]).T
    except:
        salt = np.zeros(temp.shape)

    nc.close()

    return temp,salt


###################################################################################
# Compute T,S at each particle position
###################################################################################
#@profile


def map_ts(simul,temp,salt,px,py,pz,**kwargs):  
    '''
    T,S on horizontal and vertical rho-grid
    '''

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    [j0,j1,i0,i1]=coord; k0=0; nqmx=px.shape[0]
    
    [ptemp,psalt] = partF.interp_3d_ts(px,py,pz,temp,salt,nqmx,i0,j0,k0)

    return [ptemp,psalt]


###################################################################################
# Compute a 3D variable at each particle position
###################################################################################
#@profile


def map_var(simul,var0,px,py,pz,**kwargs):  
    '''
    var on horizontal and vertical rho-grid
    '''

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    [j0,j1,i0,i1]=coord; k0=0; nqmx=px.shape[0]
   
    print(var0.shape)
    print([i1-i0,j1-j0,len(simul.coordmax[4])])
   
    if var0.shape==(i1-i0,j1-j0,len(simul.coordmax[4])):
        pvar0 = partF.interp_3d(px,py,pz,var0,nqmx,i0,j0,k0)
    elif var0.shape==(i1-i0-1,j1-j0-1,len(simul.coordmax[4])):
        pvar0 = partF.interp_3d_psi(px,py,pz,var0,nqmx,i0,j0,k0)
    elif var0.shape==(i1-i0,j1-j0,len(simul.coordmax[4])+1):
        pvar0 = partF.interp_3d_w(px,py,pz,var0,nqmx,i0,j0,k0)
    elif var0.shape==(i1-i0-1,j1-j0-1,len(simul.coordmax[4])+1):
        pvar0 = partF.interp_3d_psiw(px,py,pz,var0,nqmx,i0,j0,k0)        
        

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


def map_topo(simul,px,py,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    [j0,j1,i0,i1]=coord;
    
    ptopo = partF.interp_2d(px,py,simul.topo,px.shape[0],i0,j0)

    return ptopo
    
    
    
###################################################################################
# Compute lon,lat at each particle position
###################################################################################
#@profile

def map_lonlat(simul,px,py,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    [j0,j1,i0,i1]=coord;
    
    plon = partF.interp_2d(px,py,simul.x,px.shape[0],i0,j0)
    plat = partF.interp_2d(px,py,simul.y,px.shape[0],i0,j0)

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

    if timing: print('create  u,v analytically....', tm.time()-tstart2)
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


    if np.rank(var_v)<3:
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

    if np.rank(var_u)<3:
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


def get_depths(simul,**kwargs):


    if 'coord' in  kwargs: 
        coord = kwargs['coord']
    else: 
        coord = simul.coord
        
    [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
    [ny1,ny2,nx1,nx2] = simul.coord[0:4]
    topo = np.asfortranarray(simul.topo[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        
    if hasattr(simul, 'zeta'): 
        zeta=simul.zeta[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1]
    else: 
        nc = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
        zeta = simul.Forder(np.squeeze(nc.variables['zeta'][simul.infiletime,ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1]))
        nc.close()

    (z_r,z_w) = partF.zlevs(topo, zeta, simul.hc, simul.Cs_r, simul.Cs_w)
        
    return [z_r,z_w]
    
       
#################################################
# get_depths (from setdepth.F in romsucla)
#################################################


def get_depths_w(simul,**kwargs):


    if 'coord' in  kwargs: 
        coord = kwargs['coord']
    else: 
        coord = simul.coord
        
    [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
    [ny1,ny2,nx1,nx2] = simul.coord[0:4]
    topo = np.asfortranarray(simul.topo[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        
    if hasattr(simul, 'zeta'): 
        zeta=simul.zeta[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1]
    else: 
        nc = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
        zeta = simul.Forder(np.squeeze(nc.variables['zeta'][simul.infiletime,ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1]))
        nc.close()

    (z_w) = partF.zlevs_w(topo, zeta, simul.hc,  simul.Cs_w)
        
    return z_w


#######################################################
#interpolate a 3D variable on horizontal levels of constant depths
#######################################################


def vinterp(var, depths, z_r, z_w=None, mask=None,imin=0,jmin=0,kmin=1, floattype=np.float64,interp=0):


    if mask==None:  mask = np.ones((z_r.shape[0],z_r.shape[1]), order='F', dtype=floattype); mask[z_r[:,:,-1]==0] = 0

    if z_w==None: 
        print('no z_w specified')
        z_w=np.zeros((z_r.shape[0],z_r.shape[1],z_r.shape[2]+1), order='F')

    if np.rank(depths)==1: newz = np.asfortranarray(np.zeros((z_r.shape[0],z_r.shape[1],len(depths))) + depths, dtype=floattype)
    else: newz = depths

    if interp==1:
        #print "data will be interpolated below ground up to topo-below"
        below=100 #depths[1]-depths[0]
        vnew=partF.sigma_to_z_intr_bot(z_r, z_w,mask,var,newz,below,imin,jmin,kmin,9999.)
    else:
        #print "no interpolation below ground"   
        vnew=partF.sigma_to_z_intr_sfc(z_r, z_w,mask,var,newz,imin,jmin,kmin,9999.)

    vnew[np.abs(vnew)==9999.]=np.nan

    return vnew


####################################################### 
  
    
    
