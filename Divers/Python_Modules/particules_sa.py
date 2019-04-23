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
import particules_3d_sa as partF


###################################################################################






###################################################################################
#@profile
def subsection(px,py,dx=1,maxvel=1,delt=1,nx=2000,ny=2000,**kwargs):
    """ Finds index subrange to move particles around"""

    if 'offset' in kwargs:
       offset = kwargs['offset']
    else:
       offset = np.abs(maxvel*delt/dx)
       
    i0 = int(np.floor(np.nanmin(px) - offset))
    i0 = max(0,i0)
    i1 = int(np.ceil(np.nanmax(px)  + offset))
    i1 = min(nx,i1)
    j0 = int(np.floor(np.nanmin(py) - offset))
    j0 = max(0,j0)
    j1 = int(np.ceil(np.nanmax(py)  + offset))
    j1 = min(ny,j1)

    return [j0,j1,i0,i1]
   
###################################################################################

   
#@profile   
def cull(px,py,pz,nx,ny,nz):
    """ Set particle positions that are out of range to nan"""

    px[px<=0] = np.nan; px[px>=nx-2] = np.nan; 
    py[py<=0] = np.nan; py[py>=ny-2] = np.nan; 
    pz[pz<=0] = np.nan; pz[pz>=nz-1.01] = nz-1.01; 

    py[np.isnan(px)] = np.nan
    pz[np.isnan(px)] = np.nan
    px[np.isnan(py)] = np.nan
    pz[np.isnan(py)] = np.nan
    px[np.isnan(pz)] = np.nan
    py[np.isnan(pz)] = np.nan

    return [px,py,pz]
   
   
###################################################################################

   
#@profile   
def cull_topo(simul,px,py,pz,**kwargs):
    """ Set particle positions that are out of range to nan"""

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    if 'depths' in kwargs:
        depths =  kwargs['depths']
    else:
        depths = simul.coord[4]

    ptopo = map_topo(simul,px,py,coord=coord)   
    pdepth = map_depth(simul,pz,depths=depths)   
    
    N=np.nansum(pdepth<-1*ptopo)
    pz[pdepth<-1*ptopo] = np.nan
    px[np.isnan(pz)] = np.nan
    py[np.isnan(pz)] = np.nan

    return [px,py,pz,N]
   
   
   
   
###################################################################################

    
#@profile   
def get_vel(simul,timing=False,**kwargs):   

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]   

    if 'depths' in kwargs:
        depths =  kwargs['depths']
    else:
        depths = simul.coord[4]

    u0,v0 = get_vel_io(simul,timing,coord=coord)
    
    u,v,w = get_vel_comp(simul,u0,v0,timing,coord=coord, depths=depths)
    
    return u,v,w
   
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

    u = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1:ny2,nx1:nx2-1]))
    v = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny1:ny2-1,nx1:nx2]))
    
    nc.close()
    
    if timing: print('get u,v from file....', tm.time()-tstart2)
    if timing: tstart2 = tm.time()    

    return u,v


###################################################################################
  
#@profile   
def get_vel_comp(simul,u,v,timing=False,**kwargs):  

    
    if 'coord' in  kwargs:
        coord = kwargs['coord']
        [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
        [ny1,ny2,nx1,nx2] = simul.coord[0:4]
        pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
    else: 
        coord = simul.coord[0:4]
        pm = simul.pm; pn = simul.pn
        
    if 'depths' in kwargs:
        depths =  kwargs['depths']
    else:
        depths = simul.coord[4]
        
    if timing: tstart2 = tm.time()     

    [z_r,z_w] = get_depths(simul,coord=coord)

    if timing: print('compute depth........', tm.time()-tstart2)
    if timing: tstart2 = tm.time()       

    w = partF.get_wvlcty(u,v,z_r,z_w,pm,pn)
    w[0,:,:] = np.nan; w[-1,:,:] = np.nan
    w[:,0,:] = np.nan; w[:,-1,:] = np.nan
    
    if timing: print('compute w............', tm.time()-tstart2)
    if timing: tstart2 = tm.time()      
    
    u0 = np.asfortranarray(u2rho(u),dtype= simul.floattype)
    v0 = np.asfortranarray(v2rho(v),dtype = simul.floattype)
    w0 = np.asfortranarray(w,dtype = simul.floattype)

    #if timing: print 'vertical interpolat..', tm.time()-tstart2   

    u = vinterp(u0,depths,z_r,z_w,interp=1)
    v = vinterp(v0,depths,z_r,z_w,interp=1) 
    w = vinterp(w0,depths,z_r,z_w,interp=1)
    
    if timing: print('vertical interpolat..', tm.time()-tstart2, w.shape)
    
    return u,v,w

    
    
###################################################################################
    
#@profile   
def get_ts(simul,**kwargs):   

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]   
        
    temp0,salt0 = get_ts_io(simul,coord=coord)
    
    temp,salt = get_ts_comp(simul,temp0,salt0,coord=coord)
    
    return temp,salt
    
    
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
    try:
        salt = np.squeeze(simul.Forder(nc.variables['salt'][simul.infiletime,:,ny1:ny2,nx1:nx2])) 
    except:
        salt = np.zeros(temp.shape)

    nc.close()

    return temp,salt

###################################################################################

    
#@profile
def get_ts_comp(simul,temp0,salt0,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    depths = simul.coord[4]
    
    [z_r,z_w] = get_depths(simul,coord=coord)
    
    temp = np.asfortranarray(vinterp(temp0,depths,z_r,z_w,interp=1))   
    
    salt = np.asfortranarray(vinterp(salt0,depths,z_r,z_w,interp=1))

    return temp,salt
    

###################################################################################
#@profile
def map_ts(simul,temp,salt,px,py,pz,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    [j0,j1,i0,i1]=coord; k0=0; nqmx=px.shape[0]
    
    [ptemp,psalt] = partF.interp_3d(px,py,pz,temp,salt,nqmx,i0,j0,k0)

    return [ptemp,psalt]

###################################################################################
#@profile
def map_var(simul,var0,px,py,pz,**kwargs):  

    if 'coord' in  kwargs:
        coord = kwargs['coord']
    else: 
        coord = simul.coord[0:4]
        
    [j0,j1,i0,i1]=coord; k0=0; nqmx=px.shape[0]
    
    pvar0 = partF.oneterp_3d(px,py,pz,var0,nqmx,i0,j0,k0)

    return pvar0
    
    
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
    
    
##################################################################################
#@profile
def map_depth(simul,pz,**kwargs):  
  
    
    if 'depths' in kwargs:
        depths =  np.asfortranarray(kwargs['depths'],dtype=np.float64)
    else:
        depths = np.asfortranarray(simul.coord[4],dtype=np.float64)
        
    pdepth = partF.interp_1d(pz,depths,pz.shape[0])

    return pdepth
    
    
    
    
    
    
    
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
  
    
    
