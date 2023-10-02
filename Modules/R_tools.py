###################################################################################
# R_TOOLS
###################################################################################
"""





"""

###################################################################################
#Load modules
###################################################################################

#for numeric functions
import numpy as np

#for netcdf files
#from Scientific.IO.NetCDF import *
from netCDF4 import Dataset

#copy data
from copy import copy

#ROMSTOOLS
import R_tools_fort as toolsF

#Simulations (path, data...)
import R_vars as va

#for plotting
import matplotlib.pyplot as py

import time as tm

###################################################################################


import R_smooth as sm


#################################################
# get_depths (from setdepth.F in romsucla)
#################################################


def get_depths(simul,**kwargs):


    if 'coord' in  kwargs:
        coord = kwargs['coord']
        [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
        [ny1,ny2,nx1,nx2] = simul.coord[0:4]
        print(ny1i,ny2i,nx1i,nx2i)
        print(ny1,ny2,nx1,nx2)
        topo = np.asfortranarray(simul.topo[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        #topo = np.asfortranarray(simul.topo[nx1i:nx2i,ny1i:ny2i])
    else: 
        coord = simul.coord
        topo = simul.topo
        
    
    if hasattr(simul, 'zeta'): 
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]      
            zeta=simul.zeta[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1]
            #zeta=simul.zeta[nx1i:nx2i,ny1i:ny2i]
        else:
            zeta=simul.zeta
        print('using simul.zeta')
        #print 'simul.zeta[5,82]',zeta[5,82]
    else: 
        zeta = va.var('zeta',simul,depths=[0],coord=coord[0:4],masked=0.).data
        #print 'zeta[5,82]',zeta[5,82]
    
    if 'agrif' in kwargs:
        (z_r,z_w) = toolsF.zlevs_agrif(topo, zeta, simul.hc, simul.Cs_r, simul.Cs_w, simul.sc_r, simul.sc_w)
        print('using AGRIF version')
    elif simul.VertCoordType == 'KAU':
        (z_r,z_w) = toolsF.zlevs_kau(topo, simul.hc, simul.Cs_r, simul.Cs_w)
        print('using KAU version')
    else:
        (z_r,z_w) = toolsF.zlevs(topo, zeta, simul.hc, simul.Cs_r, simul.Cs_w)
        #print 'zw[5,82]',z_w[5,82,:]

    if 'sub' in  kwargs: 
        z_r = np.asfortranarray(z_r[:,:,simul.coord[4]-1])
        z_w = np.asfortranarray(z_w[:,:,np.arange(np.min(simul.coord[4])-1,np.max(simul.coord[4])+1)])
    elif 'depths' in kwargs:
        depths = np.array(kwargs['depths'])-1
        z_r = np.asfortranarray(z_r[:,:,depths])
        z_w = np.asfortranarray(z_w[:,:,np.arange(np.min(depths),np.max(depths)+2)])

     
    return [z_r,z_w]



#################################################
# get_depth (from setdepth.F in romsucla)
#################################################


def get_depth(simul,**kwargs):


    if 'coord' in  kwargs: 
        coord = kwargs['coord']
        [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
        [ny1,ny2,nx1,nx2] = simul.coord[0:4]
        topo = np.asfortranarray(simul.topo[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        #topo = np.asfortranarray(simul.topo[nx1i:nx2i,ny1i:ny2i])

    else: 
        coord = simul.coord
        topo = simul.topo

    if hasattr(simul, 'zeta'): 
        #zeta=simul.zeta[nx1i:nx2i,ny1i:ny2i] 
        zeta=simul.zeta[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1]
    else: 
        zeta = va.var('zeta',simul,depths=[0],coord=coord[0:4]).data

    (z_r) = toolsF.zlev(topo, zeta, simul.hc, simul.Cs_r, simul.Cs_w)

    return z_r



#######################################################
#interpolate a 3D variable on horizontal levels of constant depths (FORTRAN version, much faster)
#######################################################


def vinterp(var, depths, z_r, z_w=None, mask=None,imin=0,jmin=0,kmin=1, floattype=np.float64,interp_sfc=1,interp_bot=0,below=None,verbo=False,bounded=False,**kwargs):

    '''
    bounded == True,  means that data above surface take surface value and data below bottom take bottom value
    '''
    if mask is None:  mask = np.ones((z_r.shape[0],z_r.shape[1]), order='F', dtype=floattype); mask[z_r[:,:,-1]==0] = 0

    if z_w is None: 
        print('no z_w specified')
        z_w=np.zeros((z_r.shape[0],z_r.shape[1],z_r.shape[2]+1), order='F')
        z_w[:,:,1:-1] = 0.5*(z_r[:,:,1:] + z_r[:,:,:-1])
        z_w[:,:,0] = z_r[:,:,0] - (z_r[:,:,1]-z_r[:,:,0])
        z_w[:,:,-1] = z_r[:,:,-1] + (z_r[:,:,-1]-z_r[:,:,-2])
        
    if np.ndim(depths)==1: newz = np.asfortranarray(np.zeros((z_r.shape[0],z_r.shape[1],len(depths))) + depths, dtype=floattype)
    else: newz = depths

    if bounded:
        if verbo: print("data bounded")
        vnew=toolsF.sigma_to_z_intr_bounded(z_r, z_w,mask,var,newz,imin,jmin,kmin,9999.)
    elif interp_bot==1:
        if verbo: print("data will be interpolated below ground")
        below=10000.
        vnew=toolsF.sigma_to_z_intr_bot(z_r, z_w,mask,var,newz,below,imin,jmin,kmin,9999.)
    elif interp_sfc==1:
        if verbo: print("no interpolation below ground")
        vnew=toolsF.sigma_to_z_intr_sfc(z_r, z_w,mask,var,newz,imin,jmin,kmin,9999.)
    else:
        if verbo: print("no interpolation below ground or above surface")
        vnew=toolsF.sigma_to_z_intr(z_r, z_w,mask,var,newz,imin,jmin,kmin,9999.)   	

    vnew[np.abs(vnew)==9999.]=np.nan

    return vnew



#######################################################
#interpolate a 3D variable on horizontal levels of constant depths (FORTRAN version, much faster)
#######################################################


def vinterp_2d(var, depths, z_r, z_w=None, mask=None,imin=0,jmin=0,kmin=1, floattype=np.float64,below=None,verbo=False,**kwargs):

    
    if mask is None:  mask = np.ones((z_r.shape[0]), order='F', dtype=floattype); mask[z_r[:,-1]==0] = 0

    if z_w is None: 
        print('no z_w specified')
        z_w=np.zeros((z_r.shape[0],z_r.shape[1]+1), order='F')
        z_w[:,1:-1] = 0.5*(z_r[:,1:] + z_r[:,:-1])
        z_w[:,0] = z_r[:,0] - (z_r[:,1]-z_r[:,0])
        z_w[:,-1] = z_r[:,-1] + (z_r[:,-1]-z_r[:,-2])
        
    if np.ndim(depths)==1: newz = np.asfortranarray(np.zeros((z_r.shape[0],len(depths))) + depths, dtype=floattype)
    else: newz = depths


    below=0.
    vnew=toolsF.sigma_to_z_intr_bot_2d(z_r, z_w,mask,var,newz,below,imin,kmin,9999.)

    vnew[np.abs(vnew)==9999.]=np.nan

    return vnew








#######################################################
#Transfert a field at psi points to rho points
#######################################################

def psi2rho(var_psi):

    if np.ndim(var_psi)<3:
        var_rho = psi2rho_2d(var_psi)
    else:
        var_rho = psi2rho_3d(var_psi)

    return var_rho


##############################

def psi2rho_2d(var_psi):

    [M,L]=var_psi.shape
    Mp=M+1
    Lp=L+1
    Mm=M-1
    Lm=L-1

    var_rho=np.zeros((Mp,Lp))
    var_rho[1:M,1:L]=0.25*(var_psi[0:Mm,0:Lm]+var_psi[0:Mm,1:L]+var_psi[1:M,0:Lm]+var_psi[1:M,1:L])
    var_rho[0,:]=var_rho[1,:]
    var_rho[Mp-1,:]=var_rho[M-1,:]
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Lp-1]=var_rho[:,L-1]

    return var_rho

#############################

def psi2rho_3d(var_psi):

    [Mz,Lz,Nz]=var_psi.shape
    var_rho=np.zeros((Mz+1,Lz+1,Nz))

    for iz in range(0, Nz, 1):    
        var_rho[:,:,iz]=psi2rho_2d(var_psi[:,:,iz])

    return var_rho



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





#######################################################
#Transfert a field at rho points to u points
#######################################################

def rho2u(var_rho):

    if np.ndim(var_rho)==1:
        var_u = 0.5*(var_rho[1:]+var_rho[:-1])
    elif np.ndim(var_rho)==2:       
        var_u = rho2u_2d(var_rho)
    else:
        var_u = rho2u_3d(var_rho)

    return var_u


##############################

def rho2u_2d(var_rho):

    var_u = 0.5*(var_rho[1:,:]+var_rho[:-1,:])

    return var_u

#############################

def rho2u_3d(var_rho):

    var_u = 0.5*(var_rho[1:,:,:]+var_rho[:-1,:,:])

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

def rho2v_2d(var_rho):

    var_v = 0.5*(var_rho[:,1:]+var_rho[:,:-1])

    return var_v

#############################

def rho2v_3d(var_rho):

    var_v = 0.5*(var_rho[:,1:,:]+var_rho[:,:-1,:])

    return var_v





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




#######################################################
#Transfert a 3-D field from verical w points to vertical rho-points
#######################################################

def w2rho(var_w):

    [M,L,N]=var_w.shape
    print('[M,L,N]',[M,L,N])
    
    var_rho = np.zeros((M,L,N-1))
    
    for iz in range(1,N-2):
        var_rho[:,:,iz]  = 0.5625*(var_w[:,:,iz+1] + var_w[:,:,iz]) -0.0625*(var_w[:,:,iz+2] + var_w[:,:,iz-1])
    
    var_rho[:,:,0]  = -0.125*var_w[:,:,2] + 0.75*var_w[:,:,1] +0.375*var_w[:,:,0] 
    var_rho[:,:,N-2]  = -0.125*var_w[:,:,N-3] + 0.75*var_w[:,:,N-2] +0.375*var_w[:,:,N-1] 
    
    return var_rho









#################################################
# rho_eos (from rho_eos.F in romsucla)
#################################################


def rho_eos(T,S,z_r,z_w,rho0):
    '''
    Compute rho(T,S,z)
    '''

    (rho) = toolsF.rho_eos(T,S,z_r,z_w,rho0)

    return rho

#################################################
# rho_eos (from rho_eos.F in romsucla)
#################################################


def rho1_eos(T,S,z_r,z_w,rho0):
    '''
    Compute rho1(T,S) the sea-water density perturbation[kg/m^3] at
    standard pressure of 1 Atm (sea surface)
    '''
    
    (rho1) = toolsF.rho1_eos(T,S,z_r,rho0)

    return rho1


#################################################
# rho_grad (from rho_eos.F and prsgrd.F in romsucla)
#################################################


def rho_grad(T,S,z_r,z_w,rho0,pm,pn):
    '''
    Compute neutral density gradients
    '''
    
    (drdz,drdx,drdy) = toolsF.rho_grad(T,S,z_r,z_w,rho0)

    return [drdz,drdx,drdy]














#################################################
# Rotationnel (on vertical levels)
#################################################


def rot(u,v,pm,pn,mask=None):
    
    if np.ndim(u)==3:
        for iz in range(u.shape[2]):
    
            rot = rot_2d(u[:,:,iz],v[:,:,iz],pm,pn)
    else:

        rot = rot_2d(u,v,pm,pn)

    return rot


################

def rot_2d(u,v,pm,pn,mask=None):
    
    if u.shape == v.shape:
        u=rho2u(u); v=rho2v(v)
    
    (rot) = toolsF.get_rot(u,v,pm,pn)

    return rot







#################################################
# Gradient amplitude (on vertical levels)
#################################################


def grad(psi,pm,pn):

    (grad) = toolsF.get_grad(psi,pm,pn)

    return grad








#################################################
# divergence (on horizontal levels)
#################################################


def div(u,v,pm,pn):

    # compute div
    div = np.zeros(pm.shape)

    div[1:-1,:] = div[1:-1,:] + diffx(u,rho2u(pm))
    div[:,1:-1] = div[:,1:-1] + diffy(v,rho2v(pn))
    div[np.isnan(div)] =0

    return div




#######################################################
#x-derivative from rho-grid to u-grid
#######################################################

def diffx(var,pm,dn=1):

    if np.ndim(var)<3:
        dvardx = diffx_2d(var,pm,dn)
    else:
        dvardx = diffx_3d(var,pm,dn)

    return dvardx

###########################

def diffx_3d(var,pm,dn=1):

    [N,M,L]=var.shape

    dvardx = np.zeros((N-dn,M,L))

    for iz in range(0, L):    
        dvardx[:,:,iz]=diffx_2d(var[:,:,iz],pm,dn)

    return dvardx

###########################

def diffx_2d(var,pm,dn=1):

    if (np.ndim(pm)==2) and (var.shape[0]==pm.shape[0]): 
        dvardx = (var[dn:,:]-var[:-dn,:])*0.5*(pm[dn:,:]+pm[:-dn,:])/dn
    else: 
        dvardx = (var[dn:,:]-var[:-dn,:])*pm/dn

    return dvardx




#######################################################
#y-derivative from rho-grid to v-grid
#######################################################

def diffy(var,pn,dn=1):

    if np.ndim(var)<3: dvardy = diffy_2d(var,pn,dn)
    else: dvardy = diffy_3d(var,pn,dn)

    return dvardy

    #######################

def diffy_3d(var,pn,dn=1):

    [N,M,L]=var.shape
    dvardy = np.zeros((N,M-dn,L))
    for iz in range(0, L): dvardy[:,:,iz]=diffy_2d(var[:,:,iz],pn,dn)

    return dvardy

    #######################


def diffy_2d(var,pn,dn=1):

    if (np.ndim(pn)==2) and (var.shape[1]==pn.shape[1]):
        dvardy = (var[:,dn:]-var[:,:-dn])*0.5*(pn[:,dn:]+pn[:,:-dn])/dn
    else: 
        dvardy = (var[:,dn:]-var[:,:-dn])*pn/dn

    return dvardy



    
    
    

#######################################################
#Compute horizontal derivatives on sigma-levels (1st order)
#######################################################
'''
var on rho-rho grid
dvardxi on psi-rho grid
'''

def diffxi(var,pm,z_r,z_w=None,newz=None,mask=None):

    if np.ndim(var)>2:
        if z_r.shape[2]<=2:
            dvardxi = diffxi_3d_2z(var,pm,z_r,z_w,newz,mask)
        else:
            dvardxi = diffxi_3d(var,pm,z_r,z_w,newz,mask)
    else:
        dvardxi = diffxi_2d(var,pm,z_r,z_w,newz,mask)

    ##############################################

    return dvardxi

#######################################################
#######################################################


def diffxi_3d(var,pm,z_r,z_w=None,newz=None,mask=None):


    if newz is None: newz = 0.5*(z_r[1:,:,:] + z_r[:-1,:,:])
    else: newz = rho2u(newz)

    dvardxi = np.zeros((var.shape[0]-1,var.shape[1],var.shape[2]))

    ##############################################

    varzp = vinterp(var[1:,:,:],newz,z_r[1:,:,:],z_w[1:,:,:],interp_bot=1)
    varzm = vinterp(var[:-1,:,:],newz,z_r[:-1,:,:],z_w[:-1,:,:],interp_bot=1)

    dvardxi = ((varzp - varzm ).T*0.5*(pm[1:,:]+pm[:-1,:]).T ).T

    ##############################################

    return dvardxi


#######################################################
#######################################################


def diffxi_3d_2z(var,pm,z_r,z_w=None,newz=None,mask=None):

    dvardxi = np.zeros((z_r.shape[0]-1,z_r.shape[1]))

    ##############################################

    if newz is None: newz = 0.5*(z_r[:-1,:,0] + z_r[1:,:,0])
    else: newz = rho2u(newz)

    dz0 = (z_r[1:,:,0]-newz)
    dz1 = (newz-z_r[1:,:,1])
    varzp = (dz1*var[1:,:,0] + dz0*var[1:,:,1])/(z_r[1:,:,0]-z_r[1:,:,1])

    dz0 = (z_r[:-1,:,0]-newz)
    dz1 = (newz-z_r[:-1,:,1])
    varzm = (dz1*var[:-1,:,0] + dz0*var[:-1,:,1])/(z_r[:-1,:,0]-z_r[:-1,:,1])

    dvardxi = (varzp - varzm )*0.5*(pm[1:,:]+pm[:-1,:])
    ##############################################

    return dvardxi


#######################################################
#######################################################


def diffxi_2d(var,pm,z_r,z_w=None,newz=None,mask=None):


    if newz is None: newz = 0.5*(z_r[1:,:] + z_r[:-1,:])
    else: newz = rho2u(newz)

    dvardxi = np.zeros((var.shape[0]-1,var.shape[1]))

    ##############################################

    varzp = vinterp_2d(var[1:,:],newz,z_r[1:,:],z_w[1:,:],interp_bot=1)
    varzm = vinterp_2d(var[:-1,:],newz,z_r[:-1,:],z_w[:-1,:],interp_bot=1)

    dvardxi = ((varzp - varzm ).T*0.5*(pm[1:]+pm[:-1]).T ).T

    ##############################################

    return dvardxi



#######################################################
#Compute horizontal derivatives on sigma-levels (1st order)
#######################################################

'''
var on rho-rho grid
dvardxi on psi-rho grid
'''

def diffeta(var,pn,z_r,z_w=None,newz=None,mask=None):

    if np.ndim(var)>2:
        if z_r.shape[2]<=2:
            dvardeta = diffeta_3d_2z(var,pn,z_r,z_w,newz,mask)
        else:
            dvardeta = diffeta_3d(var,pn,z_r,z_w,newz,mask)
    else:
        dvardeta = diffeta_2d(var,pn,z_r,z_w,newz,mask)


    ##############################################

    return dvardeta


#######################################################
#######################################################


def diffeta_3d(var,pn,z_r,z_w=None,newz=None,mask=None):


    if newz is None: newz = 0.5*(z_r[:,:-1,:] + z_r[:,1:,:])
    else: newz = rho2v(newz)

    dvardeta = np.zeros((var.shape[0],var.shape[1]-1,var.shape[2]))

    ##############################################

    varzp = vinterp(var[:,1:,:],newz,z_r[:,1:,:],z_w[:,1:,:],interp_bot=1)
    varzm = vinterp(var[:,:-1,:],newz,z_r[:,:-1,:],z_w[:,:-1,:],interp_bot=1)

    dvardeta = ((varzp - varzm).T*0.5*(pn[:,:-1]+pn[:,1:]).T).T

    ##############################################


    return dvardeta



#######################################################
#Compute horizontal derivatives on sigma-levels (1st order)
#######################################################



def diffeta_3d_2z(var,pn,z_r,z_w=None,newz=None,mask=None):

    dvardeta = np.zeros((z_r.shape[0],z_r.shape[1]-1))

    ##############################################

    if newz is None: newz = 0.5*(z_r[:,:-1,0] + z_r[:,1:,0])
    else: newz = rho2v(newz)

    dz0 = (z_r[:,1:,0]-newz)
    dz1 = (newz-z_r[:,1:,1])
    varzp = (dz1*var[:,1:,0] + dz0*var[:,1:,1])/(z_r[:,1:,0]-z_r[:,1:,1])

    dz0 = (z_r[:,:-1,0]-newz)
    dz1 = (newz-z_r[:,:-1,1])
    varzm = (dz1*var[:,:-1,0] + dz0*var[:,:-1,1])/(z_r[:,:-1,0]-z_r[:,:-1,1])

    dvardeta = (varzp - varzm )*0.5*(pn[:,:-1]+pn[:,1:])

    ##############################################


    return dvardeta

#######################################################
#######################################################


def diffeta_2d(var,pn,z_r,z_w=None,newz=None,mask=None):


    if newz is None: newz = 0.5*(z_r[:-1,:] + z_r[1:,:])
    else: newz = rho2v(newz)

    dvardeta = np.zeros((var.shape[1]-1,var.shape[2]))

    ##############################################

    varzp = vinterp(var[1:,:],newz,z_r[1:,:],z_w[1:,:],interp_bot=1)
    varzm = vinterp(var[:-1,:],newz,z_r[:-1,:],z_w[:-1,:],interp_bot=1)

    dvardeta = ((varzp - varzm).T*0.5*(pn[:,:-1]+pn[:,1:]).T).T

    ##############################################


    return dvardeta





#######################################################
#Compute Jacobian on sigma-levels (1st order)
#######################################################


def jacob_sig(var1,var2,pm,pn,z_r,z_w=None,newz=None,mask=None):

    print('jacob ,var, var2', var1.shape, var2.shape)

    var = rho2v(diffxi(var1,pm,z_r,z_w,newz,mask)) * rho2u(diffeta(var2,pn,z_r,z_w,newz,mask))\
        - rho2v(diffxi(var2,pm,z_r,z_w,newz,mask)) * rho2u(diffeta(var1,pn,z_r,z_w,newz,mask))

    print('jacob ,final', var.shape)

    return var
    




#######################################################
#Compute Jacobian on vertical-levels 
#######################################################


def jacob(var1,var2,pm,pn):

    var = rho2v(diffx(var1,pm)) * rho2u(diffy(var2,pn))\
        - rho2v(diffx(var2,pm)) * rho2u(diffy(var1,pn))

    return var   
    
    
    
    
    
    
    
    
#######################################################
#Compute Laplacien
#######################################################


def laplacien(var,pm,pn):
    
    if np.ndim(var)==3:
        
        dvar = np.zeros(var.shape)

        for iz in range(var.shape[2]):
            dvar[:,:,iz] = laplacien_2d(var[:,:,iz],pm,pn)
            
    else:
        
        dvar = laplacien_2d(var,pm,pn)

    return dvar   

##########################


def laplacien_2d(var,pm,pn):

    dvar = np.zeros(var.shape)
    
    #print 'no boundary conditions yet'
    
    dvar[1:-1,:] = (var[2:,:]- var[1:-1,:]) * 0.5 * (pm[2:,:]+pm[1:-1,:]) * pm[1:-1,:]\
                 + (var[:-2,:] - var[1:-1,:])* 0.5 * (pm[:-2,:]+pm[1:-1,:]) * pm[1:-1,:]
    
    dvar[:,1:-1] =  dvar[:,1:-1]\
                 + (var[:,2:]- var[:,1:-1]) * 0.5 * (pn[:,2:]+pn[:,1:-1]) * pn[:,1:-1]\
                 + (var[:,:-2] - var[:,1:-1])* 0.5 * (pn[:,:-2]+pn[:,1:-1]) * pn[:,1:-1]
                 
    dvar = nanbnd(dvar)

    return dvar   
    
    


#######################################################
#Compute verticale derivative on vertical levels
#######################################################

def diffz(var,depths):  

    if len(var.shape)==3:
        z_depths = copy(var[:,:,:-1]);
        for i in range(var.shape[0]): 
            for j in range(var.shape[1]):
                z_depths[i,j,:] = 0.5*(depths[1:]+depths[:-1])
                
        dvardz = vinterp((var[:,:,1:] - var[:,:,:-1])/(depths[1:]-depths[:-1]),depths,z_depths)
        
    elif len(var.shape)==2:
        z_depths = copy(var[:,:-1]);
        for i in range(var.shape[0]): 
            z_depths[i,:] = 0.5*(depths[1:]+depths[:-1])
                
        dvardz = vinterp((var[:,1:] - var[:,:-1])/(depths[1:]-depths[:-1]),depths,z_depths)     
        
    return dvardz
      

#######################################################
#Compute verticale derivative on sigma levels
#######################################################

def diffz_sig(var,z_r,z_w):  

    dvardz = vinterp((var[:,:,1:]-var[:,:,:-1])/(z_r[:,:,1:]-z_r[:,:,:-1]),z_r,z_w[:,:,1:-1],z_r)

    return dvardz  
        


#######################################################
#Compute 2nd order verticale derivative on vertical levels
#######################################################

def diffzz(var,depths):


        dvardz = np.zeros(var.shape)
        dvardz[:,:,1:-1] = ((var[:,:,2:]- var[:,:,1:-1])/(depths[2:]-depths[1:-1])\
                         -  (var[:,:,1:-1] - var[:,:,:-2])/(depths[1:-1]-depths[:-2]))\
                         / (0.5*(depths[2:]-depths[:-2]))
                         
        z_depths = copy(var[:,:,1:-1]);
        for i in range(var.shape[0]): 
            for j in range(var.shape[1]):
                z_depths[i,j,:] = depths[1:-1]
                
        dvardz[:,:,0] = vinterp(dvardz[:,:,1:-1],[depths[0]],z_depths)[:,:,0]
        dvardz[:,:,-1] = vinterp(dvardz[:,:,1:-1],[depths[-1]],z_depths)[:,:,0]
        
        return dvardz
        
        
              
#######################################################
#Compute mean
#######################################################

def nanmean(data,axis=None, *args):
    '''
    Compute the mean of an array along a given axe or not, nan are masked
    '''
    
    dataout = np.ma.filled(np.ma.masked_array(data,np.isnan(data)).mean(*args,axis=axis), fill_value=np.nan)
    if dataout.shape==(): dataout = float(dataout)

    return dataout


              
#######################################################
#Compute max
#######################################################


def nanmax_n(var2d,N,axis=1):
    '''
    Compute the mean of the N max values of an array along a given axe 
    (N=1 gives the same result than the regular max function: np.nanmax)
    '''
    
    if np.ndim(var2d)==2 and axis>0:
        a = np.ma.array(var2d, mask=False)
        if axis==1:
            amax=np.zeros(a.shape[0])
            for ivar in range(a.shape[0]):
                imax=np.zeros(N,int)
                for i in range(N):
                    try:
                        imax[i] = np.nanargmax(a[ivar,:])
                        amax[ivar] =  amax[ivar]+ a[ivar,imax[i]]/N
                        a.mask[ivar,imax[i]] = True
                    except:
                        print('not enough values')
        elif axis==0:
            amax=np.zeros(a.shape[1])
            for ivar in range(a.shape[1]):
                imax=np.zeros(N,int)
                for i in range(N):
                    try:
                        imax[i] = np.nanargmax(a[:,ivar])
                        amax[ivar] =  amax[ivar]+ a[imax[i],ivar]/N
                        a.mask[imax[i],ivar] = True      
                    except:
                        print('not enough values')

    else:
        a = np.ma.array(var2d.ravel(), mask=False)
        amax=0.; imax=np.zeros(N,int)
        for i in range(N):
            try:
                imax[i] = np.nanargmax(a)
                amax =  amax+ a[imax[i]]/N
                a.mask[imax[i]] = True      
            except:
                print('not enough values')
                
    return amax



#######################################################
#Rotate winds or u,v to lat,lon coord -> result on psi grid
#######################################################



def rotuv(simul,u,v,psi=True,**kwargs):


    if isinstance(simul, float):
        angle = simul 
    else:
        if 'coord' in  kwargs: 
            [ny1,ny2,nx1,nx2]= kwargs['coord']
        else:
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]
        
        ncfile = Dataset(simul.ncname.grd, 'r', format='NETCDF3_CLASSIC')
        if psi:
            angle = rho2psi(simul.Forder( np.array(ncfile.variables['angle'][ny1:ny2,nx1:nx2]) ))
        else:
            angle = simul.Forder( np.array(ncfile.variables['angle'][ny1:ny2,nx1:nx2] ))

    if u.shape!=v.shape:
        u=rho2v(u)
        v=rho2u(v)

    'rotate vectors by geometric angle'
    urot = u*np.cos(angle) - v*np.sin(angle)
    vrot = u*np.sin(angle) + v*np.cos(angle)

    return [urot,vrot]

#######################################################
#Rotate winds or u,v to lat,lon coord -> result on rho grid
#######################################################



def rotuv_rho(simul,u,v,**kwargs):


    if isinstance(simul, float):
        angle = simul 
    else:
        if 'coord' in  kwargs: 
            [ny1,ny2,nx1,nx2]= kwargs['coord']
        else:
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]
        
        try:
            angle = simul.angle[nx1:nx2,ny1:ny2]
        except:
            ncfile = Dataset(simul.ncname.grd, 'r', format='NETCDF3_CLASSIC')
            angle = simul.Forder( np.array(ncfile.variables['angle'][ny1:ny2,nx1:nx2]) )

    if u.shape!=v.shape:
        u=u2rho(u)
        v=v2rho(v)
    
    'rotate vectors by geometric angle'
    urot = u*np.cos(angle) - v*np.sin(angle)
    vrot = u*np.sin(angle) + v*np.cos(angle)

    return [urot,vrot]


#######################################################
#Get distance from lon-lat grid
#######################################################

def lonlat_to_m(lon,lat):
    lon = lon*2*np.pi/360.
    lat = lat*2*np.pi/360.
    if len(lat.shape)==2:
        dx = np.arccos(np.sin(lat[1:,:])*np.sin(lat[:-1,:]) + np.cos(lat[1:,:])*np.cos(lat[:-1,:])*np.cos(lon[1:,:]-lon[:-1,:]))*6371000.
        dy = np.arccos(np.sin(lat[:,1:])*np.sin(lat[:,:-1]) + np.cos(lat[:,1:])*np.cos(lat[:,:-1])*np.cos(lon[:,1:]-lon[:,:-1]))*6371000.
    else:
        dx = np.arccos(np.sin(lat[1:])*np.sin(lat[:-1]) + np.cos(lat[1:])*np.cos(lat[:-1])*np.cos(lon[1:]-lon[:-1]))*6371000.
        dy = np.arccos(np.sin(lat[1:])*np.sin(lat[:-1]) + np.cos(lat[1:])*np.cos(lat[:-1])*np.cos(lon[1:]-lon[:-1]))*6371000.        
    return dx,dy





#######################################################
#Put nan at boundaries
#######################################################

def nanbnd(var,nbp=1):


    if len(var.shape)==1:
        var[:nbp] = np.nan
        var[-nbp:] = np.nan       

    elif len(var.shape)==2:
        var[:nbp,:] = np.nan
        var[-nbp:,:] = np.nan       
        var[:,:nbp] = np.nan       
        var[:,-nbp:] = np.nan
        
    elif len(var.shape)==3:
        var[:nbp,:,:] = np.nan
        var[-nbp:,:,:] = np.nan       
        var[:,:nbp,:] = np.nan       
        var[:,-nbp:,:] = np.nan    
         
    elif len(var.shape)==4:
        var[:nbp,:,:,:] = np.nan
        var[-nbp:,:,:,:] = np.nan       
        var[:,:nbp,:,:] = np.nan       
        var[:,-nbp:,:,:] = np.nan           
        
        
    return var

    
    
#######################################################
#Put nan at boundaries
#######################################################

def zerobnd(var,nbp=1):


    if len(var.shape)==1:
        var[:nbp] = 0
        var[-nbp:] = 0       

    elif len(var.shape)==2:
        var[:nbp,:] = 0
        var[-nbp:,:] = 0       
        var[:,:nbp] = 0       
        var[:,-nbp:] = 0
        
    elif len(var.shape)==3:
        var[:nbp,:,:] = 0
        var[-nbp:,:,:] = 0       
        var[:,:nbp,:] = 0       
        var[:,-nbp:,:] = 0    
         
    elif len(var.shape)==4:
        var[:nbp,:,:,:] = 0
        var[-nbp:,:,:,:] = 0       
        var[:,:nbp,:,:] = 0       
        var[:,-nbp:,:,:] = 0           
        
        
    return var
    
    
    
    
    

#######################################################
# Find indices of a sub region
#######################################################

#from matplotlib.nxutils import points_inside_poly (nxutils not included in matplotlib since 1.3 version)
from matplotlib.path import Path

def nansub_bool(x,y,xsub,ysub):

    
    # get contour of the sub region
    polygon=np.array((np.hstack((xsub[0,:],xsub[:,-1],xsub[-1,::-1],xsub[::-1,0])),np.hstack((ysub[0,:],ysub[:,-1],ysub[-1,::-1],ysub[::-1,0])))).T  
    
    # get all points from the region
    points = np.array((x.ravel(),y.ravel())).T
    
    # test wether points are inside the sub region -> return a boolean array
    #var_bool = nxutils.points_inside_poly(points, polygon).reshape(x.shape)
    var_bool = Path(polygon).contains_points(points).reshape(x.shape)

    return var_bool  
    
    
     
    
#######################################################
#Put nan in a sub region
#######################################################


def nansub(var,x,y,xsub,ysub):

    var_bool = nansub_bool(x,y,xsub,ysub)
    
    if isinstance(var,int) or isinstance(var,float): var[var_bool] = np.nan
    else: var.data[var_bool] = np.nan
    
    return var   
    
   
   
      
      
      
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

    return i,j

    
#######################################################
#Compute Potential Vorticity of a 3-D field on psi-w grid
#######################################################
'''

Set of functions used to compute ertel potential vorticity using buoyancy (b=-g rho/rho0)

T and S on horizontal rho grids and vertical rho-grid (specified by z_r)
U and V on horizontal u- and v- grids and vertical rho-grid (specified by z_r)

PV is computed on horizontal psi-grid and vertical w-grid (specified by z_w)

Computation directly on the sigma grid (no interpolation needed)

updated: 17/06/15

'''

############################################

def get_buoy(T,S,z_r,z_w,rho0,g,userho1=False,simul=None):
    if simul!=None:
        if 'NONLIN_EOS' not in simul.cpp:
            print('using LIN_EOS')
            buoy = get_buoy_lineos(simul,T,S)
        else:
            if userho1: buoy = -g*toolsF.rho1_eos(T,S,z_r,rho0)/rho0
            else: buoy = toolsF.get_buoy(T,S,z_r,z_w,rho0)
    else:
        if userho1: buoy = -g*toolsF.rho1_eos(T,S,z_r,rho0)/rho0
        else: buoy = toolsF.get_buoy(T,S,z_r,z_w,rho0)
    return buoy

############################################

def get_buoy_lineos(simul,T,S):
    '''
    ! R0         Coefficients for linear Equation of State (EOS)
    ! T0,Tcoef  
    ! S0,Scoef          rho = R0 - Tcoef*(T-T0) + Scoef*(S-S0)

    lin_EOS_cff:  R0 [kg/m3], T0 [Celsius], S0 [PSU], TCOEF [1/Celsius], SCOEF [1/PSU]
                  30.         15.           35.       0.28d0             0.78d0

    '''
    T0 = 15.; S0=35.
    rho1 = simul.R0 - simul.Tcoef*(T - T0) + simul.Scoef*(S - S0)
    rho = rho1 + 1000. - simul.rho0
    buoy = -simul.g * rho/simul.rho0
    ################
    return buoy

############################################
def rho_grad_buoy_sig(T,S,z_r,z_w,rho0,g,pm,pn,simul=None):
    buoy = get_buoy(T,S,z_r,z_w,rho0,g,simul=simul)
    #########################################
    dbdx = diffx(buoy,pm)
    dbdy = diffy(buoy,pn)
    dbdz = (buoy[:,:,1:] - buoy[:,:,:-1])/(z_r[:,:,1:] - z_r[:,:,:-1])
    #########################################
    return [dbdx,dbdy,dbdz]

############################################
def PV(temp,salt,u,v,z_r,z_w,f,g,rho0,pm,pn,mask=None,simul=None):
    #rho on rho-rho grid      bvf on rho-w grid
    [dbdx,dbdy,dbdz] = rho_grad_buoy_sig(temp,salt,z_r,z_w,rho0,g,pm,pn,simul=simul)
    dz_r = z_r[:,:,1:]- z_r[:,:,:-1]
    dz_r[dz_r==0] = np.nan
    pv=np.zeros((z_w.shape[0]-1,z_w.shape[1]-1,z_w.shape[-1]))*np.nan
    ##########################
    #Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*db/dz
    #dudy and dvdx on psi-rho grid
    dvdx = diffx(v, rho2v(pm))
    dudy = diffy(u, rho2u(pn))
    dbdz = rho2psi(dbdz)
    #vrt on psi-rho grid
    vrt = dvdx - dudy
    #PV1 on psi-w grid
    pv[:,:,1:-1] =  (( rho2psi(f).T + 0.5*(vrt[:,:,1:] + vrt[:,:,:-1]).T).T * dbdz)
    del vrt,dbdz
    ##########################
    #'Ertel potential vorticity, term 2: (dv/dz)*(db/dx)'
    #'dvdz on psi-w grid'
    dvdz = rho2u((v[:,:,1:]-v[:,:,:-1])/(0.5*(dz_r[:,1:,:]+ dz_r[:,:-1,:])))
    #'dbdx on psi-rho grid'
    dbdx = rho2v(dbdx)
    #PV1 on psi-w grid
    pv[:,:,1:-1] = pv[:,:,1:-1] -1*dvdz*0.5*(dbdx[:,:,1:] + dbdx[:,:,:-1])
    del dbdx,dvdz
    ##########################
    #'Ertel potential vorticity, term 3: (du/dz)*(db/dy)'
    #'dudz on psi-w grid'
    dudz = rho2v((u[:,:,1:]-u[:,:,:-1])/(0.5*(dz_r[1:,:,:]+ dz_r[:-1,:,:])))
    #'dbdy on psi-rho grid'
    dbdy = rho2u(dbdy)
    #PV3 on psi-w grid
    pv[:,:,1:-1] = pv[:,:,1:-1] + dudz*0.5*(dbdy[:,:,1:] + dbdy[:,:,:-1])
    del dbdy,dudz
    return pv

############################################
def PVr(temp,salt,u,v,z_r,z_w,f,g,rho0,pm,pn,mask=None,simul=None):
    #rho on rho-rho grid      bvf on rho-psi grid
    [dbdx,dbdy,dbdz] = rho_grad_buoy_sig(temp,salt,z_r,z_w,rho0,pm,pn,simul=simul)
    dz_r = z_r[:,:,1:]- z_r[:,:,:-1]
    dz_r[dz_r==0] = np.nan
    pv=np.zeros((z_w.shape[0]-1,z_w.shape[1]-1,z_r.shape[-1]))*np.nan
    ##########################
    #Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*db/dz
    #dudy and dvdx on psi-rho grid
    dvdx = diffx(v, rho2v(pm))
    dudy = diffy(u, rho2u(pn))
    dbdz = vinterp(dbdz,z_r,z_w[:,:,1:-1],z_r)
    dbdz = rho2psi(dbdz)
    #vrt on psi-rho grid
    vrt = dvdx - dudy
    #PV1 on psi-rho grid
    pv =  (( rho2psi(f).T + vrt.T).T * dbdz)
    del vrt,dbdz
    ##########################
    #'Ertel potential vorticity, term 2: (dv/dz)*(db/dx)'
    #'dvdz on psi-rho grid'
    dvdz = rho2u((v[:,:,1:]-v[:,:,:-1])/(0.5*(dz_r[:,1:,:]+ dz_r[:,:-1,:])))
    dvdz = vinterp(dvdz, rho2psi(z_r), rho2psi(z_w[:,:,1:-1]), rho2psi(z_r))
    #'dbdx on psi-rho grid'
    dbdx = rho2v(dbdx)
    #PV1 on psi-w grid
    pv = pv -1*dvdz*dbdx
    del dbdx,dvdz
    ##########################
    #'Ertel potential vorticity, term 3: (du/dz)*(db/dy)'
    #'dudz on psi-w grid'
    dudz = rho2v((u[:,:,1:]-u[:,:,:-1])/(0.5*(dz_r[1:,:,:]+ dz_r[:-1,:,:])))
    dudz = vinterp(dudz, rho2psi(z_r), rho2psi(z_w[:,:,1:-1]), rho2psi(z_r))
    #'dbdy on psi-rho grid'
    dbdy = rho2u(dbdy)
    #PV3 on psi-w grid
    pv = pv + dudz*dbdy
    del dbdy,dudz
    ##########################
    return pv


    
#######################################################
'''
DEPRECATED : no need for vertical interpolation...

def PV_old(temp,salt,u,v,z_r,z_w,f,g,rho0,pm,pn,mask=None):

    #print 'we are using python version for PV'

    #rho on rho-rho grid      bvf on rho-w grid
    [dbdx,dbdy,dbdz] = toolsF.rho_grad(temp,salt,z_r,z_w,rho0,pm,pn)
    dz_r = z_r[:,:,1:]- z_r[:,:,:-1]
    dz_r[dz_r==0] = np.nan
    pv=np.zeros((z_w.shape[0]-1,z_w.shape[1]-1,z_w.shape[-1]))*np.nan

##########################
#Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*db/dz


    #dudy and dvdx on psi-rho grid
    dvdx = diffxi(v,rho2v(pm),rho2v(z_r),rho2v(z_w),mask=mask)
    dudy = diffeta(u,rho2u(pn),rho2u(z_r),rho2u(z_w),mask=mask)

    dbdz = rho2psi(dbdz)[:,:,1:-1]

    #vrt on psi-rho grid
    vrt = dvdx - dudy

    #PV1 on psi-w grid
    pv[:,:,1:-1] =  ((rho2psi(f).T + 0.5*(vrt[:,:,1:] + vrt[:,:,:-1]).T).T * dbdz)
    del vrt,dbdz

##########################
#'Ertel potential vorticity, term 2: (dv/dz)*(db/dx)'

    #'dvdz on psi-w grid'
    dvdz = rho2u((v[:,:,1:]-v[:,:,:-1])/(0.5*(dz_r[:,1:,:]+ dz_r[:,:-1,:])))

    #'dbdx on psi-rho grid'
    dbdx = rho2v(dbdx)

    #PV1 on psi-w grid
    pv[:,:,1:-1] = pv[:,:,1:-1] -1*dvdz*0.5*(dbdx[:,:,1:] + dbdx[:,:,:-1])
    del dbdx,dvdz

##########################
#'Ertel potential vorticity, term 3: (du/dz)*(db/dy)'

    #'dudz on psi-w grid'
    dudz = rho2v((u[:,:,1:]-u[:,:,:-1])/(0.5*(dz_r[1:,:,:]+ dz_r[:-1,:,:])))

    #'dbdy on psi-rho grid'
    dbdy = rho2u(dbdy)

    #PV3 on psi-w grid
    pv[:,:,1:-1] = pv[:,:,1:-1] + dudz*0.5*(dbdy[:,:,1:] + dbdy[:,:,:-1])
    
    del dbdy,dudz

##########################

    return pv


#######################################################
#Compute Potential Vorticity of a 3-D field on psi-rho grid
#######################################################

def PVr(temp,salt,u,v,z_r,z_w,f,g,rho0,pm,pn,mask=None):

    #print 'we are using python version for PV'

    #rho on rho-rho grid      bvf on rho-psi grid
    [dbdx,dbdy,dbdz] = toolsF.rho_grad(temp,salt,z_r,z_w,rho0,pm,pn)
    dz_r = z_r[:,:,1:]- z_r[:,:,:-1]
    dz_r[dz_r==0] = np.nan
    pv=np.zeros((z_w.shape[0]-1,z_w.shape[1]-1,z_r.shape[-1]))*np.nan


##########################
#Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*db/dz


    #dudy and dvdx on psi-rho grid
    dvdx = diffxi(v,rho2v(pm),rho2v(z_r),rho2v(z_w),mask=mask)
    dudy = diffeta(u,rho2u(pn),rho2u(z_r),rho2u(z_w),mask=mask)

    dbdz = vinterp(dbdz[:,:,1:-1],z_r,z_w[:,:,1:-1],z_r)
    dbdz = rho2psi(dbdz)

    #vrt on psi-rho grid
    vrt = dvdx - dudy

    #PV1 on psi-rho grid
    pv =  ((rho2psi(f).T + vrt.T).T * dbdz)
    del vrt,dbdz


##########################
#'Ertel potential vorticity, term 2: (dv/dz)*(db/dx)'

    #'dvdz on psi-rho grid'
    dvdz = rho2u((v[:,:,1:]-v[:,:,:-1])/(0.5*(dz_r[:,1:,:]+ dz_r[:,:-1,:])))
    dvdz = vinterp(dvdz,rho2psi(z_r),rho2psi(z_w[:,:,1:-1]),rho2psi(z_r))

    #'dbdx on psi-rho grid'
    dbdx = rho2v(dbdx)

    #PV1 on psi-w grid
    pv = pv -1*dvdz*dbdx
    del dbdx,dvdz

##########################
#'Ertel potential vorticity, term 3: (du/dz)*(db/dy)'

    #'dudz on psi-w grid'
    dudz = rho2v((u[:,:,1:]-u[:,:,:-1])/(0.5*(dz_r[1:,:,:]+ dz_r[:-1,:,:])))
    dudz = vinterp(dudz,rho2psi(z_r),rho2psi(z_w[:,:,1:-1]),rho2psi(z_r))

    #'dbdy on psi-rho grid'
    dbdy = rho2u(dbdy)

    #PV3 on psi-w grid
    pv = pv + dudz*dbdy

    del dbdy,dudz

##########################

    return pv

'''



    
#######################################################
#Compute Potential Vorticity Stetching of a 3-D field on rho-w grid
#######################################################
'''

S = -F d(Dz)/dz

'''

def stretching(temp,salt,z_r,z_w,f,g,rho0,pm,pn,mask=None):


##########################

    rho = toolsF.rho_eos(temp,salt,z_r,z_w,rho0)
    depths = np.linspace(z_w.min(),z_w.max(),5000)
    rho_z = vinterp(rho,depths,z_r,z_w)
    rho_z_mean = nanmean(nanmean(rho_z,0),0)
    del rho_z


    z_rho_mean = copy(z_r)

    for i in range(z_r.shape[0]):
      for j in range(z_r.shape[1]):
        for k in range(z_r.shape[2]):
          if not np.isnan(rho[i,j,k]):
            z_rho_mean[i,j,k] = depths[np.nanargmin(np.abs(rho[i,j,k] - rho_z_mean))]

    delta_z = z_r - z_rho_mean
    del rho,z_rho_mean

    ddelta_z = (delta_z[:,:,1:] - delta_z[:,:,:-1])/(z_r[:,:,1:]- z_r[:,:,:-1])

    return -ddelta_z

##########################





    
    
#######################################################
#Compute absolute vorticity of a 3-D field on psi grid
#######################################################

def get_absvrt(u,v,z_r,z_w,f,pm,pn,mask=None):

##########################
#Absolute vorticity,  [f + (dv/dx - du/dy)]

    vrt = get_vrt(u,v,z_r,z_w,pm,pn,mask)
    
    var =  (rho2psi(f).T + vrt.T).T 
    
    return var



#######################################################
#Compute relative vorticity of a 3-D field on psi grid
#######################################################

def get_vrt(u,v,z_r,z_w,pm,pn,mask=None):

    if len(u.shape)==3:
        #dudy and dvdx on psi grid
        dvdx = diffxi(v,rho2v(pm),rho2v(z_r),rho2v(z_w),mask)
        dudy = diffeta(u,rho2u(pn),rho2u(z_r),rho2u(z_w),mask)       
    else:      
        dvdx = diffx(v,rho2v(pm))
        dudy = diffy(u,rho2u(pn))  
        
    #vrt on psi grid
    vrt = dvdx - dudy    
    
    return vrt





#######################################################
# Vertical integration of a 3D variable between depth1 and depth2
#######################################################

def vert_int(var,z_w,depth1,depth2):

    cff2 = np.min([depth1,depth2])
    cff1 = np.max([depth1,depth2])

    Hz = z_w[:,:,1:] - z_w[:,:,:-1]
    
    cff = z_w[:,:,:-1] - cff1
    Hz[cff>0] = 0.
    Hz[np.logical_and(cff<Hz,cff>0)] = cff[np.logical_and(cff<Hz,cff>0)]
    
    cff = z_w[:,:,1:] - cff2
    Hz[cff<0] = 0.
    Hz[np.logical_and(cff<Hz,cff>0)] = cff[np.logical_and(cff<Hz,cff>0)]
    
    varint = np.nansum(Hz * var,2)
    #varint = np.nansum( var,2)    
    
    return varint
    


    
#######################################################
# Integrate function along vertical for a 3D variable
#######################################################

import scipy.integrate as integrate
 
def cumtrapz(var,z,inv=False):
    
    varint = np.zeros((var.shape[0],var.shape[1],var.shape[2]))
    
    if np.ndim(z)==0:
    
        if inv:
            varint[:,:,1:] = integrate.cumtrapz(var*z,axis=2)
        else:
            varint[:,:,:-1] = integrate.cumtrapz(var[:,:,::-1]*z[:,:,::-1],axis=2)[:,:,::-1]
            
    elif np.ndim(z)==1:
    
        for i in range(var.shape[0]):
            for j in range(var.shape[1]):
                var[i,j,:] = var[i,j,:]*z
                
        if inv:
            varint[:,:,1:] = integrate.cumtrapz(var,axis=2)
        else:
            varint[:,:,:-1] = integrate.cumtrapz(var[:,:,::-1],axis=2)[:,:,::-1]           
            
    elif np.ndim(z)==3:
    
        if inv:
            varint[:,:,1:] = integrate.cumtrapz(var*z,axis=2)
        else:
            varint[:,:,:-1] = integrate.cumtrapz(var[:,:,::-1]*z[:,:,::-1],axis=2)[:,:,::-1]
            
    return varint
    



#######################################################
# Compute circumcircle radius given a npoints*2 array containing [x,y(x)]
#######################################################

 
def radius(var,nsmooth=1):
    '''
    Compute curvature and radius of curvature
    
    http://mathworld.wolfram.com/Circle.html
    
    '''
    
    r = np.zeros(var.shape[0])*np.nan
    center = np.zeros((var.shape[0],2))*np.nan
     
    for ix in range(1,var.shape[0]-1):

        xmin1=np.max([ix-nsmooth,0]); xmax2=np.min([ix+nsmooth,var.shape[0]])
        xmin2=np.max([ix-nsmooth/2,0]); xmax1=np.min([ix+nsmooth/2,var.shape[0]])
            
        if nsmooth<20: xmin2=xmin1; xmax1=xmax2

        [x1,y1] = [nanmean(var[xmin1:xmin2,0]),nanmean(var[xmin1:xmin2,1])]
        [x2,y2] = [var[ix,0],var[ix,1]]
        [x3,y3] = [nanmean(var[xmax1:xmax2,0]),nanmean(var[xmax1:xmax2,1])]


        xs = np.array([x1, x2, x3])
        ys = np.array([y1, y2, y3])

        os = np.ones(xs.shape[0])
        ss = xs**2+ys**2

        a = np.linalg.det(np.column_stack((xs, ys, os))); 
        d = -np.linalg.det(np.column_stack((ss, ys, os))); 
        e = np.linalg.det(np.column_stack((ss, xs, os))); 
        f = -np.linalg.det(np.column_stack((ss, xs, ys))); 

        mysign = np.sign(np.arctan2(x2-x3,y2-y3) - np.arctan2(x1-x2,y1-y2) )

        r[ix] = mysign*np.sqrt((d**2+e**2)/(4*a**2)-(f/a)); 
        

        #k[ix] = 1/r; 

        x0,y0 = -d/(2*a), -e/(2*a)

        '''
        #old version (works too)
        a = np.sqrt((x2-x1)**2+(y2-y1)**2)
        b = np.sqrt((x2-x3)**2+(y2-y3)**2)
        c = np.sqrt((x3-x1)**2+(y3-y1)**2)

        mysign = np.sign(np.arccos(((y2-y1)*(y3-y2)+(x3-x2)*(x2-x1))/(a*b)))
        r[ix] = mysign*a*b*c/np.sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c)) 
        
        '''
        
        center[ix,0] = x0
        center[ix,1] = y0
        
         
    return r, center
    

    
    

#######################################################
# Compute season
#######################################################

 
def season(simul,winter='DJF'):
    day = (simul.oceantime%(360.*24*3600))/(24*3600.)
    if winter=='DJF':
        if 330 <= day < 360 or 0.0 <= day < 60: season = 'winter'
        elif 60 <= day < 150: season = 'spring'    
        elif 150 <= day < 240: season = 'summer'
        elif 240 <= day < 330: season = 'autumn'
    elif winter=='JFM':
        if 0.0 <= day < 90: season = 'winter'
        elif 90 <= day < 180: season = 'spring'    
        elif 180 <= day < 270: season = 'summer'
        elif 270 <= day < 360: season = 'autumn'
    return season



#######################################################
# Compute season
#######################################################

 
def season_int(simul,winter='DJF'):
    day = (simul.oceantime%(360.*24*3600))/(24*3600.)
    if winter=='DJF':
        if 330 <= day < 360 or 0.0 <= day < 60: season = 0
        elif 60 <= day < 150: season = 1   
        elif 150 <= day < 240: season = 2
        elif 240 <= day < 330: season = 3
    elif winter=='JFM':
        if 0.0 <= day < 90: season = 0
        elif 90 <= day < 180: season = 1    
        elif 180 <= day < 270: season = 2
        elif 270 <= day < 360: season = 3
    return season


###################################################################################
#define levels
###################################################################################

def levels(vavar, var=None , nlev=100):

    minvar=np.nanmin(var); maxvar=np.nanmax(var)
    levelsvar=np.arange(minvar,maxvar+(maxvar-minvar)/nlev,(maxvar-minvar)/nlev)

    return levelsvar



###################################################################################
#define colorabar levels using levelsvar
###################################################################################

def clabels(levelsvar,nblab=4.,sym=0):

    tot=levelsvar.max()-levelsvar.min()
    test=0; i=0
    eps = tot*1e-6
    while test==0:
        test=round(tot,i)
        i=i+1

    dlab=round(tot/nblab,i)
    labmin=round(levelsvar.min(),i)
    labmax=round(levelsvar.max(),i)

    labels=np.arange(labmin,labmax+dlab,dlab)
    
    # Impose that 0 is among labels if sym=1
    if sym==1:
        labels = np.arange(-2*dlab,3*dlab,dlab)

 
    return labels[np.logical_and(labels<=levelsvar.max()+eps,labels>=levelsvar.min()-eps)]



######################################################
def rhop(T,S):
######################################################
    """Density of seawater at zero pressure"""

    # --- Define constants ---
    a0 = 999.842594
    a1 =   6.793952e-2
    a2 =  -9.095290e-3
    a3 =   1.001685e-4
    a4 =  -1.120083e-6
    a5 =   6.536332e-9
    b0 =   8.24493e-1
    b1 =  -4.0899e-3
    b2 =   7.6438e-5
    b3 =  -8.2467e-7
    b4 =   5.3875e-9
    c0 =  -5.72466e-3
    c1 =   1.0227e-4
    c2 =  -1.6546e-6
    d0 =   4.8314e-4
    # --- Computations ---
    # Density of pure water
    SMOW = a0 + (a1 + (a2 + (a3 + (a4 + a5*T)*T)*T)*T)*T
    # More temperature polynomials
    RB = b0 + (b1 + (b2 + (b3 + b4*T)*T)*T)*T
    RC = c0 + (c1 + c2*T)*T
    return SMOW + RB*S + RC*(S**1.5) + d0*S*S




  

