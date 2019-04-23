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
        topo = np.asfortranarray(simul.topo[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])

    else: 
        coord = simul.coord
        topo = simul.topo
        

    if hasattr(simul, 'zeta'): 
        zeta=simul.zeta[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1]
    else: 
        zeta = va.var('zeta',simul,depths=[0],coord=coord[0:4]).data

    (z_r,z_w) = toolsF.zlevs(topo, zeta, simul.hc, simul.Cs_r, simul.Cs_w)

    if 'sub' in  kwargs: 
        z_r = np.asfortranarray(z_r[:,:,simul.coord[4]-1])
        z_w = np.asfortranarray(z_w[:,:,np.arange(np.min(simul.coord[4])-1,np.max(simul.coord[4])+1)])
        
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

    else: 
        coord = simul.coord
        topo = simul.topo

    if hasattr(simul, 'zeta'): zeta=simul.zeta[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1]
    else: zeta = va.var('zeta',simul,depths=[0],coord=coord[0:4]).data

    (z_r) = toolsF.zlev(topo, zeta, simul.hc, simul.Cs_r, simul.Cs_w)

    return z_r


#################################################
# rho_eos (from rho_eos.F in romsucla)
#################################################


def rho_eos(T,S,z_r,z_w,rho0):

    (rho) = toolsF.rho_eos(T,S,z_r,z_w,rho0)

    return rho

#################################################
# rho_eos (from rho_eos.F in romsucla)
#################################################


def rho1_eos(T,S,z_r,z_w,rho0):

    (rho1) = toolsF.rho1_eos(T,S,z_r,rho0)

    return rho1


#################################################
# rho_grad (from rho_eos.F and prsgrd.F in romsucla)
#################################################


def rho_grad(T,S,z_r,z_w,rho0,pm,pn):

    (drdz,drdx,drdy) = toolsF.rho_grad(T,S,z_r,z_w,rho0)

    return [drdz,drdx,drdy]


#######################################################
#interpolate a 3D variable on horizontal levels of constant depths (FORTRAN version, much faster)
#######################################################


def vinterp(var, depths, z_r, z_w=None, mask=None,imin=0,jmin=0,kmin=1, floattype=np.float64,interp_sfc=1,interp_bot=0,below=None,**kwargs):


    if mask==None:  mask = np.ones((z_r.shape[0],z_r.shape[1]), order='F', dtype=floattype); mask[z_r[:,:,-1]==0] = 0

    if z_w==None: 
        print('no z_w specified')
        z_w=np.zeros((z_r.shape[0],z_r.shape[1],z_r.shape[2]+1), order='F')
        z_w[:,:,1:-1] = 0.5*(z_r[:,:,1:] + z_r[:,:,:-1])
        z_w[:,:,0] = z_r[:,:,0] - (z_r[:,:,1]-z_r[:,:,0])
        z_w[:,:,-1] = z_r[:,:,-1] + (z_r[:,:,-1]-z_r[:,:,-2])
        
    if np.rank(depths)==1: newz = np.asfortranarray(np.zeros((z_r.shape[0],z_r.shape[1],len(depths))) + depths, dtype=floattype)
    else: newz = depths

    if interp_bot==1:
        print("data will be interpolated below ground")
        below=1000.
    	vnew=toolsF.sigma_to_z_intr_bot(z_r, z_w,mask,var,newz,below,imin,jmin,kmin,9999.)
    elif interp_sfc==1:
        print("no interpolation below ground")
        print(z_r.shape, z_w.shape,mask.shape,var.shape,newz.shape)
    	vnew=toolsF.sigma_to_z_intr_sfc(z_r, z_w,mask,var,newz,imin,jmin,kmin,9999.)
    else:
        print("no interpolation below ground")
        vnew=toolsF.sigma_to_z_intr(z_r, z_w,mask,var,newz,imin,jmin,kmin,9999.)   	

    
    vnew[np.abs(vnew)==9999.]=np.nan

    return vnew










#######################################################
#Transfert a field at psi points to rho points
#######################################################

def psi2rho(var_psi):

    if np.rank(var_psi)<3:
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

    if np.rank(var_rho)<3:
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

    if np.rank(var_rho)==1:
        var_u = 0.5*(var_rho[1:]+var_rho[:-1])
    elif np.rank(var_rho)==2:       
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

    if np.rank(var_rho)==1:
        var_v = 0.5*(var_rho[1:]+var_rho[:-1])
    elif np.rank(var_rho)==2:
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











#######################################################
#Transfert a 3-D field from verical w points to vertical rho-points
#######################################################

def w2rho(var_w):


    [M,L,N]=var_w.shape[2]
    
    var_rho = np.zeros((M,L,N-1))
    
    for iz in range(1,N-1):
        var_rho[:,:,iz]  = 0.5625*(var_w[:,:,iz+1] + var_w[:,:,iz]) -0.0625*(var_w[:,:,iz+2] + var_w[:,:,iz-1])
    
    var_rho[:,:,0]  = -0.125*var_w[:,:,2] + 0.75*var_w[:,:,1] +0.375*var_w[:,:,0] 
    var_rho[:,:,N-1]  = -0.125*var_w[:,:,N-2] + 0.75*var_w[:,:,N-1] +0.375*var_w[:,:,N] 
    

    return var_rho











'''
####################################################################################################################################
#Load variables
###################################################################################


    def load(self,varname,ncfile,simul,**kwargs):

        [ny1,ny2,nx1,nx2,depths] = self.coord

        if 'coord' in  kwargs: [ny1,ny2,nx1,nx2] = kwargs['coord'][0:4]
        if 'depths' in  kwargs: depths = kwargs['depths']

        [imin,jmin,kmin] = self.dico.get(varname)[2]; depth = np.array(depths)-1
        if len(depth)==1: depth = depth[0]


        try:
            data = np.squeeze(simul.Forder(ncfile.variables[varname][simul.infiletime,depth,ny1:ny2-jmin,nx1:nx2-imin]))
        except:
            data = np.squeeze(simul.Forder(ncfile.variables[varname][simul.infiletime,ny1:ny2-jmin,nx1:nx2-imin]))


        return data
# get_diagsPV 
#################################################


def get_pv_sol1(T,S,u,v,z_r,z_w,rho0,pm,pn,f):

    #print toolsF.get_diagspv.__doc__

    (pv) = toolsF.get_diagspv_sol1(T,S,u,v,z_r,z_w,rho0,pm,pn,f)

    return pv




#################################################
# get_diagsPV 
#################################################


def get_pv_sol2(T,S,u,v,z_r,z_w,rho0,pm,pn,f):

    #print toolsF.get_diagspv.__doc__

    (pv) = toolsF.get_diagspv_sol2(T,S,u,v,z_r,z_w,rho0,pm,pn,f)
    pv[pv==-9999.] = np.nan

    return pv

'''

#################################################
# get_J1
#################################################


def get_j1_sol1(stflx,ssflx,u,v,z_r,z_w,rho0,pm,pn,hbls,f):

    (J1) = toolsF.get_j1_sol1(stflx,ssflx,u,v,z_r,z_w,rho0,pm,pn,hbls,f)

    return J1

#################################################
# get_J1
#################################################


def get_j1_sol2(stflx,ssflx,u,v,z_r,z_w,rho0,pm,pn,hbls,f):

    (J1) = toolsF.get_j1_sol2(stflx,ssflx,u,v,z_r,z_w,rho0,pm,pn,hbls,f)

    return J1

#################################################
# get_J2
#################################################


def get_j2_sol1(T,S,u,v,z_r,z_w,rho0,pm,pn,hbls):

    (J2) = toolsF.get_j2_sol1(T,S,u,v,z_r,z_w,rho0,pm,pn,hbls)

    return J2

#################################################


def get_j2_sol2(T,S,u,v,z_r,z_w,rho0,pm,pn,hbls):

    (J2) = toolsF.get_j2_sol2(T,S,u,v,z_r,z_w,rho0,pm,pn,hbls)

    return J2



#################################################
# get_Jbot
#################################################


def get_jbot_sol1(T,S,u,v,z_r,z_w,rho0,pm,pn,hbbls,rdrg):


    (Jbot) = toolsF.get_jbot_sol1(T,S,u,v,z_r,z_w,rho0,pm,pn,hbbls,rdrg)

    return Jbot


#################################################


def get_jbot_sol2(T,S,u,v,z_r,z_w,rho0,pm,pn,hbbls,rdrg):


    (Jbot) = toolsF.get_jbot_sol2(T,S,u,v,z_r,z_w,rho0,pm,pn,hbbls,rdrg)

    return Jbot


'''

#################################################
# get bottom drag
#################################################


def get_bottom_drag(u,v,Hz,rdrg):

    (ubot,vbot) = toolsF.get_bot(u,v,Hz,rdrg)

    return [ubot,vbot]



#################################################
# get bottom pressure torque
#################################################


def get_bpt(T,S, z_r,z_w,rho0,pm,pn):


    (bpt) = toolsF.get_bpt(T,S, z_r,z_w,rho0,pm,pn)


    #joe = np.zeros((T.shape[0]-1,T.shape[1]-1))
    #joe[1:-1,1:-1] = np.asfortranarray(bpt[1:-1,1:-1])

    return bpt


#################################################
# get planetary term from vorticity balance
#################################################


def get_vortplantot(u,v,H,pm,pn,f):

    (vrtp) = toolsF.get_vortplantot(u,v,H,pm,pn,f)

    return vrtp


#################################################
# get planetary term from vorticity balance
#################################################


def get_vortplanet(u,v,H,pm,pn,f):

    (vrtp) = toolsF.get_vortplanet(u,v,H,pm,pn,f)

    return vrtp

#################################################
# get planetary stretching term from vorticity balance
#################################################


def get_vortstretch(u,v,H,pm,pn,f):

    (vrts) = toolsF.get_vortstretch(u,v,H,pm,pn,f)

    return vrts

#################################################
# get NL advective term from vort. balance equ.
#################################################


#def get_vortadv(u,v, z_r,z_w,pm,pn):

    #NOT READY YET
    #(vrta) = toolsF.get_vortadv(u,v, z_r,z_w,pm,pn)


    #return vrta

'''













#################################################
# Rotationnel
#################################################


def rot(u,v,pm,pn):
    
    if u.shape == v.shape:
        u=rho2u(u); v=rho2v(v)
    
    (rot) = toolsF.get_rot(u,v,pm,pn)

    return rot



#################################################
# Gradient (amplitude)
#################################################


def grad(psi,pm,pn):

    (grad) = toolsF.get_grad(psi,pm,pn)

    return grad















#######################################################
#x-derivative from rho-grid to u-grid
#######################################################

def diffx(var,pm,dn=1):

    if np.rank(var)<3:
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

    if (np.rank(pm)==2) and (var.shape[0]==pm.shape[0]): 
        dvardx = (var[dn:,:]-var[:-dn,:])*0.5*(pm[dn:,:]+pm[:-dn,:])/dn
    else: 
        dvardx = (var[dn:,:]-var[:-dn,:])*pm/dn

    return dvardx




#######################################################
#y-derivative from rho-grid to v-grid
#######################################################

def diffy(var,pn,dn=1):

    if np.rank(var)<3: dvardy = diffy_2d(var,pn,dn)
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

    if (np.rank(pn)==2) and (var.shape[1]==pn.shape[1]):
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


    if z_r.shape[2]<=2:
        dvardxi = diffxi_2d(var,pm,z_r,z_w,newz,mask)
    else:
        dvardxi = diffxi_3d(var,pm,z_r,z_w,newz,mask)

    ##############################################

    return dvardxi

#######################################################
#######################################################


def diffxi_3d(var,pm,z_r,z_w=None,newz=None,mask=None):


    if newz==None: newz = 0.5*(z_r[1:,:,:] + z_r[:-1,:,:])
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


def diffxi_2d(var,pm,z_r,z_w=None,newz=None,mask=None):

    dvardxi = np.zeros((z_r.shape[0]-1,z_r.shape[1]))

    ##############################################

    if newz==None: newz = 0.5*(z_r[:-1,:,0] + z_r[1:,:,0])
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
#Compute horizontal derivatives on sigma-levels (1st order)
#######################################################

'''
var on rho-rho grid
dvardxi on psi-rho grid
'''

def diffeta(var,pn,z_r,z_w=None,newz=None,mask=None):


    if z_r.shape[2]<=2:
        dvardeta = diffeta_2d(var,pn,z_r,z_w,newz,mask)
    else:
        dvardeta = diffeta_3d(var,pn,z_r,z_w,newz,mask)

    ##############################################

    return dvardeta


#######################################################
#######################################################


def diffeta_3d(var,pn,z_r,z_w=None,newz=None,mask=None):


    if newz==None: newz = 0.5*(z_r[:,:-1,:] + z_r[:,1:,:])
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



def diffeta_2d(var,pn,z_r,z_w=None,newz=None,mask=None):

    dvardeta = np.zeros((z_r.shape[0],z_r.shape[1]-1))

    ##############################################

    if newz==None: newz = 0.5*(z_r[:,:-1,0] + z_r[:,1:,0])
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
#Compute Jacobian on sigma-levels (1st order)
#######################################################


def jacob_sig(var1,var2,pm,pn,z_r,z_w=None,newz=None,mask=None):

    print('jacob ,var, var2', var1.shape, var2.shape)

    var = rho2v(diffxi(var1,pm,z_r,z_w,newz,mask)) * rho2u(diffeta(var2,pn,z_r,z_w,newz,mask))\
        - rho2v(diffxi(var2,pm,z_r,z_w,newz,mask)) * rho2u(diffeta(var1,pn,z_r,z_w,newz,mask))

    print('jacob ,final', var.shape)

    return var
    

 #######################################################
#Compute Jacobian on sigma-levels (1st order)
#######################################################


def jacob(var1,var2,pm,pn):

    var = rho2v(diffx(var1,pm)) * rho2u(diffy(var2,pn))\
        - rho2v(diffx(var2,pm)) * rho2u(diffy(var1,pn))

    return var   
    
    
    
   
#######################################################
#Compute mean
#######################################################



def nanmean(data,axis=None, *args):

    dataout = np.ma.filled(np.ma.masked_array(data,np.isnan(data)).mean(*args,axis=axis), fill_value=np.nan)
    if dataout.shape==(): dataout = float(dataout)

    return dataout






    
#######################################################
#Strain (direction + amplitude)
#######################################################


def straindir(u,v,pm=1,pn=1):
    
    if u.shape[0]==v.shape[0]: 
        ux=u2rho( diffx(u,pm))
        vy=v2rho( diffy(v,pn))
        uy=v2rho( diffy(u,pn))
        vx=u2rho( diffx(v,pm))   
        
    else:    
        if pm==1:
            ux=u2rho( u2rho( diffx(u,pm)))
            vy=v2rho( v2rho( diffy(v,pn)))
            uy=psi2rho( diffy(u,pn))
            vx=psi2rho( diffx(v,pm))           
        else:
            ux = np.zeros(pm.shape)*np.nan
            vy,uy,vx = copy(ux),copy(ux),copy(ux)
            
            ux[1:-1,:]=diffx(u,rho2u(pm))
            vy[:,1:-1]=diffy(v,rho2v(pn))
            uy=psi2rho( diffy(u,rho2u(pn)))
            vx=psi2rho( diffx(v,rho2v(pm)))
        
    s1 = ux-vy; s2= vx+uy;
    thetas = np.arctan(s2/s1)/2; thetaps = np.arctan(-1*s1/s2)/2;
    #see if division by 0
    eps = 1e-15; 
    thetas[np.abs(s1)<eps] = np.sign(s2[np.abs(s1)<eps])*np.pi/4
    #check if s1'>0 (s1<0 means that you are on the perpendicular axis)
    s1bis = s1 * np.cos(2*thetas) + s2*np.sin(2*thetas)
    thetas[s1bis<0] = thetas[s1bis<0]+np.pi/2
    s1bis = s1 * np.cos(2*thetas) + s2*np.sin(2*thetas)
    return thetas,s1bis
    

    
    
#######################################################
#Shear (direction + amplitude)
#######################################################


def sheardir(u,v,pm=1,pn=1):
    
    if u.shape[0]==v.shape[0]: 
        ux=u2rho( diffx(u,pm))
        vy=v2rho( diffy(v,pn))
        uy=v2rho( diffy(u,pn))
        vx=u2rho( diffx(v,pm))   
        
    else:    
        if pm==1:
            ux=u2rho( u2rho( diffx(u,pm)))
            vy=v2rho( v2rho( diffy(v,pn)))
            uy=psi2rho( diffy(u,pn))
            vx=psi2rho( diffx(v,pm))           
        else:   
            ux = np.zeros(pm.shape)*np.nan
            vy,uy,vx = copy(ux),copy(ux),copy(ux)
            
            ux[1:-1,:]=diffx(u,rho2u(pm))
            vy[:,1:-1]=diffy(v,rho2v(pn))
            uy=psi2rho( diffy(u,rho2u(pn)))
            vx=psi2rho( diffx(v,rho2v(pm)))
        
    s1 = ux-vy; s2= vx+uy; div=ux+vy
    #thetas = np.arctan(s2/s1)/2; 
    thetas = np.arctan(-1*s1/s2)/2;
    #see if division by 0
    eps = 1e-15; 
    thetas[np.abs(s2)<eps] = 0.
    #check if s2'>0 (s2'>0 means that you are on the perpendicular axis)
    s2bis = -s1 * np.sin(2*thetas) + s2*np.cos(2*thetas)
    thetas[s2bis>0] = thetas[s2bis>0]+np.pi/2
    s2bis = -s1 * np.sin(2*thetas) + s2*np.cos(2*thetas)


    return thetas,s2bis
    
    
    
    


#######################################################
#Rotate winds or u,v to lat,lon coord -> result on psi grid
#######################################################



def rotuv(simul,u,v,**kwargs):


    if isinstance(simul, float):
        angle = simul 
    else:
        if 'coord' in  kwargs: 
            [ny1,ny2,nx1,nx2]= kwargs['coord']
        else:
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]
        
        ncfile = Dataset(simul.ncname.grd, 'r', format='NETCDF3_CLASSIC')
        angle = rho2psi(simul.Forder( np.array(ncfile.variables['angle'][ny1:ny2,nx1:nx2]) ))

    if u.shape!=v.shape:
        u=rho2v(u)
        v=rho2u(v)

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

    dx = np.arccos(np.sin(lat[1:,:])*np.sin(lat[:-1,:]) + np.cos(lat[1:,:])*np.cos(lat[:-1,:])*np.cos(lon[1:,:]-lon[:-1,:]))*6371000.

    dy = np.arccos(np.sin(lat[:,1:])*np.sin(lat[:,:-1]) + np.cos(lat[:,1:])*np.cos(lat[:,:-1])*np.cos(lon[:,1:]-lon[:,:-1]))*6371000.

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
#Put nan in a sub region
#######################################################

import matplotlib.nxutils as nxutils 

def nansub(var,x,y,xsub,ysub):

    var_bool = nansub_bool(x,y,xsub,ysub)
    
    if isinstance(var,int) or isinstance(var,float): var[var_bool] = np.nan
    else: var.data[var_bool] = np.nan
    
    return var   
    
    

#######################################################


def nansub_bool(x,y,xsub,ysub):

    #polygon=np.array([(xsub[0,0],ysub[0,0]),(xsub[0,-1],ysub[0,-1]),(xsub[-1,-1],ysub[-1,-1]),(xsub[-1,0],ysub[-1,0])]) 
    
    # get contour of the sub region
    polygon=np.array((np.hstack((xsub[0,:],xsub[:,-1],xsub[-1,::-1],xsub[::-1,0])),np.hstack((ysub[0,:],ysub[:,-1],ysub[-1,::-1],ysub[::-1,0])))).T  
    
    # get all points from the region
    points = np.array((x.ravel(),y.ravel())).T
    
    # test wether points are inside the sub region -> return a boolean array
    var_bool = nxutils.points_inside_poly(points, polygon).reshape(x.shape)
    
    return var_bool  
    
    
 
   
   
   
   
    
    
    
#######################################################
#Compute Potential Vorticity of a 3-D field on psi-w grid
#######################################################
'''

Compute ertel potential vorticity using buoyancy (b=-g rho/rho0)

T and S on horizontal rho grids and vertical rho-grid (specified by z_r)
U and V on horizontal u- and v- grids and vertical rho-grid (specified by z_r)

PV is computed on horizontal psi-grid and vertical w-grid (specified by z_w)

'''

def PV(temp,salt,u,v,z_r,z_w,f,g,rho0,pm,pn,mask=None):

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
#Compute Potential Vorticity of a 3-D field on psi-w grid
#######################################################
'''

Compute ertel potential vorticity using buoyancy (b=-g rho/rho0)

T and S on horizontal rho grids and vertical rho-grid (specified by z_r)
U and V on horizontal u- and v- grids and vertical rho-grid (specified by z_r)

PV is computed on horizontal psi-grid and vertical w-grid (specified by z_w)

'''

def PV_terms(temp,salt,u,v,z_r,z_w,f,g,rho0,pm,pn,mask=None):

    #print 'we are using python version for PV'

    #rho on rho-rho grid      bvf on rho-w grid
    [dbdx,dbdy,dbdz] = toolsF.rho_grad(temp,salt,z_r,z_w,rho0,pm,pn)
    dz_r = z_r[:,:,1:]- z_r[:,:,:-1]
    dz_r[dz_r==0] = np.nan

    pv1=np.zeros((z_w.shape[0]-1,z_w.shape[1]-1,z_w.shape[-1]))*np.nan
    pv2=np.zeros((z_w.shape[0]-1,z_w.shape[1]-1,z_w.shape[-1]))*np.nan
    pv3=np.zeros((z_w.shape[0]-1,z_w.shape[1]-1,z_w.shape[-1]))*np.nan

    
##########################
#Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*db/dz


    #dudy and dvdx on psi-rho grid
    dvdx = diffxi(v,rho2v(pm),rho2v(z_r),rho2v(z_w),mask=mask)
    dudy = diffeta(u,rho2u(pn),rho2u(z_r),rho2u(z_w),mask=mask)

    dbdz = rho2psi(dbdz)[:,:,1:-1]

    #vrt on psi-rho grid
    vrt = dvdx - dudy

    #PV1 on psi-w grid
    pv1[:,:,1:-1] =  ((rho2psi(f).T + 0.5*(vrt[:,:,1:] + vrt[:,:,:-1]).T).T * dbdz)
    del vrt,dbdz

##########################
#'Ertel potential vorticity, term 2: (dv/dz)*(db/dx)'

    #'dvdz on psi-w grid'
    dvdz = rho2u((v[:,:,1:]-v[:,:,:-1])/(0.5*(dz_r[:,1:,:]+ dz_r[:,:-1,:])))

    #'dbdx on psi-rho grid'
    dbdx = rho2v(dbdx)

    #PV1 on psi-w grid
    pv2[:,:,1:-1] =  -1*dvdz*0.5*(dbdx[:,:,1:] + dbdx[:,:,:-1])
    del dbdx,dvdz

##########################
#'Ertel potential vorticity, term 3: (du/dz)*(db/dy)'

    #'dudz on psi-w grid'
    dudz = rho2v((u[:,:,1:]-u[:,:,:-1])/(0.5*(dz_r[1:,:,:]+ dz_r[:-1,:,:])))

    #'dbdy on psi-rho grid'
    dbdy = rho2u(dbdy)

    #PV3 on psi-w grid
    pv3[:,:,1:-1] =  dudz*0.5*(dbdy[:,:,1:] + dbdy[:,:,:-1])
    
    del dbdy,dudz

##########################

    return [pv1,pv2,pv3]
    
    
    
    
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
    
    varint = nansum(Hz * var,2)

    print(z_w[10,10,:])
    print(Hz[10,10,:])   
    
    
    return varint
    
#######################################################
# Compute solution of TTW equation on sigma levels
# 
#######################################################


import scipy.integrate as integrate
 
def cumtrapz(var,z,inv=False):
    
    varint = np.zeros((var.shape[0],var.shape[1],var.shape[2]))
    
    if np.rank(z)==0:
    
        if inv:
            varint[:,:,1:] = integrate.cumtrapz(var*z,axis=2)
        else:
            varint[:,:,:-1] = integrate.cumtrapz(var[:,:,::-1]*z[:,:,::-1],axis=2)[:,:,::-1]
            
    elif np.rank(z)==1:
    
        for i in range(var.shape[0]):
            for j in range(var.shape[1]):
                var[i,j,:] = var[i,j,:]*z
                
        if inv:
            varint[:,:,1:] = integrate.cumtrapz(var,axis=2)
        else:
            varint[:,:,:-1] = integrate.cumtrapz(var[:,:,::-1],axis=2)[:,:,::-1]           
            
    elif np.rank(z)==3:
    
        if inv:
            varint[:,:,1:] = integrate.cumtrapz(var*z,axis=2)
        else:
            varint[:,:,:-1] = integrate.cumtrapz(var[:,:,::-1]*z[:,:,::-1],axis=2)[:,:,::-1]
            
    return varint
    
    
#######################################################
#Compute tendency given u,v,buoy
#######################################################

def get_tendency(u,v,buoy,pm,pn):
    
    if u.shape==v.shape:
        
        vx = u2rho(diffx(v ,pm))
        uy = v2rho(diffy(u ,pn))
        ux = u2rho(diffx(u ,pm))
        vy = v2rho(diffy(v ,pn))

    else:

        vx = psi2rho(diffx(v ,rho2v(pm)))
        uy = psi2rho(diffy(u ,rho2u(pn)))
        
        if len(u.shape)==3:
            vy = np.zeros((pm.shape[0],pm.shape[1],u.shape[2]))*np.nan
            vy[:,1:-1,:] = diffy(v ,rho2v(pn))
            ux = np.zeros((pm.shape[0],pm.shape[1],u.shape[2]))*np.nan
            ux[1:-1,:,:] = diffx(u ,rho2u(pm))
        else:
            vy = np.zeros(pm.shape)*np.nan
            vy[:,1:-1] = diffy(v ,rho2v(pn))
            ux = np.zeros(pm.shape)*np.nan
            ux[1:-1,:] = diffx(u ,rho2u(pm))
    
    ##############
    
    bx = u2rho(diffx(buoy,pm))
    by = v2rho(diffy(buoy,pn))
    
    tend = -1*(bx * ux * bx + by * uy * bx + bx * vx * by + by * vy * by)
    
    return tend
    
#######################################################
# Compute divergent part of the flow 
# by solving Poisson equation for velocity potential
#######################################################
'''

Note that we are using uniform grid spacing, which is fine for small
scale grids but not for very large grids such as Pacific.

We are using Dirichlet boundary conditions

'''
from pyamg import *
from pyamg.gallery import *
from scipy import *
from scipy.linalg import *  


def div2uvs(u,v,pm,pn):
    

    if len(u.shape)>2:
        
        udiv = np.zeros(u.shape)*np.nan
        vdiv = np.zeros(v.shape)*np.nan
        
        for iz in range(u.shape[2]):
            udiv[:,:,iz],vdiv[:,:,iz] = div2uv(u[:,:,iz],v[:,:,iz],pm,pn)
            
    else:
        
        udiv,vdiv = div2uv(u,v,pm,pn)

    return udiv,vdiv

##################


def div2uv(u,v,pm,pn):


    pm = np.ones(pm.shape)*np.mean(pm)
    pn = np.ones(pm.shape)*np.mean(pn)
  
  
    # compute div
    div = np.zeros(pm.shape)
    div[1:-1,:] = div[1:-1,:] + diffx(u,rho2u(pm)) 
    div[:,1:-1] = div[:,1:-1] + diffy(v,rho2v(pn))
    div[isnan(div)] =0
    
    # solve poisson
    A = poisson(div.shape, format='csr')     # 2D Poisson problem 
    ml =ruge_stuben_solver(A)                # construct the multigrid hierarchy
    print(ml)                                 # print hierarchy information
    b = -1*div.flatten()*1/np.mean(pm)**2    # right hand side
    x = ml.solve(b, tol=1e-10)               # solve Ax=b to a tolerance of 1e-8
    print("residual norm is", norm(b - A*x))  # compute norm of residual vector
    
    udiv = diffx(x.reshape(div.shape),pm)
    vdiv = diffy(x.reshape(div.shape),pn)
    
    return udiv,vdiv
    

    
#######################################################
# Compute solution of omega equation
# 
#######################################################
'''


'''
#from pyamg import *
#from pyamg.gallery import *
from scipy import *
from scipy.linalg import *  
from scipy.sparse import *
from scipy.sparse.linalg import spsolve
    
def solve_omega(buoy,pm,pn,f,N2,depths,ur=None,vr=None,nh=1,forcing=0.,mixrotuv=True):
    
    
    #######################################################
    #Create forcing (Q vector divergence)
    #######################################################
    rotuv=True; 
    if ur==None: rotuv=False

    new = np.zeros(buoy.shape)

    nz=len(depths); [nx,ny] = pm.shape; ndim = (nz-1)*ny*nx;
    print('number of points is ',nx,ny,nz)
    dz = depths[1:]-depths[:-1]

    #get gradients
    bx,by = copy(new),copy(new)
    bx[1:-1,:,:]= diffx(buoy,pm,2); by[:,1:-1,:] = diffy(buoy,pn,2)  
    
    #if periodicity in y (used for jet example)
    #by[:,0,:] = (buoy[:,1,:]-buoy[:,-1,:])*0.5*(pn[:,1]+pn[:,-1])/2
    #by[:,-1,:] = (buoy[:,0,:]-buoy[:,-2,:])*0.5*(pn[:,1]+pn[:,-1])/2   
    
    ux,uy,uz = copy(new),copy(new),copy(new)
    vx,vy,vz = copy(new),copy(new),copy(new)

    if rotuv:
        u = ur; v=vr
        ux[1:-1,:,:] = diffx(u,pm,2); uy[:,1:-1,:] = diffy(u,pn,2)   
        vx[1:-1,:,:] = diffx(v,pm,2); vy[:,1:-1,:] = diffy(v,pn,2)       
        
        if mixrotuv:
            print('using velocity and buoyancy field to define Q vector')
            uz = -(by.T/f.T).T; vz =  (bx.T/f.T).T;
        
        else:
            #using only non-divergent velocity field
            print('using only velocity field to define Q vector')
            z_depths = copy(uz[:,:,:-1]);
            for i in range(nx): 
                for j in range(ny):
                    z_depths[i,j,:] = 0.5*(depths[1:]+depths[:-1])
            uz = vinterp((u[:,:,1:] - u[:,:,:-1])/dz,depths,z_depths)
            vz = vinterp((v[:,:,1:] - v[:,:,:-1])/dz,depths,z_depths)
        
    else:
        #or thermal wind      
        print('using buoyancy gradient to define Q vector')
        uz = -(by.T/f.T).T; vz =  (bx.T/f.T).T;
        #u = np.cumsum(uz,2)*dz;
        #v = np.cumsum(vz,2)*dz;
        u = cumtrapz(uz,depths)
        v = cumtrapz(vz,depths)
        ux[1:-1,:,:] = diffx(u,pm,2); uy[:,1:-1,:] = diffy(u,pn,2)   
        vx[1:-1,:,:] = diffx(v,pm,2); vy[:,1:-1,:] = diffy(v,pn,2)
        
        #if periodicity in y (used for jet example)
        #uy[:,0,:] = (u[:,1,:]-u[:,-1,:])*0.5*(pn[:,1]+pn[:,-1])/2
        #uy[:,-1,:] = (u[:,0,:]-u[:,-2,:])*0.5*(pn[:,1]+pn[:,-1])/2      
        #vy[:,0,:] = (v[:,1,:]-v[:,-1,:])*0.5*(pn[:,1]+pn[:,-1])/2
        #vy[:,-1,:] = (v[:,0,:]-v[:,-2,:])*0.5*(pn[:,1]+pn[:,-1])/2
        
    #Components of Q vector = (Qx,Qy)
    Qx = 2*(f.T*(vx*uz + vy*vz).T).T;
    Qy =-2*(f.T*(ux*uz + uy*vz).T).T;

    Qxx, Qyy =  copy(new),copy(new)
    Qxx[1:-1,:,:]= diffx(Qx,pm,2); Qyy[:,1:-1,:] = diffy(Qy,pn,2)
    #if periodicity in y (used for jet example)    
    #Qyy[:,0,:] = (Qy[:,1,:]-Qy[:,-1,:])*0.5*(pn[:,1]+pn[:,-1])/2
    #Qyy[:,-1,:] = (Qy[:,0,:]-Qy[:,-2,:])*0.5*(pn[:,1]+pn[:,-1])/2
    
    print('N2.shape', N2.shape)
    divQ = (Qxx + Qyy)/N2
    
    # smooothing...
    if nh>1: 
        divQ = sm.moy_h(divQ,nh)
        f=  sm.moy_h(f,nh); 
        pm = sm.moy_h(pm,nh)/nh
        pn = sm.moy_h(pn,nh)/nh
        if len(N2.shape)==3: N2=sm.moy_h(N2,nh);
        new = np.zeros(divQ.shape)
        [nx,ny] = pm.shape; ndim = (nz-1)*ny*nx;
        print('nh is', nh)
        print('number of points is now ',nx,ny,nz)


    #######################################################
    # reorder forcings from (i,j,k) to vector
    R = np.zeros(ndim);
    for i in range(nx): 
        for j in range(ny): 
            for k in range(nz-1):
                idx = i*ny*(nz-1) + k*ny + j;
                R[idx] = divQ[i,j,k]

    R[np.isnan(R)]=0

    #######################################################
    #Create matrix A
    #######################################################    
    
    print('creating matrix A')
    A = omega_matrix(pm,pn,depths,f,N2)

    
    #######################################################
    #Solve matrix A
    #######################################################   

    A = A.tocsr()
    
    #print 'solving equation'
    
    tstart = tm.time() 

    
    # Method 1 
    ml =ruge_stuben_solver(A)                # construct the multigrid hierarchy
    print(ml)                                 # print hierarchy information
    X = ml.solve(R, tol=1e-8)               # solve Ax=b to a tolerance of 1e-8
    print("residual norm is", norm(R - A*X))  # c
    print('Using ruge_stuben_solver.........', tm.time()-tstart)

    
    # Method 2 
    #X = spsolve(A,R)
    #print "residual norm is", norm(R - A*X)  # c
    #print 'Using spsolve.........', tm.time()-tstart
    #tstart = tm.time()  
    
    
    #######################################################  
    # reorder results in (i,j,k)
    
    w = np.zeros((nx,ny,nz))
    for i in range(nx): 
        for j in range(ny): 
            for k in range(nz-1):
                idx = i*ny*(nz-1) + k*ny + j; 
                w[i,j,k] = X[idx];
    
    
    return w
    
####################################################### 

######################################################  

def omega_matrix(pm,pn,depths,f,N2):    
    
    
    
    # elliptic equation matrix divided by N^2: (f/N)^2 d_zz + d_xx + d_yy
    dx =1/np.mean(pm)
    dy =1/np.mean(pn)   
    dz = depths[1]-depths[0]
    
    nz=len(depths)
    [nx,ny] = pm.shape

    ############################
    

    
    dx2i = 1./(dx*dx);
    dy2i = 1./(dy*dy);
    dz2i = 1./(dz*dz);


    ndim = (nz-1)*ny*nx;
    i_s  = ny*(nz-1);
    js  = 1;
    bjs = ny-1;
    ks  = ny;
    #A = np.zeros((ndim,ndim));
    #A=csc_matrix((ndim,ndim))
    A=lil_matrix((ndim,ndim))
    
    ############################
    
    for i in range(nx): 
    
        for j in range(ny): 
        
            for k in range(nz-1):
                
                if len(N2.shape)==3: 
                    f2N2 = f[i,j]**2/N2[i,j,k]
                else: 
                    f2N2 = f[i,j]**2/N2[k]
                    
                if len(pm.shape)==2: 
                    dx2i= pm[i,j]**2
                    dy2i= pn[i,j]**2
      
                
                idx = i*ny*(nz-1) + k*ny + j;
                diag = 0.;

                if j>0:
                    A[idx,idx-js] = dy2i;
                    diag = diag - dy2i;
                #else:
                    #A[idx,idx+bjs] = dy2i;
                    #diag = diag - dy2i;

                if k>0:
                    dz2m = 1./((depths[k]-depths[k-1])*0.5*(depths[k+1]- depths[k-1]))
                    A[idx,idx-ks] = f2N2*dz2m;
                    diag = diag - f2N2*dz2m;
                else:
                    dz2m = 1./((depths[1]-depths[0])**2)
                    diag = diag - f2N2*dz2m;

                if i>0:
                    A[idx,idx-i_s] = dx2i;
                    diag = diag - dx2i;

                if i<nx-1:
                    A[idx,idx+i_s] = dx2i;
                    diag = diag - dx2i;

                if k==0:
                    dz2p = 1./((depths[k+1]-depths[k])*(depths[k+1]- depths[k]))
                    A[idx,idx+ks] = f2N2*dz2p;
                    diag = diag - f2N2*dz2p;
                elif k<nz-2:                   
                    dz2p = 1./((depths[k+1]-depths[k])*0.5*(depths[k+1]- depths[k-1]))
                    A[idx,idx+ks] = f2N2*dz2p;
                    diag = diag - f2N2*dz2p;                   
                else:
                    dz2p = 1./((depths[k+1]-depths[k])*0.5*(depths[k+1]- depths[k-1]))
                    diag = diag - f2N2*dz2p;
                    

                if j<ny-1:
                    A[idx,idx+js] = dy2i;
                    diag = diag - dy2i;
                #else:
                    #A[idx,idx-bjs] = dy2i;
                    #diag = diag - dy2i;

                A[idx,idx] = diag;

                

    return A
    
    
    
'''   
#######################################################
# Compute solution of TTW equation
# 
#######################################################

#from pyamg import *
#from pyamg.gallery import *
from scipy import *
from scipy.linalg import *  
from scipy.sparse import *
from scipy.sparse.linalg import spsolve
import time as tm
    
def solve_ttw(bx,by,AKv0,sustr,svstr,f,pm,pn,depths,timing=False):
    

    if timing: tstart = tm.time() 
    #######################################################
    #Create forcing (bx,by)
    #######################################################
    new = np.zeros(AKv0.shape)

    nz=len(depths); [nx,ny] = pm.shape; 
    print 'number of points is ',nx,ny,nz
    dz = depths[1]-depths[0]
    dz2 = dz**2
    ks = 2;
    
    #get gradients
    #bx,by = copy(new),copy(new)
    #bx[1:-1,:,:]= diffx(buoy,pm,2); by[:,1:-1,:] = diffy(buoy,pn,2)
    
    AKv = np.zeros((AKv0.shape[0],AKv0.shape[1],AKv0.shape[2]+1))
    AKv[:,:,:-1] = AKv0[:]; AKv[:,:,-1] = AKv0[:,:,-2]; AKv[:,:,-2] = AKv0[:,:,-2]
    del AKv0
    
    # Solutions
    uz,vz = copy(new)*np.nan,copy(new)

    #######################################################
    # Create Matrix
    #######################################################

    ndim = 2*(nz+1);
    A = lil_matrix((ndim,ndim))
    R = np.zeros(ndim);

    for i in range(nx): 
    
        if i%20==0: print 'solving equation:', round(100.* i/(nx-1)) , ' %'
        
        for j in range(ny): 
        
            A = lil_matrix((ndim,ndim))
            
            #idx = 0
            #A[idx,idx+1] =  f[i,j];
            #A[idx+1,idx] = -f[i,j];

            #for k in range(1,nz):
                #idx = 2*k;
                #A[idx,idx+ks] = AKv[i,j,k+1]/dz2;
                #A[idx,idx] =-AKv[i,j,k]/dz2 - AKv[i,j,k]/dz2;
                #A[idx,idx-ks] = AKv[i,j,k-1]/dz2;
                #A[idx,idx+1] = f[i,j];
        
                #A[idx+1,idx+1+ks] = AKv[i,j,k+1]/dz2;
                #A[idx+1,idx+1] =-AKv[i,j,k]/dz2 - AKv[i,j,k]/dz2;
                #A[idx+1,idx+1-ks] = AKv[i,j,k-1]/dz2;
                #A[idx+1,idx] =-f[i,j];
    
            #idx = 2*nz;
            #A[idx,idx] = AKv[i,j,-1];
            #A[idx+1,idx+1] = AKv[i,j,-1];

            #######################################################

            idx = 0
            A[idx+1,idx] =  f[i,j];
            A[idx,idx+1] = -f[i,j];

            for k in range(1,nz):
                idx = 2*k;
                A[idx+ks,idx] = AKv[i,j,k+1]/dz2;
                A[idx,idx] =-AKv[i,j,k]/dz2 - AKv[i,j,k]/dz2;
                A[idx-ks,idx] = AKv[i,j,k-1]/dz2;
                A[idx+1,idx] = f[i,j];
        
                A[idx+1+ks,idx+1] = AKv[i,j,k+1]/dz2;
                A[idx+1,idx+1] =-AKv[i,j,k]/dz2 - AKv[i,j,k]/dz2;
                A[idx+1-ks,idx+1] = AKv[i,j,k-1]/dz2;
                A[idx,idx+1] =-f[i,j];
    
            idx = 2*nz;
            A[idx+1,idx] = AKv[i,j,-1];
            A[idx,idx+1] = AKv[i,j,-1];
            A = A*-1
            
            #######################################################

            for k in range(nz):   
                idx = 2*k
                R[idx] = bx[i,j,k]
                R[idx+1] = by[i,j,k]
                
   
            # add winds (for Ekman effect)
            #idx = 2*nz-4    
            #R[idx] = R[idx]-sustr[i,j]
            #R[idx+1] = R[idx+1]-svstr[i,j]            
            
            #idx = 2*nz-2    
            #R[idx] = R[idx]+2*sustr[i,j]
            #R[idx+1] = R[idx+1]+2*svstr[i,j]       
            
            #idx = 2*nz-2    
            #R[idx] = R[idx]-sustr[i,j]
            #R[idx+1] = R[idx+1]-svstr[i,j]       
                  
            idx = 2*nz    
            R[idx] = sustr[i,j]
            R[idx+1] = svstr[i,j]
            
            R[np.isnan(R)]=0
            if timing: print 'Matrix definition OK.........', tm.time()-tstart               
            if timing: tstart = tm.time()         

            #######################################################
            #Solve matrix A
            #######################################################   

            A = A.tocsr() 
            
            if timing: print 'Starting computation.........', tm.time()-tstart
            if timing: tstart = tm.time()   
  

            X = spsolve(A,R)
            del A
            
            if timing: print 'computation OK.........', tm.time()-tstart
            if timing: tstart = tm.time()         
            
            #######################################################
            
            # reorder results in (i,j,k)
            for k in range(nz):
                idx = 2*k
                uz[i,j,k] = X[idx];
                vz[i,j,k] = X[idx+1];
                
            if timing: print 'allocation OK.........', tm.time()-tstart               
            if timing: tstart = tm.time()  
            
    u = np.cumsum(uz,2)*dz
    v = np.cumsum(vz,2)*dz  
    
    return u,v
    
    
'''
    
  
    
#######################################################
# Compute solution of TTW equation
# 
#######################################################

#from pyamg import *
#from pyamg.gallery import *
from scipy import *
from scipy.linalg import *  
from scipy.sparse import *
from scipy.sparse.linalg import spsolve
import time as tm
import scipy.integrate as integrate
 
def solve_ttw(bx,by,AKv0,sustr,svstr,f,pm,pn,depths,timing=False):
    

    if timing: tstart = tm.time() 
    #######################################################
    #Create forcing (bx,by)
    #######################################################
    new = np.zeros(AKv0.shape)

    nz=len(depths); [nx,ny] = pm.shape; 
    print('number of points is ',nx,ny,nz)
    dz = depths[1]-depths[0]
    dz2 = dz**2
    ks = 2;
    
    #get gradients
    #bx,by = copy(new),copy(new)
    #bx[1:-1,:,:]= diffx(buoy,pm,2); by[:,1:-1,:] = diffy(buoy,pn,2)
    
    AKv = np.zeros((AKv0.shape[0],AKv0.shape[1],AKv0.shape[2]+1))
    #AKv[:,:,:-1] = AKv0[:]; AKv[:,:,-1] = AKv0[:,:,-2]; AKv[:,:,-2] = AKv0[:,:,-2]
    AKv[:,:,:-1] = AKv0[:]; AKv[:,:,-1] = AKv0[:,:,-3]; AKv[:,:,-2] = AKv0[:,:,-2]  
    #AKv[:,:,:-1] = AKv0[:,:,:]; AKv[:,:,-1] = AKv0[:,:,-1]; del AKv0  


    
    # Solutions
    uz,vz = copy(new)*np.nan,copy(new)

    #######################################################
    # Create Matrix
    #######################################################

    ndim = 2*(nz+1);
    #A = lil_matrix((ndim,ndim))
    #R = np.zeros(ndim);

    for i in range(nx): 
    
        if i%20==0: print('solving equation:', round(100.* i/(nx-1)) , ' %')
        
        for j in range(ny): 
        
            A = lil_matrix((ndim,ndim))
            R = np.zeros(ndim);
            
            idx = 0
            A[idx,idx+1] =  f[i,j];
            A[idx+1,idx] = -f[i,j];

            for k in range(1,nz):
                idx = 2*k;
                A[idx,idx+ks] = AKv[i,j,k+1]/dz2;
                A[idx,idx] =-AKv[i,j,k]/dz2 - AKv[i,j,k]/dz2;
                A[idx,idx-ks] = AKv[i,j,k-1]/dz2;
                A[idx,idx+1] = f[i,j];
        
                A[idx+1,idx+1+ks] = AKv[i,j,k+1]/dz2;
                A[idx+1,idx+1] =-AKv[i,j,k]/dz2 - AKv[i,j,k]/dz2;
                A[idx+1,idx+1-ks] = AKv[i,j,k-1]/dz2;
                A[idx+1,idx] =-f[i,j];
    
            idx = 2*nz;
            A[idx,idx+1] = AKv[i,j,-1];
            A[idx+1,idx] = AKv[i,j,-1];


            
            #######################################################

            for k in range(nz):   
                idx = 2*k
                R[idx] = bx[i,j,k]
                R[idx+1] = by[i,j,k]

            idx = 2*nz    
            R[idx] = svstr[i,j]
            R[idx+1] = sustr[i,j]

            
            if i==20 and j==80: 
                print(R)
                print(AKv[i,j,:])         
            if timing: print('Matrix definition OK.........', tm.time()-tstart)               
            if timing: tstart = tm.time()         

            #######################################################
            #Solve matrix A
            #######################################################   

            A = A.tocsr() 
            
            if timing: print('Starting computation.........', tm.time()-tstart)
            if timing: tstart = tm.time()   
  

            X = spsolve(A,R)
            del A
            
            if timing: print('computation OK.........', tm.time()-tstart)
            if timing: tstart = tm.time()         
            
            #######################################################
            
            # reorder results in (i,j,k)
            for k in range(nz):
                idx = 2*k
                uz[i,j,k] = X[idx];
                vz[i,j,k] = X[idx+1];
                
            if timing: print('allocation OK.........', tm.time()-tstart)               
            if timing: tstart = tm.time()  
            
    #u = np.cumsum(uz,2)*dz
    #v = np.cumsum(vz,2)*dz  
    
    u,v = copy(new), copy(new)
    print(uz.shape,len(depths))
    u[:,:,1:] = integrate.cumtrapz(uz,dx=dz, axis=2 )
    v[:,:,1:] = integrate.cumtrapz(vz,dx=dz, axis=2 )
    
    #ut[:,:,:-1] = integrate.cumtrapz(uz[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1]
    #vt[:,:,:-1] = integrate.cumtrapz(vz[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1]
    #ug[:,:,:-1] = (integrate.cumtrapz(-by[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1].T/f.T).T
    #vg[:,:,:-1] = (integrate.cumtrapz(bx[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1].T/f.T).T
    
    return u,v
    
      
 
    
#######################################################
# Compute solution of TTW equation on sigma levels
# 
#######################################################


from scipy.sparse import *
from scipy.sparse.linalg import spsolve
import scipy.integrate as integrate
import time as tm
 
def solve_ttw_sig(bx,by,AKv,sustr,svstr,f,pm,pn,z_w,timing=False,debug=0):
    
    '''
    AKv and by,by are on vertical w levels
    
    uz,vz (solutions of TTW) also

    '''

    
    if timing: tstart = tm.time() 
    #######################################################
    #Create forcing (bx,by)
    #######################################################
    new = np.zeros(bx.shape)

    nz=AKv.shape[2]-1;  [nx,ny] = pm.shape; 
    print('number of points is ',nx,ny,nz)
    dz = z_w[:,:,1:] - z_w[:,:,:-1]
    ks = 2;

    # Solutions
    uz,vz = copy(new),copy(new)
    
    #######################################################
    # Create Matrix
    #######################################################

    ndim = 2*(nz+1);
    #A = lil_matrix((ndim,ndim))
    #R = np.zeros(ndim);

    for i in range(nx): 
    
        if i%20==0: print('solving equation:', round(100.* i/(nx-1)) , ' %')
        
        for j in range(ny): 
        
            A = lil_matrix((ndim,ndim))
            R = np.zeros(ndim);
            
            idx = 0
            A[idx,idx+1] =  f[i,j];
            A[idx+1,idx] = -f[i,j];

            for k in range(1,nz):
                idx = 2*k;
                dz2 = 0.5*(dz[i,j,k]+dz[i,j,k-1])
                
                A[idx,idx+ks] = AKv[i,j,k+1]/dz[i,j,k]/dz2;
                A[idx,idx] =-AKv[i,j,k]/dz[i,j,k]/dz2 - AKv[i,j,k]/dz[i,j,k-1]/dz2;
                A[idx,idx-ks] = AKv[i,j,k-1]/dz[i,j,k-1]/dz2;
                A[idx,idx+1] = f[i,j];
                
                A[idx+1,idx+1+ks] = AKv[i,j,k+1]/dz[i,j,k]/dz2;
                A[idx+1,idx+1] =-AKv[i,j,k]/dz[i,j,k]/dz2 - AKv[i,j,k]/dz[i,j,k-1]/dz2;
                A[idx+1,idx+1-ks] = AKv[i,j,k-1]/dz[i,j,k-1]/dz2;
                A[idx+1,idx] =-f[i,j];
    
            idx = 2*nz;
            A[idx,idx+1] = AKv[i,j,-1];
            A[idx+1,idx] = AKv[i,j,-1];
     
            #######################################################

            for k in range(nz):   
                idx = 2*k
                R[idx] = bx[i,j,k]
                R[idx+1] = by[i,j,k]

            idx = 2*nz    
            R[idx] = svstr[i,j]
            R[idx+1] = sustr[i,j]

                
            if timing: print('Matrix definition OK.........', tm.time()-tstart)               
            if timing: tstart = tm.time()         

            #######################################################
            #Solve matrix A
            #######################################################   

            A = A.tocsr() 
            
            if timing: print('Starting computation.........', tm.time()-tstart)
            if timing: tstart = tm.time()   
  

            X = spsolve(A,R)
            del A
            
            if timing: print('computation OK.........', tm.time()-tstart)
            if timing: tstart = tm.time()         
            
            #######################################################
            
            # reorder results in (i,j,k)
            for k in range(nz+1):
                idx = 2*k
                uz[i,j,k] = X[idx];
                vz[i,j,k] = X[idx+1];
                
            if timing: print('allocation OK.........', tm.time()-tstart)               
            if timing: tstart = tm.time()  
            
            
    #######################################################
    # Integrate vertically to get u,v,ug,vg
    print(uz.shape, z_w.shape)
    #print uz[10,10,10],bx[10,10,10], by[10,10,10], AKv[10,10,10]
    
    ut,vt,ug,vg = copy(new), copy(new), copy(new), copy(new)
    ut[:,:,1:] = integrate.cumtrapz(uz,z_w, axis=2 )
    vt[:,:,1:] = integrate.cumtrapz(vz,z_w, axis=2 )
    ug[:,:,1:] = (integrate.cumtrapz(-by,z_w, axis=2 ).T/f.T).T
    vg[:,:,1:] = (integrate.cumtrapz(bx,z_w, axis=2 ).T/f.T).T
    
    #ut[:,:,:-1] = integrate.cumtrapz(uz[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1]
    #vt[:,:,:-1] = integrate.cumtrapz(vz[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1]
    #ug[:,:,:-1] = (integrate.cumtrapz(-by[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1].T/f.T).T
    #vg[:,:,:-1] = (integrate.cumtrapz(bx[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1].T/f.T).T
    
    if debug==1:
        return ut,vt,ug,vg,uz,vz
    else:
        return ut,vt,ug,vg  
    
    
    
    
    
    
   