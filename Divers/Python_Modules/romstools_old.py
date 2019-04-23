#!/usr/bin/python
# Filename: romstools.py


#Netcdf IO module
#from Scientific.IO.NetCDF import *
from netCDF4 import Dataset

#module for numerics
import numpy as np

#module for copy
#from copy import copy

#module for plots and colormaps
#import matplotlib.pyplot as py
#import matplotlib.colors as col

#load fortran routines
import romstoolsfort_old as romsF

####################
from itertools import product


#################################################
# get_depths NEW VERSION (similar to z.slice.F)
#################################################


def get_depths(fname,gname,tindex,type):

    #UCLA version 2
    if isinstance(gname, str):
        nc = Dataset(gname, 'r', format='NETCDF3_CLASSIC')
        h = np.array(nc.variables['h'])
        nc.close()
    else:
        h = gname

    if isinstance(fname, str):
        nc = Dataset(fname, 'r', format='NETCDF3_CLASSIC'); opened = True
    else:
        nc = fname; opened = False

    hc=nc.hc
    N=len(nc.dimensions['s_rho'])

    if type=='w':
        Cs = nc.Cs_w
        N=N+1
    else:
        Cs = nc.Cs_r   
 
    z=zlevs(h,hc,N,Cs,type)

    if opened==True: nc.close()

    return z

#################################################
# get_depths NEW VERSION (similar to z.slice.F)
#################################################


def get_depths_zoom(fname,gname,tindex,type,ny1,ny2,nx1,nx2):

    #UCLA version 2

    if isinstance(gname, str):
        nc = Dataset(gname, 'r', format='NETCDF3_CLASSIC')
        h = np.array(nc.variables['h'][ny1:ny2,nx1:nx2])
        nc.close()
    else:
        h = gname
    

    if isinstance(fname, str):
        nc = Dataset(fname, 'r', format='NETCDF3_CLASSIC'); opened = True
    else:
        nc = fname; opened = False

    zeta = np.array(nc.variables['zeta'][tindex,ny1:ny2,nx1:nx2])

    hc=nc.hc
    N=len(nc.dimensions['s_rho'])

    if type=='w':
        Cs = nc.Cs_w
        N=N+1
    else:
        Cs = nc.Cs_r   
 
    if zeta.shape==h.shape:
        z=zlevs_new(h,zeta,hc,N,Cs,type) 
    else:
        z=zlevs(h,hc,N,Cs,type)


    if opened==True: nc.close()

    return z

##################################################

def zlevs(h,hc,N,Cs,type):

    [M,L]=h.shape
    z=np.zeros((N,M,L))

    if type=='w': 
        z[N-1,:,:] = 0.

        for k in range(N-1):

            cff= hc*((k+1-N))/N
            z[k,:,:]=h*(cff+Cs[k]*h)/(hc+h)


    else:
        cff= -0.5*hc/N; z[N-1,:,:]=h*(cff+Cs[N-1]*h)/(hc+h)

        for k in range(N-1):

            cff= hc*((k+1-N)-0.5)/N
            z[k,:,:]=h*(cff+Cs[k]*h)/(hc+h)

    return z
    
    

##################################################
'''
def zlevs_new(h,zeta,hc,N,Cs,type):

    [M,L]=h.shape

    z=np.zeros((N,M,L))

    if type=='w': 
        z[N-1,:,:] = 0.

        for k in range(N-1):

            cff= hc*((k+1-N))/N
            z[k,:,:]=zeta+(zeta+h)*(cff+Cs[k]*h)/(hc+h)

    else:
        cff= -0.5*hc/N; z[N-1,:,:]=zeta+(zeta+h)*(cff+Cs[N-1]*h)/(hc+h)

        for k in range(N-1):

            cff= hc*((k+1-N)-0.5)/N
            z[k,:,:]=zeta+(zeta+h)*(cff+Cs[k]*h)/(hc+h)

    return z
'''

    

##################################################

def zlevs_new(h,zeta,hc,N,Cs,type):

    [M,L]=h.shape

    z=np.zeros((N,M,L))

    if type=='w': 
    
        z[0,:,:] = - h

        for k in range(1,N):

            cff= hc*((k+1-N))/(N-1)
            z[k,:,:]=zeta+(zeta+h)*(cff+Cs[k]*h)/(hc+h)
        
    else:

        for k in range(N):

            cff= hc*((k+1-N)-0.5)/N
            z[k,:,:]=zeta+(zeta+h)*(cff+Cs[k]*h)/(hc+h)
   
            
    return z

    

  
##################################################
# Exact version (identical to zlevs.F)
#################################################



def zlevs_latest(h,zeta, hc, Cs_r, Cs_w):
    
    N=Cs_r.shape[0]
    z_r=np.zeros((N,h.shape[0],h.shape[1]))
    
    for k in range(N):
        cff= hc*((k+1-N)-0.5)/N
        z_r[k,:,:]=zeta+(zeta+h)*(cff+Cs_r[k]*h)/(hc+h)
        
    z_w=np.zeros((N+1,h.shape[0],h.shape[1]))
    z_w[0,:,:] = - h
    z_w[-1,:,:] = zeta
    
    for k in range(1,N+1):
        cff= hc*((k-N))/N
        z_w[k,:,:]=zeta+(zeta+h)*(cff+Cs_w[k]*h)/(hc+h)
        
    return z_r,z_w

      
#################################################
# get_depths OLD VERSION (similar to romstools)
#################################################



def get_depths_zoomold(fname,gname,tindex,type,ny1,ny2,nx1,nx2):

    #UCLA version 2
    if isinstance(gname, str):
        nc = Dataset(gname, 'r', format='NETCDF3_CLASSIC'); 
        h = np.array(nc.variables['h'][ny1:ny2,nx1:nx2])
        nc.close()
    else:
        h = gname;


    if isinstance(fname, str):
        nc = Dataset(fname, 'r', format='NETCDF3_CLASSIC'); opened=True
    else:
        nc = fname; opened=False

    zeta=np.array(nc.variables['zeta'][tindex,ny1:ny2,nx1:nx2])
    theta_s=nc.theta_s
    theta_b=nc.theta_b
    hc=nc.hc

    N=len(nc.dimensions['s_rho'])

    try:
        if type=='w':
            Cs = nc.Cs_w
        else:
            Cs = nc.Cs_r   
    except:
        Cs = None       

    z=zlevsold(h,zeta,theta_s,theta_b,hc,N,type,'new2008')

    if opened==True: nc.close()

    return z



##################################################


def CSF(sc,theta_s,theta_b):

    if theta_s>0.:
        csrf=(1.-np.cosh(theta_s*sc))/(np.cosh(theta_s)-1.)
    else:
        csrf=-sc*sc

    sc1=csrf+1.

    if theta_b>0.:
        Cs =(np.exp(theta_b*sc1)-1.)/(np.exp(theta_b)-1.) -1.
    else:
        Cs = csrf

    return Cs

##################################################


def zlevsold(h,zeta,theta_s,theta_b,hc,N,type,scoord):

	[M,L]=h.shape

# Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
	cff1=1./np.sinh(theta_s)
	cff2=0.5/np.tanh(0.5*theta_s)

	if type=='w':
  		sc=(np.arange(0.,N+1)-N)/N
  		N=N+1
	else:
		sc=(np.arange(1.,N+1)-N-0.5)/N

	if scoord=='new2008':
		Cs = CSF(sc,theta_s,theta_b)
	else:
		Cs=(1.-theta_b)*cff1*np.sinh(theta_s*sc)+theta_b*(cff2*np.tanh(theta_s*(sc+0.5))-0.5)

	hinv=1./(h+hc)
	cff=hc*sc
	cff1=Cs

	z=np.zeros((N,M,L))
	for k in range(N):
		z0=cff[k]+cff1[k]*h
		z[k,:,:]=zeta+(zeta+h)*(cff[k]+cff1[k]*h)*hinv

	return z




#######################################################
#interpolate a 3D variable on a horizontal level of constant depth
#######################################################

def vinterp(var,z,depth,topo=None,cubic=0):

    [N,Mp,Lp]=z.shape

    #depth=0, just take surface field
    if depth>0: 
        varz = np.nan

    #Simple linear interpolation
    elif cubic==0:

        levs2=sum(z<depth)-1
        levs2[levs2==N-1]=N-2
        levs2[levs2==-1]=0
        levs1=levs2+1

        X,Y=np.meshgrid(np.arange(0,Lp),np.arange(0,Mp))

        pos1=levs1,Y,X
        pos2=levs2,Y,X

        z1=z[pos1]
        z2=z[pos2]

        v1=var[pos1]
        v2=var[pos2]
	
        varz = (((v1-v2)*depth+v2*z1-v1*z2)/(z1-z2))
        if topo is not None: varz[depth<-1*topo]=np.nan

    #Cubic interpolation (see ShchepetkinMcWilliams08.pdf)
    elif cubic==1:


        print('cubic interpolation')

        #find the closest level BELOW depth
        levs2=sum(z<depth)-1
        levs1=copy(levs2)
        #levs2[levs1==N-1]=N-2

        #cubic interpolation will use 4 values of var and z in the vertical (2 below and 2 above depth)
        Nlev = 4 

        #prepare arrays for intermediate variables:
        X,Y=np.meshgrid(np.arange(0,Lp),np.arange(0,Mp))
        levs=np.zeros((Nlev,Mp,Lp),int); Xlev=np.zeros((Nlev,Mp,Lp),int); Ylev=np.zeros((Nlev,Mp,Lp),int)
        for ilev in range(Nlev):
            levs[ilev,:,:]=levs2+ilev-1
            Xlev[ilev,:,:]=X
            Ylev[ilev,:,:]=Y


        levs[levs>N-1]=N-1
        levs[levs<0]=0

        pos=levs,Y,X; zz=z[pos]; vark=var[pos]


        #######################################################

        test0=np.zeros((Mp,Lp)); test0[levs2==-1]=1; 
        test1=np.zeros((Mp,Lp)); test1[levs2==0]=1;
        testN1=np.zeros((Mp,Lp)); testN1[levs2==N-2]=1; 
        testN=np.zeros((Mp,Lp)); testN[levs2==N-1]=1;
 
        #######################################################

        zz[1:-1,:,:] = testN * zz[:-2,:,:] + test0 * zz[2:,:,:] + (1 - test0 - testN) * zz[1:-1,:,:]

        dzz = zz[1:,:,:]- zz[:-1,:,:]; 
        dzz[-1,:,:] = testN1 * dzz[1,:,:] + (1-testN1)* dzz[-1,:,:]
        dzz[0,:,:] = test1 * dzz[1,:,:] + (1-test1)* dzz[0,:,:]

        vark[1:-1,:,:] = testN * vark[:-2,:,:] + test0 * vark[2:,:,:] + (1 - test0 - testN) * vark[1:-1,:,:]

        dvark = vark[1:,:,:]-vark[:-1,:,:]; 
        dvark[-1,:,:] = testN1 * dvark[1,:,:] + (1-testN1)* dvark[-1,:,:]
        dvark[0,:,:] = test1 * dvark[1,:,:] + (1-test1)* dvark[0,:,:]


        FC0 = (dvark[1:,:,:]+dvark[:-1,:,:])*dzz[1:,:,:]*dzz[:-1,:,:]
        FC1 = (dzz[1:,:,:]+dzz[:-1,:,:])*dvark[1:,:,:]*dvark[:-1,:,:]
        val=dvark[1:,:,:]*dvark[:-1,:,:]

        FC0[val<=0]=1; FC1[val<=0]=0; FC = FC1/FC0


        #######################################################
        cff = 1/dzz[1,:,:]; p=depth-zz[1,:,:]; q=zz[2,:,:]-depth

        varz = cff*(q*vark[1,:,:]+p*vark[2,:,:]- (1-test0-testN) * cff*p*q*(cff*(q-p)*dvark[1,:,:]+p*FC[1,:,:]-q*FC[0,:,:]))     


        #######################################################

        #varz[depth<-1*topo]=np.nan

   
    return varz





#######################################################
#interpolate a 3D variable on horizontal levels of constant depths 
#######################################################

def vinterps(var,z,depths,topo, cubic=0):

    [N,Mp,Lp]=var.shape
    Nz=len(depths)

    #if var not on rho-grid: interpolate z and topo to the same grid than var (u,v,or psi)
    if var.shape!=z.shape:
        if (var.shape[1]==z.shape[1]-1) and (var.shape[2]==z.shape[2]-1):
            z = rho2psi(z); topo = rho2psi(topo)
        elif (var.shape[1]==z.shape[1]-1):
            z = rho2v(z); topo = rho2v(topo)
        elif (var.shape[2]==z.shape[2]-1):
            z = rho2u(z); topo = rho2u(topo)

    if len(depths)==1:
        vnew=vinterp(var,z,depths[0],topo,cubic)

    else:
        [N,Mp,Lp]=var.shape; Nz=len(depths); vnew=np.zeros((Nz, Mp,Lp))
        for iz in range(0, Nz, 1):
            vnew[iz,:,:]=vinterp(var,z,depths[iz],topo,cubic)

    return vnew




#######################################################
#interpolate a 3D variable on horizontal levels of constant depths (FORTRAN version, much faster)
#######################################################

def vinterpsFold(var,depths, z_r,z_w=None,mask=None):

    if mask is None:  mask = np.ones((z_r.shape[2],z_r.shape[1])); mask[z_r[-1,:,:]==0] = 0

    if z_w is None: 
        print('no z_w specified')
        z_w=np.zeros((z_r.shape[0]+1,z_r.shape[1],z_r.shape[2]))

    if var.shape[0]==z_r.shape[0]: kmin=1
    else: kmin=0

    if var.shape[1]==z_r.shape[1]-1: jmin=1
    else: jmin=0

    if var.shape[2]==z_r.shape[2]-1: imin=1
    else: imin=0

    vnew=romsF.sigma_to_z_introld(z_r.T, z_w.T,mask.T,var.T,depths,imin,jmin,kmin,9999.)

    vnew[vnew==9999.]=np.nan
    vnew[vnew==-9999.]=np.nan

    vnew = vnew.T

    return vnew




#######################################################
#interpolate a 3D variable on horizontal levels of constant depths (FORTRAN version, much faster)
#######################################################

def vinterpsF(var, newz, z_r, z_w=None,mask=None):

    if mask is None:  mask = np.ones((z_r.shape[1],z_r.shape[2])); mask[z_r[-1,:,:]==0] = 0

    if z_w is None: 
        print('no z_w specified')
        z_w=np.zeros((z_r.shape[0]+1,z_r.shape[1],z_r.shape[2]))


    if var.shape[0]==z_r.shape[0]: kmin=1
    else: kmin=0

    if var.shape[1]==z_r.shape[1]-1: jmin=1
    else: jmin=0

    if var.shape[2]==z_r.shape[2]-1: imin=1
    else: imin=0

    if np.rank(newz)==1:
        vnew=romsF.sigma_to_z_introld(z_r.T, z_w.T,mask.T,var.T,newz,imin,jmin,kmin,9999.)
    else:
        vnew=romsF.sigma_to_z_intr(z_r.T, z_w.T,mask.T,var.T,newz.T,imin,jmin,kmin,9999.)

    vnew[vnew==9999.]=np.nan
    vnew[vnew==-9999.]=np.nan

    vnew = vnew.T

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


    [Nz,Mz,Lz]=var_psi.shape
    var_rho=np.zeros((Nz,Mz+1,Lz+1))

    for iz in range(0, Nz, 1):    
        var_rho[iz,:,:]=psi2rho_2d(var_psi[iz,:,:])


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

    var_psi = 0.25*(var_rho[:,1:,1:]+var_rho[:,1:,:-1]+var_rho[:,:-1,:-1]+var_rho[:,:-1,1:])

    return var_psi


#######################################################
#Transfert a 2 or 3-D field at rho points to u points
#######################################################

def rho2u(var_rho):

    if np.rank(var_rho)<3:
        var_u = rho2u_2d(var_rho)
    else:
        var_u = rho2u_3d(var_rho)

    return var_u

def rho2u_2d(var_rho):

    [Mp,Lp]=var_rho.shape
    L=Lp-1
    var_u=0.5*(var_rho[:,0:L]+var_rho[:,1:Lp])

    return var_u


def rho2u_3d(var_rho):

    [N,Mp,Lp]=var_rho.shape
    L=Lp-1
    var_u=0.5*(var_rho[:,:,0:L]+var_rho[:,:,1:Lp])

    return var_u

#######################################################
#Transfert a 3-D field at rho points to v points
#######################################################

def rho2v(var_rho):

    if np.rank(var_rho)<3:
        var_v = rho2v_2d(var_rho)
    else:
        var_v = rho2v_3d(var_rho)

    return var_v

#######################################################

def rho2v_2d(var_rho):

    [Mp,Lp]=var_rho.shape
    M=Mp-1
    var_v=0.5*(var_rho[0:M,:]+var_rho[1:Mp,:]);

    return var_v

#######################################################

def rho2v_3d(var_rho):

    [N,Mp,Lp]=var_rho.shape
    M=Mp-1
    var_v=0.5*(var_rho[:,0:M,:]+var_rho[:,1:Mp,:]);

    return var_v

#######################################################
#Transfert a 2-D field at u points to the rho points
#######################################################

def u2rho(var_u):


    if np.rank(var_u)<3:
        var_rho = u2rho_2d(var_u)
    else:
        var_rho = u2rho_3d(var_u)

    return var_rho

#######################################################

def u2rho_2d(var_u):

    [Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[:,1:L]=0.5*(var_u[:,0:Lm]+var_u[:,1:L])
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Lp-1]=var_rho[:,L-1]
    return var_rho

#######################################################

def u2rho_3d(var_u):

    [N,Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((N,Mp,Lp))
    var_rho[:,:,1:L]=0.5*(var_u[:,:,0:Lm]+var_u[:,:,1:L])
    var_rho[:,:,0]=var_rho[:,:,1]
    var_rho[:,:,Lp-1]=var_rho[:,:,L-1]
    return var_rho


#######################################################
#Transfert a 2 or 2-D field at v points to the rho points
#######################################################

def v2rho(var_v):

    if np.rank(var_v)<3:
        var_rho = v2rho_2d(var_v)
    else:
        var_rho = v2rho_3d(var_v)

    return var_rho

#######################################################

def v2rho_2d(var_v):

    [M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[1:M,:]=0.5*(var_v[0:Mm,:]+var_v[1:M,:])
    var_rho[0,:]=var_rho[1,:]
    var_rho[Mp-1,:]=var_rho[M-1,:]

    return var_rho

#######################################################

def v2rho_3d(var_v):

    [N,M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((N,Mp,Lp))
    var_rho[:,1:M,:]=0.5*(var_v[:,0:Mm,:]+var_v[:,1:M,:])
    var_rho[:,0,:]=var_rho[:,1,:]
    var_rho[:,Mp-1,:]=var_rho[:,M-1,:]

    return var_rho

#######################################################
#Compute vorticity of a 2-D field
#######################################################
"""
    vrt = dv/dx - du/dy

    u,v on originals u and v grids, respectively
    vrt is outputed on psi grid 

"""

def vort(u,v,pm,pn):

    if np.rank(u)<3:
        vrt = vort_2d(u,v,pm,pn)
    else:
        vrt = vort_3d(u,v,pm,pn)

    return vrt

#############################

def vort_2d(u,v,pm,pn):

    [Mp,Lp]=pm.shape
    L=Lp-1
    M=Mp-1

    dm_u=2*u/(pm[:,0:L]+pm[:,1:Lp])
    dn_v=2*v/(pn[0:M,:]+pn[1:Mp,:])

    iA_q=0.0625*(pm[0:M,0:L]+pm[0:M,1:Lp]+pm[1:Mp,1:Lp]+pm[1:Mp,0:L])\
        *(pn[0:M,0:L]+pn[0:M,1:Lp]+pn[1:Mp,1:Lp]+pn[1:Mp,0:L])
    
    vrt=iA_q*(dn_v[:,1:Lp]-dn_v[:,0:L]-dm_u[1:Mp,:]+dm_u[0:M,:])

    return vrt


#############################

def vort_3d(u,v,pm,pn):

    [Nz,Mz,Lz]=u.shape
    [Mp,Lp]=pm.shape

    vrt = np.zeros((Nz,Mp-1,Lp-1))

    for iz in range(0, Nz, 1):    
        vrt[iz,:,:]=vort_2d(u[iz,:,:],v[iz,:,:],pm,pn)

    return vrt


#######################################################
#Compute divergence of a 2-D field
#######################################################
'''
div = dudx + dvdy

div is computed at rho points (first and last points as nan)
'''



def div(u,v,pm,pn):

    [Mp,Lp]=pm.shape
    L=Lp-1
    M=Mp-1
    Lm=L-1
    Mm=M-1


    dudx = np.zeros((Mp,Lp))*np.nan
    dudx[:,1:L] = (u[:,1:L]-u[:,0:Lm])*pm[:,1:L]
    #dudx[:,0] = dudx[:,1]
    #dudx[:,L] = dudx[:,Lm]

    dvdy = np.zeros((Mp,Lp))*np.nan
    dvdy[1:M,:] = (v[1:M,:]-v[0:Mm,:])*pn[1:M,:]
    #dvdy[0,:] = dvdy[1,:]
    #dvdy[M,:] = dvdy[Mm,:]

    var = dudx + dvdy

    return var


def divs(u,v,pm,pn):

    if np.rank(u)>2:

        [Nz,Mz,Lz]=u.shape
        [Mp,Lp]=pm.shape

        divs = np.zeros((Nz,Mp,Lp))*np.nan

        for iz in range(0, Nz, 1):
            divs[iz,:,:]=div(u[iz,:,:],v[iz,:,:],pm,pn)

    else:
    
        divs=div(u,v,pm,pn)

    return divs


#######################################################
#Compute strain of a 2 or 3-D field
#######################################################
'''
 S = sqrt((ux-vy)^2+(vx+uy)^2)

Strain is computed at rho points
'''


def strain(u,v,pm,pn):

    [Mp,Lp]=pm.shape
    L=Lp-1; M=Mp-1; Lm=L-1; Mm=M-1

    'dudx on rho points'
    dudx = pm*0
    dudx[:,1:L] = (u[:,1:L]-u[:,0:Lm])*pm[:,1:L]
    dudx[:,0] = dudx[:,1]
    dudx[:,L] = dudx[:,Lm]

    'dvdx on rho points'
    dvdy = pn*0
    dvdy[1:M,:] = (v[1:M,:]-v[0:Mm,:])*pn[1:M,:]
    dvdy[0,:] = dvdy[1,:]
    dvdy[M,:] = dvdy[Mm,:]

    'dvdx on psi points'
    dvdx = (v[:,1:Lp]-v[:,0:L])*0.25*(pm[1:Mp,1:Lp]+pm[1:Mp,0:L]+pm[0:M,1:Lp]+pm[0:M,0:L])
    'dvdx on rho points'
    dvdx = psi2rho(dvdx)

    'dudy on psi points'
    dudy = (u[1:Mp,:]-u[0:M,:])*0.25*(pn[1:Mp,1:Lp]+pn[1:Mp,0:L]+pn[0:M,1:Lp]+pn[0:M,0:L])
    'dudy on rho points'
    dudy = psi2rho(dudy)

    var=np.sqrt((dudx-dvdy)**2 + (dudy+dvdx)**2)

    return var

#######################################################

def strainold(u,v,pm,pn):



    [Mp,Lp]=pm.shape
    L=Lp-1; M=Mp-1; Lm=L-1; Mm=M-1

    if u.shape!=v.shape:
        u = u2rho(u)
        v = v2rho(v)

    dx = 0.5*sum(sum(pm))/((Mp+1)*(Lp+1))
    dy = 0.5*sum(sum(pn))/((Mp+1)*(Lp+1))

    'dudx on rho points'
    dudx = pm*0
    dudx[:,1-1:L-1] = (u[:,2:Lp]-u[:,0:Lm])*dx
    dudx[:,0] = dudx[:,1]
    dudx[:,L] = dudx[:,Lm]

    'dvdx on rho points'
    dvdy = pn*0
    dvdy[1-1:M-1,:] = (v[2:Mp,:]-v[0:Mm,:])*dy
    dvdy[0,:] = dvdy[1,:]
    dvdy[M,:] = dvdy[Mm,:]

    dvdx = pm*0
    dvdx[:,1-1:L-1] = (v[:,2:Lp]-v[:,0:Lm])*dx
    dvdx[:,0] = dvdx[:,1]
    dvdx[:,L] = dvdx[:,Lm]

    dudy = pn*0
    dudy[1-1:M-1,:] = (u[2:Mp,:]-u[0:Mm,:])*dy
    dudy[0,:] = dudy[1,:]
    dudy[M,:] = dudy[Mm,:]

    var=np.sqrt((dudx-dvdy)**2 + (dudy+dvdx)**2)

    return var

#######################################################

def strains(u,v,pm,pn):


    if np.rank(u)>2:

        [Nz,Mz,Lz]=u.shape
        [Mp,Lp]=pm.shape

        var = np.zeros((Nz,Mp,Lp))

        for iz in range(0, Nz, 1):
            var[iz,:,:]=strain(u[iz,:,:],v[iz,:,:],pm,pn)

    else:
    
        var=strain(u,v,pm,pn)

    return var






#######################################################
#Compute back ground strain...
#######################################################
'''
 S = sqrt((ux-vy)^2+(vx+uy)^2)

Strain is computed at rho points
'''


def bgstrain(u,v,pm,pn):

    [M,L]=pm.shape

    [x,y]=np.meshgrid(np.arange(L),np.arange(M))  
    x=(x-np.min(x))/pm; y=(y-np.min(y))/pn


    vrt =  -psi2rho( vort(u,v,pm,pn))
    
    ubg = np.zeros((pm.shape))
    vbg = np.zeros((pm.shape))

    for i,j in product(list(range(L)),list(range(M))):

        ubg[j,i] = 1/(2*np.pi) *np.nansum(-(y-y[j,i])/((x-x[j,i])**2+(y-y[j,i])**2)*vrt/(pm*pn))
        vbg[j,i] = 1/(2*np.pi) *np.nansum((x-x[j,i])/((x-x[j,i])**2+(y-y[j,i])**2)*vrt/(pm*pn))


    return ubg,vbg


#######################################################
#Compute A-S of a 2 or 3-D field
#######################################################

def ageo(u,v,pm,pn,f):

    xi = f + psi2rho(vort(u,v,pm,pn)) - abs(strain(u,v,pm,pn))

    return xi

#######################################################

def ageos(u,v,pm,pn,f):

    if np.rank(u)>2:

        [Nz,Mz,Lz]=u.shape
        [Mp,Lp]=pm.shape

        xi = np.zeros((Nz,Mp,Lp))

        for iz in range(0, Nz, 1):    
            xi[iz,:,:]=ageo(u[iz,:,:],v[iz,:,:],pm,pn,f)

    else:

        xi=ageo(u,v,pm,pn,f)

    return xi

#######################################################
#Compute A-S of a 2 or 3-D field
#######################################################

def ow(u,v,pm,pn):

    xi =  np.abs(strain(u,v,pm,pn))**2 -  np.abs(psi2rho(vort(u,v,pm,pn)))**2

    return xi

#######################################################

def ows(u,v,pm,pn):

    if np.rank(u)>2:

        [Nz,Mz,Lz]=u.shape
        [Mp,Lp]=pm.shape

        xi = np.zeros((Nz,Mp,Lp))

        for iz in range(0, Nz, 1):    
            xi[iz,:,:]=ow(u[iz,:,:],v[iz,:,:],pm,pn)

    else:

        xi=ow(u,v,pm,pn)

    return xi

#######################################################
#Compute A-S of a 2 or 3-D field
#######################################################

def ow_alter(u,v,pm,pn,part=0):

    ux = np.zeros(pm.shape)*np.nan
    vy = np.zeros(pm.shape)*np.nan
    
    ux[:,1:-1] = diffx(u,rho2u(pm))
    vy[1:-1,:] = diffy(v,rho2v(pn))

    vx = diffx(v,rho2v(pm))
    uy = diffy(u,rho2u(pn))

    if part==0:
        xi =  (ux - vy)**2 + 4 * psi2rho(vx * uy)
    elif part==1:
        xi =  (ux - vy)**2
    elif part==2:   
        xi =  4 * psi2rho(vx * uy)

    return xi

#######################################################

def ows_alter(u,v,pm,pn,part=0):

    if np.rank(u)>2:

        [Nz,Mz,Lz]=u.shape
        [Mp,Lp]=pm.shape

        xi = np.zeros((Nz,Mp,Lp))

        for iz in range(0, Nz, 1):    
            xi[iz,:,:]=ow_alter(u[iz,:,:],v[iz,:,:],pm,pn,part)

    else:

        xi=ow_alter(u,v,pm,pn,part)

    return xi
    
    
#######################################################
#Compute Potential Vorticity of a 3-D field
#######################################################


def PVold(tempz,saltz,uz,vz,z_r,z_w,f,g,rho0,pm,pn):

    uz=u2rho(uz)
    vz=v2rho(vz)
     
    [rho,bvf] = rho_eos(tempz,saltz,z_r,g,rho0,z_w)
    
    dx  = 1/np.mean(pn)
    dy  = 1/np.mean(pm)              
    dxi = 0.5/dx
    dyi = 0.5/dy
    
    dzb = bvf[1:-1,:,:]
    
    dzu = (uz[1:,:,:]-uz[:-1,:,:])/(z_r[1:,:,:]- z_r[:-1,:,:])
    dzv = (vz[1:,:,:]-vz[:-1,:,:])/(z_r[1:,:,:]- z_r[:-1,:,:])
    
    umean = 0.5*(uz[1:,:,:]+uz[:-1,:,:])
    vmean = 0.5*(vz[1:,:,:]+vz[:-1,:,:])
    bmean = -0.5*g*(rho[1:,:,:] + rho[:-1,:,:])/rho0
      
    dyu = 0*bmean
    dxv = 0*bmean
    dxb = 0*bmean
    dyb = 0*bmean

    dyu[:,1:-1,:] = (umean[:,2:,:] - umean[:,:-2,:] )*dyi
    dxv[:,:,1:-1] = (vmean[:,:,2:] - vmean[:,:,:-2] )*dxi
    dxb[:,:,1:-1] = (bmean[:,:,2:] - bmean[:,:,:-2] )*dxi
    dyb[:,1:-1,:] = (bmean[:,2:,:] - bmean[:,:-2,:] )*dyi
    
    pv1 =  (f + dxv - dyu) * dzb
    pv2 =  -1*dzv*dxb 
    pv3 =  dzu*dyb   
    
    pv =  pv1 + pv2 + pv3

    return [pv,pv1,pv2,pv3]




#######################################################
#Compute Potential Vorticity of a 3-D field on psi-w grid
#######################################################
'''

Compute ertel potential vorticity using buoyancy (b=-g rho/rho0)

T and S on horizontal rho grids and vertical rho-grid (specified by z_r)
U and V on horizontal u- and v- grids and vertical rho-grid (specified by z_r)

PV is computed on horizontal psi-grid and vertical w-grid (specified by z_w)

'''

def PV(tempz,saltz,uz,vz,z_r,z_w,f,g,rho0,pm,pn):


    #rho on rho-rho grid      bvf on rho-w grid
    [rho,bvf] = rho_eos(tempz,saltz,z_r,g,rho0,z_w)

    buoy = -g*rho/rho0
    #buoy=rho
    dz_r = z_r[1:,:,:]- z_r[:-1,:,:]

##########################
#Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*db/dz

    #dudy and dvdx on psi-rho grid
    #dvdx = (v[:,:,1:Lp]-v[:,:,0:L])*0.25*(pm[1:Mp,1:Lp]+pm[1:Mp,0:L]+pm[0:M,1:Lp]+pm[0:M,0:L])
    #dudy = (u[:,1:Mp,:]-u[:,0:M,:])*0.25*(pn[1:Mp,1:Lp]+pn[1:Mp,0:L]+pn[0:M,1:Lp]+pn[0:M,0:L])

    #vrt on psi-rho grid
    #vrt = dvdx - dudy
    vrt = vort(uz,vz,pm,pn)

    #dbdz on psi-w grid
    dbdz = rho2psi(bvf)
    #dbdz = rho2psi((buoy[1:,:,:]-buoy[:-1,:,:])/(dz_r))
    #print dbdz.shape
    #print rho2psi(bvf).shape

    #PV1 on psi-w grid
    pv1=z_w[1:-1,1:,1:]*np.nan
    pv1 =  (rho2psi(f) + 0.5*(vrt[1:,:,:]+vrt[:-1,:,:])) * dbdz[1:-1,:,:]
    #pv1 =  (rho2psi(f) + 0.5*(vrt[1:,:,:]+vrt[:-1,:,:])) * dbdz

##########################
#'Ertel potential vorticity, term 2: (dv/dz)*(db/dx)'

    #'dvdz on psi-w grid'
    dvdz = rho2u((vz[1:,:,:]-vz[:-1,:,:])/(0.5*(dz_r[:,1:,:]+ dz_r[:,:-1,:])))

    #'dbdx on psi-rho grid'
    dbdx = rho2v((buoy[:,:,1:]-buoy[:,:,:-1])*0.5*(pm[:,1:]+pm[:,:-1]))

    
    #PV1 on psi-w grid
    pv2 = z_w[1:-1,1:,1:]*np.nan
    pv2 =  -1*dvdz*0.5*(dbdx[1:,:,:]+dbdx[:-1,:,:])

##########################
#'Ertel potential vorticity, term 3: (du/dz)*(db/dy)'

    #'dudz on psi-w grid'
    dudz = rho2v((uz[1:,:,:]-uz[:-1,:,:])/(0.5*(dz_r[:,:,1:]+ dz_r[:,:,:-1])))


    #'dbdy on psi-rho grid'
    dbdy = rho2u((buoy[:,1:,:]-buoy[:,:-1,:])*0.5*(pn[1:,:]+pm[:-1,:]))

    #PV3 on psi-w grid
    pv3 = z_w[1:-1,1:,1:]*np.nan
    pv3 =  dudz*0.5*(dbdy[1:,:,:]+dbdy[:-1,:,:])

##########################

    return pv1 + pv2 + pv3








#######################################################
#Compute Potential Vorticity of a 3-D field on psi-w grid
#######################################################
'''

Compute ertel potential vorticity using buoyancy (b=-g rho/rho0)

T and S on horizontal rho grids and vertical rho-grid (specified by z_r)
U and V on horizontal u- and v- grids and vertical rho-grid (specified by z_r)

PV is computed on horizontal psi-grid and vertical w-grid (specified by z_w)

'''

def PV_sig_bot(temp,salt,u,v,z_r,z_w,f,g,rho0,pm,pn):

    print('bottom PV?')
    if np.nanmin(z_r[0,:,:])==np.nanmax(z_r[0,:,:]): print('you should use PV')

    #rho on rho-rho grid      bvf on rho-w grid
    [rho,bvf] = rho_eosF(temp,salt,z_r,g,rho0,z_w)

    buoy = -g*rho/rho0
    #buoy=rho
    dz_r = z_r[1:,:,:]- z_r[:-1,:,:]

##########################
#Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*db/dz


    #dudy and dvdx on psi-rho grid
    dvdx = diffxi(v,rho2v(pm),rho2v(z_r),rho2v(z_w[1,:,:]))
    dudy = diffeta(u,rho2u(pn),rho2u(z_r),rho2u(z_w[1,:,:]))

    #vrt on psi-rho grid
    vrt = dvdx - dudy

    #dbdz on psi-w grid
    #dbdz = rho2psi(bvf[1,:,:])
    dbdz = rho2psi((buoy[1:,:,:]-buoy[:-1,:,:])/(dz_r))

    #PV1 on psi-w grid
    pv1 =  (rho2psi(f) + vrt) * dbdz


##########################
#'Ertel potential vorticity, term 2: (dv/dz)*(db/dx)'

    #'dvdz on psi-w grid'
    dvdz = rho2u((v[1:,:,:]-v[:-1,:,:])/(0.5*(dz_r[:,1:,:]+ dz_r[:,:-1,:])))

    #'dbdx on psi-rho grid'
    dbdx = rho2v(diffxi(buoy,pm,z_r,z_w[1,:,:]))

    #PV1 on psi-w grid
    pv2 =  -1*dvdz*dbdx

##########################
#'Ertel potential vorticity, term 3: (du/dz)*(db/dy)'

    #'dudz on psi-w grid'
    dudz = rho2v((u[1:,:,:]-u[:-1,:,:])/(0.5*(dz_r[:,:,1:]+ dz_r[:,:,:-1])))


    #'dbdy on psi-rho grid'
    dbdy = rho2u(diffeta(buoy,pn,z_r,z_w[1,:,:]))

    #PV3 on psi-w grid
    pv3 =  dudz*dbdy

##########################

    return pv1 + pv2 + pv3





#######################################################
#Compute Potential Vorticity of a 3-D field on psi-w grid
#######################################################
'''

Compute ertel potential vorticity using buoyancy (b=-g rho/rho0)

T and S on horizontal rho grids and vertical rho-grid (specified by z_r)
U and V on horizontal u- and v- grids and vertical rho-grid (specified by z_r)

PV is computed on horizontal psi-grid and vertical w-grid (specified by z_w)

'''

def PV_sig(temp,salt,u,v,z_r,z_w,f,g,rho0,pm,pn,mask=None):

    print('we are using PV_sig')
    #rho on rho-rho grid      bvf on rho-w grid
    [rho,bvf] = rho_eosF(temp,salt,z_r,g,rho0,z_w)

    buoy = -g*rho/rho0
    #buoy=rho
    dz_r = z_r[1:,:,:]- z_r[:-1,:,:]

##########################
#Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*db/dz


    #dudy and dvdx on psi-rho grid
    dvdx = diffxi(v,rho2v(pm),rho2v(z_r),rho2v(z_w),mask)
    dudy = diffeta(u,rho2u(pn),rho2u(z_r),rho2u(z_w),mask)

    #vrt on psi-rho grid
    vrt = dvdx - dudy

    #dbdz on psi-w grid
    dbdz = rho2psi(bvf[1:-1,:,:])
    #dbdz = rho2psi((buoy[1:,:,:]-buoy[:-1,:,:])/(dz_r))

    #PV1 on psi-w grid
    pv1 =  (rho2psi(f) + 0.5*(vrt[1:,:,:] + vrt[:-1,:,:])) * dbdz
    del vrt,dbdz

##########################
#'Ertel potential vorticity, term 2: (dv/dz)*(db/dx)'

    #'dvdz on psi-w grid'
    dvdz = rho2u((v[1:,:,:]-v[:-1,:,:])/(0.5*(dz_r[:,1:,:]+ dz_r[:,:-1,:])))

    #'dbdx on psi-rho grid'
    dbdx = rho2v(diffxi(buoy,pm,z_r,z_w,mask))

    #PV1 on psi-w grid
    pv2 =  -1*dvdz*0.5*(dbdx[1:,:,:] + dbdx[:-1,:,:])
    del dbdx,dvdz

##########################
#'Ertel potential vorticity, term 3: (du/dz)*(db/dy)'

    #'dudz on psi-w grid'
    dudz = rho2v((u[1:,:,:]-u[:-1,:,:])/(0.5*(dz_r[:,:,1:]+ dz_r[:,:,:-1])))

    #'dbdy on psi-rho grid'
    dbdy = rho2u(diffeta(buoy,pn,z_r,z_w,mask))

    #PV3 on psi-w grid
    pv3 =  dudz*0.5*(dbdy[1:,:,:] + dbdy[:-1,:,:])
    del dbdy,dudz

##########################

    return [pv1+pv2+pv3,pv1,pv2,pv3]




#######################################################
#Compute density and Brunt-Vaissala frequency
#######################################################
'''

compute rho from equation of state

rho on rho (vert and hor) grid

bvf is computed only if z_w is not None
bvf computed on rho-w grid (first and last levels set to 0) 


'''

def rho_eos(Tt,Ts,z_r,g,rho0,z_w=None):

    if np.rank(Tt)==2:
        [M,L]=Tt.shape
    else:
        [N,M,L]=Tt.shape

    A00=+19092.56;A01=+209.8925;
    A02=-3.041638;A03=-1.852732e-3;A04=-1.361629e-5;A10=104.4077;
    A11=-6.500517;A12=+0.1553190;A13=2.326469e-4;AS0=-5.587545;
    AS1=+0.7390729;AS2=-1.909078e-2;B00=+4.721788e-1;B01=+1.028859e-2;
    B02=-2.512549e-4;B03=-5.939910e-7;B10=-1.571896e-2;B11=-2.598241e-4;
    B12=+7.267926e-6;BS1=+2.042967e-3;E00=+1.045941e-5;E01=-5.782165e-10;
    E02=+1.296821e-7;E10=-2.595994e-7;E11=-1.248266e-9;E12=-3.508914e-9;

    QR=+999.842594;Q01=+6.793952e-2;Q02=-9.095290e-3;
    Q03=+1.001685e-4;Q04=-1.120083e-6;Q05=+6.536332e-9;Q10=+0.824493;
    Q11=-4.08990e-3;Q12=+7.64380e-5;Q13=-8.24670e-7;Q14=+5.38750e-9;
    QS0=-5.72466e-3;QS1=+1.02270e-4;QS2=-1.65460e-6;Q20=+4.8314e-4;

    sqrtTs=Ts ** 0.5;
    
    K0=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))\
    +Ts*(A10+Tt*(A11+Tt*(A12+Tt*A13))\
    +sqrtTs*(AS0+Tt*(AS1+Tt*AS2)));
    
    K1=B00+Tt*(B01+Tt*(B02+Tt*B03))\
    +Ts*(B10+Tt*(B11+Tt*B12)+sqrtTs*BS1);
    
    K2=E00+Tt*(E01+Tt*E02)\
    +Ts*(E10+Tt*(E11+Tt*E12));
    
    rho1=QR+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))\
    +Ts*(Q10+Tt*(Q11+Tt*(Q12+Tt*(Q13+Tt*Q14)))\
    +sqrtTs*(QS0+Tt*(QS1+Tt*QS2))+Ts*Q20);

    rho=rho1/(1+0.1*z_r/(K0-z_r*(K1-z_r*K2)));

    #######################################################

    if z_w is not None:

        bvf=0.*z_w;
        cff=g/rho0;

        bvf[1:N,:]=-cff*(rho1[1:N,:,:]/\
        (1.+0.1*z_w[1:N,:,:]/\
        ( K0[1:N,:,:]-z_w[1:N,:,:]*(K1[1:N,:,:]-z_w[1:N,:,:]*K2[1:N,:,:])))\
        -rho1[0:N-1,:,:]/( 1.+0.1*z_w[1:N,:,:]/\
        ( K0[0:N-1,:,:]-z_w[1:N,:,:]*(K1[0:N-1,:,:]-z_w[1:N,:,:]*K2[0:N-1,:,:]))))\
        /(z_r[1:N,:,:]-z_r[0:N-1,:,:]);


        return [rho,bvf]

    else:

        return rho

#######################################################
#Compute density and Brunt-Vaissala frequency
#######################################################
'''

compute rho1 from equation of state

rho1 on rho (vert and hor) grid

bvf is computed only if z_w is not None
bvf computed on rho-w grid (first and last levels set to 0) 


'''

def rho1_eos(Tt,Ts,z_r,g,rho0,z_w=None):

    if np.rank(Tt)==2:
        [M,L]=Tt.shape
    else:
        [N,M,L]=Tt.shape

    A00=+19092.56;A01=+209.8925;
    A02=-3.041638;A03=-1.852732e-3;A04=-1.361629e-5;A10=104.4077;
    A11=-6.500517;A12=+0.1553190;A13=2.326469e-4;AS0=-5.587545;
    AS1=+0.7390729;AS2=-1.909078e-2;B00=+4.721788e-1;B01=+1.028859e-2;
    B02=-2.512549e-4;B03=-5.939910e-7;B10=-1.571896e-2;B11=-2.598241e-4;
    B12=+7.267926e-6;BS1=+2.042967e-3;E00=+1.045941e-5;E01=-5.782165e-10;
    E02=+1.296821e-7;E10=-2.595994e-7;E11=-1.248266e-9;E12=-3.508914e-9;

    QR=+999.842594;Q01=+6.793952e-2;Q02=-9.095290e-3;
    Q03=+1.001685e-4;Q04=-1.120083e-6;Q05=+6.536332e-9;Q10=+0.824493;
    Q11=-4.08990e-3;Q12=+7.64380e-5;Q13=-8.24670e-7;Q14=+5.38750e-9;
    QS0=-5.72466e-3;QS1=+1.02270e-4;QS2=-1.65460e-6;Q20=+4.8314e-4;

    sqrtTs=Ts ** 0.5;
    
    
    rho1=QR+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))\
    +Ts*(Q10+Tt*(Q11+Tt*(Q12+Tt*(Q13+Tt*Q14)))\
    +sqrtTs*(QS0+Tt*(QS1+Tt*QS2))+Ts*Q20);

    return rho1

#######################################################
#Compute Potential Vorticity of a 3-D field on psi-w grid
#######################################################
'''

Compute ertel potential vorticity using buoyancy (b=-g rho/rho0)

T and S on horizontal rho grids and vertical rho-grid (specified by z_r)
U and V on horizontal u- and v- grids and vertical rho-grid (specified by z_r)

PV is computed on horizontal psi-grid and vertical w-grid (specified by z_w)

'''

def PVF(tempz,saltz,uz,vz,z_r,z_w,f,g,rho0,pm,pn):

    #rho on rho-rho grid      bvf on rho-w grid
    rp_r = -1*g*rho_eosF(tempz,saltz,z_r,g,rho0)/rho0
    #[rp_r,rp_w] = rho_eosF(tempz,saltz,z_r,g,rho0)


    #buoy = -g*rho/rho0
    #buoy=rho
    dz_r = z_r[1:,:,:]- z_r[:-1,:,:]

##########################
#Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*db/dz

    #dudy and dvdx on psi-rho grid
    #dvdx = (v[:,:,1:Lp]-v[:,:,0:L])*0.25*(pm[1:Mp,1:Lp]+pm[1:Mp,0:L]+pm[0:M,1:Lp]+pm[0:M,0:L])
    #dudy = (u[:,1:Mp,:]-u[:,0:M,:])*0.25*(pn[1:Mp,1:Lp]+pn[1:Mp,0:L]+pn[0:M,1:Lp]+pn[0:M,0:L])

    #vrt on psi-rho grid
    #vrt = dvdx - dudy
    vrt = vort(uz,vz,pm,pn)

    #dbdz on psi-w grid
    #dbdz = rho2psi(rp_w)
    dbdz = rho2psi((rp_r[1:,:,:]-rp_r[:-1,:,:])/(dz_r))


    #PV1 on psi-w grid
    pv1=z_w[1:-1,1:,1:]*np.nan
    pv1 =  (rho2psi(f) + 0*0.5*(vrt[1:,:,:]+vrt[:-1,:,:])) * dbdz

##########################
#'Ertel potential vorticity, term 2: (dv/dz)*(db/dx)'

    #'dvdz on psi-w grid'
    dvdz = rho2u((vz[1:,:,:]-vz[:-1,:,:])/(0.5*(dz_r[:,1:,:]+ dz_r[:,:-1,:])))

    #'dbdx on psi-rho grid'
    dbdx = rho2v((rp_r[:,:,1:]-rp_r[:,:,:-1])*0.5*(pm[:,1:]+pm[:,:-1]))
    
    #PV1 on psi-w grid
    pv2=z_w[1:-1,1:,1:]*np.nan
    pv2 =  -1*dvdz*0.5*(dbdx[1:,:,:]+dbdx[:-1,:,:])
 
##########################
#'Ertel potential vorticity, term 3: (du/dz)*(db/dy)'

    #'dudz on psi-w grid'
    dudz = rho2v((uz[1:,:,:]-uz[:-1,:,:])/(0.5*(dz_r[:,:,1:]+ dz_r[:,:,:-1])))


    #'dbdy on psi-rho grid'
    dbdy = rho2u((rp_r[:,1:,:]-rp_r[:,:-1,:])*0.5*(pn[1:,:]+pm[:-1,:]))

    #PV3 on psi-w grid
    pv3=z_w[1:-1,1:,1:]*np.nan
    pv3 =  dudz*0.5*(dbdy[1:,:,:]+dbdy[:-1,:,:])

##########################

    pv =  pv1 + pv2 + pv3


    return [pv,pv1,pv2,pv3]


#######################################################
#Compute density and Brunt-Vaissala frequency using ROMS subroutine
#######################################################
'''

compute rho from equation of state

rho on rho (vert and hor) grid

bvf is computed only if z_w is not None
bvf computed on rho-w grid (first and last levels set to 0) 


'''


def rho_eosF(Tt,Ts,z_r,g,rho0,z_w=None):
    #subroutine rho_eos(Lm,Mm,N, T,S, z_r,z_w, rho1,qp1,rho,bvf)

    if z_w is None: 

        z_w=np.zeros((z_r.shape[0]+1,z_r.shape[1],z_r.shape[2]))
        (rho1,qp1,rho,bvf) = romsF.rho_eos(Tt.T,Ts.T, z_r.T,z_w.T,rho0)

        return rho.T

    else:

        print('computing bvf using Fortran routine')
        print(z_w.shape)
        print('in romstools', rho0)
        (rho1,qp1,rho,bvf) = romsF.rho_eos(Tt.T,Ts.T, z_r.T,z_w.T,rho0)
        #(rho1,qp1,rho,bvf) = romsF.rho_eos_alex(Tt.T,Ts.T, z_r.T,z_w.T)

        return [rho.T,bvf.T]

#######################################################
#Compute density and Brunt-Vaissala frequency using ROMS subroutine
#######################################################
'''

compute rho from equation of state

rho on rho (vert and hor) grid

bvf is computed only if z_w is not None
bvf computed on rho-w grid (first and last levels set to 0) 


'''


def rho1_eosF(Tt,Ts,z_r,g,rho0,z_w=None):
    #subroutine rho_eos(Lm,Mm,N, T,S, z_r,z_w, rho1,qp1,rho,bvf)

    if z_w is None: 

        z_w=np.zeros((z_r.shape[0]+1,z_r.shape[1],z_r.shape[2]))
        (rho1) = romsF.rho1_eos(Tt.T,Ts.T, z_r.T,z_w.T,rho0)

        return rho1.T

    else:

        print('computing bvf using Fortran routine')
        print(z_w.shape)
        print('in romstools', rho0)
        (rho1,qp1,rho,bvf) = romsF.rho_eos(Tt.T,Ts.T, z_r.T,z_w.T,rho0)
        #(rho1,qp1,rho,bvf) = romsF.rho_eos_alex(Tt.T,Ts.T, z_r.T,z_w.T)

        return [rho1.T,bvf.T]

#######################################################
#Compute adiabatic density gradients
#######################################################
'''

compute rho from equation of state

rho on rho (vert and hor) grid

drdz computed on rho-w grid (first and last levels set to 0) 
drdx computed on u-rho grid 
drdy computed on v-rho grid
'''


def rho_grad(Tt,Ts,z_r,g,rho0,z_w,pm,pn):
    #subroutine rho_grad(Lm,Mm,N, T,S, z_r,z_w,rho0,pm,pn,
    #  & rho1,qp1,drdz,drdx,drdy)

        (rho1,qp1,drdz,drdx,drdy) = romsF.rho_grad(Tt.T,Ts.T, z_r.T,z_w.T,rho0,pm.T,pn.T)

        return [drdx.T,drdy.T,drdz.T]





#######################################################
#Compute Potential Vorticity using neutral density (rho_grad)
#######################################################
'''

Compute ertel potential vorticity using adiabatic gradients of buoyancy

T and S on horizontal rho grids and vertical rho-grid (specified by z_r)
U and V on horizontal u- and v- grids and vertical rho-grid (specified by z_r)

PV is computed on horizontal psi-grid and vertical w-grid (specified by z_w)

'''

def PV_grad(temp,salt,u,v,z_r,z_w,f,g,rho0,pm,pn,mask=None):


    #rho on rho-rho grid      bvf on rho-w grid
    [dbdx,dbdy,dbdz] = rho_grad(temp,salt,z_r,g,rho0,z_w,pm,pn)


    dz_r = z_r[1:,:,:]- z_r[:-1,:,:]

##########################
#Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*db/dz


    #dudy and dvdx on psi-rho grid
    dvdx = diffxi(v,rho2v(pm),rho2v(z_r),rho2v(z_w),mask)
    dudy = diffeta(u,rho2u(pn),rho2u(z_r),rho2u(z_w),mask)

    dbdz = rho2psi(dbdz)[1:-1,:,:]

    #vrt on psi-rho grid
    vrt = dvdx - dudy


    #PV1 on psi-w grid
    pv1 =  (rho2psi(f) + 0.5*(vrt[1:,:,:] + vrt[:-1,:,:])) * dbdz
    del vrt,dbdz

##########################
#'Ertel potential vorticity, term 2: (dv/dz)*(db/dx)'

    #'dvdz on psi-w grid'
    dvdz = rho2u((v[1:,:,:]-v[:-1,:,:])/(0.5*(dz_r[:,1:,:]+ dz_r[:,:-1,:])))

    #'dbdx on psi-rho grid'
    dbdx = rho2v(dbdx)

    #PV1 on psi-w grid
    pv2 =  -1*dvdz*0.5*(dbdx[1:,:,:] + dbdx[:-1,:,:])
    del dbdx,dvdz

##########################
#'Ertel potential vorticity, term 3: (du/dz)*(db/dy)'

    #'dudz on psi-w grid'
    dudz = rho2v((u[1:,:,:]-u[:-1,:,:])/(0.5*(dz_r[:,:,1:]+ dz_r[:,:,:-1])))

    #'dbdy on psi-rho grid'
    dbdy = rho2u(dbdy)

    #PV3 on psi-w grid
    pv3 =  dudz*0.5*(dbdy[1:,:,:] + dbdy[:-1,:,:])
    del dbdy,dudz

##########################

    return pv1 + pv2 + pv3








#######################################################
#Compute density and Brunt-Vaissala frequency using Fortran subroutine fron alex
#######################################################
'''

compute rho from equation of state

rho on rho (vert and hor) grid

bvf is computed only if z_w is not None
bvf computed on rho-w grid (first and last levels set to 0) 


'''


def rho_eosF_alex(Tt,Ts,z_r,g,rho0,z_w=None):
    #subroutine rho_eos(Lm,Mm,N, T,S, z_r,z_w, rho1,qp1,rho,bvf)

    if z_w is None: 

        z_w=np.zeros((z_r.shape[0]+1,z_r.shape[1],z_r.shape[2]))
        (rho1,qp1,rho,bvf) = romsF.rho_eos(Tt.T,Ts.T, z_r.T,z_w.T,rho0)

        return rho1.T

    else:

        print('computing bvf using Fortran routine')
        print(z_w.shape)
        print('in romstools', rho0)
        (rho1,qp1,rho,bvf) = romsF.rho_eos_alex(Tt.T,Ts.T, z_r.T,z_w.T)

        return [rho1.T,qp1.T,rho.T,bvf.T]

#######################################################
#Get surface wind stress from forcing files interpolated to current time-step
#######################################################

def get_winds(ncname,infiletime,ncnamewind,ny1,ny2,nx1,nx2):


    if isinstance(ncname, str):
        ncfile = Dataset(ncname, 'r', format='NETCDF3_CLASSIC'); opened = True
    else:
        ncfile = ncname; opened = False

    ncfilewind = Dataset(ncnamewind, 'r', format='NETCDF3_CLASSIC')

    oceantime = int(np.array(ncfile.variables['ocean_time'][infiletime]))%(360*24*3600)
    oceanday=oceantime/(24*3600.)


    #if type=='d':
    if ncfilewind.variables['sustr'].shape[0]==360: #daily winds

        datewind1=int(np.floor(oceanday-0.5))%360
        datewind2=int(np.ceil(oceanday-0.5))%360

        if datewind1==datewind2:
            coef1=0.5
            coef2=0.5
        else:
            coef1=abs(oceanday-0.5 - np.ceil(oceanday-0.5))
            coef2=abs(oceanday-0.5 - np.floor(oceanday-0.5))

    elif ncfilewind.variables['sustr'].shape[0]==12: #monthly winds

        datewind1=int(np.floor((oceanday-14.5)/30))%12
        datewind2=int(np.ceil((oceanday-14.5)/30))%12
        
        if datewind1==datewind2:
            coef1=0.5
            coef2=0.5
        else:
            coef1=abs((oceanday-14.5)/30-np.ceil((oceanday-14.5)/30))
            coef2=abs((oceanday-14.5)/30-np.floor((oceanday-14.5)/30))

    else:

        print('sure about your wind forcing file?')


    uwind1=np.array(ncfilewind.variables['sustr'][datewind1,ny1:ny2,nx1:nx2-1]) 
    vwind1=np.array(ncfilewind.variables['svstr'][datewind1,ny1:ny2-1,nx1:nx2]) 

    uwind2=np.array(ncfilewind.variables['sustr'][datewind2,ny1:ny2,nx1:nx2-1]) 
    vwind2=np.array(ncfilewind.variables['svstr'][datewind2,ny1:ny2-1,nx1:nx2]) 

    uwind=coef1*uwind1+coef2*uwind2
    vwind=coef1*vwind1+coef2*vwind2

    if opened==True: ncfile.close()
    ncfilewind.close()

    return [uwind,vwind]

#######################################################
#Get surface wind stress from forcing files interpolated to current time-step
#######################################################

def get_buoy_flux(ncname,infiletime,ncnamewind,ny1,ny2,nx1,nx2,type,Hz,temp,salt):


    if isinstance(ncname, str):
        ncfile = Dataset(ncname, 'r', format='NETCDF3_CLASSIC'); opened = True
    else:
        ncfile = ncname; opened = False

    ncfilewind = Dataset(ncnamewind, 'r', format='NETCDF3_CLASSIC')

    oceantime = int(np.array(ncfile.variables['ocean_time'][infiletime]))%(360*24*3600)
    oceanday=oceantime/(24*3600.)


    if type=='d':

        datewind1=int(np.floor(oceanday-0.5))%360
        datewind2=int(np.ceil(oceanday-0.5))%360

        if datewind1==datewind2:
            coef1=0.5
            coef2=0.5
        else:
            coef1=abs(oceanday-0.5 - np.ceil(oceanday-0.5))
            coef2=abs(oceanday-0.5 - np.floor(oceanday-0.5))

    else:

        datewind1=int(np.floor((oceanday-14.5)/30))%12
        datewind2=int(np.ceil((oceanday-14.5)/30))%12
        
        if datewind1==datewind2:
            coef1=0.5
            coef2=0.5
        else:
            coef1=abs((oceanday-14.5)/30-np.ceil((oceanday-14.5)/30))
            coef2=abs((oceanday-14.5)/30-np.floor((oceanday-14.5)/30))


    #######################################################

    shflx1=np.array(ncfilewind.variables['shflux'][datewind1,ny1:ny2,nx1:nx2]) 
    shflx2=np.array(ncfilewind.variables['shflux'][datewind2,ny1:ny2,nx1:nx2])
    swflx1=np.array(ncfilewind.variables['swflux'][datewind1,ny1:ny2,nx1:nx2]) 
    swflx2=np.array(ncfilewind.variables['swflux'][datewind2,ny1:ny2,nx1:nx2])
    dqdt1=np.array(ncfilewind.variables['dQdSST'][datewind1,ny1:ny2,nx1:nx2]) 
    dqdt2=np.array(ncfilewind.variables['dQdSST'][datewind2,ny1:ny2,nx1:nx2])
    sst1=np.array(ncfilewind.variables['SST'][datewind1,ny1:ny2,nx1:nx2]) 
    sst2=np.array(ncfilewind.variables['SST'][datewind2,ny1:ny2,nx1:nx2])
    sss1=np.array(ncfilewind.variables['SSS'][datewind1,ny1:ny2,nx1:nx2]) 
    sss2=np.array(ncfilewind.variables['SSS'][datewind2,ny1:ny2,nx1:nx2])
    swrad1=np.array(ncfilewind.variables['swrad'][datewind1,ny1:ny2,nx1:nx2]) 
    swrad2=np.array(ncfilewind.variables['swrad'][datewind2,ny1:ny2,nx1:nx2])

    #######################################################

    shflx = coef1*shflx1 + coef2*shflx2
    swflx = coef1*swflx1 + coef2*swflx2
    dqdt = coef1*dqdt1 + coef2*dqdt2
    sst = coef1*sst1 + coef2*sst2
    sss = coef1*sss1 + coef2*sss2
    swrad = coef1*swrad1 + coef2*swrad2




    #######################################################

    rho0=ncfile.rho0
    Cp=3985.
    stflx = shflx/(rho0*Cp) + dqdt/(rho0*Cp)*(temp - sst)

    dsdt = 1./(90.*86400)
    ssflx = swflx*0.01/86400*salt - dsdt*abs(Hz)*(salt - sss)


    if opened==True: ncfile.close()
    ncfilewind.close()

    return [stflx, ssflx]


#######################################################
#Compute thermal expansion and saline contraction coefficients
#######################################################

def alphabeta(Tt,Ts,rho0):

    Q01=6.793952E-2; Q02=-9.095290E-3;
    Q03=+1.001685E-4; Q04=-1.120083E-6; Q05=+6.536332E-9;
    U00=+0.824493; U01=-4.08990E-3; U02=+7.64380E-5 ;
    U03=-8.24670E-7; U04=+5.38750E-9; V00=-5.72466E-3 ;
    V01=+1.02270E-4; V02=-1.65460E-6; W00=+4.8314E-4;


    sqrtTs=Ts ** 0.5;
    cff=1/rho0

    alpha=-cff*( Q01+Tt*( 2.*Q02+Tt*( 3.*Q03+Tt*(4.*Q04 +Tt*5.*Q05 )))\
    +Ts*( U01+Tt*( 2.*U02+Tt*(3.*U03 +Tt*4.*U04 ))+sqrtTs*( V01+Tt*2.*V02)))

    beta= cff*( U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))\
    +1.5*(V00+Tt*(V01+Tt*V02))*sqrtTs+2.*W00*Ts )


    return [alpha, beta]








#######################################################
#z-derivative (messy, needs to be rewritten)
#######################################################

def diffz(var,z):

    if np.rank(var)<3:
        dvardz = diffz_2d(var,z)
    else:
        dvardz = diffz_3d(var,z)

    return dvardz


    #######################

def diffz_2d(var,z):

    [N,M]=var.shape
    dvardz = var*np.nan

    if np.rank(z)==2:
        dvardz[1:,:] = (var[1:,:]-var[:-1,:])/(z[1:,:]-z[:-1,:])
        dvardz[0,:] =  dvardz[1,:]
    else:
        for iz in range(1,N):
            dvardz[iz,:] = (var[iz,:]-var[iz-1,:])/(z[iz]-z[iz-1])
        dvardz[0,:] =  dvardz[1,:]

    return dvardz

    #######################

def diffz_3d(var,z):

    [N,M,L]=var.shape
    dvardz = var*np.nan

    if np.rank(z)==3:
        dvardz[1:,:,:] = (var[1:,:,:]-var[:-1,:,:])/(z[1:,:,:]-z[:-1,:,:])
        dvardz[0,:,:] =  dvardz[1,:,:]
    else:
        for iz in range(1,N):
            dvardz[iz,:,:] = (var[iz,:,:]-var[iz-1,:,:])/(z[iz]-z[iz-1])
        dvardz[0,:,:] =  dvardz[1,:,:]

    return dvardz

#######################################################
#x-derivative from rho-grid to u-grid
#######################################################

def diffx(var,pm):

    if np.rank(var)<3:
        dvardx = diffx_2d(var,pm)
    else:
        dvardx = diffx_3d(var,pm)

    return dvardx

###########################

def diffx_3d(var,pm):

    [N,M,L]=var.shape

    dvardx = np.zeros((N,M,L-1))

    for iz in range(0, N):    
        dvardx[iz,:,:]=diffx_2d(var[iz,:,:],pm)

    return dvardx

###########################

def diffx_2d(var,pm):


    if np.rank(pm)==2: dvardx = (var[:,1:]-var[:,:-1])*0.5*(pm[:,1:]+pm[:,:-1])

    else: dvardx = (var[:,1:]-var[:,:-1])*pm

    return dvardx


#######################################################
#y-derivative from rho-grid to v-grid
#######################################################

def diffy(var,pn):

    if np.rank(var)<3: dvardy = diffy_2d(var,pn)
    else: dvardy = diffy_3d(var,pn)

    return dvardy

    #######################

def diffy_3d(var,pn):

    [N,M,L]=var.shape
    dvardy = np.zeros((N,M-1,L))
    for iz in range(0, N): dvardy[iz,:,:]=diffy_2d(var[iz,:,:],pn)

    return dvardy

    #######################


def diffy_2d(var,pn):

    if np.rank(pn)==2: dvardy = (var[1:,:]-var[:-1,:])*0.5*(pn[1:,:]+pn[:-1,:])
    else: dvardy = (var[1:,:]-var[:-1,:])*pn

    return dvardy

######################################################

def diffysmooth(var,pn):

    if np.rank(var)<3: dvardy = diffy_2dsmooth(var,pn)
    else: dvardy = diffy_3dsmooth(var,pn)

    return dvardy

    #######################

def diffy_3dsmooth(var,pn):

    [N,M,L]=var.shape
    dvardy = np.zeros((N,M,L))

    for iz in range(0, N): dvardy[iz,:,:]=diffy_2dsmooth(var[iz,:,:],pn)

    return dvardy

    #######################

def diffy_2dsmooth(var,pn):

    [M,L]=var.shape
    dvardy = np.zeros((M,L))

    dvardy[1:-1,:] = (var[2:,:]-var[:-2,:])/2*(pn[1:-1,:])
    dvardy[0,:] = dvardy[1,:]
    dvardy[-1,:] = dvardy[-2,:]

    return dvardy





#######################################################
#Compute horizontal derivatives on sigma-levels (1st order)
#######################################################
'''
var on rho-rho grid
dvardxi on psi-rho grid
'''

def diffxi(var,pm,z_r,z_w=None,newz=None,mask=None):


    if z_r.shape[0]<=2:
        dvardxi = diffxi_2d(var,pm,z_r,z_w,newz,mask)
    else:
        dvardxi = diffxi_3d(var,pm,z_r,z_w,newz,mask)

    ##############################################

    return dvardxi

#######################################################
#######################################################


def diffxi_3d(var,pm,z_r,z_w=None,newz=None,mask=None):


    if newz is None: newz = (z_r[:,:,:-1] + z_r[:,:,1:])/2
    else: newz = rho2u(newz)

    dvardxi = np.zeros((var.shape[0],var.shape[1],var.shape[2]-1))

    ##############################################

    varzp = vinterpsF(var[:,:,1:],newz,z_r[:,:,1:],z_w[:,:,1:])
    varzm = vinterpsF(var[:,:,:-1],newz,z_r[:,:,:-1],z_w[:,:,:-1])

    dvardxi = (varzp - varzm )*0.5*(pm[:,:-1]+pm[:,1:])  

    ##############################################

    return dvardxi


#######################################################
#######################################################


def diffxi_2d(var,pm,z_r,z_w=None,newz=None,mask=None):

    dvardxi = np.zeros((z_r.shape[1],z_r.shape[2]-1))

    ##############################################

    if newz is None: newz = (z_r[0,:,:-1] + z_r[0,:,1:])/2
    else: newz = rho2u(newz)

    dz0 = (z_r[0,:,1:]-newz)
    dz1 = (newz-z_r[1,:,1:])
    varzp = (dz1*var[0,:,1:] + dz0*var[1,:,1:])/(z_r[0,:,1:]-z_r[1,:,1:])

    dz0 = (z_r[0,:,:-1]-newz)
    dz1 = (newz-z_r[1,:,:-1])
    varzm = (dz1*var[0,:,:-1] + dz0*var[1,:,:-1])/(z_r[0,:,:-1]-z_r[1,:,:-1])

    dvardxi = (varzp - varzm )*0.5*(pm[:,:-1]+pm[:,1:])
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


    if z_r.shape[0]<=2:
        dvardeta = diffeta_2d(var,pn,z_r,z_w,newz,mask)
    else:
        dvardeta = diffeta_3d(var,pn,z_r,z_w,newz,mask)

    ##############################################

    return dvardeta


#######################################################
#######################################################


def diffeta_3d(var,pn,z_r,z_w=None,newz=None,mask=None):


    if newz is None: newz = (z_r[:,:-1,:] + z_r[:,1:,:])/2
    else: newz = rho2v(newz)

    dvardeta = np.zeros((var.shape[0],var.shape[1]-1,var.shape[2]))

    ##############################################

    varzp = vinterpsF(var[:,1:,:],newz,z_r[:,1:,:],z_w[:,1:,:])
    varzm = vinterpsF(var[:,:-1,:],newz,z_r[:,:-1,:],z_w[:,:-1,:])


    dvardeta = (varzp - varzm )*0.5*(pn[:-1,:]+pn[1:,:])

    ##############################################


    return dvardeta



#######################################################
#Compute horizontal derivatives on sigma-levels (1st order)
#######################################################



def diffeta_2d(var,pn,z_r,z_w=None,newz=None,mask=None):

    dvardeta = np.zeros((z_r.shape[1]-1,z_r.shape[2]))

    ##############################################

    if newz is None: newz = (z_r[0,:-1,:] + z_r[0,1:,:])/2
    else: newz = rho2v(newz)

    dz0 = (z_r[0,1:,:]-newz)
    dz1 = (newz-z_r[1,1:,:])
    varzp = (dz1*var[0,1:,:] + dz0*var[1,1:,:])/(z_r[0,1:,:]-z_r[1,1:,:])

    dz0 = (z_r[0,:-1,:]-newz)
    dz1 = (newz-z_r[1,:-1,:])
    varzm = (dz1*var[0,:-1,:] + dz0*var[1,:-1,:])/(z_r[0,:-1,:]-z_r[1,:-1,:])

    dvardeta = (varzp - varzm )*0.5*(pn[:-1,:]+pn[1:,:])

    ##############################################


    return dvardeta





#######################################################
#operator grad.
#######################################################

def grad(var,pm=1.,pn=1.):


    amplitude = u2rho(diffx(var,pm)**2) + v2rho(diffy(var,pn)**2)

    ##############################################


    return amplitude






#######################################################
#Compute vertical velocity (see Wvlcty.F)
#######################################################

def getw(u,v,pm,pn,z_r,z_w):

    
    #rho grid = Nrho,Mrho,Lrho
    #u grid = Nrho,Mrho,Lu
    #v grid = Nrho,Mv,Lrho
    #w grid = Nw,Mrho,Lrho

    [Mrho,Lrho]=pm.shape
    Mv=Mrho-1
    Lu=Lrho-1

    if len(u.shape)==3: Nrho=u.shape[0]
    else: Nrho = 1

    Nw=Nrho+1

    ###########################

    flxu = (rho2u(z_w[1:,:,:]-z_w[:-1,:,:]))/(0.5*(pn[:,1:]+pn[:,:-1]))*u
    flxv = (rho2v(z_w[1:,:,:]-z_w[:-1,:,:]))/(0.5*(pm[1:,:]+pm[:-1,:]))*v

    ###########################

    wrk = np.zeros((Nw,Mrho,Lrho))*np.nan
    wvlc = np.zeros((Nw,Mrho,Lrho))*np.nan

    wrk[0,1:-1,1:-1]=0
    wvlc[0,1:-1,1:-1]=0
    for iz in range(1,Nw):
        wrk[iz,1:-1,1:-1] = wrk[iz-1,1:-1,1:-1] - pm[1:-1,1:-1] * pn[1:-1,1:-1] * (flxu[iz-1,1:-1,1:] - flxu[iz-1,1:-1,:-1] + flxv[iz-1,1:,1:-1] - flxv[iz-1,:-1,1:-1])

    if Nw==2:
        wvlc[1,:,:] = wrk[1,:,:]

    else:
        #move to vertical rho points
        wvlc[1,:,:] = -0.125*wrk[2,:,:] + 0.75*wrk[1,:,:] + 0.375 * wrk[0,:,:]
        wvlc[-1,:,:] = 0.375*wrk[-1,:,:] + 0.75*wrk[-2,:,:] - 0.125 * wrk[-3,:,:]
        for iz in range(2,Nw-1):
            wvlc[iz,:,:] = 0.5625*(wrk[iz-1,:,:]+wrk[iz,:,:])-0.0625*(wrk[iz-2,:,:]+wrk[iz+1,:,:])


    #add contributions due to S-coord slopes (u*dz/dx and v*dz/dy)
    Wx = u*(z_r[:,:,1:]-z_r[:,:,:-1])*(pm[:,1:]+pm[:,:-1])
    Wy = v*(z_r[:,1:,:]-z_r[:,:-1,:])*(pn[1:,:]+pn[:-1,:])


    wvlc[1:,1:-1,1:-1] = wvlc[1:,1:-1,1:-1] + 0.25 * (Wy[:,1:,1:-1]+Wy[:,:-1,1:-1]+Wx[:,1:-1,1:]+Wx[:,1:-1,:-1])

    return [wvlc]



   
#######################################################
#Compute mean
#######################################################


def nanmean(data, *args):

    dataout = np.ma.filled(np.ma.masked_array(data,np.isnan(data)).mean(*args), fill_value=np.nan)
    if dataout.shape==(): dataout = float(dataout)

    return dataout



#######################################################
#Lagrange Interpolation
#######################################################


def lag_intrp(x1,x2,x3,x):

    a1 = (x-x2)*(x-x3)/((x1-x2)*(x1-x3))
    a2 = (x-x1)*(x-x3)/((x2-x1)*(x2-x3))
    a3 = (x-x1)*(x-x2)/((x3-x1)*(x3-x2))

    return [a1,a2,a3]





#######################################################
#Rotate winds or u,v to lat,lon coord
#######################################################

def rotuv(gname,u,v,ny1=None,ny2=None,nx1=None,nx2=None):

    if isinstance(gname, str):
        print('are you sure you want to use angle?')
        nc = Dataset(gname, 'r', format='NETCDF3_CLASSIC')
        angle = np.array(nc.variables['angle'][ny1:ny2,nx1:nx2])
        nc.close()
    else:
        angle = gname 

    if u.shape!=v.shape:
        u=u2rho(u)
        v=v2rho(v)

    'rotate vectors by geometric angle'
    urot = u*np.cos(angle) - v*np.sin(angle)
    vrot = u*np.sin(angle) + v*np.cos(angle)

    return [urot,vrot]

#######################################################
#Rotation
#######################################################

#the function
def rotate_2d(pts,cnt,ang=np.pi/4):
    '''pts = {} Rotates points(nx2) about center cnt(2) by angle ang(1) in radian'''
    return np.dot(pts-cnt,np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]]))+cnt


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
        ux=u2rho( u2rho( diffx(u,pm)))
        vy=v2rho( v2rho( diffy(v,pn)))
        uy=psi2rho( diffy(u,pn))
        vx=psi2rho( diffx(v,pm))
    s1 = ux-vy; s2= vx+uy; 
    thetas = np.arctan(s2/s1)/2; thetaps = np.arctan(-1*s1/s2)/2;
    #see if division by 0
    eps = 1e-15; 
    thetas[np.abs(s1)<eps] = np.sign(s2[np.abs(s1)<eps])*np.pi/4
    #check if s1'>0 (s1<0 means that you are on the perpendicular axis)
    s1bis = s1 * np.cos(2*thetas) + s2*np.sin(2*thetas)
    thetas[s1bis<0] = thetas[s1bis<0]+np.pi/2
    s1bis = s1 * np.cos(2*thetas) + s2*np.sin(2*thetas)

    s2bis = -s1 * np.sin(2*thetas) + s2*np.cos(2*thetas)

    return thetas,s1bis,s2bis






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
            ux=u2rho( u2rho( diffx(u,rho2u(pm))))
            vy=v2rho( v2rho( diffy(v,rho2v(pn))))
            uy=psi2rho( diffy(u,rho2u(pn)))
            vx=psi2rho( diffx(v,rho2v(pm)))
        
    s1 = ux-vy; s2= vx+uy; div=ux+vy
    #thetas = np.arctan(s2/s1)/2; 
    thetas = np.arctan(-1*s1/s2)/2;
    #see if division by 0
    eps = 1e-15; 
    thetas[np.abs(s2)<eps] = 0.
    #check if s2'>0 (s2'>0 means that you are on the perpendicular axis)
    s2bis = -s1/2. * np.sin(2*thetas) + s2*np.cos(2*thetas)
    thetas[s2bis>0] = thetas[s2bis>0]+np.pi/2
    s2bis = -s1/2. * np.sin(2*thetas) + s2*np.cos(2*thetas)

    s1bis = s1 * np.cos(2*thetas) + s2*np.sin(2*thetas)

    return thetas,s2bis,s1bis
    


















# End of romstools.py
