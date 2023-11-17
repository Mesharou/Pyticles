#!/usr/bin/python
# Filename: romstools.py


#Netcdf IO module
from netCDF4 import Dataset

#module for numerics
import numpy as np

#module for copy
from copy import copy

####################
from itertools import product



####################

def zlevs(h,zeta, hc, Cs_r, Cs_w):
    
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

#######################################################
#interpolate a 3D variable on a horizontal level of constant depth
#######################################################

def vinterp(var,z,depth,topo=None,cubic=0):

    [N,Mp,Lp]=z.shape

    #depth=0, just take surface field
    if np.all(depth>0):
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
        #if np.all(topo!=None): varz[depth<-1*topo]=np.nan

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

        varz[depth<-1*topo]=np.nan

   
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


    [Nz,Mz,Lz]=var_psi.shape
    var_rho=np.zeros((Nz,Mz+1,Lz+1))

    for iz in range(0, Nz, 1):    
        var_rho[iz,:,:]=psi2rho_2d(var_psi[iz,:,:])


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

    var_psi = 0.25*(var_rho[:,1:,1:]+var_rho[:,1:,:-1]+var_rho[:,:-1,:-1]+var_rho[:,:-1,1:])

    return var_psi


#######################################################
#Transfert a 2 or 3-D field at rho points to u points
#######################################################

def rho2u(var_rho):

    if np.ndim(var_rho)<3:
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

    if np.ndim(var_rho)<3:
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


    if np.ndim(var_u)<3:
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

    if np.ndim(var_v)<3:
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

    if np.ndim(u)<3:
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

    if np.ndim(u)>2:

        [Nz,Mz,Lz]=u.shape
        [Mp,Lp]=pm.shape

        divs = np.zeros((Nz,Mp,Lp))*np.nan

        for iz in range(0, Nz, 1):
            divs[iz,:,:]=div(u[iz,:,:],v[iz,:,:],pm,pn)

    else:
    
        divs=div(u,v,pm,pn)

    return divs




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

    if np.ndim(Tt)==2:
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

    if z_w!=None:

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

    if np.ndim(var)<3:
        dvardz = diffz_2d(var,z)
    else:
        dvardz = diffz_3d(var,z)

    return dvardz


    #######################

def diffz_2d(var,z):

    [N,M]=var.shape
    dvardz = var*np.nan

    if np.ndim(z)==2:
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

    if np.ndim(z)==3:
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

    if np.ndim(var)<3:
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


    if np.ndim(pm)==2: dvardx = (var[:,1:]-var[:,:-1])*0.5*(pm[:,1:]+pm[:,:-1])

    else: dvardx = (var[:,1:]-var[:,:-1])*pm

    return dvardx


#######################################################
#y-derivative from rho-grid to v-grid
#######################################################

def diffy(var,pn):

    if np.ndim(var)<3: dvardy = diffy_2d(var,pn)
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

    if np.ndim(pn)==2: dvardy = (var[1:,:]-var[:-1,:])*0.5*(pn[1:,:]+pn[:-1,:])
    else: dvardy = (var[1:,:]-var[:-1,:])*pn

    return dvardy

######################################################

def diffysmooth(var,pn):

    if np.ndim(var)<3: dvardy = diffy_2dsmooth(var,pn)
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


    if newz==None: newz = (z_r[:,:,:-1] + z_r[:,:,1:])/2
    else: newz = rho2u(newz)

    dvardxi = np.zeros((var.shape[0],var.shape[1],var.shape[2]-1))

    ##############################################

    varzp = vinterps(var[:,:,1:],newz,z_r[:,:,1:],z_w[:,:,1:])
    varzm = vinterps(var[:,:,:-1],newz,z_r[:,:,:-1],z_w[:,:,:-1])

    dvardxi = (varzp - varzm )*0.5*(pm[:,:-1]+pm[:,1:])  

    ##############################################

    return dvardxi


#######################################################
#######################################################


def diffxi_2d(var,pm,z_r,z_w=None,newz=None,mask=None):

    dvardxi = np.zeros((z_r.shape[1],z_r.shape[2]-1))

    ##############################################

    if newz==None: newz = (z_r[0,:,:-1] + z_r[0,:,1:])/2
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


    if newz==None: newz = (z_r[:,:-1,:] + z_r[:,1:,:])/2
    else: newz = rho2v(newz)

    dvardeta = np.zeros((var.shape[0],var.shape[1]-1,var.shape[2]))

    ##############################################

    varzp = vinterps(var[:,1:,:],newz,z_r[:,1:,:],z_w[:,1:,:])
    varzm = vinterps(var[:,:-1,:],newz,z_r[:,:-1,:],z_w[:,:-1,:])


    dvardeta = (varzp - varzm )*0.5*(pn[:-1,:]+pn[1:,:])

    ##############################################


    return dvardeta



#######################################################
#Compute horizontal derivatives on sigma-levels (1st order)
#######################################################



def diffeta_2d(var,pn,z_r,z_w=None,newz=None,mask=None):

    dvardeta = np.zeros((z_r.shape[1]-1,z_r.shape[2]))

    ##############################################

    if newz==None: newz = (z_r[0,:-1,:] + z_r[0,1:,:])/2
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














# End of romstools.py
