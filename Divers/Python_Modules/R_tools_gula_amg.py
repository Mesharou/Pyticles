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
from R_tools import *
import R_tools_fort as toolsF
import R_tools_fort_gula as toolsF_g

#Simulations (path, data...)
import R_vars_gula as va

#for plotting
import matplotlib.pyplot as py

import time as tm

import R_smooth as sm







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


    (Jbot) = toolsF_g.get_jbot_sol1(T,S,u,v,z_r,z_w,rho0,pm,pn,hbbls,rdrg)

    return Jbot


#################################################


def get_jbot_sol2(T,S,u,v,z_r,z_w,rho0,pm,pn,hbbls,rdrg):


    (Jbot) = toolsF_g.get_jbot_sol2(T,S,u,v,z_r,z_w,rho0,pm,pn,hbbls,rdrg)

    return Jbot















    
#######################################################
#Strain (direction + amplitude)
#######################################################


def straindir(u,v,pm=1.,pn=1.):
    
    if u.shape[0]==v.shape[0]: 
        ux=u2rho( diffx(u,pm))
        vy=v2rho( diffy(v,pn))
        uy=v2rho( diffy(u,pn))
        vx=u2rho( diffx(v,pm))   
        
    else:    
        if isinstance(pm,float):
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

def div2uvs(u,v,pm,pn,verbo=False,variable=False,fast=False):
    
    if len(u.shape)>2:
        
        udiv = np.zeros(u.shape)*np.nan
        vdiv = np.zeros(v.shape)*np.nan
        
        for iz in range(u.shape[2]):
            if variable:
                udiv[:,:,iz],vdiv[:,:,iz] = div2uv_variable(u[:,:,iz],v[:,:,iz],pm,pn,verbo=verbo,fast=fast)
            else:
                udiv[:,:,iz],vdiv[:,:,iz] = div2uv(u[:,:,iz],v[:,:,iz],pm,pn,verbo=verbo)           
    else:
        if variable:    
            udiv,vdiv = div2uv_variable(u,v,pm,pn,verbo=verbo,fast=fast)
        else:
            udiv,vdiv = div2uv(u,v,pm,pn,verbo=verbo)
   
    return udiv,vdiv

##################

# if pm,pn is constant

def div2uv(u,v,pm,pn,verbo=False):

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
    if verbo: print(ml)                                 # print hierarchy information
    b = -1*div.flatten()*1/(np.mean(pm)*np.mean(pn))  # right hand side
    x = ml.solve(b, tol=1e-10)               # solve Ax=b to a tolerance of 1e-8
    if verbo: print("residual norm is", norm(b - A*x))  # compute norm of residual vector
    
    udiv = diffx(x.reshape(div.shape),pm)
    vdiv = diffy(x.reshape(div.shape),pn)
    
    return udiv,vdiv


##################
# test with variable pm,pn

def div2uv_variable(u,v,pm,pn,verbo=False,fast=False):

    # compute div
    div = np.zeros(pm.shape)
    div[1:-1,:] = div[1:-1,:] + diffx(u,rho2u(pm))
    div[:,1:-1] = div[:,1:-1] + diffy(v,rho2v(pn))
    div[isnan(div)] =0
    
    # solve poisson
    #A = poisson(div.shape, format='csr')     # 2D Poisson problem 
    #######################################################
    #Create matrix A
    #######################################################    
    print('creating matrix A')
    if fast:
        A = poisson_matrix_fast(pm,pn)
    else:
        A = poisson_matrix(pm,pn)    
    #######################################################
    #Solve matrix A
    A = A.tocsr()
    #######################################################

    ml =ruge_stuben_solver(A)                # construct the multigrid hierarchy
    if verbo: print(ml)                                 # print hierarchy information
    b = -1*div.flatten() # right hand side
    x = ml.solve(b, tol=1e-10)               # solve Ax=b to a tolerance of 1e-8
    if verbo: print("residual norm is", norm(b - A*x))  # compute norm of residual vector
    
    udiv = diffx(x.reshape(div.shape),pm)
    vdiv = diffy(x.reshape(div.shape),pn)
    
    return udiv,vdiv

##################
# test with variable pm,pn

def div2uv_noamg(u,v,pm,pn,verbo=False):

    print('using  div2uv_noamg (debug)')
    # compute div
    div = np.zeros(pm.shape)
    div[1:-1,:] = div[1:-1,:] + diffx(u,rho2u(pm))
    div[:,1:-1] = div[:,1:-1] + diffy(v,rho2v(pn))
    div[isnan(div)] =0

    # solve poisson
    #A = poisson(div.shape, format='csr')     # 2D Poisson problem
    #######################################################
    #Create matrix A
    #######################################################
    print('creating matrix A')
    #A = poisson_matrix(pm,pn)
    A = poisson_matrix_fast(pm,pn)
    #######################################################
    #Solve matrix A
    A = A.tocsr()
    #######################################################

    #ml =ruge_stuben_solver(A)                # construct the multigrid hierarchy
    if verbo: print(ml)                                 # print hierarchy information
    b = -1*div.flatten() # right hand side
    #x = ml.solve(b, tol=1e-10)               # solve Ax=b to a tolerance of 1e-8
    x = spsolve(A,b)
    if verbo: print("residual norm is", norm(b - A*x))  # compute norm of residual vector

    udiv = diffx(x.reshape(div.shape),pm)
    vdiv = diffy(x.reshape(div.shape),pn)

    return udiv,vdiv



#######################################################
# Compute solution of omega equation
# 
#######################################################
'''


'''
from pyamg import *
from pyamg.gallery import *
from scipy import *
from scipy.linalg import *  
from scipy.sparse import *
from scipy.sparse.linalg import spsolve
    
def solve_omega(buoy,pm,pn,f,N2,depths,ur=None,vr=None,nh=1,forcing=0.,mixrotuv=True,bbc=0,wbot=None):
    
    print('in solve_omega mixrotuv is', mixrotuv)

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
            print(' ')
            print('using velocity and buoyancy field to define Q vector')
            print(' ')            
            uz = -(by.T/f.T).T; vz =  (bx.T/f.T).T;
        
        else:
            #using only non-divergent velocity field
            print(' ')
            print('using only velocity field to define Q vector')
            print(' ')
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
    #Qx = 2*(f.T*(vx*uz + vy*vz).T).T;
    #Qy =-2*(f.T*(ux*uz + uy*vz).T).T;
    
    Qx = 2*(f.T*(vx*uz - ux*vz).T).T;
    Qy = 2*(f.T*(vy*uz - uy*vz).T).T;
    
    Qxx, Qyy =  copy(new),copy(new)
    Qxx[1:-1,:,:]= diffx(Qx,pm,2); Qyy[:,1:-1,:] = diffy(Qy,pn,2)
    #if periodicity in y (used for jet example)    
    #Qyy[:,0,:] = (Qy[:,1,:]-Qy[:,-1,:])*0.5*(pn[:,1]+pn[:,-1])/2
    #Qyy[:,-1,:] = (Qy[:,0,:]-Qy[:,-2,:])*0.5*(pn[:,1]+pn[:,-1])/2
    
    #print 'N2.shape', N2.shape
    divQ = (Qxx + Qyy) #/N2
    
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
    
    #######################################################
    if bbc==2:
      for i in range(nx):
        for j in range(ny): 
          k=0
          idx = i*ny*(nz-1) + k*ny + j;
          R[idx] = wbot[i,j]

    #######################################################

    R[np.isnan(R)]=0

    #######################################################
    #Create matrix A
    #######################################################    
    
    print('creating matrix A')
    if bbc==2:
        A = omega_matrix_bottom(pm,pn,depths,f,N2,bbc=bbc)
    else:
        A = omega_matrix(pm,pn,depths,f,N2,bbc=bbc)
    
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

def omega_matrix(pm,pn,depths,f,N2,bbc=0):
    
    #bbc = bottom boundary condition.
    # bbc = 0 <> w = 0
    # bbc = 1 <> dw/dz = 0
    
    # elliptic equation matrix divided by N^2: (f/N)^2 d_zz + d_xx + d_yy
    dx =1/np.mean(pm)
    dy =1/np.mean(pn)   
    #dz = depths[1]-depths[0]
    
    nz=len(depths)
    [nx,ny] = pm.shape

    ############################
    
    print('bottom boundary condition is ', bbc)
    
    dx2i = 1./(dx*dx);
    dy2i = 1./(dy*dy);
    #dz2i = 1./(dz*dz);


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
                
                f2N2 = f[i,j]**2
                    
                if len(pm.shape)==2: 
                    dx2i= pm[i,j]**2
                    dy2i= pn[i,j]**2
      
                
                idx = i*ny*(nz-1) + k*ny + j;
                diag = 0.;

                if j>0:
                    A[idx,idx-js] = dy2i*N2[k];
                    diag = diag - dy2i*N2[k];
                #else:
                    #A[idx,idx+bjs] = dy2i;
                    #diag = diag - dy2i;

                if k>0:
                    dz2m = 1./((depths[k]-depths[k-1])*0.5*(depths[k+1]- depths[k-1]))
                    A[idx,idx-ks] = f2N2*dz2m;
                    diag = diag - f2N2*dz2m;
                elif bbc==0:
                    dz2m = 1./((depths[1]-depths[0])**2)
                    diag = diag - f2N2*dz2m;

                if i>0:
                    A[idx,idx-i_s] = dx2i*N2[k];
                    diag = diag - dx2i*N2[k];

                if i<nx-1:
                    A[idx,idx+i_s] = dx2i*N2[k];
                    diag = diag - dx2i*N2[k];

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
                    A[idx,idx+js] = dy2i*N2[k];
                    diag = diag - dy2i*N2[k];
                #else:
                    #A[idx,idx-bjs] = dy2i;
                    #diag = diag - dy2i;

                A[idx,idx] = diag;

                

    return A
    

######################################################  

def omega_matrix_bottom(pm,pn,depths,f,N2,bbc=0):
    
    #bbc = bottom boundary condition.
    # bbc = 0 <> w = 0
    # bbc = 1 <> dw/dz = 0
    
    # elliptic equation matrix divided by N^2: (f/N)^2 d_zz + d_xx + d_yy
    dx =1/np.mean(pm)
    dy =1/np.mean(pn)   
    #dz = depths[1]-depths[0]
    
    nz=len(depths)
    [nx,ny] = pm.shape

    ############################
    
    print('bottom boundary condition is ', bbc)
    
    dx2i = 1./(dx*dx);
    dy2i = 1./(dy*dy);
    #dz2i = 1./(dz*dz);


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
                
                f2N2 = f[i,j]**2
                    
                if len(pm.shape)==2: 
                    dx2i= pm[i,j]**2
                    dy2i= pn[i,j]**2
      
                
                idx = i*ny*(nz-1) + k*ny + j;
                diag = 0.;
                
                if k==0:
                    diag = 1.
                else:
                    if j>0:
                        A[idx,idx-js] = dy2i*N2[k];
                        diag = diag - dy2i*N2[k];
                    #else:
                        #A[idx,idx+bjs] = dy2i;
                        #diag = diag - dy2i;

                    if k>0:
                        dz2m = 1./((depths[k]-depths[k-1])*0.5*(depths[k+1]- depths[k-1]))
                        A[idx,idx-ks] = f2N2*dz2m;
                        diag = diag - f2N2*dz2m;
                    elif bbc==0:
                        dz2m = 1./((depths[1]-depths[0])**2)
                        diag = diag - f2N2*dz2m;

                    if i>0:
                        A[idx,idx-i_s] = dx2i*N2[k];
                        diag = diag - dx2i*N2[k];

                    if i<nx-1:
                        A[idx,idx+i_s] = dx2i*N2[k];
                        diag = diag - dx2i*N2[k];

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
                        A[idx,idx+js] = dy2i*N2[k];
                        diag = diag - dy2i*N2[k];
                    #else:
                        #A[idx,idx-bjs] = dy2i;
                        #diag = diag - dy2i;

                A[idx,idx] = diag;

                

    return A
    

  
#######################################################
# Compute solution of omega equation
# 
#######################################################
'''


'''
from pyamg import *
from pyamg.gallery import *
from scipy import *
from scipy.linalg import *  
from scipy.sparse import *
from scipy.sparse.linalg import spsolve
    
def solve_genomega(buoy, pm, pn, f, N2, depths, ur=None, vr=None, u=None, v=None, nh=1, forcing=0., terms=[0], terms_data=[0], debug=False, smoothing=False, nsmooth=2):
    
    
    
    #######################################################
    #Create forcing (Q vector divergence)
    #######################################################

    if ur==None: ur,vr = u,v
    
    new = np.zeros((pm.shape[0],pm.shape[1],len(depths)))

    nz=len(depths); [nx,ny] = pm.shape; ndim = (nz-1)*ny*nx;
    print('number of points is ',nx,ny,nz)
    dz = depths[1:]-depths[:-1]

    if terms['ROMS']==0 or debug:
        
        #get gradients
        bx,by = copy(new),copy(new)
        bx[1:-1,:,:]= diffx(buoy,pm,2); by[:,1:-1,:] = diffy(buoy,pn,2)  
        
        #if periodicity in y (used for jet example)
        #by[:,0,:] = (buoy[:,1,:]-buoy[:,-1,:])*0.5*(pn[:,1]+pn[:,-1])/2
        #by[:,-1,:] = (buoy[:,0,:]-buoy[:,-2,:])*0.5*(pn[:,1]+pn[:,-1])/2   
        
        ux,uy,uz = copy(new),copy(new),copy(new)
        vx,vy,vz = copy(new),copy(new),copy(new)
    
    
        if terms['VECQ']==1:   

                
            ux[1:-1,:,:] = diffx(ur,pm,2); uy[:,1:-1,:] = diffy(ur,pn,2)   
            vx[1:-1,:,:] = diffx(vr,pm,2); vy[:,1:-1,:] = diffy(vr,pn,2)       
                        
            print(' ')
            print('using velocity and buoyancy field to define Q vector')
            print(' ')
            
            uz = -(by.T/f.T).T; vz =  (bx.T/f.T).T;



            #################################################################
                
            #Components of Q vector = (Qx,Qy)
            #Qx = 2*(f.T*(vx*uz + vy*vz).T).T;
            #Qy =-2*(f.T*(ux*uz + uy*vz).T).T;
            
            Qx = 2*(f.T*(vx*uz - ux*vz).T).T;
            Qy = 2*(f.T*(vy*uz - uy*vz).T).T;       

            Qxx, Qyy =  copy(new),copy(new)
            Qxx[1:-1,:,:]= diffx(Qx,pm,2); Qyy[:,1:-1,:] = diffy(Qy,pn,2)
            #if periodicity in y (used for jet example)    
            #Qyy[:,0,:] = (Qy[:,1,:]-Qy[:,-1,:])*0.5*(pn[:,1]+pn[:,-1])/2
            #Qyy[:,-1,:] = (Qy[:,0,:]-Qy[:,-2,:])*0.5*(pn[:,1]+pn[:,-1])/2
            
            #print 'N2.shape', N2.shape
            divQ = (Qxx + Qyy) + forcing #/N2

            
            if debug: terms_data['VECQ'] = (Qxx + Qyy)
            
            #################################################################
            
            
            
        #################################################################
        # Horizontal vorticity of the flow

        
        ux[1:-1,:,:] = diffx(u,pm,2); uy[:,1:-1,:] = diffy(u,pn,2)   
        vx[1:-1,:,:] = diffx(v,pm,2); vy[:,1:-1,:] = diffy(v,pn,2)     
        
        if terms['VA']==1 or terms['LHS1']==1 or terms['DVAD']==1:
            vrt = (vx - uy)
            vrt[np.isnan(vrt)]=0
        else:
            vrt=None
            
        #################################################################
 
            
        
        if terms['LTAD']==1 or terms['DVAD']==1:
            ''' 
            alternate version of the r.h.s
            
            '''

            print(' ')       
            print('using alternate version of the r.h.s')
            print('terms[LTAD], terms[DVAD]', terms['LTAD'], terms['DVAD'])
            print(' ')   
  
            
            #################################################################
            if terms['DVAD']==1:
                
                absvrt = (vrt.T + f.T).T
                vrtadv = copy(new)        
                vrtadv[1:-1,:,:] = u[1:-1,:,:] * diffx(absvrt,pm,2) 
                vrtadv[:,1:-1,:] = vrtadv[:,1:-1,:]  + v[:,1:-1,:] * diffy(absvrt,pn,2)
                
                del absvrt
                vrtadv = diffz(vrtadv,depths)
            
            #################################################################
            
            divQ = forcing
            
            if terms['LTAD']==1: 
                divQ = divQ - laplacien(u*bx+v*by,pm,pn)
                if debug: terms_data['LTAD'] = - laplacien(u*bx+v*by,pm,pn)
                
            if terms['DVAD']==1: 
                divQ = divQ + ( vrtadv.T * f.T ).T 
                if debug: terms_data['DVAD'] = ( vrtadv.T * f.T ).T 


            del vrtadv, forcing
            
            #################################################################
            
        if terms['LTAD']==0 and terms['DVAD']==0 and terms['VECQ']==0:   
        
            divQ = forcing
            
            
            
            #################################################################            

        if terms['T1']==1:
            
            print(' ')                  
            print('adding T1 term to r.h.s')
            print(' ')                
            uz = diffz(u,depths); vz = diffz(v,depths)
            
            divx,divy = copy(new),copy(new)
            divx[1:-1,:,:] = diffx(ux + vy,pm,2)
            divy[:,1:-1,:] = diffy(ux + vy,pn,2)

            T1 = divx*vz - divy*uz
            
            divQ = divQ - ( T1.T * f.T ).T
            
            if debug: terms_data['T1'] = - ( T1.T * f.T ).T
            
            #################################################################

                    
    else:
        
        divQ = forcing

        vrt=None    
            
    #################################################################
    #################################################################


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
        
        
        
    if smoothing:  
        print('Smoothing divQ using nsmooth=', nsmooth)
        for iz in range(len(depths)):
            divQ[:,:,iz] =  sm.smooth_2d(divQ[:,:,iz],nsmooth)
            N2[:,:,iz] = sm.smooth_2d(N2[:,:,iz],nsmooth)
            
            if terms['VA']==1 or terms['LHS1']==1 or terms['DVAD']==1:
                vrt[:,:,iz] = sm.smooth_2d(vrt[:,:,iz],nsmooth)
                


    #######################################################
    # reorder forcings from (i,j,k) to vector
    R = np.zeros(ndim);
    for i in range(nx): 
        for j in range(ny): 
            for k in range(nz-1):
                idx = i*ny*(nz-1) + k*ny + j;
                R[idx] = divQ[i,j,k]

    R[np.isnan(R)]=0
    
    
    
    
    
    if not debug:
        
        #######################################################
        #Create matrix A
        #######################################################    
        print('creating matrix A')
        A = genomega_matrix(pm,pn,depths,f,N2,vrt,u,v,terms=terms)
        
        #######################################################
        #Solve matrix A
        A = A.tocsr()
        #print 'solving equation'
        tstart = tm.time() 
        
        # Method 1 
        try:
            ml =ruge_stuben_solver(A)                # construct the multigrid hierarchy
            print(ml)                                 # print hierarchy information
            X = ml.solve(R, tol=1e-8)               # solve Ax=b to a tolerance of 1e-8
            print("residual norm is", norm(R - A*X))  # c
            print('Using ruge_stuben_solver.........', tm.time()-tstart)

        
        # Method 2 
        except:
            print('   ')     
            print('no convergence with the ruge_stuben_solver')
            #print '   '     
            #X = spsolve(A,R)
            #print "residual norm is", norm(R - A*X)  # c
            #print 'Using spsolve.........', tm.time()-tstart
            #tstart = tm.time()  
            raise
        
        #######################################################  
        # reorder results in (i,j,k)
        
        w = np.zeros((nx,ny,nz))
        for i in range(nx): 
            for j in range(ny): 
                for k in range(nz-1):
                    idx = i*ny*(nz-1) + k*ny + j; 
                    w[i,j,k] = X[idx];

        #######################################################  
        

    if debug:
        w = np.zeros((nx,ny,nz))
        return w, divQ, terms_data
    else:
        return w
        
    
####################################################### 

######################################################  

def genomega_matrix(pm,pn,depths,f,N2,vrt,u,v,terms):    
    '''
    Differences with omega_matrix:
    
    add horizontal variations of N2
    add vorticity of the flow
    
    
    '''
    
    
    # elliptic equation matrix divided by N^2: (f/N)^2 d_zz + d_xx + d_yy
    
    nz=len(depths)
    [nx,ny] = pm.shape

    ############################


    ndim = (nz-1)*ny*nx;
    i_s  = ny*(nz-1);
    js  = 1;
    bjs = ny-1;
    ks  = ny;
    #A = np.zeros((ndim,ndim));
    #A=csc_matrix((ndim,ndim))
    A=lil_matrix((ndim,ndim))
    
    ############################

    if terms['VA']==1: 
        print('VA term is added to the matrix')
        vrtzz = diffzz(vrt,depths)

    if terms['T2']==1:         
        print('T2 term is added to the matrix')
        uzz = diffzz(u,depths)
        vzz = diffzz(v,depths)
        
    if terms['LHS2']==0 and terms['LHS2QG']==0: LHS2 = False
    else: LHS2 = True
    ############################
    
    for i in range(nx): 
    
        for j in range(ny): 
        
            for k in range(nz-1):
                

                f2N2 = f[i,j]**2
                
                if terms['LHS1']==1: f2N2 = np.nanmax([f2N2 + f[i,j]*vrt[i,j,k],0])
                if terms['LHS1']==0 and terms['LHS1QG']==0: f2N2 = 0.

                
                #if len(pm.shape)==2: 
                    #dx2i= pm[i,j]**2
                    #dy2i= pn[i,j]**2
                
                idx = i*ny*(nz-1) + k*ny + j;
                diag = 0.;

                if j>0 and LHS2:
                    dy2i = 0.5*(pn[i,j]+pn[i,j-1])*pn[i,j]
                    A[idx,idx-js] = dy2i*N2[i,j-1,k];
                    diag = diag - dy2i*N2[i,j,k];
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

                if i>0 and LHS2:
                    dx2i = 0.5*(pm[i,j]+pm[i-1,j])*pm[i,j]
                    A[idx,idx-i_s] = dx2i*N2[i-1,j,k];
                    diag = diag - dx2i*N2[i,j,k];

                if i<nx-1 and LHS2:
                    dx2i = 0.5*(pm[i,j]+pm[i+1,j])*pm[i,j]
                    A[idx,idx+i_s] = dx2i*N2[i+1,j,k];
                    diag = diag - dx2i*N2[i,j,k];

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

                if j<ny-1 and LHS2:
                    dy2i = 0.5*(pn[i,j]+pn[i,j+1])*pn[i,j]
                    A[idx,idx+js] = dy2i*N2[i,j+1,k];
                    diag = diag - dy2i*N2[i,j,k];
                #else:
                    #A[idx,idx-bjs] = dy2i;
                    #diag = diag - dy2i;

                A[idx,idx] = diag
                
                if terms['VA']==1: 
                    #print 'VA_term is added to the matrix'
                    A[idx,idx] = A[idx,idx]  - f[i,j] * vrtzz[i,j,k];
                    
                    
                if terms['T2']==1:
                    print('work in proress')
                    
                    

                

    return A
       
'''   
#######################################################
# Compute solution of TTW equation
# 
#######################################################

from pyamg import *
from pyamg.gallery import *
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
#################################################
# Compute solution of TTW equation on sigma levels
# 
#######################################################


from scipy.sparse import *
from scipy.sparse.linalg import spsolve
import scipy.integrate as integrate
import time as tm
 
def solve_ttw(bx,by,AKv,sustr,svstr,f,pm,pn,z_w,timing=False,debug=0):
    
    '''
    AKv and bx,by are on vertical w levels
    
    uz,vz (solutions of TTW) also

    '''
    
    print('welcome to solve_ttw_z')
    
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

    ndim = 2*(nz+1)*(ny-2)*(nx-2);
    A = lil_matrix((ndim,ndim))
    R = np.zeros(ndim);

    for i in range(1,nx-1): 
        for j in range(1,ny-1): 
        
            idx0 = (i-1)*(ny-2)*2*(nz+1) + (j-1)*2*(nz+1)
            
            A[idx0,idx0+1] =  f[i,j];
            A[idx0+1,idx0] = -f[i,j];

            for k in range(1,nz):
                idx = idx0 + 2*k;
                dz2 = 0.5*(dz[i,j,k]+dz[i,j,k-1])
                
                A[idx,idx+ks] = AKv[i,j,k+1]/dz[i,j,k]/dz2;
                A[idx,idx] =-AKv[i,j,k]/dz[i,j,k]/dz2 - AKv[i,j,k]/dz[i,j,k-1]/dz2;
                A[idx,idx-ks] = AKv[i,j,k-1]/dz[i,j,k-1]/dz2;
                A[idx,idx+1] = f[i,j];
                
                A[idx+1,idx+1+ks] = AKv[i,j,k+1]/dz[i,j,k]/dz2;
                A[idx+1,idx+1] =-AKv[i,j,k]/dz[i,j,k]/dz2 - AKv[i,j,k]/dz[i,j,k-1]/dz2;
                A[idx+1,idx+1-ks] = AKv[i,j,k-1]/dz[i,j,k-1]/dz2;
                A[idx+1,idx] =-f[i,j];
    
            idx = idx0 +2*nz;
            A[idx,idx+1] = AKv[i,j,-1];
            A[idx+1,idx] = AKv[i,j,-1];
     
            #######################################################

            for k in range(nz):   
                idx = idx0 +2*k
                R[idx] = bx[i,j,k]
                R[idx+1] = by[i,j,k]


            idx =idx0 + 2*nz    
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


    #ml =ruge_stuben_solver(A)                # construct the multigrid hierarchy
    #print ml                                 # print hierarchy information
    #X = ml.solve(R, tol=1e-8)               # solve Ax=b to a tolerance of 1e-8
    #print "residual norm is", norm(R - A*X)  # c
    #print 'Using ruge_stuben_solver.........'



    if timing: print('computation OK.........', tm.time()-tstart)
    if timing: tstart = tm.time()         
            
            #######################################################
            
    for i in range(1,nx-1): 
        for j in range(1,ny-1):         
            # reorder results in (i,j,k)
            for k in range(nz+1):
                idx = (i-1)*(ny-2)*2*(nz+1) + (j-1)*2*(nz+1) +2*k
                uz[i,j,k] = X[idx];
                vz[i,j,k] = X[idx+1];
        
    if timing: print('allocation OK.........', tm.time()-tstart)               
    if timing: tstart = tm.time()  
    
 

    #######################################################
    # Integrate vertically to get u,v,ug,vg
    
    ut,vt,ug,vg = copy(new), copy(new), copy(new), copy(new)
    ut[:,:,1:] = integrate.cumtrapz(uz,z_w, axis=2 )
    vt[:,:,1:] = integrate.cumtrapz(vz,z_w, axis=2 )
    ug[:,:,1:] = (integrate.cumtrapz(-by,z_w, axis=2 ).T/f.T).T
    vg[:,:,1:] = (integrate.cumtrapz(bx,z_w, axis=2 ).T/f.T).T

    #######################################################

    
    if debug==1:
        return ut,vt,ug,vg,uz,vz  
    else:
        return ut,vt,ug,vg  
    
   
'''
from Modules import *
from Modules_gula import *

from pyamg import *
from pyamg.gallery import *
from scipy import *
from scipy.linalg import *  

pm = np.ones((10,10))
pn = np.ones((10,10))


A = poisson(pm.shape, format='csr') 

B = tools_g.poisson_matrix(pm,pn)
B = B.tocsr()

'''



'''
from pyamg import *
from pyamg.gallery import *
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
    print 'number of points is ',nx,ny,nz
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
    
        if i%20==0: print 'solving equation:', round(100.* i/(nx-1)) , ' %'
        
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
                print R
                print AKv[i,j,:]         
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
            
    #u = np.cumsum(uz,2)*dz
    #v = np.cumsum(vz,2)*dz  
    
    u,v = copy(new), copy(new)
    print uz.shape,len(depths)
    u[:,:,1:] = integrate.cumtrapz(uz,dx=dz, axis=2 )
    v[:,:,1:] = integrate.cumtrapz(vz,dx=dz, axis=2 )
    
    #ut[:,:,:-1] = integrate.cumtrapz(uz[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1]
    #vt[:,:,:-1] = integrate.cumtrapz(vz[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1]
    #ug[:,:,:-1] = (integrate.cumtrapz(-by[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1].T/f.T).T
    #vg[:,:,:-1] = (integrate.cumtrapz(bx[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1].T/f.T).T
    
    return u,v
'''
    
#######################################################
# Compute solution of TTW equation on sigma levels
# 
#######################################################


from scipy.sparse import *
from scipy.sparse.linalg import spsolve
import scipy.integrate as integrate
import time as tm
 
def solve_ttw_sig(bx,by,AKv,sustr,svstr,f,pm,pn,z_w,timing=False,debug=0,ekman=0):
    
    '''
    AKv and bx,by are on vertical w levels
    
    uz,vz (solutions of TTW) also

    '''
    
    print('welcome to solve_ttw_sig')
    
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
    uze,vze = copy(new),copy(new)    
    
    #######################################################
    # Create Matrix
    #######################################################

    ndim = 2*(nz+1);
    #A = lil_matrix((ndim,ndim))
    #R = np.zeros(ndim);

    for i in range(nx): 
    
        if i%20==0: print('solving TTW equation:', round(100.* i/(nx-1)) , ' %')
        
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

            if ekman==1:
                # Compute separately the Ekman part
                idx = 2*nz    
                R[idx] = 0.
                R[idx+1] = 0.
            else:
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
            
            
            #ml =ruge_stuben_solver(A)                # construct the multigrid hierarchy
            #print ml                                 # print hierarchy information
            #X = ml.solve(R, tol=1e-8)               # solve Ax=b to a tolerance of 1e-8
            #print "residual norm is", norm(R - A*X)  # c
            #print 'Using ruge_stuben_solver.........'



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

            if ekman==1: 
                
                R = np.zeros(ndim);
                        
                idx = 2*nz    
                R[idx] = svstr[i,j]
                R[idx+1] = sustr[i,j]
                
                X = spsolve(A,R)
                
                #######################################################
            
                # reorder results in (i,j,k)
                for k in range(nz+1):
                    idx = 2*k
                    uze[i,j,k] = X[idx];
                    vze[i,j,k] = X[idx+1];      
                    
                #######################################################
            

    #######################################################
    # Integrate vertically to get u,v,ug,vg
    
    ut,vt,ug,vg = copy(new), copy(new), copy(new), copy(new)
    ut[:,:,1:] = integrate.cumtrapz(uz,z_w, axis=2 )
    vt[:,:,1:] = integrate.cumtrapz(vz,z_w, axis=2 )
    ug[:,:,1:] = (integrate.cumtrapz(-by,z_w, axis=2 ).T/f.T).T
    vg[:,:,1:] = (integrate.cumtrapz(bx,z_w, axis=2 ).T/f.T).T
    
    
    # Put on rho-level
    ut = 0.5*(ut[:,:,1:] + ut[:,:,:-1])
    vt = 0.5*(vt[:,:,1:] + vt[:,:,:-1])
    ug = 0.5*(ug[:,:,1:] + ug[:,:,:-1])
    vg = 0.5*(vg[:,:,1:] + vg[:,:,:-1])


    #######################################################

    uek,vek= copy(new[:,:,:]), copy(new[:,:,:])
    uek[:,:,1:] = integrate.cumtrapz(uze,z_w, axis=2 )
    vek[:,:,1:] = integrate.cumtrapz(vze,z_w, axis=2 )   
    
    uek = 0.5*(uek[:,:,1:] + uek[:,:,:-1])
    vek = 0.5*(vek[:,:,1:] + vek[:,:,:-1])
    
    #######################################################
    
    #ut[:,:,:-1] = integrate.cumtrapz(uz[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1]
    #vt[:,:,:-1] = integrate.cumtrapz(vz[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1]
    #ug[:,:,:-1] = (integrate.cumtrapz(-by[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1].T/f.T).T
    #vg[:,:,:-1] = (integrate.cumtrapz(bx[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1].T/f.T).T
    
    
    #######################################################

    
    if debug==1 and ekman==1:
        return ut,vt,ug,vg,uz,vz,uek,vek
    elif debug==1:
        return ut,vt,ug,vg,uz,vz  
    else:
        return ut,vt,ug,vg  
    
    

#######################################################
# Compute solution of TTW equation on sigma levels
# 
#######################################################


from scipy.sparse import *
from scipy.sparse.linalg import spsolve
import scipy.integrate as integrate
import time as tm
 
def solve_ttw_sig_fast(bx,by,AKv,sustr,svstr,f,pm,pn,z_w,timing=False,debug=0):
    
    '''
    AKv and bx,by are on vertical w levels
    
    uz,vz (solutions of TTW) also

    '''
    
    print('welcome to solve_ttw_sig')
    
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

    ndim = 2*(nz+1)*(ny-2)*(nx-2);
    A = lil_matrix((ndim,ndim))
    R = np.zeros(ndim);

    for i in range(1,nx-1): 
        for j in range(1,ny-1): 
        
            idx0 = (i-1)*(ny-2)*2*(nz+1) + (j-1)*2*(nz+1)
            
            A[idx0,idx0+1] =  f[i,j];
            A[idx0+1,idx0] = -f[i,j];

            for k in range(1,nz):
                idx = idx0 + 2*k;
                dz2 = 0.5*(dz[i,j,k]+dz[i,j,k-1])
                
                A[idx,idx+ks] = AKv[i,j,k+1]/dz[i,j,k]/dz2;
                A[idx,idx] =-AKv[i,j,k]/dz[i,j,k]/dz2 - AKv[i,j,k]/dz[i,j,k-1]/dz2;
                A[idx,idx-ks] = AKv[i,j,k-1]/dz[i,j,k-1]/dz2;
                A[idx,idx+1] = f[i,j];
                
                A[idx+1,idx+1+ks] = AKv[i,j,k+1]/dz[i,j,k]/dz2;
                A[idx+1,idx+1] =-AKv[i,j,k]/dz[i,j,k]/dz2 - AKv[i,j,k]/dz[i,j,k-1]/dz2;
                A[idx+1,idx+1-ks] = AKv[i,j,k-1]/dz[i,j,k-1]/dz2;
                A[idx+1,idx] =-f[i,j];
    
            idx = idx0 +2*nz;
            A[idx,idx+1] = AKv[i,j,-1];
            A[idx+1,idx] = AKv[i,j,-1];
     
            #######################################################

            for k in range(nz):   
                idx = idx0 +2*k
                R[idx] = bx[i,j,k]
                R[idx+1] = by[i,j,k]


            idx =idx0 + 2*nz    
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


    #ml =ruge_stuben_solver(A)                # construct the multigrid hierarchy
    #print ml                                 # print hierarchy information
    #X = ml.solve(R, tol=1e-8)               # solve Ax=b to a tolerance of 1e-8
    #print "residual norm is", norm(R - A*X)  # c
    #print 'Using ruge_stuben_solver.........'



    if timing: print('computation OK.........', tm.time()-tstart)
    if timing: tstart = tm.time()         
            
            #######################################################
            
    for i in range(1,nx-1): 
        for j in range(1,ny-1):         
            # reorder results in (i,j,k)
            for k in range(nz+1):
                idx = (i-1)*(ny-2)*2*(nz+1) + (j-1)*2*(nz+1) +2*k
                uz[i,j,k] = X[idx];
                vz[i,j,k] = X[idx+1];
        
    if timing: print('allocation OK.........', tm.time()-tstart)               
    if timing: tstart = tm.time()  
    
 

    #######################################################
    # Integrate vertically to get u,v,ug,vg
    
    ut,vt,ug,vg = copy(new), copy(new), copy(new), copy(new)
    ut[:,:,1:] = integrate.cumtrapz(uz,z_w, axis=2 )
    vt[:,:,1:] = integrate.cumtrapz(vz,z_w, axis=2 )
    ug[:,:,1:] = (integrate.cumtrapz(-by,z_w, axis=2 ).T/f.T).T
    vg[:,:,1:] = (integrate.cumtrapz(bx,z_w, axis=2 ).T/f.T).T
    
    
    # Put on rho-level
    ut = 0.5*(ut[:,:,1:] + ut[:,:,:-1])
    vt = 0.5*(vt[:,:,1:] + vt[:,:,:-1])
    ug = 0.5*(ug[:,:,1:] + ug[:,:,:-1])
    vg = 0.5*(vg[:,:,1:] + vg[:,:,:-1])

    #######################################################
    
    #ut[:,:,:-1] = integrate.cumtrapz(uz[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1]
    #vt[:,:,:-1] = integrate.cumtrapz(vz[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1]
    #ug[:,:,:-1] = (integrate.cumtrapz(-by[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1].T/f.T).T
    #vg[:,:,:-1] = (integrate.cumtrapz(bx[:,:,::-1],z_w[:,:,::-1], axis=2 )[:,:,::-1].T/f.T).T
    
    
    #######################################################

    
    if debug==1:
        return ut,vt,ug,vg,uz,vz  
    else:
        return ut,vt,ug,vg  
    
   
'''
from Modules import *
from Modules_gula import *

from pyamg import *
from pyamg.gallery import *
from scipy import *
from scipy.linalg import *  

pm = np.ones((10,10))
pn = np.ones((10,10))


A = poisson(pm.shape, format='csr') 

B = tools_g.poisson_matrix(pm,pn)
B = B.tocsr()

'''







######################################################  
def poisson_matrix(pm,pn):    
######################################################  

    # elliptic equation matrix:  d_xx + d_yy
    [nx,ny] = pm.shape
    
    ############################

    ndim = ny*nx;
    i_s  = ny;
    js  = 1;
    A=lil_matrix((ndim,ndim))
    #A=np.zeros((ndim,ndim))
    ############################
    
    for i in range(nx): 
        for j in range(ny): 
                idx = i*ny + j;
                diag = 0.;
                if j>0:
                    dy2i = -0.5*(pn[i,j]+pn[i,j-1])*pn[i,j]
                    A[idx,idx-js] = dy2i;
                    diag -= dy2i;
                if i>0:
                    dx2i = -0.5*(pm[i,j]+pm[i-1,j])*pm[i,j]
                    A[idx,idx-i_s] = dx2i;
                    diag -= dx2i;
                if i<nx-1:
                    dx2i = -0.5*(pm[i,j]+pm[i+1,j])*pm[i,j]
                    A[idx,idx+i_s] = dx2i;
                    diag -= dx2i;
                if j<ny-1:
                    dy2i = -0.5*(pn[i,j]+pn[i,j+1])*pn[i,j]
                    A[idx,idx+js] = dy2i;
                    diag -= dy2i;
                A[idx,idx] = diag


    # normalize A; otherwise it is singular
    #i = nx/2;
    #j = ny/2;
    #idx0 = (i-1)*ny + j-1;
    #A[idx0,:] = 0;
    #A[idx0,idx0] = np.sqrt(dx2i*dy2i);
    
    return A




from scipy.sparse import *

######################################################  
def poisson_matrix_fast(pm,pn):    
######################################################  

    print('using poisson_matrix_fast')
    # elliptic equation matrix:  d_xx + d_yy
    [nx,ny] = pm.shape
    
    ############################

    Dxx = spdiags([-pm[:,ny/2]**2, 2*pm[:,ny/2]**2, -pm[:,ny/2]**2], [-1, 0, 1], nx, nx)
    Dyy = spdiags([-pn[nx/2,:]**2, 2*pn[nx/2,:]**2, -pn[nx/2,:]**2], [-1, 0, 1], ny, ny)
    A = kronsum(Dxx, Dyy)

    return A











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





