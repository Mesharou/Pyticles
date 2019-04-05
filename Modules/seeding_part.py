# seeding_part.py
# module 
# J.Collin 14-03-2019
#=============================================================================
# Need to check whether it crashes or not when grid_box goes out of the domain
# 


import numpy as np
from scipy.interpolate import interp1d
import pyticles_sig_sa as part
from copy import *

def seed_box(ic=0, jc=0, lev0=0, lev1=0, iwd=0, jwd=0, nx=1, ny=1, nnx=1,
             nny=1, nnlev=1):
    '''
    Spatial box for seeding partciles in sigma coordinates located at
    psi_w points
    If initial_depth : z will be redefined at depths0
    
    returns z, y, x

    parameters:
        ic, jc : center location 
        lev0, lev1 : first and final vertical levels
        iwd, jwd : box's horizontal half width
        nx, ny : domain's width
        nnx, nny, nnlev : seeding density; ex nnx = 2 (particles every 2 x_grid points)
    '''
    z, y, x = np.mgrid[lev0:lev1+1:nnlev,
              np.max([jc-jwd,1]):np.min([jc+jwd+np.min([1.,jwd]),ny]):nny,
              np.max([ic-iwd,1]):np.min([ic+iwd+np.min([1.,iwd]),nx]):nnx]

    return z, y, x

def ini_depth(maskrho,simul,depths0,x,y,z,z_w,ng=0):
    '''
Cubic interpolation of sigma levels at depths0 (seeding depths in meters)
z : sigma level for particles released at depths0

return z list

parameters: simul
            maskrho : land mask
            depths0 : vector with seeding depths
            x,y,z : 3 dimensionnal grid point for seeded particles
            z_w : depths in meters of sigma levels
            ng : nunber of ghost point for horizontal interpolation

z_part vector of size nz+1: is z_w interpolated at particles horizontal location
    '''
    z_part = np.arange(len(simul.coord[4])+1, dtype='float')

    for k in range(len(depths0)):
        for i in range(x.shape[2]):
            for j in range(x.shape[1]):

                ix = np.int(np.floor(x[k,j,i])) + ng
                iy = np.int(np.floor(y[k,j,i])) + ng
                cfx = x[k,j,i] - ix + 0.5 + ng
                cfy = y[k,j,i] - iy + 0.5 + ng
                
                zxp = z_w[ix+1,iy]
                zyp = z_w[ix,iy+1]
                zxyp = z_w[ix+1,iy+1]

               # zxp = z_wpsi[ix+1,iy]
               # zyp = z_wpsi[ix,iy+1]
               # zxyp = z_wpsi[ix+1,iy+1]
                
                z_part[:] = (1-cfx)*(1-cfy)*z_w[ix,iy] + (1-cfx)*cfy*zyp \
                        + cfx*(1-cfy)*zxp + cfx*cfy*zxyp
                
                if maskrho[ix,iy]==1:
                    f = interp1d(z_part, list(range(z_w.shape[2])), kind='cubic')
                    z[k,j,i] = f(depths0[k])

                else:
                    z[k,j,i] = 0.

                debug_pdepth = False
                if debug_pdepth:
                    print('---------- SEEDING----------')
                    print(f'cfx = {cfx}')
                    print(f'cfy = {cfy}')
                    print('----------------------------')
                    print(f'zxp = {zxp};   zyp = {zyp};  zxyp = {zxyp}')
                    print(f' z_part = {z_part}')

    return z

##############################################################################
def remove_mask(simul,topolim,x,y,z,px0,py0,pz0,nq,ng=0,pcond=np.array(False)):
    '''
    Remove particles found in land mask and particles below sea-floor if ADV2D
    Modify in place px0, py0, pz0 with particles coordinates
    Returns ipcount to keep count of seeded particles 
    '''
    # Recomputing some data to help readablility (less argument to remove_mask)
    # Carefull of Side Effects
    

    mask = simul.mask
    maskrho = copy(mask)
    maskrho[np.isnan(maskrho)] = 0.
    ptopo = part.map_topo(simul,x.reshape(-1),y.reshape(-1))
    pmask = part.map_var2d(simul,maskrho,x.reshape(-1),y.reshape(-1))
    
    debug_seed = False
    if debug_seed:
        print(f'mask = {mask}')
        print(f'maskrho = {maskrho}')
        print(f'ptopo = {ptopo}')
        print(f'pmask = {pmask}')
        print(f"topolim = {topolim} ")
        print(f'nq = {nq}')

   # print(f'pcond = {pcond}') 
    ipcount = 0
    if pcond.any():

         for ip in range(len(x.reshape(-1))):
            if (ptopo[ip]>topolim) and (pmask[ip]>=1.) and (ipcount<nq) and (pcond[ip]==1.):
                px0.append(x.reshape(-1)[ip])
                py0.append(y.reshape(-1)[ip])
                pz0.append(z.reshape(-1)[ip])
                ipcount +=1
    else:

        for ip in range(len(x.reshape(-1))):
            if (ptopo[ip]>topolim) and (pmask[ip]>=1.) and (ipcount<nq):
                px0.append(x.reshape(-1)[ip])
                py0.append(y.reshape(-1)[ip])
                pz0.append(z.reshape(-1)[ip])
                ipcount +=1
    
    return ipcount



##############################################################################
def ini_cond(simul,x,y,z,px0,py0,pz0,nq,cond):
    '''
    '''
    pcond = partF.interp3_d(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0)
    debug_cond = False
    if debug_seed:
        print(f'mask = {mask}')
        print(f'maskrho = {maskrho}')
        print(f'ptopo = {ptopo}')
        print(f'pmask = {pmask}')
        print(f"topolim = {topolim} ")
        print(f'nq = {nq}')


    for ip in range(len(x.reshape(-1))):
        if (ptopo[ip]>topolim) and (pmask[ip]>=1.) and (ipcount<nq):
            px0.append(x.reshape(-1)[ip])
            py0.append(y.reshape(-1)[ip])
            pz0.append(z.reshape(-1)[ip])
            ipcount +=1
    return ipcount

