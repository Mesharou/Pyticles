# seeding_part.py
# module 
# J.Collin 14-03-2019
#=================================================================
# parameters :
# particles location:
# [ic, jc] seeding center location 
# iwd, jwd : half width of seeding patch [in grid points] (x, y) direction
# depths0 : depths vector for seeding patch
# nx, ny, nz partciles shape 
# nnx, nny seeding horrizontal density
# ex: nnx, nny = 1, 2 => every x grid point, every 2 y grid-point:  
# lev0, lev1 !!! To see might be an issue here ?
# nnlev : vertical density
# mashrho : 
# simul 
# ng : number of ghost points 

import numpy as np
from scipy.interpolate import interp1d
import pyticles_sig_sa as part
from copy import *

def ini_depth(maskrho, simul, depths0, x, y, z, z_w):

    for k in range(len(depths0)):
        for i in range(x.shape[2]):
            for j in range(x.shape[1]):
                ix,iy = np.int(np.floor(x[k,j,i])),np.int(np.floor(y[k,j,i]))
                if maskrho[ix,iy]==1:
                    f = interp1d(z_w[ix,iy], list(range(z_w.shape[2])), kind='cubic')
                    z[k,j,i] = f(depths0[k])
                else:
                    z[k,j,i] = 0.

    return z
##############################################################################
def remove_mask(simul,topolim,x,y,z,px0,py0,pz0,nq):
    
    # Recomputing some data to help readablility (less argument to remove_mask)
    # Carefull of Side Effects
    mask = simul.mask
    maskrho = copy(mask)
    maskrho[np.isnan(maskrho)] = 0.
    ptopo = part.map_topo(simul,x.reshape(-1),y.reshape(-1))
    pmask = part.map_var2d(simul,maskrho,x.reshape(-1),y.reshape(-1))
    print(topolim)
    ipcount = 0
    for ip in range(len(x.reshape(-1))):
        if (ptopo[ip]>topolim) and (pmask[ip]>=1.) and (ipcount<nq):
            px0.append(x.reshape(-1)[ip])
            py0.append(y.reshape(-1)[ip])
            pz0.append(z.reshape(-1)[ip])
            ipcount +=1
    return ipcount
