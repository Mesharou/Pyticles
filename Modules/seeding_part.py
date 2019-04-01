# seeding_part.py
# module 
# J.Collin 14-03-2019
#=================================================================

import numpy as np
from scipy.interpolate import interp1d
import pyticles_sig_sa as part
from copy import *

def ini_depth(maskrho,simul,depths0,x,y,z,z_wpsi):
    '''
Cubic interpolation of sigma levels at depths0 (seeding depths in meters)
z : sigma level for particles released at depths0

return z list

parameters: simul
            maskrho : land mask
            depths0 : vector with seeding depths
            x,y,z : 3 dimensionnal grid point for seeded particles
            z_w : depths in meters of sigma levels 
    '''
    z_part = np.arange(len(simul.coord[4])+1, dtype='float')

    for k in range(len(depths0)):
        for i in range(x.shape[2]):
            for j in range(x.shape[1]):
                ix,iy = np.int(np.floor(x[k,j,i])),np.int(np.floor(y[k,j,i]))
                cfx = x[k,j,i] - ix
                cfy = y[k,j,i] - iy
                
                zxp = z_wpsi[ix+1, iy]
                zyp = z_wpsi[ix, iy+1]
                zxyp = z_wpsi[ix+1, iy+1]
                
                z_part[:] = (1-cfx)*(1-cfy)*z_wpsi[ix, iy]+ (1-cfx)*cfy*zyp + cfx*(1-cfy)*zxp + cfx*cfy*zxyp

                if debug_pdeth:
                    print('---------- SEEDING----------')
                    print(f'cfx = {cfx}')
                    print(f'cfy = {cfy}')
                    print('----------------------------')
                    print(f'zxp = {zxp};   zyp = {zyp};  zxyp = {zxyp}')
                    print(f' z_part = {z_part}')

                # JC Need to check also surrounding points
                if maskrho[ix,iy]==1:
                #    f = interp1d(z_psi[ix,iy], list(range(z_wpsi.shape[2])), kind='cubic')
                 
                 #   fxp = interp1d(z_wpsi[ix+1,iy], list(range(z_wpsi.shape[2])), kind='cubic')
                 #   fyp = interp1d(z_wpsi[ix,iy+1], list(range(z_wpsi.shape[2])), kind='cubic')
                 #   fxyp = interp1d(z_wpsi[ix+1,iy+1], list(range(z_wpsi.shape[2])), kind='cubic')
                 #   #f = interp1d(z_wpsi[5,5], list(range(z_wpsi.shape[2])), kind='cubic')
                 #   zpsi = f(depths0[k])
                 #   zxp = fxp(depths0[k])
                 #   zyp = fyp(depths0[k])
                 #   zxyp = fxyp(depths0[k])
                    
                    print('-----------------')
                 #   print(f'zpsi = {zpsi};    zxp = {zxp};   zyp = {zyp};  zxyp = {zxyp}')
                 #   z[k,j,i] = (1-cfx)*(1-cfy)*zpsi + (1-cfx)*cfy*zxp + cfx*(1-cfy)*zyp + cfx*cfy*zxyp
                    
                    f = interp1d(z_part, list(range(z_wpsi.shape[2])), kind='cubic')
                    z[k,j,i] = f(depths0[k])

                    print(f'z[{k}{j}{i}] = {z[k,j,i]}' )
                else:
                    z[k,j,i] = 0.

    return z

##############################################################################
def remove_mask(simul,topolim,x,y,z,px0,py0,pz0,nq):
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

  
    ipcount = 0
    for ip in range(len(x.reshape(-1))):
        if (ptopo[ip]>topolim) and (pmask[ip]>=1.) and (ipcount<nq):
            px0.append(x.reshape(-1)[ip])
            py0.append(y.reshape(-1)[ip])
            pz0.append(z.reshape(-1)[ip])
            ipcount +=1
    return ipcount
