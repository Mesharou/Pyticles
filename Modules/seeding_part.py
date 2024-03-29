# seeding_part.py
# module 
# J.Collin 14-03-2019
#=============================================================================
# Need to check whether it crashes or not when grid_box goes out of the domain
# 
# Create a file with parameters 
# 
# Bug found: When pcond = False (due to an error into initialisation)
# Then particles are released at every vertical grid-point

import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
import pyticles_sig_sa as part
from copy import copy
from netCDF4 import Dataset


##############################################################################
def prho(ptemp=0, psalt=0, pdepth=0):
    '''
    Computes density via Equation Of State (EOS) for seawater at particules.
    If so prescribed, non-linear EOS of Jackett and McDougall (1995)
    is used.
    
    ptemp potential temperature [deg Celsius].
    psalt salinity [PSU].
    
    K0, K1 and K2 are the pressure polynomial coefficients for secant
    bulk modulus, so that
    
    bulk = K0 - K1 * z + K2 * z**2 ;
    
    while rho1 is sea-water density [kg/m^3] at standard pressure
    of 1 Atm, so that the density anomaly at in-sity pressure is
    
                  rho = rho1 / (1 + z / bulk) - 1000
    
    Check Values: T=3 C S=35.5 PSU Z=-5000 m rho=1050.3639165364
    Adapted from Roms_tools
    '''
    A00 = 19092.56
    A01 = 209.8925
    A02 = -3.041638
    A03 = -1.852732*10**-3 
    A04 = -1.361629*10**-5
    A10 = 104.4077
    A11 = -6.500517
    A12 = 0.1553190 
    A13 = 2.326469*10**-4
    AS0 = -5.587545
    AS1 = 0.7390729 
    AS2 = -1.909078*10**-2 
    B00 = 4.721788*10**-1
    B01 = 1.028859*10**-2
    B02 = -2.512549*10**-4
    B03 = -5.939910*10**-7
    B10 = -1.571896*10**-2
    B11 = -2.598241*10**-4
    B12 = 7.267926*10**-6 
    BS1 = 2.042967*10**-3 
    E00 = 1.045941*10**-5 
    E01 = -5.782165*10**-10
    E02 = 1.296821*10**-7 
    E10 = -2.595994*10**-7 
    E11 = -1.248266*10**-9 
    E12 = -3.508914*10**-9

    QR = 999.842594
    Q01 = 6.793952*10**-2
    Q02 = -9.095290*10**-3
    Q03 = 1.001685*10**-4
    Q04 = -1.120083*10**-6
    Q05 = 6.536332*10**-9
    Q10 = 0.824493
    Q11 = -4.08990*10**-3 
    Q12 = 7.64380*10**-5 
    Q13 = -8.24670*10**-7 
    Q14 = 5.38750*10**-9
    QS0 = -5.72466*10**-3
    QS1 = 1.02270*10**-4
    QS2 = -1.65460*10**-6
    Q20 = 4.8314*10**-4

    sqrtpsalt = np.sqrt(psalt)

    K0 = A00 + ptemp * (A01 + ptemp * (A02 + ptemp * (A03 + ptemp * A04))) \
    + psalt * (A10 + ptemp * (A11 + ptemp * (A12 + ptemp * A13)) \
    + sqrtpsalt * (AS0 + ptemp * (AS1 + ptemp * AS2)))
    
    K1 = B00 + ptemp * (B01 + ptemp * (B02 + ptemp * B03)) \
    + psalt * (B10 + ptemp * (B11 + ptemp * B12) + sqrtpsalt * BS1);
    
    K2 = E00 + ptemp * (E01 + ptemp * E02) \
    + psalt * (E10 + ptemp * (E11 + ptemp * E12))
    
    rho1 = QR  \
    + ptemp * (Q01 + ptemp * (Q02 + ptemp * (Q03 + ptemp * (Q04 + ptemp * Q05)))) \
    + psalt * (Q10 + ptemp * (Q11 + ptemp * (Q12 + ptemp * (Q13 + ptemp * Q14))) \
    + sqrtpsalt * (QS0 + ptemp * (QS1 + ptemp * QS2)) + psalt * Q20)
    
    prho = rho1/(1 + 0.1 * pdepth / (K0 - pdepth *(K1 - pdepth * K2)))
    
    return prho


##############################################################################
def seed_box(ic=10, jc=10, lev0=0, lev1=0, iwd=2, jwd=2, nx=1, ny=1, nnx=1,
             nny=1, nnlev=1):
    '''

    Spatial box for seeding particles in sigma coordinates located at
    psi_w points
    If initial_depth : z will be redefined at depths0
    
    returns z, y, x

    parameters:
    ic, jc : center location 
    lev0, lev1 : first and final vertical levels
    iwd, jwd : box's horizontal half width
    nx, ny : domain's width
    nnx, nny, nnlev : seeding density;
    Ex nnx = 2 (particles every 2 x_grid points)

    '''
    z, y, x = np.mgrid[lev0:lev1+1:nnlev,
              np.max([jc-jwd,0]):np.min([jc+jwd+nny,ny]):nny,
              np.max([ic-iwd,0]):np.min([ic+iwd+nnx,nx]):nnx]

    return z, y, x

##############################################################################
def get_dx(ic=0, jc=0., simul=None):
    '''
    returns zonal dx in meters at [jc, ic] location

    parameters:
    float jc sigma meridionnal coordinate to compute dx
    R_files.load simul 

    '''
    dx = 1/simul.pm[int(jc), int(ic)]
    
    return dx


##############################################################################
def seed_meter(ic=10, jc=10, lev0=0, lev1=1, nnlev=1, nx_box=10, ny_box=10,
                 dx_box=2000, simul=None, ng=1, debug=False):
    '''

    Spatial box for seeding particles in sigma coordinates located at
    psi_w points
    If initial_depth : z will be redefined at depths0

    returns z, y, x
    
    Size of the patch and distance between particles in meters are conserved
    even when box's center moves during simulation

    parameters:
    ic, jc : center location (floats or int)
    lev0, lev1 : first and final vertical levels (ints)
    nnlev : seeding vertical density (int)
    nx_box, ny_box : number of intervals  in both x, y directions
    dx_box (meters):  particle's spacing in meters

    '''
    nx, ny = simul.pm.shape
    dx = get_dx(ic=ic, jc=jc, simul=simul)
    # index for mesh
    j0 = np.max([jc - (ny_box)/2*dx_box/dx, 1])
    j1 = np.min([jc + (ny_box+1)/2*dx_box/dx, ny + ng])
    
    i0 = np.max([ic - nx_box/2*dx_box/dx, 1])
    i1 = np.min([ic + (nx_box+1)/2*dx_box/dx, nx + ng])

    z, y, x = np.mgrid[lev0:lev1+1:nnlev, j0:j1:dx_box/dx, i0:i1:dx_box/dx]
    
    if debug:
        print(f"nx, ny {nx, ny}")
        print(f"dx {dx}")
        print(f"i0, i1, j0, j1 {i0, i1, j0 ,j1}")
        print(f"dx_box {dx_box}")

    return z, y, x
    
    
    
##############################################################################

def particles_on_a_sphere(d, lonmin=-180,lonmax=180,latmin=-90,latmax = -90):
    '''
    gives lon,lat position for particles on a sphere
    with a distance of d [ in meters ] between particles
    '''
    
    #earth radius [in m]
    r = 6371e3

    #####
    # initialization
    i = 0
    plon,plat = [],[]
    
    #####
    Mtheta = int(np.pi * r / d)
    dtheta = np.pi * r / Mtheta
    dphi = d**2 / dtheta

    for m in range(Mtheta):
        theta = np.pi * (m + 0.5) / Mtheta
        Mphi = int( 2 * np.pi * r * np.sin(theta) / dtheta)

        for n in range(Mphi):
            phi = 2 * np.pi * n / Mphi
            ####
            plon_tmp = phi * 180 / np.pi - 180.
            plat_tmp = theta * 180 / np.pi - 90.
            
            if lonmax>plon_tmp>lonmin and latmax>plat_tmp>latmin:
                plon.append(plon_tmp)
                plat.append(plat_tmp)

            ####
            i +=1

    return plon,plat

    
##############################################################################
def seed_meter_sphere(lev0=0, lev1=1, nnlev=1,
                     dx_box=2000,
                     box_coords = [-180,180,-90,90],
                     simul=None, sig=None, ng=1, debug=False):
    '''

    Seeding particles with (approximate) constant distance on a spherical grid:

    returns z, y, x
    
    parameters:
    ic, jc : center location (floats or int)
    lev0, lev1 : first and final vertical levels (ints)
    nnlev : seeding vertical density (int)
    dx_box (meters):  particle's spacing in meters
    sig: sigma layer to initialize particles
    
    '''
    
    # define the domain
    [lonmin,lonmax,latmin,latmax] = box_coords
    
    lonmin = np.min([lonmin,simul.x.min()])
    lonmax = np.min([lonmax,simul.x.max()])
    latmin = np.min([latmin,simul.y.min()])
    latmax = np.min([latmax,simul.y.max()])
    
    # define particles positions on lat/lon
    plon,plat = particles_on_a_sphere(d=dx_box,\
                                  lonmin=simul.x.min(),lonmax=simul.x.max(),\
                                  latmin=simul.y.min(),latmax=simul.y.max())

    # project on grid indices
    y,x = np.meshgrid(range(simul.x.shape[1]),range(simul.x.shape[0]))
    px = griddata((simul.x.ravel(),simul.y.ravel()), x.ravel(), (plon,plat), method='linear')
    py = griddata((simul.x.ravel(),simul.y.ravel()), y.ravel(), (plon,plat), method='linear')

    #remove particles outside domain
    py = py[np.isfinite(px)]; px = px[np.isfinite(px)]
    px = px[np.isfinite(py)]; py = py[np.isfinite(py)]
    py = py[np.logical_and(simul.x.shape[0]-1>px , px>0)]
    px = px[np.logical_and(simul.x.shape[0]-1>px , px>0)]
    px = px[np.logical_and(simul.x.shape[1]-1>py , py>0)]
    py = py[np.logical_and(simul.x.shape[1]-1>py , py>0)]
    
    z, x = np.meshgrid(range(lev0,lev1+1,nnlev), px)
    z, y = np.meshgrid(range(lev0,lev1+1,nnlev), py)
    
    return z, y, x
    

##############################################################################

def ini_surf(simul, rho0, x, y, z, rho, ng=0):
    ''' Interpolate sigma-levels onto rho0-isosurface 

    simul : 
    rho0 : vector of value for iso surface
    rho : variable at psi_w points (same as particles)
    x, y : horizontal location where particles will be seeded
    z : sigma levels over the whole water-column
    ng: number of ghost points 

    '''
    mask = simul.mask
    for k in range(len(rho0)):
        for i in range(x.shape[2]):
            for j in range(x.shape[1]):
                if mask[i,j]==1:
                    f = interp1d(rho[:, j, i], list(range(rho.shape[0])),
                                 kind='cubic')
                    z[k, j, i] = f(rho0[k])
                else:
                    z[k, j, i] = 0.
    return z

##############################################################################
def ini_depth(maskrho, simul, depths0, x, y, z, z_w, ng=0):
    ''' 
    depths0 < 0 : depth in meter
    depths0 = 0 : ocean surface
    depths0 > 0 : sigma layer
    
    else:
    Cubic interpolation of sigma levels at depths0 (seeding depths in meters)
    z : sigma level for particles released at depths0

    return z list

    parameters: simul
                maskrho : land mask
                depths0 : vector with seeding depths
                x, y, z : 3 dimensionnal coordinates for seeded particles
                z_w : depths in meters of sigma levels
                ng : nunber of ghost point for horizontal interpolation

    z_part vector of size nz+1: is z_w interpolated at particles horizontal location
    
    careful: depths0 = [0] means sea-s

    '''
    
    z_part = np.arange(len(simul.coord[4])+1, dtype='float')

    for k in range(len(depths0)):
        # surface layer
        if depths0[k] == 0:
            z[k, :, :] = simul.coord[4].max()
            
        # sigma layer
        elif depths0[k] > 0:
            z[k, :, :] = depths0[k]
            
        # vertical depth  
        else:
            for i in range(x.shape[2]):
                for j in range(x.shape[1]):
                    ix = int(np.floor(x[k, j, i])) + ng
                    iy = int(np.floor(y[k, j, i])) + ng
                    if maskrho[ix,iy]==1:
                        cfx = x[k, j, i] - ix + 0.5 + ng
                        cfy = y[k, j, i] - iy + 0.5 + ng
                
                        zxp = z_w[ix+1, iy]
                        zyp = z_w[ix, iy+1]
                        zxyp = z_w[ix+1, iy+1]

                        z_part[:] = (1-cfx) * (1-cfy) * z_w[ix, iy] \
                        + (1-cfx)*cfy*zyp + cfx*(1-cfy)*zxp + cfx*cfy*zxyp
                
                        f = interp1d(z_part, list(range(z_w.shape[2])), kind='cubic')
                        try:
                            z[k,j,i] = f(depths0[k])
                        except ValueError:
                            z[k,j,i] = -1.
                      
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
def ini_trap_depth(maskrho, simul, depths0, x, y, z, z_w, ng=0):
    '''
    retrieve pz from px, py (sigman coordinate), and depths0 in meter being vectors
    '''
    z_part = np.arange(len(simul.coord[4])+1, dtype='float')
    for ip in range(len(x)):
        ix = int(np.floor(x[ip])) + ng
        iy = int(np.floor(y[ip])) + ng
        if maskrho[ix,iy]==1:
            cfx = x[ip] - ix + 0.5 + ng
            cfy = y[ip] - iy + 0.5 + ng

            zxp = z_w[ix+1, iy]
            zyp = z_w[ix, iy+1]
            zxyp = z_w[ix+1, iy+1]

            z_part[:] = (1-cfx) * (1-cfy) * z_w[ix, iy] \
               + (1-cfx)*cfy*zyp + cfx*(1-cfy)*zxp + cfx*cfy*zxyp

            f = interp1d(z_part, list(range(z_w.shape[2])), kind='cubic')
            try:
                z[ip] = f(depths0)
            except ValueError:
                print('')
                print(np.min(z_part))
                print(np.max(z_part))
                z[ip] = -1.
            '''
            Need to handle exeption where 
                 
            raise ValueError("A value in x_new is below the interpolation "
            ValueError: A value in x_new is below the interpolation range.
                
            just need to put z=0
            '''
        else:
            z[ip] = 0.

    return z



##############################################################################
def ini_trap(trap_file, simul, maskrho, itime_fwd=0, x_periodic=False,
             y_periodic=False, ng=0):
    '''
    Release particles from netcdf pyticles file 'trap_file'
    returns nq_1save, ipmx, px0, py0, pz0

    parameters: trap_file: netcdf file with particles position
                simul: simulname
                maskrho: mask land at rho points
                itime_fwd : time index to start in trap_file
                ng: number of ghost points for linear interpolation
                    see input_file
    
    retrieve particles position from trap_file
    filter NANs
    
    keep number of particles at itime (including NANs) in order to allocate
    sufficient memory for px, py, pz

    retrieve sigma level pz0 corresponding to advection depth form trap_file
    '''
    ######
    nc = Dataset(trap_file, 'r')
    ipmx = len(nc.variables['px'][itime_fwd, :])
    depths0 = [getattr(nc, 'depth')]
    indx = ~np.isnan(nc.variables['px'][itime_fwd, :])
    indy = ~np.isnan(nc.variables['py'][itime_fwd, :])
    px0 = nc.variables['px'][itime_fwd, indx]
    py0 = nc.variables['py'][itime_fwd, indx]
    pz0 = 0 * px0
    nc.close()
    ######
    #maskrho = simul.maskrho
    z_w = part.get_depths_w(simul, x_periodic=x_periodic,
                            y_periodic=y_periodic, ng=ng)
    pz0 = ini_trap_depth(maskrho, simul, depths0, px0, py0, pz0, z_w, ng=ng)
    nq_1save = len(pz0)
    return nq_1save, ipmx, px0, py0, pz0
##############################################################################

def remove_mask(simul,maskrho, topolim, x, y, z, px0, py0, pz0, nq, ng=0,
                pcond=np.array(False)):
    '''
    Remove particles found in land mask and particles below sea-floor if ADV2D
    Also chech for particles below floor in advd3d and(z.reshape(-1)[ip]>=0.)
    Modify in place px0, py0, pz0 with particles coordinates
    Returns ipcount to keep count of seeded particles
    
    Found an issue due to machine precision
    After mask interpolation on particules you might get 0.999...999
    And therefore remove particles from patch area
    Therefore we introduce eps to fiw this numerical issue
    
    '''
    
    eps = 1e-8

    ptopo = part.map_topo(simul,x.reshape(-1),y.reshape(-1))
    pmask = part.map_var2d(simul,maskrho,x.reshape(-1),y.reshape(-1)) + eps
    
    debug_seed = False
    if debug_seed:
        print(f'mask = {mask}')
        print(f'maskrho = {maskrho}')
        print(f'ptopo = {ptopo}')
        print(f'pmask = {pmask}')
        print(f"topolim = {topolim} ")
        print(f'nq = {nq}')

    ipcount = 0
    if pcond.any():

         for ip in range(len(x.reshape(-1))):
            if (ptopo[ip]>topolim) and (pmask[ip]>=1.) and (ipcount<nq)\
            and (pcond[ip]==1.) and(z.reshape(-1)[ip]>=0.):
                px0.append(x.reshape(-1)[ip])
                py0.append(y.reshape(-1)[ip])
                pz0.append(z.reshape(-1)[ip])
                ipcount +=1
    else:

        for ip in range(len(x.reshape(-1))):
            if (ptopo[ip]>topolim) and (pmask[ip]>=1.) and (ipcount<nq) \
            and(z.reshape(-1)[ip]>=0.):
                px0.append(x.reshape(-1)[ip])
                py0.append(y.reshape(-1)[ip])
                pz0.append(z.reshape(-1)[ip])
                ipcount +=1
            else:
                if debug_seed:
                    print(f'ipcount = {ipcount}')
                    print(f'ptopo[ip] {ptopo[ip]}')
                    print(f'pmask[ip] {pmask[ip]}')
    return ipcount


##############################################################################

def remove_mask_surf(simul,maskrho, x, y, px0, py0, nq, ng=0,
                pcond=np.array(False)):
    '''
    Remove particles found in land mask
    '''
    
    eps = 1e-8
    pmask = part.map_var2d(simul,maskrho,x.reshape(-1),y.reshape(-1)) + eps

    ipcount = 0
    if pcond.any():

         for ip in range(len(x.reshape(-1))):
            if (pmask[ip]>=1.) and (ipcount<nq)\
            and (pcond[ip]==1.):
                px0.append(x.reshape(-1)[ip])
                py0.append(y.reshape(-1)[ip])

                ipcount +=1
    else:

        for ip in range(len(x.reshape(-1))):
            if (pmask[ip]>=1.) and (ipcount<nq):
                px0.append(x.reshape(-1)[ip])
                py0.append(y.reshape(-1)[ip])
                ipcount +=1

    return ipcount



##############################################################################
#def from_nc(ncname, itime)
#    '''
#    
#    '''
