# coding=utf-8



###################################################################################
# VARIABLES 
###################################################################################

"""

16/01/20: Modif for periodic files of JC [add parameter iper,jper], because xi-u has the same size than xi-rho
          So if periodic we shift indices 1 to the right


"""
from __future__ import print_function


###################################################################################
#Load modules
###################################################################################

import os

#for numeric functions
#from builtins import range
#from builtins import object
import numpy as np

#for netcdf files
#from Scientific.IO.NetCDF import *
from netCDF4 import Dataset

#For nested loop
from itertools import product

#copy data
from copy import copy

#convert seconds to date
import time as tm

###################################################################################


#ROMSTOOLS
import R_tools as tools
import R_tools_fort as toolsF

import R_tools_gula as tools_g
import R_tools_fort_gula as toolsF_g
#import R_tools_fort_cuc as toolsF_cuc
import R_tools_fort_gigatl as toolsF_gi

#import streamfunction as st

#Simulations (path, data...)
import simulations_old as oldsim

import R_smooth as sm


##############

from gc import *
from datetime import datetime, timedelta

###################################################################################


class var(object):

    ''' 
        This class is used to load/compute/interpolate any variable from the ROMS outputs

        and is used as follows: 
            myvar = R_vars.var(variable_name, simul)

        where variable_name is a string 
            (to see all valid variable names one should type: __???__ )

        where simul is a R_files object which has to include:
            simul.ncfile (path to netcdf file)
            simul.infiletime (time step (0,1,2...) in the netcdf file)
            simul.coord (subdomain coordinates [ny1,ny2,nx1,nx2,depths])
            

        NOTE that depths is used as follows:

        depth <= 0 means depths in meters
        depth = 0 means surface
        depth > 0 means sigma-level (1 = bottom [0 in netcdf file],...
                             ..., Nz = surface [Nz-1 in netcdf file])

                                    
    '''


    dico =\
    {'temp': ['Temperature', r'$T\,(^{\circ}C)$', [0,0,1]],\
    'salt': ['Salinity', 'PSU', [0,0,1]],\
    'u': ['u', 'm/s', [1,0,1]],\
    'v': ['v', 'm/s', [0,1,1]],\
    'ubar': ['ubar', 'm/s', [1,0,-1]],\
    'vbar': ['vbar', 'm/s', [0,1,-1]],\
    'zeta': ['SSH', r'$\eta\,(m)$' ,[0,0,-1]],\
    'hbls': ['Thickness of KPP surface boundary layer', 'm', [0,0,-1]],\
    'hbbl': ['Thickness of KPP bottom boundary layer', 'm', [0,0,-1]],\
    'hbls_rho': ['Surface mixed-layer (based on rho = rhos+ 0.03)', 'm', [0,0,-1]],\
    'hbls_t': ['Surface mixed-layer (based on t = ts - 0.2)', 'm', [0,0,-1]],\
    'hbls_akv': ['Surface mixed-layer (based on Akv<2e-4)', 'm', [0,0,-1]],\
    'Ekman': ['Ekman number (based on hbls_akv)', ' ', [0,0,-1]],\
    'AKt': ['Temperature vertical diffusion coef', 'm2/s', [0,0,0]],\
    'AKv': ['Momentum vertical diffusion coef', 'm2/s', [0,0,0]],\
    'omega': ['S-coordinate vertical velocity', 'm/s ?', [0,0,0]],\
    \
    'sustr': ['sustr', 'm/s', [1,0,-1]],\
    'svstr': ['svstr', 'm/s', [0,1,-1]],\
    \
    'psi': ['psi', 'Streamfunction' ,[1,1,-1]],\
    'psi_v2': ['psi', 'Streamfunction' ,[1,1,-1]],\
    'zonal_psi': ['zonal_psi', 'Streamfunction' ,[1,1,-1]],\
    'meridional_psi': ['meridional_psi', 'Streamfunction' ,[1,1,-1]],\
    'int_psi': ['psi', 'Streamfunction' ,[1,1,-1]],\
    'psi_surf': ['psi_surf', 'Streamfunction' ,[1,1,-1]],\
    'psir': ['psir', 'Streamfunction of rotational part' ,[1,1,-1]],\
    \
    'rho': ['in-situ density', 'kg.m-3', [0,0,1]],\
    'rho1': ['in-situ density1', 'kg.m-3', [0,0,1]],\
    'rho1_sol2': ['in-situ density', 'kg.m-3', [0,0,1]],\
    'rhop': ['potential density', 'kg.m-3', [0,0,1]],\
    'bvf': ['Brunt-Vaisala Frequency squared: N2', 's-2', [0,0,0]],\
    'buoy': ['buoyancy', 'm/s-2', [0,0,1]],\
    'buoy1': ['buoyancy1', 'm/s-2', [0,0,1]],\
    \
    'w': ['Vertical velocity', r'$w\,(m\,s^{-1})$', [0,0,1]],\
    'w_vmix': ['Vartical mixing contribution to vertical velocity', 'm/s', [0,0,1]],\
    'w_ttw': ['TTW vertical velocity', 'm/s', [0,0,1]],\
    \
    'absvrt': ['Absolute Vorticity', 's-1' ,[1,1,1]],\
    'vrt': ['Relative Vorticity', r'$\frac{\zeta}{f}$' ,[1,1,1]],\
    'str': ['Strain', 'S' ,[0,0,1]],\
    'ow': ['Okubo-Weiss parameter', 'ow' ,[0,0,1]],\
    'pv': ['Potential Vorticity', 'PVU' ,[1,1,0]],\
    'pv1': ['[f + (dv/dx - du/dy)]*db/dz', 'PVU' ,[1,1,2]],\
    'pv2': ['-(dv/dz)*(db/dx)', 'PVU' ,[1,1,2]],\
    'pv3': ['(du/dz)*(db/dy)', 'PVU' ,[1,1,2]],\
    'stretching': ['stretching', ' ' ,[0,0,2]],\
    'pvr': ['Potential Vorticity on rho levels', 'PVU' ,[1,1,1]],\
    'Ri': ['Richardson Number', ' ', [0,0,0]],\
    'dzu': ['Vertical shear of U', ' ', [1,0,0]],\
    'dzv': ['Vertical shear of V', ' ', [0,1,0]],\
    'dzu2': ['Vertical shear amplitude', ' ', [0,0,0]],\
    'kdzu': ['AKv*(du/dz)', ' ', [0,0,0]],\
    'kdzv': ['AKv*(dv/dz)', ' ', [0,0,0]],\
    'kdzu2': ['Interior Energy dissipation due to Vmix', ' ', [0,0,0]],\
    \
    'vortadv': ['NL advective term of vort balance', ' ' ,[1,1,-1]],\
    'vortadv_sol2': ['NL advective term of vort balance (as computed in ROMS)', ' ' ,[1,1,-1]],\
    'int_vortadv_sol2': ['NL advective term of vort balance (as computed in ROMS)', ' ' ,[1,1,-1]],\
    'vortadv_sol3': ['NL advective term of vort balance (as computed in ROMS) with curvgrid', ' ' ,[1,1,-1]],\
    'vortadv_mix': ['mix advective term of vort balance (as computed in ROMS)', ' ' ,[1,1,-1]],\
    'vortadv_centered': ['centered advective term of vort balance (as computed in ROMS)', ' ' ,[1,1,-1]],\
    'vortadv_uvgrid': ['curvgrid term', ' ' ,[1,1,-1]],\
    'rotwind': ['Rot. of Wind Stress', ' ' ,[1,1,-1]],\
    'rotbot': ['Rot of bottom Stress', ' ' ,[1,1,-1]],\
    'bpt': ['bottom presssure torque', ' ' ,[1,1,-1]],\
    'bpts': ['bottom presssure torque term 1: J(zeta,h)', ' ' ,[1,1,-1]],\
    'fwb': ['bottom vortex stretching', ' ' ,[0,0,-1]],\
    'fws': ['surface vortex stretching', ' ' ,[0,0,-1]],\
    'fdivub': ['bottom vortex stretching', ' ' ,[0,0,-1]],\
    'fwdivub': ['bottom vortex stretching', ' ' ,[0,0,-1]],\
    'vortbar_int': ['Baro. vorticity (int)', ' ', [1,1,-1]],\
    'intvrt': ['integrated vorticity', ' ', [1,1,-1]],\
    'vortplanet': ['planetary vorticity term of vort balance', ' ' ,[1,1,-1]],\
    'vortstretch': ['planetary stretching term of vort balance', ' ' ,[1,1,-1]],\
    'vortstretch2': ['planetary stretching term of vort balance', ' ' ,[1,1,-1]],\
    'int_vortplanet': ['planetary vorticity term of vort balance integrated only from 0 to depth', ' ' ,[1,1,-1]],\
    'vortstretch_sol2': ['planetary stretching term of vort balance', ' ' ,[1,1,-1]],\
    'vortplantot': ['planetary vorticity and stretching term of vort balance', ' ' ,[1,1,-1]],\
    'vortplantot_sol2': ['Cor term of vorticity balance (as computed in ROMS)', ' ' ,[1,1,-1]],\
    'u_Prsgrd': ['bottom presssure torque', ' ' ,[1,0,1]],\
    \
    'vortadv_sol2_mean': ['NL advective term of vort balance (as computed in ROMS)', ' ' ,[1,1,-1]],\
    'vortadv_sol3_mean': ['NL advective term of vort balance (as computed in ROMS)', ' ' ,[1,1,-1]],\
    'rotwind_mean': ['Rot. of Wind Stress', ' ' ,[1,1,-1]],\
    'rotbot_mean': ['Rot of bottom Stress', ' ' ,[1,1,-1]],\
    'bpt_mean': ['bottom presssure torque', ' ' ,[1,1,-1]],\
    'vortbar_mean': ['Baro. vorticity (mean)', ' ', [1,1,-1]],\
    'vortplanet_mean': ['planetary vorticity term of vort balance', ' ' ,[1,1,-1]],\
    'vorttopo_mean': ['planetary vorticity term of vort balance', ' ' ,[1,1,-1]],\
    'vortf_mean': ['planetary vorticity term of vort balance', ' ' ,[1,1,-1]],\
    'vortstretch_mean': ['planetary stretching term of vort balance', ' ' ,[1,1,-1]],\
    'vortstretch_sol2_mean': ['planetary stretching term of vort balance', ' ' ,[1,1,-1]],\
    'vortplantot_mean': ['planetary vorticity and stretching term of vort balance', ' ' ,[1,1,-1]],\
    'vortplantot_sol2_mean': ['Cor term of vorticity balance (as computed in ROMS)', ' ' ,[1,1,-1]],\
    \
    'J1_sol1': ['Surface Buoy PV flux', ' ' ,[1,1,-1]],\
    'J1_sol2': ['Surface Buoy PV flux', ' ' ,[1,1,-1]],\
    'J2_sol1': ['Surface Wind PV flux', ' ' ,[1,1,-1]],\
    'J2_sol2': ['Surface Wind PV flux', ' ' ,[1,1,-1]],\
    'Jbot_sol1': ['Bottom PV flux', ' ' ,[1,1,-1]],\
    'Jbot_sol2': ['Bottom PV flux', ' ' ,[1,1,-1]],\
    'Jbot_sol1_nohbbls': ['Bottom PV flux', ' ' ,[1,1,-1]],\
     \
    'Tadv': ['Temp. advection', 'C.s-1' ,[0,0,1]],\
    'THdiff': ['Temp. implicit horizontal diff.', 'C.s-1' ,[0,0,1]],\
    'TVmix': ['Temp. vertical diffusion', 'C.s-1' ,[0,0,1]],\
    'TForc': ['Temp. boundary forcing', 'C.s-1' ,[0,0,1]],\
    'Trate': ['Temp. rate of change', 'C.s-1' ,[0,0,1]],\
    'Sadv': ['Salt advection', 'PSU.s-1' ,[0,0,1]],\
    'SHmix': ['Salt implicit horizontal diff.', 'PSU.s-1' ,[0,0,1]],\
    'SVmix': ['Salt vertical diffusion', 'PSU.s-1' ,[0,0,1]],\
    'SForc': ['Salt boundary forcing', 'PSU.s-1' ,[0,0,1]],\
    'Srate': ['Salt rate of change', 'PSU.s-1' ,[0,0,1]],\
    'TXadv': ['Temp. advection', 'C.s-1' ,[0,0,1]],\
    'TYadv': ['Temp. advection', 'C.s-1' ,[0,0,1]],\
    'TVadv': ['Temp. advection', 'C.s-1' ,[0,0,1]],\
    \
    'Tadv_pert': ['Temp. advection', 'C.s-1' ,[0,0,1]],\
    'THdiff_pert': ['Temp. implicit horizontal diff.', 'C.s-1' ,[0,0,1]],\
    'TVmix_pert': ['Temp. vertical diffusion', 'C.s-1' ,[0,0,1]],\
    'TForc_pert': ['Temp. boundary forcing', 'C.s-1' ,[0,0,1]],\
    'Trate_pert': ['Temp. rate of change', 'C.s-1' ,[0,0,1]],\
    'Sadv_pert': ['Salt advection', 'PSU.s-1' ,[0,0,1]],\
    'SHmix_pert': ['Salt implicit horizontal diff.', 'PSU.s-1' ,[0,0,1]],\
    'SVmix_pert': ['Salt vertical diffusion', 'PSU.s-1' ,[0,0,1]],\
    'SForc_pert': ['Salt boundary forcing', 'PSU.s-1' ,[0,0,1]],\
    'Srate_pert': ['Salt rate of change', 'PSU.s-1' ,[0,0,1]],\
    \
    'Madv': ['Momentum advection', 'm.s-2' ,[0,0,1]],\
    'MHmix': ['Momentum explicit horizontal diff.', 'm.s-2' ,[0,0,1]],\
    'MHdiss': ['Momentum implicit horizontal diff.', 'm.s-2' ,[0,0,1]],\
    'MVmix': ['Momentum vertical diffusion', 'm.s-2' ,[0,0,1]],\
    'MCor': ['Momentum Coriolis', 'm.s-2' ,[0,0,1]],\
    'Mrate': ['Momentum rate of change', 'm.s-2' ,[0,0,1]],\
    'MXadv': ['Momentum advection', 'm.s-2' ,[0,0,1]],\
    'MYadv': ['Momentum advection', 'm.s-2' ,[0,0,1]],\
    'MPrsgrd': ['Momentum pressure gradient', 'm.s-2' ,[0,0,1]],\
    \
    \
    'u_prsgrd': ['online Momentum pressure gradient', 'm.s-2' ,[1,0,1]],\
    'v_prsgrd': ['online Momentum pressure gradient', 'm.s-2' ,[0,1,1]],\
    \
    'kediss': ['dissipation', ' ' ,[0,0,1]],\
    \
    'dxbuoy2': ['buoyancy gradients', ' ' ,[0,0,1]],\
    'dxu2': ['velocity gradients', ' ' ,[0,0,1]],\
    'tendency': ['Horizontal tendency for buoyancy', ' ' ,[0,0,1]],\
    'tendency_3d': ['Advective tendency for buoyancy', ' ' ,[0,0,1]],\
    'tendency_u': ['Horizontal tendency for buoyancy', ' ' ,[0,0,1]],\
    'tendency_3d_u': ['Advective tendency for buoyancy', ' ' ,[0,0,1]],\
    'tendency_full_u': ['Tendency for buoyancy', ' ' ,[0,0,1]]\
    }
    
    
    
    
    '''
    
    
#    'vrt': ['relative vorticity', 's-1', [1,1,1]],\

#    'bvf': ['Brunt-Vaisala Frequency squared: N2', 's-2', [0,0,0]],\



    'J1_sol1': ['Surface Buoy PV flux', ' ' ,[1,1,-1]],\
    'J1_sol2': ['Surface Buoy PV flux', ' ' ,[1,1,-1]],\
    'J2_sol1': ['Surface Wind PV flux', ' ' ,[1,1,-1]],\
    'J2_sol2': ['Surface Wind PV flux', ' ' ,[1,1,-1]],\
    'Jbot_sol1': ['Bottom PV flux', ' ' ,[1,1,-1]],\
    'Jbot_sol2': ['Bottom PV flux', ' ' ,[1,1,-1]],\
    'Jbot_sol1_nohbbls': ['Bottom PV flux', ' ' ,[1,1,-1]],\
    'pv_sol2': ['Potential Vorticity', 'PVU' ,[1,1,0]],\
    'pv_sol1': ['Potential Vorticity', 'PVU' ,[1,1,0]],\
    'vortbar': ['Baro. vorticity (moy)', ' ', [1,1,-1]],\
    'vortbar_int': ['Baro. vorticity (int)', ' ', [1,1,-1]],\
    'vortplanet': ['planetary vorticity term of vort balance', ' ' ,[1,1,-1]],\
    'vortstretch': ['planetary stretching term of vort balance', ' ' ,[1,1,-1]],\
    'vortplantot': ['planetary stretching term of vort balance', ' ' ,[1,1,-1]],\
    'vortadv': ['NL advective term of vort balance', ' ' ,[1,1,-1]],\
    'vortadv_sol1': ['NL advective term of vort balance', ' ' ,[1,1,-1]],\
    'vortadv_sol2': ['NL advective term of vort balance', ' ' ,[1,1,-1]],\
    'bpt': ['bottom presssure torque', ' ' ,[1,1,-1]],\
    'u_Prsgrd': ['bottom presssure torque', ' ' ,[1,1,-1]],\
    'bpt_sol2': ['bottom presssure torque', ' ' ,[1,1,-1]],\
    'rotwind': ['Rot. of Wind Stress', ' ' ,[1,1,-1]],\
    'rotbot': ['Rot of bottom Stress', ' ' ,[1,1,-1]]}
    '''



###################################################################################
#Load variables
###################################################################################

    def __init__(self,varname,simul,n2max=50000,method='new',debug=False,**kwargs):

    
        if debug: print('using R_vars gula version')
    

        if debug: print('#######################################')
        if debug: print('load var', varname)
        if debug: print('#######################################')


        self.name = copy(varname)
        self.dictionnary()     
        self.coord = copy(simul.coord)

        #self.coord are coordinates of variable (not necessarily the same than simul.coord, can be a subdomain)
        # if coordinates are given directly as an argument, use it instead of simul.coord[0:4]
        if 'coord' in  kwargs: self.coord[0:4] = kwargs['coord']
        [ny1,ny2,nx1,nx2,depths] = self.coord

        # if depths is given directly as an argument, use it instead of simul.coord[4]
        if 'depths' in  kwargs: depths = kwargs['depths']

        # Load NETCDF file
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
                
        if 'u' in  kwargs: 
            u = kwargs['u']; v = kwargs['v']
        else:
            u = None; v= None
        
        
        if 'masked' in  kwargs: 
            masked_value =  kwargs['masked']
        else: 
            masked_value =  np.nan
            
     
        ############################# 
        # test if periodic
        #############################
        xperiodic = False; iper=0
        yperiodic = False; jper=0

        try:
            if ncfile.dimensions['xi_u'].size==ncfile.dimensions['xi_rho'].size:
                xperiodic = True; iper=1
        except:
            pass

        try:
            if ncfile.dimensions['eta_v'].size==ncfile.dimensions['eta_rho'].size:
                yperiodic = True; jper=1
        except:
            pass

       
 
        '''
        Different scenarii:
            1. var is already in output
                1.1...  We want a 2D variable (nothing to do)
                1.2...  We want a 3D variable on one or more sigmal-levels (nothing to do)
                1.3...  We want a 3D variable on z-levels (z interpolation needed)
            2. var is not in file (computations needed)
                2.1...  We want it on sigma level (computation on sigma)
                2.2...  We want it on z-levels
                    2.2.1 Computation on sigma + z-interpolation
                    2.2.2 Computation on z-levels    
        '''
        
        ####################################################################
        #1. var is already in netcdf output 
        ####################################################################
        if (self.name in list(ncfile.variables.keys())) and (method=='new'): # and (self.name not in ['w']):
        #check if variable is already in output:  



            ##################################################################
            # load attributes
            try:
                self.unit = ncfile.variables[varname].units
            except:
                print('no units in file')         

            try:
                self.longname = ncfile.variables[varname].long_name
            except:
                print('no long_name in file')
            ##################################################################


        
            if debug: print('self.imin, self.jmin ,self.kmin')
            if debug: print(self.imin, self.jmin ,self.kmin)
            if debug: print(' ')
            
            if debug: print('min(depths)' , np.min(depths))

            #Check type of vertical and horizontal grid
            self.iper = 0; self.jper = 0

            if 'xi_u' in ncfile.variables[varname].dimensions: 
                self.imin = 1
                self.iper = iper
            if 'eta_v' in ncfile.variables[varname].dimensions: 
                self.jmin = 1
                self.jper = jper
            if 's_w' in ncfile.variables[varname].dimensions: self.kmin = 0



            ####################################################################
            #1.1 It is a 2D variable
            ####################################################################
            if 's_rho' not in ncfile.variables[varname].dimensions and 's_w' not in ncfile.variables[varname].dimensions:
                #this is a 2d variable -> just load it
                self.data = np.squeeze( simul.Forder( ncfile.variables[varname][simul.infiletime,ny1+self.jper:ny2-self.jmin+self.jper,nx1+self.iper:nx2-self.imin+self.iper]) )


            ####################################################################
            #1.2 It is a 3D variable on one or more sigmal-levels
            ####################################################################
            elif np.min(depths)>=0:
                
                if debug: print('3D variable on one or more sigmal-levels')
                if debug: print(depths)
                
                if 's_w' in ncfile.variables[varname].dimensions:

                    if (len(depths)==1):   depth = depths[0] -1 
                    else:   depth = np.append(np.array(depths)-1,np.max(depths))
                    if debug: print(depth)
                    self.data = np.squeeze(simul.Forder(ncfile.variables[varname][simul.infiletime,depth,ny1+self.jper:ny2-self.jmin+self.jper,nx1+self.iper:nx2-self.imin+self.iper]))

                elif 's_rho' in ncfile.variables[varname].dimensions:

                    if (len(depths)==1):   depth = depths[0] -1 
                    else:   depth = np.array(depths) - 1
                    self.data = np.squeeze(simul.Forder(ncfile.variables[varname][simul.infiletime,depth,ny1+self.jper:ny2-self.jmin+self.jper,nx1+self.iper:nx2-self.imin+self.iper]))
                    
                if debug: print('var shape', self.data.shape)
                if debug: print('simul.infiletime,depth,ny1,ny2-self.jmin,nx1,nx2-self.imin')
                if debug: print(simul.infiletime,depth,ny1,ny2-self.jmin,nx1,nx2-self.imin)
                
            ####################################################################
            #1.3 It is a 3D variable on z-levels (interpolation needed)
            ####################################################################
            else:
                #if max(depths)>0: raise NameError('Depths are ill-defined. Check again please.')

                #Check how large is the domain___ Divide computations in chunk if $n_x*n_y > n2max$ 
                nchunk = int(np.max([np.sqrt((ny2-ny1)*(nx2-nx1)/n2max),1]))

                #depths can be a scalar or a 2d/3d array.
                if np.ndim(depths)==1:
                    zsize = len(depths)
                elif np.ndim(depths)==2:
                    zsize = 1
                elif np.ndim(depths)==3:
                    zsize = depths.shape[2]
                    

                self.data = np.zeros((nx2-nx1-self.imin,ny2-ny1-self.jmin,np.max([zsize,1])))*np.nan

                for i,j in product(list(range(nchunk)),list(range(nchunk))):
            
                    dx1=4; dx2=4; dy1=4; dy2=4; imin=0; jmin=0
                    if i==0: dx1=0
                    if i==nchunk-1: dx2=0; imin=self.imin
                    if j==0: dy1=0 
                    if j==nchunk-1: dy2=0; jmin=self.jmin

                    nx1i = int(nx1+i*(nx2-nx1)/nchunk-2*dx1)
                    nx2i = int(nx1+(i+1)*(nx2-nx1)/nchunk+2*dx2)
                    ny1i = int(ny1+j*(ny2-ny1)/nchunk-2*dy1)
                    ny2i = int(ny1+(j+1)*(ny2-ny1)/nchunk+2*dy2)

                    #we need to perform some vertical interpolation_ compute z_r,z_w for subdomain only
                    [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

                    # Check of depths is a 2-d array:
                    if np.ndim(depths)==2:
                        subdepths = depths[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1]
                    elif np.ndim(depths)==3:
                        subdepths = depths[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1,:]
                    else:
                        subdepths = depths
                    
                    #Compute variable in subdomain
                    chunk = tools.vinterp(simul.Forder(ncfile.variables[varname][simul.infiletime,:,ny1i+self.jper:ny2i-self.jmin+self.jper,nx1i+self.iper:nx2i-self.imin+self.iper]),\
                            subdepths,z_r,z_w,imin=self.imin,jmin=self.jmin,kmin=self.kmin,floattype = simul.floattype)

                    #Include in full variable                
                    self.data[nx1i+dx1-nx1:nx2i-dx2-nx1-imin,ny1i+dy1-ny1:ny2i-dy2-ny1-jmin,:] = \
                             chunk[dx1:nx2i-nx1i-dx2-imin,dy1:ny2i-ny1i-dy2-jmin,:]

                    #print 'chunck',i,j,nx1i+dx1-nx1,nx2i-dx2-nx1-self.imin,ny1i+dy1-ny1,ny2i-dy2-ny1-self.jmin
                if zsize==1: self.data=self.data[:,:,0]


        ####################################################################


            try: 
                self.data[self.data==ncfile.variables[varname]._FillValue] = masked_value
            except:
                try: 
                    self.data[self.data==ncfile.variables[varname].fill_value] = masked_value          
                except:       
                    try:
                        self.data[self.data==ncfile.variables['vbar']._FillValue] = masked_value
                    except:              
                        #print 'no FillValue in file for', varname
                        pass
                
            
        ####################################################################
        #2. var is not in netcdf output but is defined in dictionnary
        ####################################################################
        elif (self.longname != 'unknown') and (method=='new'):

            if debug: print('Variable not in ROMS outputs _ will be computed using new Fortran tools')

            ####################################################################
            # Create an array to store the variable (self.data)
            ####################################################################
            if debug: print(varname, self.kmin, len(depths))

            #depths can be a scalar or a 2d/3d array.
            if np.ndim(depths)==1:
                zsize = len(depths)
            elif np.ndim(depths)==2:
                zsize = 1
            elif np.ndim(depths)==3:
                zsize = depths.shape[2]
                    
            if (self.kmin>=0) and (np.min(depths)>0) and (zsize>1):
                self.data = np.zeros((nx2-nx1-self.imin,ny2-ny1-self.jmin,np.max([zsize+1-self.kmin,1])))*np.nan
            elif  (self.kmin>=0) and (zsize>1):
                self.data = np.zeros((nx2-nx1-self.imin,ny2-ny1-self.jmin,np.max([zsize,1])))*np.nan
            else:
                self.data = np.zeros((nx2-nx1-self.imin,ny2-ny1-self.jmin))*np.nan

            #print 'self.data', self.data.shape
            # Number of tiles used for computation
            nchunk = int(np.max([np.sqrt((ny2-ny1)*(nx2-nx1)/n2max),1]))

            #You cannot compute psi in chunks: 
            if self.name in ['psi','psi_v2','psir','psi_surf','int_psi','meridional_psi']: nchunk=1

            if debug: print('Domain will be divided in ', nchunk**2 , ' chunks')

            for i,j in product(list(range(nchunk)),list(range(nchunk))):
                
                dx1=4; dx2=4; dy1=4; dy2=4; # extra pts around the tile
                if i==0: dx1=0
                if i==nchunk-1: dx2=0 
                if j==0: dy1=0 
                if j==nchunk-1: dy2=0

                nx1i = int(nx1+i*(nx2-nx1)/nchunk-dx1)
                nx2i = int(nx1+(i+1)*(nx2-nx1)/nchunk+dx2)
                ny1i = int(ny1+j*(ny2-ny1)/nchunk-dy1)
                ny2i = int(ny1+(j+1)*(ny2-ny1)/nchunk+dy2)
                
                #print '[ny1i,ny2i,nx1i,nx2i]',[ny1i,ny2i,nx1i,nx2i]

                ####################################################################
                # Compute variable
                ####################################################################   


                if (self.name[:4] in ['int_']):
                    #Variables with names beginning by int_ are variables vertically integrated between 2 depths (so 2D variables but zsize!=1)
                    if debug: print('depths for integration is', depths)
                    chunk = self.get_sig(ncfile,simul,depths=depths,coord=[ny1i,ny2i,nx1i,nx2i],subcoord=[ny1i-ny1,ny2i-ny1,nx1i-nx1,nx2i-nx1],debug=debug)
                    
                elif (np.min(depths)>0) or (self.kmin<0):
                    chunk = self.get_sig(ncfile,simul,depths=depths,coord=[ny1i,ny2i,nx1i,nx2i],subcoord=[ny1i-ny1,ny2i-ny1,nx1i-nx1,nx2i-nx1],debug=debug)

                elif (zsize==1) and (np.min(depths)==0) and (self.name in ['vrt','str','absvrt','ow',\
                                                                           'dxbuoy2','tendency']):
                    
                    chunk = self.get_sig(ncfile,simul,depths=depths,coord=[ny1i,ny2i,nx1i,nx2i],subcoord=[ny1i-ny1,ny2i-ny1,nx1i-nx1,nx2i-nx1],debug=debug)

                else:

                    [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])
                    
                    # Check of depths is a 2-d array:
                    if np.ndim(depths)==2:
                        subdepths = depths[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1]
                    elif np.ndim(depths)==3:
                        subdepths = depths[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1,:]
                    else:
                        subdepths = depths
                    
                    
                    #Compute variable in subdomain
                    chunk = tools.vinterp(self.get_sig(ncfile,simul,depths=simul.coordmax[4],\
                                                        coord=[ny1i,ny2i,nx1i,nx2i],\
                                                        subcoord=[ny1i-ny1,ny2i-ny1,nx1i-nx1,nx2i-nx1]),\
                                                        subdepths,\
                                                        z_r,z_w,\
                                                        imin=self.imin,jmin=self.jmin,kmin=self.kmin,\
                                                        floattype = simul.floattype,debug=debug)
                                        
                ####################################################################
                # Write the chunk into the self.data
                ####################################################################


                if (self.kmin>=0) and (zsize>1):
                    # 3D variable
                    self.data[nx1i+dx1-nx1:nx2i-dx2-nx1-self.imin+1,ny1i+dy1-ny1:ny2i-dy2-ny1-self.jmin+1,:] = \
                        chunk[dx1:nx2i-nx1i-dx2-self.imin+1,dy1:ny2i-ny1i-dy2-self.jmin+1,:]                   
                else:
                    # 2D variable                 
                    self.data[nx1i+dx1-nx1:nx2i-dx2-nx1-self.imin+1,ny1i+dy1-ny1:ny2i-dy2-ny1-self.jmin+1] = \
                         np.squeeze(chunk)[dx1:nx2i-nx1i-dx2-self.imin+1,dy1:ny2i-ny1i-dy2-self.jmin+1]


            if zsize==1 and len(self.data.shape)>=3: self.data=self.data[:,:,0]

        ####################################################################

        else:
            
            if debug: print('We will use older script version')
            self.oldvar(varname,simul,depths = depths, u=u,v=v)
            
            #raise NameError('Sorry. I don t know how to compute '  + varname + '.')


        ###################################################################################



        ncfile.close()





###################################################################################
#long name for each variable (used as plot title)
###################################################################################



    def dictionnary(self):
       
         [self.longname, self.unit, [self.imin, self.jmin, self.kmin]] =\
                         self.dico.get(self.name,['unknown','unknown',[0,0,1]])


###################################################################################
#Load variables
###################################################################################


    def load(self,varname,ncfile,simul,debug=False,**kwargs):

        if debug: print(varname)
        [ny1,ny2,nx1,nx2,depths] = self.coord

        if 'coord' in  kwargs: [ny1,ny2,nx1,nx2] = kwargs['coord'][0:4]
        if 'depths' in  kwargs: depths = kwargs['depths']

        [imin,jmin,kmin] = self.dico.get(varname)[2]; depth = np.array(depths)-1
        if len(depth)==1: depth = depth[0]

        if 'masked' in  kwargs: 
            masked_value =  kwargs['masked']
        else: 
            masked_value =  np.nan
            
        ############################# 
        # test if periodic or tile
        #############################
        xperiodic = False; iper=0*imin
        yperiodic = False; jper=0*jmin

        try:
            if ncfile.dimensions['xi_u'].size==ncfile.dimensions['xi_rho'].size:
                xperiodic = True; iper=1*imin
        except:
            pass

        try:
            if ncfile.dimensions['eta_v'].size==ncfile.dimensions['eta_rho'].size:
                yperiodic = True; jper=1*jmin
        except:
            pass


        #############################
        if debug: print(varname, ncfile.variables[varname].shape)

        try:
            data = np.squeeze(simul.Forder(ncfile.variables[varname][simul.infiletime,depth,ny1+jper:ny2-jmin+jper,nx1+iper:nx2-imin+iper]))
        except:
            data = np.squeeze(simul.Forder(ncfile.variables[varname][simul.infiletime,ny1+jper:ny2-jmin+jper,nx1+iper:nx2-imin+iper]))



        try: 
            data[data==ncfile.variables[varname]._FillValue] = masked_value
        except:
            try: 
                data[data==ncfile.variables[varname].fill_value] = masked_value        
            except:       
                try:
                    data[self.data==ncfile.variables['vbar']._FillValue] = masked_value
                except:              
                    #print 'no FillValue in file for', varname
                    pass
                    
        data[data>1e30] = masked_value



        return data

    

###################################################################################
#Compute variables on sigma-levels
###################################################################################


    def get_sig(self,ncfile,simul,debug=False,**kwargs):


        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = self.coord[0:4]
            [ny1s,ny2s,nx1s,nx2s] = simul.coord[0:4]          
            pm = np.asfortranarray(simul.pm[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])
            pn = np.asfortranarray(simul.pn[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])
            f = np.asfortranarray(simul.f[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])
            mask = np.asfortranarray(simul.mask[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])
            topo = np.asfortranarray(simul.topo[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])            
        else: 
            coord = self.coord[0:4]
            [ny1,ny2,nx1,nx2] = self.coord[0:4]
            [ny1s,ny2s,nx1s,nx2s] = simul.coord[0:4]    
            pm = np.asfortranarray(simul.pm[nx1s:nx2-nx1s,ny1s:ny2-ny1s])
            pn = np.asfortranarray(simul.pn[nx1s:nx2-nx1s,ny1s:ny2-ny1s])
            f = np.asfortranarray(simul.f[nx1s:nx2-nx1s,ny1s:ny2-ny1s])
            mask = np.asfortranarray(simul.mask[nx1s:nx2-nx1s,ny1s:ny2-ny1s])
            topo = np.asfortranarray(simul.topo[nx1s:nx2-nx1s,ny1s:ny2-ny1s])
            
        if 'depths' in  kwargs: 
            depths = kwargs['depths']
        else: 
            depths = self.coord[4]
        

        #mask for fortran routines
        rmask = copy(mask)
        rmask[np.isnan(mask)] = 0

        # We use all sigma levels
        #depths = simul.coordmax[4]

        ################################################


        if self.name in ['rho','rho1','rhop','bvf','buoy','buoy1']:


            if self.name not in ['rho1','buoy1']:
                [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            T = self.load('temp',ncfile,simul,coord=coord, depths=depths)
            try:
                S = self.load('salt',ncfile,simul,coord=coord, depths=depths)
            except:
                print('no S in file')
                S = T*0.
              
            if self.name in ['rho']: var = toolsF.rho_eos(T,S,z_r,z_w,simul.rho0)
            elif self.name in ['rho1']: var = toolsF.rho1_eos(T,S,simul.rho0)    
            elif self.name in ['bvf']: var = toolsF.bvf_eos(T,S,z_r,z_w,simul.rho0)   
            elif self.name in ['buoy']: var = toolsF.get_buoy(T,S,z_r,z_w,simul.rho0)
            elif self.name in ['buoy1']: var = toolsF.rho1_eos(T,S,simul.rho0)*(-simul.g/simul.rho0)
            elif self.name in ['rhop']: var = tools_g.rhop(T,S)

        ################################################
        
        elif self.name in ['hbls_rho']:
            '''
            mixed-layer depth computed as depth where rho>rho_surface+0.03
            '''
            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            T = self.load('temp',ncfile,simul,coord=coord, depths=simul.coordmax[4])
            try:
                S = self.load('salt',ncfile,simul,coord=coord, depths=simul.coordmax[4])
            except:
                print('no S in file')
                S = T*0.
              
            rho = toolsF.rho_eos(T,S,z_r,z_w,simul.rho0)
            var = copy(rho[:,:,-1])*np.nan

            for i in range(rho.shape[0]):
                for j in range(rho.shape[1]):
                    try:
                        var[i,j] = np.min([-z_r[i,j,np.min([1+np.nanargmax(np.where((rho[i,j,:] - rho[i,j,-1])>0.03)),rho.shape[2]-1])],topo[i,j]])
                    except:
                        var[i,j] = topo[i,j]

        ################################################
        
        elif self.name in ['hbls_akv']:
            '''
            Mixed-layer depth based on AKv (Akv > 1e-4)
            '''
            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            try:
                AKv = self.load('AKv',ncfile,simul,coord=coord, depths=simul.coordmax[4])
            except:
                AKv = self.load('AKt',ncfile,simul,coord=coord, depths=simul.coordmax[4])

            var = copy(z_r[:,:,-1])*np.nan

            for i in range(var.shape[0]):
                for j in range(var.shape[1]):
                    try:
                        k_hbl = np.min([np.nanmax((AKv[i,j,:]<2e-4).nonzero()),z_w.shape[2]-1])
                        var[i,j] = np.min([-z_w[i,j,k_hbl],topo[i,j]])
                    except:
                        var[i,j] = topo[i,j]

        ################################################
        
        elif self.name in ['Ekman']:
            '''
            Ekman number based on AKv and Active mixing ML (Akv > 1e-4)
            '''
            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            try:
                AKv = self.load('AKv',ncfile,simul,coord=coord, depths=simul.coordmax[4])
            except:
                AKv = self.load('AKt',ncfile,simul,coord=coord, depths=simul.coordmax[4])
            var = copy(z_r[:,:,-1])*np.nan

            for i in range(var.shape[0]):
                for j in range(var.shape[1]):
                    k_hbl = np.min([np.nanmax((AKv[i,j,:]<2e-4).nonzero()),z_w.shape[2]-1])
                    hbl = np.min([-z_w[i,j,k_hbl],topo[i,j]])
                    var[i,j] = np.nanmean(AKv[i,j,k_hbl:]) / (f[i,j] * hbl**2)

         ################################################
        
        elif self.name in ['hbls_t']:
            '''
            mixed-layer depth computed as depth where temp<temp_surface-0.2
            '''
            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            T = self.load('temp',ncfile,simul,coord=coord, depths=simul.coordmax[4])
            
            var = copy(T[:,:,-1])*np.nan

            for i in range(T.shape[0]):
                for j in range(T.shape[1]):
                    try:
                        ind = 1+np.nanargmax(np.where(np.abs((T[i,j,-1] - T[i,j,:]))>0.2))
                    except:
                        ind = 0
                    var[i,j] = np.min([-z_r[i,j,np.min([ind,T.shape[2]-1])],topo[i,j]])
            

         ################################################
        
        elif self.name in ['hbbl']:
            '''
            mixed-layer depth computed as depth where temp<temp_surface-0.2
            '''

            print('evaluating hbbl using AKt')

            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            AKt = self.load('AKt',ncfile,simul,coord=coord,depths=simul.coordmax[4])
            var = toolsF_g.get_hbbls_from_akt(AKt,z_w)

         ################################################

        elif self.name in ['psi']:

            u = self.load('ubar',ncfile,simul,coord=coord)
            v = self.load('vbar',ncfile,simul,coord=coord)
    
            u[np.isnan(u)]=0.
            v[np.isnan(v)]=0.   
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); #Hz[Hz==0] = np.nan
            del z_r, z_w

            u = 0.5*(Hz[1:,:] + Hz[:-1,:])* u * 2 /(pn[1:,:]+pn[:-1,:])
            v = 0.5*(Hz[:,1:] + Hz[:,:-1])* v * 2 /(pm[:,1:]+pm[:,:-1])
            
            lx,ly = pm.shape[0]-1,pm.shape[1]-1

            psi = np.zeros((lx,ly))*np.nan

            psi[0,:] = -1*np.cumsum(u[0,:-1])
            psi[:,0] = np.cumsum(v[:-1,0])

            for i in range(1,lx):
                for j in range(1,ly):
                    psi[i,j]= psi[i-1,j-1] + (v[i,j-1] - u[i,j])

            var = psi*tools.rho2psi(mask)

         ################################################
        

        elif self.name in ['psi_v2']:

            u = self.load('ubar',ncfile,simul,coord=coord)
            v = self.load('vbar',ncfile,simul,coord=coord)
    
            u[np.isnan(u)]=0.
            v[np.isnan(v)]=0.   
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); #Hz[Hz==0] = np.nan
            del z_r, z_w

            u = 0.5*(Hz[1:,:] + Hz[:-1,:])* u * 2 /(pn[1:,:]+pn[:-1,:])
            v = 0.5*(Hz[:,1:] + Hz[:,:-1])* v * 2 /(pm[:,1:]+pm[:,:-1])
            
            lx,ly = pm.shape[0]-1,pm.shape[1]-1

            psi = np.zeros((lx,ly))*np.nan

            psi[lx-1,:] = -1*np.cumsum(u[lx-1,1:])
            psi[:,ly-1] = np.cumsum(v[1:,ly-1])

            for i in range(1,lx):
                for j in range(1,ly):
                    psi[lx-1-i,ly-1-j]= psi[lx-i,ly-j] - (v[lx-i,ly-1-j] - u[lx-i,ly-j])

            var = psi*tools.rho2psi(mask)


        ################################################
        

        elif self.name in ['psir']:

            u0 = self.load('u',ncfile,simul,depths=[0],coord=coord)
            v0 = self.load('v',ncfile,simul,depths=[0],coord=coord)

            ud,vd = tools_g.div2uvs(u0,v0,pm,pn)
            u,v = u0-ud, v0-vd
            del u0,v0,ud,vd

            u[np.isnan(u)]=0.
            v[np.isnan(v)]=0.   
            
            lx,ly = pm.shape[0]-1,pm.shape[1]-1

            psi = np.zeros((lx,ly))*np.nan

            psi[0,:] = -1*np.cumsum(u[0,:-1]/(0.5*(pn[0,:-1]+pn[1,:-1])))
            psi[:,0] = np.cumsum(v[:-1,0]/(0.5*(pm[:-1,0]+pm[:-1,1])))

            for i in range(1,lx):
                for j in range(1,ly):
                    psi[i,j]= psi[i-1,j-1] + (v[i,j-1]/(0.5*(pm[i,j]+pm[i,j-1])) - u[i,j]/(0.5*(pn[i,j]+pn[i-1,j])))

            var = psi*tools.rho2psi(mask)* tools.rho2psi(simul.f)/simul.g
            
         
        ################################################

        elif self.name in ['zonal_psi']:

            u = self.load('ubar',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); #Hz[Hz==0] = np.nan
            del z_r, z_w

            u = 0.5*(Hz[1:,:] + Hz[:-1,:])* u * 2 /(pn[1:,:]+pn[:-1,:])
            
            lx,ly = pm.shape[0]-1,pm.shape[1]-1

            psi = np.zeros((lx,ly))*np.nan
            psi[:,ly-1] = 0.
            
            for j in range(1,ly):
                psi[:,ly-1-j]= psi[:,ly-j] + u[:,ly-j]

            var = psi*tools.rho2psi(mask)/1e6

        ################################################

        elif self.name in ['meridional_psi']:

            v = self.load('vbar',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); #Hz[Hz==0] = np.nan
            del z_r, z_w

            v = 0.5*(Hz[:,1:] + Hz[:,:-1])* v * 2 /(pm[:,1:]+pm[:,:-1])
            
            lx,ly = pm.shape[0]-1,pm.shape[1]-1

            psi = np.zeros((lx,ly))*np.nan
            psi[lx-1,:] = 0.
            
            for i in range(1,lx):
                psi[lx-1-i,:]= psi[lx-i,:] +v[lx-i,:]

            var = psi*tools.rho2psi(mask)/1e6


        ################################################
        

        elif self.name in ['int_psi']:

            u = self.load('u',ncfile,simul,coord=coord, depths=simul.coordmax[4])
            v = self.load('v',ncfile,simul,coord=coord, depths=simul.coordmax[4])

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            
            u = tools.vert_int(u, 0.5*(z_w[1:,:,:] + z_w[:-1,:,:]),depths[0],depths[1])* 2 /(pn[1:,:]+pn[:-1,:])
            v = tools.vert_int(v, 0.5*(z_w[:,1:,:] + z_w[:,:-1,:]),depths[0],depths[1])* 2 /(pm[:,1:]+pm[:,:-1])

            
            lx,ly = pm.shape[0]-1,pm.shape[1]-1
            psi = np.zeros((lx,ly))*np.nan

            psi[0,:] = -1*np.cumsum(u[0,:-1])
            psi[:,0] = np.cumsum(v[:-1,0])

            for i in range(1,lx):
                for j in range(1,ly):
                    psi[i,j]= psi[i-1,j-1] + (v[i,j-1] - u[i,j])

            var = psi*tools.rho2psi(mask)
            
                        
         ################################################
        

        elif self.name in ['psi_surf']:

            u = self.load('u',ncfile,simul,depths=[0],coord=coord)
            v = self.load('v',ncfile,simul,depths=[0],coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,-2]); #Hz[Hz==0] = np.nan
            del z_r, z_w

            u = 0.5*(Hz[1:,:] + Hz[:-1,:])* u * 2 /(pn[1:,:]+pn[:-1,:])
            v = 0.5*(Hz[:,1:] + Hz[:,:-1])* v * 2/(pm[:,1:]+pm[:,:-1])
            
            lx,ly = pm.shape[0]-1,pm.shape[1]-1

            psi = np.zeros((lx,ly))*np.nan

            psi[0,:] = -1*np.cumsum(u[0,:-1])
            psi[:,0] = np.cumsum(v[:-1,0])

            for i in range(1,lx):
                for j in range(1,ly):
                    psi[i,j]= psi[i-1,j-1] + (v[i,j-1] - u[i,j])

            var = psi*tools.rho2psi(mask)
                       
            
        ################################################


        elif self.name in ['w']:
            
            if debug: print('computing w', coord)
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var= tools.nanbnd(toolsF.get_wvlcty(u,v,z_r,z_w,pm,pn))
            #var= toolsF.get_wvlcty(u,v,z_r,z_w,pm,pn)

        ################################################

        
        elif self.name in ['w_vmix']:
            
            depths_r = simul.coordmax[4]
            depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            T = self.load('temp',ncfile,simul,coord=coord, depths=depths_r)
            S = self.load('salt',ncfile,simul,coord=coord, depths=depths_r)
            
            #buoy = toolsF.get_buoy(T,S,z_r,z_w,simul.rho0)
            buoy = toolsF.rho1_eos(T,S,simul.rho0)*(-simul.g/simul.rho0)
            del T,S
            
            AKv = self.load('AKv',ncfile,simul,coord=coord, depths=depths_w)
            AKv = 0.5*(AKv[:,:,1:] + AKv[:,:,:-1])

            ################################################
            
            byy = tools.diffeta(tools.rho2v(AKv) * tools.diffeta(buoy,pn,z_r,z_w),tools.rho2v(pn),tools.rho2v(z_r),tools.rho2v(z_w))
            bxx = tools.diffxi(tools.rho2u(AKv) * tools.diffxi(buoy,pm,z_r,z_w),tools.rho2u(pm),tools.rho2u(z_r),tools.rho2u(z_w))
            
            ################################################ 
            
            var = np.zeros(buoy.shape)*np.nan
            var[1:-1,:,:] = bxx 
            var[:,1:-1,:] = var[:,1:-1,:] + byy        
            var = -1*(var.T/(simul.f**2).T).T
            
            return var

        ################################################

        
        elif self.name in ['w_ttw']:
            
            var = self.get_w_ttw_sig(simul,coord=coord)
            
            return var
            
        ################################################
        
        elif self.name in ['absvrt']:  


            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            #this is the python version
            var = tools.get_absvrt(u,v,z_r,z_w,f,pm,pn)
            
            
        ################################################
        
        elif self.name in ['vrt']:  

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            if debug: print(u.shape,depths)
            
            if len(depths)==1:
                z_r,z_w = None,None
            else:
                [z_r,z_w] = tools.get_depths(simul,coord=coord)

            #this is the python version
            var = tools.get_vrt(u,v,pm,pn,z_r,z_w)
            
        ################################################
        
        elif self.name in ['str']:

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            if debug: print(u.shape,depths)
            
            if len(depths)==1:
                z_r,z_w = None,None
            else:
                [z_r,z_w] = tools.get_depths(simul,coord=coord)

            #this is the python version
            var = tools.get_strain(u,v,pm,pn,z_r,z_w)
            
        ################################################
        
        elif self.name in ['ow']:


            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            if debug: print(u.shape,depths)
            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            #this is the python version
            var = tools.get_ow(u,v,z_r,z_w,pm,pn)
            
        ################################################
        
        elif self.name in ['Ri']:

            
            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
            u = tools.u2rho(self.load('u',ncfile,simul,coord=coord,depths=depths))
            v = tools.v2rho(self.load('v',ncfile,simul,coord=coord,depths=depths))

            [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)

            bvf = toolsF.bvf_eos(T,S,z_r,z_w,simul.rho0)
            
            dzu = (u[:,:,1:] - u[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1])
            dzv = (v[:,:,1:] - v[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1])

            var = bvf * 0.
            var[:,:,1:-1] = bvf[:,:,1:-1]/(dzu**2 + dzv**2);
            
        ################################################
        
        elif self.name in ['dzu']:

            u = tools.u2rho(self.load('u',ncfile,simul,coord=coord,depths=depths))

            [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)
            
            dzu = (u[:,:,1:] - u[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1])

            var = z_w * np.nan
            var[:,:,1:-1] = dzu
            
        ################################################
        
        elif self.name in ['dzv']:

            v = tools.v2rho(self.load('v',ncfile,simul,coord=coord,depths=depths))

            [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)

            dzv = (v[:,:,1:] - v[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1])

            var = z_w * np.nan
            var[:,:,1:-1] = dzv
            
        ################################################
        
        elif self.name in ['dzu2']:

            u = tools.u2rho(self.load('u',ncfile,simul,coord=coord,depths=depths))
            v = tools.v2rho(self.load('v',ncfile,simul,coord=coord,depths=depths))

            [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)

            
            dzu = (u[:,:,1:] - u[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1])
            dzv = (v[:,:,1:] - v[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1])

            var = z_w * np.nan
            var[:,:,1:-1] = np.sqrt(dzu**2 + dzv**2);
            
        ################################################
        
        elif self.name in ['kdzu']:
        
            depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

            u = tools.u2rho(self.load('u',ncfile,simul,coord=coord,depths=depths))

            [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)

            dzu = (u[:,:,1:] - u[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1])
            
            AKv = self.load('AKv',ncfile,simul,coord=coord, depths=depths_w)
            
            var = z_w * np.nan
            var[:,:,1:-1] = AKv[:,:,1:-1]*dzu
            
        ################################################
            
        elif self.name in ['kdzv']:
        
            depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

            v = tools.v2rho(self.load('v',ncfile,simul,coord=coord,depths=depths))

            [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)

            
            dzv = (v[:,:,1:] - v[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1])
            
            AKv = self.load('AKv',ncfile,simul,coord=coord, depths=depths_w)
            
            var = z_w * np.nan
            var[:,:,1:-1] = AKv[:,:,1:-1]* dzv
            
        ################################################
            
        elif self.name in ['kdzu2']:

            depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

            u = tools.u2rho(self.load('u',ncfile,simul,coord=coord,depths=depths))
            v = tools.v2rho(self.load('v',ncfile,simul,coord=coord,depths=depths))

            [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)

            
            dzu = (u[:,:,1:] - u[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1])
            dzv = (v[:,:,1:] - v[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1])
            
            AKv = self.load('AKv',ncfile,simul,coord=coord, depths=depths_w)
            
            var = z_w * np.nan
            var[:,:,1:-1] = AKv[:,:,1:-1]*(dzu**2 + dzv**2);
            
        ################################################
        elif self.name in ['pv']:  

            
            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)

            if debug: print('for PV using levels', depths)
            if debug: print(T.shape, S.shape)
            if debug: print(z_r.shape, z_w.shape)


            #this is the python version
            var = tools.PV(T,S,u,v,z_r,z_w,f,simul.g,simul.rho0,pm,pn)

        ################################################
        elif self.name in ['pvr']:


            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)

            if debug: print('for PV using levels', depths)
            if debug: print(T.shape, S.shape)
            if debug: print(z_r.shape, z_w.shape)


            #this is the python version
            var = tools.PVr(T,S,u,v,z_r,z_w,f,simul.g,simul.rho0,pm,pn)

            
        ################################################
        elif self.name in ['stretching']:  

            
            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            #this is the python version
            var = tools.stretching(T,S,z_r,z_w,f,simul.g,simul.rho0,pm,pn)
            
        ################################################
        elif self.name in ['dxbuoy2']:
            
            if debug: print('-------------------------------------------------')
            if debug:  print('Computing dxbuoy2 using new R_tools_gula version')
            if debug: print('-------------------------------------------------')
            

            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
            
            if len(T.shape)==2:
                if debug:  print('Computing dxbuoy2 at the surface')
                buoy = toolsF.rho1_eos_2d(T,S,simul.rho0)*(-simul.g/simul.rho0)
                bx = tools.u2rho(tools.diffx(buoy,pm))
                by = tools.v2rho(tools.diffy(buoy,pn))
            else:
                buoy = toolsF.rho1_eos(T,S,simul.rho0)*(-simul.g/simul.rho0)
                [z_r,z_w] = tools.get_depths(simul,coord=coord)
                bx = tools.u2rho(tools.diffxi(buoy,pm,z_r,z_w))
                by = tools.v2rho(tools.diffeta(buoy,pn,z_r,z_w))
            del buoy,T,S

            var = 0.5 * (bx**2 + by**2)
        
        ################################################
        elif self.name in ['dxu2']:
            if debug: print('-------------------------------------------------')
            if debug: print('Computing dxu2 using new R_tools_gula version')
            if debug: print('-------------------------------------------------')
            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            ux = tools.u2rho(tools.diffxi(tools.u2rho(u),pm,z_r,z_w))
            uy = tools.v2rho(tools.diffeta(tools.u2rho(u),pn,z_r,z_w))
            vx = tools.u2rho(tools.diffxi(tools.v2rho(v),pm,z_r,z_w))
            vy = tools.v2rho(tools.diffeta(tools.v2rho(v),pn,z_r,z_w))

            var = 0.5 * (ux**2 + uy**2 + vx**2 + vy**2)
        
        ################################################
        elif self.name in ['tendency','tendency_3d']:

            if debug: print('-------------------------------------------------')
            if debug: print('Computing tendency using new R_tools_gula version')
            if debug: print('-------------------------------------------------')

            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            if len(T.shape)==2:
                if debug:  print('Computing tendency at the surface')
                ux = tools.u2rho(tools.diffx(tools.u2rho(u),pm))
                uy = tools.v2rho(tools.diffy(tools.u2rho(u),pn))
                vx = tools.u2rho(tools.diffx(tools.v2rho(v),pm))
                vy = tools.v2rho(tools.diffy(tools.v2rho(v),pn))
                del u,v

                buoy = toolsF.rho1_eos_2d(T,S,simul.rho0)*(-simul.g/simul.rho0)
                del T,S

                bx = tools.u2rho(tools.diffx(buoy,pm))
                by = tools.v2rho(tools.diffy(buoy,pn))
                del buoy

            else:
                [z_r,z_w] = tools.get_depths(simul,coord=coord)

                ux = tools.u2rho(tools.diffxi(tools.u2rho(u),pm,z_r,z_w))
                uy = tools.v2rho(tools.diffeta(tools.u2rho(u),pn,z_r,z_w))
                vx = tools.u2rho(tools.diffxi(tools.v2rho(v),pm,z_r,z_w))
                vy = tools.v2rho(tools.diffeta(tools.v2rho(v),pn,z_r,z_w))
                del u,v

                buoy = toolsF.rho1_eos(T,S,simul.rho0)*(-simul.g/simul.rho0)
                del T,S

                bx = tools.u2rho(tools.diffxi(buoy,pm,z_r,z_w))
                by = tools.v2rho(tools.diffeta(buoy,pn,z_r,z_w))
                del buoy

            var = -1*(bx * ux * bx + by * uy * bx + bx * vx * by + by * vy * by)
            del ux,uy,vx,vy
 
            ##############
            if self.name == 'tendency_3d':
                w = tools.nanbnd(toolsF.get_wvlcty(u,v,z_r,z_w,pm,pn))
                del u,v
     
                wx = tools.u2rho(tools.diffxi(w,pm,z_r,z_w))
                wy = tools.v2rho(tools.diffeta(w,pn,z_r,z_w))

                bz = tools.vinterp((buoy[:,:,1:] - buoy[:,:,:-1])/ (z_r[:,:,1:] -z_r[:,:,:-1]),z_r,z_w[:,:,1:-1],z_r)

                var += -1*(bz * wx * bx + bz * wy * by)


        ################################################
        elif self.name in ['tendency_u','tendency_3d_u']:
            if debug: print('-------------------------------------------------')
            if debug: print('Computing tendency using new R_tools_gula version')
            if debug: print('-------------------------------------------------')
            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            ux = tools.u2rho(tools.diffxi(tools.u2rho(u),pm,z_r,z_w))
            uy = tools.v2rho(tools.diffeta(tools.u2rho(u),pn,z_r,z_w))
            vx = tools.u2rho(tools.diffxi(tools.v2rho(v),pm,z_r,z_w))
            vy = tools.v2rho(tools.diffeta(tools.v2rho(v),pn,z_r,z_w))

            var = -1*(ux * ux * ux + uy * uy * ux + ux * vx * uy + uy * vy * uy)\
                  -1*(vx * ux * vx + vy * uy * vx + vx * vx * vy + vy * vy * vy)

            ##############
            if self.name == 'tendency_3d_u':
                w = tools.nanbnd(toolsF.get_wvlcty(u,v,z_r,z_w,pm,pn))
     
                wx = tools.u2rho(tools.diffxi(w,pm,z_r,z_w))
                wy = tools.v2rho(tools.diffeta(w,pn,z_r,z_w))

                uz = tools.vinterp(tools.u2rho(u[:,:,1:] - u[:,:,:-1])/ (z_r[:,:,1:] -z_r[:,:,:-1]),z_r,z_w[:,:,1:-1],z_r)
                vz = tools.vinterp(tools.v2rho(v[:,:,1:] - v[:,:,:-1])/ (z_r[:,:,1:] -z_r[:,:,:-1]),z_r,z_w[:,:,1:-1],z_r)
                
                var += -1*(uz * wx * ux + uz * wy * uy)
                var += -1*(vz * wx * vx + vz * wy * vy)


        ################################################
        elif self.name in ['tendency_full_u']:
        
            if debug: print('-------------------------------------------------')
            if debug: print('Computing full tendency using new R_tools_gula version')
            if debug: print('-------------------------------------------------')
            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            ux = tools.u2rho(tools.diffxi(tools.u2rho(u),pm,z_r,z_w))
            uy = tools.v2rho(tools.diffeta(tools.u2rho(u),pn,z_r,z_w))
            vx = tools.u2rho(tools.diffxi(tools.v2rho(v),pm,z_r,z_w))
            vy = tools.v2rho(tools.diffeta(tools.v2rho(v),pn,z_r,z_w))

            var = -1*(ux * ux * ux + uy * uy * ux + ux * vx * uy + uy * vy * uy)\
                  -1*(vx * ux * vx + vy * uy * vx + vx * vx * vy + vy * vy * vy)

            ################################################

            w = tools.nanbnd(toolsF.get_wvlcty(u,v,z_r,z_w,pm,pn))
     
            wx = tools.u2rho(tools.diffxi(w,pm,z_r,z_w))
            wy = tools.v2rho(tools.diffeta(w,pn,z_r,z_w))
            del w
            
            uz = tools.vinterp(tools.u2rho(u[:,:,1:] - u[:,:,:-1])/ (z_r[:,:,1:] -z_r[:,:,:-1]),z_r,z_w[:,:,1:-1],z_r)
            vz = tools.vinterp(tools.v2rho(v[:,:,1:] - v[:,:,:-1])/ (z_r[:,:,1:] -z_r[:,:,:-1]),z_r,z_w[:,:,1:-1],z_r)
                
            var += -1*(uz * wx * ux + uz * wy * uy)
            var += -1*(vz * wx * vx + vz * wy * vy)
            
            del uz,vz,wx,wy
            
            ################################################
            '''
            [MXadv,MYadv,MVadv,MHdiss,MHmix,MVmix,MCor,MPrsgrd] = self.get_uv_evolution(simul,coord=coord,rhs=True)
            del MXadv,MYadv,MVadv
            
            rhs = MHdiss+MHmix+MVmix+MCor+MPrsgrd
            del MHdiss,MHmix,MVmix,MCor,MPrsgrd
            '''
            
            rhs =  np.ones((z_r.shape[0],z_r.shape[1],z_r.shape[2],2))
            rhs[:] =  self.get_uv_evolution(simul,coord=coord,rhs=True)
            
            rhs_ux = tools.u2rho(tools.diffxi(tools.u2rho(rhs[1:,:,:,0]),pm,z_r,z_w))
            rhs_uy = tools.v2rho(tools.diffeta(tools.u2rho(rhs[1:,:,:,0]),pn,z_r,z_w))
            rhs_vx = tools.u2rho(tools.diffxi(tools.v2rho(rhs[:,1:,:,1]),pm,z_r,z_w))
            rhs_vy = tools.v2rho(tools.diffeta(tools.v2rho(rhs[:,1:,:,1]),pn,z_r,z_w))
            del rhs,z_r,z_w
            
            var += ux * rhs_ux + uy * rhs_uy + vx * rhs_vx + vy * rhs_vy
            del ux,rhs_ux,uy,rhs_uy,vx,rhs_vx,vy,rhs_vy

        ################################################
        elif self.name in ['pv1','pv2','pv3']:  

            
            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            #this is the python version
            [pv1,pv2,pv3] = tools_g.PV_terms(T,S,u,v,z_r,z_w,f,simul.g,simul.rho0,pm,pn)       
            
            
            if self.name in ['pv1']: var = pv1[:,:,1:-1]
            elif self.name in ['pv2']: var = pv2[:,:,1:-1]
            elif self.name in ['pv3']: var = pv3[:,:,1:-1]

        ################################################

        elif self.name in ['kediss']:

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = toolsF_g.get_kediss(u,v, z_r,z_w,pm,pn,rmask)
            

            
        ################################################
        # BAROTROPIC VORTICITY EQU. TERMS
        # vortbar on L.H.S   
        # all other terms are computed as if R.H.S (so beware of sign)
        ################################################


        elif self.name in ['vortbar_int']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = z_w[:,:,-1] - z_w[:,:,0]
            del z_r, z_w

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])*tools.rho2u(Hz)
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])*tools.rho2v(Hz)
            del Hz

            var = tools.rot(ubar,vbar,pm,pn) * tools.rho2psi(mask)
            del ubar,vbar

        ################################################
            

        elif self.name in ['intvrt']:

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            vrt = tools.get_vrt(u,v,pm,pn,z_r,z_w)
            del u,v
            
            Hz = tools.rho2psi(z_w[:,:,1:] - z_w[:,:,:-1])
            
            var = np.nansum(vrt*Hz,axis=2)
            
                        

        ################################################
        # updated 16/08/18

        elif self.name in ['bpt']:

            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = tools.nanbnd(toolsF_g.get_bpt(T,S, z_r,z_w,simul.rho0,pm,pn,rmask))

        ################################################


        elif self.name in ['bpts']:

            zeta = self.load('zeta',ncfile,simul,coord=coord)
            
            var = tools.jacob(zeta,topo,pm,pn)
            
        ################################################
        
   

        elif self.name in ['fwb']:     
            
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            var = toolsF_g.get_fwb(u,v,z_r,z_w,pm,pn,f)

            var= tools.nanbnd(var)
        
         ################################################
        
   

        elif self.name in ['fws']:     
            
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            var = toolsF_g.get_fws(u,v,z_r,z_w,pm,pn,f)

            var= tools.nanbnd(var)

        ################################################
   
        
             
   

        elif self.name in ['fdivub']:     
            
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            var = toolsF_g.get_fdivub(u,v,z_r,z_w,pm,pn,f)

            var= tools.nanbnd(var)
        
         ################################################    
       
            
   

        elif self.name in ['fwdivub']:     
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            var = toolsF_g.get_fwdivub(u,v,z_r,z_w,pm,pn,f)

            var= tools.nanbnd(var)
        
         ################################################    


        elif self.name in ['u_Prsgrd']:

            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = tools.nanbnd(toolsF_g.get_u_prsgrd(T,S, z_r,z_w,simul.rho0,pm,pn))

        ################################################

        elif self.name in ['vortplanet']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])
 
            var = -1*toolsF_g.get_vortplanet(ubar,vbar,Hz,pm,pn,f) #* tools.rho2psi(mask)

        ################################################

        elif self.name in ['vortstretch']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = -1*toolsF_g.get_vortstretch(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)
            
        ################################################

        #elif self.name in ['vortstretch2']:

            #[z_r,z_w] = tools.get_depths(simul,coord=coord)
            #Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            #del z_r, z_w

            #ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            #vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            #var = -1*toolsF.get_vortstretch2(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)            
            
        ################################################

        elif self.name in ['int_vortplanet']:
            '''
            int_vortplanet in vortplanet integrated from surface to depth
            
            '''
            if debug: print('integrating transport down to', depths[0])
            u = self.load('u',ncfile,simul,coord=coord,depths=simul.coordmax[4])
            v = self.load('v',ncfile,simul,coord=coord,depths=simul.coordmax[4])

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            if debug: print('integrate between:', depths[0],depths[1])
            var = -1*toolsF_g.get_intvortplanet(u,v, z_r,z_w,pm,pn,f,depths[0],depths[1]) * tools.rho2psi(mask)


        ################################################

        #elif self.name in ['vortstretch_sol1']:

            #[z_r,z_w] = tools.get_depths(simul,coord=coord)
            #Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            #del z_r, z_w

            #ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            #vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            #var = toolsF.get_vortstretch_sol1(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)
            
         ################################################

        elif self.name in ['vortstretch_sol2']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = toolsF_g.get_vortstretch_sol2(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)
            
         ################################################

        #elif self.name in ['vortstretch_sol3']:

            #[z_r,z_w] = tools.get_depths(simul,coord=coord)
            #Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            #del z_r, z_w

            #ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            #vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            #var = toolsF.get_vortstretch_sol3(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)


   
            
        ################################################
        # updated 16/08/18
        elif self.name in ['vortplantot_sol2']:

            #vortplantot_sol2 is Cor + CURVGRID

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = z_w[:,:,1:] - z_w[:,:,:-1]
            del z_r, z_w

            var = toolsF_g.get_vortplantot_sol2(u,v, Hz,pm,pn,f,rmask)

        ################################################

        elif self.name in ['vortplantot_sol3']:

            # vortplantot_sol3 is only Cor (should be same than vortplantot)
            # except computed with ROMS routines

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = toolsF_g.get_vortplantot_sol3(u,v, z_r,z_w,pm,pn,f) * tools.rho2psi(mask)

            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan   

        ################################################

        elif self.name in ['vortplantot']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = -1*toolsF_g.get_vortplantot(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)         


        ################################################


        elif self.name in ['vortadv']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,1:] - z_w[:,:,:-1]); Hz[Hz==0] = np.nan
            del z_r, z_w

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            u = tools.u2rho(u); v = tools.v2rho(v)

            u2 = tools.diffy(tools.diffx(np.nansum(Hz * u**2 ,axis=2),pm),tools.rho2u(pn))
            v2 = tools.diffy(tools.diffx(np.nansum(Hz * v**2 ,axis=2),pm),tools.rho2u(pn))

            uv = np.nansum(Hz * u * v ,axis=2)
            del u,v,Hz

            uvy = tools.v2rho(tools.rho2u(tools.diffy(tools.diffy(  uv     ,pn),tools.rho2v(pn))))
            uvx = tools.u2rho(tools.rho2v(tools.diffx(tools.diffx(  uv     ,pm),tools.rho2u(pm))))
            del uv

            var = -1*(v2 - u2 + uvx - uvy)
            del u2, v2, uvx, uvy

        ################################################


        elif self.name in ['int_vortadv_sol2']:

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)


            var = toolsF_g.get_int_adv_sol2(u,v, z_r,z_w,pm,pn) * tools.rho2psi(mask)
            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan      
     
        ################################################
        
        elif self.name in ['vortadv_sol1']:

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            depths  = self.coord[4]          
            depthsw = np.append(np.array(depths),np.max(depths)+1)
            w = self.load('omega',ncfile,simul,coord=coord,depths = depthsw)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)


            var = toolsF_g.get_adv_sol1(u,v,w, z_r,z_w,pm,pn) * tools.rho2psi(mask)
            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan


        ################################################
        # updated 16/08/18

        elif self.name in ['vortadv_sol2']:

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = tools.nanbnd(toolsF.get_adv_sol2(u,v, z_r,z_w,pm,pn,rmask),2)

        ################################################
        # updated 16/08/18

        elif self.name in ['vortadv_sol3']:

            # vortadv_sol3 is vortadv_sol2 + adv due to rotation of grid

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = toolsF_g.get_adv_sol2(u,v, z_r,z_w,pm,pn,rmask)
            var = tools.nanbnd(var + toolsF_g.get_uvgrid(u,v, z_r,z_w,pm,pn,f,rmask))

        
        ################################################
        # updated 16/08/18

        elif self.name in ['vortadv_mix']:

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = tools.nanbnd(toolsF.get_adv_mix(u,v, z_r,z_w,pm,pn,rmask),2)
            
        ################################################
        # updated 16/08/18

        elif self.name in ['vortadv_centered']:

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = tools.nanbnd(toolsF_g.get_adv_4th(u,v, z_r,z_w,pm,pn,rmask),2)
            
        ################################################        
        # updated 16/08/18
        
        elif self.name in ['vortadv_uvgrid']:

            # adv due to rotation of grid

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = tools.nanbnd(toolsF_g.get_uvgrid(u,v, z_r,z_w,pm,pn,f,rmask))
            
        ###################################################################################
        # updated 16/08/18

        elif self.name in ['rotwind','rotbot']:


            if self.name=='rotwind':
                
                [ut,vt] = self.get_winds(simul,coord=coord)

                umask = rmask[1:,:] * rmask[:-1,:]
                vmask = rmask[:,1:] * rmask[:,:-1]
                ut*=umask/simul.rho0; vt*=vmask/simul.rho0
                
                var = toolsF_g.get_rotv(ut,vt,pm,pn)


            elif self.name=='rotbot':
                
                u = self.load('u',ncfile,simul,coord=coord,depths=[1])
                v = self.load('v',ncfile,simul,coord=coord,depths=[1])
                [z_r,z_w] = tools.get_depths(simul,coord=coord)
                Hz = z_w[:,:,1] - z_w[:,:,0]
                del z_r,z_w

                (ut,vt) = toolsF_g.get_bot(u,v,Hz,simul.rdrg,simul.Zob)

                umask = rmask[1:,:] * rmask[:-1,:]
                vmask = rmask[:,1:] * rmask[:,:-1]
                ut*=umask; vt*=vmask

                var = -1*tools.nanbnd(toolsF.get_rotv(ut,vt,pm,pn))

        ###################################################################################

        ################################################
        # BAROTROPIC VORTICITY EQU. TERMS FOR DEPTH AVERAGED EQUATION
        # vortbar_mean on L.H.S   
        # all other terms are computed as if R.H.S (so beware of sign)
        ################################################

            
        elif self.name in ['vortbar_mean']:

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = tools.rot(ubar,vbar,pm,pn) * tools.rho2psi(mask)
            del ubar,vbar            
            
                        

        ################################################


        elif self.name in ['bpt_mean']:

            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = tools.nanbnd(toolsF_g.get_bpt_mean(T,S, z_r,z_w,simul.rho0,pm,pn) * tools.rho2psi(mask))

            
        ################################################

        elif self.name in ['vorttopo_mean']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w
            
            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])
 
            var = -1*toolsF_g.get_vortplanet(ubar,vbar,Hz,pm,pn,1/Hz) * tools.rho2psi(mask) * tools.rho2psi(f)
            
        ################################################

        elif self.name in ['vortf_mean']:


            Hz = np.ones(pm.shape)
            
            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])
 
            var = -1*toolsF_g.get_vortplanet(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)
            
        ################################################

        elif self.name in ['vortplanet_mean']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w
            
            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])
 
            var = -1*toolsF_g.get_vortplanet(ubar,vbar,Hz,pm,pn,f/Hz) * tools.rho2psi(mask)

        ################################################

        elif self.name in ['vortstretch_mean']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w
            
            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = -1*toolsF_g.get_vortstretch(ubar,vbar,Hz,pm,pn,f/Hz) * tools.rho2psi(mask)
            


         ################################################

        elif self.name in ['vortstretch_sol2_mean']:

            Hz = np.ones(pm.shape)

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = toolsF_g.get_vortstretch_sol2(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)
            

 
            
        ################################################

        elif self.name in ['vortplantot_sol2_mean']:

            #vortplantot_sol2 is Cor + adv dur to rotation of grid

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = toolsF_g.get_vortplantot_sol2_mean(u,v, z_r,z_w,pm,pn,f) * tools.rho2psi(mask)

            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan   



        ################################################

        elif self.name in ['vortplantot_mean']:

            Hz = np.ones(pm.shape)

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = -1*toolsF_g.get_vortplantot(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)         


        ################################################

        elif self.name in ['vortadv_sol2_mean']:
            ''' 
            In vortadv_sol2_mean we do the same computations than in 
            vortadv_sol2 but with a Hz(k)/H
               
            This is wrong and vortadv_sol3_mean should
            be used instead
            '''
               
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)


            var = toolsF_g.get_adv_sol2_mean(u,v, z_r,z_w,pm,pn) * tools.rho2psi(mask)
            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan      

 
        ################################################

        elif self.name in ['vortadv_sol3_mean']:

            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = toolsF_g.get_adv_sol3_mean(u,v, z_r,z_w,pm,pn) * tools.rho2psi(mask)
            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan      
           
        ###################################################################################

        elif self.name in ['rotwind_mean','rotbot_mean']:
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan

            if self.name=='rotwind_mean':
                
                [ut,vt] = self.get_winds(simul,coord=coord)

                ut=ut/simul.rho0/tools.rho2u(Hz); vt=vt/simul.rho0/tools.rho2v(Hz)
                
                var = tools.rot(ut,vt,pm,pn) * tools.rho2psi(mask)     


            elif self.name=='rotbot_mean':
                
                u = self.load('u',ncfile,simul,coord=coord,depths=[1])
                v = self.load('v',ncfile,simul,coord=coord,depths=[1])

                (ut,vt) = toolsF_g.get_bot(u,v,Hz,simul.rdrg,simul.Zob)
                
                ut=ut/tools.rho2u(Hz); vt=vt/tools.rho2v(Hz)

                var = -1*tools.rot(ut,vt,pm,pn) * tools.rho2psi(mask)     

        ###################################################################################
        ###################################################################################
        ###################################################################################

        
        elif self.name in ['J2_sol1','J2_sol2']:

            N=10;

            T = self.load('temp',ncfile,simul,coord=coord,depths=simul.coordmax[4][-1*N:])
            S = self.load('salt',ncfile,simul,coord=coord,depths=simul.coordmax[4][-1*N:])
            hbls = self.load('hbls',ncfile,simul,coord=coord,depths=[0])

            [u,v] = self.get_winds(simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            # if hbls = 0 (or very small), use depth of first cell instead.
            Hz = z_w[:,:,-1] - z_w[:,:,-2]
            hbls[hbls<Hz] = Hz[hbls<Hz]

            z_r = np.asfortranarray(z_r[:,:,simul.coordmax[4][-1*N:]-1])
            z_w = np.asfortranarray(z_w[:,:,simul.coordmax[4][-1*N-1:]])

            if self.name=='J2_sol1':
                var = tools_g.get_j2_sol1(T,S,u,v,z_r,z_w,simul.rho0,pm,pn,hbls)
            elif self.name=='J2_sol2':
                var = tools_g.get_j2_sol2(T,S,u,v,z_r,z_w,simul.rho0,pm,pn,hbls)

            var = var*tools.rho2psi(mask)




        ################################################


        elif self.name in ['Jbot_sol1','Jbot_sol2']:

            N=10; depths = np.arange(1,N+1,1)
            
            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)

            u = self.load('u',ncfile,simul,coord=coord,depths=[1])
            v = self.load('v',ncfile,simul,coord=coord,depths=[1])

            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            z_r = np.asfortranarray(z_r[:,:,:N])
            z_w = np.asfortranarray(z_w[:,:,:N+1])

            try:
                hbbls = self.load('hbbls',ncfile,simul,coord=coord,depths=[0])
                print('using actual hbbls')
            except:
                AKt = self.load('AKt',ncfile,simul,coord=coord,depths=list(range(50)))
                hbbls = toolsF_g.get_hbbls_from_akt(AKt,z_w)
                print('evaluating hbbls using AKt')

            if self.name=='Jbot_sol1':
                var = tools_g.get_jbot_sol1(T,S,u,v,z_r,z_w,simul.rho0,pm,pn,hbbls,simul.rdrg)
            elif self.name=='Jbot_sol2':
                var = tools_g.get_jbot_sol2(T,S,u,v,z_r,z_w,simul.rho0,pm,pn,hbbls,simul.rdrg)

            #var = self.others[0]

            var = var*tools.rho2psi(mask)

        ################################################


        elif self.name in ['Jbot_sol1_nohbbls']:

            N=10; depths = np.arange(1,N+1,1)
            
            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)

            u = self.load('u',ncfile,simul,coord=coord,depths=[1])
            v = self.load('v',ncfile,simul,coord=coord,depths=[1])

            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            z_r = np.asfortranarray(z_r[:,:,:N])
            z_w = np.asfortranarray(z_w[:,:,:N+1])


            hbbls = z_w[:,:,1]-z_w[:,:,0]

            self.others = tools.get_jbot_sol1(T,S,u,v,z_r,z_w,simul.rho0,pm,pn,hbbls,simul.rdrg)

            var = self.others[0]

            var = var*tools.rho2psi(mask)



        ################################################


        elif self.name in ['J1_sol1','J1_sol2']:


            N= 4
            u = self.load('u',ncfile,simul,coord=coord,depths=simul.coordmax[4][-1*N:])
            v = self.load('v',ncfile,simul,coord=coord,depths=simul.coordmax[4][-1*N:])


            hbls = self.load('hbls',ncfile,simul,coord=coord,depths=[0])
            [stflx, ssflx] = self.get_buoy_flux(simul,coord=coord,alphabeta=True)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            Hz = z_w[:,:,-1] - z_w[:,:,-2]

            # if hbls = 0 (or very small), use depth of first cell instead.
            hbls[hbls<Hz] = Hz[hbls<Hz]

            z_r = np.asfortranarray(z_r[:,:,simul.coordmax[4][-1*N:]-1])
            z_w = np.asfortranarray(z_w[:,:,simul.coordmax[4][-1*N-1:]])

            if self.name=='J1_sol1':
                var = tools.get_j1_sol1(stflx, ssflx, u,v,z_r,z_w,simul.rho0,pm,pn,hbls,f)
            elif self.name=='J1_sol2':
                var = tools.get_j1_sol2(stflx, ssflx, u,v,z_r,z_w,simul.rho0,pm,pn,hbls,f)

            var = var*tools.rho2psi(mask)

        ################################################
        
        
        elif self.name in ['Tadv','TXadv','TYadv','TVadv','THdiff','TVmix','TForc','Sadv','SHmix','SVmix','SForc','Trate','Srate']:

            [TXadv,TYadv,TVadv,THdiff,TVmix,TForc] = self.get_tracer_evolution(simul,coord=coord)
            
            if self.name=='Tadv':
              var = TXadv[:,:,:,0] + TYadv[:,:,:,0] + TVadv[:,:,:,0]  
            elif self.name=='TXadv':
              var = TXadv[:,:,:,0]
            elif self.name=='TYadv':
              var = TYadv[:,:,:,0]
            elif self.name=='TVadv':
              var = TVadv[:,:,:,0]                 
            elif self.name=='THdiff':
              var = THdiff[:,:,:,0]
            elif self.name=='TVmix':
              var = TVmix[:,:,:,0]+TForc[:,:,:,0]
            elif self.name=='TForc':
              var = TForc[:,:,:,0]
            elif self.name=='Trate':
              var = TXadv[:,:,:,0] + TYadv[:,:,:,0] + TVadv[:,:,:,0]  + TForc[:,:,:,0]  +  TVmix[:,:,:,0]        
            elif self.name=='Sadv':
              var = TXadv[:,:,:,1] + TYadv[:,:,:,1] + TVadv[:,:,:,1]
            elif self.name=='SHmix':
              var = THdiff[:,:,:,1]
            elif self.name=='SVmix':
              var = TVmix[:,:,:,1]+TForc[:,:,:,1]
            elif self.name=='SForc':
              var = TForc[:,:,:,1]
            elif self.name=='Srate':
              var = SXadv[:,:,:,0] + SYadv[:,:,:,0] + SVadv[:,:,:,0]  + SForc[:,:,:,0]  +  SVmix[:,:,:,0] 

              
        ################################################
        
        
        
        elif self.name in ['Tadv_pert','THdiff_pert','TVmix_pert','TForc_pert','Sadv_pert','SHmix_pert','SVmix_pert','SForc_pert','Trate_pert','Srate_pert']:

            [TXadv,TYadv,TVadv,THdiff,TVmix,TForc] = self.get_tracer_evolution(simul,coord=coord,pert=1)
            
            if self.name=='Tadv_pert':
              var = TXadv[:,:,:,0] + TYadv[:,:,:,0] + TVadv[:,:,:,0]  
            elif self.name=='THdiff_pert':
              var = THdiff[:,:,:,0]
            elif self.name=='TVmix_pert':
              var = TVmix[:,:,:,0]+TForc[:,:,:,0]
            elif self.name=='TForc_pert':
              var = TForc[:,:,:,0]
            elif self.name=='Trate_pert':
              var = TXadv[:,:,:,0] + TYadv[:,:,:,0] + TVadv[:,:,:,0]  + TForc[:,:,:,0]  +  TVmix[:,:,:,0]        
            elif self.name=='Sadv_pert':
              var = TXadv[:,:,:,1] + TYadv[:,:,:,1] + TVadv[:,:,:,1]
            elif self.name=='SHmix_pert':
              var = THdiff[:,:,:,1]
            elif self.name=='SVmix_pert':
              var = TVmix[:,:,:,1]+TForc[:,:,:,1]
            elif self.name=='SForc_pert':
              var = TForc[:,:,:,1]
            elif self.name=='Srate_pert':
              var = SXadv[:,:,:,0] + SYadv[:,:,:,0] + SVadv[:,:,:,0]  + SForc[:,:,:,0]  +  SVmix[:,:,:,0] 
              

        ################################################
        
        
        elif self.name in ['Madv','MXadv','MYadv','MVadv','MHmix','MHdiss','MVmix','MCor','MPrsgrd']:

            [MXadv,MYadv,MVadv,MHdiss,MHmix,MVmix,MCor,MPrsgrd] = self.get_uv_evolution(simul,coord=coord)
            
            if self.name=='Madv':
              var = MXadv[:,:,:,:] + MYadv[:,:,:,:] + MVadv[:,:,:,:]  
            elif self.name=='MXadv':
              var = MXadv[:,:,:,:]
            elif self.name=='MYadv':
              var = MYadv[:,:,:,:]
            elif self.name=='MVadv':
              var = MVadv[:,:,:,:]                 
            elif self.name=='MHmix':
              var = MHmix[:,:,:,:]
            elif self.name=='MHdiss':
              var = MHdiss[:,:,:,:]             
            elif self.name=='MVmix':
              var = MVmix[:,:,:,:]
            elif self.name=='MCor':
              var = MCor[:,:,:,:]
            elif self.name=='MPrsgrd':
              var = MPrsgrd[:,:,:,:]

              
        ################################################
        
        return var



            
        
        '''

        ################################################

        elif self.name in ['pv_sol1']:

            T = self.load('temp',ncfile,simul,coord=coord)
            S = self.load('salt',ncfile,simul,coord=coord)
            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            var = tools.get_pv_sol1(T,S,u,v,z_r,z_w,simul.rho0,pm,pn,f)



        ################################################

        elif self.name in ['pv_sol2']:

            T = self.load('temp',ncfile,simul,coord=coord)
            S = self.load('salt',ncfile,simul,coord=coord)
            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            var = tools.get_pv_sol2(T,S,u,v,z_r,z_w,simul.rho0,pm,pn,f)


        ################################################



        ###################################################################################
        #TERMS OF VORTICITY BALANCE
        ###################################################################################

        elif self.name in ['rotwind','rotbot']:


            if self.name=='rotwind':
                [ut,vt] = self.get_winds(simul,coord=coord)
                ut=ut/simul.rho0; vt=vt/simul.rho0

            elif self.name=='rotbot':
                u = self.load('u',ncfile,simul,coord=coord,depths=[1])
                v = self.load('v',ncfile,simul,coord=coord,depths=[1])
                [z_r,z_w] = tools.get_depths(simul,coord=coord)
                Hz = z_w[:,:,1] - z_w[:,:,0]

                [ut,vt] = tools.get_bottom_drag(u,v,Hz,simul.rdrg,simul.Zob)


            var = tools.rot(ut,vt,pm,pn) * tools.rho2psi(mask)


        ################################################

        elif self.name in ['vortbar']:


            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = tools.rot(ubar,vbar,pm,pn) * tools.rho2psi(mask)
            del ubar,vbar

        '''

        ################################################






###################################################################################
#Use old routines (temporary)
###################################################################################


    def oldvar(self,varname,simul,debug=False,**kwargs):

        ######################################
        #transformation to old system (temporary)
        ######################################

        ncfile = simul.ncfile
        #simul.ncname.his+'{0:04}'.format(simul.filetime)+ '.nc'
        
        if 'u' in  kwargs: 
            u = kwargs['u']; v = kwargs['v']
        else:
            u = None; v= None
            
            
            
        '''
        self.coord = simul.coord
         # if coordinates are given directly as an argument, use it instead of simul.coord[0:4]
        if 'coord' in  kwargs: self.coord[0:4] = kwargs['coord']
        '''
        
        [ny1,ny2,nx1,nx2,depths] = self.coord    
        if debug: print('[ny1,ny2,nx1,nx2,depths]',[ny1,ny2,nx1,nx2,depths])
        
        [ny0,nx0] = simul.coord[0],simul.coord[2]
        


        varnames = [varname]
        infiletime = simul.infiletime

        '''  
        Important:
        
        Error in the computation of grd variables (wrong subsetting)  corrected on June 02 2013
        Affected are all computation using computations in the form of:
        
        variables(varname, coord=coord, method='old)
        
        where old method is used with a coord smaller than the simul.coord, and where topo is used for the computation)
        '''
        '''
        [topo,pm,pn,f] = [simul.topo.T,
                            simul.pm.T,
                            simul.pn.T,
                            simul.f.T]        
        print topo.shape , topo[3,3]             
             
        
        [topo,pm,pn,f,lat,lon] = oldsim.variables_grd(simul.ncname.grd,ny1,ny2,nx1,nx2)
        print topo.shape , topo[3,3]      
        '''
        
        [topo,pm,pn,f] = [simul.topo.T[ny1-ny0:ny2-ny0,nx1-nx0:nx2-nx0],
                            simul.pm.T[ny1-ny0:ny2-ny0,nx1-nx0:nx2-nx0],
                            simul.pn.T[ny1-ny0:ny2-ny0,nx1-nx0:nx2-nx0],
                            simul.f.T[ny1-ny0:ny2-ny0,nx1-nx0:nx2-nx0]]     
                            
        if u!= None: u,v = u.T[ny1-ny0:ny2-ny0,nx1-nx0:nx2-nx0], v.T[ny1-ny0:ny2-ny0,nx1-nx0:nx2-nx0]
        if u!= None: print('in oldvar shapes are:', u.shape, v.shape)


        
        ncname = simul.ncname
        
        if 'depths' in  kwargs: depths = kwargs['depths']

        oldvar = oldsim.variables_fast(ncfile,varnames[0],ny1,ny2,nx1,nx2,depths,infiletime,\
                    topo=topo,pm=pm,pn=pn,f=f,u=u,v=v,cubic=1,ncname=ncname)

       
        self.data = np.squeeze(np.asfortranarray(oldvar.T))

        ######################################

        self.name = varname
        [self.imin, self.jmin, self.kmin] = [0,0,1]

        if varname in ['pv']:
            [self.imin, self.jmin, self.kmin] = [1,1,2]
        elif varname in ['vrt']:
            [self.imin, self.jmin, self.kmin] = [1,1,1]

#######################################################
#Get surface wind stress from forcing files interpolated to current time-step
#######################################################



    @staticmethod
    def get_winds(simul,**kwargs):
        ''' Load wind stress from forcing files '''

        if 'coord' in  kwargs: 
            [ny1,ny2,nx1,nx2]= kwargs['coord']
        else:
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]


        # Load NETCDF files
        ncfile = Dataset(simul.ncfile, 'r')

        try:
            print('check if sustr in ', simul.ncfile)
            
            #############################
            # test if periodic or tile
            #############################
            iper=0; jper=0

            try:
                if ncfile.dimensions['xi_u'].size==ncfile.dimensions['xi_rho'].size: iper=1
            except:
                pass

            try:
                if ncfile.dimensions['eta_v'].size==ncfile.dimensions['eta_rho'].size: jper=1
            except:
                pass
                
            print('iper,jper',iper,jper)
            #uwind = self.load('sustr',ncfile,simul,coord=[ny1,ny2,nx1,nx2])
            #vwind = self.load('svstr',ncfile,simul,coord=[ny1,ny2,nx1,nx2])
            uwind = simul.Forder( np.array(ncfile.variables['sustr'][simul.infiletime, ny1:ny2, nx1+iper:nx2-1+iper]) )
            vwind = simul.Forder( np.array(ncfile.variables['svstr'][simul.infiletime, ny1+jper:ny2-1+jper, nx1:nx2]) )
            print('loading sustr,svstr from his/avg file')
        
        except:
            print('computing sustr,svstr from frc file')
            ncfilewind = Dataset(simul.ncname.wind, 'r')

            try:
                oceantime = int(np.array(ncfile.variables['ocean_time'][simul.infiletime]))%(360*24*3600)
            except:
                oceantime = int(np.array(ncfile.variables['scrum_time'][simul.infiletime]))%(360*24*3600)

            oceanhour=oceantime/(3600.)
            oceanday=oceantime/(24*3600.)


            #if type=='d':

            print('Wind timesteps', ncfilewind.variables['sustr'].shape[0])
            if 360>=ncfilewind.variables['sustr'].shape[0]>12: #daily winds

                datewind1=int(np.floor(oceanday-0.5))%360
                datewind2=int(np.ceil(oceanday-0.5))%360

                if datewind1==datewind2:
                    coef1=0.5
                    coef2=0.5
                else:
                    coef1=abs(oceanday-0.5 - np.ceil(oceanday-0.5))
                    coef2=abs(oceanday-0.5 - np.floor(oceanday-0.5))


            elif ncfilewind.variables['sustr'].shape[0]>360: #hourly winds

                datewind1=int(np.floor(oceanhour-0.5))%(360*24)
                datewind2=int(np.ceil(oceanhour-0.5))%(360*24)

                if datewind1==datewind2:
                    coef1=0.5
                    coef2=0.5
                else:
                    coef1=abs(oceanhour-0.5 - np.ceil(oceanhour-0.5))
                    coef2=abs(oceanhour-0.5 - np.floor(oceanhour-0.5))

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


            uwind1=simul.Forder( np.array(ncfilewind.variables['sustr'][datewind1,ny1:ny2,nx1:nx2-1]) )
            vwind1=simul.Forder( np.array(ncfilewind.variables['svstr'][datewind1,ny1:ny2-1,nx1:nx2]) )

            uwind2=simul.Forder( np.array(ncfilewind.variables['sustr'][datewind2,ny1:ny2,nx1:nx2-1]) )
            vwind2=simul.Forder( np.array(ncfilewind.variables['svstr'][datewind2,ny1:ny2-1,nx1:nx2]) )

            uwind=coef1*uwind1+coef2*uwind2
            vwind=coef1*vwind1+coef2*vwind2


            ncfile.close()
            ncfilewind.close()


        return [uwind,vwind]





#######################################################
#Get surface wind stress from forcing files interpolated to current time-step
#######################################################


    #@staticmethod
    def get_buoy_flux(self,simul,alphabeta=False,solar=False,**kwargs):
        '''
        if solar is True the solar part of radiation will be outupted
        
        if alphabeta is True the flux will be outputed in the form alpha * stflx and beta*ssfls
        which is useful to compute buoyancy flux
        
        '''
        if 'coord' in  kwargs: 
            [ny1,ny2,nx1,nx2]= kwargs['coord']
        else:
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]


        # Load NETCDF files
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
        ncfilewind = Dataset(simul.ncname.frc, 'r', format='NETCDF3_CLASSIC')

        oceantime = int(np.array(ncfile.variables['ocean_time'][simul.infiletime]))%(360*24*3600)
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

        shflx1=simul.Forder( np.array(ncfilewind.variables['shflux'][datewind1,ny1:ny2,nx1:nx2]) )
        shflx2=simul.Forder( np.array(ncfilewind.variables['shflux'][datewind2,ny1:ny2,nx1:nx2]) )
        swflx1=simul.Forder( np.array(ncfilewind.variables['swflux'][datewind1,ny1:ny2,nx1:nx2]) ) 
        swflx2=simul.Forder( np.array(ncfilewind.variables['swflux'][datewind2,ny1:ny2,nx1:nx2]) )
        dqdt1=simul.Forder( np.array(ncfilewind.variables['dQdSST'][datewind1,ny1:ny2,nx1:nx2])  )
        dqdt2=simul.Forder( np.array(ncfilewind.variables['dQdSST'][datewind2,ny1:ny2,nx1:nx2]) )
        sst1=simul.Forder( np.array(ncfilewind.variables['SST'][datewind1,ny1:ny2,nx1:nx2])  )
        sst2=simul.Forder( np.array(ncfilewind.variables['SST'][datewind2,ny1:ny2,nx1:nx2]) )
        sss1=simul.Forder( np.array(ncfilewind.variables['SSS'][datewind1,ny1:ny2,nx1:nx2]) ) 
        sss2=simul.Forder( np.array(ncfilewind.variables['SSS'][datewind2,ny1:ny2,nx1:nx2]) )
        swrad1=simul.Forder( np.array(ncfilewind.variables['swrad'][datewind1,ny1:ny2,nx1:nx2])  )
        swrad2=simul.Forder( np.array(ncfilewind.variables['swrad'][datewind2,ny1:ny2,nx1:nx2]) )

        #######################################################

        shflx = coef1*shflx1 + coef2*shflx2
        swflx = coef1*swflx1 + coef2*swflx2
        dqdt = coef1*dqdt1 + coef2*dqdt2
        sst = coef1*sst1 + coef2*sst2
        sss = coef1*sss1 + coef2*sss2
        swrad = coef1*swrad1 + coef2*swrad2

        #######################################################

        
        
        temp = self.load('temp',ncfile,simul,depths=[0],coord=[ny1,ny2,nx1,nx2])
        salt = self.load('salt',ncfile,simul,depths=[0],coord=[ny1,ny2,nx1,nx2])

        
        [alpha,beta] = self.alphabeta(temp,salt,simul.rho0)

        '''
        #          Bo(i,j)=g*( alpha(i,j)*(stflx(i,j,itemp)-srflx(i,j))
        # ifdef SALINITY
        #     &                              -beta(i,j)*stflx(i,j,isalt)
        # endif
        #     &                                                        )
        '''

        [z_r,z_w] = tools.get_depths(simul,coord=[ny1,ny2,nx1,nx2])
        Hz = np.asfortranarray(z_w[:,:,-1] - z_w[:,:,-2])
        del z_r, z_w

        #swr_frac = 0.58 * np.exp(-Hz/0.35) + (1-0.58) * np.exp(-Hz/23.)

        #######################################################

        rho0=simul.rho0
        Cp=3985.
        
        #######################################################

        
        
        if 'QCORRECTION' in simul.cpp: clim = 1.
        else: clim = 0.

        stflx = shflx/(rho0*Cp) + clim * dqdt/(rho0*Cp)*(temp - sst)

        #stflx2 = stflx - swr_frac * swrad/(rho0*Cp)

        dsdt = 1./(90.*86400.)
        ssflx = swflx*0.01/86400*salt - clim * dsdt*Hz*(salt - sss)
        
        
        #######################################################

        if 'DIURNAL_SRFLUX' in simul.cpp or 'ROBUST_DIURNAL_SRFLUX' in simul.cpp:


            sec2day=1./86400.
            tdays = simul.oceantime*sec2day

            if 'ROBUST_DIURNAL_SRFLUX' in simul.cpp:
                print('ROBUST_DIURNAL_SRFLUX has not been validated yet')

                phase = copy(stflx)*0.

                deg2rad = np.pi/180.
                cff = 1/360.
                dec=-0.406*np.cos(deg2rad*(tdays-np.floor(cff*tdays)*360.))
                cos_d = np.cos(dec)
                sin_d = np.sin(dec)
                tan_d = np.tan(dec)
                for i in range(stflx.shape[0]):
                    for j in range(stflx.shape[1]):
                        cos_h = np.cos(2*np.pi*((tdays-np.floor(tdays)) + cff *simul.x[i,j]))
                        phi = simul.y[i,j]*deg2rad
                        h0 = np.acos(-np.tan(phi)*tan_d)
                        cff1 = cos_d*np.cos(phi)
                        cff2 = sin_d*np.sin(phi)
                        phase[i,j] = np.pi*(cos_h*cff1+cff2)/(np.sin(h0)*cff1+h0*cff2)
                        phase[i,j] = np.max([0.,phase[i,j]])

            else:

                cff=2.*simul.dt_model*sec2day
                phase=4.*(tdays-np.floor(tdays))-2.
                cff1=np.max([-1., np.min([1., np.min(phase-cff)])])
                cff2=np.max([-1., np.min([1., np.min(phase+cff)])])
                phase=(cff2-cff1)/cff + (np.sin(np.pi*cff2)-np.sin(np.pi*cff1))/(np.pi*cff)

            print('phase is ', phase)


            srflx =  phase * swrad/(rho0*Cp)
            stflx = stflx + (phase - 1.) * swrad/(rho0*Cp)


        else:

            srflx =  swrad/(rho0*Cp)

        
        #######################################################


        ncfile.close()
        ncfilewind.close()

        if alphabeta: 
            print('outputing alpha*stflx and beta*ssflx ')
            return [np.asfortranarray(alpha*stflx,dtype=simul.floattype), np.asfortranarray(beta*ssflx,dtype=simul.floattype)]
        else:
            if solar:
                return [stflx, ssflx, srflx]
            else:
                return [np.asfortranarray(stflx,dtype=simul.floattype), np.asfortranarray(ssflx,dtype=simul.floattype)]




    #######################################################
    # Get wind stress using ONLINE interpolation
    # ------------------------------------------
    # Main function is get_online_bulk_winds
    #
    # it uses subfunctions:
    #       get_bulk_coords
    #       get_bulk_winds
    #       get_bulk_vars
    #       interpolate_bulk_winds
    #       interpolate_bulk_vars
    #
    #######################################################
    


    def get_bulk_coords(self, simul, time=None, **kwargs):
    
        '''
        Get Bulk files location and space-time coordinates corresponding to simul
        
        Designed for hourly CFSR
        
        '''
        
        try:
            if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                folder='/ccc/scratch/cont003/gen7638/gulaj/CFSR1h/'
            else:
                folder= '/home/datawork-lops-osi/jgula/CFSR1h/'
        except:
            folder = '/net/omega/local/tmp/1/gula/CFSR/'
        
        if 'coord' in  kwargs:
            coord = kwargs['coord'][0:4]
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1s,ny2s,nx1s,nx2s] = simul.coord[0:4]
            x = np.asfortranarray(simul.x[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])
            y = np.asfortranarray(simul.y[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])
        else:
            coord = simul.coord[0:4]
            x = simul.x; y = simul.y
            
        ####################
        # time of simulation
        if time is None: time = simul.oceantime
        time_in_days = np.float(time)/3600./24.
        
        #############################
        # closest date for bulk files
        date1 = simul.ncname.realyear_origin + timedelta(days=time_in_days)
        year1 = date1.year; month1 = date1.month
        file_U1 = folder + '/U-component_of_wind_Y' + format(year1) + 'M' + format(month1) + '.nc'
        #file_V1 = folder + '/V-component_of_wind_Y' + format(year1) + 'M' + format(month1) + '.nc'
        
        # one hour before
        date0 = simul.ncname.realyear_origin + timedelta(days=time_in_days-1./24)
        year0 = date0.year; month0 = date0.month

        ###
        # one hour after
        date2 = simul.ncname.realyear_origin + timedelta(days=time_in_days+1./24.)
        year2 = date2.year; month2 = date2.month

        ##########
        # get indices
        ncfile_U = Dataset(file_U1,'r')
        
        # find time indices
        time_bulk = ncfile_U.variables['time'][:]; nt = len(time_bulk)
        
        if time_bulk.min()<time_in_days<time_bulk.max():
            itime2 = np.nanargmax(time_bulk >= time_in_days)
            itime1 = itime2-1
        else:
            print('between 2 files [to be checked]')
            itime1 = nt-1
            itime2 = 0
            print('file 1 : ',file_U1)
            print('year1,month1,itime1,year2,month2,itime2 :\
                    ',year1,month1,itime1,year2,month2,itime2)
            
        ##########
        
        # find domain
        #ncfile_U.variables['lon'].valid_range = [-180.,180.] # fix because format issue
        lon_bulk = ncfile_U.variables['lon'][:]
        lat_bulk = ncfile_U.variables['lat'][:]
        
        # min/max domain coordinates
        lonmin = x.min()
        lonmax = x.max()
        latmin = y.min()
        latmax = y.max()
        
        # CFSR indices
        imin = np.max([0,np.nanargmax(lon_bulk >= lonmin)-2])
        imax = np.min([len(lon_bulk),np.nanargmin(lon_bulk <= lonmax)+2])
        jmin = np.max([0,np.nanargmax(lat_bulk >= latmin)-2])
        jmax = np.min([len(lat_bulk),np.nanargmin(lat_bulk <= latmax)+2])
        print('imin,imax,jmin,jmax',imin,imax,jmin,jmax)
        ncfile_U.close()
        
        lon_bulk = lon_bulk[imin:imax]
        lat_bulk = lat_bulk[jmin:jmax]
        
        return year1,  month1, itime1, year2, month2, itime2, imin, imax, jmin, jmax, lon_bulk, lat_bulk, folder

    ########################################

    def get_bulk_winds(self, simul, time=None, **kwargs):
            
        '''
        
        Load winds from bulk files
        
        Designed for hourly CFSR
        
        '''
        
        if 'coord' in  kwargs:
            coord = kwargs['coord'][0:4]
        else:
            coord = simul.coord[0:4]


        if time is None: time = simul.oceantime
        
        [year1,month1,itime1,year2,month2,itime2,imin,imax,jmin,jmax,lon_bulk,lat_bulk,folder] =\
                                                    self.get_bulk_coords(simul, coord = coord, time = time)
        
        print('get_bulk_coords OK')
        print('year1,month1,itime1,year2,month2,itime2,imin,imax,jmin,jmax',\
               year1,month1,itime1,year2,month2,itime2,imin,imax,jmin,jmax)
        
        file_U1 = folder + '/U-component_of_wind_Y' + format(year1) + 'M' + format(month1) + '.nc'
        file_V1 = folder + '/V-component_of_wind_Y' + format(year1) + 'M' + format(month1) + '.nc'
        file_U2 = folder + '/U-component_of_wind_Y' + format(year2) + 'M' + format(month2) + '.nc'
        file_V2 = folder + '/V-component_of_wind_Y' + format(year2) + 'M' + format(month2) + '.nc'
        
        # for U at time 1
        ncfile_U = Dataset(file_U1,'r')
        time_bulk = ncfile_U.variables['time'][:]
        uwnd1 = ncfile_U.variables['U-component_of_wind'][itime1,jmin:jmax,imin:imax]
        ncfile_U.close()
        print('get uwnd1 OK')

        # for U  at time 2
        ncfile_U = Dataset(file_U2,'r')
        uwnd2 = ncfile_U.variables['U-component_of_wind'][itime2,jmin:jmax,imin:imax]
        ncfile_U.close()
        print('get uwnd2 OK')
        
        # for V  at time 1
        ncfile_U = Dataset(file_V1,'r')
        vwnd1 = ncfile_U.variables['V-component_of_wind'][itime1,jmin:jmax,imin:imax]
        ncfile_U.close()
        print('get vwnd1 OK')
        
        # for V  at time 2
        ncfile_U = Dataset(file_V2,'r')
        vwnd2 = ncfile_U.variables['V-component_of_wind'][itime2,jmin:jmax,imin:imax]
        ncfile_U.close()
        print('get vwnd2 OK')
        
        #####################
        # time interpolation
        
        time_in_days = np.float(time)/3600./24.
        
        cff1 = time_bulk[itime2] - time_in_days
        cff2 = time_in_days - time_bulk[itime1]
        cff = 1./(cff1+cff2)
        cff1 *= cff; cff2 *= cff
        print('time coef OK')
        #####################
        
        return lon_bulk, lat_bulk, cff1,uwnd1,vwnd1, cff2,uwnd2,vwnd2


    ########################################

    def get_bulk_vars(self, simul, time=None, **kwargs):
    
        '''
        
        Load tair,rhum from bulk files
        
        Designed for hourly CFSR
        
        '''
        
        if 'coord' in  kwargs:
            coord = kwargs['coord'][0:4]
        else:
            coord = simul.coord[0:4]
        
        if time is None: time = simul.oceantime
        
        [year1,month1,itime1,year2,month2,itime2,imin,imax,jmin,jmax,lon_bulk,lat_bulk,folder] =\
                                                    self.get_bulk_coords(simul, coord = coord, time = time)
                                                    
        file_t1 = folder + '/Temperature_height_above_ground_Y' + format(year1) + 'M' + format(month1) + '.nc'
        file_r1 = folder + '/Specific_humidity_Y' + format(year1) + 'M' + format(month1) + '.nc'
        
        file_t2 = folder + '/Temperature_height_above_ground_Y' + format(year2) + 'M' + format(month2) + '.nc'
        file_r2 = folder + '/Specific_humidity_Y' + format(year2) + 'M' + format(month2) + '.nc'
        
        # for tair at time 1
        ncfile_t = Dataset(file_t1,'r')
        time_bulk = ncfile_t.variables['time'][:]
        
        #ncfile_t.variables['lon'].valid_range = [-180.,180.] # fix because format issue
        
        tair1 = ncfile_t.variables['Temperature_height_above_ground'][itime1,jmin:jmax,imin:imax]
        ncfile_t.close()
        
        # for tair  at time 2
        ncfile_t = Dataset(file_t2,'r')
        tair2 = ncfile_t.variables['Temperature_height_above_ground'][itime2,jmin:jmax,imin:imax]
        ncfile_t.close()
        
        # for rhum  at time 1
        ncfile_r = Dataset(file_r1,'r')
        rhum1 = ncfile_r.variables['Specific_humidity'][itime1,jmin:jmax,imin:imax]
        ncfile_r.close()
        
        # for rhum  at time 2
        ncfile_r = Dataset(file_r2,'r')
        rhum2 = ncfile_r.variables['Specific_humidity'][itime2,jmin:jmax,imin:imax]
        ncfile_r.close()
        
        #####################
        # time interpolation
        
        time_in_days = np.float(time)/3600./24.
        
        cff1 = time_bulk[itime2] - time_in_days
        cff2 = time_in_days - time_bulk[itime1]
        cff = 1./(cff1+cff2)
        cff1 *= cff; cff2 *= cff
        
        tair = cff1 * tair1 + cff2 * tair2
        rhum = cff1 * rhum1 + cff2 * rhum2
        
        #####################
        
        return lon_bulk, lat_bulk, tair, rhum
        

    ########################################

    def interpolate_bulk_winds(self, simul, time=None, **kwargs):
    
        '''
        
        Interpolate winds on model grid
        
        Designed for hourly CFSR
        
        '''
        
        if 'coord' in  kwargs:
            coord = kwargs['coord'][0:4]
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1s,ny2s,nx1s,nx2s] = simul.coord[0:4]
            x = np.asfortranarray(simul.x[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])
            y = np.asfortranarray(simul.y[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])
            angle = np.asfortranarray(simul.angle[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])
        else:
            coord = simul.coord[0:4]
            x = simul.x; y = simul.y; angle = simul.angle
            
        if time is None: time = simul.oceantime
        
        [lon_cfsr,lat_cfsr, cff1,uwnd1,vwnd1, cff2,uwnd2,vwnd2] = \
                               self.get_bulk_winds(simul,coord=coord, time=time)
        print('get bulk winds OK')
        
        NX_cfsr = len(lon_cfsr); NY_cfsr = len(lat_cfsr)
        [nx,ny] = x.shape
        Istr, Iend, Jstr, Jend = 0, nx-1, 0, ny-1

        '''   call MYINTERP(1, NXref, 1, NYref,
         &                lonref, latref, varref,
         &                I_RANGE,J_RANGE,
         &                lonr(GLOBAL_2D_ARRAY),
         &                latr(GLOBAL_2D_ARRAY),
         &                uwndg(GLOBAL_2D_ARRAY,iblkrec))'''


        #print(1, NX_cfsr, 1, NY_cfsr,lon_cfsr, lat_cfsr,Istr, Iend, Jstr, Jend)
        #print(x.min(),x.max(),y.min(),y.max())

        
        u_blk1 = toolsF_gi.cinterp2d(1, NX_cfsr, 1, NY_cfsr,\
                                 lon_cfsr, lat_cfsr, uwnd1.T,\
                                 Istr, Iend, Jstr, Jend,\
                                 x,y)
        print('interp uwnd1 OK')
        
        v_blk1 = toolsF_gi.cinterp2d(1, NX_cfsr, 1, NY_cfsr,\
                                 lon_cfsr,lat_cfsr, vwnd1.T,\
                                 Istr, Iend, Jstr, Jend,\
                                 x,y)
        print('interp vwnd1 OK')
        
        u_blk2 = toolsF_gi.cinterp2d(1, NX_cfsr, 1, NY_cfsr,\
                                 lon_cfsr,lat_cfsr, uwnd2.T,\
                                 Istr, Iend, Jstr, Jend,\
                                 x,y)
        print('interp uwnd2 OK')
        
        v_blk2 = toolsF_gi.cinterp2d(1, NX_cfsr, 1, NY_cfsr,\
                                 lon_cfsr,lat_cfsr, vwnd2.T,\
                                 Istr, Iend, Jstr, Jend,\
                                 x,y)
        print('interp vwnd2 OK')
        
        # compute wspd:
        wspd = cff1 * np.sqrt(u_blk1**2 + v_blk1**2) + cff2 * np.sqrt(u_blk2**2 + v_blk2**2)
        print('wspd OK')
        
        # rotate to model directions:
        u_blk_rot1 = u_blk1*np.cos(angle) + v_blk1*np.sin(angle)
        v_blk_rot1 = v_blk1*np.cos(angle) - u_blk1*np.sin(angle)
        print('rot 1 OK')
        del u_blk1,v_blk1
     
        # rotate to model directions:
        u_blk_rot2 = u_blk2*np.cos(angle) + v_blk2*np.sin(angle)
        v_blk_rot2 = v_blk2*np.cos(angle) - u_blk2*np.sin(angle)
        print('rot 2 OK')
        del u_blk2,v_blk2

        # time interpolation
        u_blk_rot = cff1 * u_blk_rot1 + cff2 * u_blk_rot2
        v_blk_rot = cff1 * v_blk_rot1 + cff2 * v_blk_rot2
        print('time interp OK')
        
        u_blk_rot = tools.rho2u(u_blk_rot)
        v_blk_rot = tools.rho2v(v_blk_rot)
        print('to u,v grid OK')
        
        return u_blk_rot,v_blk_rot, wspd

    ####################


    def interpolate_bulk_vars(self, simul, time=None, **kwargs):
    
        '''
        
        Interpolate tair,rhum on model grid
        
        Designed for hourly CFSR
        
        '''
        
        if 'coord' in  kwargs:
            coord = kwargs['coord'][0:4]
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1s,ny2s,nx1s,nx2s] = simul.coord[0:4]
            x = np.asfortranarray(simul.x[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])
            y = np.asfortranarray(simul.y[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])
            angle = np.asfortranarray(simul.angle[nx1i-nx1s:nx2i-nx1s,ny1i-ny1s:ny2i-ny1s])
        else:
            coord = simul.coord[0:4]
            x = simul.x; y = simul.y
            
        if time is None: time = simul.oceantime
        
        [lon_cfsr,lat_cfsr, tair, rhum] = self.get_bulk_vars(simul, coord = coord, time = time)

        NX_cfsr = len(lon_cfsr); NY_cfsr = len(lat_cfsr)
        [nx,ny] = x.shape
        Istr, Iend, Jstr, Jend = 0, nx-1, 0, ny-1

        tair = toolsF_gi.cinterp2d(1, NX_cfsr, 1, NY_cfsr,\
                                 lon_cfsr,lat_cfsr, tair.T,
                                 Istr, Iend, Jstr, Jend,
                                 x,y)
                                 
        ###########
        # Temperature: Convert from Kelvin to Celsius
        tair -= 273.15
            
        ###############
        
        rhum = toolsF_gi.cinterp2d(1, NX_cfsr, 1, NY_cfsr,\
                                 lon_cfsr,lat_cfsr, rhum.T,
                                 Istr, Iend, Jstr, Jend,
                                 x,y)
        
        ##########
        # Convert specific to relative humidity
        Pref=1020.         # default air pressure [mbars]
        ew=6.1121*(1.0007+3.46e-6*Pref)* np.exp( (17.502*tair) / (240.97+tair))
        Qsat=0.62197*(ew/(Pref-0.378*ew))
        rhum /= Qsat
        print('interp other vars OK')
        
        return tair, rhum

    #################

    ###########

    def cfb_stress(self,sustr_nocfb,svstr_nocfb,u,v,wspd,rho0):
        
        '''
        
        Compute wind-stress correction
        
        '''
        
        #wind-stress correction using wind speed:  rho0*sustr + s_tau*Uo
        #   s_tau = cfb_slope * wspd + cfb_offset [N.m^-3.s]
        #  (recommendended and default if BULK_FLUX - needs wspd data)
        
        cfb_slope=-0.0029
        cfb_offset=0.008
        Wspd_min=3.        # [m/s]
        Wstr_min=0.045     # [N/m2]
        stau_ref=-0.0027   # [N.m^-3.s]
         
        ########
        
        cff1=tools.rho2u(wspd)
        cff = cfb_slope*cff1 + cfb_offset
        cff[cff1<=Wspd_min] = stau_ref

        sustr = sustr_nocfb + cff*u/rho0

        ########
                
        cff1=tools.rho2v(wspd)
        cff = cfb_slope*cff1 + cfb_offset
        cff[cff1<=Wspd_min] = stau_ref
            
        svstr = svstr_nocfb + cff*v/rho0
        
        ########
        
        return sustr,svstr
        

    #######################################################
    # Get bulk wind stress
    #######################################################
    
    def get_online_bulk_winds(self, simul, time=None, **kwargs):
    
        ''' Main function'''
        
        if 'coord' in  kwargs:
            coord = kwargs['coord'][0:4]
        else:
            coord = simul.coord[0:4]
            
        if time is None:
            if 'xios' in simul.ncname.model:
                time = simul.oceantime - 0.5 * simul.dt_model
            else:
                time = simul.oceantime + 0.5 * simul.dt_model
                
        [u_blk_rot,v_blk_rot, wspd] = \
        self.interpolate_bulk_winds(simul, coord = coord, time=time)
        [tair, rhum] = self.interpolate_bulk_vars(simul, coord = coord, time=time)
        print('get all vars OK')
        
        ncfile = Dataset(simul.ncfile, 'r')
        temp = self.load('temp',ncfile,simul,coord=coord,depths=[0])
        ncfile.close()
        
        print('get temp OK')
        
        [sustr_nocfb,svstr_nocfb] = toolsF_gi.bulk_stress(temp, tair, rhum, \
                                       u_blk_rot, v_blk_rot, wspd, simul.rho0)
                                       
        del temp, tair, rhum, u_blk_rot, v_blk_rot
        
        print('bulk stress OK')
        
        # Load u,v
        ncfile = Dataset(simul.ncfile, 'r')
        u = self.load('u',ncfile,simul,coord=coord,depths=[0])
        v = self.load('v',ncfile,simul,coord=coord,depths=[0])
        ncfile.close()
        print('get u,v OK')
        
        [sustr,svstr] = self.cfb_stress(sustr_nocfb,svstr_nocfb,u,v,wspd,simul.rho0)
        print('get cfb_stress OK')
        del sustr_nocfb,svstr_nocfb,u,v,wspd
        
        return sustr,svstr
        
        


#######################################################
#Compute thermal expansion and saline contraction coefficients
#######################################################

    @staticmethod
    def alphabeta(Tt,Ts,rho0):
        '''
! Compute thermal expansion and saline contraction coefficients
! as functions of potential temperature, salinity from a polynomial
! expression (Jackett & McDougall, 1992). The coefficients are
! evaluated at the surface.
!
!  alpha(Ts,Tt,0)=-d(rho1(Ts,Tt,0))/d(Tt) / rho0
!  beta(Ts,Tt,0) = d(rho1(Ts,Tt,0))/d(Ts) / rho0
!
!  Adapted from original "rati" and "beta" routines.
        '''
        
        '''
 In order to get linearized buoyancy:
    [alpha,beta] = var.alphabeta(T[:,:,iz],S[:,:,iz],simul.rho0)
    buoy_lin = simul.g * (alpha.T * (T.T - T[:,:,iz].T) - beta.T * (S.T - S[:,:,iz].T) ).T + buoy[:,:,iz]
        
        '''

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
#Compute  tendency terms fron the tracer equation
#######################################################



    def get_tracer_evolution(self,simul,**kwargs):
    

        #######################################################
        if 'pert' in  kwargs: 
            pert = kwargs['pert']
        else:
            pert=None
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            print('original coord', coord)
            
            # We need at leat 4 pts to compute advection:
            addx1,addy1,addx2,addy2 =0,0,None,None
            if (ny2i-ny1i)<6: 
                addy=6-(ny2i-ny1i); addy1= (addy/2); addy2=-(addy/2+addy%2)
                ny1i=ny1i-addy1; ny2i=ny2i-addy2
            if (nx2i-nx1i)<6: 
                addx=6-(nx2i-nx1i); addx1= (addx/2); addx2=-(addx/2+addx%2)
                nx1i=nx1i-addx1; nx2i=nx2i-addx2      
            coord = [ny1i,ny2i,nx1i,nx2i]        
            print('temporary coord', coord)

            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        else:
             # We need at leat 4 pts to compute advection:
             # But we can't enlarge the domain...
            addx1,addy1,addx2,addy2 =0,0,None,None
            
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            

           
        #if 'depths' in  kwargs: 
            #depths = kwargs['depths']
        #else: 
            #depths = simul.coord[4]
            
        depths = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

        #######################################################
 
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
            
        #######################################################
        [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)
        
        u = self.load('u',ncfile,simul,coord=coord,depths=depths)
        #if pert==1: u = u - np.mean(u); print 'removing mean velocity'
        #elif pert==0: u = u * np.nan + np.mean(u); print 'using only mean velocity'
        v = self.load('v',ncfile,simul,coord=coord,depths=depths)
        #if pert==1: v = v - np.mean(v)        
        #elif pert==0: v = v * np.nan + np.mean(v); print 'using only mean velocity'
        
        try:
            omega = self.load('omega',ncfile,simul,coord=coord,depths=depths_w)
        except:
            print('no omega in file, computing')
            omega = tools.nanbnd(toolsF.get_omega (u,v,z_r,z_w,pm,pn))
        
        AKt = self.load('AKt',ncfile,simul,coord=coord,depths=depths_w)
        hbls = self.load('hbls',ncfile,simul,coord=coord)

        t=np.zeros((u.shape[0]+1,u.shape[1],u.shape[2],2))
        t[:,:,:,0] = self.load('temp',ncfile,simul,coord=coord,depths=depths)
        t[:,:,:,1] = self.load('salt',ncfile,simul,coord=coord,depths=depths)
        
        Hz =  z_w[:,:,1:] - z_w[:,:,:-1]
        
        stflx=np.zeros((u.shape[0]+1,u.shape[1],2))
        [stflx[:,:,0], stflx[:,:,1],srflx] = self.get_buoy_flux(simul,coord=coord,solar=True)
        
        [alpha,beta] = self.alphabeta(t[:,:,-1,0],t[:,:,-1,1],simul.rho0)

        

        (swr_frac) = \
                toolsF.get_swr_frac (Hz)    
        del Hz
       
        (ghat) = \
                toolsF.get_ghat (alpha,beta, z_r,z_w\
                             ,stflx, srflx,  hbls,  swr_frac) 
        #ghat= z_r*0.              
     
    
        (TXadv,TYadv,TVadv,THdiff,TVmix,TForc) = \
                toolsF_g.get_tracer_evolution (u,v, z_r,z_w,pm,pn\
                             ,simul.dt_model,t,stflx, srflx, ghat, swr_frac, omega, AKt) 
            
        [TXadv,TYadv,TVadv,THdiff,TVmix,TForc] = \
             [tools.nanbnd(TXadv,2),
             tools.nanbnd(TYadv,2),     
             tools.nanbnd(TVadv,2),
             tools.nanbnd(THdiff,2),      
             tools.nanbnd(TVmix,2),
             tools.nanbnd(TForc,2)]        
             
             
             
             
        #print 'computing old '
        #(TXadv,TYadv,TVadv,THdiff,TVmix,TForc) = \
                #toolsF.get_tracer_evolution_old (u,v, z_r,z_w,pm,pn\
                             #,simul.dt_model,t,stflx, srflx, ghat, omega, AKt)
        if len(TXadv.shape)==4:
            return [TXadv[addx1:addx2,addy1:addy2,:,:],TYadv[addx1:addx2,addy1:addy2,:,:],TVadv[addx1:addx2,addy1:addy2,:,:],\
                THdiff[addx1:addx2,addy1:addy2,:,:],TVmix[addx1:addx2,addy1:addy2,:,:],TForc[addx1:addx2,addy1:addy2,:,:]]
        elif len(TXadv.shape)==3:
            return [TXadv [addx1:addx2,addy1:addy2,:],TYadv [addx1:addx2,addy1:addy2,:],TVadv [addx1:addx2,addy1:addy2,:],\
                THdiff [addx1:addx2,addy1:addy2,:],TVmix [addx1:addx2,addy1:addy2,:],TForc [addx1:addx2,addy1:addy2,:]]
        
        




#######################################################
#Compute tendency terms fron the tracer equation
#######################################################

    def get_tend(self,simul,**kwargs):

        #######################################################
        if 'pert' in  kwargs: 
            pert = kwargs['pert']
        else:
            pert=None
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            print('original coord', coord)
            
            # We need at leat 4 pts to compute advection:
            addx1,addy1,addx2,addy2 =0,0,None,None
            if (ny2i-ny1i)<6: 
                addy=6-(ny2i-ny1i); addy1= (addy/2); addy2=-(addy/2+addy%2)
                ny1i=ny1i-addy1; ny2i=ny2i-addy2
            if (nx2i-nx1i)<6: 
                addx=6-(nx2i-nx1i); addx1= (addx/2); addx2=-(addx/2+addx%2)
                nx1i=nx1i-addx1; nx2i=nx2i-addx2      
            coord = [ny1i,ny2i,nx1i,nx2i]        
            print('temporary coord', coord)

            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        else:
             # We need at leat 4 pts to compute advection:
             # But we can't enlarge the domain...
            addx1,addy1,addx2,addy2 =0,0,None,None
            
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            

           
        #if 'depths' in  kwargs: 
            #depths = kwargs['depths']
        #else: 
            #depths = simul.coord[4]
            
        depths = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

        #######################################################
 
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
            
        #######################################################
        [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)
        
        u = self.load('u',ncfile,simul,coord=coord,depths=depths)
        if pert==1: u = u - np.mean(u); print('removing mean velocity')
        elif pert==0: u = u * np.nan + np.mean(u); print('using only mean velocity')
        v = self.load('v',ncfile,simul,coord=coord,depths=depths)
        if pert==1: v = v - np.mean(v)        
        elif pert==0: v = v * np.nan + np.mean(v); print('using only mean velocity')
        try:
            omega = self.load('omega',ncfile,simul,coord=coord,depths=depths_w)
        except:
            print('no omega in file: ', simul.ncfile)
            omega = toolsF.get_omega (u,v,z_r,z_w,pm,pn)    
            
        AKt = self.load('AKt',ncfile,simul,coord=coord,depths=depths_w)
        hbls = self.load('hbls',ncfile,simul,coord=coord)

        t=np.zeros((u.shape[0]+1,u.shape[1],u.shape[2],2))
        t[:,:,:,0] = self.load('temp',ncfile,simul,coord=coord,depths=depths)
        t[:,:,:,1] = self.load('salt',ncfile,simul,coord=coord,depths=depths)
        
        Hz =  z_w[:,:,1:] - z_w[:,:,:-1]
        
        stflx=np.zeros((u.shape[0]+1,u.shape[1],2))
        [stflx[:,:,0], stflx[:,:,1],srflx] = self.get_buoy_flux(simul,coord=coord,solar=True)
        
        [alpha,beta] = self.alphabeta(t[:,:,-1,0],t[:,:,-1,1],simul.rho0)

        print('swr_frac')
        (swr_frac) = \
                toolsF.get_swr_frac (Hz)    
        del Hz
        
        print('ghat')
        (ghat) = \
                toolsF.get_ghat (alpha,beta, z_r,z_w\
                             ,stflx, srflx,  hbls, swr_frac) 
        #ghat= z_r*0.              
     
        print('get_tracer_evolution')   
        (TXadv,TYadv,TVadv,THdiff,TVmix,TForc) = \
                toolsF_g.get_tracer_evolution (u,v, z_r,z_w,pm,pn\
                             ,simul.dt_model,t,stflx, srflx, ghat, swr_frac, omega, AKt) 
        
        
        TVmix = tools.nanbnd(TVmix + TForc,2)
        THdiff = tools.nanbnd(THdiff,2)
        
        del TXadv,TYadv,TVadv,TForc
        
        
        #######################################################

        buoy = -simul.g * toolsF.rho_eos(t[:,:,:,0],t[:,:,:,1],z_r,z_w,simul.rho0)/simul.rho0

        #######################################################
        # Everything is similar to get_tracer_evolution till now
        # Now we are computing
        # tend = \nabla_h (alpha TXxxx - beta SXxxx) .  \nabla_h (b)
        #######################################################     
        
        if 'depths' in  kwargs: 
            depths = kwargs['depths']
        
        print('depths',  depths)
        
        if np.min(depths)<=0:
            
            #interpolate Tracer terms to depths
            
            ########################################################   
            # We first prepare w, bz cause we need full 3D buoy and u,v
            print('1')                           
            w = tools.vinterp(toolsF.get_wvlcty(u,v,z_r,z_w,pm,pn),depths,z_r,z_w)            
            bz = tools.vinterp((buoy[:,:,1:] - buoy[:,:,:-1])/ (z_r[:,:,1:] -z_r[:,:,:-1]),depths,z_w[:,:,1:-1],z_r)
            print('2')             
            bz2 = (buoy[:,:,1:] - buoy[:,:,:-1])/ (z_r[:,:,1:] -z_r[:,:,:-1])           
            bz2 = (bz2[:,:,1:] - bz2[:,:,:-1])/ (z_w[:,:,2:-1] -z_w[:,:,1:-2])
            bzz = tools.vinterp(bz2,depths,z_r[:,:,1:-1],z_w[:,:,1:-1])
            print('3')           
            AKtz = (AKt[:,:,1:] - AKt[:,:,:-1])/ (z_w[:,:,1:] -z_w[:,:,:-1])
            print('AKt.shape',AKt.shape)
            print('AKtz.shape',AKtz.shape)            
            print('z_r.shape',z_r.shape)    
            print('z_w.shape',z_w.shape)
            AKtz =  tools.vinterp(AKtz,depths,z_r,z_w)
            AKt =   tools.vinterp(AKt,depths,z_w)  
            print('4')
            
            ########################################################    

            if len(depths)>len(simul.coordmax[4]):
                
                THdiff0 =np.zeros((TVmix.shape[0],TVmix.shape[1],len(depths),2))*np.nan
                TVmix0 = np.zeros((TVmix.shape[0],TVmix.shape[1],len(depths),2))*np.nan 
                
                for it in [0,1]:
                    [THdiff0[:,:,:,it],TVmix0[:,:,:,it]]=\
                            [tools.vinterp(THdiff[:,:,:,it],depths,z_r,z_w),\
                             tools.vinterp(TVmix[:,:,:,it],depths,z_r,z_w)] 
                             
                [THdiff, TVmix] = [THdiff0, TVmix0]
                
            else:
           
                for it in [0,1]:
                    [THdiff[:,:,:len(depths),it],TVmix[:,:,:len(depths),it]]=\
                [tools.vinterp(THdiff[:,:,:,it],depths,z_r,z_w),\
                tools.vinterp(TVmix[:,:,:,it],depths,z_r,z_w)]
            
            ########################################################    
            
            
            buoy = tools.vinterp(buoy,depths,z_r,z_w)[:,:,:len(depths)]

            bHdiff =np.zeros((TVmix.shape[0],TVmix.shape[1],len(depths)))*np.nan
            bVmix = np.zeros((TVmix.shape[0],TVmix.shape[1],len(depths)))*np.nan
            
            #print 'temporary modification alpha =1/g, beta=0'
            #alpha = 1./simul.g; beta=0.

            for iz in range(len(depths)):
                bHdiff[:,:,iz] = simul.g * (alpha*THdiff[:,:,iz,0] - beta*THdiff[:,:,iz,1])
                bVmix[:,:,iz] = simul.g * (alpha*TVmix[:,:,iz,0] - beta*TVmix[:,:,iz,1])
 
            #print 'temporary modification alpha =1/g, beta=0'          
            #buoy = t[:,:,:len(depths),0]
            
            del THdiff,TVmix            
            
            ########################################################       
            
            tenddv = tools.u2rho(tools.diffx(bVmix,pm)*tools.diffx(buoy,pm)) + tools.v2rho(tools.diffy(bVmix,pn)*tools.diffy(buoy,pn))
            tenddh = tools.u2rho(tools.diffx(bHdiff,pm)*tools.diffx(buoy,pm)) + tools.v2rho(tools.diffy(bHdiff,pn)*tools.diffy(buoy,pn))

            tendw = -1*tools.u2rho(tools.diffx(w,pm)*tools.diffx(buoy,pm))*bz - tools.v2rho(tools.diffy(w,pn)*tools.diffy(buoy,pn))*bz
            tendw = tools.nanbnd(tendw,2)
                  
            ########################################################       

            Vmix1 = tools.u2rho(tools.diffx(AKtz,pm)*tools.diffx(buoy,pm))*bz + tools.v2rho(tools.diffy(AKtz,pn)*tools.diffy(buoy,pn))*bz
            Vmix2 = tools.u2rho(tools.diffx(bz,pm)*tools.diffx(buoy,pm))*AKtz + tools.v2rho(tools.diffy(bz,pn)*tools.diffy(buoy,pn))*AKtz
            Vmix3 = tools.u2rho(tools.diffx(AKt,pm)*tools.diffx(buoy,pm))*bzz + tools.v2rho(tools.diffy(AKt,pn)*tools.diffy(buoy,pn))*bzz
            Vmix4 = tools.u2rho(tools.diffx(bzz,pm)*tools.diffx(buoy,pm))*AKt + tools.v2rho(tools.diffy(bzz,pn)*tools.diffy(buoy,pn))*AKt
            
            #Vmix1 = bz
            #Vmix2 = tools.u2rho(tools.diffx(bz,pm))
            #Vmix3 = bzz
            #Vmix4 = tools.u2rho(tools.diffx(bzz,pm)*tools.diffx(buoy,pm))
            
            ########################################################       
            # Compute other terms
            ########################################################      
            
            u = tools.vinterp(tools.u2rho(u),depths,z_r,z_w)
            v = tools.vinterp(tools.v2rho(v),depths,z_r,z_w)

            ux = tools.u2rho(tools.diffx(u,pm))
            uy = tools.v2rho(tools.diffy(u,pn))
            vx = tools.u2rho(tools.diffx(v,pm))
            vy = tools.v2rho(tools.diffy(v,pn))

            tendD = -0.5 * (ux+vy) * ( tools.u2rho(tools.diffx(buoy,pm)**2) + tools.v2rho(tools.diffy(buoy,pn)**2) )
            tendE = -0.5 * (ux-vy) * ( tools.u2rho(tools.diffx(buoy,pm)**2) - tools.v2rho(tools.diffy(buoy,pn)**2) )
            tendF = -1. * (vx+uy)  *  tools.u2rho(tools.diffx(buoy,pm)) * tools.v2rho(tools.diffy(buoy,pn))
             
             
            tendD = tools.nanbnd(tendD,2)
            tendE = tools.nanbnd(tendE,2)
            tendF = tools.nanbnd(tendF,2)
            

            ########################################################    

        
        elif np.min(depths)>0:
        
            #Compuite on sigma-levels
            print('not tested yet!!!')
            ########################################################   
            # We first prepare w, bz cause we need full 3D buoy and u,v
                        
            w = toolsF.get_wvlcty(u,v,z_r,z_w,pm,pn)
            
            bz = (buoy[:,:,1:] - buoy[:,:,:-1])/ (z_r[:,:,1:] -z_r[:,:,:-1])
            bz = tools.vinterp(bz,z_r,z_w[1:-1],z_r)           

            
            ########################################################   

            bHdiff =np.zeros((TVmix.shape[0],TVmix.shape[1],TVmix.shape[2]))*np.nan
            bVmix = np.zeros((TVmix.shape[0],TVmix.shape[1],TVmix.shape[2]))*np.nan
                
            for iz in range(TVmix[:,:,:,0].shape[2]):
                bHdiff[:,:,iz] = simul.g * (alpha*THdiff[:,:,iz,0] - beta*THdiff[:,:,iz,1])
                bVmix[:,:,iz] = simul.g * (alpha*TVmix[:,:,iz,0] - beta*TVmix[:,:,iz,1])
                
            del THdiff,TVmix
            
            ########################################################

            tenddv = tools.u2rho(tools.diffxi(bVmix,pm,z_r,z_w)*tools.diffxi(buoy,pm,z_r,z_w)) + tools.v2rho(tools.diffeta(bVmix,pn,z_r,z_w)*tools.diffeta(buoy,pn,z_r,z_w))
            tenddh = tools.u2rho(tools.diffxi(bHdiff,pm,z_r,z_w)*tools.diffxi(buoy,pm,z_r,z_w)) + tools.v2rho(tools.diffeta(bHdiff,pn,z_r,z_w)*tools.diffeta(buoy,pn,z_r,z_w))
            
            tendw = -1*tools.u2rho(tools.diffxi(w,pm,z_r,z_w)*tools.diffxi(buoy,pm,z_r,z_w))*bz \
                     - tools.v2rho(tools.diffeta(w,pn,z_r,z_w)*tools.diffeta(buoy,pn,z_r,z_w))*bz

            ########################################################
            
            
            ########################################################       
            # Compute other terms
            ########################################################      

            ux = tools.u2rho(tools.diffxi(tools.u2rho(u),pm,z_r,z_w))
            uy = tools.v2rho(tools.diffeta(tools.u2rho(u),pn,z_r,z_w))
            vx = tools.u2rho(tools.diffxi(tools.v2rho(v),pm,z_r,z_w))
            vy = tools.v2rho(tools.diffeta(tools.v2rho(v),pn,z_r,z_w))

            tendD = -0.5 * (ux+vy) * ( tools.u2rho(tools.diffxi(buoy,pm,z_r,z_w)**2) + tools.v2rho(tools.diffeta(buoy,pn,z_r,z_w)**2) )
            tendE = -0.5 * (ux-vy) * ( tools.u2rho(tools.diffxi(buoy,pm,z_r,z_w)**2) - tools.v2rho(tools.diffeta(buoy,pn,z_r,z_w)**2) )
            tendF = -1. * (vx+uy)  *  tools.u2rho(tools.diffxi(buoy,pm,z_r,z_w)) * tools.v2rho(tools.diffeta(buoy,pn,z_r,z_w))       
            
            ########################################################
            
                                
        return [tendD[addx1:addx2,addy1:addy2,:],\
                tendE[addx1:addx2,addy1:addy2,:],\
                tendF[addx1:addx2,addy1:addy2,:],\
                tendw[addx1:addx2,addy1:addy2,:],\
                tenddv[addx1:addx2,addy1:addy2,:],\
                tenddh[addx1:addx2,addy1:addy2,:],\
                Vmix1[addx1:addx2,addy1:addy2,:],\
                Vmix2[addx1:addx2,addy1:addy2,:],\
                Vmix3[addx1:addx2,addy1:addy2,:],\
                Vmix4[addx1:addx2,addy1:addy2,:],\
                ]
                    
           
           
           

#######################################################
#Compute tendency terms fron the tracer equation
#######################################################

    def get_advtend(self,simul,**kwargs):

        #######################################################
        if 'pert' in  kwargs: 
            pert = kwargs['pert']
        else:
            pert=None
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            print('original coord', coord)
            
            # We need at leat 4 pts to compute advection:
            addx1,addy1,addx2,addy2 =0,0,None,None
            if (ny2i-ny1i)<20: 
                addy=20-(ny2i-ny1i); addy1= (addy/2); addy2=-(addy/2+addy%2)
                ny1i=ny1i-addy1; ny2i=ny2i-addy2
            if (nx2i-nx1i)<20: 
                addx=20-(nx2i-nx1i); addx1= (addx/2); addx2=-(addx/2+addx%2)
                nx1i=nx1i-addx1; nx2i=nx2i-addx2      
            coord = [ny1i,ny2i,nx1i,nx2i]        
            print('temporary coord', coord)

            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        else:
             # We can't enlarge the domain...
            addx1,addy1,addx2,addy2 =0,0,None,None
            
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            
        ########################################################    
         
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
            
        #######################################################
        
        depths = simul.coordmax[4]

        [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)
        
        temp = self.load('temp',ncfile,simul,coord=coord,depths=depths)
        salt = self.load('salt',ncfile,simul,coord=coord,depths=depths)
        buoy = -simul.g * toolsF.rho_eos(temp,salt,z_r,z_w,simul.rho0)/simul.rho0
        
        u = self.load('u',ncfile,simul,coord=coord,depths=depths)
        v = self.load('v',ncfile,simul,coord=coord,depths=depths)
        
        #######################################################

        if 'depths' in  kwargs: 
            depths = kwargs['depths']
            
        ########################################################            
        
        buoy = tools.vinterp(buoy,depths,z_r,z_w)[:,:,:len(depths)]        
        u = tools.vinterp(u,depths,tools.rho2u(z_r),tools.rho2u(z_w))
        v = tools.vinterp(v,depths,tools.rho2v(z_r),tools.rho2v(z_w))
        
        ud,vd = tools_g.div2uvs(u,v,pm,pn)
        ur,vr = u-ud, v-vd
        
        ########################################################      

        tenddiv = tools_g.get_tendency(ud,vd,buoy,pm,pn)        
        tendrot = tools_g.get_tendency(ur,vr,buoy,pm,pn)       

   
        ########################################################      
  
                                
        return [tenddiv[addx1:addx2,addy1:addy2,:],\
                tendrot[addx1:addx2,addy1:addy2,:]]
                    
           
           
           
           
                
#######################################################
#Comput
#######################################################



    def get_uv_evolution(self,simul,rhs=False,**kwargs):
    

        #######################################################
        if 'pert' in  kwargs: 
            pert = kwargs['pert']
        else:
            pert=None
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            print('original coord', coord)
            
            # We need at leat 4 pts to compute advection:
            addx1,addy1,addx2,addy2 =0,0,None,None
            #if (ny2i-ny1i)<6: 
                #addy=6-(ny2i-ny1i); addy1= (addy/2); addy2=-(addy/2+addy%2)
                #ny1i=ny1i-addy1; ny2i=ny2i-addy2
            #if (nx2i-nx1i)<6: 
                #addx=6-(nx2i-nx1i); addx1= (addx/2); addx2=-(addx/2+addx%2)
                #nx1i=nx1i-addx1; nx2i=nx2i-addx2      
            #coord = [ny1i,ny2i,nx1i,nx2i]        
            #print 'temporary coord', coord

            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        else:
             # We need at leat 4 pts to compute advection:
             # But we can't enlarge the domain...
            addx1,addy1,addx2,addy2 =0,0,None,None
            
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            
#        if 'recomputeakv' in  kwargs: 
#            recomputeakv = kwargs['recomputeakv']

        if 'myakv' in  kwargs: 
            myakv = kwargs['myakv']
        else: 
            myakv=False


        rmask = copy(mask)
        rmask[np.isnan(mask)] = 0
           
        #if 'depths' in  kwargs: 
            #depths = kwargs['depths']
        #else: 
            #depths = simul.coord[4]
            
        depths = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

        #######################################################
 

        print('collect 0'); collect(); del garbage[:];

        ncfile = Dataset(simul.ncfile, 'r')

        print('collect 1'); collect(); del garbage[:];

        #######################################################

        [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)
        
        u = self.load('u',ncfile,simul,coord=coord,depths=depths,masked=0.)
        if pert==1: u = u - np.mean(u); print('removing mean velocity')
        elif pert==0: u = u * np.nan + np.mean(u); print('using only mean velocity')
        v = self.load('v',ncfile,simul,coord=coord,depths=depths,masked=0.)
        if pert==1: v = v - np.mean(v)        
        elif pert==0: v = v * np.nan + np.mean(v); print('using only mean velocity')

        print('collect 2'); collect(); del garbage[:];

        try:
            omega = self.load('omega',ncfile,simul,coord=coord,depths=depths_w,masked=0.)
        except:
            print('no omega in file, computing')
            omega = toolsF.get_omega (u,v,z_r,z_w,pm,pn) 

        print('collect 3'); collect(); del garbage[:];

            
        try:          
            AKv = self.load('AKv',ncfile,simul,coord=coord,depths=depths_w)
        except:          
            print('no AKv!!!! ....')
            if myakv=='mlconvec':
                print('using MLCONVEC')
                (AKv,AKt) = self.get_akt(simul,coord=coord)     
            elif myakv=='akt':
                print('using AKt for AKv')              
                (AKt,AKv) = self.get_akt(simul,coord=coord)    
            elif myakv=='realakt':
                print('using real AKt')              
                AKv = self.load('AKt',ncfile,simul,coord=coord,depths=depths_w)                 
            else :
                print('NOT using MLCONVEC')              
                AKv = self.get_akv(simul,coord=coord)

        #hbls = self.load('hbls',ncfile,simul,coord=coord)

#        if recomputeakv:
#            print 'force recompute of AKv'
#            AKv = self.get_akv(simul,coord=coord)

        #t=np.zeros((u.shape[0]+1,u.shape[1],u.shape[2],2))
        T = self.load('temp',ncfile,simul,coord=coord,depths=depths,masked=0.)
        try:
            S = self.load('salt',ncfile,simul,coord=coord, depths=depths,masked=0.)
        except:
            print('no S in file')
            S = T*0.        
        
        #Hz =  z_w[:,:,1:] - z_w[:,:,:-1]
        #stflx=np.zeros((u.shape[0]+1,u.shape[1],2))
        #[stflx[:,:,0], stflx[:,:,1],srflx] = self.get_buoy_flux(simul,coord=coord,solar=True)
        #[alpha,beta] = self.alphabeta(t[:,:,-1,0],t[:,:,-1,1],simul.rho0)

        print('collect 4'); collect(); del garbage[:];    

        try:
            #print "TEST NEGATIVE WINDS!!!!"
            [sustr,svstr] = self.get_winds(simul,coord=coord)
            sustr=sustr/simul.rho0; svstr=svstr/simul.rho0
        except:
            print('Cannot find the winds')
            sustr, svstr = np.zeros(u[:,:,-1].shape), np.zeros(v[:,:,-1].shape)
            
        visc2 = simul.visc2
        if not simul.visc2>0: visc2 = 0.

        try:
            v_sponge = simul.v_sponge
            visctype = simul.visctype #Old version (=0) or new version (=1)?
        except:
            print('no v_sponge in file')
            v_sponge = 0.
            visctype = 0 
            
        print('collect 5'); collect(); del garbage[:];

        print('v_sponge is ', v_sponge, ' visctype is ', visctype)
        (MXadv,MYadv,MVadv,MHdiss,MHmix,MVmix,MCor,MPrsgrd) = \
                                toolsF_g.get_uv_evolution (u,v,T,S, z_r,z_w,pm,pn,f\
                                                         ,simul.dt_model,rmask,simul.rdrg,simul.rho0\
                                                         ,omega,AKv,sustr,svstr ,visc2, v_sponge\
                                                         ,visctype\
                                                         , coord[:4], simul.coordmax[:4]   ) 
                            
        #(MXadv,MYadv,MVadv,MHdiss,MHmix,MVmix,MCor,MPrsgrd) = \
                                #toolsF_g.get_uv_evolution_old (u,v,T,S, z_r,z_w,pm,pn,f\
                                                         #,simul.dt_model,rmask,simul.rdrg,simul.rho0\
                                                         #,omega,AKv,sustr,svstr ,visc2, v_sponge\
                                                         #, simul.coord[:4], simul.coordmax[:4]   ) 
        print('MXadv.shape', MXadv.shape)
        print('collect 6'); collect(); del garbage[:];  

        
        
        if len(MXadv.shape)==4:
            [MXadv0,MYadv0,MVadv0,MHdiss0,MHmix0,MVmix0,MCor0,MPrsgrd0] = \
                                [MXadv[addx1:addx2,addy1:addy2,:,:],\
                                 MYadv[addx1:addx2,addy1:addy2,:,:],\
                                 MVadv[addx1:addx2,addy1:addy2,:,:],\
                                 MHdiss[addx1:addx2,addy1:addy2,:,:],\
                                 MHmix[addx1:addx2,addy1:addy2,:,:],\
                                 MVmix[addx1:addx2,addy1:addy2,:,:],\
                                 MCor[addx1:addx2,addy1:addy2,:,:],\
                                 MPrsgrd[addx1:addx2,addy1:addy2,:,:]]
        elif len(MXadv.shape)==3:
            [MXadv0,MYadv0,MVadv0,MHdiss0,MHmix0,MVmix0,MCor0,MPrsgrd0] = \
                                [MXadv[addx1:addx2,addy1:addy2,:],\
                                 MYadv[addx1:addx2,addy1:addy2,:],\
                                 MVadv[addx1:addx2,addy1:addy2,:],\
                                 MHdiss[addx1:addx2,addy1:addy2,:],\
                                 MHmix[addx1:addx2,addy1:addy2,:],\
                                 MVmix[addx1:addx2,addy1:addy2,:],\
                                 MCor[addx1:addx2,addy1:addy2,:],\
                                 MPrsgrd[addx1:addx2,addy1:addy2,:]]

        if rhs:
            return MHdiss0 + MHmix0 + MVmix0 + MCor0 + MPrsgrd0
        else:
            return [MXadv0,MYadv0,MVadv0,MHdiss0,MHmix0,MVmix0,MCor0,MPrsgrd0]



#######################################################
#
#######################################################



    def get_pv(self,simul,**kwargs):
    

        #######################################################

            
        if 'coord' in  kwargs: 
        
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   

            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            
        else:
            
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            

        depths = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))     
            
        #######################################################
        
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')

        #######################################################

            
        T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
        S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
        u = self.load('u',ncfile,simul,coord=coord,depths=depths)
        v = self.load('v',ncfile,simul,coord=coord,depths=depths)

        [z_r,z_w] = tools.get_depths(simul,coord=coord)

        [pv1,pv2,pv3] = tools.PV_terms(T,S,u,v,z_r,z_w,f,simul.g,simul.rho0,pm,pn) 


        #######################################################
            
        return [pv1,pv2,pv3]

        #######################################################
    
            
        
        
        
        
        
                 
                
#######################################################
#Compute Mixed layer depth using lmd.kpp.F
#######################################################



    def get_hbl(self,simul,**kwargs):
    

        #######################################################
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   

            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        else:
            
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            
        if 'idealized' in  kwargs:
            idealized= kwargs['idealized']
        else: 
            idealized= False
            
        #######################################################
        
        rmask = copy(mask)
        rmask[np.isnan(mask)] = 0

        #######################################################

        depths = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

        #######################################################
 
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
            
        #######################################################
        
        u = self.load('u',ncfile,simul,coord=coord,depths=depths)
        v = self.load('v',ncfile,simul,coord=coord,depths=depths)
        print(coord,coord)
        
        
        
        if idealized:               
        #######################################################
        # Idealized forcings...
        #######################################################
            print('using idealized values')

            sustr=np.zeros((u.shape[0],u.shape[1]))
            svstr=np.zeros((v.shape[0],v.shape[1]))          
            stflx=np.zeros((u.shape[0]+1,u.shape[1],2))
            srflx=np.zeros((u.shape[0]+1,u.shape[1]))
            
            
            sustr=sustr*0.-0.0315/simul.rho0; svstr=svstr*0.-0.0124/simul.rho0     
             
            stflx[:,:,0] = stflx[:,:,0] * 0. - 7.5e-5
            stflx[:,:,1] = stflx[:,:,1] * 0.      
            srflx = srflx * 0.     

        else:   
            [sustr,svstr] = self.get_winds(simul,coord=coord)
            sustr=sustr/simul.rho0; svstr=svstr/simul.rho0
            
            stflx=np.zeros((u.shape[0]+1,u.shape[1],2))
            [stflx[:,:,0], stflx[:,:,1],srflx] = self.get_buoy_flux(simul,coord=coord,solar=True)

        print('u.shape, v.shape', u.shape, v.shape)
        
        #######################################################
        #t=np.zeros((u.shape[0]+1,u.shape[1],u.shape[2],2))
        T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
        
        if idealized:
            S = T*0.+35.         
            [alpha,beta] = self.alphabeta(T[:,:,-1],S[:,:,-1],simul.rho0)
            Tcoef = 0.2
            Scoef = 0.
            alpha = alpha*0.+ Tcoef
            beta = beta*0. + Scoef
            R0 = 0. # R0 = alpha * T0
        else:
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
            [alpha,beta] = self.alphabeta(T[:,:,-1],S[:,:,-1],simul.rho0) 
        
        
        [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)
        Hz =  z_w[:,:,1:] - z_w[:,:,:-1]
        (swr_frac) = \
                toolsF.get_swr_frac (Hz)    
        del Hz
        
        if idealized:
            bvf = toolsF.bvf_lineos(T,S,Tcoef,Scoef,z_r,R0,simul.rho0)
        else:
            bvf = toolsF.bvf_eos(T,S,z_r,z_w,simul.rho0)      
        
        hbls = self.load('hbls',ncfile,simul,coord=coord)
        
        #print 'previous hbl temporary set to 0!!!!!!'
        #f = 0.*f
        #print 'f temporary set to 0!!!!!!'       
        #######################################################

        

         
        #######################################################
        
        try:
            Ricr = simul.Ricr
            lmdkpp = simul.lmdkpp
        except:
            if simul.simul in ['atlbig','gulfz','nesea','neseb','nefro']:
                Ricr=0.15
                lmdkpp = 1 #new version 
            else:
                Ricr=0.45 # for SHING family
                lmdkpp = 0 #old version (no bkpp, molemaker patch)
                
        #######################################################

        if  lmdkpp==0:
            print('using old lmd_kppp version for CUC')
      
            (hbl) = toolsF_cuc.get_hbl_cuc (alpha, beta, z_r, z_w, stflx, srflx, swr_frac, sustr, svstr, Ricr, hbls, f,\
                                u, v, bvf,  rmask)
        else:
            
            (hbl, out1, out2, out3, out4) = toolsF_g.get_hbl (alpha, beta, z_r, z_w, stflx, srflx, swr_frac, sustr, svstr, Ricr, hbls, f,\
                                u, v, bvf)
    
    
             
        #######################################################

        #return hbl,hbls,out1,out2,out3,out4
        return hbl

        

       
                
#######################################################
#Compute Mixed layer depth using lmd.kpp.F
#######################################################



    def get_drag(self,simul,**kwargs):
    

        #######################################################
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   

            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        else:
            
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask

        #######################################################
        
        rmask = copy(mask)
        rmask[np.isnan(mask)] = 0

        #######################################################
 
        ncfile = Dataset(simul.ncfile, 'r')
            
        #######################################################
        
        
        u = self.load('u',ncfile,simul,coord=coord,depths=[1])
        v = self.load('v',ncfile,simul,coord=coord,depths=[1])
        [z_r,z_w] = tools.get_depths(simul,coord=coord)
        Hz = z_w[:,:,1] - z_w[:,:,0]

        if simul.ncname.model in ['croco'] and False:
            # Needs some test
            print('using croco version for bottom drag')
            Hz2 = z_r[:,:,0] - z_w[:,:,0]
            (ut,vt) = toolsF_g.get_bot_croco(u,v,Hz,Hz2,simul.rdrg,simul.rdrg2,simul.Zob,\
                                             simul.Cdb_min,simul.Cdb_max,simul.dt_model)

        else:
            (ut,vt) = toolsF_g.get_bot(u,v,Hz,simul.rdrg,simul.Zob)

        #######################################################


        return ut,vt

       
                 
                
#######################################################
#Compute AKv coefficient using lmd.kpp.F
#######################################################



    def get_akv(self,simul,**kwargs):
    

        #######################################################
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   

            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        else:
            
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            
        #######################################################
        
        rmask = copy(mask)
        rmask[np.isnan(mask)] = 0

        #######################################################

        depths = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

        #######################################################
 
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
            
        #######################################################
        
        u = self.load('u',ncfile,simul,coord=coord,depths=depths)
        v = self.load('v',ncfile,simul,coord=coord,depths=depths)
        print(coord,coord)
        
        
        #######################################################


        [sustr,svstr] = self.get_winds(simul,coord=coord)
        sustr=sustr/simul.rho0; svstr=svstr/simul.rho0
        
        stflx=np.zeros((u.shape[0]+1,u.shape[1],2))
        [stflx[:,:,0], stflx[:,:,1],srflx] = self.get_buoy_flux(simul,coord=coord,solar=True)

        
        #######################################################
        #t=np.zeros((u.shape[0]+1,u.shape[1],u.shape[2],2))
        
        T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
        S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
        [alpha,beta] = self.alphabeta(T[:,:,-1],S[:,:,-1],simul.rho0) 
        
        
        [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)
        Hz =  z_w[:,:,1:] - z_w[:,:,:-1]
        (swr_frac) = \
                toolsF.get_swr_frac (Hz)    
        del Hz
        
        bvf = toolsF.bvf_eos(T,S,z_r,z_w,simul.rho0)      
        
        hbls = self.load('hbls',ncfile,simul,coord=coord)
             
        #######################################################

        

         
        #######################################################
        
        try:
            Ricr = simul.Ricr
            lmdkpp = simul.lmdkpp
        except:
            if simul.simul in ['atlbig','gulfz','nesea','neseb','nefro']:
                Ricr=0.15
                lmdkpp = 1 #new version 
            else:
                Ricr=0.45 # for SHING family
                lmdkpp = 0 #old version (no bkpp, molemaker patch)

        rd = rmask*0.
        print('no rd yet!!!!!!!!!')

        if  lmdkpp==0:
            print('using old lmd_kppp version')
            #(akv) = toolsF_g.get_akv_old (alpha, beta, z_r, z_w, stflx, srflx, swr_frac, sustr, svstr, Ricr, hbls, f,\
                                #u, v, bvf,  rmask, rd)
            (akv) = toolsF_cuc.get_akv_cuc (alpha, beta, z_r, z_w, stflx, srflx, swr_frac, sustr, svstr, Ricr, hbls, f,\
                                u, v, bvf,  rmask, rd)            
        else:
            print('using latest lmd_kpp version')
            (akv) = toolsF_g.get_akv (alpha, beta, z_r, z_w, stflx, srflx, swr_frac, sustr, svstr, Ricr, hbls, f,\
                                u, v, bvf,  rmask, rd)
            
        #######################################################

        return akv

        
#######################################################
#Compute AKv coefficient using lmd.kpp.F
#######################################################



    def get_akt(self,simul,**kwargs):
    

        #######################################################
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   

            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        else:
            
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            
        #######################################################
        
        rmask = copy(mask)
        rmask[np.isnan(mask)] = 0

        #######################################################

        depths = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

        #######################################################
 
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
            
        #######################################################
        
        u = self.load('u',ncfile,simul,coord=coord,depths=depths)
        v = self.load('v',ncfile,simul,coord=coord,depths=depths)
        print(coord,coord)
        
        
        #######################################################


        [sustr,svstr] = self.get_winds(simul,coord=coord)
        sustr=sustr/simul.rho0; svstr=svstr/simul.rho0
        
        stflx=np.zeros((u.shape[0]+1,u.shape[1],2))
        [stflx[:,:,0], stflx[:,:,1],srflx] = self.get_buoy_flux(simul,coord=coord,solar=True)

        
        #######################################################
        #t=np.zeros((u.shape[0]+1,u.shape[1],u.shape[2],2))
        
        T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
        S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
        [alpha,beta] = self.alphabeta(T[:,:,-1],S[:,:,-1],simul.rho0) 
        
        
        [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)
        Hz =  z_w[:,:,1:] - z_w[:,:,:-1]
        (swr_frac) = \
                toolsF.get_swr_frac (Hz)    
        del Hz
        
        bvf = toolsF.bvf_eos(T,S,z_r,z_w,simul.rho0)      
        
        hbls = self.load('hbls',ncfile,simul,coord=coord)
             
        #######################################################

        

         
        #######################################################
        
        try:
            Ricr = simul.Ricr
            lmdkpp = simul.lmdkpp
        except:
            if simul.simul in ['atlbig','gulfz','nesea','neseb','nefro']:
                Ricr=0.15
                lmdkpp = 1 #new version 
            else:
                Ricr=0.45 # for SHING family
                lmdkpp = 0 #old version (no bkpp, molemaker patch)

        rd = rmask*0.
        print('no rd yet!!!!!!!!!')

        if  lmdkpp==0:
            (akv,akt) = toolsF_cuc.get_akt_cuc (alpha, beta, z_r, z_w, stflx, srflx, swr_frac, sustr, svstr, Ricr, hbls, f,\
                                u, v, bvf,  rmask, rd)
        else:
            print('nonononono')
            
        #######################################################

        return akv,akt

        
#######################################################
#Compute vertical velocity due to mixing (in the GaLo81 way)
#######################################################

    def get_w_vmix(self,simul,frontal=False, thermal = True, **kwargs):
        '''
        Solving the equation :
        
        Frontal= True means we output u,v,wu and wv
        Frontal = False mean we just oputput w

        '''

        #######################################################
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        else:           
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            
        #######################################################
            
        depths_r = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

        #######################################################
 
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
            
        [z_r,z_w] = tools.get_depths(simul,coord=coord)
        T = self.load('temp',ncfile,simul,coord=coord, depths=depths_r)
        S = self.load('salt',ncfile,simul,coord=coord, depths=depths_r)
        
        #######################################################        
        
        depths = kwargs['depths']   

        #######################################################        

        depths = kwargs['depths']
        depths0 = np.array(depths)-5.
        depths1 = np.array(depths)+5.
        
        buoy = tools.vinterp(-simul.g * toolsF.rho1_eos(T,S,simul.rho0)/simul.rho0,depths,z_r,z_w)
        AKv = tools.vinterp(self.load('AKv',ncfile,simul,coord=coord, depths=depths_w),depths,z_r,z_w,kmin=0) #*0.+4e-2
        
        buoy0 = tools.vinterp(-simul.g * toolsF.rho1_eos(T,S,simul.rho0)/simul.rho0,depths0,z_r,z_w)
        AKv0 = tools.vinterp(self.load('AKv',ncfile,simul,coord=coord, depths=depths_w),depths0,z_r,z_w,kmin=0) #*0.+4e-2
        
        buoy1 = tools.vinterp(-simul.g * toolsF.rho1_eos(T,S,simul.rho0)/simul.rho0,depths1,z_r,z_w)
        AKv1 = tools.vinterp(self.load('AKv',ncfile,simul,coord=coord, depths=depths_w),depths1,z_r,z_w,kmin=0) #*0.+4e-2
        
        #print 'AKv is', AKv
        hbls = self.load('hbls',ncfile,simul,coord=coord)
        
        
        for iz in range(len(depths)):
            AKv[:,:,iz][depths[iz]<-hbls] = 1e-4
            AKv0[:,:,iz][depths0[iz]<-hbls] = 1e-4       
            AKv1[:,:,iz][depths1[iz]<-hbls] = 1e-4
            
        ################################### 
        
                        
        wu = np.zeros(buoy.shape)*np.nan
        wv = np.zeros(buoy.shape)*np.nan
        
        if thermal:
        
            uz = -1*(tools.diffy(buoy,pn)); vz = (tools.diffx(buoy,pm))
            
            if frontal:
                uz0 = -1*(tools.diffy(buoy0,pn)); vz0 = (tools.diffx(buoy0,pm))
                uz1 = -1*(tools.diffy(buoy1,pn)); vz1 = (tools.diffx(buoy1,pm))
                        
            byy = tools.diffy(tools.rho2v(AKv) * uz,tools.rho2v(pn))
            bxx = -1*tools.diffx(tools.rho2u(AKv) * vz,tools.rho2u(pm))     
                                        
            wu = np.zeros(buoy.shape)*np.nan
            wv = np.zeros(buoy.shape)*np.nan
            
            wu[1:-1,:,:] = bxx 
            wv[:,1:-1,:] = byy   
            
            wu = (wu.T/(simul.f**2).T).T
            wv = (wv.T/(simul.f**2).T).T
        else:
            
            u = self.load('u',ncfile,simul,coord=coord,depths=depths_r)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths_r)
            uz = (u[:,:,1:] - u[:,:,:-1])/(tools.rho2u(z_r[:,:,1:])-tools.rho2u(z_r[:,:,:-1]))
            vz = (v[:,:,1:] - v[:,:,:-1])/(tools.rho2v(z_r[:,:,1:])-tools.rho2v(z_r[:,:,:-1]))
            uz = tools.vinterp(uz,depths,tools.rho2u(z_w[:,:,1:-1]),tools.rho2u(z_r))
            vz = tools.vinterp(vz,depths,tools.rho2v(z_w[:,:,1:-1]),tools.rho2v(z_r))
            
            byy = tools.diffy(tools.rho2u(AKv) * uz,tools.rho2u(pn))
            bxx = -1*tools.diffx(tools.rho2v(AKv) * vz,tools.rho2v(pm))
            ##var = np.zeros((u.shape[0],v.shape[1],buoy.shape[2])*np.nan
            wu = tools.psi2rho(bxx + byy)      
            wv = 0.*wu
            
            wu = (wu.T/(simul.f).T).T
            
        if frontal: 
        
            um = (AKv1 * tools.u2rho(vz1) - AKv0 * tools.u2rho(vz0))/10.
            vm = -1*(AKv1 * tools.v2rho(uz1) - AKv0 * tools.v2rho(uz0))/10.

            um = (um.T/(simul.f**2).T).T
            vm = (vm.T/(simul.f**2).T).T

        #else:
            
            #u = self.load('u',ncfile,simul,coord=coord,depths=depths_r)
            #v = self.load('v',ncfile,simul,coord=coord,depths=depths_r)
            #uz = (u[:,:,1:] - u[:,:,:-1])/(tools.rho2u(z_r[:,:,1:])-tools.rho2u(z_r[:,:,:-1]))
            #vz = (v[:,:,1:] - v[:,:,:-1])/(tools.rho2v(z_r[:,:,1:])-tools.rho2v(z_r[:,:,:-1]))
            ##uz[:,:,1:-1] = (u[:,:,2:] - u[:,:,:-2])/(2*dz)
            ##uz[:,:,1] = (u[:,:,1] - u[:,:,0])/(dz)*0; uz[:,:,-1] = (u[:,:,-1] - u[:,:,-2])/(dz)*0
            ##vz[:,:,1:-1] = (v[:,:,2:] - v[:,:,:-2])/(2*dz)
            ##vz[:,:,1] = (v[:,:,1] - v[:,:,0])/(dz)*0; vz[:,:,-1] = (v[:,:,-1] - v[:,:,-2])/(dz)*0       
            #print uz.shape
            #print tools.rho2u(z_w[:,:,1:-1]).shape
            #uz = tools.vinterp(uz,depths,tools.rho2u(z_w[:,:,1:-1]),tools.rho2u(z_r)) 
            #vz = tools.vinterp(vz,depths,tools.rho2v(z_w[:,:,1:-1]),tools.rho2v(z_r) )          
            #byy = tools.diffy(tools.rho2u(AKv) * uz,tools.rho2u(pn))
            #bxx = -1*tools.diffx(tools.rho2v(AKv) * vz,tools.rho2v(pm))
            ##var = np.zeros((u.shape[0],v.shape[1],buoy.shape[2])*np.nan
            #var = tools.psi2rho(bxx + byy)           
            #if len(depths)==1: var = var[:,:,0]    

        ################################################        
        #add winds
        
        [sustr,svstr] = self.get_winds(simul,coord=coord)     
        print(sustr.shape, svstr.shape)
        wu = (wu.T - tools.psi2rho((- tools.diffx(svstr,tools.rho2v(pm)) )).T/simul.rho0).T
        wv = (wv.T - tools.psi2rho((tools.diffy(sustr,tools.rho2u(pn)) )).T/simul.rho0).T   
        
        ################################################        

        
        if frontal:
          return np.squeeze(wu),np.squeeze(wv),np.squeeze(um),np.squeeze(vm)
        else:
          return np.squeeze(wu+wv)       
                
        ################################################        
        return var
        
        
        
        
        
        
        
        
#######################################################
#Compute horizomtal velocities due to TTW
#######################################################

    def get_uv_ttw(self,simul,**kwargs):
        '''
        
        Solving the TTW equation to get u_ttw, v_ttw (and w_ttw by div)
        
        Given a b fieds we solve:
        
        f.du/dz - d2/dz2(AKv.dvdz) = - d/dy(b)
        -f.dv/dz - d2/dz2(AKv.dudz) = - d/dx(b)
        
        
        marche pas!!!
        
        '''

        #######################################################
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            topo = np.asfortranarray(simul.topo[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])           
        else:           
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            topo = simul.topo
            
        #######################################################
            
        depths_r = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))
        
        ############################################################
        # Load and compute known variables (bx, by, AKv) on z-levels
        ############################################################
         
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
        [z_r,z_w] = tools.get_depths(simul,coord=coord)
        
        T = self.load('temp',ncfile,simul,coord=coord, depths=depths_r)
        S = self.load('salt',ncfile,simul,coord=coord, depths=depths_r)

        depths = kwargs['depths']
        buoy = tools.vinterp(-simul.g * toolsF.rho_eos(T,S,z_r,z_w,simul.rho0)/simul.rho0,depths,z_r,z_w)
        del T,S
        
        bx,by =  np.zeros(buoy.shape)*np.nan,np.zeros(buoy.shape)*np.nan
        bx[1:-1,:,:]= tools.diffx(buoy,pm,2); by[:,1:-1,:] = tools.diffy(buoy,pn,2)
        #bx[1:,:,:]= tools.diffx(buoy,pm); by[:,1:,:] = tools.diffy(buoy,pn)
        
        del buoy
        #xnan = np.isnan(bx)
        #bx[xnan]=0
        #bx[np.isnan(by)]=0
        #by[xnan]=0
        #by[np.isnan(by)]=0
        
        ################################################        

        AKv = self.load('AKv',ncfile,simul,coord=coord, depths=depths_w)
        AKv = tools.vinterp(AKv,depths,z_r,z_w,kmin=0); 
        
        ################################################      
        # Create mask at each depth (we don;t want nan for the solver)
        ################################################      

        mask_depth = np.zeros(bx.shape)*np.nan
        #mask_depth0 = np.zeros(bx.shape)*np.nan        
        for iz in range(len(depths)):
            mask_depth[:,:,iz][-1*topo<depths[iz]]=1.
            #mask_depth0[:,:,iz][-1*topo+100<depths[iz]]=1.   
            
        #mask_depth=tools.psi2rho(tools.rho2psi(mask_depth))
        AKv = AKv*mask_depth
        bx = bx*mask_depth
        by = by*mask_depth
        
        AKv[np.isnan(AKv)]=0
        bx[np.isnan(bx)]=0
        by[np.isnan(by)]=0
        
        ################################################        
        #add winds
        
        [sustr,svstr] = self.get_winds(simul,coord=coord)
        sustr,svstr = tools.u2rho(sustr)/simul.rho0, tools.v2rho(svstr)/simul.rho0
      
            
            
        ################################################        
        # Solve Equation:
        ################################################ 
        timing=True
        
        if timing: tstart = tm.time()   
        #u,v = tools_g.solve_ttwall(bx,by,AKv,sustr,svstr,f,pm,pn,depths)
        
        if timing: print('computation all.........', tm.time()-tstart)
        if timing: tstart = tm.time()   
        
        zw = copy(AKv)
        for i in range(AKv.shape[0]):
          for j in range(AKv.shape[1]):           
            zw[i,j,:] = depths
        
        ut,vt,ug,vg = tools_g.solve_ttw(bx,by,AKv,sustr,svstr,f,pm,pn,zw)
        if timing: print('computation old.........', tm.time()-tstart)
        
        
        ut = ut*mask_depth
        vt = vt*mask_depth        
        ################################################     
        
        return ut,vt,ug,vg
        
        
        
#######################################################
#Compute vertical velocity due to TTW
#######################################################

    def get_w_ttw(self,simul,**kwargs):
        '''
        marche pas!!!
        '''
        #######################################################   
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        else:           
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            
            
            
        #######################################################   
        depths = kwargs['depths']  
        
        ut,vt,ug,vg = self.get_uv_ttw(simul,depths=depths)
        
        w = -1*tools.nanbnd(np.cumsum(tools.u2rho(tools.diffx(ut-ug,simul.pm)) + tools.v2rho(tools.diffy(vt-vg,simul.pn)),2)*(depths[1] - depths[0]),1)
        
        ################################################     
        
        return w        
   
   
   
    
#######################################################
#Compute horizomtal velocities due to TTW
#######################################################

    def get_uv_ttw_sig(self,simul,ekman=1,fast=1,**kwargs):
        '''
        
        Solving the TTW equation to get u_ttw, v_ttw (and w_ttw by div)
        
        Given a b fieds we solve:
        
        f.du/dz - d2/dz2(AKv.dvdz) = - d/dy(b)
        -f.dv/dz - d2/dz2(AKv.dudz) = - d/dx(b)
        '''

        #######################################################
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            topo = np.asfortranarray(simul.topo[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])           
        else:           
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            topo = simul.topo
            
        print('coord', coord)
        print('topo.shape in get_uv_ttw_sig', topo.shape)
        #######################################################
            
        depths_r = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))
        
        ############################################################
        # Load and compute known variables (bx, by, AKv) on z-levels
        ############################################################
         
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
        [z_r,z_w] = tools.get_depths(simul,coord=coord)
        
        T = self.load('temp',ncfile,simul,coord=coord, depths=depths_r)
        S = self.load('salt',ncfile,simul,coord=coord, depths=depths_r)

        buoy = -simul.g * toolsF.rho_eos(T,S,z_r,z_w,simul.rho0)/simul.rho0
        del T,S
        
        print('what the fuck?')

        bx = tools.u2rho(tools.diffxi(buoy,pm,z_r,z_w,z_w[:,:,:]));
        by = tools.v2rho(tools.diffeta(buoy,pn,z_r,z_w,z_w[:,:,:]))
        #bx[1:,:,:]= tools.diffx(buoy,pm); by[:,1:,:] = tools.diffy(buoy,pn)


        del buoy
        #xnan = np.isnan(bx)
        #bx[xnan]=0
        #bx[np.isnan(by)]=0
        #by[xnan]=0
        #by[np.isnan(by)]=0
        
        ################################################        
        #zeta = self.load('zeta',ncfile,simul,coord=coord)      
        AKv = self.load('AKv',ncfile,simul,coord=coord, depths=depths_w)      
        
        #Interpolation to AKv to the free surface give spurious values in computaiton of ttw ( Uz = sustr / AKv)
        #AKv[:,:,-1] = tools.vinterp(AKv[:,:,:-1],zeta,z_w[:,:,:-1])[:,:,-1]; 
        
        #Try interpolation at z_r instead
        AKv[:,:,-1] = tools.vinterp(AKv[:,:,:-1],z_r[:,:,-1],z_w[:,:,:-1])[:,:,-1]; 
        #AKv[:,:,-1] = tools.vinterp(AKv[:,:,:-1],zeta - 0.5*(zeta-z_r[:,:,-1]),z_w[:,:,:-1])[:,:,-1]; 
        ################################################        

        
        AKv[np.isnan(AKv)]=0
        bx[np.isnan(bx)]=0
        by[np.isnan(by)]=0
        
        ################################################        
        #add winds
        
        [sustr,svstr] = self.get_winds(simul,coord=coord)
        sustr,svstr = tools.u2rho(sustr)/simul.rho0, tools.v2rho(svstr)/simul.rho0
      
            
            
        ################################################        
        # Solve Equation:
        ################################################ 
        timing=True
        
        #if timing: tstart = tm.time()   
        #u,v = tools_g.solve_ttwall(bx,by,AKv,sustr,svstr,f,pm,pn,depths)
        #if timing: print 'computation all.........', tm.time()-tstart
        
        if timing: tstart = tm.time()   
        
        if fast==1:
            if ekman==1:
                u,v,ug,vg,uek,vek = tools_g.solve_ttw_sig_fast(bx,by,AKv,sustr,svstr,f,pm,pn,z_w,ekman=1)
            else:
                u,v,ug,vg = tools_g.solve_ttw_sig_fast(bx,by,AKv,sustr,svstr,f,pm,pn,z_w)
        else:
            if ekman==1:
                u,v,ug,vg,uek,vek = tools_g.solve_ttw_sig(bx,by,AKv,sustr,svstr,f,pm,pn,z_w,ekman=1)
            else:
                u,v,ug,vg = tools_g.solve_ttw_sig(bx,by,AKv,sustr,svstr,f,pm,pn,z_w)
            
        if timing: print('computation .........', tm.time()-tstart)
        
        # Return result on vertical rho-grid
        #u = tools.vinterp(u,z_r,z_r,z_w,kmin=0)
        #v = tools.vinterp(v,z_r,z_r,z_w,kmin=0)
        #ug = tools.vinterp(ug,z_r,z_r,z_w,kmin=0)
        #vg = tools.vinterp(vg,z_r,z_r,z_w,kmin=0)
        
        print('get_uv_ttw_sig OK')   

        
        ################################################     
        if ekman==1:
            return u,v,ug,vg,uek,vek
        else:  
            return u,v,ug,vg
    
    
#######################################################
#Compute horizomtal velocities due to TTW
#######################################################

    def get_uv_ttw_sig_fort(self,simul,debug=0,**kwargs):
        '''
        
        Solving the TTW equation to get u_ttw, v_ttw (and w_ttw by div)
        
        Given a b fieds we solve:
        
        f.du/dz - d2/dz2(AKv.dvdz) = - d/dy(b)
        -f.dv/dz - d2/dz2(AKv.dudz) = - d/dx(b)
        '''

        #######################################################
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            topo = np.asfortranarray(simul.topo[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])           
        else:           
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            topo = simul.topo
            
        #######################################################
            
        depths_r = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))
        
        ############################################################
        # Load and compute known variables (bx, by, AKv) on z-levels
        ############################################################
         
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
        [z_r,z_w] = tools.get_depths(simul,coord=coord)
        
        T = self.load('temp',ncfile,simul,coord=coord, depths=depths_r)
        S = self.load('salt',ncfile,simul,coord=coord, depths=depths_r)
    

        buoy = -simul.g * toolsF.rho_eos(T,S,z_r,z_w,simul.rho0)/simul.rho0
        del T,S
        
        buoyz = copy(buoy); N = buoy.shape[-1]
        buoyz[:,:,-1] = simul.g * z_w[:,:,-1] - buoy[:,:,-1] * (z_w[:,:,-1]-z_r[:,:,-1])
        
        for k in range(N-2,-1,-1):
            buoyz[:,:,k] = buoyz[:,:,k+1]\
               - 0.5*(buoy[:,:,k]+buoy[:,:,k+1]) * (z_r[:,:,k+1]-z_r[:,:,k])

        #py.plot(buoyz[10,10,:],z_r[10,10,:]); py.show()
        bx= tools.u2rho(tools.diffxi(buoyz,pm,z_r,z_w,)); 
        by = tools.v2rho(tools.diffeta(buoyz,pn,z_r,z_w))
        bx[np.isnan(bx)]=0
        by[np.isnan(by)]=0
        del buoy

        '''
        rho = toolsF.rho_eos(T,S,z_r,z_w,simul.rho0)/simul.rho0
        rho[np.isnan(rho)]=0
        '''
    
        ################################################        

        AKv = self.load('AKv',ncfile,simul,coord=coord, depths=depths_w)      

        ################################################        

        
        AKv[np.isnan(AKv)]=0

        
        ################################################        
        #add winds
        
        [sustr,svstr] = self.get_winds(simul,coord=coord)
        sustr,svstr = tools.u2rho(sustr)/simul.rho0, tools.v2rho(svstr)/simul.rho0
    
        ################################################        
        # Solve Equation:
        ################################################ 
        timing=True

        if timing: tstart = tm.time()   


        u,v,ug,vg  = toolsF_g.solve_ttw(bx,by,AKv,sustr,svstr,f,pm,pn,z_r,z_w)
            
        if timing: print('computation .........', tm.time()-tstart)

        
        ################################################     
        return u,v,ug,vg
              
#######################################################
#Compute vertical velocity due to TTW
#######################################################

    def get_w_ttw_sig(self,simul,debug=0,fast=1,**kwargs):
        
        #######################################################   
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        else:           
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            
            
        #######################################################   
        if fast==1:
            debug=0;
            ut,vt,ug,vg = self.get_uv_ttw_sig(simul,fast=1,coord=coord)
        elif debug==1:
            ut,vt,ug,vg,uz,vz = self.get_uv_ttw_sig(simul,debug=1,coord=coord)
        else:
            ut,vt,ug,vg = self.get_uv_ttw_sig(simul,coord=coord)

        [z_r,z_w] = tools.get_depths(simul,coord=coord)
        w = toolsF.get_wvlcty(tools.rho2u(tools.nanbnd(ut-ug)),tools.rho2v(tools.nanbnd(vt-vg)),z_r,z_w,pm,pn)
        
        ################################################     

        if debug==1:
            return w,uz,vz
        else:  
            return w

        
#######################################################
#Compute vertical velocity due to TTW
#######################################################

    def get_w_ttw_sig_fort(self,simul,**kwargs):
        
        #######################################################   
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
        else:           
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            
        #######################################################   
        
        ut,vt,ug,vg = self.get_uv_ttw_sig_fort(simul)

        [z_r,z_w] = tools.get_depths(simul,coord=coord)
        w = toolsF.get_wvlcty(ut-ug,vt-vg,z_r,z_w,pm,pn)
        
        ################################################     
        
        return w      
        
#######################################################
#Compute QG omega equation
#######################################################


    def get_w_omega(self,simul,mixrotuv=True,u_ttw=None,v_ttw=None,variable=False,fast=False,bbc=0,wbot=None,**kwargs):

        #######################################################
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            topo = np.asfortranarray(simul.topo[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])           
        else:           
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            topo = simul.topo


        #######################################################
        
        depths_r = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))
        depths = kwargs['depths']

        #######################################################
        # get u,v,buoy on sigmal levels

        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
            
        [z_r,z_w] = tools.get_depths(simul,coord=coord)

        T = self.load('temp',ncfile,simul,coord=coord, depths=depths_r)
        
        if 'NONLIN_EOS' in simul.cpp:
            S = self.load('salt',ncfile,simul,coord=coord, depths=depths_r)
            rho = toolsF.rho_eos(T,S,z_r,z_w,simul.rho0)
        else:
            print('using Linear EOS')
            try:
                S = self.load('salt',ncfile,simul,coord=coord, depths=depths_r)
            except:
                S = 0.
            T0,S0 = 0.,0.
            rho = simul.R0 - simul.Tcoef*(T - T0) + simul.Scoef*(S - S0)

        ########################################################       
        # add mixing tracer term: 
        ########################################################       

        forcing = 0
        
        ########################################################       
        # Interpolate to vertical levels

        buoy = tools.vinterp(-simul.g *rho/simul.rho0,depths,z_r,z_w)

        ########################################################       
   
        bvf = (rho[:,:,1:] - rho[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1] )
        bvf = tools.vinterp(-simul.g *bvf/simul.rho0,depths,z_w[:,:,1:-1],z_r)
        
        #bvf = tools.vinterp(toolsF.bvf_eos(T,S,z_r,z_w,simul.rho0),depths,z_w)
        
        del T,S    
        
        N2 = tools.nanmean(tools.nanmean(bvf,0),0); 
        N2[N2<=0] = 1e-10  
        del bvf
        
        N2[np.isnan(N2)] = 1e-10

        print('--------------------------------------------------------------')        
        print('real N2', N2)
        print('--------------------------------------------------------------')               
        #N2 = 1e-3*np.exp(2*depths/np.max(np.abs(depths)))
        #print 'used N2', N2        

        
        
        if 'field' in  kwargs:
            field = kwargs['field']
        else:
            field='buoy'
            
        print('field is', field)
        print(kwargs)
  
        if field in ['urot','utot','udiv']:
            u = tools.vinterp(self.load('u',ncfile,simul,coord=coord,depths=depths_r),depths,tools.rho2u(z_r),tools.rho2u(z_w))
            v = tools.vinterp(self.load('v',ncfile,simul,coord=coord,depths=depths_r),depths,tools.rho2v(z_r),tools.rho2v(z_w))
            
        del z_r,z_w
        
        #######################################################        
        # We use only the non divergent field
        if field in ['urot']:
            print('using non divergent velocity field')
            ud,vd = tools_g.div2uvs(u,v,pm,pn,variable=variable,fast=fast)
            u,v = u-ud, v-vd
            del ud,vd
            
        #######################################################        
        # We use only the  divergent field
        if field in ['udiv']:
            print('using divergent velocity field')
            u,v = tools_g.div2uvs(u,v,pm,pn,variable=variable,fast=fast)      

        ########################################################        
        ## degrade resolution horizontally by a factor nn
        
        #if 'nh' in  kwargs: 
            #nn=kwargs['nh']
            #buoy =  sm.moy_h(buoy,nn); 
            #f=  sm.moy_h(f,nn); 
            #pm = sm.moy_h(pm,nn)/nn
            #pn = sm.moy_h(pn,nn)/nn
            #print 'nh is', nn
            
            #if len(N2.shape)==3: N2=sm.moy_h(N2,nn);
            
            #if field in ['urot','utot','udiv']:
                #u =  sm.moy_h(tools.u2rho(u),nn);
                #v =  sm.moy_h(tools.v2rho(v),nn);
                
                #print '--------------------------------------------------------------'        
                #print u.shape, v.shape, buoy.shape
                #print '--------------------------------------------------------------'        
                
        #done inside solve_omega instead

        
        ########################################################       
        
        mask_depth = np.zeros(buoy.shape)*np.nan
        for iz in range(len(depths)):
            mask_depth[:,:,iz][-1*topo<depths[iz]]=1.
            mask_depth[:,:,iz]=mask_depth[:,:,iz]*mask
            
        buoy = buoy * mask_depth

        ########################################################       
        print('in get_w_omega mixrotuv is', mixrotuv)
        
        if 'nh' in  kwargs: nh=kwargs['nh']
        
        if field in ['urot','utot','udiv']: 
            u,v = tools.u2rho(u)* mask_depth, tools.v2rho(v)* mask_depth
            print(u.shape, v.shape, buoy.shape)
            w = tools_g.solve_omega(buoy,pm,pn,f,N2,depths,u,v,nh=nh,forcing = forcing,mixrotuv=mixrotuv,bbc=bbc,wbot=wbot)
        elif field in ['uttw']:
            w = tools_g.solve_omega(buoy,pm,pn,f,N2,depths,u_ttw,v_ttw,nh=nh,forcing = forcing,mixrotuv=mixrotuv,bbc=bbc,wbot=wbot)
        else:
            w = tools_g.solve_omega(buoy,pm,pn,f,N2,depths,nh=nh,forcing = forcing,bbc=bbc,wbot=wbot)

        #######################################################        


        return w







      
        
#######################################################
#Compute QG omega equation
#######################################################


    def get_w_genomega(self, simul, mixing = False,  terms=None, terms_data=None, stab=True, debug=False, smoothing=False, nsmooth=1, **kwargs):

        #######################################################
            
        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]   
            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            topo = np.asfortranarray(simul.topo[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])           
        else:           
            coord = simul.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            topo = simul.topo


        #######################################################
        
        depths_r = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))
        depths = kwargs['depths']

        #######################################################
        # get u,v,buoy on sigmal levels

        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
            
        [z_r,z_w] = tools.get_depths(simul,coord=coord)

        T = self.load('temp',ncfile,simul,coord=coord, depths=depths_r)
        S = self.load('salt',ncfile,simul,coord=coord, depths=depths_r)
        rho = toolsF.rho_eos(T,S,z_r,z_w,simul.rho0)
        #rho = toolsF.rho1_eos(T,S,simul.rho0)
        
        
        ########################################################       
        # Interpolate to vertical levels

        buoy = tools.vinterp(-simul.g *rho/simul.rho0,depths,z_r,z_w)
        
        forcing = buoy*0.
        
        ########################################################       
        # Create forcing from momentum and tracer ROMS equation for test
        # if terms['ROMS']==1 all aother terms are pot to 0
        ########################################################       

        if terms['ROMS']==1: 
        
            [alpha,beta] = self.alphabeta(T[:,:,:],S[:,:,:],simul.rho0)
            
            [TXadvtot,TYadvtot,TVadvtot,THdifftot,TVmixtot,TForctot] = self.get_tracer_evolution(simul)
            del TXadvtot,TYadvtot,THdifftot,TVmixtot,TForctot
                        
            TVadvz = tools.vinterp(TVadvtot[:,:,:,0],depths,z_r,z_w)
            SVadvz = tools.vinterp(TVadvtot[:,:,:,1],depths,z_r,z_w)
            del TVadvtot
            
            temp = var('temp',simul,depths=depths).data
            salt = var('salt',simul,depths=depths).data
            u = var('u',simul,depths=depths).data
            v = var('v',simul,depths=depths).data
            div = np.zeros(buoy.shape)
            div[1:-1,:,:] = tools.diffx(u,tools.rho2u(pm))
            div[:,1:-1,:] = div[:,1:-1,:]+tools.diffy(v,tools.rho2v(pn))
            del u,v
            
            [alpha,beta] = self.alphabeta(temp,salt,simul.rho0)
            
            TVadv = tools.laplacien(TVadvz - div * temp,pm,pn)
            SVadv = tools.laplacien(SVadvz - div * salt,pm,pn)         
            del TVadvz, SVadvz
            
            Tbadv = simul.g * (TVadv * alpha - SVadv * beta)
                       
            #TVadv = simul.g * (TVadvtot[:,:,:,0] * alpha - TVadvtot[:,:,:,1] * beta)
            #TVadv = tools.vinterp(TVadv,depths,z_r,z_w)  
            #del TVadvtot 
            
            #forcing = simul.g * (TTmix.T * alpha.T - TSmix.T * beta.T).T
            #div = np.zeros(buoy.shape)
            #u = tools.vinterp(self.load('u',ncfile,simul,coord=coord,depths=depths_r),depths,tools.rho2u(z_r),tools.rho2u(z_w))
            #v = tools.vinterp(self.load('v',ncfile,simul,coord=coord,depths=depths_r),depths,tools.rho2v(z_r),tools.rho2v(z_w))
            #div[1:-1,:,:] = tools.diffx(u,tools.rho2u(pm))
            #div[:,1:-1,:] = div[:,1:-1,:]+tools.diffy(v,tools.rho2v(pn))
            #del u,v

            #buoylin = simul.g * (T.T * alpha.T - S.T * beta.T).T
            #buoylin = tools.vinterp(buoylin,depths,z_r,z_w)
            
            #forcing = -1*tools.laplacien(TVadv - div * buoylin,pm,pn)
            #del div
            
            forcing = -1*Tbadv
            
            #######################################
        
            [MXadv,MYadv,MVadv,MHdiss,MHmix,MVmix,MCor,MPrsgrd] = self.get_uv_evolution(simul)
            del MXadv,MYadv,MVadv,MHdiss,MHmix,MVmix,MPrsgrd
            
            MCoru, MCorv, MCorz = np.zeros(buoy.shape), np.zeros(buoy.shape), np.zeros(buoy.shape)
            
            MCoru[1:-1,1:-1,:] = tools.vinterp(MCor[:,:,:,0] ,depths,z_r[1:-1,1:-1,:],z_w[1:-1,1:-1,:])
            MCorv[1:-1,1:-1,:] = tools.vinterp(MCor[:,:,:,1] ,depths,z_r[1:-1,1:-1,:],z_w[1:-1,1:-1,:])

            MCorz[1:-1,:,:] = tools.diffx(MCorv,pm,2) 
            MCorz[:,1:-1,:] =  MCorz[:,1:-1,:] - tools.diffy(MCoru,pn,2)
            del MCoru, MCorv
            
            MCorz = tools.diffz(MCorz,depths)
            
            forcing = forcing + ( MCorz.T * f.T).T

        
            print('using MCor and TVadv to compute forcing')
            print('all other terms put to zero:', terms)
            
            terms_data['ROMS'] = ( MCorz.T * f.T).T - Tbadv
            
            del MCorz, MCor
        ########################################################       
        # add mixing tracer term: 
        ########################################################          
        
        
        
        if terms['MIXB']==1:

            [aaa,aaa,aaa,THdifftot,TVmixtot,TForctot] = self.get_tracer_evolution(simul)
            del aaa
            
            TTmix = TVmixtot[:,:,:,0] + TForctot[:,:,:,0] + THdifftot[:,:,:,0]
            TSmix = TVmixtot[:,:,:,1] + TForctot[:,:,:,1] + THdifftot[:,:,:,1]
            del TVmixtot, TForctot, THdifftot
            
            [alpha,beta] = self.alphabeta(T[:,:,-1],S[:,:,-1],simul.rho0)
            
            TTmixz = tools.vinterp(TTmix,depths,z_r,z_w)
            TSmixz = tools.vinterp(TSmix,depths,z_r,z_w)

            TTmix = tools.laplacien(TTmixz,pm,pn)
            TSmix = tools.laplacien(TSmixz,pm,pn)         
            del TTmixz, TSmixz
            
            forcing = simul.g * (TTmix.T * alpha.T - TSmix.T * beta.T).T
            
            terms_data['MIXB'] = simul.g * (TTmix.T * alpha.T - TSmix.T * beta.T).T
            
            del TTmix, TSmix

            
            
        ########################################################       
        # add mixing momentum term: 
        ########################################################       

        if terms['MIXUV']==1:
            
            
            [MXadv,MYadv,MVadv,MHdiss,MHmix,MVmix,MCor,MPrsgrd] = self.get_uv_evolution(simul)
            del MXadv,MYadv,MVadv,MCor,MPrsgrd,MHdiss
            
            Mmixu, Mmixv, Mmix = np.zeros(buoy.shape), np.zeros(buoy.shape), np.zeros(buoy.shape)
            
            Mmixu[1:-1,1:-1,:] = tools.vinterp(MHmix[:,:,:,0] + MVmix[:,:,:,0],depths,z_r[1:-1,1:-1,:],z_w[1:-1,1:-1,:])
            Mmixv[1:-1,1:-1,:] = tools.vinterp(MHmix[:,:,:,1] + MVmix[:,:,:,1],depths,z_r[1:-1,1:-1,:],z_w[1:-1,1:-1,:])

            Mmix[1:-1,:,:] = tools.diffx(Mmixv,pm,2) 
            Mmix[:,1:-1,:] =  Mmix[:,1:-1,:] - tools.diffy(Mmixu,pn,2)
            del Mmixu, Mmixv
            
            Mmixz = tools.diffz(Mmix,depths)
            
            forcing = forcing - ( Mmixz.T * f.T).T
            
            terms_data['MIXUV'] = - ( Mmixz.T * f.T).T
            
            del Mmixz, Mmix, MHmix,  MVmix

            
        ########################################################       
        # 
        ########################################################       

        if terms['AVTN']==1: 

            [TXadvtot,TYadvtot,TVadvtot,THdifftot,TVmixtot,TForctot] = self.get_tracer_evolution(simul)
            
            Ttot = TXadvtot+TYadvtot+TVadvtot+TForctot+TVmixtot
            del TXadvtot,TYadvtot,TVadvtot,THdifftot,TVmixtot,TForctot
            
            Tt = tools.vinterp(Ttot[:,:,:,0],depths,z_r,z_w)  
            St = tools.vinterp(Ttot[:,:,:,1],depths,z_r,z_w)
            
            #[alpha,beta] = self.alphabeta(T[:,:,-1],S[:,:,-1],simul.rho0)
            temp = var('temp',simul,depths=depths).data
            salt = var('salt',simul,depths=depths).data
            
            [alpha,beta] = self.alphabeta(temp ,salt,simul.rho0)
            del temp,salt

            Tt = tools.laplacien(Tt,simul.pm,simul.pn)
            St = tools.laplacien(St,simul.pm,simul.pn)

            Ttot =  simul.g * (Tt.T * alpha.T - St.T * beta.T).T

            #Tt = simul.g * (Ttot[:,:,:,0] * alpha - Ttot[:,:,:,1] * beta)
            #Tt = tools.vinterp(Ttot,depths,z_r,z_w)  

            del Tt,St
            
           
            forcing = forcing - Ttot
            
            #######################################
        
            [MXadv,MYadv,MVadv,MHdiss,MHmix,MVmix,MCor,MPrsgrd] = self.get_uv_evolution(simul)
            Mt = MXadv+MYadv+MVadv+MCor+MPrsgrd+MVmix
            del MXadv,MYadv,MVadv,MHdiss,MHmix,MVmix,MPrsgrd,MCor
            
            Mtu, Mtv, Mtz = np.zeros(buoy.shape), np.zeros(buoy.shape), np.zeros(buoy.shape)
            
            Mtu[1:-1,1:-1,:] = tools.vinterp(Mt[:,:,:,0] ,depths,z_r[1:-1,1:-1,:],z_w[1:-1,1:-1,:])
            Mtv[1:-1,1:-1,:] = tools.vinterp(Mt[:,:,:,1] ,depths,z_r[1:-1,1:-1,:],z_w[1:-1,1:-1,:])

            Mtz[1:-1,:,:] = tools.diffx(Mtv,pm,2) 
            Mtz[:,1:-1,:] =  Mtz[:,1:-1,:] - tools.diffy(Mtu,pn,2)
            del Mtu, Mtv
            
            Mtz = tools.diffz(Mtz,depths)
            
            forcing = forcing + ( Mtz.T * f.T).T
            
            terms_data['AVTN'] = - Ttot + ( Mtz.T * f.T).T
            
            
            del Mtz, Mt
 


        ########################################################       
   
        bvf = (rho[:,:,1:] - rho[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1] )
        bvf = tools.vinterp(-simul.g *bvf/simul.rho0,depths,z_w[:,:,1:-1],z_r)
        
        #bvf = tools.vinterp(toolsF.bvf_eos(T,S,z_r,z_w,simul.rho0),depths,z_w)
        
        del T,S    

        N2 = copy(bvf); N2mean = tools.nanmean(tools.nanmean(bvf,0),0); 
        
        if terms['LHS2']==0 or debug:
        #Use N2 mean
            for i in range(N2.shape[0]):
                for j in range(N2.shape[1]):
                    N2[i,j,:] = N2mean; 

        #del bvf
        
        N2[np.isnan(N2)] = 0.
        if stab: N2[N2<=0] = 1e-10   
        if stab: bvf[bvf<=0] = 1e-10  
            
        ########################################################       
 
        if terms['NWE']==1: 
        
            w_ttw=kwargs['w_ttw']
            
            print('N2.shape, w_ttw.shape', N2.shape, w_ttw.shape)
            forcing = forcing + tools.laplacien( N2 * w_ttw,pm,pn)
            terms_data['NWE'] = tools.laplacien( N2 * w_ttw,pm,pn)
            
        if terms['ROMS']==0 or debug:   
        
        ########################################################       
     
            
            if 'field' in  kwargs:
                field = kwargs['field']
            else:
                field='buoy'
                
            print('field is', field)
            print(kwargs)
    
            #if field in ['urot','utot','udiv']:
            u = tools.vinterp(self.load('u',ncfile,simul,coord=coord,depths=depths_r),depths,tools.rho2u(z_r),tools.rho2u(z_w))
            v = tools.vinterp(self.load('v',ncfile,simul,coord=coord,depths=depths_r),depths,tools.rho2v(z_r),tools.rho2v(z_w))
            ur,vr = None, None
            
            #######################################################        
            # We use only the non divergent field
            if field in ['urot']:
                print('using non divergent velocity field')
                ud,vd = tools_g.div2uvs(u,v,pm,pn)
                ur,vr = u-ud, v-vd
                del ud,vd
                
            #######################################################        
            # We use only the  divergent field
            if field in ['udiv']:
                print('using divergent velocity field')
                ur,vr = tools_g.div2uvs(u,v,pm,pn)      

            
            #########################################################       
            
            mask_depth = np.zeros(buoy.shape)*np.nan
            for iz in range(len(depths)):
                mask_depth[:,:,iz][-1*topo<depths[iz]]=1.
                mask_depth[:,:,iz]=mask_depth[:,:,iz]*mask
                
            buoy = buoy * mask_depth
            
            
            u,v = tools.u2rho(u)* mask_depth, tools.v2rho(v)* mask_depth
            
            if field in ['urot','udiv']: 
                ur,vr = tools.u2rho(ur)* mask_depth, tools.v2rho(vr)* mask_depth
            elif field in ['uttw','ugeo','urest']:
                [ur,vr] = kwargs['uv']

        else:

            # if ROMS==1 we don;t need anything except forcing
            u,v,buoy = 0,0,0            
            
            
            
            
        ########################################################       
        
        if 'nh' in  kwargs: nh=kwargs['nh']
        
        ########################################################  
        
        
        if debug:
        ########################################################

        
            new = np.zeros(buoy.shape)*np.nan
            ux,uy,uz = copy(new),copy(new),copy(new)
            vx,vy,vz = copy(new),copy(new),copy(new)  
            ux[1:-1,:,:] = tools.diffx(u,pm,2); uy[:,1:-1,:] = tools.diffy(u,pn,2)   
            vx[1:-1,:,:] = tools.diffx(v,pm,2); vy[:,1:-1,:] = tools.diffy(v,pn,2)         
            vrt = (vx - uy)
            vrt[np.isnan(vrt)]=0        
            del ux,vx,uy,vy
            
            #################################################################

            w = var('w',simul,depths=depths).data   
            
            ########################################################  
  
            terms_data['LHS1']= ( tools.diffzz(w,depths).T * ( f.T + vrt.T) * f.T).T
            
            terms_data['LHS1QG']= ( tools.diffzz(w,depths).T * (f**2).T).T
                
            ########################################################  
            #N2[N2<=0] = 1e-10   
            terms_data['LHS2QG']= tools.laplacien(N2*w,pm,pn)
            
            #bvf[bvf<=0] = 1e-10    
            terms_data['LHS2']= tools.laplacien(bvf*w,pm,pn)
            
            ########################################################  
            
            vrtzz = tools.diffzz(vrt,depths)
            del vrt
            
            terms_data['VA']= - (f.T * vrtzz.T).T * w
            
            wx,wy = copy(new),copy(new)
            wx[1:-1,:,:] = tools.diffx(w,pm,2); wy[:,1:-1,:] = tools.diffy(w,pn,2)   
            vrtzz = (wx * tools.diffzz(v,depths) - wy*tools.diffzz(u,depths))
            terms_data['T2']=  - (f.T * vrtzz.T).T
            del vrtzz         
            
        ########################################################  

        if terms['LHS2']==1: N2 = bvf # Use variable N2:
        
        #if terms['LHS2']==0 and terms['LHS2QG']==0: N2 = 0.*N2
        

        ########################################################  

        if debug:

            w, divQ, terms_data = tools_g.solve_genomega(buoy,pm,pn,f,N2,depths,ur=ur,vr=vr,u=u,v=v,nh=nh,forcing = forcing,terms=terms,terms_data=terms_data,debug=True)  
            return w, divQ, terms_data

        else:
            
            w = tools_g.solve_genomega(buoy,pm,pn,f,N2,depths,ur=ur,vr=vr,u=u,v=v,nh=nh,forcing = forcing,terms=terms,smoothing=smoothing,nsmooth=nsmooth)   
            return w

        
        #######################################################        






















