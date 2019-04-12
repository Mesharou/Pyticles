



###################################################################################
# VARIABLES 
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
#import streamfunction as st

#Simulations (path, data...)
import simulations_old as oldsim

import R_smooth as sm

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
    {'temp': ['Temperature', 'C', [0,0,1]],\
    'salt': ['Salinity', 'PSU', [0,0,1]],\
    'u': ['u', 'm/s', [1,0,1]],\
    'v': ['v', 'm/s', [0,1,1]],\
    'ubar': ['ubar', 'm/s', [1,0,-1]],\
    'vbar': ['vbar', 'm/s', [0,1,-1]],\
    'zeta': ['SSH', 'm' ,[0,0,-1]],\
    'hbls': ['Thickness of KPP surface boundary layer', 'm', [0,0,-1]],\
    'hbbls': ['Thickness of KPP bottom boundary layer', 'm', [0,0,-1]],\
    'AKt': ['Temperature vertical diffusion coef', 'm2/s', [0,0,0]],\
    'AKv': ['Momentum vertical diffusion coef', 'm2/s', [0,0,0]],\
    'omega': ['S-coordinate vertical velocity', 'm/s ?', [0,0,0]],\
    \
    'psi': ['psi', 'Streamfunction' ,[1,1,-1]],\
    'int_psi': ['psi', 'Streamfunction' ,[1,1,-1]],\
    'psi_surf': ['psi_surf', 'Streamfunction' ,[1,1,-1]],\
    \
    'rho': ['in-situ density', 'kg.m-3', [0,0,1]],\
    'rho1': ['in-situ density', 'kg.m-3', [0,0,1]],\
    'rho1_sol2': ['in-situ density', 'kg.m-3', [0,0,1]],\
    'bvf': ['Brunt-Vaisala Frequency squared: N2', 's-2', [0,0,0]],\
    'buoy': ['buoyancy', 'm/s-2', [0,0,1]],\
    \
    'w': ['Vertical velocity', 'm/s', [0,0,1]],\
    'w_vmix': ['Vartical mixing contribution to vertical velocity', 'm/s', [0,0,1]],\
    \
    'absvrt': ['Absolute Vorticity', 's-1' ,[1,1,1]],\
    'vrt': ['Relative Vorticity', 's-1' ,[1,1,1]],\
    'pv': ['Potential Vorticity', 'PVU' ,[1,1,0]],\
    'pv1': ['[f + (dv/dx - du/dy)]*db/dz', 'PVU' ,[1,1,2]],\
    'pv2': ['-(dv/dz)*(db/dx)', 'PVU' ,[1,1,2]],\
    'pv3': ['(du/dz)*(db/dy)', 'PVU' ,[1,1,2]],\
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
    'MHdiff': ['Momentum implicit horizontal diff.', 'm.s-2' ,[0,0,1]],\
    'MVmix': ['Momentum vertical diffusion', 'm.s-2' ,[0,0,1]],\
    'MCor': ['Momentum Coriolis', 'm.s-2' ,[0,0,1]],\
    'Mrate': ['Momentum rate of change', 'm.s-2' ,[0,0,1]],\
    'MXadv': ['Momentum advection', 'm.s-2' ,[0,0,1]],\
    'MYadv': ['Momentum advection', 'm.s-2' ,[0,0,1]],\
    'MPrsgrd': ['Momentum pressure gradient', 'm.s-2' ,[0,0,1]],\
    \
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

    def __init__(self,varname,simul,n2max=100000,method='new',verbo=False,**kwargs):

        if varname!='zeta':
            print(' ')
            print('computing ', varname)
    

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
            
        '''
        3 different scenarii:
            1. var is already in output & we want surface or sigma-level (nothing to do)
            2. var is already in output BUT we want it on z-level (vert. interpolation needed)
            3. var is not in file (computations needed)
        '''

        ####################################################################
        if (self.name in list(ncfile.variables.keys())) and (method=='new'):
        #check if variable is already in output:  
        
            if verbo: print('self.imin, self.jmin ,self.kmin')
            if verbo: print(self.imin, self.jmin ,self.kmin) 
            if verbo: print(' ')
            
            if verbo: print('min(depths)' , min(depths))

            #Check type of vertical and horizontal grid
            if 'xi_u' in ncfile.variables[varname].dimensions: self.imin = 1
            if 'eta_v' in ncfile.variables[varname].dimensions: self.jmin = 1
            if 's_w' in ncfile.variables[varname].dimensions: self.kmin = 0

            #2D variable
            if 's_rho' not in ncfile.variables[varname].dimensions and 's_w' not in ncfile.variables[varname].dimensions:
                #this is a 2d variable -> just load it
                self.data = np.squeeze(simul.Forder(ncfile.variables[varname][simul.infiletime,ny1:ny2-self.jmin,nx1:nx2-self.imin]))


            #3D variable on one or more sigmal-levels
            elif min(depths)>=0:
                
                if verbo: print('3D variable on one or more sigmal-levels')
                if verbo: print(depths)
                
                if 's_w' in ncfile.variables[varname].dimensions:

                    if (len(depths)==1):   depth = depths[0] -1 
                    else:   depth = np.append(np.array(depths)-1,np.max(depths))
                    if verbo: print(depth)
                    self.data = np.squeeze(simul.Forder(ncfile.variables[varname][simul.infiletime,depth,ny1:ny2-self.jmin,nx1:nx2-self.imin]))

                elif 's_rho' in ncfile.variables[varname].dimensions:

                    if (len(depths)==1):   depth = depths[0] -1 
                    else:   depth = np.array(depths) - 1
                    self.data = np.squeeze(simul.Forder(ncfile.variables[varname][simul.infiletime,depth,ny1:ny2-self.jmin,nx1:nx2-self.imin]))
                    
                if verbo: print('var shape', self.data.shape)
                if verbo: print('simul.infiletime,depth,ny1,ny2-self.jmin,nx1,nx2-self.imin')  
                if verbo: print(simul.infiletime,depth,ny1,ny2-self.jmin,nx1,nx2-self.imin)     

            #3D variable on z-levels (interpolation needed)
            else:
                #if max(depths)>0: raise NameError('Depths are ill-defined. Check again please.')

                #Check how large is the domain___ Divide computations in chunk if $n_x*n_y > n2max$ 
                nchunk = int(np.max([np.sqrt((ny2-ny1)*(nx2-nx1)/n2max),1]))

                self.data = np.zeros((nx2-nx1-self.imin,ny2-ny1-self.jmin,np.max([len(depths),1])))*np.nan

                for i,j in product(list(range(nchunk)),list(range(nchunk))):

                    dx1=1; dx2=1; dy1=1; dy2=1;
                    if i==0: dx1=0
                    if i==nchunk-1: dx2=0 
                    if j==0: dy1=0 
                    if j==nchunk-1: dy2=0

                    nx1i,nx2i = nx1+i*(nx2-nx1)/nchunk-2*dx1,nx1+(i+1)*(nx2-nx1)/nchunk+2*dx2
                    ny1i,ny2i = ny1+j*(ny2-ny1)/nchunk-2*dy1,ny1+(j+1)*(ny2-ny1)/nchunk+2*dy2  

                    #we need to perform some vertical interpolation_ compute z_r,z_w for subdomain only
                    [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

                    #Compute variable in subdomain
                    chunk = tools.vinterp(simul.Forder(ncfile.variables[varname][simul.infiletime,:,ny1i:ny2i-self.jmin,nx1i:nx2i-self.imin]),\
                            depths,z_r,z_w,imin=self.imin,jmin=self.jmin,kmin=self.kmin,floattype = simul.floattype)

                    #Include in full variable                
                    self.data[nx1i+dx1-nx1:nx2i-dx2-nx1-self.imin,ny1i+dy1-ny1:ny2i-dy2-ny1-self.jmin,:] = \
                             chunk[dx1:nx2i-nx1i-dx2-self.imin,dy1:ny2i-ny1i-dy2-self.jmin,:]

                if len(depths)==1: self.data=self.data[:,:,0]


        ####################################################################
        elif (self.longname != 'unknown') and (method=='new'):

            if verbo: print('Variable not in ROMS outputs _ will be computed using new Fortran tools')  

            if (self.kmin>=0) and (min(depths)>0) and (len(depths)>1):
                self.data = np.zeros((nx2-nx1-self.imin,ny2-ny1-self.jmin,np.max([len(depths)+1-self.kmin,1])))*np.nan           
            elif  (self.kmin>=0) and (len(depths)>1):
                self.data = np.zeros((nx2-nx1-self.imin,ny2-ny1-self.jmin,np.max([len(depths),1])))*np.nan
            else:
                self.data = np.zeros((nx2-nx1-self.imin,ny2-ny1-self.jmin))*np.nan

            nchunk = int(np.max([np.sqrt((ny2-ny1)*(nx2-nx1)/n2max),1]))

            #You cannot compute psi in chunks
            if self.name in ['psi','psi_surf','int_psi']: nchunk=1

            if verbo: print('Domain will be divided in ', nchunk**2 , ' chunks')

            for i,j in product(list(range(nchunk)),list(range(nchunk))):

                dx1=2; dx2=2; dy1=2; dy2=2;
                if i==0: dx1=0
                if i==nchunk-1: dx2=0 
                if j==0: dy1=0 
                if j==nchunk-1: dy2=0

                nx1i,nx2i = nx1+i*(nx2-nx1)/nchunk-2*dx1,nx1+(i+1)*(nx2-nx1)/nchunk+2*dx2
                ny1i,ny2i = ny1+j*(ny2-ny1)/nchunk-2*dy1,ny1+(j+1)*(ny2-ny1)/nchunk+2*dy2  

                ####################################################################
                if (min(depths)==0) and (len(depths)==1) and (self.name not in ['rho1','pv']):
                    #problem here for variables computes using z_r (ex: rho1), False is temporary 
                    if i+j==0 and verbo: print("computing data at surface: ", depths)
                    
                    chunk = self.get(ncfile,simul,depths=[-1],coord=[ny1i,ny2i,nx1i,nx2i],subcoord=[ny1i-ny1,ny2i-ny1,nx1i-nx1,nx2i-nx1])
                    
                elif (min(depths)>0) or (self.name[:4]=='int_'):

                    if i+j==0 and verbo: print("computing data at sigma-level: ", depths)
                    
                    chunk = self.get(ncfile,simul,depths=depths,coord=[ny1i,ny2i,nx1i,nx2i],subcoord=[ny1i-ny1,ny2i-ny1,nx1i-nx1,nx2i-nx1])               

                else:
                    
                    #if i+j==0: print "computing data at z-levels: ", depths
                    #we need to perform some vertical interpolation_ compute z_r,z_w for subdomain only
                    [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])
                    
                    #Compute variable in subdomain
                    chunk = tools.vinterp(self.get(ncfile,simul,depths=simul.coordmax[4],coord=[ny1i,ny2i,nx1i,nx2i],subcoord=[ny1i-ny1,ny2i-ny1,nx1i-nx1,nx2i-nx1]),\
                            depths,z_r,z_w,imin=self.imin,jmin=self.jmin,kmin=self.kmin,floattype = simul.floattype)
                         
                if (self.kmin>=0) and (len(depths)>1):
                    self.data[nx1i+dx1-nx1:nx2i-dx2-nx1-self.imin,ny1i+dy1-ny1:ny2i-dy2-ny1-self.jmin,:] = \
                        chunk[dx1:nx2i-nx1i-dx2-self.imin,dy1:ny2i-ny1i-dy2-self.jmin,:]                    
                else:         
                
                    if len(chunk.shape)==3: chunk=chunk[:,:,0]
                
                    self.data[nx1i+dx1-nx1:nx2i-dx2-nx1-self.imin,ny1i+dy1-ny1:ny2i-dy2-ny1-self.jmin] = \
                         chunk[dx1:nx2i-nx1i-dx2-self.imin,dy1:ny2i-ny1i-dy2-self.jmin]

                #else:
                    #self.data[nx1i+dx1-nx1:nx2i-dx2-nx1-self.imin,ny1i+dy1-ny1:ny2i-dy2-ny1-self.jmin] = \
                         #chunk[dx1:nx2i-nx1i-dx2-self.imin,dy1:ny2i-ny1i-dy2-self.jmin]


            if len(depths)==1 and len(self.data.shape)>=3: self.data=self.data[:,:,0]

        ####################################################################

        else:

            if verbo: print('We will use older script version')  
            self.oldvar(varname,simul,depths = depths, u=u,v=v)
            #raise NameError('Sorry. I don t know how to compute '  + varname + '.')



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

    

###################################################################################
#Compute some variables 
###################################################################################



    def get(self,ncfile,simul,**kwargs):


        if 'coord' in  kwargs: 
            coord = kwargs['coord']
            [ny1i,ny2i,nx1i,nx2i] = coord[0:4]
            [ny1,ny2,nx1,nx2] = self.coord[0:4]
            pm = np.asfortranarray(simul.pm[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            pn = np.asfortranarray(simul.pn[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            f = np.asfortranarray(simul.f[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            mask = np.asfortranarray(simul.mask[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])
            topo = np.asfortranarray(simul.topo[nx1i-nx1:nx2i-nx1,ny1i-ny1:ny2i-ny1])            
        else: 
            coord = self.coord[0:4]
            pm = simul.pm
            pn = simul.pn
            f = simul.f
            mask = simul.mask
            topo = simul.topo
            
        if 'depths' in  kwargs: 
            depths = kwargs['depths']
        else: 
            depths = self.coord[4]
            
        #print "depths used in the get routine:", depths
        
        ################################################
        
        

        if self.name in ['rho']:

            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            T = self.load('temp',ncfile,simul,coord=coord, depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord, depths=depths)

            var = toolsF.rho_eos(T,S,z_r,z_w,simul.rho0)       
            ################################################
        
        

        elif self.name in ['rho1']:
            
            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])
            
            T = self.load('temp',ncfile,simul,coord=coord, depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord, depths=depths)
            
            var = toolsF.rho1_eos(T,S,z_r,simul.rho0)           

        ################################################

        elif self.name in ['bvf']:

            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            T = self.load('temp',ncfile,simul,coord=coord, depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord, depths=depths)

            var = toolsF.bvf_eos(T,S,z_r,z_w,simul.rho0)   
        
        
       ################################################


        elif self.name in ['buoy']:

            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            T = self.load('temp',ncfile,simul,coord=coord, depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord, depths=depths)

            var = toolsF.get_buoy(T,S,z_r,z_w,simul.rho0)

        ################################################
        

        elif self.name in ['psi']:

            u = self.load('ubar',ncfile,simul,coord=coord)
            v = self.load('vbar',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); #Hz[Hz==0] = np.nan
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
    
            u = self.load('u',ncfile,simul,coord=coord, depths=simul.coordmax[4])
            v = self.load('v',ncfile,simul,coord=coord, depths=simul.coordmax[4])
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            var = toolsF.get_wvlcty(u,v,z_r,z_w,pm,pn)

            var[0,:,:] = np.nan
            var[-1,:,:] = np.nan
            var[:,0,:] = np.nan
            var[:,-1,:] = np.nan
            
            print('depths', depths)
            print(var.shape, z_w.shape)
            
            if len(depths)==1:
                if depths==[-1]: 
                    var=var[:,:,-1] #should be 0 anyway


        ################################################


        
        elif self.name in ['w_vmix']:
            
            depths_r = simul.coordmax[4]
            depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            T = self.load('temp',ncfile,simul,coord=coord, depths=depths_r)
            S = self.load('salt',ncfile,simul,coord=coord, depths=depths_r)
            
            #buoy = toolsF.get_buoy(T,S,z_r,z_w,simul.rho0)
            buoy = toolsF.rho1_eos(T,S,z_r,simul.rho0)*(-simul.g/simul.rho0)
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
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            #this is the python version
            var = tools.get_vrt(u,v,z_r,z_w,pm,pn)
            
            
            
        ################################################
        elif self.name in ['pv']:  

            
            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            #this is the python version
            var = tools.PV(T,S,u,v,z_r,z_w,f,simul.g,simul.rho0,pm,pn)


               
            
        ################################################
        elif self.name in ['pv1','pv2','pv3']:  

            
            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            #this is the python version
            [pv1,pv2,pv3] = tools.PV_terms(T,S,u,v,z_r,z_w,f,simul.g,simul.rho0,pm,pn)       
            
            
            if self.name in ['pv1']: var = pv1[:,:,1:-1]
            elif self.name in ['pv2']: var = pv2[:,:,1:-1]
            elif self.name in ['pv3']: var = pv3[:,:,1:-1]

        
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

            print(depths)
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            vrt = tools.get_vrt(u,v,z_r,z_w,pm,pn)
            print('vrt.shape',vrt.shape) 
            del u,v
            
            Hz = tools.rho2psi(z_w[:,:,1:] - z_w[:,:,:-1])
            
            var = np.nansum(vrt*Hz,axis=2)
            
                        

        ################################################


        elif self.name in ['bpt']:

            T = self.load('temp',ncfile,simul,coord=coord)
            S = self.load('salt',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = tools.nanbnd(toolsF.get_bpt(T,S, z_r,z_w,simul.rho0,pm,pn) * tools.rho2psi(mask))

        ################################################


        elif self.name in ['bpts']:

            zeta = self.load('zeta',ncfile,simul,coord=coord)
            
            var = tools.jacob(zeta,topo,pm,pn)
            
        ################################################


        elif self.name in ['u_Prsgrd']:

            T = self.load('temp',ncfile,simul,coord=coord)
            S = self.load('salt',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = tools.nanbnd(toolsF.get_u_prsgrd(T,S, z_r,z_w,simul.rho0,pm,pn))

        ################################################

        elif self.name in ['vortplanet']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])
 
            var = -1*toolsF.get_vortplanet(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)

        ################################################

        elif self.name in ['vortstretch']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = -1*toolsF.get_vortstretch(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)
            
        ################################################

        elif self.name in ['vortstretch2']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = -1*toolsF.get_vortstretch2(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)            
            
        ################################################

        elif self.name in ['int_vortplanet']:
            '''
            int_vortplanet in vortplanet integrated from surface to depth
            
            '''
            print('integrating transrpot down to', depths[0])
            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            print('integrate between:', depths[0],depths[1])
            var = -1*toolsF.get_intvortplanet(u,v, z_r,z_w,pm,pn,f,depths[0],depths[1]) * tools.rho2psi(mask)

        ################################################

        elif self.name in ['vortstretch_sol1']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = toolsF.get_vortstretch_sol1(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)
            
         ################################################

        elif self.name in ['vortstretch_sol2']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = toolsF.get_vortstretch_sol2(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)
            
         ################################################

        elif self.name in ['vortstretch_sol3']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = toolsF.get_vortstretch_sol3(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)


   
            
        ################################################

        elif self.name in ['vortplantot_sol2']:

            #vortplantot_sol2 is Cor + adv dur to rotation of grid

            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = toolsF.get_vortplantot_sol2(u,v, z_r,z_w,pm,pn,f) * tools.rho2psi(mask)

            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan   

        ################################################

        elif self.name in ['vortplantot_sol3']:

            # vortplantot_sol3 is only Cor (should be same than vortplantot)
            # except computed with ROMS routines

            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = toolsF.get_vortplantot_sol3(u,v, z_r,z_w,pm,pn,f) * tools.rho2psi(mask)

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

            var = -1*toolsF.get_vortplantot(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)         


        ################################################


        elif self.name in ['vortadv']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,1:] - z_w[:,:,:-1]); Hz[Hz==0] = np.nan
            del z_r, z_w

            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

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

            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)


            var = toolsF.get_int_adv_sol2(u,v, z_r,z_w,pm,pn) * tools.rho2psi(mask)
            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan      
     
        ################################################
        
        elif self.name in ['vortadv_sol1']:

            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)
            depths  = self.coord[4]          
            depthsw = np.append(np.array(depths),np.max(depths)+1)
            w = self.load('omega',ncfile,simul,coord=coord,depths = depthsw)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)


            var = toolsF.get_adv_sol1(u,v,w, z_r,z_w,pm,pn) * tools.rho2psi(mask)
            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan


        ################################################

        elif self.name in ['vortadv_sol2']:

            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)


            var = toolsF.get_adv_sol2(u,v, z_r,z_w,pm,pn) * tools.rho2psi(mask)
            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan      
     
        ################################################
        #USELESS
        elif self.name in ['vortadv_sol3']:

            # vortadv_sol3 is vortadv_sol2 + adv due to rotation of grid


            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = toolsF.get_adv_sol2(u,v, z_r,z_w,pm,pn) * tools.rho2psi(mask)
            var = var + toolsF.get_uvgrid(u,v, z_r,z_w,pm,pn,f) * tools.rho2psi(mask)

            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan 
        
         ################################################

        elif self.name in ['vortadv_mix']:

            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)


            var = toolsF.get_adv_mix(u,v, z_r,z_w,pm,pn) * tools.rho2psi(mask)
            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan
            
        ################################################

        elif self.name in ['vortadv_centered']:

            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)


            var = toolsF.get_adv_4th(u,v, z_r,z_w,pm,pn) * tools.rho2psi(mask)
            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan      
            
         ################################################        

        #USELESS
        
        elif self.name in ['vortadv_uvgrid']:

            # vortadv_sol3 is vortadv_sol2 + adv due to rotation of grid


            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            #var = toolsF.get_adv_sol2(u,v, z_r,z_w,pm,pn) * tools.rho2psi(mask)
            var = toolsF.get_uvgrid(u,v, z_r,z_w,pm,pn,f) * tools.rho2psi(mask)

            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan 
            
        ###################################################################################

        elif self.name in ['rotwind','rotbot']:


            if self.name=='rotwind':
                
                [ut,vt] = self.get_winds(simul,coord=coord)
                
                ut=ut/simul.rho0; vt=vt/simul.rho0
                
                var = tools.rot(ut,vt,pm,pn) * tools.rho2psi(mask)     


            elif self.name=='rotbot':
                
                u = self.load('u',ncfile,simul,coord=coord,depths=[1])
                v = self.load('v',ncfile,simul,coord=coord,depths=[1])
                [z_r,z_w] = tools.get_depths(simul,coord=coord)
                Hz = z_w[:,:,1] - z_w[:,:,0]

                (ut,vt) = toolsF.get_bot(u,v,Hz,simul.rdrg)
                
                var = -1*tools.rot(ut,vt,pm,pn) * tools.rho2psi(mask)     

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

            T = self.load('temp',ncfile,simul,coord=coord)
            S = self.load('salt',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = tools.nanbnd(toolsF.get_bpt_mean(T,S, z_r,z_w,simul.rho0,pm,pn) * tools.rho2psi(mask))

            
        ################################################

        elif self.name in ['vorttopo_mean']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w
            
            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])
 
            var = -1*toolsF.get_vortplanet(ubar,vbar,Hz,pm,pn,1/Hz) * tools.rho2psi(mask) * tools.rho2psi(f)
            
        ################################################

        elif self.name in ['vortf_mean']:


            Hz = np.ones(pm.shape)
            
            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])
 
            var = -1*toolsF.get_vortplanet(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)
            
        ################################################

        elif self.name in ['vortplanet_mean']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w
            
            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])
 
            var = -1*toolsF.get_vortplanet(ubar,vbar,Hz,pm,pn,f/Hz) * tools.rho2psi(mask)

        ################################################

        elif self.name in ['vortstretch_mean']:

            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            Hz = (z_w[:,:,-1] - z_w[:,:,0]); Hz[Hz==0] = np.nan
            del z_r, z_w
            
            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = -1*toolsF.get_vortstretch(ubar,vbar,Hz,pm,pn,f/Hz) * tools.rho2psi(mask)
            


         ################################################

        elif self.name in ['vortstretch_sol2_mean']:

            Hz = np.ones(pm.shape)

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = toolsF.get_vortstretch_sol2(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)
            

 
            
        ################################################

        elif self.name in ['vortplantot_sol2_mean']:

            #vortplantot_sol2 is Cor + adv dur to rotation of grid

            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = toolsF.get_vortplantot_sol2_mean(u,v, z_r,z_w,pm,pn,f) * tools.rho2psi(mask)

            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan   



        ################################################

        elif self.name in ['vortplantot_mean']:

            Hz = np.ones(pm.shape)

            ubar = self.load('ubar',ncfile,simul,coord=coord,depths=[0])
            vbar = self.load('vbar',ncfile,simul,coord=coord,depths=[0])

            var = -1*toolsF.get_vortplantot(ubar,vbar,Hz,pm,pn,f) * tools.rho2psi(mask)         


        ################################################

        elif self.name in ['vortadv_sol2_mean']:
            ''' 
            In vortadv_sol2_mean we do the same computations than in 
            vortadv_sol2 but with a Hz(k)/H
               
            This is wrong and vortadv_sol3_mean should
            be used instead
            '''
               
            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)


            var = toolsF.get_adv_sol2_mean(u,v, z_r,z_w,pm,pn) * tools.rho2psi(mask)
            var[1,:] = np.nan
            var[-1,:] = np.nan
            var[:,1] = np.nan
            var[:,-1] = np.nan      

 
        ################################################

        elif self.name in ['vortadv_sol3_mean']:

            u = self.load('u',ncfile,simul,coord=coord)
            v = self.load('v',ncfile,simul,coord=coord)

            [z_r,z_w] = tools.get_depths(simul,coord=coord)

            var = toolsF.get_adv_sol3_mean(u,v, z_r,z_w,pm,pn) * tools.rho2psi(mask)
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

                (ut,vt) = toolsF.get_bot(u,v,Hz,simul.rdrg)
                
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
                var = tools.get_j2_sol1(T,S,u,v,z_r,z_w,simul.rho0,pm,pn,hbls)
            elif self.name=='J2_sol2':
                var = tools.get_j2_sol2(T,S,u,v,z_r,z_w,simul.rho0,pm,pn,hbls)

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
                hbbls = toolsF.get_hbbls_from_akt(AKt,z_w)
                print('evaluating hbbls using AKt')

            if self.name=='Jbot_sol1':
                var = tools.get_jbot_sol1(T,S,u,v,z_r,z_w,simul.rho0,pm,pn,hbbls,simul.rdrg)
            elif self.name=='Jbot_sol2':
                var = tools.get_jbot_sol2(T,S,u,v,z_r,z_w,simul.rho0,pm,pn,hbbls,simul.rdrg)

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
            [stflx, ssflx] = self.get_buoy_flux(simul,coord=coord)

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
        
        
        elif self.name in ['Madv','MXadv','MYadv','MVadv','MHdiff','MVmix','MCor','MPrsgrd']:

            [MXadv,MYadv,MVadv,MHdiff,MVmix,MCor,MPrsgrd] = self.get_tracer_evolution(simul,coord=coord)
            
            if self.name=='Madv':
              var = MXadv[:,:,:,:] + MYadv[:,:,:,:] + MVadv[:,:,:,:]  
            elif self.name=='MXadv':
              var = MXadv[:,:,:,:]
            elif self.name=='MYadv':
              var = MYadv[:,:,:,:]
            elif self.name=='MVadv':
              var = MVadv[:,:,:,:]                 
            elif self.name=='MHdiff':
              var = MHdiff[:,:,:,:]
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

                [ut,vt] = tools.get_bottom_drag(u,v,Hz,simul.rdrg)


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


    def oldvar(self,varname,simul,verbo=False,**kwargs):

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
        if verbo: print('[ny1,ny2,nx1,nx2,depths]',[ny1,ny2,nx1,nx2,depths])
        
        [ny0,nx0] = simul.coord[0],simul.coord[2]
        
        print('[ny1,ny2,nx1,nx2,depths]',[ny1,ny2,nx1,nx2])
        print('[ny1,ny2,nx1,nx2,depths]',[ny0,nx0])

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


#######################################################
#Get surface wind stress from forcing files interpolated to current time-step
#######################################################






    @staticmethod
    def get_winds(simul,**kwargs):


        if 'coord' in  kwargs: 
            [ny1,ny2,nx1,nx2]= kwargs['coord']
        else:
            [ny1,ny2,nx1,nx2] = simul.coord[0:4]


        # Load NETCDF files
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
        ncfilewind = Dataset(simul.ncname.wind, 'r', format='NETCDF3_CLASSIC')

        oceantime = int(np.array(ncfile.variables['ocean_time'][simul.infiletime]))%(360*24*3600)
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

        dsdt = 1./(90.*86400)
        ssflx = swflx*0.01/86400*salt - clim * dsdt*Hz*(salt - sss)
        
        
        #######################################################

        if 'DIURNAL_SRFLUX' in simul.cpp:
            
            sec2day=1./86400.
            tdays = simul.oceantime*sec2day
            cff=2.*simul.dt_model*sec2day

            phase=4.*(tdays-np.floor(tdays))-2.
            cff1=np.max([-1., np.min([1., np.min(phase-cff)])])
            cff2=np.max([-1., np.min([1., np.min(phase+cff)])])
            phase=(cff2-cff1)/cff + (np.sin(np.pi*cff2)-np.sin(np.pi*cff1))/(np.pi*cff)

            print('phase is ', phase)
            #######################################################

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
    b = g (alpha T - beta S)
        
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
#Compute thermal expansion and saline contraction coefficients
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
        
        u = self.load('u',ncfile,simul,coord=coord,depths=depths)
        if pert==1: u = u - np.mean(u); print('removing mean velocity')
        elif pert==0: u = u * np.nan + np.mean(u); print('using only mean velocity')
        v = self.load('v',ncfile,simul,coord=coord,depths=depths)
        if pert==1: v = v - np.mean(v)        
        elif pert==0: v = v * np.nan + np.mean(v); print('using only mean velocity')
        omega = self.load('omega',ncfile,simul,coord=coord,depths=depths_w)
        AKt = self.load('AKt',ncfile,simul,coord=coord,depths=depths_w)
        hbls = self.load('hbls',ncfile,simul,coord=coord)

        t=np.zeros((u.shape[0]+1,u.shape[1],u.shape[2],2))
        t[:,:,:,0] = self.load('temp',ncfile,simul,coord=coord,depths=depths)
        t[:,:,:,1] = self.load('salt',ncfile,simul,coord=coord,depths=depths)
        
        [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)
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
                toolsF.get_tracer_evolution (u,v, z_r,z_w,pm,pn\
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
        
        u = self.load('u',ncfile,simul,coord=coord,depths=depths)
        if pert==1: u = u - np.mean(u); print('removing mean velocity')
        elif pert==0: u = u * np.nan + np.mean(u); print('using only mean velocity')
        v = self.load('v',ncfile,simul,coord=coord,depths=depths)
        if pert==1: v = v - np.mean(v)        
        elif pert==0: v = v * np.nan + np.mean(v); print('using only mean velocity')
        omega = self.load('omega',ncfile,simul,coord=coord,depths=depths_w)
        AKt = self.load('AKt',ncfile,simul,coord=coord,depths=depths_w)
        hbls = self.load('hbls',ncfile,simul,coord=coord)

        t=np.zeros((u.shape[0]+1,u.shape[1],u.shape[2],2))
        t[:,:,:,0] = self.load('temp',ncfile,simul,coord=coord,depths=depths)
        t[:,:,:,1] = self.load('salt',ncfile,simul,coord=coord,depths=depths)
        
        [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)
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
                toolsF.get_tracer_evolution (u,v, z_r,z_w,pm,pn\
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
        
        if min(depths)<=0:
            
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

        
        elif min(depths)>0:     
        
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
        
        ud,vd = tools.div2uvs(u,v,pm,pn)
        ur,vr = u-ud, v-vd
        
        ########################################################      

        tenddiv = tools.get_tendency(ud,vd,buoy,pm,pn)        
        tendrot = tools.get_tendency(ur,vr,buoy,pm,pn)       

   
        ########################################################      
  
                                
        return [tenddiv[addx1:addx2,addy1:addy2,:],\
                tendrot[addx1:addx2,addy1:addy2,:]]
                    
           
           
           
           
                
#######################################################
#Comput
#######################################################



    def get_uv_evolution(self,simul,**kwargs):
    

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
        
        u = self.load('u',ncfile,simul,coord=coord,depths=depths)
        if pert==1: u = u - np.mean(u); print('removing mean velocity')
        elif pert==0: u = u * np.nan + np.mean(u); print('using only mean velocity')
        v = self.load('v',ncfile,simul,coord=coord,depths=depths)
        if pert==1: v = v - np.mean(v)        
        elif pert==0: v = v * np.nan + np.mean(v); print('using only mean velocity')
        
        omega = self.load('omega',ncfile,simul,coord=coord,depths=depths_w)
        AKv = self.load('AKv',ncfile,simul,coord=coord,depths=depths_w)
        #hbls = self.load('hbls',ncfile,simul,coord=coord)

        #t=np.zeros((u.shape[0]+1,u.shape[1],u.shape[2],2))
        T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
        S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
        
        [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)
        
        #Hz =  z_w[:,:,1:] - z_w[:,:,:-1]
        #stflx=np.zeros((u.shape[0]+1,u.shape[1],2))
        #[stflx[:,:,0], stflx[:,:,1],srflx] = self.get_buoy_flux(simul,coord=coord,solar=True)
        #[alpha,beta] = self.alphabeta(t[:,:,-1,0],t[:,:,-1,1],simul.rho0)
        
        [sustr,svstr] = self.get_winds(simul,coord=coord)
        sustr=sustr/simul.rho0; svstr=svstr/simul.rho0
    
        (MXadv,MYadv,MVadv,MHdiff,MVmix,MCor,MPrsgrd) = \
                                toolsF.get_uv_evolution (u,v,T,S, z_r,z_w,pm,pn,f\
                                                         ,simul.dt_model,simul.rdrg,simul.rho0,omega,AKv,sustr,svstr) 

        print('MXadv.shape', MXadv.shape)
        
        [MXadv,MYadv,MVadv,MHdiff,MVmix,MCor,MPrsgrd] = \
             [tools.nanbnd(MXadv,2),
             tools.nanbnd(MYadv,2),     
             tools.nanbnd(MVadv,2),
             tools.nanbnd(MHdiff,2),      
             tools.nanbnd(MVmix,2),
             tools.nanbnd(MCor,2),          
             tools.nanbnd(MPrsgrd,2)]        
             

        if len(MXadv.shape)==4:
            return [MXadv[addx1:addx2,addy1:addy2,:,:],MYadv[addx1:addx2,addy1:addy2,:,:],MVadv[addx1:addx2,addy1:addy2,:,:],\
                MHdiff[addx1:addx2,addy1:addy2,:,:],MVmix[addx1:addx2,addy1:addy2,:,:],MCor[addx1:addx2,addy1:addy2,:,:],\
                MPrsgrd [addx1:addx2,addy1:addy2,:,:]]
        elif len(MXadv.shape)==3:
            return [MXadv [addx1:addx2,addy1:addy2,:],MYadv [addx1:addx2,addy1:addy2,:],MVadv [addx1:addx2,addy1:addy2,:],\
                MHdiff [addx1:addx2,addy1:addy2,:],MVmix [addx1:addx2,addy1:addy2,:],MCor [addx1:addx2,addy1:addy2,:],\
                MPrsgrd [addx1:addx2,addy1:addy2,:]]



#######################################################
#Compute thermal expansion and saline contraction coefficients
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
            
            
        #######################################################

        depths = simul.coordmax[4]
        depths_w = np.concatenate((simul.coordmax[4],[simul.coordmax[4][-1]+1]))

        #######################################################
 
        ncfile = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
            
        #######################################################
        
        u = self.load('u',ncfile,simul,coord=coord,depths=depths)
        v = self.load('v',ncfile,simul,coord=coord,depths=depths)

        [sustr,svstr] = self.get_winds(simul,coord=coord)
        sustr=sustr/simul.rho0; svstr=svstr/simul.rho0    

        print('u.shape, v.shape', u.shape, v.shape)
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

        stflx=np.zeros((u.shape[0]+1,u.shape[1],2))
        
        [stflx[:,:,0], stflx[:,:,1],srflx] = self.get_buoy_flux(simul,coord=coord,solar=True)
        
        #######################################################
        if simul.simul in ['atlbig','gulfz','nesea','neseb','nefro']:
            Ricr=0.15 # for SHING family      
        else:
            Ricr=0.45 # for SHING family
    
        (hbl, out1, out2, out3, out4) = toolsF.get_hbl (alpha, beta, z_r, z_w, stflx, srflx, swr_frac, sustr, svstr, Ricr, hbls, f,\
                                u, v, bvf)
    
             
        #######################################################

        return hbl,hbls,out1,out2,out3,out4

        

        
        
        
#######################################################
#Compute vertical velocity due to mixing (in the GaLo81 way)
#######################################################

    def get_w_vmix(self,simul,**kwargs):
        '''
        Solving the equation 

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
        buoy = tools.vinterp(-simul.g * toolsF.rho_eos(T,S,z_r,z_w,simul.rho0)/simul.rho0,depths,z_r,z_w)
        del T,S

        print('depths_w', depths_w)
        AKv = self.load('AKv',ncfile,simul,coord=coord, depths=depths_w)
        print('AKv.shape', AKv.shape)
        
        #AKv = tools.vinterp(AKv,depths,z_w)     
        AKv = tools.vinterp(0.5*(AKv[:,:,1:]+AKv[:,:,:-1]),depths,z_r,z_w)     
        
        ################################################dsdsdsd
        thermal = True
        
        if thermal:
            
            uz = -1*(tools.diffy(buoy,pn))
            vz = (tools.diffx(buoy,pm))
            
            byy = tools.diffy(tools.rho2v(AKv) * uz,tools.rho2v(pn))
            bxx = -1*tools.diffx(tools.rho2u(AKv) * vz,tools.rho2u(pm))     
                  
            var = np.zeros(buoy.shape)*np.nan
            var[1:-1,:,:] = bxx 
            var[:,1:-1,:] = var[:,1:-1,:] + byy        
            var = var[:,:,0]    
            var = (var.T/simul.f.T).T
            
        else:
            
            u = self.load('u',ncfile,simul,coord=coord,depths=depths_r)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths_r)

            uz = (u[:,:,1:] - u[:,:,:-1])/(tools.rho2u(z_r[:,:,1:])-tools.rho2u(z_r[:,:,:-1]))
            vz = (v[:,:,1:] - v[:,:,:-1])/(tools.rho2v(z_r[:,:,1:])-tools.rho2v(z_r[:,:,:-1]))

            #uz[:,:,1:-1] = (u[:,:,2:] - u[:,:,:-2])/(2*dz)
            #uz[:,:,1] = (u[:,:,1] - u[:,:,0])/(dz)*0; uz[:,:,-1] = (u[:,:,-1] - u[:,:,-2])/(dz)*0
            #vz[:,:,1:-1] = (v[:,:,2:] - v[:,:,:-2])/(2*dz)
            #vz[:,:,1] = (v[:,:,1] - v[:,:,0])/(dz)*0; vz[:,:,-1] = (v[:,:,-1] - v[:,:,-2])/(dz)*0
            
            print(uz.shape)
            print(tools.rho2u(z_w[:,:,1:-1]).shape)
            
            uz = tools.vinterp(uz,depths,tools.rho2u(z_w[:,:,1:-1]),tools.rho2u(z_r)) 
            vz = tools.vinterp(vz,depths,tools.rho2v(z_w[:,:,1:-1]),tools.rho2v(z_r) )
            
            
            byy = tools.diffy(tools.rho2u(AKv) * uz,tools.rho2u(pn))
            bxx = -1*tools.diffx(tools.rho2v(AKv) * vz,tools.rho2v(pm))
            
            #var = np.zeros((u.shape[0],v.shape[1],buoy.shape[2])*np.nan
            var = tools.psi2rho(bxx + byy)           
            var = var[:,:,0]
            
        ################################################        
        #add winds
        
        [sustr,svstr] = self.get_winds(simul,coord=coord)
        
        print(sustr.shape, svstr.shape)
        var = var - tools.psi2rho((tools.diffy(sustr,tools.rho2u(pn))- tools.diffx(svstr,tools.rho2v(pm)) ))/simul.rho0

        var = (var.T/simul.f.T).T
                
                
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
        #u,v = tools.solve_ttwall(bx,by,AKv,sustr,svstr,f,pm,pn,depths)
        
        if timing: print('computation all.........', tm.time()-tstart)
        if timing: tstart = tm.time()   
        
        u,v = tools.solve_ttw(bx,by,AKv,sustr,svstr,f,pm,pn,depths)
        if timing: print('computation old.........', tm.time()-tstart)
        
        
        u = u*mask_depth
        v = v*mask_depth        
        ################################################     
        
        return u,v
        
        
        
#######################################################
#Compute vertical velocity due to TTW
#######################################################

    def get_w_ttw(self,simul,**kwargs):
        
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
        
        u,v = self.get_uv_ttw(simul,depths=depths)
        
        w = -1*tools.nanbnd(np.cumsum(tools.u2rho(tools.diffx(u_ttw,simul.pm)) + tools.v2rho(tools.diffy(v_ttw,simul.pn)),2)*(depths[1] - depths[0]),1)
        
        ################################################     
        
        return w        
   
   
   
        
#######################################################
#Compute horizomtal velocities due to TTW
#######################################################

    def get_uv_ttw_sig(self,simul,debug=0,**kwargs):
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
        

        bx= tools.u2rho(tools.diffxi(buoy,pm,z_r,z_w,z_w[:,:,:])); 
        by = tools.v2rho(tools.diffeta(buoy,pn,z_r,z_w,z_w[:,:,:]))
            #bx[1:,:,:]= tools.diffx(buoy,pm); by[:,1:,:] = tools.diffy(buoy,pn)


        del buoy
        #xnan = np.isnan(bx)
        #bx[xnan]=0
        #bx[np.isnan(by)]=0
        #by[xnan]=0
        #by[np.isnan(by)]=0
        
        ################################################        
        zeta = self.load('zeta',ncfile,simul,coord=coord)      
        AKv = self.load('AKv',ncfile,simul,coord=coord, depths=depths_w)      
        AKv[:,:,-1] = tools.vinterp(AKv[:,:,:-1],zeta,z_w[:,:,:-1])[:,:,-1]; 

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
        
        if timing: tstart = tm.time()   
        #u,v = tools.solve_ttwall(bx,by,AKv,sustr,svstr,f,pm,pn,depths)
        
        if timing: print('computation all.........', tm.time()-tstart)
        if timing: tstart = tm.time()   
        
        if debug==1:        
            u,v,ug,vg,uz,vz = tools.solve_ttw_sig(bx,by,AKv,sustr,svstr,f,pm,pn,z_w,debug=1)
        else:
            u,v,ug,vg = tools.solve_ttw_sig(bx,by,AKv,sustr,svstr,f,pm,pn,z_w)          
        if timing: print('computation old.........', tm.time()-tstart)
        
        # Return result on vertical rho-grid
        u = tools.vinterp(u,z_r,z_r,z_w,kmin=0)
        v = tools.vinterp(v,z_r,z_r,z_w,kmin=0)
        ug = tools.vinterp(ug,z_r,z_r,z_w,kmin=0)
        vg = tools.vinterp(vg,z_r,z_r,z_w,kmin=0)
     
        ################################################     
        if debug==1:
            return u,v,ug,vg,uz,vz
        else:  
            return u,v,ug,vg
        
        
#######################################################
#Compute vertical velocity due to TTW
#######################################################

    def get_w_ttw_sig(self,simul,debug=0,**kwargs):
        
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
        
        if debug==1:
            ut,vt,ug,vg,uz,vz = self.get_uv_ttw_sig(simul,debug=1)
        else:
            ut,vt,ug,vg = self.get_uv_ttw_sig(simul)

        [z_r,z_w] = tools.get_depths(simul,coord=coord)
        w = toolsF.get_wvlcty(tools.rho2u(ut-ug),tools.rho2v(vt-vg),z_r,z_w,pm,pn)
        
        ################################################     
        
        if debug==1:
            return w,uz,vz
        else:  
            return w

        
        
        
#######################################################
#Compute QG omega equation
#######################################################


    def get_w_omega(self,simul,mixrotuv=True,u_ttw=None,v_ttw=None,**kwargs):

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
            
        mixing = False

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
        #rho = toolsF.rho1_eos(T,S,z_r,simul.rho0)

        ########################################################       
        # add mixing tracer term: 
        ########################################################       

        if mixing:
            
            AKt = self.load('AKt',ncfile,simul,coord=coord, depths=depths_w)
            buoy = -simul.g *rho/simul.rho0
            
            dbdz = np.zeros(AKt.shape)
            dbdz[:,:,1:-1] = (buoy[:,:,1:] - buoy[:,:,:-1]) / (z_r[:,:,1:] - z_r[:,:,:-1] )
            
            #no flux for now: we should add surface forcings!!!
            dbdz[:,:,0] = dbdz[:,:,1]
            dbdz[:,:,-1] =dbdz[:,:,-2]
            
            dbdz2 = (AKt[:,:,1:]*dbdz[:,:,1:] - AKt[:,:,:-1]*dbdz[:,:,:-1]) / (z_w[:,:,1:] - z_w[:,:,:-1] )
            
            TVmix = copy(dbdz2)
            TVmix[1:-1,:,:] = tools.diffx(tools.diffx(dbdz2,pm),tools.rho2u(pm))
            TVmix[:,1:-1,:] = TVmix[:,1:-1,:] + tools.diffy(tools.diffy(dbdz2,pn),tools.rho2v(pn))

            forcing = tools.vinterp(TVmix,depths,z_r,z_w)
            del TVmix, dbdz2, dbdz
            
        else:
            
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
        #N2 = bvf
        del bvf
        
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
            ud,vd = tools.div2uvs(u,v,pm,pn)
            u,v = u-ud, v-vd
            del ud,vd
            
        #######################################################        
        # We use only the  divergent field
        if field in ['udiv']:
            print('using divergent velocity field')
            u,v = tools.div2uvs(u,v,pm,pn)      

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

        
        if 'nh' in  kwargs: nh=kwargs['nh']
        
        if field in ['urot','utot','udiv']: 
            u,v = tools.u2rho(u)* mask_depth, tools.v2rho(v)* mask_depth
            w = tools.solve_omega(buoy,pm,pn,f,N2,depths,u,v,nh=nh,forcing = forcing,mixrotuv=mixrotuv)
        elif field in ['uttw']:
            w = tools.solve_omega(buoy,pm,pn,f,N2,depths,u_ttw,v_ttw,nh=nh,forcing = forcing,mixrotuv=mixrotuv)            
        else:
            w = tools.solve_omega(buoy,pm,pn,f,N2,depths,nh=nh,forcing = forcing)

        #######################################################        


        return w






























