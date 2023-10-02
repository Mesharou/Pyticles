



###################################################################################
# VARIABLES 
###################################################################################

"""
16/06/08: Corrected bug in checking for fillvalue in the netcdf file leading to extensive memory usage
16/01/20: Modif for periodic files of JC [add parameter iper,jper], because xi-u has the same size than xi-rho
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
    {'temp': ['Temperature', r'$T\,(^{\circ}C)$', [0,0,1]],\
    'salt': ['Salinity', 'PSU', [0,0,1]],\
    'u': ['u', 'm/s', [1,0,1]],\
    'v': ['v', 'm/s', [0,1,1]],\
    'ubar': ['ubar', 'm/s', [1,0,-1]],\
    'vbar': ['vbar', 'm/s', [0,1,-1]],\
    'zeta': ['SSH', r'$\eta\,(m)$' ,[0,0,-1]],\
    'hbls': ['Thickness of KPP surface boundary layer', 'm', [0,0,-1]],\
    'hbbls': ['Thickness of KPP bottom boundary layer', 'm', [0,0,-1]],\
    'hbls_rho': ['Surface mixed-layer (based on rho = rhos+ 0.03)', 'm', [0,0,-1]],\
    'hbls_t': ['Surface mixed-layer (based on t = ts - 0.2)', 'm', [0,0,-1]],\
    'AKt': ['Temperature vertical diffusion coef', 'm2/s', [0,0,0]],\
    'AKv': ['Momentum vertical diffusion coef', 'm2/s', [0,0,0]],\
    'omega': ['S-coordinate vertical velocity', 'm/s ?', [0,0,0]],\
    \
    'psi': ['psi', 'Streamfunction for depth integrated flow' ,[1,1,-1]],\
    'psi_surf': ['psi_surf', 'Streamfunction of surface flow' ,[1,1,-1]],\
    \
    'rho': ['in-situ density', 'kg.m-3', [0,0,1]],\
    'rho1': ['in-situ density standard pressure', 'kg.m-3', [0,0,1]],\
    'rhop': ['potential density', 'kg.m-3', [0,0,1]],\
    'bvf': ['Brunt-Vaisala Frequency squared: N2', 's-2', [0,0,0]],\
    'buoy': ['buoyancy', 'm/s-2', [0,0,1]],\
    \
    'w': ['Vertical velocity', r'$w\,(m\,s^{-1})$', [0,0,1]],\
    \
    'absvrt': ['Absolute Vorticity', 's-1' ,[1,1,1]],\
    'vrt': ['Relative Vorticity', r'$\frac{\zeta}{f}$' ,[1,1,1]],\
    'pv': ['Potential Vorticity', 'PVU' ,[1,1,0]],\
    'pvr': ['Potential Vorticity on rho levels', 'PVU' ,[1,1,1]],\
    }



###################################################################################
#Load variables
###################################################################################

    def __init__(self,varname,simul,n2max=100000,method='new',verbo=False,**kwargs):

    


        #print '#######################################'
        #print 'load var', varname
        #print '#######################################'


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
            if ncfile.variables['xi_u'].shape[0]==ncfile.variables['xi_rho'].shape[0]:
                xperiodic = True; iper=1
        except:
            pass

        try:
            if ncfile.variables['eta_v'].shape[0]==ncfile.variables['eta_rho'].shape[0]:
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
        if (self.name in list(ncfile.variables.keys())) and (self.name not in ['w']):
        # check if variable is already in output:  (add test for w because problem with w in avg files (if produced before 05/05/14))
        # test for ['rho'] removed 17/11/24

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


            if verbo: print('self.imin, self.jmin ,self.kmin')
            if verbo: print(self.imin, self.jmin ,self.kmin) 
            if verbo: print(' ')
            
            if verbo: print('min(depths)' , min(depths))

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
            elif min(depths)>=0:
                
                if verbo: print('3D variable on one or more sigmal-levels')
                if verbo: print(depths)
                
                if 's_w' in ncfile.variables[varname].dimensions:

                    if (len(depths)==1):   depth = depths[0] -1 
                    else:   depth = np.append(np.array(depths)-1,np.max(depths))
                    if verbo: print(depth)
                    self.data = np.squeeze(simul.Forder(ncfile.variables[varname][simul.infiletime,depth,ny1+self.jper:ny2-self.jmin+self.jper,nx1+self.iper:nx2-self.imin+self.iper]))

                elif 's_rho' in ncfile.variables[varname].dimensions:

                    if (len(depths)==1):   depth = depths[0] -1 
                    else:   depth = np.array(depths) - 1
                    self.data = np.squeeze(simul.Forder(ncfile.variables[varname][simul.infiletime,depth,ny1+self.jper:ny2-self.jmin+self.jper,nx1+self.iper:nx2-self.imin+self.iper]))
                    
                if verbo: print('var shape', self.data.shape)
                if verbo: print('simul.infiletime,depth,ny1,ny2-self.jmin,nx1,nx2-self.imin')  
                if verbo: print(simul.infiletime,depth,ny1,ny2-self.jmin,nx1,nx2-self.imin)     
                
            ####################################################################
            #1.3 It is a 3D variable on z-levels (interpolation needed)
            ####################################################################
            else:
                #if max(depths)>0: raise NameError('Depths are ill-defined. Check again please.')

                #Check how large is the domain___ Divide computations in chunk if $n_x*n_y > n2max$ 
                nchunk = int(np.max([np.sqrt((ny2-ny1)*(nx2-nx1)//n2max),1]))

                self.data = np.zeros((nx2-nx1-self.imin,ny2-ny1-self.jmin,np.max([len(depths),1])))*np.nan

                for i,j in product(list(range(nchunk)),list(range(nchunk))):

                    dx1=4; dx2=4; dy1=4; dy2=4;
                    if i==0: dx1=0
                    if i==nchunk-1: dx2=0 
                    if j==0: dy1=0 
                    if j==nchunk-1: dy2=0

                    nx1i,nx2i = nx1+i*(nx2-nx1)//nchunk-2*dx1,nx1+(i+1)*(nx2-nx1)//nchunk+2*dx2
                    ny1i,ny2i = ny1+j*(ny2-ny1)//nchunk-2*dy1,ny1+(j+1)*(ny2-ny1)//nchunk+2*dy2  

                    #we need to perform some vertical interpolation_ compute z_r,z_w for subdomain only
                    [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

                    #Compute variable in subdomain
                    chunk = tools.vinterp(simul.Forder(ncfile.variables[varname][simul.infiletime,:,ny1i+self.jper:ny2i-self.jmin+self.jper,nx1i+self.iper:nx2i-self.imin+self.iper]),\
                            depths,z_r,z_w,imin=self.imin,jmin=self.jmin,kmin=self.kmin,floattype = simul.floattype)

                    #Include in full variable                
                    self.data[nx1i+dx1-nx1:nx2i-dx2-nx1-self.imin,ny1i+dy1-ny1:ny2i-dy2-ny1-self.jmin,:] = \
                             chunk[dx1:nx2i-nx1i-dx2-self.imin,dy1:ny2i-ny1i-dy2-self.jmin,:]

                if len(depths)==1: self.data=self.data[:,:,0]

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
        elif (self.longname != 'unknown'):

            if verbo: print('Variable not in ROMS outputs _ will be computed using new Fortran tools')  

            ####################################################################
            # Create an array to store the variable (self.data)
            ####################################################################

            if (self.kmin>=0) and (min(depths)>0) and (len(depths)>1):
                self.data = np.zeros((nx2-nx1-self.imin,ny2-ny1-self.jmin,np.max([len(depths)+1-self.kmin,1])))*np.nan           
            elif  (self.kmin>=0) and (len(depths)>1):
                self.data = np.zeros((nx2-nx1-self.imin,ny2-ny1-self.jmin,np.max([len(depths),1])))*np.nan
            else:
                self.data = np.zeros((nx2-nx1-self.imin,ny2-ny1-self.jmin))*np.nan

            # Number of tiles used for computation
            nchunk = int(np.max([np.sqrt((ny2-ny1)*(nx2-nx1)//n2max),1]))

            #You cannot compute psi in chunks: 
            if self.name in ['psi','psi_surf']: nchunk=1

            if verbo: print('Domain will be divided in ', nchunk**2 , ' chunks')

            for i,j in product(list(range(nchunk)),list(range(nchunk))):
                
                dx1=2; dx2=2; dy1=2; dy2=2; # extra pts around the tile
                if i==0: dx1=0
                if i==nchunk-1: dx2=0 
                if j==0: dy1=0 
                if j==nchunk-1: dy2=0

                nx1i,nx2i = nx1+i*(nx2-nx1)//nchunk-2*dx1,nx1+(i+1)*(nx2-nx1)//nchunk+2*dx2
                ny1i,ny2i = ny1+j*(ny2-ny1)//nchunk-2*dy1,ny1+(j+1)*(ny2-ny1)//nchunk+2*dy2  

                ####################################################################
                # Compute variable
                ####################################################################         
                print('depths',depths,min(depths),self.name)

                if (min(depths)>0) or (self.kmin<0):
                    
                    chunk = self.get_sig(ncfile,simul,depths=depths,coord=[ny1i,ny2i,nx1i,nx2i],subcoord=[ny1i-ny1,ny2i-ny1,nx1i-nx1,nx2i-nx1])


                elif (len(depths)==1) and (min(depths)==0) and (self.name in ['vrt','absvrt']):
                    
                    print('computing vrt')
                    chunk = self.get_sig(ncfile,simul,depths=depths,coord=[ny1i,ny2i,nx1i,nx2i],subcoord=[ny1i-ny1,ny2i-ny1,nx1i-nx1,nx2i-nx1])
                    
                    
                else:
                    [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])
                    
                    #Compute variable in subdomain
                    chunk = tools.vinterp(self.get_sig(ncfile,simul,depths=simul.coordmax[4],coord=[ny1i,ny2i,nx1i,nx2i],subcoord=[ny1i-ny1,ny2i-ny1,nx1i-nx1,nx2i-nx1]),\
                            depths,z_r,z_w,imin=self.imin,jmin=self.jmin,kmin=self.kmin,floattype = simul.floattype)

                
                ####################################################################
                # Write the chunk into the self.data
                ####################################################################


                if (self.kmin>=0) and (len(depths)>1):
                    # 3D variable
                    self.data[nx1i+dx1-nx1:nx2i-dx2-nx1-self.imin,ny1i+dy1-ny1:ny2i-dy2-ny1-self.jmin,:] = \
                        chunk[dx1:nx2i-nx1i-dx2-self.imin,dy1:ny2i-ny1i-dy2-self.jmin,:]                   
                else:
                    # 2D variable            
                    self.data[nx1i+dx1-nx1:nx2i-dx2-nx1-self.imin,ny1i+dy1-ny1:ny2i-dy2-ny1-self.jmin] = \
                         np.squeeze(chunk)[dx1:nx2i-nx1i-dx2-self.imin,dy1:ny2i-ny1i-dy2-self.jmin]


            if len(depths)==1 and len(self.data.shape)>=3: self.data=self.data[:,:,0]

        ####################################################################

        else:
            
            #if verbo: print 'We will use older script version'/oldvar
            #self.oldvar(varname,simul,depths = depths, u=u,v=v)
            
            raise NameError('Sorry. I don t know how to compute '  + varname + '.')



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

        if 'masked' in  kwargs: 
            masked_value =  kwargs['masked']
        else: 
            masked_value =  np.nan


            
        ############################# 
        # test if periodic
        #############################
        xperiodic = False; iper=0*imin
        yperiodic = False; jper=0*jmin

        try:
            if ncfile.variables['xi_u'].shape[0]==ncfile.variables['xi_rho'].shape[0]:
                xperiodic = True; iper=1*imin
        except:
            pass

        try:
            if ncfile.variables['eta_v'].shape[0]==ncfile.variables['eta_rho'].shape[0]:
                yperiodic = True; jper=1*jmin
        except:
            pass

        #############################

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


    def get_sig(self,ncfile,simul,**kwargs):

        verbo = False
        
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


        ################################################


        if self.name in ['rho','rho1','rhop','bvf','buoy','buoy1']:

            [z_r,z_w] = tools.get_depths(simul,coord=[ny1i,ny2i,nx1i,nx2i])

            T = self.load('temp',ncfile,simul,coord=coord, depths=depths)
            try:
                S = self.load('salt',ncfile,simul,coord=coord, depths=depths)
            except:
                print('no S in file')
                S = T*0.
              
            if self.name in ['rho']: var = toolsF.rho_eos(T,S,z_r,z_w,simul.rho0)
            elif self.name in ['rho1']: var = toolsF.rho1_eos(T,S,z_r,simul.rho0)   
            elif self.name in ['rhop']: var = tools.rhop(T,S)
            elif self.name in ['bvf']: var = toolsF.bvf_eos(T,S,z_r,z_w,simul.rho0)   
            elif self.name in ['buoy']: var = toolsF.get_buoy(T,S,z_r,z_w,simul.rho0)
            elif self.name in ['buoy1']: var = toolsF.rho1_eos(T,S,z_r,simul.rho0)*(-simul.g/simul.rho0)


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
                    var[i,j] = np.min([-z_r[i,j,np.min([1+np.nanargmax(np.where((rho[i,j,:] - rho[i,j,-1])>0.03)),rho.shape[2]-1])],topo[i,j]])
       
       #######################################################

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


        elif self.name in ['w']:
    
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)
            
            [z_r,z_w] = tools.get_depths(simul,coord=coord)
            var = toolsF.get_wvlcty(u,v,z_r,z_w,pm,pn)

            #var= tools.nanbnd(var)


                    
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

            [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)
            
            print('for PV using levels', depths)

            #this is the python version
            var = tools.PV(T,S,u,v,z_r,z_w,f,simul.g,simul.rho0,pm,pn)

        ################################################
        elif self.name in ['pvr']:


            T = self.load('temp',ncfile,simul,coord=coord,depths=depths)
            S = self.load('salt',ncfile,simul,coord=coord,depths=depths)
            u = self.load('u',ncfile,simul,coord=coord,depths=depths)
            v = self.load('v',ncfile,simul,coord=coord,depths=depths)

            [z_r,z_w] = tools.get_depths(simul,coord=coord,depths=depths)

            print('for PV using levels', depths)

            #this is the python version
            var = tools.PVr(T,S,u,v,z_r,z_w,f,simul.g,simul.rho0,pm,pn)




        ################################################
        
        return var



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
        ncfile = Dataset(simul.ncfile, 'r')

        try:
            print('check if sustr in ', simul.ncfile)
            #uwind = self.load('sustr',ncfile,simul,coord=[ny1,ny2,nx1,nx2])
            #vwind = self.load('svstr',ncfile,simul,coord=[ny1,ny2,nx1,nx2])
            uwind = simul.Forder( np.array(ncfile.variables['sustr'][simul.infiletime,ny1:ny2,nx1:nx2-1]) )
            vwind = simul.Forder( np.array(ncfile.variables['svstr'][simul.infiletime,ny1:ny2-1,nx1:nx2]) )
            print('loading sustr,svstr from his/avg file')
        
        except:
            print('computing sustr,svstr from frc file')
            ncfilewind = Dataset(simul.ncname.wind, 'r')

            try:
                oceantime = int(np.array(ncfile.variables['ocean_time'][simul.infiletime]))%(360*24*3600)
            except:
                oceantime = int(np.array(ncfile.variables['scrum_time'][simul.infiletime]))%(360*24*3600)

            oceanhour=oceantime//(3600.)
            oceanday=oceantime//(24*3600.)


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
        oceanday=oceantime//(24*3600.)


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
#Compute correlation between variable and observations
#######################################################


    def correlation_dong(variable,obs):
                 
        return 0.95
    
    #######################################################


    def correlation_dong_optimale(variable,obs):
                 
        return 0.97
    
        








