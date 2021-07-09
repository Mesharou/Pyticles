
"""

"""
from __future__ import print_function



###################################################################################
#Load some common modules
###################################################################################

#from builtins import object
import sys

#for netcdf files
#from Scientific.IO.NetCDF import *
from netCDF4 import Dataset

#for numeric functions
import numpy as np

import os


###################################################################################


class ionetcdf(object):

###################################################################################
#   Main 
###################################################################################

    def __init__(self,newfile,simul,var,nctime=None,name='var',shape = [0,0,1],coord=[0,None,0,None,0,None],netcdf_format='NETCDF4_CLASSIC',zlib=False,static=False,**kwargs):

        """

        """
        
        if 'dims' in kwargs: dims=kwargs['dims']
        else: dims = [simul.x.shape[0], simul.x.shape[1], len(simul.coord[4])]
        
        if not os.path.isfile(newfile): newnc = self.create(newfile,simul,netcdf_format=netcdf_format,dims=dims)
        #else:  os.remove(newfile); newnc = self.create(newfile,simul)
        else:  newnc = Dataset(newfile, 'a', format=netcdf_format)

        #[dim0, dim3, dim2, dim1]= ['time', ['depth_w','depth','depth_w_old'], ['lat','lat_v'], ['lon','lon_u']]
        [dim0, dim3, dim2, dim1]= ['time', ['s_w','s_rho','s_w_old'], ['eta_rho','eta_v'], ['xi_rho','xi_u']]
        
        
        # If you want to write only part of the array
        [ix1,ix2,iy1,iy2,iz1,iz2] = coord
        
        try:
            name = var.name
            data = var.data
            [imin,jmin,kmin] = [var.imin,var.jmin,var.kmin]
        except (AttributeError, TypeError):
            data = var
            [imin,jmin,kmin] = shape

        if nctime==None:  nctime = len(newnc.dimensions['time'])

        if isinstance(data,int) or isinstance(data,float):

            if name not in list(newnc.variables.keys()):
                newnc.createVariable(name, 'f', (dim0,) )
                
            newnc.variables[name][nctime]=data
            
        else:

            print('                           ')
            if name not in list(newnc.variables.keys()):
                print('name not in newnc.variables.keys()')
                if netcdf_format=='NETCDF3_CLASSIC':
                    print('netcdf_format==NETCDF3_CLASSIC')
                    if static:
                        newnc.createVariable(name, 'f', (dim2[jmin],dim1[imin],) )
                    else:
                        if len(data.shape)>2:
                            print('len(data.shape)>2')
                            newnc.createVariable(name, 'f', (dim0,dim3[kmin],dim2[jmin],dim1[imin],))
                        else:
                            print('len(data.shape)=2')
                            newnc.createVariable(name, 'f', (dim0,dim2[jmin],dim1[imin],) )
                            print('ok')
                else:  
                    print('netcdf_format==NETCDF4')
                    if static:
                        newnc.createVariable(name, 'f', (dim2[jmin],dim1[imin],) ,zlib=zlib)
                    else:
                        if len(data.shape)>2:
                            print('len(data.shape)>2')
                            #newnc.createVariable(name, 'f4', (dim0,dim3[kmin],dim2[jmin],dim1[imin],) ,zlib=zlib,chunksizes = (1,1,700,600))
                            newnc.createVariable(name, 'f4', (dim0,dim3[kmin],dim2[jmin],dim1[imin],) ,zlib=zlib)
                        else:
                            print('len(data.shape)=2')
                            #newnc.createVariable(name, 'f4', (dim0,dim2[jmin],dim1[imin],) ,zlib=zlib,chunksizes = (1,700,600))
                            newnc.createVariable(name, 'f4', (dim0,dim2[jmin],dim1[imin],) ,zlib=zlib)
                        
            if len(data.shape)>2:
                 #print 'newnc.variables[name][:].shape',newnc.variables[name][:].shape
                 #print 'data.shape',data.shape
                 newnc.variables[name][nctime,iz1:iz2,iy1:iy2,ix1:ix2]=data.T
            else:
                 if static:
                     newnc.variables[name][iy1:iy2,ix1:ix2]=data.T
                 else:
                     #print 'data.T.shape',data.T.shape
                     #print '[ix1,ix2,iy1,iy2,iz1,iz2]',[ix1,ix2,iy1,iy2,iz1,iz2]
                     #print 'imin,jmin,nctime',imin,jmin,nctime
                     #print 'newnc.variables[name][:].shape',newnc.variables[name][:].shape
                     #print 'imin,jmin,nctime',imin,jmin,nctime
                     newnc.variables[name][nctime,iy1:iy2,ix1:ix2]=data.T



            ##############################################################
            # add optionnal attributes
            ##############################################################

            if 'long_name' in kwargs:
                newnc.variables[name].long_name = kwargs['long_name']
            else:
                try:
                    newnc.variables[name].long_name = var.longname
                except :
                    print('no long_name specified')


            if 'units' in kwargs:
                newnc.variables[name].units = kwargs['units']
            else:
                try:
                    newnc.variables[name].units = var.unit
                except :
                    print('no long_name specified')
 

                

            ##############################################################
            # Add ocean_time
            ##############################################################




            if 'ocean_time' not in list(newnc.variables.keys()):
                newnc.createVariable('ocean_time', 'f', (dim0,) )
                


        if not static: newnc.variables['ocean_time'][nctime]=simul.oceantime

        newnc.close()



###################################################################################
#   Create netcdf file
###################################################################################

    def create(self,newfile,simul,netcdf_format='NETCDF4_CLASSIC',**kwargs):

        """

        """
        if 'dims' in kwargs: [nx,ny,nz] = kwargs['dims']
        else: [nx,ny,nz] = [simul.x.shape[0], simul.x.shape[1], len(simul.coord[4])]
        
        newnc = Dataset(newfile, 'w', format=netcdf_format)

        #[dim0, dim3, dim2, dim1]= ['time', ['depth_w','depth','depth_w_old'], ['lat','lat_v'], ['lon','lon_u']]
        [dim0, dim3, dim2, dim1]= ['time', ['s_w','s_rho','s_w_old'], ['eta_rho','eta_v'], ['xi_rho','xi_u']]

        newnc.createDimension(dim0, None) 
        newnc.createDimension(dim1[0], nx) 
        newnc.createDimension(dim2[0], ny) 
        newnc.createDimension(dim3[0], nz+1) 
        newnc.createDimension(dim1[1], nx-1) 
        newnc.createDimension(dim2[1], ny-1) 
        newnc.createDimension(dim3[1], nz) 
        newnc.createDimension(dim3[2], np.max([nz-1,1])) 
        print('create', newfile)
        
        return newnc





###################################################################################
#   Get variable
###################################################################################


    @staticmethod
    def get(ncfile,varname,simul,netcdf_format= 'NETCDF4_CLASSIC',**kwargs):

        """

        """
        if isinstance(ncfile, str): 
            nc = Dataset(ncfile, 'r', format=netcdf_format)
            opened=True
        else: 
            nc = ncfile
            opened = False

        dims = len(nc.variables[varname].dimensions)

            
        if 'time' in kwargs:
            time = kwargs['time']

            if dims==1:
                var = simul.Forder(nc.variables[varname][time])
            elif dims==2:
                var = simul.Forder(nc.variables[varname][time,:])
            elif dims==3:
                var = simul.Forder(nc.variables[varname][time,:,:])
            elif dims==4:
                var = simul.Forder(nc.variables[varname][time,:,:,:])

        else:

            var = simul.Forder(nc.variables[varname][:])
                
                
        if ('level' in kwargs) and (len(var.shape)==3):
            var = var[:,:,kwargs['level']]

            
            
        try: 
            var[var==nc.variables[varname]._FillValue] = np.nan
        except:
            try: 
                var[var==nc.variables[varname].fill_value] = np.nan          
            except:       
                #print 'no FillValue in file'
                pass
            
            
            
        if opened==True: nc.close()

        return np.squeeze(var)



###################################################################################
#   Create netcdf file
###################################################################################

    @staticmethod
    def write(newfile,var,nctime=None,name='var',shape = [0,0,0],netcdf_format= 'NETCDF4_CLASSIC',zlib=False,**kwargs):

        """

        """

        data = var
        
        ###################################################
        #var is a float (only for existing netcdf file)
        
        if isinstance(data,int) or isinstance(data,float):
            
            
            newnc = Dataset(newfile, 'a', format=netcdf_format)
            [dim0]= ['time']
            
            if nctime==None:  nctime = len(newnc.dimensions['time'])

            if name not in list(newnc.variables.keys()):
                newnc.createVariable(name, 'f', (dim0,) )
            newnc.variables[name][nctime]=data
            
            
            
        ###################################################
        #var is an array 
        
        else:
        
            if len(var.shape)>1: 
                Ly = var.shape[1]
                if len(var.shape)>2:
                    Lz = var.shape[2]
                else:
                    Lz=1
            else:
                Ly,Lz=1,1
            
            #[dim0, dim3, dim2, dim1]= ['time', ['depth'],['lat','lat_v'], ['lon','lon_u']]
            [dim0, dim3, dim2, dim1]= ['time', ['s_rho'], ['eta_rho','eta_v'], ['xi_rho','xi_u']]

            
            if not os.path.isfile(newfile): 
                newnc = Dataset(newfile, 'w', format=netcdf_format)
                newnc.createDimension(dim0, None) 
                newnc.createDimension(dim1[0], var.shape[0]+shape[0]) 
                newnc.createDimension(dim2[0], Ly+shape[1]) 
                newnc.createDimension(dim3[0], Lz) 
                newnc.createDimension(dim1[1], var.shape[0]+shape[0]-1) 
                newnc.createDimension(dim2[1], Ly+shape[1]-1) 
            else:  
                newnc = Dataset(newfile, 'a', format=netcdf_format)

            
            #[dim0, dim3, dim2, dim1]= ['time', ['depth'],['lat','lat_v'], ['lon','lon_u']]
            [dim0, dim3, dim2, dim1]= ['time', ['s_rho'], ['eta_rho','eta_v'], ['xi_rho','xi_u']]

            [imin,jmin,kmin] = shape

            if nctime==None:  nctime = len(newnc.dimensions['time'])



            print(name)

            if name not in list(newnc.variables.keys()):
                if netcdf_format=='NETCDF3_CLASSIC':
                    if len(data.shape)>2:
                        newnc.createVariable(name, 'f', (dim0,dim3[kmin],dim2[jmin],dim1[imin],) )
                    else:
                        newnc.createVariable(name, 'f', (dim0,dim2[jmin],dim1[imin],) )
                else:  
                    if len(data.shape)>2:
                        newnc.createVariable(name, 'f4', (dim0,dim3[kmin],dim2[jmin],dim1[imin],) ,zlib=zlib)
                    else:
                        newnc.createVariable(name, 'f4', (dim0,dim2[jmin],dim1[imin],) ,zlib=zlib)                  

            if len(data.shape)>2:
                 newnc.variables[name][nctime,:,:,:]=data.T
            else:
                 newnc.variables[name][nctime,:,:]=data.T

            try:
                print(simul.oceantime)
                if 'ocean_time' not in list(newnc.variables.keys()):
                    newnc.createVariable('ocean_time', 'f', (dim0,) )

                newnc.variables['ocean_time'][nctime]=simul.oceantime
                
            except:
                print('no oceantime in simul')

                
        newnc.close()


###################################################################################
#   Create netcdf file
###################################################################################

    @staticmethod
    def create_man(newfile,shape,netcdf_format= 'NETCDF4_CLASSIC',**kwargs):

        """

        """

        newnc = Dataset(newfile, 'w', format=netcdf_format)

        #[dim0, dim3, dim2, dim1]= ['time', ['depth'], ['lat','lat_v'], ['lon','lon_u']]
        [dim0, dim3, dim2, dim1]= ['time', ['s_rho'], ['eta_rho','eta_v'], ['xi_rho','xi_u']]

        newnc.createDimension(dim0, None) 
        newnc.createDimension(dim1[0], shape[0]) 
        newnc.createDimension(dim2[0], shape[1]) 
        newnc.createDimension(dim3[0], shape[2]) 
        newnc.createDimension(dim1[1], shape[0]-1) 
        newnc.createDimension(dim2[1], shape[1]-1) 
        
        return newnc


###################################################################################
#END CLASS LOAD
###################################################################################


'''
            if len(varnames)>1: newnc.createVariable(varnames[1], 'f', (dim0,dim1,dim2,) )
            if len(varnames)>2: newnc.createVariable(varnames[2], 'f', (dim0,dim1,dim2,) )
            if len(varnames)>3: newnc.createVariable(varnames[3], 'f', (dim0,dim1,dim2,) )
            if len(varnames)>4: newnc.createVariable(varnames[4], 'f', (dim0,dim1,dim2,) )
            if len(varnames)>5: newnc.createVariable(varnames[5], 'f', (dim0,dim1,dim2,) )
            newnc.flush()


        newnc.variables[varnames[0]][((time-time0)/deltatime)%ncsize,:,:]=np.float32(roms.psi2rho(var))
'''




