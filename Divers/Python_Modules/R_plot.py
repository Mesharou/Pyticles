###################################################################################
# R_PLOT
###################################################################################
"""


  Various functions useful for plotting with matplotlib

  Copyright (c) 2014 by Jonathan Gula
  e-mail:gula@ucla.edu
  
  Updated  14-01-30 ( add ncview_colormap, find_file )

  
"""

###################################################################################
#Load modules
###################################################################################

#for numeric functions
import numpy as np

#for plotting
import matplotlib.pyplot as py
import matplotlib.colors as col


###################################################################################

import os,sys


def findfile(path):
    '''Find the file named path in the sys.path.
    Returns the full path name if found, None if not found'''
    for dirname in sys.path:
        possible = os.path.join(dirname, path)
        if os.path.isfile(possible):
            return possible
    return None


###################################################################################


def nc_colormap(name):
    '''
    name is the string used in ncview for the colormap
    colormaps.h files needs to be in somewhere in the python path
    '''
    
    ncfile = findfile('colormaps_' + name + '.h')
    if ncfile==None: ncfile = findfile('./Colormaps/colormaps_' + name + '.h')
    if ncfile==None: print('no colormap file for ' + name)       
    
    f = open(ncfile, 'r')
        
    
    numbers = f.read().split("{")[1].split("}")[0].split(",")
    
    rgb=[]

    for i in range(len(numbers)/3):
        rgb.append((np.float(np.int(numbers[3*i])/255.),np.float(np.int(numbers[3*i+1])/255.), np.float(np.int(numbers[3*i+2])/255.)))

    my_cmap = col.LinearSegmentedColormap.from_list('my_colormap',rgb,256)          

    return my_cmap

###################################################################################


def nc_colormap_asym(name,ratio=2):
    '''
    name is the string used in ncview for the colormap
    colormaps.h files needs to be in somewhere in the python path
    '''

    ncfile = findfile('colormaps_' + name + '.h')
    if ncfile==None: ncfile = findfile('./Colormaps/colormaps_' + name + '.h')
    if ncfile==None: print('no colormap file for ' + name)

    f = open(ncfile, 'r')


    numbers = f.read().split("{")[1].split("}")[0].split(",")

    rgb=[]

    for i in range(len(numbers)/3):
        rgb.append((np.float(np.int(numbers[3*i])/255.),np.float(np.int(numbers[3*i+1])/255.), np.float(np.int(numbers[3*i+2])/255.)))

    print(len(numbers)/3, len(rgb))
    newrgb = rgb[:len(rgb)/2:ratio]

    for i in range(len(rgb)/2,len(rgb)):
        newrgb.append(rgb[i])

    my_cmap = col.LinearSegmentedColormap.from_list('my_colormap',newrgb,256)

    return my_cmap



###################################################################################


def nc_colormap_r(name):
    '''
    name is the string used in ncview for the colormap
    colormaps.h files needs to be in somewhere in the python path
    '''
    
    ncfile = findfile('colormaps_' + name + '.h')
    if ncfile==None: ncfile = findfile('./Colormaps/colormaps_' + name + '.h')
    if ncfile==None: print('no colormap file for ' + name)       
    
    f = open(ncfile, 'r')
        
    
    numbers = f.read().split("{")[1].split("}")[0].split(",")
    
    rgb=[]

    for i in range(len(numbers)/3):
        rgb.append((np.float(np.int(numbers[3*i])/255.),np.float(np.int(numbers[3*i+1])/255.), np.float(np.int(numbers[3*i+2])/255.)))

    rgb = rgb[::-1]

    my_cmap = col.LinearSegmentedColormap.from_list('my_colormap',rgb,256)          

    return my_cmap


###################################################################################


def nc_colormap_for_mlab(name):
    '''
    name is the string used in ncview for the colormap
    colormaps.h files needs to be in somewhere in the python path
    
    
    Output is a 255x4 array, with the columns representing RGBA
    (red, green, blue, alpha) coded witintegers going from 0 to 255.
    
    '''
    
    ncfile = findfile('colormaps_' + name + '.h')
    if ncfile==None: ncfile = findfile('./Colormaps/colormaps_' + name + '.h')
    if ncfile==None: print('no colormap file for ' + name)       
    
    f = open(ncfile, 'r')
        
    numbers = f.read().split("{")[1].split("}")[0].split(",")
    
    rgb=np.zeros((len(numbers)/3,4))

    for i in range(len(numbers)/3):
        rgb[i,0]=np.int(numbers[3*i])
        rgb[i,1]=np.int(numbers[3*i+1])
        rgb[i,2]=np.int(numbers[3*i+2])
        rgb[i,3]=255 #alpha Channel
        
    #my_cmap = col.LinearSegmentedColormap.from_list('my_colormap',rgb,256)          

    return rgb

###################################################################################


def fsu_colormap():
    '''


    '''
    myfile = findfile('colormap_fsu.txt')
    if myfile==None: myfile = findfile('./Colormaps/colormap_fsu.txt')

    f = open(myfile, 'r'); joe = f.read()  

    numbers = joe.split("\n")

    rgb=[]

    for i in range(len(numbers)-1):
        val = numbers[i].split(" ");
        while val.count('') > 0:
            val.remove('')
        rgb.append((np.float(val[0]),np.float(val[1]), np.float(val[2])))

    my_cmap = col.LinearSegmentedColormap.from_list('my_colormap',rgb,256)          

    return my_cmap


###################################################################################



class plot(object):

###################################################################################
#   Contour plot
###################################################################################


    def __init__(self,var,simul=None,N=50,sym=0,coef=1.,levels=[0],minmax=0,**kwargs):

        
        if type(var).__name__=='ndarray': data = var
        else: data = var.data

        if 'name' in kwargs: 
            name = kwargs['name']
        else:
            try:
                name = var.name
            except:
                name = 'unknown'


        if 'cmap' in kwargs: cmap=self.colormap(kwargs['cmap'] )
        else: cmap=self.colormap(name)
        
        if np.nanmin(levels)==np.nanmax(levels):
            if sym==1:
                if minmax==0: minmax = np.nanmax([np.abs(np.nanmin(data)),np.abs(np.nanmax(data))])*coef
                levels=np.arange(-1*minmax,minmax+minmax/(2*N+1),minmax/(2*N+1))
            else:
                levels=np.arange(np.nanmin(data),np.nanmax(data)+(np.nanmax(data)-np.nanmin(data))/N,(np.nanmax(data)-np.nanmin(data))/N)
        
        #self.plot = py.contourf(data,levels,
        #                        cmap=self.colormap('pv'),
        #                        extend='both'); 

        self.plot = py.pcolormesh(data,vmin=min(levels), vmax=max(levels),
                                cmap=self.colormap('pv'))

        py.colorbar(); py.show()



###################################################################################
# Fast plot different fields with the same scale
###################################################################################


    @staticmethod
    def p4(u,u2,u3,u4,x=None,y=None,N=10,samelev=False,coef=1., sym=0.,minmax=0):
        
        
        
        if sym==1:
            if minmax==0: minmax = np.nanmax([np.abs(np.nanmin(u)),np.abs(np.nanmax(u))])*coef
            levels=np.arange(-1*minmax,minmax+minmax/(2*N+1),minmax/(2*N+1))
        else:
            levels=np.arange(np.nanmin(u),np.nanmax(u)+(np.nanmax(u)-np.nanmin(u))/N,(np.nanmax(u)-np.nanmin(u))/N)

        py.subplot(2,2,1)
        if x==None: py.contourf(u,levels,extend='both'); py.colorbar();
        else: py.contourf(x,y,u,levels,extend='both'); py.colorbar();

        if samelev: 
            levels=np.arange(np.nanmin(u2),np.nanmax(u2)+(np.nanmax(u2)-np.nanmin(u2))/N,(np.nanmax(u2)-np.nanmin(u2))/N)
        py.subplot(2,2,2)
        if x==None: py.contourf(u2,levels,extend='both'); py.colorbar(); 
        else: py.contourf(x,y,u2,levels,extend='both'); py.colorbar();

        if samelev: 
            levels=np.arange(np.nanmin(u3),np.nanmax(u3)+(np.nanmax(u3)-np.nanmin(u3))/N,(np.nanmax(u3)-np.nanmin(u3))/N)
        py.subplot(2,2,3)
        if x==None: py.contourf(u3,levels,extend='both'); py.colorbar();
        else: py.contourf(x,y,u3,levels,extend='both'); py.colorbar(); 

        if samelev: 
            levels=np.arange(np.nanmin(u4),np.nanmax(u4)+(np.nanmax(u4)-np.nanmin(u4))/N,(np.nanmax(u4)-np.nanmin(u4))/N)
        py.subplot(2,2,4)
        if x==None: py.contourf(u4,levels,extend='both'); py.colorbar(); py.show()
        else: py.contourf(x,y,u4,levels,extend='both'); py.colorbar(); py.show()


###################################################################################
# Fast plot different fields with the same scale
###################################################################################


    @staticmethod
    def p2(u,u2,x=None,y=None,N=10,samelev=True,coef=1., sym=0.,minmax=0.):
        
        if minmax>0:
            levels=np.arange(-1*minmax,minmax+minmax/(2*N+1),minmax/(2*N+1))           
        elif sym==1:
            minmax = np.nanmax([np.abs(np.nanmin(u)),np.abs(np.nanmax(u))])*coef
            levels=np.arange(-1*minmax,minmax+minmax/(2*N+1),minmax/(2*N+1))
        else:
            levels=np.arange(np.nanmin(u),np.nanmax(u)+(np.nanmax(u)-np.nanmin(u))/N,(np.nanmax(u)-np.nanmin(u))/N)
            
        py.subplot(1,2,1)
        if x==None: py.contourf(u,levels,extend='both'); py.colorbar();
        else: py.contourf(x,y,u,levels,extend='both'); py.colorbar();


        if not samelev:
            py.colorbar();
            levels=np.arange(np.nanmin(u2),np.nanmax(u2)+(np.nanmax(u2)-np.nanmin(u2))/N,(np.nanmax(u2)-np.nanmin(u2))/N)
        py.subplot(1,2,2)
        if x==None: py.contourf(u2,levels,extend='both'); py.colorbar(); py.show()
        else: py.contourf(x,y,u2,levels,extend='both'); py.colorbar(); py.show()

        
        
###################################################################################
#define colorabar levels using levelsvar
###################################################################################

    @staticmethod
    def clabels(levelsvar,nblab=4.,sym=0):

        tot=levelsvar.max()-levelsvar.min()
        test=0; i=0
        eps = tot*1e-6
        while test==0:
            test=round(tot,i)
            i=i+1

        dlab=round(tot/nblab,i)
        labmin=round(levelsvar.min(),i)
        labmax=round(levelsvar.max(),i)

        labels=np.arange(labmin,labmax+dlab,dlab)
        
        # Impose that 0 is among labels if sym=1
        if sym==1:
            labels = np.arange(-2*dlab,3*dlab,dlab)
    
        return labels[np.logical_and(labels<=levelsvar.max()+eps,labels>=levelsvar.min()-eps)]
    
    

###################################################################################
#Load some custom colormaps defined for each variable
###################################################################################

    @staticmethod
    def colormap(vavar):


        cdict5 = {'blue':  ((0.0, 0.5, 0.5),
                           (0.2,1.0, 1.0),
                           (0.45, 1.0, 1.0),  
                           (0.55, 1.0, 1.0),                    
                           (0.8,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'green': ((0.0, 0.0, 0.0),
                           (0.2, 0.0, 0.0),
                           (0.45, 1.0, 1.0),
                           (0.55, 1.0, 1.0),          
                           (0.8,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'red':  ((0.0, 0.0, 0.0),
                           (0.2,0.0, 0.0),
                           (0.45, 1.0, 1.0),  
                           (0.55, 1.0, 1.0),                   
                           (0.8,1.0, 1.0),
                           (1.0, 0.5, 0.5))
        }                  
        cdict_pv_asym = {'blue':  ((0.0, 0.5, 0.5),
                           (0.1,1.0, 1.0),
                           (0.225, 1.0, 1.0),
                           (0.275, 1.0, 1.0),
                           (0.4,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'green': ((0.0, 0.0, 0.0),
                           (0.1, 0.0, 0.0),
                           (0.225, 1.0, 1.0),
                           (0.275, 1.0, 1.0),
                           (0.4,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'red':  ((0.0, 0.0, 0.0),
                           (0.1,0.0, 0.0),
                           (0.225, 1.0, 1.0),
                           (0.275, 1.0, 1.0),
                           (0.4,1.0, 1.0),
                           (1.0, 0.5, 0.5))
        }
  
        cdict6 = {'blue':  ((0.0, 0.5, 0.5),
                           (0.2,1.0, 1.0),
                           (0.48, 1.0, 1.0),  
                           (0.52, 0.0, 0.0),                  
                           (0.8,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'green': ((0.0, 0.0, 0.0),
                           (0.2, 0.0, 0.0),
                           (0.48, 1.0, 1.0), 
                           (0.52, 1.0, 1.0),                   
                           (0.8,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'red':  ((0.0, 0.0, 0.0),
                           (0.2,0.0, 0.0),
                           (0.48, 0.0, 0.0), 
                           (0.52, 1.0, 1.0),                   
                           (0.8,1.0, 1.0),
                           (1.0, 0.5, 0.5))
        }
        cdict7 = {'blue':  ((0.0, 0.5, 0.5),
                           (0.2,1.0, 1.0),
                           (0.49, 1.0, 1.0),  
                           (0.51, 1.0, 1.0),                    
                           (0.8,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'green': ((0.0, 0.0, 0.0),
                           (0.2, 0.0, 0.0),
                           (0.49, 1.0, 1.0),
                           (0.51, 1.0, 1.0),          
                           (0.8,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'red':  ((0.0, 0.0, 0.0),
                           (0.2,0.0, 0.0),
                           (0.49, 1.0, 1.0),  
                           (0.51, 1.0, 1.0),                   
                           (0.8,1.0, 1.0),
                           (1.0, 0.5, 0.5))
        }                     
        cdict8 = {'blue':  ((0.0,0.1,0.1),
                           (0.6, 1.0, 1.0),  
                           (0.8, 0.0, 0.0),                  
                           (0.92,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'green': ((0.0, 0.0, 0.0),
                           (0.1, 0.0, 0.0),
                           (0.6, 1.0, 1.0), 
                           (0.8, 1.0, 1.0),                   
                           (0.92,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'red':  ((0.0, 0.0, 0.0),
                           (0.1,0.0, 0.0),
                           (0.6, 0.0, 0.0), 
                           (0.8, 1.0, 1.0),                   
                           (0.92,1.0, 1.0),
                           (1.0, 0.5, 0.5))
        }
        cdict9 = {'blue':  ((0.0,0.1,0.1),
                           (0.6, 1.0, 1.0),  
                           (0.8, 0.0, 0.0),                  
                           (0.8,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'green': ((0.0, 0.0, 0.0),
                           (0.0, 0.0, 0.0),
                           (0.6, 1.0, 1.0), 
                           (0.8, 1.0, 1.0),                   
                           (1.0, 0.1, 0.1)),

                 'red':  ((0.0, 0.0, 0.0),
                           (0.0,0.0, 0.0),
                           (0.6, 0.0, 0.0), 
                           (0.8, 1.0, 1.0),                   
                           (0.95,1.0, 1.0),
                           (1.0, 0.8, 0.8))
        }
        cdict9bis = {'blue':  ((0.0,0.2,0.2),
                           (0.6, 1.0, 1.0),  
                           (0.8, 0.0, 0.0),                  
                           (0.8,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'green': ((0.0, 0.0, 0.0),
                           (0.0, 0.0, 0.0),
                           (0.6, 1.0, 1.0), 
                           (0.8, 1.0, 1.0),                   
                           (1.0, 0.1, 0.1)),

                 'red':  ((0.0, 0.0, 0.0),
                           (0.0,0.0, 0.0),
                           (0.6, 0.0, 0.0), 
                           (0.8, 1.0, 1.0),                   
                           (0.95,1.0, 1.0),
                           (1.0, 0.8, 0.8))
        }
        cdict10 = {'blue':  ((0.0, 1.0, 1.0),
                           (0.0,1.0, 1.0),
                           (0.5, 1.0, 1.0),                    
                           (1.0,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'green': ((0.0, 0.0, 0.0),
                           (0.0, 0.0, 0.0),
                           (0.5, 1.0, 1.0),          
                           (1.0,0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'red':  ((0.0, 0.0, 0.0),
                           (0.0,0.0, 0.0),
                           (0.5, 1.0, 1.0),                   
                           (1.0,1.0, 1.0),
                           (1.0, 1.0, 1.0))
        }                
        cdict_meanvel = {'blue':  ((0.0,1.0,1.0),
                           (0.2, 1.0, 1.0),    
                           (0.4, 1.0, 1.0),   
                           (0.5, 0.0, 0.0),
                           (0.6, 0.0, 0.0),
                           (0.8, 0.0, 0.0),    
                           (1.0, 0.0, 0.0)),

                 'green': ((0.0, 1.0, 1.0),
                           (0.3, 1.0, 1.0),
                           (0.4, 0.5, 0.5),
                           (0.5, 0.5, 0.5),
                           (0.6, 1.0, 1.0),                 
                           (0.7, 1.0, 1.0), 
                           (1.0, 0.0, 0.0)),

                 'red':  ((0.0, 1.0, 1.0),
                           (0.2, 1.0, 1.0), 
                           (0.4, 0.0, 0.0), 
                           (0.5, 0.0, 0.0),  
                           (0.6, 0.0, 0.0),  
                           (0.7, 1.0, 1.0),
                           (0.8, 0.5, 0.5),
                           (1.0, 1.0, 1.0))
        }
        cdict_meanvel_du = {'blue':  ((0.0, 1.0, 1.0),
                           (0.4, 1.0, 1.0),  
                           (0.6, 1.0, 1.0),                    
                           (1.0, 0.0, 0.0)),

                 'green': ((0.0, 0.0, 0.0),
                           (0.4, 1.0, 1.0),
                           (0.6, 1.0, 1.0),          
                           (1.0, 0.0, 0.0)),

                 'red':  ((0.0, 0.0, 0.0),
                           (0.4, 1.0, 1.0),  
                           (0.6, 1.0, 1.0),                   
                           (1.0, 1.0, 1.0))
        }
        cdict_meanu = {'blue':  ((0.0,1.0,1.0),
                           (0.3, 1.0, 1.0),    
                           (0.4, 1.0, 1.0),   
                           (0.5, 0.0, 0.0),
                           (0.6, 0.0, 0.0),
                           (0.8, 0.0, 0.0),    
                           (1.0, 0.0, 0.0)),

                 'green': ((0.0, 1.0, 1.0),
                           (0.3, 1.0, 1.0),
                           (0.4, 0.5, 0.5),
                           (0.5, 0.5, 0.5),
                           (0.6, 1.0, 1.0),                 
                           (0.7, 1.0, 1.0), 
                           (1.0, 0.0, 0.0)),

                 'red':  ((0.0, 1.0, 1.0),
                           (0.3, 1.0, 1.0), 
                           (0.4, 0.0, 0.0), 
                           (0.5, 0.0, 0.0),  
                           (0.6, 0.0, 0.0),  
                           (0.7, 1.0, 1.0),
                           (0.8, 0.5, 0.5),
                           (1.0, 1.0, 1.0))
        }
                 


        if vavar in ['str','salt','dsalt']:
            my_cmap = py.cm.jet
        # Classic Blue-white-red colorscale
        elif vavar in ['div', 'dvdz', 'pv', 'pvF', 'pv1', 'pv2', 'pv3', 'J1', 'J2', 'Jbot', 'J1old', 'J2old', 'omega', 'absvrt', 'vrt', 'pvold', 'dxpv', 'wb', 'u', 'uper', 'v']:
            my_cmap = col.LinearSegmentedColormap('my_colormap',cdict5,256)
        # Blue-yellow-red colorscale
        elif vavar in ['vrt', 'ageo', 'ow']:
            my_cmap = col.LinearSegmentedColormap('my_colormap',cdict6,256)
        # Blue-white-red colorscale with less white (to see pos-neg)
        elif vavar in ['test']:
            my_cmap = col.LinearSegmentedColormap('my_colormap',cdict7,256)
        # blask-yellow-red 
        elif vavar in ['temp']:
            my_cmap = col.LinearSegmentedColormap('my_colormap',cdict9,256)
        # blue-yellow-red 
        elif vavar in ['temp_blue']:
            my_cmap = col.LinearSegmentedColormap('my_colormap',cdict9bis,256)
        elif vavar in ['w']:
            my_cmap = col.LinearSegmentedColormap('my_colormap',cdict10,256)
        elif vavar in ['meanvel']:
            my_cmap = col.LinearSegmentedColormap('my_colormap',cdict_meanvel,256)
        elif vavar in ['meanvel_du']:
            my_cmap = col.LinearSegmentedColormap('my_colormap',cdict_meanvel_du,256)
        elif vavar in ['meanu']:
            my_cmap = col.LinearSegmentedColormap('my_colormap',cdict_meanu,256)
        elif vavar in ['bwyr']:
            my_cmap = col.rgb.from_list('my_colormap',['blue','white','yellow','red'],256)
        elif vavar in ['pv_asym']:
	    my_cmap = col.LinearSegmentedColormap('my_colormap',cdict_pv_asym,256)

        else:
            # Load a NCVIEW colormap
            try:
                my_cmap = nc_colormap(vavar)  
            except:        
                my_cmap = py.cm.spectral 

        return my_cmap



###################################################################################
# Load ncview colormaps
###################################################################################


    ########################
   
    @staticmethod
    def ncview_colormap(name):
        '''
        name is the string used in ncview for the colormap
        colormaps.h files needs to be in somewhere in the python path
        '''
        
        my_cmap = nc_colormap(name)  

        return my_cmap


    ########################
   
    @staticmethod
    def ncview_colormap_r(name):
        '''
        name is the string used in ncview for the colormap
        colormaps.h files needs to be in somewhere in the python path
        '''
        
        my_cmap = nc_colormap_r(name)  

        return my_cmap

    ########################

    @staticmethod
    def ncview_colormap_asym(name,ratio=2):
        '''
        name is the string used in ncview for the colormap
        colormaps.h files needs to be in somewhere in the python path
        '''

        my_cmap = nc_colormap_asym(name,ratio=ratio)

        return my_cmap

    
###################################################################################
# Load fsu colormap
###################################################################################


    ########################
   
    @staticmethod
    def fsu_colormap():
        '''

        '''
        
        my_cmap = fsu_colormap()  

        return my_cmap
    
    
    
    
    
    
