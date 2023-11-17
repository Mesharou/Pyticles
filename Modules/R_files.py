
"""
Routines used to load ROMS simulations parameters, use it in a python script as follows:

#######################

From R_files import load

my_simul = load(simul = 'filam [0,1000,0,1600,[-100,0,10]] 190')

or

my_simul = load('filam') #default being all domain and time=0

#######################


Where Arguments are in this order:

1 - *not used here*
2 - name of the simulation (should be defined in the class files, see below)
3 - domain size defined as [ymin,ymax,xmin,xmax,[zmin,zmax,deltaz]] or [ymin,ymax,xmin,xmax,[z]]
    if z >0 : z is considered as a sigma level
    if z <0 : z is depth in meter
    if z =0 : z = surface
4 - *not used here*
5 - initial time (corresponding to ROMS output files number)
6 - delta time (optional, default is 1)
5 - end time (optional)

#######################

All parameters such as name, date, static fields, files path ...etc..
are loaded in the form of "my_simul" attributes

Check "dir(my_simul)" to see what has been loaded:

['Cs_r', 'Cs_w', 'Forder',  'coord', 'coordmax', 'cst', 'date', 'day', 'domain', 'dt', 'dtime', 'f', 'filetime', 'floattype', 'g', 'get_domain', 'hc', 'hour', 'infiletime', 'load_file', 'mask', 'min', 'month', 'ncfile', 'ncname', 'oceandate', 'oceantime', 'pm', 'pn', 'rdrg', 'rho0', 'simul', 'time', 'time0', 'topo', 'update', 'variables_grd', 'x', 'y', 'year']

#######################

16/10/11 : Changed format from 
           my_simul = load(simul = 'whatever filam [0,1000,0,1600,[-100,0,10]] whatever 190 1 250')
           to
           my_simul = load(simul = 'filam [0,1000,0,1600,[-100,0,10]] 190') or my_simul = load('filam') 
           

17/11/01 : add option realyear and modify oceandate function
           
"""
from __future__ import print_function



###################################################################################
#Load some common modules
###################################################################################

#from builtins import object
import sys

#for netcdf files
from netCDF4 import Dataset

#for numeric functions
import numpy as np

import time as tm

from datetime import datetime, timedelta

import socket


###################################################################################


class load(object):
    """
    Load ROMS/CROCO files and variables

    Parameters:
    ======================================================================
        - str, simulname: name of ROMS/CROCO simulation defined in Class
            files
        - float, time: time to load files
        - bool, light: if True load only surface variables
            (used in GIGATL_surf) to incease performance
        - bool, touchfile: if False it doesnt matter if first file is on
            disk or not
        - bool, output: increase verbosity

    """

###################################################################################
#   Main 
###################################################################################
    def __init__(self, simulname=None, time=None, floattype=float, light=False,
                touchfile=True, output=True, **kwargs):

        if output: print('simulname is',simulname)

        #define type of variables
        self.floattype = floattype

        #get simulation parameters from arguments
        if 'simul' in  kwargs: 
            self.load_file(kwargs['simul'].split(' '), output=output)
        #elif len(sys.argv)>0:
        #    self.load_file(sys.argv)
        else: 
            self.load_file([simulname], output=output)
        
        #for time in range(self.time0, self.ncname.tend, self.dtime ):
        if time==None: time = self.time0
        self.time = int(time)
        
        self.ncname = files(self.simul, time=self.time, output = output)
        
        if self.ncname.model in ['ucla', 'croco']:
            if output: print('time of simulation is:', self.time)
            self.infiletime = self.time%self.ncname.tfile
            self.filetime = self.time-self.infiletime
            if self.ncname.digits == 4:
                self.ncfile = self.ncname.his+'{0:04}'.format(self.filetime) + self.ncname.fileformat
            elif self.ncname.digits == 5:
                self.ncfile = self.ncname.his+'{0:05}'.format(self.filetime) + self.ncname.fileformat
            elif self.ncname.digits == 6:
                self.ncfile = self.ncname.his+'{0:06}'.format(self.filetime) + self.ncname.fileformat
            elif self.ncname.digits == 0:
                self.ncfile = self.ncname.his + self.ncname.fileformat
        elif self.ncname.model == 'agrif':
            dsec = 30*24*3600 //self.ncname.tfile #dt in seconds between 2 outputs
            month = (self.ncname.Mstart + self.time * dsec// (30*24*3600) - 1 )%12 + 1
            year = self.ncname.Ystart + ((self.ncname.Mstart-1) * 30 * 24 * 3600 + self.time * dsec ) // (360*24*3600)
            self.infiletime = (self.time * dsec % (30*24*3600))//dsec
            #self.infiletime = (self.time * dhour % (30*24))/dhour
            self.filetime=self.time
            self.ncfile = self.ncname.his+'Y' +format(year)+'M'+format(month) + self.ncname.fileformat
        elif self.ncname.model == 'croco_lionel':
            date = self.ncname.realyear_origin + timedelta(days=self.time-self.ncname.tstart)
            year = date.year
            month = date.month
            self.infiletime = self.time%self.ncname.tfile
            self.filetime = self.time
            self.ncfile = self.ncname.his+'Y' +'{0:04}'.format(year)+'M'+'{0:02}'.format(month) + '.' +'{0:04}'.format(self.filetime) + self.ncname.fileformat
            if output: print('file is ',self.ncfile)
        elif self.ncname.model == 'croco_xios' or 'croco_gigatl1' in self.ncname.model:
            # find first and last date in file to reconstruct name, ex: 1999-01-25-1999-01-29
            time1 = self.time - self.time%self.ncname.tfile
            date1 = self.ncname.realyear_tstart + timedelta(seconds=time1*self.ncname.dtfile)
            year1 = date1.year
            month1 = date1.month
            day1 = date1.day
            hour1 = date1.hour
            
            time2 = time1 + self.ncname.tfile #- 1

            date2 = self.ncname.realyear_tstart + timedelta(seconds=time2*self.ncname.dtfile)

            ###############
            if 'rrexnum' not in self.simul and\
                'brazil' not in self.simul: date2 = date2 -timedelta(days=1)
            ###############
            
            year2 = date2.year
            month2 = date2.month
            day2 = date2.day

 
            #fix for GIGATL3_6h_knl! shifted by one day.
            if (self.simul == 'gigatl3_6h' and self.time>=1550)\
                 or 'gigatl6_1h_tides' in self.simul:
                if output: print('fix for GIGATL3_6h_knl and GIGATL6_1h_tides')
                date1 = date1 - timedelta(days = 1)
                year1 = date1.year
                month1 = date1.month
                day1 = date1.day

            

            self.infiletime = self.time%self.ncname.tfile
            self.filetime = self.time
            
            if self.ncname.model == 'croco_xios':
                self.ncfile = self.ncname.his\
                        + '{0:04}'.format(year1)+'-'+'{0:02}'.format(month1) + '-' +'{0:02}'.format(day1)\
                   + '-'+ '{0:04}'.format(year2)+'-'+'{0:02}'.format(month2) + '-' +'{0:02}'.format(day2)\
                        + self.ncname.fileformat
            elif self.ncname.model == 'croco_gigatl1':
                self.ncfile = self.ncname.his\
                        + '{0:04}'.format(year1)+'-'+'{0:02}'.format(month1) + '-' +'{0:02}'.format(day1)+ self.ncname.fileformat
            elif self.ncname.model == 'croco_gigatl1_hours':
                self.ncfile = self.ncname.his\
                        + '{0:04}'.format(year1)+'-'+'{0:02}'.format(month1)\
                        + '-' +'{0:02}'.format(day1)+ '-' +'{0:02}'.format(hour1)\
                        + self.ncname.fileformat

        if touchfile:
            try:
                self.oceandate() #define self.date and self.oceantime
            except:
                if output: print("no time in file")

        try:
            if touchfile:
                if output: print('coord')
                self.coord=self.get_domain(self.ncfile,self.ncname.grd,self.domain\
                            ,time,output = output)
                if output: print('coordmax')
                self.coordmax=self.get_domain(self.ncfile, self.ncname.grd,\
                            '[0,1e9,0,1e9,[1,1e9,1]]', time ,output = output)
                if output: print('cst')
                self.cst(output=output);
                if output: print('dt')
                self.dt(output=output)
            else:
                self.coord=[0,-1,0,-1,[0]]
                if output: print(self.ncfile)
                if output: print("not touching _his file, loading _grd anyway")
        except:
            self.coord=[0,-1,0,-1,[0]]
            if output: print(self.ncfile)
            if output: print("no _his file, loading _grd anyway")
            
        self.getin() #load some data from .in file
     
        #self.grd = [self.topo,self.pm,self.pn,self.f,self.y,self.x]
        if light:
            [self.mask,self.pm,self.pn,self.x,self.y] = self.variables_grd(output = output, light=True)
        else:
            [self.topo,self.mask,self.pm,self.pn,self.f,self.x,self.y,self.angle] = self.variables_grd(output = output)

        if self.simul in ['shing','reshing','filam','dilam','filam_avg','dilam_avg','slope','hatt2','gulfs']: 
            self.topo_corrected = self.variables_grd_corrected()

###################################################################################
#   Update time 
###################################################################################
    def update(self, time=None, output =True):

        """
        
        Update 'my_simul' to the next timestep
        (path to files, infiletime, filetime, ncfile, oceandate, date ...etc...)

        """

        #for time in range(self.time0, self.ncname.tend, self.dtime ):
        if time==None: self.time += 1
        else: self.time = int(time)
        
        self.ncname=files(self.simul, time=self.time,output = output)

        if self.ncname.model in ['ucla','croco']:
            if output: print('time of simulation is:', self.time)
            self.infiletime=self.time%self.ncname.tfile
            self.filetime=self.time-self.infiletime
            if self.ncname.digits==4:
                self.ncfile = self.ncname.his+'{0:04}'.format(self.filetime) + self.ncname.fileformat
            elif self.ncname.digits==5:
                self.ncfile = self.ncname.his+'{0:05}'.format(self.filetime) + self.ncname.fileformat
            elif self.ncname.digits==6:
                self.ncfile = self.ncname.his+'{0:06}'.format(self.filetime) + self.ncname.fileformat
            elif self.ncname.digits==0:
                self.ncfile = self.ncname.his + self.ncname.fileformat
        elif self.ncname.model == 'agrif_jc':
            dsec = 30*24*3600 /self.ncname.tfile #dt in seconds between 2 outputs
            month = (self.ncname.Mstart + self.time * dsec/ (30*24*3600) - 1 )%12 + 1
            year = self.ncname.Ystart + ((self.ncname.Mstart-1) * 30 * 24 * 3600 + self.time * dsec ) / (360*24*3600)
            self.infiletime = (self.time * dsec % (30*24*3600))/dsec
            #self.infiletime = (self.time * dhour % (30*24))/dhour
            self.filetime=self.time
            self.ncfile = self.ncname.his+'Y' +format(year)+'M'+format(month)+ self.ncname.fileformat
            if output: print('file is ',self.ncfile)
        elif self.ncname.model == 'croco_lionel':
            date = self.ncname.realyear_origin + timedelta(days=self.time-self.ncname.tstart)
            year = date.year
            month = date.month
            self.infiletime = self.time%self.ncname.tfile
            self.filetime = self.time
            self.ncfile = self.ncname.his+'Y' +'{0:04}'.format(year)+'M'+'{0:02}'.format(month) + '.' +'{0:04}'.format(self.filetime) + self.ncname.fileformat
            if output: print('file is ',self.ncfile)
        elif self.ncname.model == 'croco_xios' or 'croco_gigatl1' in self.ncname.model:
            # find first and last date in file to reconstruct name, ex: 1999-01-25-1999-01-29
            time1 = self.time - self.time%self.ncname.tfile
            date1 = self.ncname.realyear_tstart + timedelta(seconds=time1*self.ncname.dtfile)
            year1 = date1.year
            month1 = date1.month
            day1 = date1.day
            hour1 = date1.hour
            
            time2 = time1 + self.ncname.tfile #- 1
            date2 = self.ncname.realyear_tstart + timedelta(seconds=time2*self.ncname.dtfile)
            
            
            ###############
            if 'rrexnum' not in self.simul and\
                'brazil' not in self.simul: date2 = date2 -timedelta(days=1)
            ###############
                
            year2 = date2.year
            month2 = date2.month
            day2 = date2.day
            
            self.infiletime = self.time%self.ncname.tfile
            self.filetime = self.time-self.infiletime
            
            if output: print((self.simul,format(self.time)))
            #fix for GIGATL3_6h_knl! shifted by one day.
            if (self.simul == 'gigatl3_6h' and self.time>=1550)\
                 or 'gigatl6_1h_tides' in self.simul:
                if output: print('fix for GIGATL3_6h_knl and GIGATL6_1h_tides')
                date1 = date1 - timedelta(days = 1)
                year1 = date1.year
                month1 = date1.month
                day1 = date1.day

            if self.ncname.model == 'croco_xios':
                self.ncfile = self.ncname.his\
                        + '{0:04}'.format(year1)+'-'+'{0:02}'.format(month1) + '-' +'{0:02}'.format(day1)\
                   + '-'+ '{0:04}'.format(year2)+'-'+'{0:02}'.format(month2) + '-' +'{0:02}'.format(day2)\
                        + self.ncname.fileformat
            elif self.ncname.model == 'croco_gigatl1':
                self.ncfile = self.ncname.his\
                        + '{0:04}'.format(year1)+'-'+'{0:02}'.format(month1) + '-' +'{0:02}'.format(day1)+ self.ncname.fileformat
            elif self.ncname.model == 'croco_gigatl1_hours':
                self.ncfile = self.ncname.his\
                        + '{0:04}'.format(year1)+'-'+'{0:02}'.format(month1)\
                        + '-' +'{0:02}'.format(day1)+ '-' +'{0:02}'.format(hour1)\
                        + self.ncname.fileformat
            
            if output: print('self.time is ',self.time)

            ############################
            
        #print('cst not updated anymore')
        #self.cst();

        if output: print(self.ncfile)
        
        try:
            self.oceandate()
        except: 
            if output: print("no oceantime in file")
            

###################################################################################
# GET SIMULATION PARAMETERS
###################################################################################


    def load_file(self, *args, **kwargs):

        if 'output' in  kwargs:
            output = kwargs['output']
        else:
            output = True

        if output:
            print('args',args)
            print('args[0]',args[0])
            print('len(args[0])',len(args[0]))

        ########################
        if len(args[0])==0:

            if output: print("""
Try again with:
            1. Simulation name 
            2. domain location
            3. time
     optional:
            4. time-step
            5. final step 

example: for an interactive session:
         load(simul = 'filam [0,1000,0,1600,[-100,0,10]] 190')

                   """)
            #sys.exit()
            raise Exception("No args provided")

        elif len(args[0])==1:

            self.simul=args[0][0]
            self.domain='[0,10000,0,10000,[1,1000,1]]'
            self.ncname=files(self.simul, time=0, output = output)
            try:
                self.time0=self.ncname.tstart
            except:
                self.time0=0

        else:
            self.simul=args[0][0]
            self.domain=args[0][1]
            self.time0=int(args[0][2])
            self.ncname=files(self.simul, time=self.time0, output = output)


        if len(args[0])>3: 
            self.dtime=int(args[0][3])
            if len(args[0])>4: 
                self.ncname.tend=int(args[0][4])

        else: self.dtime=1000


    
###################################################################################
# GET DOMAIN
###################################################################################

    @staticmethod
    def get_domain(ncname,ncname0,domainname,time, output =True,*args):

        if output: print('loading', ncname0)
        ncfile0 = Dataset(ncname0, 'r')
        if output: print('loading', ncname)
        ncfile = Dataset(ncname, 'r')
        
        if output: print('get domain', domainname, domainname[:5])
        
        if domainname[:5] in ['filam','shing','slope']:
            
            import simulations_old as oldsim

            if output: print('loading custom domain using old scripts')
            [ny1,ny2,nx1,nx2,depths]=oldsim.domain(ncfile0,domainname,time)

        else:
            
            [ny1,ny2,nx1,nx2] = eval(domainname)[0:4]

            if len(eval(domainname)[4])==1:
                depths=eval(domainname)[4]
            else:
              try:
                depths=np.arange(eval(domainname)[4][0],int(np.min([len(ncfile.dimensions['s_rho']),eval(domainname)[4][1]]))+1,eval(domainname)[4][2])
              except:
                depths=[0]

        if ny1==ny2: ny1=ny1-1; ny2=ny2+2
        if nx1==nx2: nx1=nx1-1; nx2=nx2+2

        nx2 = int(np.min([nx2,len(ncfile0.dimensions['xi_rho'])]))
        ny2 = int(np.min([ny2,len(ncfile0.dimensions['eta_rho'])]))

        ncfile0.close()
        ncfile.close()


        return [ny1,ny2,nx1,nx2,depths]


###################################################################################
#Load grd variables
###################################################################################

    def variables_grd(self, light = False, output =True):

        if output: print(self.coord)

        [ny1,ny2,nx1,nx2] = self.coord[0:4]
        #ncfile0 = NetCDFFile(ncname,'r')
        ncfile0 = Dataset(self.ncname.grd, 'r')
    
        if output: print('ncname0,ny1,ny2,nx1,nx2')
        if output: print(self.ncname.grd,ny1,ny2,nx1,nx2)


        if ny2==-1 and nx2==-1:
            [ny2,nx2] = ncfile0.variables['h'].shape
            self.coord[0:4] = [ny1,ny2,nx1,nx2]

        pm = self.Forder(ncfile0.variables['pm'][ny1:ny2,nx1:nx2])
        pn = self.Forder(ncfile0.variables['pn'][ny1:ny2,nx1:nx2])


        try:
            mask = self.Forder(ncfile0.variables['mask_rho'][ny1:ny2,nx1:nx2])
        except:
            mask =  pm*0.+1.; 

        mask[mask==0] = np.nan

        
        if 'lon_rho' in list(ncfile0.variables.keys()):
            lon = self.Forder(ncfile0.variables['lon_rho'][ny1:ny2,nx1:nx2])
            lat = self.Forder(ncfile0.variables['lat_rho'][ny1:ny2,nx1:nx2])
        else:
            lon = self.Forder(ncfile0.variables['x_rho'][ny1:ny2,nx1:nx2])
            lat = self.Forder(ncfile0.variables['y_rho'][ny1:ny2,nx1:nx2])

        if not light: 
            topo = self.Forder(ncfile0.variables['h'][ny1:ny2,nx1:nx2])
            f = self.Forder(ncfile0.variables['f'][ny1:ny2,nx1:nx2])
            try:
                angle = self.Forder(ncfile0.variables['angle'][ny1:ny2,nx1:nx2])
            except:
                angle = f*0.
            
       
        ncfile0.close()

        if output: print('[topo,pm,pn,f,lat,lon] have just been loaded')
        if output: print('----------------------------------------------------------')
        if output: print('All arrays are now Fortran ordered and indices are [i,j,k]')
        if output: print('----------------------------------------------------------')

        if output and not light: print(topo.shape)
        
        if not light:
            return [topo,mask,pm,pn,f,lon,lat,angle]
        else:
            return [mask,pm,pn,lon,lat]

###################################################################################
#Load grd variables
###################################################################################

    def variables_grd_corrected(self, output =True):

        [ny1,ny2,nx1,nx2] = self.coord[0:4]
        ncfile0 = Dataset(self.ncname.grd_corrected, 'r')
        topo = self.Forder(ncfile0.variables['h'][ny1:ny2,nx1:nx2])
        ncfile0.close()

        return topo
    
    
    
###################################################################################
#Load some constant
###################################################################################


    def cst(self, output =True):

        try:
            ncfile = Dataset(self.ncfile, 'r')
        except:
            if output: print('cannot find: ', self.ncfile)

        try:
            self.rho0 = ncfile.rho0 
            self.g = 9.81
            self.hc = ncfile.hc
            self.dt_model = ncfile.dt
        except:
            pass


        try:
            self.Cs_r = self.Forder(ncfile.variables['Cs_r'][:])
            self.Cs_w = self.Forder(ncfile.variables['Cs_w'][:])
            self.sc_r = self.Forder(ncfile.variables['sc_r'][:])
            self.sc_w = self.Forder(ncfile.variables['sc_w'][:])
            if output: print('read Cs_r in ncfile.variables')
        except:
            try:
                self.Cs_r = self.Forder(ncfile.Cs_r)
                self.Cs_w = self.Forder(ncfile.Cs_w)
                self.sc_r = self.Forder(ncfile.sc_r)
                self.sc_w = self.Forder(ncfile.sc_w)
                if output: print('read Cs_r in ncfile.Cs_r')

            except:
                try:
                    grdfile = Dataset(self.ncname.grd, 'r')
                    self.Cs_r = self.Forder(grdfile.variables['Cs_r'][:])
                    self.Cs_w = self.Forder(grdfile.variables['Cs_w'][:])
                    grdfile.close()
                    if output: print('read Cs_r in grdfile.variables')
                except:
                    try:
                        grdfile = Dataset(self.ncname.grd, 'r')
                        self.Cs_r = self.Forder(grdfile.Cs_r)
                        self.Cs_w = self.Forder(grdfile.Cs_w)
                        grdfile.close()
                        if output: print('read Cs_r in grdfile.Cs_r')
                    except:
                        pass
        try:
            if np.ndim(self.Cs_r)==2:
                self.Cs_r = self.Cs_r[:,0]
                self.Cs_w = self.Cs_w[:,0]

            
            #ugly fix for now because of wrong xios files
            if self.Cs_r.max()>1:
                if output: print('really??')
                try:
                    grdfile = Dataset(self.ncname.grd, 'r')
                    self.Cs_r = self.Forder(grdfile.Cs_r)
                    self.Cs_w = self.Forder(grdfile.Cs_w)
                    grdfile.close()
                except:
                    grdfile = Dataset(self.ncname.grd, 'r')
                    self.Cs_r = self.Forder(grdfile.variables['Cs_r'][:])
                    self.Cs_w = self.Forder(grdfile.variables['Cs_w'][:])
                    grdfile.close()
        except:
            pass
            

        try:
            self.rdrg = ncfile.rdrg
        except:
            self.rdrg = 0.

        try:
            self.rdrg2 = ncfile.rdrg2
        except:
            self.rdrg2 = 0.


        try:
            self.Cdb_max = ncfile.Cdb_max
            self.Cdb_min = ncfile.Cdb_min
        except:
            self.Cdb_max = None
            self.Cdb_min = None

        try:
            self.cpp = ncfile.CPPS
        except:
            try:
                self.cpp = ncfile.__getattribute__('CPP-options')
            except:
                pass

        try:
            self.visc2 = ncfile.visc2       
        except:
            self.visc2 = None

        try:
            if 'NONLIN_EOS' not in self.cpp:
                self.Tcoef = ncfile.Tcoef
                self.Scoef = ncfile.Scoef
                self.R0 = ncfile.R0
        except:
            pass
        
        try:
            self.Zob = ncfile.Zob
        except:
            if output: print('no Zob in job ... using Zob = 0.01')
            self.Zob = 0.01
        
        try:
            self.VertCoordType = ncfile.VertCoordType
        except:
            self.VertCoordType = 'OLD'

        ncfile.close()

###################################################################################
#Load some data from the .in file not outputed
###################################################################################


    '''
    if (visctype.eq.0) then
      Old version with no east sponge
      CALL visc3d_GP  (Lm,Mm,N,u,v,z_r,z_w
            
    elseif (visctype.eq.1) then
      New version with  east sponge + new sponge sixes (1/12)
      CALL visc3d_S  (Lm,Mm,N,u,v,z_r,z_w

    elseif (visctype.eq.2) then
      New version with  east sponge + smaller sponge sixes (1/20)


    #############################################################

    


    '''


    def getin(self):

        if self.simul in ['cuc']:
            self.v_sponge=5.
            self.visctype=0
            self.lmdkpp=0
            self.Ricr=0.45
        elif self.simul in ['filam','shing']:
            self.v_sponge=3.
            self.visctype=0
            self.lmdkpp=0
            self.Ricr=0.45   
        elif self.simul in ['atlbig']:
            self.v_sponge=40.
            self.visctype=0
            self.lmdkpp=0
            self.Ricr=0.45   
        elif self.simul in ['chabu','nwat']:
            self.v_sponge=20.
            self.visctype=0
            self.lmdkpp=1
            self.Ricr=0.15
        elif self.simul in ['chabz']:
            self.v_sponge=20.
            self.visctype=1
            self.lmdkpp=1
            self.Ricr=0.15
        elif self.simul in ['chabu_sasha']:
            # A completer
            self.v_sponge=10.
            self.visctype=1
            self.lmdkpp=1
            self.Ricr=0.15   
        elif self.simul in ['bahan']:
            self.v_sponge=10.
            self.visctype=2
            self.lmdkpp=1
            self.Ricr=0.15
        elif self.simul in ['bahar']:
            self.v_sponge=10.
            self.visctype=1
            self.lmdkpp=1
            self.Ricr=0.15
        elif self.simul in ['bahas','bahat']:
            self.v_sponge=10.
            self.visctype=1
            self.lmdkpp=2
            self.Ricr=0.45
        elif self.simul in ['baham']:
            self.v_sponge=10.
            self.visctype=1
            self.lmdkpp=2
            self.Ricr=0.45
        elif self.simul in ['bahav']:
            self.v_sponge=10.
            self.visctype=2
            self.lmdkpp=2
            self.Ricr=0.45
        elif self.simul in ['bahaw']:
            self.v_sponge=10.
            self.visctype=1
            self.lmdkpp=2
            self.Ricr=0.45
        elif self.simul in ['gulfz']:
            self.v_sponge=40.
            self.visctype=0
            self.lmdkpp=1
            self.Ricr=0.15        
        elif self.simul in ['nesea','neseb']:
            self.v_sponge=10.
            self.visctype=0
            self.lmdkpp=1
            self.Ricr=0.15
        elif self.simul[:5]=='capat':
            self.v_sponge=20.         
            self.visctype=1  
            self.lmdkpp=1                
            self.Ricr=0.15
        elif self.simul in ['seamountk','seamountk160']:
            self.v_sponge=50.         
            self.visctype=1  
            self.lmdkpp=1                
            self.Ricr=0.45
        elif self.simul[:5]=='bsose':
            self.v_sponge=20.         
            self.visctype=1  
            self.lmdkpp=1                
            self.Ricr=0.15           
            
###################################################################################
#Force fortran order
###################################################################################


    def Forder(self,var):

        return np.asfortranarray(var.T,dtype=self.floattype)



###################################################################################
# get ocean date from file and convert it to '%m/%d - %H:00'
###################################################################################

    def oceandate(self):

        #print self.ncfile, self.infiletime

        ncfile = Dataset(self.ncfile, 'r')

        if self.ncname.model in ['ucla','croco_lionel']:
            self.oceantime = int(np.array(ncfile.variables['ocean_time'][self.infiletime]))
        else:
            try:
                self.oceantime = int(np.array(ncfile.variables['scrum_time'][self.infiletime]))
            except:
                try:
                    self.oceantime = int(np.array(ncfile.variables['time'][self.infiletime]))
                except:
                    self.oceantime = int(np.array(ncfile.variables['time_centered'][self.infiletime]))


        if not self.ncname.realyear:
        
            self.year = int(np.floor(self.oceantime/(360*24*3600)))
            
            self.oceantime = self.oceantime%(360*24*3600)

            self.month = self.oceantime//(24*3600)//30+1
            month_name = ["None","Jan","Feb","Mar","Apr", "May", "Jun", "Jul","Aug","Sep","Oct","Nov","Dec"] 

            self.day = self.oceantime//(24*3600) - (self.month-1) * 30 + 1

            self.hour = self.oceantime%(24*3600)//3600

            self.min = self.oceantime%(3600)//60

            self.date = month_name[self.month] + ' ' + '{0:02}'.format(self.day) + ' - ' +'{0:02}'.format(self.hour) + ':' + '{0:02}'.format(self.min) 

        else:
        
            date = self.ncname.realyear_origin + timedelta(days=self.oceantime/3600./24.)
            
            self.year = date.year
            self.month = date.month
            month_name = ["None","Jan","Feb","Mar","Apr", "May", "Jun", "Jul","Aug","Sep","Oct","Nov","Dec"] 

            self.day = date.day
            self.hour = date.hour
            self.min = date.minute

            self.date = month_name[self.month] + ' ' + '{0:02}'.format(self.day) + ' - ' +'{0:02}'.format(self.hour) + ':' + '{0:02}'.format(self.min) 

            # Could also do
            #self.date = tm.strftime("%B %d - %H:00", tm.gmtime(self.oceantime))
            
        ncfile.close()


###################################################################################
# get ocean date from file and convert it to '%m/%d - %H:00'
###################################################################################

    def dt(self, output=True):

        ncfile = Dataset(self.ncfile, 'r')
        if output: print('dt is read in ',self.ncfile)
        try:
            if self.ncname.model in ['ucla','croco_lionel']:
                self.dt = np.array(ncfile.variables['ocean_time'][1]) \
                        - np.array(ncfile.variables['ocean_time'][0])
            elif  self.ncname.model in ['croco','agrif_jc']:
                self.dt = np.array(ncfile.variables['scrum_time'][1]) \
                        - np.array(ncfile.variables['scrum_time'][0])
            elif  self.ncname.model in ['croco_xios']:
                self.dt = np.array(ncfile.variables['time_counter'][1]) \
                        - np.array(ncfile.variables['time_counter'][0])
            else:
                self.dt = np.array(ncfile.variables['time'][1]) \
                        - np.array(ncfile.variables['time'][0])
        except:
            if self.simul[:4]=='natl': self.dt = 24. * 3600
            elif self.simul =='atlbig_mean2': self.dt = 24. * 3600 * 5.
            else: self.dt = 0
        
        # Just a check
        if 'hourly' in self.simul or 'surf' in self.simul: 
            print(self.ncname.model)
            print('in dt: ',self.dt)      

        ncfile.close()


###################################################################################
#END CLASS LOAD
###################################################################################

###################################################################################
# FILES 
###################################################################################

"""
Path for ROMS simulations files.

For each simulation 
"""

###################################################################################

import sys, os

###################################################################################
class files(object):
    """
    Specify files location, file format and naming convention

    naming convention depends mainly on the model you use among

    

        
    """
###################################################################################


    def __init__(self, simul, time=0, output =True):



        #ROMSSISMS is defined in the ~/bashrc file
        ROMSSIMS=os.environ.get("ROMSSIMS")
        ROMSSIMSC=os.environ.get("ROMSSIMSC")
        ROMSSIMSK=os.environ.get("ROMSSIMSK")

        ###################################################
        # UCLA computers       
        ###################################################

        celtic='/mnt/celtic'
        avatar='/mnt/avatar'   
        shirma='/mnt/shirma'
        cherokee='/mnt/cherokee'
        comanche='/mnt/comanche'
        cheyenne='/mnt/cheyenne'
        kamiya='/mnt/kamiya'
        stratus='/mnt/stratus'
        goth='/mnt/goth'
        hun='/mnt/hun'

        if os.getenv('HOSTNAME')=='shirma.atmos.ucla.edu':
            shirma='/shirma'
        elif  os.getenv('HOSTNAME')=='celtic.atmos.ucla.edu':
            celtic='/celtic'
        elif  os.getenv('HOSTNAME')=='avatar.atmos.ucla.edu':
            avatar='/avatar'
        elif  os.getenv('HOSTNAME')=='cherokee.atmos.ucla.edu':
            cherokee='/cherokee' 
        elif  os.getenv('HOSTNAME')=='goth.atmos.ucla.edu':
            goth='/goth' 
        elif  os.getenv('HOSTNAME')=='hun.atmos.ucla.edu':
            hun='/hun'
        
        #Temporary fix while shirma unavailable
        shirma = cherokee



        ###################################################
        # UBO computers       
        ###################################################

        libra='/net/libra/local/tmp/1/'
        capella='/net/capella/local/tmp/2/'
        #print 'where is libra', libra

        krypton='/data0/project/meddle/'

        ###################################################

        if ROMSSIMS==None: 
            #ROMSSIMS= shirma + '/gula/ROMS/Simulations'
            ROMSSIMS= libra + '/gula/ROMS/Simulations'
        if ROMSSIMSC==None: 
            ROMSSIMSC= celtic + '/gula/ROMS/Simulations'
        if ROMSSIMSK==None:
            ROMSSIMSK= krypton + '/gula/ROMS/Simulations'


        ###################################################
        # Curie      
        ###################################################

        #if os.getenv('HOSTNAME')=='curie91':
        #    libra = '/ccc/store/cont003/gen7638/vicc/



        ################## 

        # default values (might be changed in simulation definition)
        self.model = 'ucla'
        self.fileformat='.nc'
        self.digits = 4
        self.realyear = False # False means using 360 days and 30 days per month

        ##################

        if simul=='gulfs':
            #self.his=ROMSSIMS+'/GULFS/HIS/gulfs_his.'
            self.his=avatar + '/nmolem/batavia/GULFS/HIS/gulfs_his.'
            self.grd=ROMSSIMS+'/GULFS/gulfs_grd.nc'  
            self.grd_corrected=ROMSSIMS+'/GULFS/gulfs_grd_corrected.nc'              
            self.Z0=avatar + '/nmolem/GULFS/HIS/Z0/gulfs_his_z0.'
            self.frc=ROMSSIMS+'/GULFS/gulfs_frc.nc'
            self.wind=ROMSSIMS+'/GULFS/gulfs_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=390

        elif simul=='gulfs_seas':
            self.his=ROMSSIMS+'/GULFS/gulfs_seas.'
            self.grd=ROMSSIMS+'/GULFS/gulfs_grd.nc'  
            self.Z0=avatar + '/nmolem/GULFS/HIS/Z0/gulfs_his_z0.'
            self.frc=ROMSSIMS+'/GULFS/gulfs_frc.nc'
            self.wind=ROMSSIMS+'/GULFS/gulfs_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=390

        elif simul=='hatt':
            self.his=avatar + '/nmolem/HATT/HIS/hatt_his.'
            #self.grd=avatar + '/nmolem/HATT/hatt_grd.nc'  
            self.grd=ROMSSIMS+'/HATT/hatt_grd.nc'  
            self.Z0=avatar + '/nmolem/HATT/HIS/Z0/hatt_his_z0.'
            self.frc=avatar + '/nmolem/HATT/hatt_frc.nc'
            self.wind=avatar + '/nmolem/HATT/hatt_frc.nc'
            self.tfile=2
            self.tstart=0
            self.tend=100

        elif simul=='hatt2':
            #self.his='/mnt/batavia/nmolem/HATT2/HIS/hatt2_his.'
            self.his='/mnt/avatar/nmolem/HATT2/HIS/hatt2_his.'            
            self.grd=ROMSSIMS+'/HATT2/hatt2_grd.nc'  
            self.grd_corrected=ROMSSIMS+'/HATT2/hatt2_grd.nc'
            if output: print(ROMSSIMS+'/HATT2/hatt2_grd.nc')
            self.Z0=ROMSSIMS+'/HATT2/Z/hatt2_his_z0.'
            self.frc=ROMSSIMS+'/HATT2/hatt2_frc.nc'
            self.wind=ROMSSIMS+'/HATT2/hatt2_frc.nc'
            self.tfile=2
            self.tstart=0
            self.tend=100

        elif simul=='hatt2_92':
            self.his=ROMSSIMS+'/HATT2/HIS/hatt2_his.'
            self.grd=ROMSSIMS+'/HATT2/hatt2_grd.nc'  
            self.grd_corrected=ROMSSIMS+'/HATT2/hatt2_grd_corrected.nc'
            self.Z0=ROMSSIMS+'/HATT2/Z/hatt2_his_z0.'
            self.frc=ROMSSIMS+'/HATT2/hatt2_frc.nc'
            self.wind=ROMSSIMS+'/HATT2/hatt2_frc.nc'
            self.tfile=2
            self.tstart=0
            self.tend=100

        elif simul=='shing':
            self.his=ROMSSIMS+'/SHING/HIS/shing_his.'
            self.grd=ROMSSIMS+'/SHING/shing_grd.nc'
            self.grd_corrected=ROMSSIMS+'/SHING/shing_grd_corrected.nc'            
            self.Z0=ROMSSIMS+'/SHING/HIS/Z0/shing_his_z0.'
            self.frc=ROMSSIMS+'/SHING/shing_frc.nc'
            self.wind=ROMSSIMS+'/SHING/shing_frc.nc'
            self.tfile=2
            self.tstart=0
            self.tend=187
            
        elif simul=='reshing':
            folder='/mnt/inca/gula/ROMS/Simulations'
            self.his=folder+'/RESHING/HIS/reshing_his.'
            self.grd=ROMSSIMS+'/SHING/shing_grd.nc'
            self.grd_corrected=ROMSSIMS+'/SHING/shing_grd_corrected.nc'            
            self.Z0=ROMSSIMS+'/SHING/HIS/Z0/shing_his_z0.'
            self.frc=ROMSSIMS+'/SHING/shing_frc.nc'
            self.wind=ROMSSIMS+'/SHING/shing_frc.nc'
            self.tfile=2
            self.tstart=0
            self.tend=187     
            

        elif simul=='filam':
            #self.his=ROMSSIMS+'/FILAM/lightfilam_his.'
            self.his=ROMSSIMS+'/FILAM/HIS/filam_his.'
            self.grd=ROMSSIMS+'/SHING/shing_grd.nc'
            self.grd_corrected=ROMSSIMS+'/SHING/shing_grd_corrected.nc'            
            self.Z0=ROMSSIMS+'/SHING/HIS/Z0/shing_his_z0.'
            self.frc=ROMSSIMS+'/SHING/shing_frc.nc'
            self.wind=ROMSSIMS+'/SHING/shing_frc.nc'
            self.tfile=2
            self.tstart=114
            self.tend=474
            
        elif simul=='filam_avg':
            #self.his=ROMSSIMS+'/FILAM/lightfilam_his.'
            self.his=ROMSSIMS+'/FILAM/HIS/filam_avg.'
            self.grd=ROMSSIMS+'/SHING/shing_grd.nc'
            self.grd_corrected=ROMSSIMS+'/SHING/shing_grd_corrected.nc'
            self.Z0=ROMSSIMS+'/SHING/HIS/Z0/shing_his_z0.'
            self.frc=ROMSSIMS+'/SHING/shing_frc.nc'
            self.wind=ROMSSIMS+'/SHING/shing_frc.nc'
            self.tfile=2
            self.tstart=114
            self.tend=474
            
        elif simul=='dilam':
            folder=cherokee+'/gula/ROMS/Simulations'
            self.his=folder + '/DILAM/HIS/dilam_his.'
            self.grd=folder + '/DILAM/shing_grd.nc'
            self.grd_corrected=ROMSSIMS+'/SHING/shing_grd_corrected.nc'
            self.Z0=ROMSSIMS+'/SHING/HIS/Z0/shing_his_z0.'
            self.frc=folder + '/DILAM/shing_frc.nc'
            self.wind=folder + '/DILAM/shing_frc.nc'
            self.tfile=2
            self.tstart=70
            self.tend=474
            
        elif simul=='dilam_avg':
            folder=cherokee+'/gula/ROMS/Simulations'
            self.his=folder + '/DILAM/HIS/dilam_avg.'
            self.grd=ROMSSIMS+'/SHING/shing_grd.nc'
            self.grd_corrected=ROMSSIMS+'/SHING/shing_grd_corrected.nc'
            self.Z0=ROMSSIMS+'/SHING/HIS/Z0/shing_his_z0.'
            self.frc=ROMSSIMS+'/SHING/shing_frc.nc'
            self.wind=ROMSSIMS+'/SHING/shing_frc.nc'
            self.tfile=2
            self.tstart=70
            self.tend=474     

        elif simul=='slope':
            folder=kamiya+'/gula/ROMS/Simulations'
            self.his=folder +'/SLOPE/HIS/slope_avg.'
            self.grd=folder +'/SLOPE/shing_grd.nc'
            self.grd_corrected=ROMSSIMS+'/SHING/shing_grd_corrected.nc'
            self.frc=folder +'/SLOPE/shing_frc.nc'
            self.wind=folder +'/SLOPE/shing_frc.nc'
            self.tfile=2
            self.tstart=60
            self.tend=349
            
            '''
            # TO MAKE dilam_avg:
            ncra -v zeta,ubar,vbar,u,v,temp,salt dilam_his.0248.nc dilam_his.0250.nc dilam_his.0252.nc dilam_avg.0250.nc
            ncra -v ocean_time,omega,AKv,AKt,hbls dilam_his.0248.nc dilam_his.0250.nc dilam_his.0252.nc dilam_avg.0250.a.nc
            ncks -A dilam_avg.0250.a.nc dilam_avg.0250.nc
            '''
        ################## FAMILY 2

        elif simul=='atlbig':
            if output: print('simul',simul)
            '''
            folder='/mnt/inca/gula'
            #folder='/shirma/gula/ROMS/Simulations'
            
            if (time<90) or (time>=204):
                self.his= celtic + '/gula/ATLBIG/DHIS/atlbig_dhis.'
            elif (time>=90) and (time<174):
                self.his=folder + '/ATLBIG/DHIS/atlbig_dhis.'
            else:
                self.his='/mnt/kamiya/gula/ATLBIG/atlbig_dhis.'
            '''

            #folder=avatar + '/nmolem'
            folder=ROMSSIMSK
            self.his=folder + '/ATLBIG/DHIS/atlbig_dhis.'
            self.grd=folder + '/ATLBIG/atlbig_grd.nc'  
            self.frc=folder + '/ATLBIG/atlbig_frc.nc'
            self.wind=folder + '/ATLBIG/atlbig_frc_wind.nc'
            self.tfile=6
            self.tstart=0
            self.tend=300

        elif simul=='atlbig_seas':
            #folder='/mnt/inca/gula'
            folder=shirma + '/gula/ROMS/Simulations'
            self.his=folder + '/ATLBIG/DHIS/atlbig_seas.'
            self.grd=folder + '/ATLBIG/atlbig_grd.nc'  
            self.frc=folder + '/ATLBIG/atlbig_frc.nc'
            self.wind=folder + '/ATLBIG/atlbig_frc_wind.nc'
            self.tfile=6
            self.tstart=0
            self.tend=200
            
        # atlbig_seas_2 is atlbig_seas minus the first year of simulation
        elif simul=='atlbig_seas_2':
            #folder='/mnt/inca/gula'
            folder=shirma + '/gula/ROMS/Simulations'
            self.his=folder + '/ATLBIG/DHIS/atlbig_seas_2.'
            self.grd=folder + '/ATLBIG/atlbig_grd.nc'  
            self.frc=folder + '/ATLBIG/atlbig_frc.nc'
            self.wind=folder + '/ATLBIG/atlbig_frc_wind.nc'
            self.tfile=6
            self.tstart=0
            self.tend=200
            
        # atlbig_seas_3 is only the 2nd year of simulation (from 0090 to 0152)
        elif simul=='atlbig_seas_3':
            #folder='/mnt/inca/gula'
            folder=shirma + '/gula/ROMS/Simulations'
            self.his=folder + '/ATLBIG/DHIS/atlbig_seas_3.'
            self.grd=folder + '/ATLBIG/atlbig_grd.nc'  
            self.frc=folder + '/ATLBIG/atlbig_frc.nc'
            self.wind=folder + '/ATLBIG/atlbig_frc_wind.nc'
            self.tfile=6
            self.tstart=0
            self.tend=200

        # atlbig_year is exactly the same year than gulfz_year (from 0102 to 0173)
        elif simul=='atlbig_year':
            #folder='/mnt/inca/gula'
            folder=shirma + '/gula/ROMS/Simulations'
            self.his=folder + '/ATLBIG/DHIS/atlbig_year.'
            self.grd=folder + '/ATLBIG/atlbig_grd.nc'  
            self.frc=folder + '/ATLBIG/atlbig_frc.nc'
            self.wind=folder + '/ATLBIG/atlbig_frc_wind.nc'
            self.tfile=6
            self.tstart=0
            self.tend=200

        # atlbig_gulfz is the total gulfz lenght - equivalent to gulfz_seas (from 0096 to 0185)
        #clim_seas atlbig_dhis.0096.nc atlbig_dhis.01[0-7]* atlbig_dhis.0180.nc
        elif simul=='atlbig_gulfz':
            #folder='/mnt/inca/gula'
            folder=shirma + '/gula/ROMS/Simulations'
            self.his=folder + '/ATLBIG/DHIS/atlbig_gulfz.'
            self.grd=folder + '/ATLBIG/atlbig_grd.nc'  
            self.frc=folder + '/ATLBIG/atlbig_frc.nc'
            self.wind=folder + '/ATLBIG/atlbig_frc_wind.nc'
            self.tfile=6
            self.tstart=0
            self.tend=200
            
            
        # atlbig_nwat the total nwat lenght - (from 0060 to 0131)
        #clim_seas atlbig_dhis.00[6-9]* atlbig_dhis.01[0-2]*
        elif simul=='atlbig_nwat_wrong':
            #folder='/mnt/inca/gula'
            folder=shirma + '/gula/ROMS/Simulations'
            self.his=folder + '/ATLBIG/DHIS/atlbig_nwat.'
            self.grd=folder + '/ATLBIG/atlbig_grd.nc'  
            self.frc=folder + '/ATLBIG/atlbig_frc.nc'
            self.wind=folder + '/ATLBIG/atlbig_frc_wind.nc'
            self.tfile=6
            self.tstart=0
            self.tend=200     
            
        # atlbig_nwat the total nwat lenght - (from 0102 to 171)
        #clim_seas atlbig_dhis.01[0-6]*
        elif simul=='atlbig_nwat':
            #folder='/mnt/inca/gula'
            folder=shirma + '/gula/ROMS/Simulations'
            self.his=folder + '/ATLBIG/DHIS/atlbig_nwat.'
            self.grd=folder + '/ATLBIG/atlbig_grd.nc'  
            self.frc=folder + '/ATLBIG/atlbig_frc.nc'
            self.wind=folder + '/ATLBIG/atlbig_frc_wind.nc'
            self.tfile=6
            self.tstart=0
            self.tend=200     
            
        #ncra atlbig_dhis.0* atlbig_dhis.mean.nc
        elif simul=='atlbig_mean':
            #folder='/mnt/inca/gula'
            #folder=shirma + '/gula/ROMS/Simulations'
            folder=ROMSSIMS
            self.his=folder + '/ATLBIG/DHIS/atlbig_mean.'
            self.grd=folder + '/ATLBIG/atlbig_grd.nc'  
            self.frc=folder + '/ATLBIG/atlbig_frc.nc'
            self.wind=folder + '/ATLBIG/atlbig_frc_wind.nc'
            self.tfile=6
            self.tstart=0
            self.tend=300
            
        #ncra atlbig_dhis.00[7-9]* atlbig_dhis.0[1-2]* atlbig_dhis.mean2.nc
        elif simul=='atlbig_mean2':
            #folder='/mnt/inca/gula'
            #folder=shirma + '/gula/ROMS/Simulations'
            folder=ROMSSIMS
            self.his=folder + '/ATLBIG/DHIS/atlbig_mean2.'
            self.grd=folder + '/ATLBIG/atlbig_grd.nc'  
            self.frc=folder + '/ATLBIG/atlbig_frc.nc'
            self.wind=folder + '/ATLBIG/atlbig_frc_wind.nc'
            self.tfile=1
            self.tstart=0
            self.tend=300              
      
      
          
        elif simul=='nwat':
            folder=ROMSSIMS
            self.his=folder + '/NWAT/DHIS/nwat_dhis.'
            self.grd=folder + '/NWAT/nwat_grd.nc'  
            self.frc=folder + '/NWAT/nwat_frc.nc'
            self.wind=folder + '/NWAT/nwat_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=1229
              
      
      
        elif simul=='nwat_sfc':
            folder='/data0/project/meddle/gula/ROMS/Simulations'
            self.his=folder + '/NWAT/sfc/nwat_sfc.'
            folder=ROMSSIMS
            self.grd=folder + '/NWAT/nwat_grd.nc'  
            self.frc=folder + '/NWAT/nwat_frc.nc'
            self.wind=folder + '/NWAT/nwat_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=1229
            

        elif simul=='nwatml':
            folder=cherokee+'/gula/ROMS/Simulations'
            #folder=shirma+'/gula/ROMS/Simulations'           
            self.his=folder + '/NWATML/DHIS/nwatmlb_dhis.'
            self.grd=folder + '/NWAT/nwat_grd.nc'  
            self.frc=folder + '/NWAT/nwat_frc.nc'
            self.wind=folder + '/NWAT/nwat_frc_wind.nc'
            self.tfile=4
            self.tstart=1100
            self.tend=1237
            
        elif simul=='nwatml_avg':
            folder=cherokee+'/gula/ROMS/Simulations'
            #folder=shirma+'/gula/ROMS/Simulations'           
            self.his=folder + '/NWATML/DHIS/nwatmlb_avg.'
            self.grd=folder + '/NWAT/nwat_grd.nc'  
            self.frc=folder + '/NWAT/nwat_frc.nc'
            self.wind=folder + '/NWAT/nwat_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=9999    
            
        # nwat_seas is nwat minus first 100 time-steps (for equilibration)
        #clim_seas nwat_dhis.0[1-5]*
        elif simul=='nwat_seas':
            folder=cherokee+'/gula/ROMS/Simulations'
            self.his=folder +'/NWAT/DHIS/nwat_dhis.seas.'
            self.grd=ROMSSIMS+'/NWAT/nwat_grd.nc'  
            self.frc=ROMSSIMS+'/NWAT/nwat_frc.nc'
            self.wind=ROMSSIMS+'/NWAT/nwat_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=450    
            
        # nwat_chabu is nwat during chabu
        #clim_seas nwat_dhis.0635.nc nwat_dhis.06[4-9]* nwat_dhis.0700.nc
        elif simul=='nwat_chabu':
            self.his=ROMSSIMS+'/NWAT/DHIS/nwat_chabu.'
            self.grd=ROMSSIMS+'/NWAT/nwat_grd.nc'  
            self.frc=ROMSSIMS+'/NWAT/nwat_frc.nc'
            self.wind=ROMSSIMS+'/NWAT/nwat_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=450    
            
            
        elif simul=='chabu':
            folder=ROMSSIMSC
            self.his=folder + '/CHABU/HIS/chabu_his.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.diags_vrt=folder + '/CHABU/HIS/chabu_diags_vrt.'
            self.tfile=2
            self.tstart=0
            self.tend=758        

        elif simul=='chabz':
            folder=avatar+'/gula/ROMS/Simulations'
            self.his=folder + '/CHABZ/HIS/chabz_his.'
            folder=celtic+'/gula/ROMS/Simulations'
            self.grd=folder + '/CHABU/chabu_grd.nc'
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc.nc'
            self.tfile=2
            self.tstart=20
            self.tend=758

            
             
        elif simul=='chabu_sasha':
            folder=celtic+'/gula/ROMS/Simulations'
            self.his=cheyenne+'/gula/Chabu/chabu_his.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=4
            self.tstart=0
            self.tend=1988  
            
        #nohup clim_seas chabu_his.00[4-9]?.nc chabu_his.0[1-9]* chabu_his.1* &
        elif simul=='chabu_sasha_his_year':
            folder=celtic+'/gula/ROMS/Simulations'
            self.his=cheyenne+'/gula/Chabu/chabu_his.seas.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=1988              
            
        #nohup clim_seas chabu_his.0
        elif simul=='chabu_sasha_his_year2':
            folder=celtic+'/gula/ROMS/Simulations'
            self.his=cheyenne+'/gula/Chabu/chabu_his.seas.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=1988        
            
        elif simul=='chabu_avg':
            folder=celtic+'/gula/ROMS/Simulations'
            self.his=folder + '/CHABU/HIS/chabu_avg.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=2
            self.tstart=0
            self.tend=759                 
            
            
        elif simul=='chabu_his_year':
            folder=ROMSSIMSC
            self.his=folder + '/CHABU/HIS/chabu_his.yearseas.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=760           

            '''
            ncra -d time,38,758 -v zeta,ubar,vbar,u,v,temp,salt chabu_avg.0*.nc chabu_avg.year.0000.nc
            ncra -d time,38,758 -v ocean_time,rho,w,hbls,hbbls chabu_avg.0*.nc chabu_avg.year.0000.a.nc
            ncks -A chabu_avg.year.0000.a.nc chabu_avg.year.0000.nc
            
            ncra -d time,38,758 -v zeta,ubar,vbar,u,v,temp,salt chabu_his.0*.nc chabu_his.year.0000.nc
            ncra -d time,38,758 -v ocean_time,omega,AKv,AKt,hbls chabu_his.0*.nc chabu_his.year.0000.a.nc
            ncks -A chabu_his.year.0000.a.nc chabu_his.year.0000.nc           
 
            
            ncra -d time,38,758 -v ocean_time chabu_avg.0*.nc chabu_avg.year.0000.nc
            
            for VAR in zeta ubar vbar hbls u v temp salt w omega AKv AKt; do
            echo $VAR
            ncra -d time,38,758 -v $VAR chabu_avg.0*.nc chabu_avg.year.0000.$VAR.nc
            done
            
            for VAR in zeta ubar vbar hbls u v temp salt w omega AKv AKt; do
            echo $VAR
            ncks -A chabu_avg.year.0000.$VAR.nc chabu_avg.year.0000.nc
            done

            clim_seas chabu_his.0038.nc chabu_his.00[4-9]?.nc chabu_his.0[1-7]??.nc
            clim_seas chabu_avg.0038.nc chabu_avg.00[4-9]?.nc chabu_avg.0[1-7]??.nc
            
            
            '''

        elif simul=='chabu_avg_year':
            folder=ROMSSIMSC
            self.his=folder + '/CHABU/HIS/chabu_avg.yearseas.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=760    
 
        elif simul=='chabu_his_meantest':
            '''
            ncra chabu_his.0138.nc chabu_his.01[4-9]*.nc chabu_his.02[0-1]* chabu_his.mean.0138-0218.nc
            ussed to test if mean is reversed duting instability events (about 40 days to mimic glideers in Aghulas)
            '''
            folder=ROMSSIMSC
            self.his=folder + '/CHABU/HIS/chabu_his.mean.0138-0218.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=760


        elif simul=='chabu_sep':
            '''
            Mean for a separated GS at the Bump
            ncra chabu_avg.052* chabu_avg.sep.0000.nc
            
            '''
            folder=celtic+'/gula/ROMS/Simulations'
            self.his=folder + '/CHABU/HIS/chabu_avg.sep.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=760
            
        elif simul=='chabu_nsep':
            '''
            Mean for a separated GS at the Bump
            ncra chabu_avg.057* chabu_avg.nsep.0000.nc
            ''' 
            folder=celtic+'/gula/ROMS/Simulations'
            self.his=folder + '/CHABU/HIS/chabu_avg.nsep.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=760    
            
        elif simul=='chabu_sepbis':
            '''
            Mean for a separated GS at the Bump
            ncra chabu_avg.052[0-4]* chabu_avg.sepbis.0000.nc
            
            '''
            folder=celtic+'/gula/ROMS/Simulations'
            self.his=folder + '/CHABU/HIS/chabu_avg.sepbis.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=760
            
        elif simul=='chabu_nsepbis':
            '''
            Mean for a separated GS at the Bump
            ncra chabu_avg.0570-4]* chabu_avg.nsepbis.0000.nc
            ''' 
            folder=celtic+'/gula/ROMS/Simulations'
            self.his=folder + '/CHABU/HIS/chabu_avg.nsepbis.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=760                
            
            
        elif simul=='chabu_sephis':
            '''
            Mean for a separated GS at the Bump
            ncra chabu_avg.052[0-4]* chabu_avg.sepbis.0000.nc
            
            '''
            folder=celtic+'/gula/ROMS/Simulations'
            self.his=folder + '/CHABU/HIS/chabu_his.sepbis.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=760
            
        elif simul=='chabu_nsephis':
            '''
            Mean for a separated GS at the Bump
            ncra chabu_avg.0570-4]* chabu_avg.nsepbis.0000.nc
            ''' 
            folder=celtic+'/gula/ROMS/Simulations'
            self.his=folder + '/CHABU/HIS/chabu_his.nsepbis.'
            self.grd=folder + '/CHABU/chabu_grd.nc'  
            self.frc=folder + '/CHABU/chabu_frc.nc'
            self.wind=folder + '/CHABU/chabu_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=760

     #################################################################
        #Same than BAHAM but with smoother topography 
        
        elif simul=='bahan':
            
            '''
            Same than original BAHAM but with smoother topography
            ''' 
            
            folder=kamiya+'/gula/ROMS/Simulations'
            self.his=folder + '/BAHAN/HIS/bahan_his.'
            self.grd=folder + '/BAHAN/bahan_grd.nc'  
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=2
            self.tstart=0
            self.tend=1000                
                       
        elif simul=='bahar':

            '''
            Same than BAHAN but with Rimix ON
            '''

            folder=kamiya+'/gula/ROMS/Simulations'
            self.his=folder + '/BAHAN/HIS/bahar_his.'
            self.grd=folder + '/BAHAN/bahan_grd.nc'
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=2
            self.tstart=0
            self.tend=1000       

        elif simul=='meanbahar':

            '''
            baham_his.007[2-8]* baham_his.00[8-9]* baham_his.01* baham_his.02[0-1]* baham_his.022[0-2]*
            '''

            folder=kamiya+'/gula/ROMS/Simulations'
            self.his=folder + '/BAHAN/HIS/baham_his.meanbahar.'
            self.grd=folder + '/BAHAN/bahan_grd.nc'
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=1
            self.tstart=0
            self.tend=10000


        elif simul=='bahas':

            '''
            Same than BAHAR but with MLCONVEC On
            '''

            folder=kamiya+'/gula/ROMS/Simulations'
            self.his=folder + '/BAHAN/HIS/bahas_his.'
            self.grd=folder + '/BAHAN/bahan_grd.nc'
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=2
            self.tstart=0
            self.tend=1000

        elif simul=='meanbahas':

            '''
            baham_his.022[6-9]* baham_his.02[3-9]* baham_his.03*
            '''

            folder=kamiya+'/gula/ROMS/Simulations'
            self.his=folder + '/BAHAN/HIS/baham_his.meanbahas.'
            self.grd=folder + '/BAHAN/bahan_grd.nc'
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=1
            self.tstart=0
            self.tend=10000

 
        elif simul=='baham':

            '''
            Wrapper of BAHAN, BAHAR and BAHAS
            '''

            folder=kamiya+'/gula/ROMS/Simulations'
            self.his=folder + '/BAHAN/HIS/baham_his.'
            self.grd=folder + '/BAHAN/bahan_grd.nc'
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=2
            self.tstart=0
            self.tend=312



        elif simul=='meanbaham':

            '''
            baham_his.007[2-8]* baham_his.00[8-9]* baham_his.0[1-3]*
            '''

            folder=kamiya+'/gula/ROMS/Simulations'
            self.his=folder + '/BAHAN/HIS/baham_his.meanbaham.'
            self.grd=folder + '/BAHAN/bahan_grd.nc'
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=1
            self.tstart=0
            self.tend=10000



        elif simul=='bahat':

            '''
            Same than BAHAR but with LMD_CONVEC Off
            '''

            folder=kamiya+'/gula/ROMS/Simulations'
            self.his=folder + '/BAHAN/HIS/bahat_his.'
            self.grd=folder + '/BAHAN/bahan_grd.nc'
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=2
            self.tstart=70
            self.tend=1000



        elif simul=='meanbahat':

            '''
            bahat_his.007[2-8]* bahat_his.00[8-9]* bahat_his.0[1-3]*
            '''

            folder=kamiya+'/gula/ROMS/Simulations'
            self.his=folder + '/BAHAN/HIS/bahat_his.meanbahat.'
            self.grd=folder + '/BAHAN/bahan_grd.nc'
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=1
            self.tstart=0
            self.tend=10000

        elif simul=='bahav':

            '''
            Different time than BAHAM series
            '''

            folder=ROMSSIMSK
            self.his=folder + '/BAHAN/HIS/bahav_his.'
            self.grd=folder + '/BAHAN/bahan_grd.nc'
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=2
            self.tstart=70
            self.tend=1000


        elif simul=='meanbahav1':


            folder=ROMSSIMSK
            self.his=folder + '/BAHAN/HIS/bahav_his.meanbahav1.'
            self.grd=folder + '/BAHAN/bahan_grd.nc'
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=1
            self.tstart=0
            self.tend=10000

        elif simul=='meanbahav2':


            folder=ROMSSIMSK
            self.his=folder + '/BAHAN/HIS/bahav_his.meanbahav2.'
            self.grd=folder + '/BAHAN/bahan_grd.nc'
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=1
            self.tstart=0
            self.tend=10000

        elif simul=='bahaw':

            '''
            rerun of BAHAM series with bahav opitons
            '''

            folder=ROMSSIMSK
            self.his=folder + '/BAHAN/HIS/bahaw_his.'
            self.grd=folder + '/BAHAN/bahan_anglecorrected_grd.nc'
            self.frc=folder + '/BAHAN/bahan_frc.nc'
            self.wind=folder + '/BAHAN/bahan_frc_wind.nc'
            self.tfile=2
            self.tstart=70
            self.tend=1000




        elif simul=='lucky0':

            '''
            nest of AZORE at 500m res (old version)
            '''

            folder=ROMSSIMS
            self.his=folder + '/LUCKY0/HIS/lucky_his.'
            self.grd=folder + '/LUCKY0/lucky_grd.nc'
            self.frc=folder + '/LUCKY0/lucky_frc.nc'
            self.wind=folder + '/LUCKY0/lucky_frc_wind.nc'
            self.tfile=2
            self.tstart=0
            self.tend=1000

        elif simul=='azore50':


            folder=ROMSSIMS
            self.his=folder + '/AZORE50/HIS/azore50_his.'
            self.grd=folder + '/AZORE50/azore_grd.nc'
            self.frc=folder + '/AZORE50/azore50_frc.nc'
            self.wind=folder + '/AZORE/azore_frc_wind.nc'
            self.tfile=2
            self.tstart=0
            self.tend=1000


        elif simul=='abaco':

            '''
            nest of NWAT at 750m res
            '''

            folder=ROMSSIMS
            self.his=folder + '/ABACO/HIS/abaco_his.'
            self.grd=folder + '/ABACO/abaco_grd.nc'
            self.frc=folder + '/ABACO/abaco_frc.nc'
            self.wind=folder + '/ABACO/abaco_frc_wind.nc'
            self.tfile=2
            self.tstat=44
            self.tend=10000


        elif simul=='gabaco':

            '''
            larger domain than abaco + code croco
            '''

            folder=ROMSSIMS
            self.model = 'croco'
            self.digits = 5
            self.his=folder + '/GABACO/HIS/gabaco_his.'
            self.grd=folder + '/GABACO/gabaco_grd.nc'
            self.frc=folder + '/GABACO/gabaco_frc.nc'
            self.wind=folder + '/GABACO/gabaco_frc_wind.nc'
            self.tfile=120
            self.tstat=0
            self.tend=100000


        elif simul=='rabaco':

            '''
            rerun of gabaco with different vertical stretching
            '''

            folder=ROMSSIMS
            self.model = 'croco'
            self.digits = 5
            self.his=folder + '/RABACO/HIS/rabaco_his.'
            self.grd=folder + '/RABACO/rabaco_grd.nc'
            self.frc=folder + '/RABACO/rabaco_frc.nc'
            self.wind=folder + '/RABACO/rabaco_frc_wind.nc'
            self.tfile=10
            self.tstat=0
            self.tend=10000


#################################################################

        ## nwat_year 
        #elif simul=='nwat_dhis_year':
            #self.his=ROMSSIMS+'/NWAT/DHIS/nwat_dhis.year.'
            #self.grd=ROMSSIMS+'/NWAT/nwat_grd.nc'  
            #self.frc=ROMSSIMS+'/NWAT/nwat_frc.nc'
            #self.wind=ROMSSIMS+'/NWAT/nwat_frc_wind.nc'
            #self.tfile=1
            #self.tstart=0
            #self.tend=9999

           # nwat_year 
        elif simul=='nwat_dhis_year':
            folder=cherokee+'/gula/ROMS/Simulations'
            self.his=folder + '/NWAT/DHIS/nwat_dhis.year.'
            self.grd=folder + '/NWAT/nwat_grd.nc'  
            self.frc=folder + '/NWAT/nwat_frc.nc'
            self.wind=folder + '/NWAT/nwat_frc_wind.nc'
            self.tfile=1
            self.tstart=0
            self.tend=9999
            
        # nwat_year 
        elif simul=='nwat_avg_year':
            self.his=ROMSSIMS+'/NWAT/DHIS/nwat_avg.year.'
            self.grd=ROMSSIMS+'/NWAT/nwat_grd.nc'  
            self.frc=ROMSSIMS+'/NWAT/nwat_frc.nc'
            self.wind=ROMSSIMS+'/NWAT/nwat_frc_wind.nc'
            self.tfile=1
            self.tstart=0
            self.tend=10   
            
        # nwat_year 
        elif simul=='nwat_diags_eddy_avg_year':
            self.his=ROMSSIMS+'/NWAT/DHIS/nwat_diags_eddy_avg.year.'
            self.grd=ROMSSIMS+'/NWAT/nwat_grd.nc'  
            self.frc=ROMSSIMS+'/NWAT/nwat_frc.nc'
            self.wind=ROMSSIMS+'/NWAT/nwat_frc_wind.nc'
            self.tfile=1
            self.tstart=0
            self.tend=10   
                      
            

#temporary################################################################

        elif simul=='gulfz':
            self.his=ROMSSIMS+'/GULFZ/HIS/gulfz_his.'
            self.grd=ROMSSIMS+'/GULFZ/gulfz_grd.nc'  
            self.frc=ROMSSIMS+'/GULFZ/gulfz_frc.nc'
            self.wind=ROMSSIMS+'/GULFZ/gulfz_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=446


        elif simul=='gulfz_surface':
            self.his=ROMSSIMS+'/GULFZ/uv_surface/gulfz_surf.'
            self.grd=ROMSSIMS+'/GULFZ/gulfz_grd.nc'
            self.frc=ROMSSIMS+'/GULFZ/gulfz_frc.nc'
            self.wind=ROMSSIMS+'/GULFZ/gulfz_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=446

        # gulfz seas is all gulfz (from 0000 to 0446) 
        # clim_seas gulfz_his.0*
        elif simul=='gulfz_seas':
            self.his=ROMSSIMS+'/GULFZ/gulfz_seas.'
            self.grd=ROMSSIMS+'/GULFZ/gulfz_grd.nc'  
            self.frc=ROMSSIMS+'/GULFZ/gulfz_frc.nc'
            self.wind=ROMSSIMS+'/GULFZ/gulfz_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=450

        # gulfz year is one year (from 0030 to 0389) 
        # clim_seas gulfz_his.00[3-9]* gulfz_his.0[1-2]* gulfz_his.03[0-8]*
        elif simul=='gulfz_year': 
            self.his=ROMSSIMS+'/GULFZ/gulfz_year.'
            self.grd=ROMSSIMS+'/GULFZ/gulfz_grd.nc'  
            self.frc=ROMSSIMS+'/GULFZ/gulfz_frc.nc'
            self.wind=ROMSSIMS+'/GULFZ/gulfz_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=450



        elif simul=='nesea':
            self.his=ROMSSIMS+'/NESEA/HIS/nesea_his.'
            self.grd=ROMSSIMS+'/NESEA/nesea_grd.nc'  
            self.frc=ROMSSIMS+'/NESEA/nesea_frc.nc'
            self.wind=ROMSSIMS+'/NESEA/nesea_frc_wind.nc'
            self.tfile=2
            self.tstart=0
            self.tend=195

        elif 'nesea_meander' in simul:
            self.model = 'croco'
            self.digits = 5
            if 'avg' in simul:
                self.his=ROMSSIMS+'/NESEA/HIS/nesea_avg.'
            else:
                self.his=ROMSSIMS+'/NESEA/HIS/nesea_his.'
            self.grd=ROMSSIMS+'/NESEA/nesea_grd.nc'
            self.frc=ROMSSIMS+'/NESEA/nesea_frc.nc'
            self.wind=ROMSSIMS+'/NESEA/nesea_frc_wind.nc'
            self.tfile=10
            self.tstart=0
            self.tend=500

        elif simul=='nesea_seas':
            self.his=ROMSSIMS+'/NESEA/HIS/nesea_his.mean.'
            self.grd=ROMSSIMS+'/NESEA/nesea_grd.nc'  
            self.frc=ROMSSIMS+'/NESEA/nesea_frc.nc'
            self.wind=ROMSSIMS+'/NESEA/nesea_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=5

        elif simul=='neseb':
            self.his=ROMSSIMS+'/NESEB/HIS/neseb_his.'
            self.grd=ROMSSIMS+'/NESEB/nesea_grd.nc'  
            self.frc=ROMSSIMS+'/NESEB/HIS/nesea_frc.nc'
            self.wind=ROMSSIMS+'/NESEB/HIS/nesea_frc_wind.nc'
            self.tfile=2
            self.tstart=80
            self.tend=500

        elif simul=='nefro':
            folder=celtic+'/gula/ROMS/Simulations'
            self.his=folder+'/NEFRO/HIS/nefro_his.'
            self.grd=folder+'/NEFRO/nefro_grd.nc'  
            self.frc=folder+'/NEFRO/nefro_frc.nc'
            self.wind=folder+'/NEFRO/nefro_frc_wind.nc'
            self.tfile=2
            self.tstart=0
            self.tend=80
            
            
            
            
            
            
            
            


        ################## Some others

        elif simul=='kuro':
            self.his='/mnt/batavia/KURO/HIS/kuro_his.'
            self.grd='/mnt/batavia/KURO/kuro_grd.nc'  
            self.frc='/mnt/batavia/KURO/kuro_frc.nc'
            self.wind='/mnt/batavia/KURO/kuro_frc.nc'
            self.tfile=2
            self.tstart=0
            self.tend=300

        elif simul=='capat':

            self.his=cherokee + '/gula/CAPAT/Run/capat_his.'
            self.grd=cherokee + '/gula/CAPAT/Run/capat_grd.nc'  
            self.frc=cherokee + '/gula/CAPAT/Run/capat_frc.nc'
            self.wind=cherokee + '/gula/CAPAT/Run/capat_frc_wind.nc'
            self.tfile=100
            self.tstart=0
            self.tend=450

        elif simul=='capat_avg':

            self.his=cherokee + '/gula/CAPAT/Run/capat_avg.'
            self.grd=cherokee + '/gula/CAPAT/Run/capat_grd.nc'
            self.frc=cherokee + '/gula/CAPAT/Run/capat_frc.nc'
            self.wind=cherokee + '/gula/CAPAT/Run/capat_frc_wind.nc'
            self.tfile=100
            self.tstart=0
            self.tend=450

        elif simul=='capat_ts':

            self.his=cherokee + '/gula/CAPAT/Run/capat_diags_ts.'
            self.grd=cherokee + '/gula/CAPAT/Run/capat_grd.nc'       
            self.frc=cherokee + '/gula/CAPAT/Run/capat_frc.nc'
            self.wind=cherokee + '/gula/CAPAT/Run/capat_frc_wind.nc'
            self.tfile=100
            self.tstart=0
            self.tend=450


        elif simul=='capat_ts_avg':

            self.his=cherokee + '/gula/CAPAT/Run/capat_diags_ts_avg.'
            self.grd=cherokee + '/gula/CAPAT/Run/capat_grd.nc'                        
            self.frc=cherokee + '/gula/CAPAT/Run/capat_frc.nc'
            self.wind=cherokee + '/gula/CAPAT/Run/capat_frc_wind.nc'
            self.tfile=100
            self.tstart=0
            self.tend=450


        ################## Some others
    
            
        elif simul[:5]=='capat':
            
            if simul[-6:]=='closed':
                bnd = 'closed'
            else:
                bnd = 'open'
                
            self.grd=shirma + '/gula/CAPAT_ENERGY/Run_'+ bnd +'/capat_grd.nc'  
            self.frc=shirma + '/gula/CAPAT_ENERGY/Run_'+ bnd +'/capat_frc.nc'
            self.wind=shirma + '/gula/CAPAT_ENERGY/Run_'+ bnd +'/capat_frc_wind.nc'
            self.tfile=100
            self.tstart=0
            self.tend=2000               

            if 'avg' in simul:
                if 'uv' in simul:
                    self.his=shirma + '/gula/CAPAT_ENERGY/Run_'+ bnd +'/capat_diags_uv_avg.'
                elif 'vrt' in simul:
                    self.his=shirma + '/gula/CAPAT_ENERGY/Run_'+ bnd +'/capat_diags_vrt_avg.'
                else:
                    self.his=shirma + '/gula/CAPAT_ENERGY/Run_'+ bnd +'/capat_avg.'
            else:
                if 'uv' in simul: 
                    self.his=shirma + '/gula/CAPAT_ENERGY/Run_'+ bnd +'/capat_diags_uv.'
                elif 'vrt' in simul:
                    self.his=shirma + '/gula/CAPAT_ENERGY/Run_'+ bnd +'/capat_diags_vrt.'
                else:   
                    self.his=shirma + '/gula/CAPAT_ENERGY/Run_'+ bnd +'/capat_his.'

            
         ################## LPO
           
            
        elif simul=='reykja':
            folder=comanche + '/gula'
            self.his=folder + '/REYKJA/HIS/reykja_his.'
            self.grd=folder + '/REYKJA/reykja_grd.nc'  
            self.frc=folder + '/REYKJA/reykja_frc.nc'
            self.wind=folder + '/REYKJA/reykja_frc_wind.nc'
            self.tfile=2
            self.tstart=0
            self.tend=9999 
            
            
        elif simul=='reykja_mean':
            '''
            
            for VAR in zeta ubar vbar hbls u v temp salt w omega AKv AKt rho; do
            echo $VAR
            ncra -v $VAR reykja_his.00[5-9]* reykja_his.0[1-9]* reykja_his.mean.$VAR.nc
            done
            
            for VAR in zeta uvbar hbls u v temp salt w omega AKv AKt rho; do    
            echo $VAR
            ncks -A reykja_his.mean.$VAR.nc reykja_his.mean.0000.nc
            done
            
            '''
            folder=comanche + '/gula'
            self.his=folder + '/REYKJA/HIS/reykja_his.mean.'
            self.grd=folder + '/REYKJA/reykja_grd.nc'  
            self.frc=folder + '/REYKJA/reykja_frc.nc'
            self.wind=folder + '/REYKJA/reykja_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=9999 
                      

        elif simul=='reykja_mean300':
            '''
            
            for VAR in zeta ubar vbar hbls u v temp salt w omega AKv AKt rho; do
            echo $VAR
            ncra -v $VAR reykja_his.00[5-9]* reykja_his.0[1-9]* reykja_his.mean.$VAR.nc
            done
            
            for VAR in zeta uvbar hbls u v temp salt w omega AKv AKt rho; do    
            echo $VAR
            ncks -A reykja_his.mean.$VAR.nc reykja_his.mean.0000.nc
            done
            
            '''
            folder= libra + '/gula/ROMS/Simulations/REYKJA/'
            self.his=folder + 'reykja_his.mean.'
            self.grd=folder + 'reykja_grd.nc'  
            self.frc=folder + 'reykja_frc.nc'
            self.wind=folder + 'reykja_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=9999 
                                
        elif simul=='reykjb':
            folder=celtic+'/gula/ROMS/Simulations'
            self.his=folder + '/REYKJB/HIS/reykjb_his.'
            self.grd=folder + '/REYKJB/reykjb_grd.nc'  
            self.frc=folder + '/REYKJB/reykjb_frc.nc'
            self.wind=folder + '/REYKJB/reykjb_frc_wind.nc'
            self.tfile=2
            self.tstart=0
            self.tend=182             
            
        ################## Some others
  
        elif simul=='msosea':
            folder='/mnt/comanche/lrenault/MSOSEA/'
            self.his=folder + 'HIS/msosea_his.'
            self.grd=folder + 'msosea_grd.nc'
            self.frc=folder + 'msosea_frc.nc'
            self.wind=folder + 'msosea_frc.nc'
            self.tfile=2
            self.tstart=0
            self.tend=210
        ################## Some others
  
        elif simul=='msosea_bug':
            folder='/mnt/comanche/lrenault/MSOSEA/'
            self.his=folder + 'HIS/Test_without_MLCONVEC_badbudg/msosea_his.'
            self.grd=folder + 'msosea_grd.nc'
            self.frc=folder + 'msosea_frc.nc'
            self.wind=folder + 'msosea_frc.nc'
            self.tfile=2
            self.tstart=0
            self.tend=14
            
        ################## Some others
                  
        elif simul=='papua':
            folder='/mnt/comanche/lrenault'
            self.his=folder + '/PAPUA/HIS/papua_his.'
            self.grd=folder + '/PAPUA/papua_grd.nc'  
            self.frc=folder + '/PAPUA/papua_frc.nc'
            self.wind=folder + '/PAPUA/papua_frc_wind.nc'
            self.tfile=6
            self.tstart=102
            self.tend=900
  
        elif simul=='sosea':
            folder='/mnt/comanche/lrenault'
            self.his=folder + '/SOSEA/HIS/sosea100_his_sp.'
            self.grd=folder + '/SOSEA/sosea_grdV4.nc'
            self.frc=folder + '/SOSEA/sosea_frc.nc'
            self.wind=folder + '/SOSEA/sosea_wind_frc.nc'
            self.tfile=2
            self.tstart=8
            self.tend=144

        elif simul=='louis_avg':
            folder='/mnt/comanche/lrenault/LOUIS'
            self.his=folder + '/AVG200/louis_avg.'
            self.grd=folder + '/louis_grd_sm.nc'
            self.frc=folder + '/louis_frc.nc'
            self.wind=folder + '/louis_wind_frc.nc'
            self.tfile=1
            self.tstart=1
            self.tend=197
            self.digits = 5

        elif simul=='solom':
            #folder='/mnt/balearic/lrenault/SOLOMON/'
            folder='/mnt/comanche/nmolem'           
            self.his=folder + '/SOLOM/HIS/solom_his.'
            self.grd=folder + '/SOLOM/solom_grd.nc'  
            self.frc=folder + '/SOLOM/solom_frc.nc'
            self.wind=folder + '/SOLOM/solom_frc_wind.nc'
            self.tfile=5
            self.tstart=10
            self.tend=130   
            
        elif simul=='bsosea':
            #folder='/mnt/balearic/lrenault/SOLOMON/'
            folder='/mnt/comanche/lrenault'           
            self.his=folder + '/BSOSEA/HISORI/bsosea_his_ori.'
            self.grd=folder + '/BSOSEA/bsosea_grd.nc'  
            self.frc=folder + '/BSOSEA/bsosea_frc.nc'
            self.wind=folder + '/BSOSEA/bsosea_wind_frc.nc'
            self.tfile=2
            self.tstart=10
            self.tend=46            
            
        elif simul=='bsosea50':
            folder='/mnt/comanche/lrenault'           
            self.his=folder + '/BSOSEA/HIS50/bsosea_his_50.'
            self.grd=folder + '/BSOSEA/bsosea_grd.nc'  
            self.frc=folder + '/BSOSEA/bsosea_frc.nc'
            self.wind=folder + '/BSOSEA/bsosea_wind_frc.nc'
            self.tfile=2
            self.tstart=10
            self.tend=46            
                            
            
        elif simul=='bsosea100':
            folder='/mnt/comanche/lrenault'           
            self.his=folder + '/BSOSEA/HIS100/bsosea_his_V100.'
            self.grd=folder + '/BSOSEA/bsosea_grd.nc'  
            self.frc=folder + '/BSOSEA/bsosea_frc.nc'
            self.wind=folder + '/BSOSEA/bsosea_wind_frc.nc'
            self.tfile=2
            self.tstart=10
            self.tend=46            
            
            
        elif simul=='swpac':
            folder='/mnt/comanche/nmolem'
            self.his=folder + '/SWPAC/HIS/swpac_his.'
            self.grd=folder + '/SWPAC/swpac_grd.nc'  
            self.frc=folder + '/SWPAC/swpac_frc.nc'
            self.wind=folder + '/SWPAC/swpac_frc_wind.nc'
            self.tfile=5
            self.tstart=20
            self.tend=437  

        elif simul=='solom_mean':
            folder='/mnt/balearic/lrenault/SOLOMON/'
            self.his=folder + '/SOLOM/HIS/solom_avg.'
            self.grd=folder + '/SOLOM/solom_grd.nc'  
            self.frc=folder + '/SOLOM/solom_frc.nc'
            self.wind=folder + '/SOLOM/solom_frc_wind.nc'
            self.tfile=1
            self.tstart=0
            self.tend=0 

        elif simul=='pacbig':
            #folder='/mnt/comanche/nmolem'
            folder='/data0/project/meddle/gula/ROMS/Simulations'
            self.his=folder + '/PACBIG/HIS/pacbig_his.'
            self.grd=folder + '/PACBIG/pacbig_grd_nw_dig.nc'  
            self.frc=folder + '/PACBIG/pacbig_frc.nc'
            self.wind=folder + '/PACBIG/pacbig_frc_wind.nc'
            self.tfile=6
            self.tstart=102
            self.tend=750      
 
        elif simul=='pacbig_seas':
            folder='/mnt/comanche/nmolem'
            self.his=folder + '/PACBIG/HIS/pacbig_seas.'
            self.grd=folder + '/PACBIG/pacbig_grd_nw_dig.nc'  
            self.frc=folder + '/PACBIG/pacbig_frc.nc'
            self.wind=folder + '/PACBIG/pacbig_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=5      
            
        elif simul=='nwatl':
            folder='/mnt/comanche/nmolem'
            self.his=folder + '/NWATL/HIS/nwatl_his.'
            self.grd=folder + '/NWATL/nwatl_grd.nc'  
            self.frc=folder + '/NWATL/nwatl_frc.nc'
            self.wind=folder + '/NWATL/nwatl_frc_wind.nc'
            self.tfile=6
            self.tstart=204
            self.tend=600       
        
        
        elif simul=='atlmed':
            #folder='/mnt/inca/gula'
            folder=shirma + '/gula/ROMS/Simulations'
            self.his=folder + '/ATLMED/HIS/atlmed_his.'
            self.grd=folder + '/ATLMED/atlmed_grd.nc'  
            self.frc=folder + '/ATLMED/atlmed_frc.nc'
            self.wind=folder + '/ATLMED/atlmed_frc_wind.nc'
            self.tfile=6
            self.tstart=0
            self.tend=200

        elif simul=='atlmed_seas':
            #folder='/mnt/inca/gula'
            folder=shirma + '/gula/ROMS/Simulations'
            self.his=folder + '/ATLMED/DHIS/atlmed_seas.'
            self.grd=folder + '/ATLMED/atlmed_grd.nc'  
            self.frc=folder + '/ATLMED/atlmed_frc.nc'
            self.wind=folder + '/ATLMED/atlmed_frc_wind.nc'
            self.tfile=5
            self.tstart=0
            self.tend=5

        #elif simul=='seamount':
            #self.his='/home/gula/Desktop/ROMS/code/seamount2temp/mount_his_res1kmu001z82.'
            #self.grd='/home/gula/Desktop/ROMS/code/seamount2temp/mount_his_res1kmu001z82.0000.nc'  
            #self.tfile=200
            #self.tstart=0
            #self.tend=200
            
        elif simul=='seamountk':
            folder=kamiya + '/gula/ROMS/Simulations'
            self.his=folder + '/SEAMOUNTK/seamount40.'
            self.grd=folder + '/SEAMOUNTK/seamount40.0000.nc'
            self.tfile=200
            self.tstart=0
            self.tend=200
            
        elif simul=='seamountk160':
            folder=kamiya + '/gula/ROMS/Simulations'
            self.his=folder + '/SEAMOUNTK/seamount160.'
            self.grd=folder + '/SEAMOUNTK/seamount160.0000.nc'
            self.tfile=200
            self.tstart=0
            self.tend=200
            
        elif simul=='roms':
            self.grd=libra + '/gula/ROMS/ToolsRoms_nc4/EGRID/roms_grd.nc'
            self.his='hehe'
            self.tfile=1
            self.tstart=0
            self.tend=0
            if output: print(self.grd)
        elif simul=='cuc':
            folder=avatar + '/nmolem'
            folder = libra + '/gula/ROMS/Simulations/'
            self.his=folder + '/CUC/RUN3/HIS/cuc_his.'
            self.grd=folder + '/CUC/cuc_grd.nc'  
            self.frc=folder + '/CUC/cuc_frc.nc'
            self.wind=folder + '/CUC/cuc_frc.nc'
            self.tfile=2
            self.tstart=0
            self.tend=999
            
        elif simul=='smcc':
            folder=avatar + '/nmolem'
            #folder = libra + '/gula/ROMS/Simulations/'
            self.his=folder + '/SMCC/RUN3/HIS/smcc_his.'
            self.grd=folder + '/SMCC/smcc_grd_mod.nc'  
            self.frc=folder + '/SMCC/smcc_frc_oro.nc'
            self.wind=folder + '/SMCC/smcc_frc_oro.nc'
            self.tfile=2
            self.tstart=0
            self.tend=812
            
        elif simul=='midcalL1':
            self.realyear = True
            self.realyear_origin = datetime(1994,1,1)
            
            #folder= '/mnt/inca/akanc'
            folder = libra + '/gula/ROMS/Simulations/'
            self.his=folder + '/MidCal/L1/04_out_1hr/usw1_avg.'
            self.grd=folder + '/MidCal/L1/04_out_1hr/usw1_grd.nc'
            self.frc=folder + '/MidCal/L1/02_input/usw1_wnd.nc'
            self.wind=folder + '/MidCal/L1/02_input/usw1_wnd.nc'
            self.tfile=24
            self.tstart=0
            self.tend=5000
            
        elif simul=='midcalL1_notides':
            self.realyear = True
            self.realyear_origin = datetime(1994,1,1)

            #folder_grd= '/mnt/inca/akanc'
            folder_grd = libra + '/gula/ROMS/Simulations/'

            folder= hun + '/roybarkan/MID_CAL_NO_Tides/'
            self.his=folder + 'usw1_avg.'
            self.grd= folder_grd + '/MidCal/L1/04_out_1hr/usw1_grd.nc'
            self.tfile=24
            self.tstart=0
            self.tend=5256

        elif simul=='midcalL1_notides_noniw':
            self.realyear = True
            self.realyear_origin = datetime(1994,1,1)

            #folder_grd= '/mnt/inca/akanc'
            folder_grd = libra + '/gula/ROMS/Simulations/'

            folder= hun + '/roybarkan/MID_CAL_NO_Tides_NO_NIWs/'
            self.his=folder + 'usw1_avg.'
            self.grd= folder_grd + '/MidCal/L1/04_out_1hr/usw1_grd.nc'
            self.tfile=24
            self.tstart=0
            self.tend=5256

        elif simul=='midcalL1_6h':
            self.realyear = True
            self.realyear_origin = datetime(1994,1,1)
            
            #folder= '/mnt/inca/akanc'
            folde = libra + '/gula/ROMS/Simulations/'
            self.his=folder + '/MidCal/L1/03_out_6hr/usw1_avg.'
            self.grd=folder + '/MidCal/L1/04_out_1hr/usw1_grd.nc'
            self.frc=folder + '/MidCal/L1/02_input/usw1_wnd.nc'
            self.wind=folder + '/MidCal/L1/02_input/usw1_wnd.nc'
            self.tfile=4
            self.tstart=0
            self.tend=1280
            
        elif simul=='midcalL2':
            self.realyear = True
            self.realyear_origin = datetime(1994,1,1)
            
            #folder= '/mnt/inca/akanc'
            folder = libra + '/gula/ROMS/Simulations/'

            self.his=folder + '/MidCal/L2/03_output/usw2_avg.'
            self.grd=folder + '/MidCal/L2/03_output/usw2_grd.nc'
            self.frc=folder + '/MidCal/L2/02_input/usw2_wnd.nc'
            self.wind=folder + '/MidCal/L2/02_input/usw2_wnd.nc'
            self.tfile=24
            self.tstart=384
            self.tend=5000

        elif simul=='midcalL3w':
            self.realyear = True
            self.realyear_origin = datetime(1994,1,1)
            
            folder= '/mnt/inca/akanc'
            folder = libra + '/gula/ROMS/Simulations/'
            self.his=folder + '/MidCal/L3/03_out_winter2006/usw3_avg.'
            self.grd=folder + '/MidCal/L3/02_input/usw3_grd.nc'
            self.frc=folder + '/MidCal/L3/02_input/usw2_wnd.nc'
            self.wind=folder + '/MidCal/L3/02_input/usw2_wnd.nc'
            self.tfile=24
            self.tstart=0
            self.tend=1000


        elif simul=='midcalL3w_daily':
            self.realyear = True
            self.realyear_origin = datetime(1994,1,1)
            
            folder= '/mnt/inca/gula/ROMS/Simulations'
            self.his=folder + '/MidCal/usw3_avg.'
            self.grd='/mnt/inca/akanc/MidCal/L3/02_input/usw3_grd.nc'
            self.frc='/mnt/inca/akanc/MidCal/L3/02_input/usw2_wnd.nc'
            self.wind='/mnt/inca/akanc/MidCal/L3/02_input/usw2_wnd.nc'
            self.tfile=1
            self.tstart=0
            self.tend=35


        elif simul=='midcalL3s':
            self.realyear = True
            self.realyear_origin = datetime(1994,1,1)
            
            folder= '/mnt/inca/akanc'
            folder = libra + '/gula/ROMS/Simulations/'
            self.his=folder + '/MidCal/L3/04_out_spring2007/usw3_avg.'
            self.grd=folder + '/MidCal/L3/02_input/usw3_grd.nc'
            self.frc=folder + '/MidCal/L3/02_input/usw2_wnd.nc'
            self.wind=folder + '/MidCal/L3/02_input/usw2_wnd.nc'
            self.tfile=24
            self.tstart=0
            self.tend=1000

        elif simul=='mnt_h1':
            folder=celtic + '/gula/MNT'
            self.his=folder + '/test_h/mount_f1_h1.'
            self.grd=folder + '/test_h/mount_grd.nc'  
            self.frc=folder + '/HIS_f0/mount_frc.nc'
            self.wind=folder + '/HIS_f0/mount_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=1000
            
        elif simul=='mnt_h01':
            folder=celtic + '/gula/MNT'
            self.his=folder + '/test_h/mount_f1_h01.'
            self.grd=folder + '/test_h/mount_grd_h0.1.nc'  
            self.frc=folder + '/HIS_f0/mount_frc.nc'
            self.wind=folder + '/HIS_f0/mount_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=1000
            
            
        elif simul=='mnt_h05':
            folder=celtic + '/gula/MNT'
            self.his=folder + '/test_h/mount_f1_h05.'
            self.grd=folder + '/test_h/mount_grd_h0.5.nc'  
            self.frc=folder + '/HIS_f0/mount_frc.nc'
            self.wind=folder + '/HIS_f0/mount_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=1000 
            
            
        elif simul=='ideal':
            folder=shirma + '/gula/ROMS/Simulations'
            self.his=folder + '/IDEAL/his.'
            self.grd=folder + '/IDEAL/his.0000.nc'  
            self.frc=folder + '/HIS_f0/mount_frc.nc'
            self.wind=folder + '/HIS_f0/mount_frc.nc'
            self.tfile=12
            self.tstart=0
            self.tend=1000


        ################## 
        # JC simulations

        elif simul=='aacc':
            self.model = 'agrif_jc'
            folder= libra + '/gula/ROMS/Simulations/AACC'
            self.his=folder + '/roms_his_Y1M3.'
            self.grd=folder + '/roms_his_Y1M3.0000.nc'
            self.frc=folder + '/roms_his_Y1M3.0000.nc'
            self.wind=folder + '/roms_his_Y1M3.0000.nc'
            self.tfile=10
            self.tstart=0
            self.tend=1000

            if output: print('AACC ', self.grd)


        elif simul=='aacc_8k':
            self.model = 'agrif_jc'
            folder= '/net/octant/local/tmp/1/collin/Gpla_8k'
            self.his=folder + '/roms_his_'
            self.grd=folder + '/roms_his_Y20M1.nc'
            self.frc=folder + '/roms_his_Y20M1.nc'
            self.wind=folder + '/roms_his_Y20M1.nc'
            self.tfile=10
            self.tstart=0
            self.Ystart=20
            self.Mstart=1
            self.tend=1000


        elif simul=='aacc_2000_8k':
            self.model = 'agrif_jc'
            folder=  '/net/octant/local/tmp/1/collin/HFreq_Pla_2000_8km'
            self.his=folder + '/roms_his_'
            self.grd=folder + '/roms_his_Y36M3.nc'
            self.frc=folder + '/roms_his_Y36M3.nc'
            self.wind=folder + '/roms_his_Y36M3.nc'
            self.tfile=50
            self.tstart=0
            self.Ystart=36
            self.Mstart=1

            self.tend=1200

        elif simul=='blkperio':
            self.model = 'agrif_jc'
            folder=  '/net/centaure/local/tmp/1/ragoasha/BUIC/Outputs_PhD/blkperio'
            self.fileformat='.nc_1'
            self.his=folder + '/roms_avg_'
            self.grd=folder + '/roms_avg_Y2006M1.nc_1'
            self.frc=folder + '/roms_avg_Y2006M1.nc_1'
            self.wind=folder + '/roms_avg_Y2006M1.nc_1'
            self.tfile=10
            self.tstart=0
            self.Ystart=2006
            self.Mstart=1

            self.tend=1200

        ##################
        
        elif simul=='natl_GS':
            '''read only first month for test, needs to increment year and month'''

            self.model = 'croco_lionel'
            self.realyear = True
            self.realyear_origin = datetime(2000,1,1)

            folder= '/data0/project/meddle/gula/ROMS/Simulations/GS/'
            self.his=folder + '/all/natl_his_'
            self.grd=folder + 'natl6_grd_agrifV2.nc'
            self.frc=folder + '/AVG_Y2000M01/natl_his_Y2000M01.0861.nc'
            self.wind=folder + '/AVG_Y2000M01/natl_his_Y2000M01.0861.nc'
            self.tfile=1
            self.tstart=862
            self.tend=3053
        
        
        ##################
        
        elif 'polgyr_natl_mean' in simul:

            
            self.realyear = True
            self.realyear_origin = datetime(1999,1,1)
            self.model = 'croco'
            self.digits = 0
            
            if os.getenv('HOSTNAME') is None:
                folder= '/net/libra/local/tmp/1/gula/ROMS/Simulations/POLGYR_NATL/'
                self.his=folder + '/mean_2002_2008'
                self.grd=folder + '/polgyr_natl_grd.nc'
            
            self.tfile=20
            self.tstart=0
            self.tend=1000
        
        ##################
        
        elif 'polgyr_natl' in simul:

            self.model = 'croco'
            self.digits = 5

            if os.getenv('HOSTNAME') is None:
                if 'mean' in simul:
                    folder= '/data0/project/meddle/lecorre/POLGYR_NATL/'
                    self.his=libra +'/gula/ROMS/Simulations/POLGYR/polgyr_natl_avg.mean'
                    self.digits = 0
                else:
                    folder= '/data0/project/meddle/lecorre/POLGYR_NATL/'
                    self.his=folder + '/HIS/polgyr_his.'
                self.grd=libra +'/gula/ROMS/Simulations/POLGYR/polgyr_natl_grd.nc'
            else:
                folder= '/ccc/store/cont003/gen7638/lecorrem/POLGYR_EDDY/'
                self.his=folder + '/HIS/polgyr_his.'
                self.grd='/ccc/store/cont003/gen7638/lecorrem/INIT_NATL/polgyr_natl_grd.nc'

            self.frc=folder + '/HIS/polgyr_his.00000.nc'
            self.wind=folder + '/HIS/polgyr_his.00000.nc'
            self.tfile=20
            self.tstart=0
            self.tend=4000

        ##################
        
        elif 'polgyr_mean' in simul:
            
            self.realyear = True
            self.realyear_origin = datetime(1999,1,1)
            self.model = 'croco'
            self.digits = 0
            
            if os.getenv('HOSTNAME') is None:
                folder= '/net/libra/local/tmp/1/gula/ROMS/Simulations/POLGYR/'
                self.his=folder + '/mean_2002_2008'
                self.grd=folder + '/polgyr_grd.nc'
            
            self.tfile=20
            self.tstart=0
            self.tend=1000

        ##################
        
        elif 'rockall_mean' in simul:
            
            self.realyear = True
            self.realyear_origin = datetime(1999,1,1)
            self.model = 'croco'
            self.digits = 0
            
            if os.getenv('HOSTNAME') is None:
                folder= '/net/libra/local/tmp/1/gula/ROMS/Simulations/POLGYR/'
                if 'winter2' in simul:
                    self.his=folder + '/mean_2011_winter2'
                elif 'winter' in simul:
                    self.his=folder + '/mean_2011_winter'
                else:
                    self.his=folder + '/mean_2011'
                    
                self.grd=folder + '/polgyr_grd.nc'
            
            self.tfile=20
            self.tstart=0
            self.tend=1000
            
        ##################
        
        elif 'polgyr' in simul:
        
            '''
            time 00000 is 2001-01-10 14:53:20
            time 06511 is 2009-12-10 02:53:20
            
            conversion to uncompressed for 2006 = from 3620 to ...
            '''
            
            self.realyear = True
            self.realyear_origin = datetime(1999,1,1)
            self.model = 'croco'
            self.digits = 5

            #if os.getenv('HOSTNAME') is None:
                #folder= ROMSSIMS + '/POLGYR/'
                #self.grd= folder + 'polgyr_grd.nc'
            #else:
            import socket
            if ('lpo' in socket.gethostname()) or ('styx' in socket.gethostname()):
                folder='/scratch/Jcollin/croco/polgyr/'
                self.grd=folder + 'polgyr_grd.nc'
            else:    
                folder= '/home/datawork-lops-osi/mlecorre/POLGYR/HIS/'
                self.grd='/home/datawork-lops-osi/mlecorre/POLGYR/INIT/polgyr_grd.nc'

            self.his=folder +'polgyr_his.'

            self.frc=folder + '/polgyr_his.03360.nc'
            self.wind=folder + '/polgyr_his.03360.nc'
            self.tfile=20
            self.tstart=0
            self.tend=1000
        
        ######################
        elif 'uncompressed' in simul:
            self.realyear = True
            self.realyear_origin = datetime(1999,1,1)

            self.model = 'croco'
            self.digits = 5

            #folder= '/home/datawork-lops-osi/mlecorre/POLGYR/HIS/'
            folder = '/home/datawork-lemar-apero/jgula/POLGYR/HIS_uncompressed/'
            self.his=folder +'polgyr_his.'

            self.grd='/home/datawork-lops-osi/mlecorre/POLGYR/INIT/polgyr_grd.nc'
            #self.grd='/home2/datawork/lwang/IDYPOP/Data/ROMS/polgyr_grd.nc'

            self.frc=folder + '/HIS/polgyr_his.00000.nc'
            self.wind=folder + '/HIS/polgyr_his.00000.nc'

            self.tfile = 20
            self.tstart = 0
            self.tend = 184



        ##################
        
        elif 'rockall_natl' in simul:

            self.realyear = True
            self.realyear_origin = datetime(1999,1,1)
            self.model = 'croco'
            self.digits = 5

            if os.getenv('HOSTNAME') is None:
                folder= ROMSSIMS + '/POLGYR_NATL/'
                self.grd=folder + 'rockall_natl_grd.nc'
            else:
                folder= '/home/datawork-lops-osi/mlecorre/ROCKALL/PRT/HIS/'
                self.grd=folder + 'rockall_his.00000.nc'
                
            if 'sig' in simul:
                self.his=folder + 'vrtz/rockall_natl_sig_uv_chunked.'
            else:
                self.his=folder + 'rockall_his.'
                
            self.frc=folder + 'rockall_his.00000.nc'
            self.wind=folder + 'rockall_his.00000.nc'
            self.tfile=20
            self.tstart=0
            self.tend=1100

        ##################
        
        elif 'rockall' in simul or 'ROCK_chd' in simul:

            self.realyear = True
            self.realyear_origin = datetime(1999,1,1)
            self.model = 'croco'
            self.digits = 5

            if os.getenv('HOSTNAME') is None:
                folder= ROMSSIMS + '/POLGYR/'
                self.grd= folder + 'polgyr_grd.nc'
            else:
                folder= '/home/datawork-lops-osi/mlecorre/ROCKALL/CHD/HIS/'
                self.grd='/home/datawork-lops-osi/jgula/POLGYR/polgyr_grd.nc'


            self.his=folder +'rockall_his.'

            self.frc=folder + '/rockall_his.00000.nc'
            self.wind=folder + '/rockall_his.00000.nc'
            self.tfile=20
            self.tstart=0
            self.tend=1460


        
        ##################
        
        elif simul=='natl':

            folder=comanche + '/lrenault'
            self.his=folder + '/NATL/HIS/natl6V2_his.'
            self.grd=folder + '/NATL/natl6_grd_agrifV2.nc'
            self.frc=folder + '/NATL/natl6_frc.nc'
            self.wind=folder + '/NATL/natl6_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=915
        
        ##################
        
        elif simul=='natl6_avg':
            '''ask Lionel which one it is'''

            folder=stratus + '/lrenault/ROMS'
            self.his=folder + '/NATL6/natl6V2_avg.'
            self.grd=folder + '/NATL6/natl6_grd_agrifV2.nc'
            self.frc=folder + '/NATL6/natl6_frc.nc'
            self.wind=folder + '/NATL6/natl6_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=910
            #self.model = 'croco'
            #code = 'croco'
        
        ##################
        
        elif simul=='natl_mean1':

            folder=cherokee +'/gula/ROMS/Simulations/'
            self.his=cherokee +'/gula/ROMS/Simulations/NATL/natl6V2_his.mean.0100_0910.'
            self.grd=folder + '/NATL/natl6_grd_agrifV2.nc'
            self.frc=folder + '/NATL/natl6_frc.nc'
            self.wind=folder + '/NATL/natl6_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=915

        ##################
        
        elif simul=='natl_mean2':

            folder=cherokee +'/gula/ROMS/Simulations/'
            self.his=cherokee +'/gula/ROMS/Simulations/NATL/natl6V2_his.mean.0265_0910.'
            self.grd=folder + '/NATL/natl6_grd_agrifV2.nc'
            self.frc=folder + '/NATL/natl6_frc.nc'
            self.wind=folder + '/NATL/natl6_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=915
            
        ##################
        
        elif simul=='natl_mean3':

            folder=cherokee +'/gula/ROMS/Simulations/'
            self.his=cherokee +'/gula/ROMS/Simulations/NATL/natl6V2_his.mean.0070_0225.'
            self.grd=folder + '/NATL/natl6_grd_agrifV2.nc'
            self.frc=folder + '/NATL/natl6_frc.nc'
            self.wind=folder + '/NATL/natl6_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=915
        
        ##################
        
        elif simul=='natlc':

            folder=stratus + '/lrenault/Coupled/GS/NATL/ROMS/'
            self.his=folder + '/ALL/natl_his.'
            self.grd=folder + 'natl6_grd_agrifV2.nc'
            self.frc=folder + 'natl6_frc.nc'
            self.wind=folder + 'natl6_frc.nc'
            self.tfile=1
            self.tstart=862
            self.tend=2322

        ##################
        
        elif simul=='natlc_avg':

            folder=stratus + '/lrenault/Coupled/GS/NATL/'
            self.his=folder + '/ROMS/ALL/natl_avg.'
            self.grd=folder + '/ROMS/natl6_grd_agrifV2.nc'
            self.frc=folder + 'natl6_frc.nc'
            self.wind=folder + 'natl6_frc.nc'
            self.tfile=1
            self.tstart=0
            self.tend=2687

        ##################
        
        elif simul=='natlu_avg':

            folder=stratus + '/lrenault/Coupled/GS/NATL/'
            self.his=folder + '/ROMS_UNCOUPLED/ALL/natl_avg.'
            self.grd=folder + '/ROMS/natl6_grd_agrifV2.nc'
            self.frc=folder + 'natl6_frc.nc'
            self.wind=folder + 'natl6_frc.nc'
            self.tfile=1
            self.tstart=0
            self.tend=2687

        ##################
        
        elif simul=='natlc_avg_mean':
            '''
            mean 01-01-2000 to 12-31-2004
            '''
            self.digits = 0
            folder=ROMSSIMS
            self.his=folder +'/NATL/natlc_avg.mean.2000_2004'
            self.grd=folder + '/NATL/natl6_grd_agrifV2.nc'
            self.frc=folder + 'natl6_frc.nc'
            self.wind=folder + 'natl6_frc.nc'
            self.tfile=1
            self.tstart=0
            self.tend=2687

        ##################
        
        elif simul=='natlu_avg_mean':
            '''
            mean 01-01-2000 to 12-31-2004
            '''
            self.digits = 0
            folder=ROMSSIMS
            self.his=folder +'/NATL/natlu_avg.mean.2000_2004'
            self.grd=folder + '/NATL/natl6_grd_agrifV2.nc'
            self.frc=folder + 'natl6_frc.nc'
            self.wind=folder + 'natl6_frc.nc'
            self.tfile=1
            self.tstart=0
            self.tend=2687

        ##################
        
        elif simul=='natluz_avg_mean':
            '''
            simu uncoupled with Zob = 0.01
            
            mean 01-01-2000 to 12-31-2004
            '''
            self.digits = 0
            folder=ROMSSIMS
            self.his=folder +'/NATL/natluz_avg.mean.2000_2004'
            self.grd=folder + '/NATL/natl6_grd_agrifV2.nc'
            self.frc=folder + 'natl6_frc.nc'
            self.wind=folder + 'natl6_frc.nc'
            self.tfile=1
            self.tstart=861
            self.tend=2687
        ##################
        
        elif simul=='natlcz_avg_mean':
            '''
            simu coupled with Zob = 0.01
            
            mean 01-01-2000 to 12-31-2004
            '''
            self.digits = 0
            folder=ROMSSIMS
            self.his=folder +'/NATL/natlcz_avg.mean.2000_2004'
            self.grd=folder + '/NATL/natl6_grd_agrifV2.nc'
            self.frc=folder + 'natl6_frc.nc'
            self.wind=folder + 'natl6_frc.nc'
            self.tfile=1
            self.tstart=861
            self.tend=2687
            
            
        ################## Test case (to test online features)
    
            
        elif simul[:4]=='test':
            
            if 'ucla' in simul:
                code = 'ucla'
            else:
                code = 'croco'
                self.model = 'croco'
                self.digits = 5

            self.grd=capella + '/gula/test_roms/Config_ucla/roms_grd.nc'
            self.frc=capella + '/gula/test_roms/Config_ucla/roms_frc.nc'
            self.wind=self.frc
            self.tfile=100
            self.tstart=0
            self.tend=2000               

            if 'without' in simul:
                withuv = 'without_diags_uv/'
            else:
                withuv = 'with_diags_uv/'

            #folder = capella + '/gula/test_roms/Config_' + code +'/'+ withuv
            folder = capella + '/gula/test_roms/Config_' + code +'/'

            if 'avg' in simul:
                if 'uv' in simul:
                    self.his = folder + 'roms_diags_uv_avg.'
                elif 'vrt' in simul:
                    self.his = folder + 'roms_diags_vrt_avg.'
                elif 'ek' in simul:
                    self.his = folder + 'roms_diags_ek_avg.'
                else:
                    self.his = folder + 'roms_avg.'
            else:
                if 'uv' in simul:
                    self.his = folder + 'roms_diags_uv.'
                elif 'vrt' in simul:
                    self.his = folder + 'roms_diags_vrt.'
                elif 'ek' in simul:
                    self.his = folder + 'roms_diags_ek.'
                else:
                    self.his = folder + 'roms_his.'
        
            
        ################## Test case (to test online features)
    
            
        elif simul[:2]=='MM':
            
            code = 'croco'
            self.model = 'croco'
            self.digits = 5
 
            self.grd=capella + '/gula/test_roms/mat_croco/input/dipole_grd.nc'
            self.frc=capella + '/gula/test_roms/mat_croco/input/dipole_frc.nc'
            self.wind=self.frc
            self.tfile=100
            self.tstart=0
            self.tend=2000               

            
            if 'fast' in simul:
                folder = capella + '/gula/test_roms/mat_croco/fast/'
            elif 'full' in simul:
                folder = capella + '/gula/test_roms/mat_croco/full/'

            if 'avg' in simul:
                if 'uv' in simul:
                    self.his = folder + 'roms_diags_uv_avg.'
                elif 'vrt' in simul:
                    self.his = folder + 'roms_diags_vrt_avg.'
                elif 'ek' in simul:
                    self.his = folder + 'roms_diags_ek_avg.'
                else:
                    self.his = folder + 'roms_avg.'
            else:
                if 'uv' in simul:
                    self.his = folder + 'roms_diags_uv.'
                elif 'vrt' in simul:
                    self.his = folder + 'roms_diags_vrt.'
                elif 'ek' in simul:
                    self.his = folder + 'dipole_ek.'
                else:
                    self.his = folder + 'dipole_his.'

        ##################   
            
        ################## Test case (to test online features)
    
            
        elif simul[:5] =='basin':
            
            code = 'croco'
            self.model = 'croco'
            self.digits = 0
            
            if 'stommel' in simul:
                folder = '/home/gula/Desktop/Work_capella/roms/1617/case2_stommel'
            elif 'topo' in simul:
                folder = libra + '/gula/ROMS/test_roms/gyre_student_topo'
            else:
                folder = libra + '/gula/ROMS/test_roms/gyre_student'


            if 'LR' in simul:
                folder = folder + '_LR'


            self.grd=folder + '/basin_his.nc'
            self.frc=folder + '/basin_his.nc'
            self.wind=self.frc
            self.tfile=10000
            self.tstart=0
            self.tend=2000

            
            if 'mean' in simul:
                if 'vrt' in simul:
                    self.his = folder + '/basin_diags_vrt_mean'
                elif 'ek' in simul:
                    self.his = folder + '/basin_diags_ek_mean'
                else:
                    self.his = folder + '/basin_his_mean'
            else:
                if 'vrt' in simul:
                    self.his = folder + '/basin_diags_vrt_avg'
                elif 'ek' in simul:
                    self.his = folder + '/basin_diags_ek_avg'
                else:
                    self.his = folder + '/basin_avg'

        ################## Test case (to test online features)
    
            
        elif simul[:6] =='budget':
            
            code = 'croco'
            self.model = 'croco'
            self.digits = 0
            
            folder = libra + '/gula/ROMS/test_roms/test_budget_ucla'


            self.grd=folder + '/basin_his.nc'
            self.frc=folder + '/basin_his.nc'
            self.wind=self.frc
            self.tfile=10000
            self.tstart=0
            self.tend=2000

            
            if 'vrt' in simul:
                self.his = folder + '/basin_diags_vrt'
            elif 'ek' in simul:
                self.his = folder + '/basin_diags_ek'
            elif 'uv' in simul:
                self.his = folder + '/basin_diaM'
            elif 'ts' in simul:
                self.his = folder + '/basin_dia'
            else:
                self.his = folder + '/basin_his'

        ################## Test case (to test online features)
    
            
        elif 'gulfstream' in simul:
            
            code = 'croco'
            self.model = 'croco'
            self.digits = 0
            
            if 'mycroco' in simul:
                folder = libra + '/gula/ROMS/test_roms/test_budget_gulfstream_mycroco'
            elif 'dev_megatl' in simul:
                folder = libra + '/gula/ROMS/test_roms/test_budget_gulfstream_dev_megatl'
            elif 'tracer' in simul:
                folder = libra + '/gula/ROMS/test_roms/test_budget_gulfstream_tracer'
            else:
                folder = libra + '/gula/ROMS/test_roms/test_budget_gulfstream'

            self.grd=folder + '/CROCO_FILES/croco_grd.nc'
            self.frc=folder + '/CROCO_FILES/croco_frc.nc'
            self.wind=self.frc
            self.tfile=10000
            self.tstart=0
            self.tend=2000

            if 'RSUP3' in simul: out = '/OUT_RSUP3'
            elif 'RSUP5' in simul: out = '/OUT_RSUP5'
            elif 'UP3' in simul and 'NOTSADV' in simul: out = '/OUT_UP3_NOTSADV'
            elif 'UP3' in simul: out = '/OUT_UP3'
            elif 'UP5' in simul: out = '/OUT_UP5'
            elif 'WENO5' in simul: out = '/OUT_WENO5'
            else: out = '/OUT'

            if 'avg' in simul:
                if 'vrt' in simul:
                    self.his = folder + out + '/croco_diags_vrt_avg'
                elif 'ek' in simul:
                    self.his = folder + out + '/croco_diags_ek_avg'
                elif 'uv' in simul:
                    self.his = folder + out + '/croco_diaM_avg'
                elif 'ts' in simul:
                    self.his = folder + out + '/croco_dia_avg'
                elif 'pv' in simul:
                    self.his = folder + out + '/croco_diags_pv_avg'
                else:
                    self.his = folder + out + '/croco_avg'
            else:
                if 'vrt' in simul:
                    self.his = folder + out + '/croco_diags_vrt'
                elif 'ek' in simul:
                    self.his = folder + out + '/croco_diags_ek'
                elif 'uv' in simul:
                    self.his = folder + out + '/croco_diaM'
                elif 'ts' in simul:
                    self.his = folder + out + '/croco_dia'
                elif 'pv' in simul:
                    self.his = folder + out + '/croco_diags_pv'
                else:
                    self.his = folder + out + '/croco_his'


        ################## Test case (to test online features)
    
            
    
            
        elif simul[:8] =='caldeira':
 
            code = 'croco'
            self.model = 'croco'
            self.digits = 5

            if os.getenv('HOSTNAME') is None:
                #folder_code = libra + '/gula/ROMS/test_roms/test_pv_budget_caldeira'
                folder_code = libra + '/gula/ROMS/Simulations/Seamount'
            else:
                folder_code = '/home/datawork-lops-osi/jgula/SEAMOUNT'
   

            if 'SW' in simul:
                if 'noS' in simul: 
                    if 'wind' in simul and 'nodiff' in simul:    
                        folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_noS'
                    else:
                        folder = folder_code +'/Config_croco_SW_GLS_noS'
                elif 'large' in simul:                 
                    folder = folder_code +'/Config_croco_SW_GLS_large'
                elif 'dz' in simul and 'GLS' in simul:                 
                    folder = folder_code +'/Config_croco_SW_GLS_dz'
                elif 'HHR' in simul:
                    if 'KPP0' in simul:
                        if '_U' in simul:
                            if 'N' in simul:
                                folder = folder_code +'/Config_croco_SW_N_KPP0_HHR_U/case1'
                elif 'VHR' in simul:
                    if 'KPP0' in simul:
                        if '_U' in simul:
                            if 'N' in simul:
                                folder = folder_code +'/Config_croco_SW_N_KPP0_VHR_U/case1'
                    elif 'GLS' in simul:    
                        folder = folder_code +'/Config_croco_SW_GLS_VHR'
                    else:               
                        folder = folder_code +'/Config_croco_SW_VHR'
                elif 'HRzz' in simul:
                    if 'KPP0' in simul:    
                        if '_U' in simul: 
                            if 'N' in simul:
                                folder = folder_code +'/Config_croco_SW_N_KPP0_HRzz_U/case1'
                                #folder = '/home2/scratch/jgula/Seamount/Config_croco_SW_N_KPP0_VHR_U'
                elif 'HRz' in simul:
                    if 'KPP0' in simul:
                        if '_U' in simul:
                            if 'N' in simul:
                                folder = folder_code +'/Config_croco_SW_N_KPP0_HRz_U/case1'
                                #folder = '/home2/scratch/jgula/Seamount/Config_croco_SW_N_KPP0_VHR_U'
                elif 'HR' in simul:
                    if 'GLSD' in simul:    
                        if '_U' in simul: 
                            if 'N' in simul:
                                folder = folder_code +'/Config_croco_SW_N_GLSD_HR_U/case1'
                                #folder = '/home2/scratch/jgula/Seamount/Config_croco_SW_N_GLSD_HR_U'
                            else:
                                folder = folder_code +'/Config_croco_SW_GLSD_HR_U/case1'
                    elif 'GLS0' in simul:    
                        if '_U' in simul: 
                            folder = folder_code +'/Config_croco_SW_GLS0_HR_U/case1'
                    elif 'KPP0D' in simul:    
                        if '_U' in simul: 
                            if 'N' in simul:
                                folder = folder_code +'/Config_croco_SW_N_KPP0D_HR_U/case1'
                                #folder = '/home2/scratch/jgula/Seamount/Config_croco_SW_N_KPP0D_HR_U'
                    elif 'KPPDNOKT' in simul:    
                        if '_U' in simul: 
                            if 'N' in simul:
                                folder = folder_code +'/Config_croco_SW_N_KPPDNOKT_HR_U/case1'
                                #folder = '/home2/scratch/jgula/Seamount/Config_croco_SW_N_KPPDNOKT_HR_U'
                    elif 'KPPD' in simul:    
                        if '_U' in simul: 
                            if 'N' in simul:
                                folder = folder_code +'/Config_croco_SW_N_KPPD_HR_U/case1'
                                #folder = '/home2/scratch/jgula/Seamount/Config_croco_SW_N_KPPD_HR_U'
                    elif 'KPP0DNOKT' in simul:    
                        if '_U' in simul: 
                            if 'N' in simul:
                                folder = folder_code +'/Config_croco_SW_N_KPP0DNOKT_HR_U/case1'
                                #folder = '/home2/scratch/jgula/Seamount/Config_croco_SW_N_KPP0DNOKT_HR_U'
                    elif 'KPP0' in simul:    
                        if '_U' in simul: 
                            if 'N' in simul:
                                folder = folder_code +'/Config_croco_SW_N_KPP0_HR_U/case1'
                                #folder = '/home2/scratch/jgula/Seamount/Config_croco_SW_N_KPP0_HR_U'
                            else:
                                folder = folder_code +'/Config_croco_SW_KPP0_HR_U/case1'
                        else:          
                            folder = folder_code +'/Config_croco_SW_KPP_HR'
                    elif 'KPPNOR' in simul:    
                        if '_U' in simul: 
                            folder = folder_code +'/Config_croco_SW_KPPNOR_HR_U/case1'
                    elif 'KPP' in simul:    
                        if '_U' in simul: 
                            folder = folder_code +'/Config_croco_SW_KPP_HR_U/case1'
                    elif 'GLS' in simul:    
                        if '_U' in simul: 
                            if 'N' in simul:
                                folder = folder_code +'/Config_croco_SW_N_GLS_HR_U/case1'
                            else:
                                folder = folder_code +'/Config_croco_SW_GLS_HR_U/case1'
                        else:          
                            folder = folder_code +'/Config_croco_SW_GLS_HR'
                    else:   
                        if '_U' in simul:    
                            folder = folder_code +'/Config_croco_SW_HR_U/case1'
                        else:          
                            folder = folder_code +'/Config_croco_SW_HR'
                   
                elif 'wind' in simul:                 
                    if 'GLS' in simul:
                        if 'nodiff' in simul:
                            if 'vadvc2' in simul:
                                if 'hadvc2' in simul:
                                    folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_vadvc2_hadvc2'
                                elif 'hadvc4' in simul:
                                    if 'test' in simul:
                                        folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_vadvc2_hadvc4_test'
                                    elif 'dif4' in simul:
                                        if 'vis4' in simul:
                                            folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_vadvc2_hadvc4_dif4_vis4'
                                        else:
                                            folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_vadvc2_hadvc4_dif4'
                                    else:
                                        folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_vadvc2_hadvc4'
                                elif 'hadvweno5' in simul:
                                    folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_vadvc2_hadvweno5'
                                elif 'vis4' in simul:
                                    folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_vadvc2_vis4'
                                elif 'testv3' in simul:
                                    folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_vadvc2_testv3'
                                elif 'testv2' in simul:
                                    folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_vadvc2_testv2'
                                elif 'test' in simul:
                                    folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_vadvc2_test'
                                else:
                                    folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_vadvc2'
                            elif 'vadvweno5' in simul:
                                folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff_vadvweno5'
                            else:
                                folder = folder_code +'/Config_croco_SW_GLS_wind_nodiff'
                        else:
                            if 'hadvup3' in simul:
                                folder = folder_code +'/Config_croco_SW_GLS_wind_hadvup3'
                            else:
                                folder = folder_code +'/Config_croco_SW_GLS_wind'
                    else:
                        folder = folder_code +'/Config_croco_SW_KPP_wind'
                elif 'GLS' in simul:             
                    folder = folder_code +'/Config_croco_SW_GLS'
                elif 'drag' in simul:             
                    folder = folder_code +'/Config_croco_SW_drag'
                elif 'clean' in simul:             
                    folder = folder_code +'/Config_croco_SW_clean'
                elif 'dz' in simul:             
                    folder = folder_code +'/Config_croco_SW_dz'
                else:
                    folder = folder_code +'/Config_croco_SW'
            else:
                if 'nodiff' in simul:
                    folder = folder_code +'/Config_croco_GLS_nodiff'
                elif 'flat' in simul:
                    folder = folder_code +'/Config_croco_GLS_flat'
                else:
                    folder = folder_code +'/Config_croco_GLS'


            self.grd=folder + '/caldeira_his.00000.nc'
            self.frc=folder + '/roms_frc.nc'
            self.wind=self.frc
            self.tfile=1000
            self.tstart=0
            self.tend=10000


            if 'avg' in simul:          
                if 'vrt' in simul:
                    self.his = folder + '/caldeira_diags_vrt_avg.'
                elif 'ek' in simul:
                    self.his = folder + '/caldeira_diags_ek_avg.'
                elif 'uv' in simul:
                    self.his = folder + '/caldeira_diaM_avg.'
                elif 'ts' in simul:
                    self.his = folder + '/caldeira_dia_avg.'
                elif 'pv' in simul:
                    self.his = folder + '/caldeira_diags_pv_avg.'
                else:
                    self.his = folder + '/caldeira_avg.'
            else:
                if 'vrt' in simul:
                    self.his = folder + '/caldeira_diags_vrt.'
                elif 'ek' in simul:
                    self.his = folder + '/caldeira_diags_ek.'
                elif 'uv' in simul:
                    self.his = folder + '/caldeira_diaM.'
                elif 'ts' in simul:
                    self.his = folder + '/caldeira_dia.'
                elif 'pv' in simul:
                    self.his = folder + '/caldeira_diags_pv.'
                else:
                    self.his = folder + '/caldeira_his.'



        ################## Test case (to test online features)
    
            
        elif simul[:5] =='nesec':
            
            code = 'croco'
            self.model = 'croco'
            self.digits = 5

            if '1h' in simul:
                if 'withdiag' in simul:
                    #folder = libra + '/gula/ROMS/Simulations/NESEA/HIS1h_withdiags'
                    folder = ROMSSIMSK + '/NESEA/HIS/HIS1h_withdiags'
                else:
                    folder = ROMSSIMSK + '/NESEA/HIS/HIS1h'
            else:
                folder = ROMSSIMSK + '/NESEA/HIS/HIS6h'


            self.grd=ROMSSIMSK + '/NESEA/nesea_grd.nc'
            self.frc=ROMSSIMSK + '/NESEA/nesea_frc.nc'
            self.wind=self.frc
            self.tfile=20
            self.tstart=0
            self.tend=2000

            
            if 'vrt' in simul:
                self.his = folder + '/nesec_diags_vrt.'
            elif 'ek' in simul:
                self.his = folder + '/nesec_diags_ek.'
            elif 'uv' in simul:
                self.his = folder + '/nesec_uv.'
            elif 'ts' in simul:
                self.his = folder + '/nesec_ts.'
            elif 'pv' in simul:
                self.his = folder + '/nesec_diags_pv.'
            else:
                self.his = folder + '/nesec_his.'

        ##################

        elif simul=='nesed_avg':
            
            code = 'croco'
            self.model = 'croco'
            self.digits = 5

            folder = ROMSSIMSK + '/NESED'

            self.grd=ROMSSIMSK + '/NESED/nesea_grd.nc'
            self.frc=ROMSSIMSK + '/NESEA/nesea_frc.nc'
            self.wind=self.frc
            self.tfile=20
            self.tstart=60
            self.tend=2000
            self.his = folder + '/nesed_avg.'

        ################## Test case (to test online features)
    

        elif simul=='WOES':

            self.model = 'agrif_jc'
            #folder=  '/media/gula/WOES_ZOOM2/SCRATCH_WOES_VTR1_ZOOM2/AVGFILES'
            folder=  '/net/capella/local/tmp/2/gula/WOES'
            self.his=folder + '/roms_avg_'
            self.fileformat='.nc.1'
            self.grd=folder + '/roms_avg_Y1991M1.nc.1'
            self.frc=folder + '/roms_avg_Y1991M1.nc.1'
            self.wind=folder + '/roms_avg_Y1991M1.nc.1'
            self.tfile=6
            self.tstart=0
            self.Ystart=1990
            self.Mstart=1

            self.tend=10000

        ##################
        
        elif 'nesev' in simul:
            
            code = 'croco'
            self.model = 'croco'
            self.digits = 5

            if 'KPP' in simul or 'kpp' in simul:

                if 'libra' in socket.gethostname() or 'lpo' in socket.gethostname():
                    folder = ROMSSIMSK + '/NESEV/NESEV_KPP'
                else:
                    folder = '/home/datawork-lops-osi/jgula/NESEV_KPP_28_16/OUT'
                    
            elif 'KEPS1' in simul or 'keps1' in simul:

                if 'libra' in socket.gethostname() or 'lpo' in socket.gethostname():
                    folder = ROMSSIMSK + '/NESEV/NESEV_KEPS1'
                else:
                    folder = '/home/datawork-lops-osi/jgula/NESEV_KEPS1_28_16/OUT'
                    
            else:

                if 'libra' in socket.gethostname() or 'lpo' in socket.gethostname():
                    folder = ROMSSIMSK + '/NESEV/NESEV'
                else:
                    folder = '/home/datawork-lops-osi/jgula/NESEV_28_16/OUT'


            self.grd = folder + '/nesea_grd.nc'
            self.frc = folder + '/nesea_frc.nc'
            self.wind = self.frc

            self.tstart=0
            self.tend=2000

            if 'mean' in simul:
            
                self.tfile=1
                self.digits = 0

                if 'vrt' in simul:
                    self.his = folder + '/nesev_diags_vrt_avg.mean'
                elif 'ek' in simul:
                    self.his = folder + '/nesev_diags_ek_avg.mean'
                elif 'uv' in simul:
                    self.his = folder + '/nesev_uv_avg.mean'
                elif 'ts' in simul:
                    self.his = folder + '/nesev_ts_avg.mean'
                elif 'pv' in simul:
                    self.his = folder + '/nesev_diags_pv_avg.mean'
                elif 'eddy' in simul:
                    self.his = folder + '/nesev_diags_eddy_avg.mean'
                elif 'his500' in simul:
                    #self.his = folder + '/nesev_his.mean.00250-00949'
                    self.his = folder + '/nesev_his.mean.00500-01399'
                elif 'his800' in simul:
                    self.his = folder + '/nesev_his.mean.00800-00899'
                elif 'avg33' in simul:
                    self.his = folder + '/nesev_avg.mean_033_132'
                elif 'his' in simul:
                    self.his = folder + '/nesev_his.mean.00250-00949'
                else:
                    self.his = folder + '/nesev_avg.mean'
                
                
            elif 'avg' in simul:
            
                self.tfile=1
                
                if 'vrt' in simul:
                    self.his = folder + '/nesec_diags_vrt_avg.'
                elif 'ek' in simul:
                    self.his = folder + '/nesev_diags_ek_avg.'
                elif 'uv' in simul:
                    self.his = folder + '/nesev_uv_avg.'
                elif 'ts' in simul:
                    self.his = folder + '/nesev_ts_avg.'
                elif 'pv' in simul:
                    self.his = folder + '/nesev_diags_pv_avg.'
                elif 'eddy' in simul:
                    self.his = folder + '/nesev_diags_eddy_avg.'
                else:
                    self.his = folder + '/nesev_avg.'
                    
            elif 'uvbot' in simul:

                self.his = folder + '/nesev_uvbot.'
                self.tfile= 2
                    
            else:
            
                if (time<224):
                    self.tfile=8
                else:
                    self.tfile=2
                    
                self.his = folder + '/nesev_his.'
                

                
        ################## AGULHAS TEST CASE
    

        elif simul=='woes5':

            self.model = 'agrif_jc'
            folder=  '/data0/project/vortex/Master_stage_tedesco2017/SCRATCH_RUN_WOES5/'
            self.his=folder + '/roms_avg_'
            self.fileformat='.nc.2'
            self.grd=folder + '/roms_grd.nc.2'
            self.frc=folder + '/roms_avg_Y1993M1.nc.2'
            self.wind=folder + '/roms_avg_Y1993M1.nc.2'
            self.tfile=31
            self.tstart=0
            self.Ystart=1993
            self.Mstart=1

            self.tend=10000
    


        elif simul=='smagu':

            self.model = 'croco'
            folder= libra +'/gula/ROMS/Simulations/SMAGU/'
            self.his=folder + '/roms_avg_'
            self.fileformat='.nc'
            self.grd=folder + '/roms_grd.nc'
            self.frc=folder + '/roms_grd.nc'
            self.wind=folder + '/roms_grd.nc'
            self.tfile=10
            self.tstart=0
            self.tend=10000


        elif simul=='smagu3':

            self.model = 'croco'
            self.digits = 5
            folder= '/data0/project/vortex/Master_stage_tedesco2017/SCRATCH_RUN_WOES5/Pauline/SMAGU_3/HIS/new_grd/'
            self.his=folder + '/smagu_his_3_newgrd.'
            self.fileformat='.nc'
            self.grd=folder + '/smagu_his_3_newgrd.00000.nc'
            self.frc=folder + '/smagu_his_3_newgrd.00000.nc'
            self.wind=folder + '/smagu_his_3_newgrd.00000.nc'
            self.tfile=121
            self.tstart=0
            self.tend=10000


        elif simul=='smagu3_mean':

            self.model = 'croco'
            self.digits = 0
            folder= '/data0/project/vortex/Master_stage_tedesco2017/SCRATCH_RUN_WOES5/Pauline/SMAGU_3/HIS/new_grd/'
            self.his=folder + '/smagu_his_3_newgrd_ymean'
            self.fileformat='.nc'
            self.grd=folder + '/smagu_his_3_newgrd.00000.nc'
            self.frc=folder + '/smagu_his_3_newgrd.00000.nc'
            self.wind=folder + '/smagu_his_3_newgrd.00000.nc'
            self.tfile=121
            self.tstart=0
            self.tend=10000


        ##################
        elif 'chabam' in simul:

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/gulaj/CHABAM/'
                    self.grd=folder + 'chaba_grd.nc'
                    self.tstart=0
            except:
                folder=krypton +'/gula/ROMS/Simulations/CHABAM/'
                self.grd=krypton +'/gula/ROMS/Simulations/CHABAM/chaba_grd.nc'
                self.tstart=250

            self.frc=folder + 'chaba_frc_monthly.nc'
            self.wind=folder + 'chaba_frc_monthly.nc'
            self.tfile=5
            self.tend=1820

            if 'surf_avg' in simul:
                self.his = folder + '/sfc/chabam_surf_avg.'
                self.tfile=30
                self.tend=9480
            elif 'surf' in simul:
                self.his = folder + '/sfc/chabam_surf.'
                self.tfile=30
                self.tend=9480
            else:
                if 'sig' in simul:
                    self.his='/data0/project/meddle/gula/ROMS/Simulations/CHABAM/iso/chabam_sig_uv.'
                    self.digits = 4
                else:
                    self.his=folder + 'chabam_his.'

        ######
        elif 'chabah' in simul:

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/gulaj/CHABAH/'
                    self.grd=folder + 'chaba_grd.nc'
                    self.tstart=0
            except:
                folder=krypton +'/gula/ROMS/Simulations/CHABAH/'
                self.grd=krypton +'/gula/ROMS/Simulations/CHABAH/chaba_grd.nc'
                self.tstart=250
            
            self.frc=folder + '../CHABAT/chabat_frc_1h.nc'
            self.wind=self.frc
            self.tfile=5
            self.tend=1820

            if 'surf_avg' in simul:
                self.his = folder + '/sfc/chabah_surf_avg.'
                self.tfile=30
                self.tend=9480
            elif 'surf' in simul:
                self.his = folder + 'chabah_surf.'
                self.tfile=30
                self.tend=9480
            else:
                if 'sig' in simul:
                    self.his='/data0/project/meddle/gula/ROMS/Simulations/CHABAH/iso/chabah_sig_uv.'
                    self.digits = 4
                else:
                    self.his=folder + 'chabah_his.'

        ######
        elif 'chabat' in simul:

            self.model = 'croco'
            self.digits = 5

            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/gulaj/CHABAT/'
                    self.grd=folder + 'chaba_grd.nc'
                    self.tstart=0

            except:
                folder=krypton +'/gula/ROMS/Simulations/CHABAT/'
                self.grd=libra +'/gula/ROMS/Simulations/CHABA/chaba_grd.nc'
                self.tstart=250
            
            self.frc=folder + 'chabat_frc_1h.nc'
            self.wind=self.frc
            self.tfile=5
            self.tend=1820

            if 'surf_avg' in simul:
                self.his = folder + '/sfc/chabat_surf_avg.'
                self.tfile=30
                self.tend=9480
            elif 'surf' in simul:
                self.his = folder + '/sfc/chabat_surf.'
                self.tfile=30
                self.tend=9480
            else:
                if 'sig' in simul:
                    self.his='/data0/project/meddle/gula/ROMS/Simulations/CHABAT/iso/chabat_sig_uv.'
                    self.digits = 4
                else:
                    self.his=folder + 'chabat_his.'


        ##################

        elif simul=='chaba':

            self.model = 'croco'
            folder=ROMSSIMSK
            self.his=folder + '/CHABA/HIS_frc_1h/chaba_his.'
            self.grd=folder + '/CHABA/chaba_grd.nc'
            self.frc=folder + '/CHABA/chaba_frc_1h.nc'
            self.wind=folder + '/CHABA/chaba_frc_1h.nc'
            self.tfile=5
            self.tstart=0
            self.tend=1555


        elif simul=='chaba_avg':

            self.model = 'croco'
            self.his=libra +'/gula/ROMS/Simulations/CHABA/HIS_frc_1h/chaba_avg.'
            self.grd=libra +'/gula/ROMS/Simulations/CHABA/chaba_grd.nc'
            self.frc=libra +'/gula/ROMS/Simulations/CHABA/chaba_frc_1h.nc'
            self.wind=libra +'/gula/ROMS/Simulations/CHABA/chaba_frc_1h.nc'
            self.tfile=5
            self.tstart=0
            self.tend=1555

        elif simul in ['sarga','SARGA']:

            folder= '/data0/project/meddle/gula/ROMS/Simulations'
            self.model = 'croco'
            self.digits = 5
            self.his=folder + '/SARGA/HIS/sarga_his.'
            self.grd=folder + '/SARGA/sarga_grd.nc'
            self.frc=folder + '/SARGA/sarga_frc_all.nc'
            self.wind=folder + '/SARGA/sarga_frc_all.nc'
            self.tfile=4
            self.tstart=0
            self.tend=384

        elif simul=='sargo':

            folder= '/data0/project/meddle/gula/ROMS/Simulations'
            self.model = 'croco'
            self.digits = 5
            self.his=folder + '/SARGA/HIS/sargo_his.'
            self.grd=folder + '/SARGA/sarga_grd.nc'
            self.frc=folder + '/SARGA/sarga_frc_all.nc'
            self.wind=folder + '/SARGA/sarga_frc_all.nc'
            self.tfile=4
            self.tstart=0
            self.tend=10000
            
        elif simul in ['sargasurf']:

            folder= '/data0/project/meddle/gula/ROMS/Simulations'
            self.model = 'croco'
            self.digits = 5
            self.his=folder + '/SARGA/HIS/sarga_surf.'
            self.grd=folder + '/SARGA/sarga_grd.nc'
            self.frc=folder + '/SARGA/sarga_frc_all.nc'
            self.wind=folder + '/SARGA/sarga_frc_all.nc'
            self.tfile=30
            self.tstart=0
            self.tend=4620

        elif simul=='sargosurf':

            folder= '/data0/project/meddle/gula/ROMS/Simulations'
            self.model = 'croco'
            self.digits = 5
            self.his=folder + '/SARGA/HIS/sargo_surf.'
            self.grd=folder + '/SARGA/sarga_grd.nc'
            self.frc=folder + '/SARGA/sarga_frc_all.nc'
            self.wind=folder + '/SARGA/sarga_frc_all.nc'
            self.tfile=30
            self.tstart=0
            self.tend=1980

        ##################

        elif simul=='lucky':

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/vicc/LUCKY/'
                    self.his=folder + 'HIS/lucky_his.'
                    self.tstart=0
            except:
                folder= libra +'/gula/ROMS/Simulations/LUCKY/'
                self.his=folder + 'lucky_his.'
                self.tstart=1000
            self.grd=folder + 'lucky_grd.nc'
            self.frc=folder + 'lucky_frc.nc'
            self.wind=folder + 'lucky_frc.nc'
            self.tfile=2
            self.tend=6118
        ##################

        elif simul=='ridge':

            self.model = 'croco'
            self.digits = 4
            folder= '/data0/project/meddle/cvic/ROMS/Simulations/RIDGE/'
            self.his=folder + 'HIS/ridge_his.'
            self.tstart=0
            self.grd=folder + 'ridge_grd.nc'
            self.frc=folder + 'ridge_frc.nc'
            self.wind=folder + 'ridge_frc.nc'
            self.tfile=5
            self.tend=515

        ##################

        elif simul=='lucky_avg':

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/vicc/LUCKY/'
            except:
                folder= libra +'/gula/ROMS/Simulations/LUCKY/'
            self.his=folder + 'AVG/lucky_avg.'
            self.grd=folder + 'lucky_grd.nc'
            self.frc=folder + 'lucky_frc.nc'
            self.wind=folder + 'lucky_frc.nc'
            self.tfile=2
            self.tstart=0
            self.tend=6118
        
        ##################

        elif simul=='luckyt':

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/vicc/LUCKYT/'
                    self.grd=folder + 'luckyt_grd.nc'
                    self.his=folder + 'HIS/luckyt_his.'
                    self.tstart=0
            except:
                folder= libra +'/gula/ROMS/Simulations/LUCKYT/'
                self.grd=folder + 'lucky_grd.nc'
                self.his=folder + 'luckyt_his.'
                self.tstart=1000
            self.frc=folder + 'luckyt_frc.nc'
            self.wind=folder + 'luckyt_frc.nc'
            self.tfile=4
            self.tend=3096

        ##################

        elif simul=='luckyt_avg':

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/vicc/LUCKYT/'
                    self.grd=folder + 'luckyt_grd.nc'
            except:
                folder= libra +'/gula/ROMS/Simulations/LUCKY/'
                self.grd=folder + 'lucky_grd.nc'
            self.his=folder + 'AVG/luckyt_avg.'
            self.frc=folder + 'luckyt_frc.nc'
            self.wind=folder + 'luckyt_frc.nc'
            self.tfile=4
            self.tstart=0
            self.tend=10000


        ##################

        elif simul=='luckyto':

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/lahayen/LUCKYTO/'
                    self.grd=folder + 'luckyto_grd.nc'
                    self.his=folder + 'HIS/luckyto_his.'
                    self.tstart=0
            except:
                folder= libra +'/gula/ROMS/Simulations/LUCKY/'
                self.grd=folder + 'lucky_grd.nc'
                self.his=folder + 'luckyt_his.'
                self.tstart=1000
            self.frc=folder + 'luckyto_frc.nc'
            self.wind=folder + 'luckyto_frc.nc'
            self.tfile=4
            self.tend=1444


        ##################

        elif simul=='luckym2w':

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/lahayen/LUCKYM2/'
                    self.grd=folder + 'luckym2_grd.nc'
                    self.his=folder + 'HIS/luckym2_his.'
                    self.tstart=0
            except:
                folder= libra +'/gula/ROMS/Simulations/LUCKY/'
                self.grd=folder + 'lucky_grd.nc'
                self.his=folder + 'luckym2_his.'
                self.tstart=1000

            self.frc=folder + 'luckym2_frc.nc'
            self.wind=folder + 'luckym2_frc.nc'
            self.tfile=4
            self.tend=1080


        ##################

        elif simul=='luckym2s':

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/lahayen/LUCKYM2s/'
                    self.grd=folder + 'luckyto_grd.nc'
                    self.his=folder + 'HIS/luckyto_his.'
                    self.tstart=0
            except:
                folder= libra +'/gula/ROMS/Simulations/LUCKY/'
                self.grd=folder + 'lucky_grd.nc'
                self.his=folder + 'luckym2_his.'
                self.tstart=1000

            self.frc=folder + 'luckyto_frc.nc'
            self.wind=folder + 'luckyto_frc.nc'
            self.tfile=4
            self.tend=720


        ##################

        ##################

        elif 'ratlbig' in simul and 'gigatl6' not in simul:

            self.model = 'croco'
            self.digits = 5

            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    if 'gls' in simul:
                        folder='/ccc/store/cont003/gen7638/gulaj/RATLBIG_GLS/'
                    else:
                        folder='/ccc/store/cont003/gen7638/gulaj/RATLBIG/'
                    self.grd='/ccc/store/cont003/gen7638/gulaj/RATLBIG/ratlbig_grd.nc'
            except:
                if 'gls' in simul:
                    folder=ROMSSIMS + '/RATLBIG_GLS/'
                else:
                    folder=ROMSSIMS + '/RATLBIG/'
                self.grd=folder + 'ratlbig_grd.nc'

            if 'mean' in simul:
                self.his=folder + 'ratlbig_avg.mean'
                self.digits = 0
            elif 'avg' in simul:
                self.his=folder + 'ratlbig_avg.'
            else:
                self.his=folder + 'ratlbig_his.'

            self.frc=folder + 'ratlbig_frc.nc'
            self.wind=folder + 'ratlbig_frc.nc'
            self.tfile=10
            self.tstart=0
            self.tend=10000


        ##################

        elif simul=='bimin':

            self.model = 'croco'
            self.digits = 5
            if os.getenv('HOSTNAME')=='curie91':
                folder='/ccc/store/cont003/gen7638/gulaj/BIMIN/'
                self.grd=folder + 'bimin_grd.nc'
            else:
                folder= libra +'/gula/ROMS/Simulations/BIMIN/'
                self.grd=folder + 'bimin_grd.nc'
            self.his=folder + 'bimin_his.'
            self.frc=folder + 'bimin_frc.nc'
            self.wind=folder + 'bimin_frc.nc'
            self.tfile=4
            self.tstart=0
            self.tend=3096

        ##################

        elif 'nbiminh' in simul:

            self.model = 'croco'
            self.digits = 4
            folder ='/data0/project/nhmg/ducousso/CROCO-NH/NBIMIN/INPUT/'
 
            self.grd=folder + 'nbimin_grd.nc'  
            self.frc=folder + 'nbimin_frc.nc'
            self.wind=folder + 'nbimin_frc.nc'

            folder= '/data0/project/nhmg/ducousso/CROCO-NH/NBIMIN2/EXP_H_dt3_nft50_OBC2DSPEC3DSPEC_SPGnuintwidthnam_5d/'

            self.his=folder + 'HIS/nbimin_his.'  

            self.tfile=4
            self.tstart=0
            self.tend=200
        ##################

        elif 'nbiminnh' in simul:

            self.model = 'croco'
            self.digits = 4
            folder ='/data0/project/nhmg/ducousso/CROCO-NH/NBIMIN/INPUT/'
 
            self.grd=folder + 'nbimin_grd.nc'  
            self.frc=folder + 'nbimin_frc.nc'
            self.wind=folder + 'nbimin_frc.nc'

            folder= '/data0/project/nhmg/ducousso/CROCO-NH/NBIMIN2/EXP_NHMG_dt3_nft50_OBC2DSPEC3DSPEC_SPGnuintwidthnam_5d/'

            self.his=folder + 'HIS/nbimin_his.'  

            self.tfile=4
            self.tstart=0
            self.tend=200

        ##################

        elif 'biminh' in simul:

            self.model = 'croco'
            self.digits = 5
            folder= libra +'/gula/ROMS/Simulations/BIMIN/'
            self.grd=folder + 'bimin_grd.nc'  

            if 'mean' in simul:
                self.his= '/data0/project/nhmg/gula/biminh_his.mean'
                self.digits = 0
            elif 'sig' in simul:
                self.his=folder + 'HIS/biminh_sig_uv.'  
            else:
                self.his=folder + 'HIS/biminh_his.'  
            
            self.frc=folder + 'bimin_frc.nc'
            self.wind=folder + 'bimin_frc.nc'
            self.tfile=1000
            self.tstart=0
            self.tend=241


        ##################

        elif 'biminnh' in simul:

            self.model = 'croco'
            self.digits = 5
            folder= libra +'/gula/ROMS/Simulations/BIMIN/'
            self.grd=folder + 'bimin_grd.nc'

            if 'mean' in simul:
                self.his= '/data0/project/nhmg/gula/biminnh_his.mean'
                self.digits = 0
            elif 'sig' in simul:
                self.his=folder + 'HIS/biminnh_sig_uv.'  
            else:
                self.his=folder + 'HIS/biminnh_his.'  


            
            self.frc=folder + 'bimin_frc.nc'
            self.wind=folder + 'bimin_frc.nc'
            self.tfile=1000
            self.tstart=0
            self.tend=241



        ##################

        elif simul=='sismi':

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/gulaj/SISMI/'
            except:
                folder= libra +'/gula/ROMS/Simulations/SISMI/'
            self.grd=folder + 'sismi_grd.nc'
            self.his=folder + 'sismi_his.'
            self.frc=folder + 'sismi_frc.nc'
            self.wind=folder + 'sismi_frc.nc'
            self.tfile=10
            self.tstart=0
            self.tend=3096

        ##################

        elif 'sismo' in simul:

            self.model = 'croco'
            self.digits = 5

            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/gulaj/SISMO/'
            except:
                folder= '/data0/project/meddle/gula/ROMS/Simulations/SISMO/'
                #folder= libra +'/gula/ROMS/Simulations/SISMO/'

            self.grd=folder + 'sismo_grd.nc'

            if 'avg' in simul:          
                if 'pv' in simul:
                    self.his = folder + 'sismo_diags_pv_avg.'
                elif 'sig' in simul:
                    self.his = folder + 'sig2/sismo_avg_sig2_uv.'
                    self.digits = 4
                elif 'zeta' in simul:
                    self.his = folder + 'HIS/sismo_avg_zeta.'
                else:
                    self.his = folder + 'sismo_avg.'
            else:
                if 'pv' in simul:
                    self.his = folder + 'sismo_diags_pv.'
                else:
                    self.his = folder + 'sismo_his.'

            self.frc=folder + 'sismo_frc.nc'
            self.wind=folder + 'sismo_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=3096


        ##################

        elif simul=='sismu':

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/gulaj/SISMU/'
            except:
                folder= libra +'/gula/ROMS/Simulations/SISMU/'
            self.grd=folder + 'sismo_grd.nc'
            self.his=folder + 'sismu_his.'
            self.frc=folder + 'sismu_frc.nc'
            self.wind=folder + 'sismu_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=3096

        ##################

        elif simul=='sismu_avg':

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/gulaj/SISMU/'
            except:
                folder= libra +'/gula/ROMS/Simulations/SISMU/'
            self.grd=folder + 'sismo_grd.nc'
            self.his=folder + 'sismu_avg.'
            self.frc=folder + 'sismu_frc.nc'
            self.wind=folder + 'sismu_frc.nc'
            self.tfile=5
            self.tstart=0
            self.tend=3096

        ##################

        elif simul in ['bahaz','bahaz_mean']:

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/gulaj/BAHAZ/'
            except:
                folder= libra +'/gula/ROMS/Simulations/BAHAZ/'

            self.grd=folder + 'bahaz_grd.nc'

            if 'mean' in simul:
                self.his=folder + 'bahaz_avg.mean.00010-00150'
                self.digits = 0
            else:
                self.his=folder + 'bahaz_his.'

            self.frc=folder + 'bahaz_frc.nc'
            self.wind=folder + 'bahaz_frc.nc'
            self.tfile=10
            self.tstart=0
            self.tend=10000

        ##################

        elif simul in ['bahazr','bahazr_mean']:

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/gulaj/BAHAZR/'
            except:
                folder= libra +'/gula/ROMS/Simulations/BAHAZR/'

            self.grd=folder + 'bahazr_grd.nc'

            if 'mean' in simul:
                self.his=folder + 'bahazr_avg.mean.00010-00150'
                self.digits = 0
            else:
                self.his=folder + 'bahazr_his.'

            self.frc=folder + 'bahazr_frc.nc'
            self.wind=folder + 'bahazr_frc.nc'
            self.tfile=10
            self.tstart=0
            self.tend=10000

        ##################

        elif simul in ['bahaz_gls','bahaz_gls_mean']:

            self.model = 'croco'
            self.digits = 5
            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/gulaj/BAHAZ_GLS/'
                    self.grd='/ccc/store/cont003/gen7638/gulaj/BAHAZ/bahaz_grd.nc'
            except:
                folder= libra +'/gula/ROMS/Simulations/BAHAZ_GLS/'

            self.grd=folder + '../BAHAZ/bahaz_grd.nc'

            if 'mean' in simul:
                self.his=folder + 'bahaz_avg.mean.00010-00150'
                self.digits = 0
            else:
                self.his=folder + 'bahaz_his.'

            self.frc=folder + 'bahaz_frc.nc'
            self.wind=folder + 'bahaz_frc.nc'
            self.tfile=10
            self.tstart=0
            self.tend=10000

        ##################        


        elif 'megatl9' in simul:

            self.realyear = True
            self.realyear_origin = datetime(1979,1,1)
            self.model = 'croco'
            self.digits = 5
            
            if 'erai' in simul:
                if 'kpp' in simul and 'nos' in simul:
                    folder_name = 'MEGATL9_ERAI_KPP_noS/'
                elif 'kpp' in simul:
                    folder_name = 'MEGATL9_ERAI_KPP/'
                elif 'v2' in simul:
                    folder_name = 'MEGATL9_ERAI_v2/'
                else:
                    folder_name = 'MEGATL9_ERAI/'
            elif 'cfsr' in simul:
                folder_name = 'MEGATL9_CFSR/'
            elif 'clim' in simul:
                if 'n50' in simul:
                    if 'newgrd' in simul:
                        folder_name = 'MEGATL9_N50_CLIM_newgrd/'
                    else:
                        folder_name = 'MEGATL9_N50_CLIM/'
                else:
                    folder_name = 'MEGATL9_CLIM/'

            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/gulaj/' + folder_name
            except:
                folder= libra +'/gula/ROMS/Simulations/MEGATL/MEGATL9/' + folder_name

            if 'newgrd' in simul:
                self.grd= libra +'/gula/ROMS/Simulations/MEGATL/MEGATL9/megatl9_grd_v2.nc'
            else:
                self.grd= libra +'/gula/ROMS/Simulations/MEGATL/MEGATL9/megatl9_grd.nc'

            if 'mean' in simul:
                self.his=folder + 'megatl9_avg.mean'
                self.digits = 0
            else:
                self.his=folder + 'megatl9_his.'

            self.frc=folder + 'megatl9_frc.nc'
            self.wind=folder + 'megatl9_frc.nc'
            self.tfile=10
            self.tstart=0
            self.tend=10000

        ##################        

        elif 'megatl6' in simul:

            self.realyear = True
            self.realyear_origin = datetime(1979,1,1)
            self.model = 'croco'
            self.digits = 5
            
            if 'erai' in simul:
                if 'kpp' in simul and 'nos' in simul:
                    folder_name = 'MEGATL6_ERAI_KPP_noS/'
                elif 'kpp' in simul:
                    folder_name = 'MEGATL6_ERAI_KPP/'
                elif 'v2' in simul:
                    folder_name = 'MEGATL6_ERAI_v2/'
                else:
                    folder_name = 'MEGATL6_ERAI/'
            elif 'cfsr6h_on' in simul:
                    folder_name = 'MEGATL6_N50_CFSR6h_on/'
            elif 'cfsr1h_on_full' in simul:
                if 'tides' in simul:
                    folder_name = 'MEGATL6_N50_CFSR1h_on_full_tides/'
                elif 'v2' in simul:
                    folder_name = 'MEGATL6_N50_CFSR1h_on_full_v2/'
                elif 'v3' in simul:
                    folder_name = 'MEGATL6_N50_CFSR1h_on_full_v3/'
                else:
                    folder_name = 'MEGATL6_N50_CFSR1h_on_full/'
            elif 'cfsr' in simul:
                if 'n50' in simul:
                    if 'v2' in simul:
                        folder_name = 'MEGATL6_N50_CFSR_v2/'
                    else:
                        folder_name = 'MEGATL6_N50_CFSR/'
                else:
                    folder_name = 'MEGATL6_CFSR/'
            elif 'clim' in simul:
                if 'n50' in simul:
                    if 'newgrd' in simul:
                        if 'kpp'  in simul:
                            if 'geo' in simul:
                                folder_name = 'MEGATL6_N50_CLIM_newgrd_kpp_geo/'
                            else:
                                folder_name = 'MEGATL6_N50_CLIM_newgrd_kpp/'
                        else:
                            if 'sig' in simul:
                                folder_name = 'MEGATL6_N50_CLIM_newgrd_sig/'
                            else:
                                folder_name = 'MEGATL6_N50_CLIM_newgrd/'
                    else:
                        folder_name = 'MEGATL6_N50_CLIM/'
                else:
                    folder_name = 'MEGATL6_CLIM/'


            if os.getenv('HOSTNAME') is None:
                #folder= libra +'/gula/ROMS/Simulations/MEGATL/MEGATL6/' + folder_name
                folder= '/home/gula/datarmor/megatl/MEGATL6/OUT/' + folder_name
            elif 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                folder='/ccc/store/cont003/gen7638/gulaj/' + folder_name
            else:
                folder= '/home/datawork-lops-osi/jgula/' + folder_name

            if 'newgrd' in simul:
                self.grd= folder + '../megatl6_newgrd.nc'
            else:
                self.grd= folder + '../megatl6_grd.nc'

            
            if 'mean_short' in simul:
                self.his=folder + 'megatl6_avg.mean_short'
                self.avg=folder + 'megatl6_avg.mean_short'
                self.vrt=folder + 'megatl6_diags_vrt_avg.mean_short'
                self.ek=folder + 'megatl6_diags_ek_avg.mean_short'
                self.digits = 0
            elif 'mean_2002' in simul:
                self.his=folder + 'megatl6_avg.mean_2002'
                self.digits = 0
            elif 'mean' in simul:
                self.his=folder + 'megatl6_avg.mean'
                self.avg=folder + 'megatl6_avg.mean'
                self.vrt=folder + 'megatl6_diags_vrt_avg.mean'
                self.ek=folder + 'megatl6_diags_ek_avg.mean'
                self.digits = 0
            else:
                if 'avg' in simul:
                    self.his=folder + 'megatl6_avg.'
                else:
                    self.his=folder + 'megatl6_his.'
                self.avg=folder + 'megatl6_avg.'
                self.vrt=folder + 'megatl6_diags_vrt_avg.'
                self.ek=folder + 'megatl6_diags_ek_avg.'

            self.frc=folder + 'megatl6_frc.nc'
            self.wind=folder + 'megatl6_frc.nc'
            self.tfile=10
            self.tstart=0
            self.tend=10000


        ##################        

        elif 'megatl3' in simul:

            self.realyear = True
            self.realyear_origin = datetime(1979,1,1)
            self.model = 'croco'
            self.digits = 5

            folder_name = 'MEGATL3_v2/'

            if os.getenv('HOSTNAME') is None:
                folder= libra +'/gula/ROMS/Simulations/MEGATL/MEGATL3/' + folder_name
            elif 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                folder='/ccc/store/cont003/gen7638/gulaj/MEGATL3/' + folder_name
            else:
                folder= '/home/datawork-lops-osi/jgula/' + folder_name

            self.grd= folder + '../megatl3_grd.nc'

            self.his=folder + 'megatl3_his.'

            self.tfile=5
            self.tstart=0
            self.tend=10000

        ##################        


        elif 'megatl' in simul:

            self.realyear = True
            self.realyear_origin = datetime(1979,1,1)
            self.model = 'croco'
            self.digits = 5
            
            folder= libra +'/gula/ROMS/Simulations/MEGATL/MEGATL/'

            self.grd= folder + 'megatl_grd.nc'

            self.his=folder + 'megatl_his.'

            self.tfile=24
            self.tstart=0
            self.tend=10000


        
        ##################


        elif 'gigatl1' in simul:
        

            self.realyear = True
            self.realyear_origin = datetime(1979,1,1)
            self.realyear_tstart = datetime(2004,1,15)

            self.model = 'croco'
            self.digits = 6
            
            if 'gigatl1_1h_tides' in simul:
                folder_name= 'GIGATL1_1h_tides/'
            elif 'gigatl1_1h' in simul:
                folder_name= 'GIGATL1_1h/'


            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen7638/gulaj/GIGATL3/' + folder_name
                    
                else:
                    folder= '/home/datawork-lops-megatl/GIGATL1/' + folder_name
                    self.grd= '/home/datawork-lops-megatl/GIGATL1/gigatl1_grd.nc'
            except:
                #folder= '/data0/project/meddle/gula/ROMS/Simulations/GIGATL/GIGATL1/' + folder_name
                #self.grd= libra +'/gula/ROMS/Simulations/GIGATL/GIGATL1/gigatl1_grd_masked.nc'
                folder= '/net/omega/local/tmp/1/gula/GIGATL1/' + folder_name
                self.grd= '/net/omega/local/tmp/1/gula/GIGATL1/gigatl1_grd.nc'


            if 'para' in  simul:
                if output: print('loading tiles ',simul[-4:])
                tile = [s for s in simul.split('_') if s.isdigit()][-1]
                #tile = simul[-4:]
                
                self.grd= libra +'/gula/ROMS/Simulations/GIGATL/GIGATL1/INIT_N100_100_100/GRD3/gigatl1_grd_masked.' +  tile + '.nc'
                
                if 'surf' in simul:
                    self.his = folder + 'gigatl1_surf.' + tile
                    self.tfile=120
                    self.dtfile=3600
                    self.tend=100000
                    
                elif 'avg' in simul:

                    self.digits = 0
                    folder= '/data0/project/meddle/gula/ROMS/Simulations/GIGATL/GIGATL1/' + folder_name
                    
                    self.tfile=1
                    self.dtfile=5*24*3600
                    self.tend=10000
                    
                    if 'ek' in simul:
                        self.his = folder + 'gigatl1_diags_ek_avg.' + tile
                    elif 'vrt' in simul:
                        self.his = folder + 'gigatl1_diags_vrt_avg.' + tile
                    elif 'uv' in simul:
                        self.his = folder + 'gigatl1_diaM_avg.' + tile
                    elif 'pv' in simul:
                        self.his = folder + 'gigatl1_diags_pv_avg.' + tile
                    elif 'eddy' in simul:
                        self.his = folder + 'gigatl1_diags_eddy_avg.' + tile
                    else:
                        self.his = folder + 'gigatl1_avg.' + tile

                elif 'hourly' in simul:
                
                    self.his = folder + 'HIS_1h/gigatl1_his.' + tile + '.'
                    self.fileformat=  '.nc'
                    self.tfile=6
                    self.dtfile=3600
                    self.tend=50000
                    
                    self.model = 'croco_gigatl1_hours'
                    self.digits = 0


                else:
                
                    self.his = folder + 'gigatl1_his.'
                    self.fileformat=  '.' + tile + '.nc'
                    self.tfile=2
                    self.dtfile=3*3600
                    self.tend=10000
                    
            else:
                if 'surf' in simul:
                
                    if 'avg' in simul:
                        self.his = folder + 'SURF/final/gigatl1_surf_avg.'
                    else:
                        self.his = folder + 'SURF/final/gigatl1_surf.' #1999-01-25-1999-01-29'
                        
                    self.tfile=120
                    self.dtfile=3600
                    self.tend=100000
                    
                    self.model = 'croco_gigatl1'
                    self.digits = 0


        ##################        



        elif 'gigatl3' in simul and 'sig' in simul:

            self.realyear = True
            self.realyear_origin = datetime(1979,1,1)
            self.realyear_tstart = datetime(2004,1,15)

            self.model = 'croco'
            self.digits = 5

            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/work/cont003/gen12051/gulaj/Post/GIGATL/iso'
                    #if 'gigatl3_1h_tides' in simul:
                    #    folder = '/ccc/scratch/cont003/gen12051/gulaj/GIGATL3/' + folder_name + 'HIS/'
                    self.grd= '/ccc/store/cont003/gen12051/gulaj/GIGATL3/gigatl3_grd.nc'
            except:
                folder= libra +'/gula/ROMS/Simulations/GIGATL/GIGATL3/'
                self.grd= libra +'/gula/ROMS/Simulations/GIGATL/GIGATL3/gigatl3_grd.nc'


            if 'gigatl3_1h' in simul:
                self.his = folder + '/gigatl3_1h_sig_uv_chunked.'
                self.dtfile=12*3600
                self.tfile=10
                self.tend=100000

       ##################

        elif 'gigatl3' in simul and 'equator' in simul:

            self.realyear = True
            self.realyear_origin = datetime(1979,1,1)
            self.realyear_tstart = datetime(2004,1,15)

            self.model = 'croco'
            self.digits = 5

            folder= '/net/omega/local/tmp/1/gula/GIGATL3/zoom_for_seamount/'
            self.grd= '/net/omega/local/tmp/1/gula/GIGATL3/zoom_for_seamount/gigatl3_1h_tides_avg_equator.00500.nc'

            if 'mean' in simul:
                self.his = folder + '/gigatl3_1h_tides_avg_equator.00500-00573'
                self.tfile=1
                self.digits = 0
            elif 'avg' in simul:
                self.his = folder + '/gigatl3_1h_tides_avg_equator.'
                self.dtfile=5*24*3600
                self.tfile=1
                self.tend=100000
       ##################
       
       
        elif 'gigatl3' in simul:

            self.realyear = True
            self.realyear_origin = datetime(1979,1,1)
            self.realyear_tstart = datetime(2004,1,15)

            self.model = 'croco_xios'
            self.digits = 0
            
            if 'gigatl3_1h_tides' in simul:
                folder_name= 'GIGATL3_1h_tides/'
            elif 'gigatl3_1h_UP5' in simul:
                folder_name= 'GIGATL3_1h_UP5/'
            elif 'gigatl3_1h' in simul:
                folder_name= 'GIGATL3_1h/'
            elif 'gigatl3_6h' in simul:
                folder_name= 'GIGATL3_6h/'


            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    folder='/ccc/store/cont003/gen12051/gulaj/GIGATL3/' + folder_name
                    #if 'gigatl3_1h_tides' in simul:
                    #    folder = '/ccc/scratch/cont003/gen12051/gulaj/GIGATL3/' + folder_name + 'HIS/'
                    self.grd= '/ccc/store/cont003/gen12051/gulaj/GIGATL3/gigatl3_grd.nc'
            except:
                folder= '/net/omega/local/tmp/1/gula/GIGATL3/' + folder_name
                #folder= libra +'/gula/ROMS/Simulations/GIGATL/GIGATL3/' + folder_name
                self.grd= libra +'/gula/ROMS/Simulations/GIGATL/GIGATL3/gigatl3_grd.nc'


            if 'surf' in simul:
                self.his = folder + 'GIGATL3_1h_inst_surf_' #1999-01-25-1999-01-29'
                self.tfile=120
                self.dtfile=3600
                self.tend=100000
            elif 'hourly' in simul:
                self.his = folder + 'GIGATL3_1h_avg_3d_' #1999-01-25-1999-01-29'
                self.tfile=120
                self.dtfile=3600
                self.tend=100000
            elif 'avg' in simul:
                self.his = folder + 'GIGATL3_5d_aver_' #1999-01-25-1999-01-29'
                self.tfile=1
                self.dtfile=5*24*3600
                self.tend=10000
            elif 'M' in simul:
                self.his = folder + 'GIGATL3_5d_M_' #1999-01-25-1999-01-29'
                self.tfile=1
                self.dtfile=5*24*3600
                self.tend=10000
            elif 'uvbot' in simul:

                if 'gigatl3_1h_tides' in simul:
                    self.his = '/net/omega/local/tmp/1/gula/GIGATL3/uvbot/gigatl3_1h_tides_uvbot.'
                    self.dtfile=3*3600
                    self.tfile=40
                elif 'gigatl3_1h' in simul:
                    if time < 2330:
                        self.his = '/net/omega/local/tmp/1/gula/GIGATL3/uvbot/gigatl3_1h_UP5_uvbot.'
                    else:
                        self.his = '/net/omega/local/tmp/1/gula/GIGATL3/uvbot/gigatl3_1h_uvbot.'
                    self.dtfile=12*3600
                    self.tfile=10
                    
                self.model = 'croco'
                self.digits = 5
                self.tend=100000

            elif 'horiz' in simul:

                if 'gigatl3_1h_tides' in simul:
                    self.his = '/net/omega/local/tmp/1/gula/GIGATL3/horizontal_section/gigatl3_1h_tides_horizontal_section.'
                    self.dtfile=3*3600
                    self.tfile=40
                elif 'gigatl3_1h' in simul:
                    if time < 2330:
                        self.his = '/net/omega/local/tmp/1/gula/GIGATL3/horizontal_section/gigatl3_1h_UP5_horizontal_section.'
                    else:
                        self.his = '/net/omega/local/tmp/1/gula/GIGATL3/horizontal_section/gigatl3_1h_horizontal_section.'
                    self.dtfile=12*3600
                    self.tfile=10

                self.model = 'croco'
                self.digits = 5
                self.tend=100000               
 
            else:
                if 'mean' in simul:
                    self.his = folder + 'GIGATL3_5d_aver_2008-2012.mean'
                    self.dtfile=1
                    self.tfile=1
                    self.digits = 0
                    self.model = 'croco'
                else:

                    if 'gigatl3_1h_tides' in simul:
                        self.his = folder + 'GIGATL3_3h_inst_' #1999-01-25-1999-01-29'
                        self.dtfile=3*3600
                        self.tfile=40
                    elif 'gigatl3_1h' in simul:
                        self.his = folder + 'GIGATL3_12h_inst_' #1999-01-25-1999-01-29'
                        self.dtfile=12*3600
                        self.tfile=10
                    self.tend=100000

       ##################        


        elif 'gigatl6' in simul:
        
            self.realyear = True
            self.realyear_origin = datetime(1979,1,1)

            #self.realyear_tstart = datetime(2004,1,15)
            self.realyear_tstart = datetime(1999,1,15)

            self.model = 'croco_xios'
            self.digits = 0

            if 'ratlbig_GLS_devweek_M3FAST_bulk' in simul:
                folder_name= 'RATLBIG_GLS_devweek_M3FAST_bulk/'
            elif 'ratlbig_GLS_devweek_M3FAST' in simul:
                folder_name= 'RATLBIG_GLS_devweek_M3FAST/'
            elif 'ratlbig_GLS_devweek' in simul:
                folder_name= 'RATLBIG_GLS_devweek/'
            elif 'gigatl6_6h_noCFB' in simul:
                folder_name= 'GIGATL6_6h_noCFB/'
            elif 'gigatl6_6h' in simul:
                self.realyear_tstart = datetime(2004,1,15)
                folder_name= 'GIGATL6_6h/'
            elif 'gigatl6_1h_GEO' in simul:
                self.realyear_tstart = datetime(2004,1,15)
                folder_name= 'GIGATL6_1h_GEO/'
            elif 'gigatl6_1h_UP3_GEO' in simul:
                self.realyear_tstart = datetime(2004,1,15)
                folder_name= 'GIGATL6_1h_UP3_GEO/'
            elif 'gigatl6_1h_UP3' in simul:
                self.realyear_tstart = datetime(2004,1,15)
                folder_name= 'GIGATL6_1h_UP3/'
            elif 'gigatl6_1h_UP5_GEO5' in simul:
                self.realyear_tstart = datetime(2004,1,15)
                folder_name= 'GIGATL6_1h_UP5_GEO5/'
            elif 'gigatl6_1h_tides' in simul:
                self.realyear_tstart = datetime(2004,1,15)
                folder_name= 'GIGATL6_1h_tides/'
            elif 'gigatl6_1h' in simul:
                self.realyear_tstart = datetime(2004,1,15)
                folder_name= 'GIGATL6_1h/'
            else:
                folder_name= 'GIGATL6/'

            try:
                if 'curie' in os.getenv('HOSTNAME') or 'irene' in os.getenv('HOSTNAME'):
                    #folder='/ccc/store/cont003/gen12051/gulaj/' + folder_name
                    folder = '/ccc/scratch/cont003/gen12051/gulaj/GIGATL6/' + folder_name + 'HIS/'
                    if 'atlbig' not in simul:
                        self.grd= '/ccc/store/cont003/gen12051/gulaj/GIGATL6/gigatl6_grd.nc'
                    else:
                        self.grd= '/ccc/store/cont003/gen12051/gulaj/GIGATL6/ratlbig_grd.nc'
                else: # 'datarmor' in os.getenv('HOSTNAME'):
                    folder= '/home/datawork-lops-megatl/GIGATL6/' + folder_name + 'HIS/'
                    if 'atlbig' not in simul:
                        self.grd= '/home/datawork-lops-megatl/GIGATL6/gigatl6_grd.nc'
            except:
                folder= libra +'/gula/ROMS/Simulations/GIGATL/GIGATL6/' + folder_name
                if 'atlbig' not in simul:
                    self.grd= libra +'/gula/ROMS/Simulations/GIGATL/GIGATL6/gigatl6_grd.nc'
                #folder= '/home/gula/datarmor/megatl/GIGATL6/' + folder_name + 'HIS/'
                #if 'atlbig' not in simul:
                #    self.grd= '/home/gula/datarmor/megatl/GIGATL6/gigatl6_grd.nc'



            if 'surf' in simul:
                self.his = folder + 'GIGATL6_1h_inst_surf_' #1999-01-25-1999-01-29'
                self.tfile=120
                self.dtfile=3600
                self.tend=100000
            elif 'mean' in simul:
                self.his=folder + 'GIGATL6_5d_aver_2002.mean'
                self.model = 'croco'
                self.digits = 0
                self.tfile=1; self.dtfile=1; self.tend=11
            elif 'hourly' in simul:
                self.his = folder + 'GIGATL6_1h_avg_3d_' 
                self.tfile=120
                self.dtfile=3600
                self.tend=100000
            elif 'avg' in simul:
                if 'vrt' in simul:
                    self.his = folder + '/GIGATL6_5d_VRT_'
                elif 'KE' in simul:
                    self.his = folder + '/GIGATL6_5d_KE_'
                elif 'uv' in simul:
                    self.his = folder + '/GIGATL6_5d_M_'
                elif 'ts' in simul:
                    self.his = folder + '/GIGATL6_5d_TS_'
                else:
                    self.his = folder + 'GIGATL6_5d_aver_' #1999-01-25-1999-01-29'
                self.tfile=1
                self.dtfile=5*24*3600
                self.tend=10000
            else:
                self.his = folder + 'GIGATL6_12h_inst_' #1999-01-25-1999-01-29'
                self.tfile=10
                self.dtfile=12*3600
                self.tend=10000


        ##################

        elif simul=='leewa':

            self.model = 'croco'
            self.digits = 5

            self.tfile=24
            self.tstart=0
            self.tend=576

            if os.getenv('HOSTNAME') is None:
                folder= libra +'/gula/ROMS/Simulations/LEEWA/'
                self.grd=folder + 'leewa_grd.nc'
                self.his=folder + 'leewa_avg.'
            else:
                folder= '/home/datawork-lops-oh/meddle/cdemarez/STAGE_M2/final_v2/'
                self.his=folder + 'leewa_avg.'
                self.grd='/home2/datawork/jgula/Post/LEEWA/leewa_grd.nc'

            self.frc=folder + 'leewa_frc.nc'
            self.wind=folder + 'leewa_frc.nc'

        ##################

        elif simul=='leewa_smoothed':

            self.model = 'croco'
            self.digits = 5

            self.tfile=24
            self.tstart=0
            self.tend=576

            if os.getenv('HOSTNAME') is None:
                folder= libra +'/gula/ROMS/Simulations/LEEWA/'
                self.grd=folder + 'leewa_grd.nc'
                self.his=folder + 'leewa_avg.'
            else:
                folder= '/home/datawork-lops-oh/meddle/cdemarez/STAGE_M2/LEEWA_smoothed/final/'
                self.his=folder + 'leewa_avg.'
                self.grd='/home2/datawork/jgula/Post/LEEWA/leewa_smoothed_grd.nc'


            self.frc=folder + 'leewa_frc.nc'
            self.wind=folder + 'leewa_frc.nc'

        ##################

        elif 'rrexnum100' in simul:
        
            self.realyear = True
            self.realyear_origin = datetime(1979,1,1)
            self.realyear_tstart = datetime(2008,8,26,12)

            self.model = 'croco_xios'
            self.digits = 0

            folder= '/home/datawork-lops-rrex/cvic/RREXNUM100'
            self.grd= '/home/datawork-lops-rrex/jgula/RREXNUM100/rrexnum_grd.nc'
            #self.grd= folder + '/HIS/@expname@_@freq@_grid.nc'

            if '1d' in simul:
                self.his = folder + '/HIS/RREXNUM100_1d_avg_' #1999-01-25-1999-01-29'
                self.tfile=5
                self.dtfile=86400
                self.tend=100000
            elif '1h' in simul:
                self.his = folder + '/HIS/RREXNUM100_1h_avg_' #1999-01-25-1999-01-29'
                self.tfile=24
                self.dtfile=3600
                self.tend=100000

        ##################

        elif 'brazil' in simul:
        
            self.realyear = True
            self.realyear_origin = datetime(1979,1,1)
            self.realyear_tstart = datetime(2008,2,23,12)

            self.model = 'croco_xios'
            self.digits = 0

            if os.getenv('HOSTNAME') is None:
                folder= libra +'/gula/ROMS/Simulations/BRAZIL/'
                self.grd=folder + 'brazil_grd.nc'
            else:
                folder= '/home/datawork-lops-osi/jgula/BRAZIL'
                self.grd= folder + '/INIT/brazil_grd.nc'


            if 'mean' in simul:
                self.his = folder + '/BRAZIL_5d_aver_mean'
                self.digits = 0
                self.model = 'croco'
                self.tfile=1; self.dtfile=1; self.tend=1
            elif 'avg' in simul:
                self.his = folder + '/HIS2/BRAZIL_5d_aver_'
                self.tfile=1
                self.dtfile=432000
                self.tend=100000
            else:
                self.his = folder + '/HIS2/BRAZIL_12h_inst_'
                self.tfile=10
                self.dtfile=12*3600
                self.tend=100000

        ##################
        
        
        
        else:

            print("""

I never heard about your simulation name, please try again...

Choices are:

'gulfs': 1.5km of Gulf Stream region (monthly forcing)
'hatt': 500m ... son of GULFS ... north of Cape Haterras
'hatt2': 500m ... son of GULFS ... south of Cape Haterras
'shing': 150m ... son of HATT2 ... close to Cape Haterras
'filam': 150m ... same than SHING with a larger output frequency (10min Vs 3h)

'atlbig': 5km of Atlantic basin (daily forcing + diurnal cycle + no Jeroen correction)
'gulfz': 1.5km ...  son of ATLBIG...  Gulf Stream region (same domain than GULFS)
'gulfz_seas': Climatology of gulfz
'nesea': 500m ... son of GULFZ ... along the Gulf Stream path above New-England seamounts
'neseb': 500m ... same than SHING with a larger output frequency (1h Vs 12h)
'nefro': 150m ... son of NESEB
'nefro0':150m ... son of NESEA (wrong due to undersampling of boundary forcing)

'seamount': idealized seamount test-case
'roms': used to visualize grid just created with the EGRID tools (on inca)
                   """)
            #sys.exit()

###################################################################################


    def findfile(path):
        '''Find the file named path in the sys.path.
        Returns the full path name if found, None if not found'''
        for dirname in sys.path:
            possible = os.path.join(dirname, path)
            if os.path.isfile(possible):
                return possible
        return None


###################################################################################
# end class files
###################################################################################













