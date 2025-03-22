
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

###################################################################################
#   Main 
###################################################################################
    def __init__(self,simulname=None,time=None, floattype=float, light = False, touchfile=True, output =True, **kwargs):

        """

        """
        if output: print('simulname is',simulname)

        #define type of variables
        self.floattype = floattype

        #get simulation parameters from arguments
        if 'simul' in  kwargs: 
            self.load_file(kwargs['simul'].split(' '),output = output)
        #elif len(sys.argv)>0:
        #    self.load_file(sys.argv)
        else: 
            self.load_file([simulname],output = output)
        
        #for time in range(self.time0, self.ncname.tend, self.dtime ):
        if time==None: time = self.time0
        self.time = int(time)
        
        self.ncname=files(self.simul, time=self.time, output = output)
        
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
            self.ncfile = self.ncname.his+'Y' +format(year)+'M'+format(month) + self.ncname.fileformat
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
    def update(self,time=None, output =True):

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
            
        self.g = 9.81
        
        try:
            self.rho0 = ncfile.rho0 
            self.hc = ncfile.hc
            self.dt_model = ncfile.dt
        except:
            try:
                self.rho0 = ncfile.variables['rho0']
                self.hc = ncfile.variables['hc']
                self.dt_model = ncfile.variables['dt']
            except:
                pass


        try:
            self.Cs_r = self.Forder(ncfile.variables['Cs_r'][:])
            self.Cs_w = self.Forder(ncfile.variables['Cs_w'][:])
            #self.sc_r = self.Forder(ncfile.variables['sc_r'][:])
            #self.sc_w = self.Forder(ncfile.variables['sc_w'][:])
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
            self.VertCoordType = 'NEW'

        ncfile.close()

            
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

    def dt(self,output=True):

        ncfile = Dataset(self.ncfile, 'r')
        if output: print('dt is read in ',self.ncfile)
        try:
            if self.ncname.model in ['ucla','croco_lionel']:
                self.dt = np.array(ncfile.variables['ocean_time'][-1]) \
                        - np.array(ncfile.variables['ocean_time'][-2])
            elif  self.ncname.model in ['croco','agrif_jc']:
                self.dt = np.array(ncfile.variables['scrum_time'][-1]) \
                        - np.array(ncfile.variables['scrum_time'][-2])
            elif  self.ncname.model in ['croco_xios']:
                self.dt = np.array(ncfile.variables['time_counter'][-1]) \
                        - np.array(ncfile.variables['time_counter'][-2])
            else:
                self.dt = np.array(ncfile.variables['time'][-1]) \
                        - np.array(ncfile.variables['time'][-2])
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
###################################################################################


    def __init__(self, simul, time=0, output =True):

        ###################################################

        # default values (might be changed in simulation definition)
        self.model = 'ucla'
        self.fileformat='.nc'
        self.digits = 4
        self.realyear = False # False means using 360 days and 30 days per month

        ##################

        if 'polgyr' in simul:
        
            '''
            time 00000 is 2001-01-10 14:53:20
            '''
            
            self.realyear = True
            self.realyear_origin = datetime(1999,1,1)
            self.model = 'croco'
            self.digits = 5

            import socket
            if 'lpo' in socket.gethostname():
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

        ##################

        elif 'rutgers' in simul:   
        
            '''
            Test with outputs from ROMS Rutgers
            '''
            
            self.realyear = True
            self.realyear_origin = datetime(1999,1,1)
            self.model = 'ucla'
            self.digits = 0

            folder= '/home/gula/libra/ROMS/Simulations/Rutgers_example/'
            self.grd= folder + 'grd.nc'

            self.his=folder +'velocity_example'
            self.tfile=2
            self.tstart=0
            self.tend=2

        ##################
        
        
        
        else:

            print("""

I never heard about your simulation name,
please add the definition in Modules/R_files.py

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













