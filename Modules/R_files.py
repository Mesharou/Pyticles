
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

21/11/23 : Jcollin
        - removed old 'agrif_jc' and 'croco_lionel
        - added croco_date 
            - realyear = True: internannual files
            - realyear = False: climatology run (1 file per month and 30 days per file)
        - added loadgrid option (if False does not load grid at __init__) this
        is handy to play around with load and update methods without having
        file on disk. Also used for testing
        - added test function 
        - regrouped filenaming and file stepping in load.filedata() method
        called in __int__() and update()

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
    """

###################################################################################
#   Main 
###################################################################################
    def __init__(self, simulname=None, time=None, floattype=float, light=False,
                touchfile=True, output=True, loadgrid=True, **kwargs):

        if output: print('simulname is',simulname)

        # define type of variables
        self.floattype = floattype

        # get simulation parameters from arguments
        if 'simul' in  kwargs: 
            self.load_file(kwargs['simul'].split(' '), output=output)
        else: 
            self.load_file([simulname], output=output)
        
        # time initialization
        if time==None:
            time = self.time0
        self.time = int(time)
        
        # get netcdf metadata and path
        self.ncname = files(self.simul, time=self.time, output = output)
        self.filetime, self.infiletime, self.ncfile = self.filedata(
            output=output)
        
        # ---- cutted here -------

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

        if loadgrid:
            if light:
                [self.mask,self.pm,self.pn,self.x,self.y] = self.variables_grd(output = output, light=True)
            else:
                [self.topo,self.mask,self.pm,self.pn,self.f,self.x,self.y,self.angle] = self.variables_grd(output = output)

            if self.simul in ['shing','reshing','filam','dilam','filam_avg','dilam_avg','slope','hatt2','gulfs']: 
                self.topo_corrected = self.variables_grd_corrected()

# ----------------------------------------------------------------------------
    def __repr__(self):
        "print data about simul"
        return f"{self.time=}\n{self.ncfile=}\n{self.infiletime=}"

# ----------------------------------------------------------------------------
    def filedata(self, output=True):
        """
        set filetime, infiletime and ncfile
        
        - 'ucla', 'croco': croco/roms_his_1200.nc
        
        - 'croco_date': classic croco formats croco_Y2000M1.nc
           depends on self.realyear
           if True: interannual, files with a file per month and an
             ouput per day
           else: climatology files with a file per month and an ouput per day
        
        - croco_xios or croco_gigatl:
          find first and last date in file to reconstruct name,
            ex: 1999-01-25-1999-01-29
        
        """
         
        if self.ncname.model in ['ucla', 'croco']:
            infiletime = self.time % self.ncname.tfile
            filetime = self.time - infiletime
            if self.ncname.digits == 4:
                ncfile = "".join([self.ncname.his, f"{filetime:04}",
                                self.ncname.fileformat])
            elif self.ncname.digits == 5:
                ncfile = "".join([self.ncname.his, f"{filetime:05}",
                                self.ncname.fileformat])
            elif self.ncname.digits == 6:
                ncfile = "".join([self.ncname.his, f"{filetime:06}",
                                self.ncname.fileformat])
            elif self.ncname.digits == 0:
                ncfile = self.ncname.his + self.ncname.fileformat

        elif self.ncname.model == 'croco_date':
            if self.ncname.realyear:
                date = self.ncname.realyear_origin + timedelta(
                    days=self.time-self.ncname.tstart)
                year = date.year
                month = int(date.month)
                day = date.day
                infiletime = day - 1
                filetime = self.time
                ncfile = "".join([self.ncname.his, 'Y', f"{year:04}",
                                'M', str(month), self.ncname.fileformat])
                
            else:
                dsec = 30*24*3600 //self.ncname.tfile 
                month = (self.ncname.Mstart \
                         + self.time * dsec//(30*24*3600)-1)%12 + 1
                year = self.ncname.Ystart \
                    + ((self.ncname.Mstart-1)*30*24*3600 + self.time*dsec)\
                        //(360*24*3600)
                infiletime = (self.time * dsec % (30*24*3600))//dsec
                filetime = self.time
                ncfile = "".join([self.ncname.his, 'Y', format(year), 'M',
                                   format(month), self.ncname.fileformat])
            
        elif self.ncname.model == 'croco_xios' \
            or 'croco_gigatl1' in self.ncname.model:
            
            time1 = self.time - self.time%self.ncname.tfile
            date1 = self.ncname.realyear_tstart \
                + timedelta(seconds=time1*self.ncname.dtfile)
            year1 = date1.year
            month1 = date1.month
            day1 = date1.day
            hour1 = date1.hour
            
            time2 = time1 + self.ncname.tfile #- 1

            date2 = self.ncname.realyear_tstart \
                  + timedelta(seconds=time2*self.ncname.dtfile)

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

            infiletime = self.time%self.ncname.tfile
            filetime = self.time
            
            if self.ncname.model == 'croco_xios':
                ncfile = "".join(
                    [self.ncname.his,
                    f"{year1:04}-{month1:02}-{day1:02}-",
                    f"{year2:04}-{month2:02}-{day2:02}{self.ncname.fileformat}"
                    ]
                )
            elif self.ncname.model == 'croco_gigatl1':
                ncfile = "".join(
                    [self.ncname.his,
                    f"{year1:04}-{month1:02}-{day1:02}", self.ncname.fileformat
                    ]
                )                
            elif self.ncname.model == 'croco_gigatl1_hours':
                ncfile = "".join(
                    [self.ncname.his,
                    f"{year1:04}-{month1:02}-{day1:02}-{hour1:02}",
                    self.ncname.fileformat
                    ]
                )      
                
        return filetime, infiletime, ncfile
    
###################################################################################
#   Update time 
###################################################################################
    def update(self, time=None, output =True):

        """
        
        Update 'my_simul' to the next timestep
        (path to files, infiletime, filetime, ncfile, oceandate, date ...etc...)

        """

        #for time in range(self.time0, self.ncname.tend, self.dtime ):
        if time==None:
            self.time += 1
        else:
            self.time = int(time)
        
        self.ncname = files(self.simul, time=self.time, output=output)
        self.filetime, self.infiletime, self.ncfile = self.filedata(
            output=output)

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

        if self.ncname.model in ['ucla']:
            self.oceantime = int(
                np.array(ncfile.variables['ocean_time'][self.infiletime])
                )
        else:
            try:
                self.oceantime = int(
                    np.array(ncfile.variables['scrum_time'][self.infiletime])
                )
            except:
                try:
                    self.oceantime = int(
                        np.array(ncfile.variables['time'][self.infiletime])
                        )
                except:
                    self.oceantime = int(
                        np.array(ncfile.variables['time_centered'][self.infiletime])
                        )


        if not self.ncname.realyear:
        
            self.year = int(np.floor(self.oceantime/(360*24*3600)))
            
            self.oceantime = self.oceantime%(360*24*3600)

            self.month = self.oceantime//(24*3600)//30+1
            month_name = ["None","Jan","Feb","Mar","Apr", "May", "Jun", "Jul",
                          "Aug","Sep","Oct","Nov","Dec"] 

            self.day = self.oceantime//(24*3600) - (self.month-1) * 30 + 1

            self.hour = self.oceantime%(24*3600)//3600

            self.min = self.oceantime%(3600)//60

            self.date = month_name[self.month] + ' ' + '{0:02}'.format(self.day) + ' - ' +'{0:02}'.format(self.hour) + ':' + '{0:02}'.format(self.min) 


        else:
        
            date = self.ncname.realyear_origin + timedelta(days=self.oceantime/3600./24.)
            
            self.year = date.year
            self.month = date.month
            month_name = ["None","Jan","Feb","Mar","Apr", "May", "Jun", "Jul",
                          "Aug","Sep","Oct","Nov","Dec"] 

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
            if self.ncname.model in ['ucla']:
                self.dt = np.array(ncfile.variables['ocean_time'][1]) \
                        - np.array(ncfile.variables['ocean_time'][0])
            # TODO add croco croco_clim, croco_inter
            elif  self.ncname.model in ['croco', 'agrif']:
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
        self.fileformat = '.nc'
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

            folder = '/home/gula/libra/ROMS/Simulations/Rutgers_example/'
            self.grd = folder + 'grd.nc'

            self.his = folder +'velocity_example'
            self.tfile = 2
            self.tstart = 0
            self.tend = 2

        ##################

        # --------------------------------------------------------------------
        # for code testing
        # >>> from R_files import test
        # >>> test()
        # --------------------------------------------------------------------
        # ucla default 'roms_0010.nc'
        elif 'code_test' in simul:
            if 'ucla' in simul:
                folder = '/ucla'
                self.his = folder + '/roms_'
                self.grd = folder + '/roms_grd.nc'
                self.frc = folder + '/roms_frc.nc'
                self.wind = folder + '/roms_frc.nc'
                self.tfile = 5
                self.tstart = 0
                self.tend = 1000

            # croco interannuel 'roms_his_Y20M3.nc'
            elif 'croco' in simul:
                if 'date_inter' in simul:
                    self.model = 'croco_date'
                    self.realyear = True
                    self.realyear_origin = datetime(2007, 1, 1)
                    folder= '/data/croco_inter'
                    self.his = folder + '/roms_his_'
                    self.grd = folder + '/roms_grid.nc'
                    self.frc = self.grd
                    self.wind = self.grd
                    self.tstart = 0
                    self.tend = 1000

                # croco climatology 'roms_his_Y20M3.nc'
                elif 'date_clim' in simul:
                    self.model = 'croco_date'
                    self.realyear = False
                    folder = '/croco_clim'
                    self.his = folder + '/roms_his_'
                    self.grd = folder + '/roms_his_Y20M1.nc'
                    self.frc = folder + '/roms_his_Y20M1.nc'
                    self.wind = folder + '/roms_his_Y20M1.nc'
                    self.tfile = 30
                    self.tstart = 0
                    self.Ystart = 20
                    self.Mstart = 1
                    self.tend = 1000

                # croco xios 'croco_his_2001-01-06-2001-01-10.nc'
                elif 'xios' in simul:
                    self.realyear = True
                    self.realyear_origin = datetime(1999, 1, 1)
                    self.realyear_tstart = datetime(2001, 1, 1)
                    self.model = 'croco_xios'

                    folder = '/croco_xios'
                    self.grd = folder + '/grd.nc'

                    self.his = folder + '/croco_his_'
                    self.dtfile = 12 * 3600
                    self.tfile = 10
                    self.tstart = 0
                    self.tend = 120    

                else:
                    # default croco croco_his_1000.nc
                    self.realyear = True
                    self.realyear_origin = datetime(1999, 1, 1)
                    self.model = 'croco'
                    self.digits = 5

                    folder = '/croco'
                    self.grd = folder +'/grd.nc'
                    self.his = folder +'/polgyr_his.'
                    self.frc = folder + '/polgyr_his.03360.nc'
                    self.wind = folder + '/polgyr_his.03360.nc'
                    self.tfile = 20
                    self.tstart = 0
                    self.tend = 1000

            elif 'fuild2d' in simul:
                pass
        
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

def test():
    """
    test code
    >>> from R_files import test
    >>> test()
    
    """
    
    def _messsage_ok(test_msg, ok_msg="OK\n"):
        print(f"{test_msg}: {ok_msg}")
        return None

    def update():
        """
        test update method
        for each velocity input there is a list of parameters defined in params
    
        to add a new input file format, specify simulation_name according to 
        file() Class, update at time and check filename and infiletime
        
        params:
            [simulation_name, time, filename, infiletime]
        """

        def _test(my_simul, time, fname, infiletime, error_msg):
            simul = load(my_simul, touchfile=False, loadgrid=False, output=False)
            simul.update(time=time, output=False)
            assert simul.ncfile == fname, error_msg
            assert simul.infiletime == infiletime, error_msg
            print(simul)
            _messsage_ok(my_simul)

        # main test.update()
        params = [
            ['code_test_croco_date_inter', 59, '/data/croco_inter/roms_his_Y2007M3.nc', 0],
            ['code_test_croco_date_clim', 60, '/croco_clim/roms_his_Y20M3.nc', 0],
            ['code_test_croco_xios', 10, '/croco_xios/croco_his_2001-01-06-2001-01-10.nc', 0],
            ['code_test_croco', 145, '/croco/polgyr_his.00140.nc', 5],
            ['code_test_ucla', 12, '/ucla/roms_0010.nc', 2],
            ['code_test_rutgers', 1, '/home/gula/libra/ROMS/Simulations/Rutgers_example/velocity_example.nc', 1]
        ]
        error_msg = "Error in R_files test.update()"
        
        print("start update \n")
        # croco inter
        for _my_simul, _time, _fname, _infiletime in params:
            _test(my_simul=_my_simul, time=_time, fname=_fname,
                   infiletime=_infiletime, error_msg=error_msg)

        print('update ok')
        print('------------------')
        return None
    
    def main():
        print(f"Start to test R_files.py methods\n")
        update()
        print(f"test OK\n")
        return None
    
    # main test()
    main()

    return None












