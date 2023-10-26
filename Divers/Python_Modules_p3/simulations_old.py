from __future__ import print_function


#from load_modules import *

#Netcdf IO module
#from Scientific.IO.NetCDF import *
from builtins import range
from builtins import object
from netCDF4 import Dataset

import os

#module for numerics
import numpy as np

#for plotting
import matplotlib.pyplot as py
import matplotlib.colors as col
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#convert seconds to date
import time as tm

#ROMSTOOLS
import romstools_old as roms


###################################################################################

if os.getenv('HOSTNAME')=='shirma.atmos.ucla.edu':
    celtic='/mnt/celtic/gula/' 
    shirma='/shirma/gula/' 
    cherokee='/mnt/cherokee/gula/'
elif  os.getenv('HOSTNAME')=='celtic.atmos.ucla.edu':
    celtic='/celtic/gula/' 
    shirma='/mnt/shirma/gula/' 
    cherokee='/mnt/cherokee/gula/'
elif  os.getenv('HOSTNAME')=='cherokee.atmos.ucla.edu':
    celtic='/mnt/celtic/gula/' 
    shirma='/mnt/shirma/gula/' 
    cherokee='/cherokee/gula/'
else:
    celtic='/mnt/celtic/gula/' 
    shirma='/mnt/shirma/gula/' 
    cherokee='/mnt/cherokee/gula/'

###################################################################################

###################################################################################
#Define files path for each simulation
###################################################################################

class cst(object):

    def __init__(self):
 
        #gravity
        self.g=9.81


###################################################################################
#Define files path for each simulation
###################################################################################


from R_files import files


###################################################################################
#Some pre-defined domains: [ny1,ny2,nx1,nx2,depth1,depth2] at each time step
###################################################################################

def is_static(domainname):

    if domainname in ['filament','filamentZ','shing','filam','filamZ1']:  
        static_domain=0
    else:
        static_domain=0  

    return static_domain

#########################################
    

def domain(ncname0,domainname,time):

    if domainname=='filament': 

        if time<=107: [nx0 ,ny0]=[900, 300]
        if time==108: [nx0 ,ny0]=[960, 300]     
        if time==109: [nx0 ,ny0]=[1000, 300]     
        if time==110: [nx0 ,ny0]=[1040, 300]
        if time==111: [nx0 ,ny0]=[1100, 300]
        if time==112: [nx0 ,ny0]=[1160, 320]
        if time==113: [nx0 ,ny0]=[1210, 340]
        if time==114: [nx0 ,ny0]=[1260, 360]
        if time==115: [nx0 ,ny0]=[1320, 380]
        if time==116: [nx0 ,ny0]=[1360, 400]
        if time==117: [nx0 ,ny0]=[1417, 400]
        if time==118: [nx0 ,ny0]=[1460, 420]
        if time==119: [nx0 ,ny0]=[1500, 450]
        if time==120: [nx0 ,ny0]=[1540, 480]
        if time==121: [nx0 ,ny0]=[1592, 467]
        if time==122: [nx0 ,ny0]=[1640, 510]
        if time==123: [nx0 ,ny0]=[1690, 544]
        if time==124: [nx0 ,ny0]=[1746, 586]
        if time==125: [nx0 ,ny0]=[1788, 628]
        if time==126: [nx0 ,ny0]=[1830, 642]
        if time>=127: [nx0 ,ny0]=[1900+40*(time-127), 677+30*(time-127)]

        dx=200
        dy=300

        #if time<132:
        #    [ny1,ny2,nx1,nx2] = [np.max([0,ny0-dy]),ny0+dy*2/3,np.max([0,nx0-dx]),nx0+dx*2/3]  
        ##else:
        [ny1,ny2,nx1,nx2] = [np.max([0,ny0-dy]),min(ny0+dy,1330),np.max([0,nx0-dx]),min(nx0+dx,2400)]  

        depths=[0]

    elif domainname=='filam': 

        time = 114 + (time-114)/18
        print('newtime')
        print(time)

        if time<115: [nx0 ,ny0]=[1320, 380]
        if time==115: [nx0 ,ny0]=[1320, 380]
        if time==116: [nx0 ,ny0]=[1360, 400]
        if time==117: [nx0 ,ny0]=[1417, 400]
        if time==118: [nx0 ,ny0]=[1460, 420]
        if time==119: [nx0 ,ny0]=[1500, 450]
        if time==120: [nx0 ,ny0]=[1540, 480]
        if time==121: [nx0 ,ny0]=[1592, 467]
        if time==122: [nx0 ,ny0]=[1640, 510]
        if time==123: [nx0 ,ny0]=[1690, 544]
        if time==124: [nx0 ,ny0]=[1746, 586]
        if time==125: [nx0 ,ny0]=[1788, 628]
        if time==126: [nx0 ,ny0]=[1830, 642]
        if time>=127: [nx0 ,ny0]=[1900+40*(time-127), 677+30*(time-127)]

        dx=200
        dy=200

        if time<132:
            [ny1,ny2,nx1,nx2] = [ny0-dy,ny0+dy*2/3,nx0-dx,nx0+dx*2/3]  
        else:
            [ny1,ny2,nx1,nx2] = [ny0-dy,min(ny0+dy,1330),nx0-dx,min(nx0+dx*2,2400)]  

        depths=[0]
        
    elif domainname=='filam_article': 

        time = 114 + (time-114)/18
        print('newtime')
        print(time)

        if time<115: [nx0 ,ny0]=[1320, 380]
        if time==115: [nx0 ,ny0]=[1320, 380]
        if time==116: [nx0 ,ny0]=[1360+20, 400+50]
        if time==117: [nx0 ,ny0]=[1417+30, 400+50]
        if time==118: [nx0 ,ny0]=[1460+40, 420+50]
        if time==119: [nx0 ,ny0]=[1500+50, 450+50]
        if time==120: [nx0 ,ny0]=[1540+50, 480+50]
        if time==121: [nx0 ,ny0]=[1592, 467]
        if time==122: [nx0 ,ny0]=[1640, 510]
        if time==123: [nx0 ,ny0]=[1690, 544]
        if time==124: [nx0 ,ny0]=[1746, 586]
        if time==125: [nx0 ,ny0]=[1788, 628]
        if time==126: [nx0 ,ny0]=[1830, 642]
        if time>=127: [nx0 ,ny0]=[1900+40*(time-127), 677+30*(time-127)]

        dx=200
        dy=250

        if time<132:
            [ny1,ny2,nx1,nx2] = [ny0-dy,ny0+dy*2/3,nx0-dx,nx0+dx*2/3]  
        else:
            [ny1,ny2,nx1,nx2] = [ny0-dy,min(ny0+dy,1330),nx0-dx,min(nx0+dx*2,2400)]  

        depths=[0]
        
        
    elif domainname=='filam_article_zoom': 

        time = 114 + (time-114)/18
        print('newtime')
        print(time)

        if time<115: [nx0 ,ny0]=[1320, 380]
        if time==115: [nx0 ,ny0]=[1320, 380]
        if time==116: [nx0 ,ny0]=[1360+20, 400+50]
        if time==117: [nx0 ,ny0]=[1417+30, 400+50]
        if time==118: [nx0 ,ny0]=[1460+40, 420+50]
        if time==119: [nx0 ,ny0]=[1500+50, 450+50]
        if time==120: [nx0 ,ny0]=[1540+50, 480+50]
        if time==121: [nx0 ,ny0]=[1592, 467]
        if time==122: [nx0 ,ny0]=[1640, 510]
        if time==123: [nx0 ,ny0]=[1690, 544]
        if time==124: [nx0 ,ny0]=[1746, 586]
        if time==125: [nx0 ,ny0]=[1788, 628]
        if time==126: [nx0 ,ny0]=[1830, 642]
        if time>=127: [nx0 ,ny0]=[1900+40*(time-127), 677+30*(time-127)]

        dx=150
        dy=250

        if time<132:
            [ny1,ny2,nx1,nx2] = [ny0-dy,ny0+dy*2/3,nx0-dx,nx0+dx*2/3]  
        else:
            [ny1,ny2,nx1,nx2] = [ny0-dy,min(ny0+dy,1330),nx0-dx,min(nx0+dx*2,2400)]  

        depths=[0]        
        
    elif domainname=='filamlarge': 

        time = 114 + (time-114)/18
        print('newtime')
        print(time)

        if time<115: [nx0 ,ny0]=[1320, 380]
        if time==115: [nx0 ,ny0]=[1320, 380]
        if time==116: [nx0 ,ny0]=[1360, 400]
        if time==117: [nx0 ,ny0]=[1417, 400]
        if time==118: [nx0 ,ny0]=[1460, 420]
        if time==119: [nx0 ,ny0]=[1500, 450]
        if time==120: [nx0 ,ny0]=[1540, 480]
        if time==121: [nx0 ,ny0]=[1592, 467]
        if time==122: [nx0 ,ny0]=[1640, 510]
        if time==123: [nx0 ,ny0]=[1690, 544]
        if time==124: [nx0 ,ny0]=[1746, 586]
        if time==125: [nx0 ,ny0]=[1788, 628]
        if time==126: [nx0 ,ny0]=[1830, 642]
        if time>=127: [nx0 ,ny0]=[1900+40*(time-127), 677+30*(time-127)]

        dx=200
        dy=300

        if time<132:
            [ny1,ny2,nx1,nx2] = [ny0-dy,ny0+dy*2/3,nx0-dx,nx0+dx*2/3]  
        else:
            [ny1,ny2,nx1,nx2] = [ny0-dy,min(ny0+dy,1330),nx0-dx,min(nx0+dx*2,2400)]  

        depths=[0]
        
        
    elif domainname=='filam3d': 

        time = 114 + (time-114)/18
        print('newtime')
        print(time)

        if time<115: [nx0 ,ny0]=[1320, 380]
        if time==115: [nx0 ,ny0]=[1320, 380]
        if time==116: [nx0 ,ny0]=[1360, 400]
        if time==117: [nx0 ,ny0]=[1417, 400]
        if time==118: [nx0 ,ny0]=[1460, 420]
        if time==119: [nx0 ,ny0]=[1500, 450]
        if time==120: [nx0 ,ny0]=[1540, 480]
        if time==121: [nx0 ,ny0]=[1592, 467]
        if time==122: [nx0 ,ny0]=[1640, 510]
        if time==123: [nx0 ,ny0]=[1690, 544]
        if time==124: [nx0 ,ny0]=[1746, 586]
        if time==125: [nx0 ,ny0]=[1788, 628]
        if time==126: [nx0 ,ny0]=[1830, 642]
        if time>=127: [nx0 ,ny0]=[1900+40*(time-127), 677+30*(time-127)]

        dx=200
        dy=200

        if time<132:
            [ny1,ny2,nx1,nx2] = [ny0-dy,ny0+dy*2/3,nx0-dx,nx0+dx*2/3]  
        else:
            [ny1,ny2,nx1,nx2] = [ny0-dy,min(ny0+dy,1330),nx0-dx,min(nx0+dx*2,2400)]  

        depths=[-200,10,5]

    elif domainname=='shing': 

        if time<=100: [nx0 ,ny0]=[850, 760]
        else: [nx0 ,ny0]=[850+20*(time-100), 760]

        dx=400
        dy=200
        [ny1,ny2,nx1,nx2] = [ny0-dy,min(ny0+dy,1330),nx0-dx,min(nx0+dx,2400)]  
        depths=[0]

    elif domainname=='filamentZ': 

        if time<115: [nx0 ,ny0]=[1320, 380]
        if time==115: [nx0 ,ny0]=[1320, 380]
        if time==116: [nx0 ,ny0]=[1360, 400]
        if time==117: [nx0 ,ny0]=[1417, 400]
        if time==118: [nx0 ,ny0]=[1460, 420]
        if time==119: [nx0 ,ny0]=[1500, 440]
        if time==120: [nx0 ,ny0]=[1540, 480]
        if time==121: [nx0 ,ny0]=[1592, 467]
        if time==122: [nx0 ,ny0]=[1640, 510]
        if time==123: [nx0 ,ny0]=[1690, 544]
        if time==124: [nx0 ,ny0]=[1746, 586]
        if time==125: [nx0 ,ny0]=[1788, 628]
        if time==126: [nx0 ,ny0]=[1830, 642]
        if time>=127: [nx0 ,ny0]=[1900+40*(time-127), 677+30*(time-127)]

        #new Z
        if time==120: [nx0 ,ny0]=[1525, 536]
        if time==121: [nx0 ,ny0]=[1575, 560]
        if time==122: [nx0 ,ny0]=[1620, 600]

        #nx1=1450
        #nx2=1550
        #ny1=440-2
        #ny2=440+2
        dx=50

        if time<132:
            [ny1,ny2,nx1,nx2] = [ny0,ny0,nx0-dx,nx0+dx]  
        else:
            [ny1,ny2,nx1,nx2] = [ny0,ny0,nx0-dx,min(nx0+dx*2,2400)]  

        depths=np.arange(-200,1,1)

    elif domainname=='shingZ': 

        nx1=1200
        nx2=1200
        ny1=500
        ny2=1100

        depths=np.arange(-500,1,10)


    elif domainname=='filamZ': 

        time = 114 + (time-114)/18
        print('newtime')
        print(time)

        if time<115: [nx0 ,ny0]=[1320, 380]
        if time==115: [nx0 ,ny0]=[1320, 380]
        if time==116: [nx0 ,ny0]=[1360, 400]
        if time==117: [nx0 ,ny0]=[1417, 400]
        if time==118: [nx0 ,ny0]=[1460, 420]
        if time==119: [nx0 ,ny0]=[1500, 440]
        if time==120: [nx0 ,ny0]=[1540, 480]
        if time==121: [nx0 ,ny0]=[1592, 467]
        if time==122: [nx0 ,ny0]=[1640, 510]
        if time==123: [nx0 ,ny0]=[1690, 544]
        if time==124: [nx0 ,ny0]=[1746, 586]
        if time==125: [nx0 ,ny0]=[1788, 628]
        if time==126: [nx0 ,ny0]=[1830, 642]
        if time>=127: [nx0 ,ny0]=[1900+40*(time-127), 677+30*(time-127)]

        #new Z
        if time==120: [nx0 ,ny0]=[1525, 536]
        if time==121: [nx0 ,ny0]=[1575, 560]
        if time==122: [nx0 ,ny0]=[1620, 600]

        dx=50

        if time<132:
            [ny1,ny2,nx1,nx2] = [ny0,ny0,nx0-dx,nx0+dx]  
        else:
            [ny1,ny2,nx1,nx2] = [ny0,ny0,nx0-dx,min(nx0+dx*2,2400)]  

        depths=np.arange(-200,1,1)



    elif domainname=='filamZ1': 

        coord=np.zeros((500,2),int)

        coord[:196,:]=[1473, 369]
        coord[200,:]=[1483, 369]
        coord[205,:]=[1493, 369]
        coord[210,:]=[1503, 369]
        coord[215,:]=[1512, 369]
        coord[220,:]=[1522, 376]
        coord[225,:]=[1532, 383]
        coord[230,:]=[1543, 390]
        coord[235,:]=[1553, 399]
        coord[240,:]=[1564, 407]
        coord[245,:]=[1574, 409]
        coord[250,:]=[1583, 413]
        coord[255,:]=[1594, 418]
        coord[260,:]=[1602, 425]
        coord[265,:]=[1612, 431]
        coord[270,:]=[1623, 438]
        coord[275,:]=[1632, 445]
        coord[280:,:]=[1640, 452]

        time1=time/5*5
        if time1==time:
            [nx0,ny0]=coord[time,:]
        else:
            time2=time1+5
            deltat=time-time1
            [nx0 ,ny0]=deltat/5.*coord[time2,:]+(5.-deltat)/5.*coord[time1,:]
        dx=50
        [ny1,ny2,nx1,nx2] = [int(ny0)-2,int(ny0)+2,int(nx0-dx),int(nx0+dx)]  
        depths=np.arange(-151,5,5)
        
        #for article
        dx=20
        [ny1,ny2,nx1,nx2] = [int(ny0)-2,int(ny0)+2,int(nx0-dx),int(nx0+dx)]  
        depths=np.arange(-101,5,5)        
        
        
        #depths=[-10]

    elif domainname=='filamZ11d': 

        coord=np.zeros((500,2),int)

        coord[:196,:]=[1473, 369]
        coord[200,:]=[1483, 369]
        coord[205,:]=[1493, 369]
        coord[210,:]=[1503, 369]
        coord[215,:]=[1512, 369]
        coord[220,:]=[1522, 376]
        coord[225,:]=[1532, 383]
        coord[230,:]=[1543, 390]
        coord[235,:]=[1553, 399]
        coord[240,:]=[1564, 407]
        coord[245,:]=[1574, 409]
        coord[250,:]=[1583, 413]
        coord[255,:]=[1594, 418]
        coord[260,:]=[1602, 425]
        coord[265,:]=[1612, 431]
        coord[270,:]=[1623, 438]
        coord[275,:]=[1632, 445]
        coord[280:,:]=[1640, 452]

        time1=time/5*5
        if time1==time:
            [nx0,ny0]=coord[time,:]
        else:
            time2=time1+5
            deltat=time-time1
            [nx0 ,ny0]=deltat/5.*coord[time2,:]+(5.-deltat)/5.*coord[time1,:]
        dx=20
        [ny1,ny2,nx1,nx2] = [int(ny0),int(ny0),int(nx0-dx),int(nx0+dx)]  

        #depths=np.arange(-101,0,5)
        depths=[0]
     
    elif domainname=='filamZ2': 

        coord=np.zeros((500,2),int)
        coord[:201,:]=[1490, 410]
        coord[205,:]=[1500, 430]
        coord[210,:]=[1510, 450]
        coord[215,:]=[1525, 460]
        coord[220,:]=[1535, 468]
        coord[225,:]=[1555, 473]
        coord[230,:]=[1565, 480]
        coord[235,:]=[1580, 490]
        coord[240,:]=[1590, 500]
        coord[245,:]=[1605, 505]
        coord[250,:]=[1620, 515]
        coord[255,:]=[1635, 525]
        coord[260,:]=[1650, 535]
        coord[265,:]=[1660, 545]
        coord[270,:]=[1675, 555]
        coord[275,:]=[1690, 565]
        coord[280:,:]=[1690, 565]

        time1=time/5*5
        if time1==time:
            [nx0,ny0]=coord[time,:]
        else:
            time2=time1+5
            deltat=time-time1
            [nx0 ,ny0]=deltat/5.*coord[time2,:]+(5.-deltat)/5.*coord[time1,:]
        ny0 = ny0
        dx=50
        [ny1,ny2,nx1,nx2] = [int(ny0)-5,int(ny0)+5,int(nx0-dx),int(nx0+dx)]  

        depths=np.arange(-201,5,5)
        
    elif domainname=='filamZ21d': 

        coord=np.zeros((500,2),int)
        coord[:201,:]=[1490, 410]
        coord[205,:]=[1500, 430]
        coord[210,:]=[1510, 450]
        coord[215,:]=[1525, 460]
        coord[220,:]=[1535, 468]
        coord[225,:]=[1555, 473]
        coord[230,:]=[1565, 480]
        coord[235,:]=[1580, 490]
        coord[240,:]=[1590, 500]
        coord[245,:]=[1605, 505]
        coord[250,:]=[1620, 515]
        coord[255,:]=[1635, 525]
        coord[260,:]=[1650, 535]
        coord[265,:]=[1660, 545]
        coord[270,:]=[1675, 555]
        coord[275,:]=[1690, 565]
        coord[280:,:]=[1690, 565]

        time1=time/5*5
        if time1==time:
            [nx0,ny0]=coord[time,:]
        else:
            time2=time1+5
            deltat=time-time1
            [nx0 ,ny0]=deltat/5.*coord[time2,:]+(5.-deltat)/5.*coord[time1,:]
        ny0 = ny0
        dx=50
        [ny1,ny2,nx1,nx2] = [int(ny0),int(ny0),int(nx0-dx),int(nx0+dx)]  

        depths=[0]
        
    elif domainname=='filamtest': 

        coord=np.zeros((500,2),int)
        coord[:191,:]=[1470, 410]
        coord[195,:]=[1480, 410]
        coord[200,:]=[1490, 410]
        coord[205,:]=[1500, 430]
        coord[210,:]=[1510, 450]
        coord[215,:]=[1525, 460]
        coord[220,:]=[1535, 468]
        coord[225,:]=[1555, 473]
        coord[230,:]=[1565, 480]
        coord[235,:]=[1580, 490]
        coord[240,:]=[1590, 500]
        coord[245,:]=[1605, 505]
        coord[250,:]=[1620, 515]
        coord[255,:]=[1635, 525]
        coord[260,:]=[1650, 535]
        coord[265,:]=[1660, 545]
        coord[270,:]=[1675, 555]
        coord[275,:]=[1690, 565]
        coord[280:,:]=[1690, 565]

        time1=time/5*5
        if time1==time:
            [nx0,ny0]=coord[time,:]
        else:
            time2=time1+5
            deltat=time-time1
            [nx0 ,ny0]=deltat/5.*coord[time2,:]+(5.-deltat)/5.*coord[time1,:]
        ny0 = ny0
        dx=50
        [ny1,ny2,nx1,nx2] = [350,700,int(nx0-100),int(nx0+50)]  

        depths=np.arange(-150,1,10)

    elif domainname=='filamtestz0': 

        coord=np.zeros((500,2),int)
        coord[:191,:]=[1470, 410]
        coord[195,:]=[1480, 410]
        coord[200,:]=[1490, 410]
        coord[205,:]=[1500, 430]
        coord[210,:]=[1510, 450]
        coord[215,:]=[1525, 460]
        coord[220,:]=[1535, 468]
        coord[225,:]=[1555, 473]
        coord[230,:]=[1565, 480]
        coord[235,:]=[1580, 490]
        coord[240,:]=[1590, 500]
        coord[245,:]=[1605, 505]
        coord[250,:]=[1620, 515]
        coord[255,:]=[1635, 525]
        coord[260,:]=[1650, 535]
        coord[265,:]=[1660, 545]
        coord[270,:]=[1675, 555]
        coord[275,:]=[1690, 565]
        coord[280:,:]=[1690, 565]

        time1=time/5*5
        if time1==time:
            [nx0,ny0]=coord[time,:]
        else:
            time2=time1+5
            deltat=time-time1
            [nx0 ,ny0]=deltat/5.*coord[time2,:]+(5.-deltat)/5.*coord[time1,:]
        ny0 = ny0
        dx=50
        [ny1,ny2,nx1,nx2] = [300,700,int(nx0-100),int(nx0+50)]  

        depths=[-50]

    elif domainname=='filamZ3': 

        coord=np.zeros((500,2),int)
        coord[:161,:]=[1400-40, 430]
        coord[165,:]=[1410-35, 440]      
        coord[170,:]=[1420-30, 450]      
        coord[175,:]=[1430-25, 460]       
        coord[180,:]=[1440-20, 470]
        coord[185,:]=[1450-15, 480]              
        coord[190,:]=[1460-10, 490]          
        coord[195,:]=[1470-5, 500]             
        coord[200,:]=[1480, 510]        
        coord[205,:]=[1490, 520]
        coord[210,:]=[1500, 530]
        coord[215,:]=[1510, 540]
        coord[220,:]=[1525, 550]
        coord[225,:]=[1535, 560]
        coord[230,:]=[1545, 570]
        coord[235,:]=[1555, 580]
        coord[240,:]=[1565, 590]
        coord[245,:]=[1575, 600]
        coord[250,:]=[1585, 610]
        coord[255,:]=[1595, 620]
        coord[260,:]=[1605, 630]
        coord[265,:]=[1605, 630]
        coord[270,:]=[1605, 630]
        coord[275,:]=[1605, 630]
        coord[280,:]=[1605, 630]
        time1=time/5*5
        if time1==time:
            [nx0,ny0]=coord[time,:]
        else:
            time2=time1+5
            deltat=time-time1
            [nx0 ,ny0]=deltat/5.*coord[time2,:]+(5.-deltat)/5.*coord[time1,:]
        ny0 = ny0
        dx=50
        [ny1,ny2,nx1,nx2] = [int(ny0)-5,int(ny0)+5,int(nx0-dx),int(nx0+dx)]  

        depths=np.arange(-200,2,2)     
        
    elif domainname=='filamZ31d': 

        coord=np.zeros((500,2),int)
        coord[:161,:]=[1400-40, 430]
        coord[165,:]=[1410-35, 440]      
        coord[170,:]=[1420-30, 450]      
        coord[175,:]=[1430-25, 460]       
        coord[180,:]=[1440-20, 470]
        coord[185,:]=[1450-15, 480]              
        coord[190,:]=[1460-10, 490]          
        coord[195,:]=[1470-5, 500]
        coord[200,:]=[1480, 510]        
        coord[205,:]=[1490, 520]
        coord[210,:]=[1500, 530]
        coord[215,:]=[1510, 540]
        coord[220,:]=[1525, 550]
        coord[225,:]=[1535, 560]
        coord[230,:]=[1545, 570]
        coord[235,:]=[1555, 580]
        coord[240,:]=[1565, 590]
        coord[245,:]=[1575, 600]
        coord[250,:]=[1585, 610]
        coord[255,:]=[1595, 620]
        coord[260,:]=[1605, 630]
        coord[265,:]=[1605, 630]
        coord[270,:]=[1605, 630]
        coord[275,:]=[1605, 630]
        coord[280,:]=[1605, 630]
        time1=time/5*5
        if time1==time:
            [nx0,ny0]=coord[time,:]
        else:
            time2=time1+5
            deltat=time-time1
            [nx0 ,ny0]=deltat/5.*coord[time2,:]+(5.-deltat)/5.*coord[time1,:]
        ny0 = ny0
        dx=50
        [ny1,ny2,nx1,nx2] = [int(ny0),int(ny0),int(nx0-dx),int(nx0+dx)]  

        depths=[-10]
        
    elif domainname=='filamZ1bis': 

        coord=np.zeros((500,2),int)
        coord[:201,:]=[1483, 369]
        coord[205,:]=[1493, 369]
        coord[210,:]=[1503, 369]
        coord[215,:]=[1512, 369]
        coord[220,:]=[1522, 376]
        coord[225,:]=[1532, 383]
        coord[230,:]=[1543, 390]
        coord[235,:]=[1553, 399]
        coord[240,:]=[1564, 407]
        coord[245,:]=[1574, 409]
        coord[250,:]=[1583, 413]
        coord[255,:]=[1594, 418]
        coord[260,:]=[1602, 425]
        coord[265,:]=[1612, 431]
        coord[270,:]=[1623, 438]
        coord[275,:]=[1632, 445]
        coord[280:,:]=[1640, 452]

        time1=time/5*5
        if time1==time:
            [nx0,ny0]=coord[time,:]
        else:
            time2=time1+5
            deltat=time-time1
            [nx0 ,ny0]=deltat/5.*coord[time2,:]+(5.-deltat)/5.*coord[time1,:]
        ny0 = ny0 +23
        dx=50
        [ny1,ny2,nx1,nx2] = [int(ny0)-5,int(ny0)+5,int(nx0-dx),int(nx0+dx)]  

        depths=np.arange(-200,1,5)     


    elif domainname=='filament3d':
        [ny1,ny2,nx1,nx2] = [1,100,1,102]  
        depths=np.arange(-200,1,5,int)


    elif domainname=='weddy':

        dx=170
        dy=120
        [ny1,ny2,nx1,nx2] = [1500-dy,1500+dy,950-dx,950+dx]  
        depths=[0]

    elif domainname=='weddyZ':

        dx=0
        dy=120
        [ny1,ny2,nx1,nx2] = [1500-dy,1500+dy,950-dx,950+dx]  
        depths=np.arange(-300,1,10)


    elif domainname=='all':

        ncfile0  = Dataset(ncname0, 'r', format='NETCDF3_CLASSIC')
        depths=[0]
        [ny1,ny2,nx1,nx2] = [0,len(ncfile0.dimensions['eta_rho']),0,len(ncfile0.dimensions['xi_rho'])]  
        ncfile0.close()

    elif domainname=='allz':

        ncfile0  = Dataset(ncname0, 'r', format='NETCDF3_CLASSIC')
        depths=np.arange(1,100,1)
        [ny1,ny2,nx1,nx2] = [0,len(ncfile0.dimensions['eta_rho']),0,len(ncfile0.dimensions['xi_rho'])]  
        ncfile0.close()

    elif domainname=='subvort': 

        coord=np.zeros((500,2),int)
        coord[:140,:]=[542, 370]
        coord[145,:]=[542, 370]
        coord[150,:]=[542, 370]
        coord[155,:]=[540, 360]
        coord[160,:]=[573, 349]
        coord[165,:]=[600, 358]
        coord[170,:]=[623, 350]
        coord[175,:]=[635, 346]
        coord[180,:]=[632, 347]
        coord[185,:]=[627, 353]
        coord[190,:]=[639, 363]
        coord[195,:]=[672, 377]
        coord[200,:]=[705, 387]
        coord[205,:]=[735, 391]
        coord[210,:]=[756, 390]
        coord[215,:]=[772, 399]
        coord[220,:]=[791, 382]
        coord[225,:]=[795, 366]
        coord[230,:]=[765, 356]
        coord[235,:]=[719, 360]
        coord[240,:]=[705, 371]
        coord[245:,:]=[716, 387]

        time1=time/5*5
        if time1==time:
            [nx0,ny0]=coord[time,:]
        else:
            time2=time1+5
            deltat=time-time1
            [nx0 ,ny0]=deltat/5.*coord[time2,:]+(5.-deltat)/5.*coord[time1,:]

        #ny0 = ny0 +23
        dx=100
        dy=100
        [ny1,ny2,nx1,nx2] = [int(ny0-dy),int(ny0+dy),int(nx0-dx),int(nx0+dx)]  

        #depths=np.arange(1,50,1)     
        depths=np.arange(-3000,-800,50)   
        
        
    elif domainname=='shing_article': 
        
        coord=np.zeros((500,2),int)
        coord[:66,:]=[295, 734]            
        coord[70,:]=[355, 750]     
        coord[75,:]=[403, 736]
        coord[80,:]=[521, 730]      
        coord[85,:]=[613, 736]
        coord[90,:]=[717, 722]
        coord[95,:]=[771, 720]
        coord[100,:]=[829, 714]
        coord[105,:]=[953, 716]
        coord[110,:]=[1041, 726]
        coord[115,:]=[1135, 740]
        coord[120,:]=[1217, 758]
        coord[125,:]=[1321, 748]
        coord[130,:]=[1423, 770]
        coord[135,:]=[1509, 782]
        coord[140:,:]=[1615, 800]

        time1=time/5*5

        if time1==time:
            [nx0,ny0]=coord[time,:]
        else:
            time2=time1+5
            deltat=time-time1
            [nx0 ,ny0]=deltat/5.*coord[time2,:]+(5.-deltat)/5.*coord[time1,:]

        
        dx=300
        dy=300
        [ny1,ny2,nx1,nx2] = [max(int(ny0-dy+100),0),min(int(ny0+dy+100),1330),max(int(nx0-dx),0),min(int(nx0+dx),2400)]  
        depths=[0]
 
###################################################################################

        
    elif domainname=='slope_article': 

        time = 60 + (time-60)/3
        print('time in SHING is', time)
        
        coord=np.zeros((500,2),int)
        coord[:66,:]=[295, 734]            
        coord[70,:]=[355, 750]     
        coord[75,:]=[403, 736]
        coord[80,:]=[521, 730]      
        coord[85,:]=[613, 736]
        coord[90,:]=[717, 722]
        coord[95,:]=[771, 720]
        coord[100,:]=[829, 714]
        coord[105,:]=[953, 716]
        coord[110,:]=[1041, 726]
        coord[115,:]=[1135, 740]
        coord[120,:]=[1217, 758]
        coord[125,:]=[1321, 748]
        coord[130,:]=[1423, 770]
        coord[135,:]=[1509, 782]
        coord[140:,:]=[1615, 800]

        time1=time/5*5

        if time1==time:
            [nx0,ny0]=coord[time,:]
        else:
            time2=time1+5
            deltat=time-time1
            [nx0 ,ny0]=deltat/5.*coord[time2,:]+(5.-deltat)/5.*coord[time1,:]
        
        dx=300
        dy=300
        
        [ny1,ny2,nx1,nx2] = [max(int(ny0-dy+100),0),min(int(ny0+dy+100),1330),max(int(nx0-dx),0),min(int(nx0+dx),2400)]  
        depths=[0]
 
###################################################################################

    else:

        ncfile0  = Dataset(ncname0, 'r', format='NETCDF3_CLASSIC')
        [ny1,ny2,nx1,nx2] = eval(domainname)[0:4]

        if len(eval(domainname)[4])==1:
            depths=eval(domainname)[4]
        else:
            depths=np.arange(eval(domainname)[4][0],eval(domainname)[4][1]+eval(domainname)[4][2],eval(domainname)[4][2])



        nx2 = int(np.min([nx2,len(ncfile0.dimensions['xi_rho'])]))
        ny2 = int(np.min([ny2,len(ncfile0.dimensions['eta_rho'])]))



        ncfile0.close()     


###################################################################################


    if ny1==ny2: ny1=ny1-1; ny2=ny2+2
    if nx1==nx2: nx1=nx1-1; nx2=nx2+2


    #depths=[0]
    return [ny1,ny2,nx1,nx2,depths]


###################################################################################
#plot paramters
###################################################################################

def plotconfig(domainname):

    subdomainname = None

    if domainname=='filament':  
        latlon=2
        subdomainname = 'filamentZ'

    elif domainname=='filam':  
        latlon=0
        subdomainname = 'filamZ'

    else:
        latlon=0
        #subdomainname = 'filamentZ' 

    return [latlon,subdomainname]


###################################################################################
#Load grd variables
###################################################################################

def variables_grd(ncname0,ny1=0,ny2=-1,nx1=0,nx2=-1):
        
    ncfile0  = Dataset(ncname0, 'r', format='NETCDF3_CLASSIC')

    mask = np.array(ncfile0.variables['mask_rho'][ny1:ny2,nx1:nx2])
    topo = np.array(ncfile0.variables['h'][ny1:ny2,nx1:nx2])

    topo[mask==0] = 0.
    mask[mask==0] = np.nan

    pm = np.array(ncfile0.variables['pm'][ny1:ny2,nx1:nx2])
    pn = np.array(ncfile0.variables['pn'][ny1:ny2,nx1:nx2])
    f = np.array(ncfile0.variables['f'][ny1:ny2,nx1:nx2])

    if 'lon_rho' in list(ncfile0.variables.keys()):
        lon = np.array(ncfile0.variables['lon_rho'][ny1:ny2,nx1:nx2])
        lat = np.array(ncfile0.variables['lat_rho'][ny1:ny2,nx1:nx2])
    else:
        lon = np.array(ncfile0.variables['x_rho'][ny1:ny2,nx1:nx2])
        lat = np.array(ncfile0.variables['y_rho'][ny1:ny2,nx1:nx2])    
   
    ncfile0.close()

    return [topo,pm,pn,f,lat,lon]


###################################################################################
#long name for each variable (used as plot title)
###################################################################################

def title(varname):

    if varname=='w': name='Vertical Velocity (m/s)'
    elif varname=='temp': name='Temperature (C)'
    elif varname=='u': name='Zonal Velocity (m/s)'
    elif varname=='uper': name=r'$u - \overline{u} \;(m s^{-1})$'   
    elif varname=='v': name=' Meridional Velocity (m/s)'
    elif varname=='pv': name='Potential Vorticity'
    elif varname=='pv1': name='[f + (dv/dx - du/dy)]*db/dz'
    elif varname=='pv2': name='-(dv/dz)*(db/dx)'
    elif varname=='pv3': name='(du/dz)*(db/dy)'
    elif varname=='vrt': name='Relative Vorticity'
    elif varname=='div': name='Horizontal Divergence'
    elif varname=='str': name='Horizontal Strain'
    elif varname=='rho': name='Density'
    elif varname=='buoy': name='Buoyancy'
    elif varname=='zeta': name='Free-surface elevation (m)'
    elif varname=='J2': name='frictional flux of PV'
    elif varname=='ow': name=' '
    elif varname=='ageo': name='A-S'
    elif varname=='stretching': name='Stretching'    
    else: name=varname

    return name



###################################################################################
#Load a variable
###################################################################################

def variables_fast(ncfile,varname,ny1,ny2,nx1,nx2,depths,infiletime,\
                topo=None,pm=None,pn=None,f=None,u=None,v=None,\
                uwinds=None,vwinds=None,cubic=0,ncname=None):

    verbo = False

    

    #if psi grid
    if varname in ['pv','pv1','pv2','pv3','upv','vpv','pvF']: corgrid=[1,1]; wgrid=1
    #if psi grid
    elif varname in ['vrt']: corgrid=[1,1]; wgrid=0
    #if u grid
    elif varname in ['u','dxbuoy','dzu']: corgrid=[0,1]; wgrid=1
    #if v grid
    elif varname in ['v','dybuoy','dzv']: corgrid=[1,0]; wgrid=1
    #if w grid     
    elif varname in ['stretching']: corgrid=[0,0]; wgrid=1    
    #if rho grid         
    else: corgrid=[0,0]; wgrid=0
    
    print('corgrid', corgrid)
    if np.max(depths)<=0: wgrid=0

    var = np.zeros((np.max([len(depths)-wgrid,1]),ny2-ny1-corgrid[0],nx2-nx1-corgrid[1]))*np.nan

    nchunk = int(np.max([np.sqrt((ny2-ny1)*(nx2-nx1)/50000),1]))
    
    if varname in ['ubgstr','vbgstr','stretching']: nchunk=1
    
    print('nchunk',nchunk)

    for i in range(nchunk):
        if i==0: dx1=0 
        else: dx1=2

        if i==nchunk-1: dx2=0 
        else: dx2=2
        nx1i,nx2i = nx1+i*(nx2-nx1)//nchunk-2*dx1,nx1+(i+1)*(nx2-nx1)//nchunk+2*dx2

        for j in range(nchunk):
            if j==0: dy1=0 
            else: dy1=10

            if j==nchunk-1: dy2=0
            else: dy2=10
            ny1i,ny2i = ny1+j*(ny2-ny1)//nchunk-2*dy1,ny1+(j+1)*(ny2-ny1)//nchunk+2*dy2

            if verbo: print('i,j',i,j,dx1,dx2,dy1,dy2)
            
            if u==None:
                vartemp = np.squeeze(variables(ncfile,varname,ny1i,ny2i,nx1i,nx2i,depths,infiletime,\
                                    topo=topo[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1],\
                                    pm=pm[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1],\
                                    pn=pn[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1],\
                                    f=f[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1],\
                                    cubic=1,ncname=ncname,verbo=verbo))
            else:
                #print 'in variables_fast shapes are:', u.shape, v.shape
                vartemp = np.squeeze(variables(ncfile,varname,ny1i,ny2i,nx1i,nx2i,depths,infiletime,\
                                    topo=topo[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1],\
                                    pm=pm[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1],\
                                    pn=pn[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1],\
                                    f=f[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1],\
                                    u=u[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1-1],\
                                    v=v[ny1i-ny1:ny2i-ny1-1,nx1i-nx1:nx2i-nx1],\
                                    cubic=1,ncname=ncname,verbo=verbo))
                #print 'in variables_fast shapes are:', u[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1-1].shape, v[ny1i-ny1:ny2i-ny1-1,nx1i-nx1:nx2i-nx1].shape
            
            #print ' '         
            #print 'i, j,dy1,dx1, vartemp[dy1,dx1]'            
            #print i, j,dy1,dx1, vartemp[dy1,dx1]
            #print ' '

            
            if (np.ndim(vartemp)==2):
                var[:,ny1i+dy1-ny1:ny2i-dy2-ny1-corgrid[0],nx1i+dx1-nx1:nx2i-dx2-nx1-corgrid[1]] = \
                                    vartemp[dy1:ny2i-ny1i-dy2-corgrid[0],dx1:nx2i-nx1i-dx2-corgrid[1]]
            else:
                var[:,ny1i+dy1-ny1:ny2i-dy2-ny1-corgrid[0],nx1i+dx1-nx1:nx2i-dx2-nx1-corgrid[1]] = \
                                    vartemp[:,dy1:ny2i-ny1i-dy2-corgrid[0],dx1:nx2i-nx1i-dx2-corgrid[1]]


    return var 




'''
            test = sim.variables(ncfile,'temp',\
                                    ny1i,ny2i,nx1i,nx2i,\
                                    depths,infiletime,\
                                    topo=topo[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1],\
                                    pm=pm[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1],\
                                    pn=pn[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1],\
                                    f=f[ny1i-ny1:ny2i-ny1,nx1i-nx1:nx2i-nx1],cubic=1,ncname=ncname)
'''




###################################################################################
#Load a variable
###################################################################################

def variables(ncfile,varname,ny1,ny2,nx1,nx2,depths,infiletime,\
                topo=None,pm=None,pn=None,f=None,u=None,v=None,\
                uwinds=None,vwinds=None,cubic=0,ncname=None,mask=None,\
                verbo=False):

    start = tm.time()

    
    #when topo<1m put nan
    depthlim = 0.
    if mask==None: mask = np.ones(topo.shape); mask[topo==depthlim] = 0.


    #check if netcdf file is already opened
    if isinstance(ncfile, str): 
        if verbo: print('open netcdf')
        ncfile  = Dataset(ncfile, 'r', format='NETCDF3_CLASSIC'); opened=True
    else: 
        opened = False




    #if u/v grid:
    nx1var=nx1; ny1var=ny1
    nx2var=nx2; ny2var=ny2

    if varname in ['u','uper']: 
        if ncfile.variables['u'].shape == ncfile.variables['zeta'].shape:
            nx1var=nx1+1
        else:
            nx2var=nx2-1

    if varname=='v': ny2var=ny2-1


    #check if variable is already in netcdffile or needs to be computed  
    if varname in list(ncfile.variables.keys()):

        if (len(depths)==1) and (depths[0]==0) and ('s_rho' not in list(ncfile.dimensions.keys())):
            var=np.array(ncfile.variables[varname][infiletime,ny1var:ny2var,nx1var:nx2var])
            var[topo[0:ny2var-ny1var,0:nx2var-nx1var]==depthlim]=np.nan
            if verbo: print('loading ' + varname + 'from Z0 file')


        #elif varname in ['zeta','pm','pn','f','angle','h','mask_rho','hbls']:
        elif len(ncfile.variables[varname].shape)==3:
            var=np.array(ncfile.variables[varname][infiletime,ny1var:ny2var,nx1var:nx2var]) 
            var[topo[ny1var:ny2var,nx1var:nx2var]==depthlim]=np.nan
            if verbo: print('loading 2D field ' + varname) 


        elif (len(depths)==1) and (depths[0]==0):
            var=np.array(ncfile.variables[varname][infiletime,-1,ny1var:ny2var,nx1var:nx2var]) 
            var[topo[0:ny2var-ny1var,0:nx2var-nx1var]==depthlim]=np.nan
            if verbo: print('loading surface ' + varname) 


        elif min(depths)>0:
            if (len(depths)==1):
                var=np.array(ncfile.variables[varname][infiletime,int(depths[0])-1,ny1var:ny2var,nx1var:nx2var]) 
            else:
                if ncfile.variables[varname].shape[1]==len(ncfile.dimensions['eta_rho']):
                    depth2=int(depths[-1])+1
                else:
                    depth2=int(depths[-1])+1
                var=np.array(ncfile.variables[varname][infiletime,int(depths[0])-1:depth2-1,ny1var:ny2var,nx1var:nx2var]) 
            if verbo: print('loading ' + varname)


        else:

            if verbo: print('interpolating ' + varname) 
            z_r=np.array(roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2,ncname=ncname))
            z_w=np.array(roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2,ncname=ncname))

            #old version
            #var=roms.vinterps(np.array(ncfile.variables[varname][infiletime,:,ny1:ny2var,nx1:nx2var]),z_r,depths,topo,cubic)


            #Fortran version (from alex)
            var=roms.vinterpsF(np.array(ncfile.variables[varname][infiletime,:,ny1var:ny2var,nx1var:nx2var]),depths,z_r,z_w,mask)
            var[:,topo[0:ny2var-ny1var,0:nx2var-nx1var]==depthlim]=np.nan


    else: 

        if verbo: print('computing ', varname) 
        var=get_var(ncfile,varname,ny1,ny2,nx1,nx2,depths,infiletime,topo,pm,pn,f,u,v,uwinds,vwinds,ncname,mask=mask,verbo=verbo)
        if verbo: print(varname, ' computed ', tm.time() - start)



    if opened==True: ncfile.close()


    
    return var 




###################################################################################
#Compute a variable
###################################################################################


def get_var(ncfile,varname,ny1,ny2,nx1,nx2,depths,infiletime,\
            topo=None,pm=None,pn=None,f=None,u=None,v=None,uwinds=None,vwinds=None,ncname=None,\
            verbo=True,mask=None):



    if isinstance(ncfile, str): ncfile  = Dataset(ncfile, 'r', format='NETCDF3_CLASSIC'); opened=True
    else: opened = False



    ################################################################################### 
    if varname in ['absvrt','vrt','div','vrtdiv','str','ageo','ow','ubgstr','vbgstr']:

        print('in get_var depths is',depths)
        
        if u==None: u = variables(ncfile,'u',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo,ncname=ncname)
        if v==None: v = variables(ncfile,'v',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo,ncname=ncname)


        if varname=='absvrt':  var = f+roms.psi2rho(roms.vort(u,v,pm,pn))
        #if varname=='vrt':  var = roms.psi2rho(roms.vort(u,v,pm,pn))/f
        if varname=='vrt':  var = roms.vort(u,v,pm,pn)
        if varname=='div':	var = roms.divs(u,v,pm,pn)/f
        if varname=='vrtdiv':  var = roms.psi2rho(roms.vort(u,v,pm,pn))/roms.divs(u,v,pm,pn)
        if varname=='str':	var = roms.strains(u,v,pm,pn)/f
        if varname=='ageo':	var = roms.ageos(u,v,pm,pn,f)/f
        if varname=='ow':	

            var = roms.ows(u,v,pm,pn) #/f**2

            '''
            s1 = (roms.strains(u,v,pm,pn))**2 
            vrt1 = (roms.psi2rho(roms.vort(u,v,pm,pn)))**2
            #den = s1 + vrt1; #den[den<1e-3]=np.nan
            var = (s1 - vrt1) #/den
            #var=vrt1
            '''

        if varname in ['ubgstr']:	

            var,v = roms.bgstrain(u,v,pm,pn)
            
        if varname in ['vbgstr']:    

            u,var = roms.bgstrain(u,v,pm,pn)
            
            
    ################################################################################### 
    elif varname in ['ug','vg']:
        
        g=cst().g
        zeta = variables(ncfile,'zeta',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo) 

        if varname=='ug':
            var = -g*roms.v2rho(roms.diffy(zeta,pn))/f

        elif varname=='vg':        
            var = g*roms.u2rho(roms.diffx(zeta,pm))/f;
            
    ################################################################################### 
    elif varname in ['ow_alter','ow_alter1','ow_alter2']:
        
        if u==None: u = variables(ncfile,'u',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo) 
        if v==None: v = variables(ncfile,'v',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo)
        
        if varname=='ow_alter': var = roms.ows_alter(u,v,pm,pn)
        if varname=='ow_alter1': var = roms.ows_alter(u,v,pm,pn,1)
        if varname=='ow_alter2': var = roms.ows_alter(u,v,pm,pn,2)   
          
    ################################################################################### 
    elif varname in ['uwinds','vwinds','uvwinds']:

        if uwinds==None: [uwinds,vwinds] = roms.get_winds(ncfile,infiletime,ncname.wind,ny1,ny2,nx1,nx2)

        if varname=='uwinds':
            var = uwinds
        elif varname=='vwinds':
            var = vwinds
        elif varname=='uvwinds':
            var = np.sqrt(roms.u2rho(uwinds**2) + roms.v2rho(vwinds**2))




    ################################################################################### 
    elif varname in ['dsalt']:	

        var0 = variables(ncfile,'salt',ny1,ny2,nx1,nx2,[-200,-100,0],infiletime,topo,verbo=verbo) 

        var = (var0[2,:,:] - var0[1,:,:])*(var0[1,:,:] - var0[0,:,:])
        
    ################################################################################### 
    elif varname in ['z']:	

        if min(depths)>0: #sigma-levels
            var,var2,var3 = np.mgrid[int(depths[0]):int(depths[-1])+1,0:ny2-ny1,0:nx2-nx1]

        else: #interpolated to z-levels   
   
            var1,var2,var3 = np.mgrid[0:len(ncfile.dimensions['s_rho']),0:ny2-ny1,0:nx2-nx1]
            z_r=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)
            var = roms.vinterps(var1,z_r,depths,topo,0)


    ################################################################################### 
    elif varname in ['zw']:	

        if min(depths)>=0: #sigma-levels
            var,var2,var3 = np.mgrid[int(depths[0]):int(depths[-1])+2,0:ny2-ny1,0:nx2-nx1]
        else: #interpolated to z-levels
            var1,var2,var3 = np.mgrid[0:len(ncfile.dimensions['s_rho'])+1,0:ny2-ny1,0:nx2-nx1]
            z_w=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)
            var = roms.vinterps(var1,z_w,depths,topo,0)       

    ################################################################################### 
    elif varname in ['w']:	

        if len(depths)==1 and depths[0]>0:

            if depths[0]==1: depths[0]=2; print('setting depths to [2] ([1] in real) cause we are on w-grid')

            z_r=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[:depths[0]-1,:,:]
            z_w=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[:depths[0],:,:]

            u=np.array(ncfile.variables['u'][infiletime,:depths[0]-1,ny1:ny2,nx1:nx2-1])
            v=np.array(ncfile.variables['v'][infiletime,:depths[0]-1,ny1:ny2-1,nx1:nx2])

            [w] = roms.getw(u,v,pm,pn,z_r,z_w)
            if verbo: print('w.shape', w.shape)
            var = w[-1,:,:]

        else:

            z_r=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)
            z_w=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)

            u=np.array(ncfile.variables['u'][infiletime,:,ny1:ny2,nx1:nx2-1])
            v=np.array(ncfile.variables['v'][infiletime,:,ny1:ny2-1,nx1:nx2])

            [w] = roms.getw(u,v,pm,pn,z_r,z_w)

            if min(depths)>0: #sigma-levels
                var=w[int(depths[0])-1:int(depths[-1]),:,:]
            elif len(depths)==1 and min(depths)==0: #surface
                var=w[-1,:,:]
            else: #interpolated to z-levels   
                var=roms.vinterps(w,z_w,depths,topo)
                


    ################################################################################### 
    elif varname in ['wb']:	

        z_r=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)
        z_w=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)

        u=np.array(ncfile.variables['u'][infiletime,:,ny1:ny2,nx1:nx2-1])
        v=np.array(ncfile.variables['v'][infiletime,:,ny1:ny2-1,nx1:nx2])

        temp = variables(ncfile,'temp',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo) 
        salt = variables(ncfile,'salt',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo)

        [w] = roms.getw(u,v,pm,pn,z_r,z_w)

        rho0=ncfile.rho0; g=cst().g     

        if min(depths)>0: #sigma-levels
            varw=w[int(depths[0])-1:int(depths[-1]),:,:]
            zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[int(depths[0])-1:int(depths[-1]),:,:]


        else: #interpolated to z-levels   
            varw=roms.vinterps(w,z_w,depths,topo)
            zr=np.zeros((len(depths),ny2-ny1,nx2-nx1))
            for iz in range(0,len(depths)): zr[iz,:,:]=depths[iz]
            if (len(depths)==1) and (np.ndim(zr)==3): zr=zr[0,:,:]

        varb = -g*(roms.rho_eos(temp,salt,zr,g,rho0)/rho0-1)
        var = varw * varb


    ################################################################################### 

    elif varname in ['bvf','pv','dxpv','dypv','pvold','pv1','pv2','pv3','dzu','dzv','dzbuoy','upv','vpv','dzt','stretching']:	

        ####################################
        #Define levels where var will be computed
        ####################################

        if (len(depths)==1) and (max(depths)<=0):

            dz=2
            if (depths[0]==0): depths=[-2.5]; depths_w=[depths[0]-dz,depths[0]+dz];
            elif (depths[0]<0): depths_w=[depths[0]-dz,depths[0]+dz];

            zr=np.zeros((len(depths_w),ny2-ny1,nx2-nx1))
            zw=np.zeros((len(depths_w)+1,ny2-ny1,nx2-nx1))

            for iz in range(0,len(depths_w)): zr[iz,:,:]=depths_w[iz]; 
            for iz in range(1,len(depths)+1): zw[iz,:,:]=depths[iz-1]; 
            zw[0,:,:]=depths[0]-2*(depths[0]-depths_w[0])
            zw[-1,:,:]=depths[-1]-2*(depths[-1]-depths_w[-1])
            if verbo: print(varname + ' will be computed at depths', depths)


        elif (len(depths)==1) and (depths[0]>0):

            depths_w=[depths[0],depths[0]+1];
            zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[[int(depths[0])-1,int(depths[0])],:,:]
            zw=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[[int(depths[0])-1,int(depths[0]),int(depths[0])+1],:,:]       
            if verbo: print(varname + ' will be computed at w-level', depths[0])

            #depths_w=[depths[0],depths[0]+1,depths[0]+2];
            #zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[[int(depths[0])-1,int(depths[0]),int(depths[0])+1],:,:]
            #zw=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[[int(depths[0])-1,int(depths[0]),int(depths[0])+1,int(depths[0])+2],:,:]       
            #if verbo: print varname + ' will be computed at w-level', depths[0]

        elif min(depths)>0:

            print ('needs to be checked...')
            depths_w=depths
            zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[int(depths[0])-1:int(depths[-1]),:,:]
            zw=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[int(depths[0])-1:int(depths[-1])+1,:,:]
            if verbo: print(varname + ' will be computed at w-levels', depths)

        else:

            depths_w=np.zeros((len(depths)+1))
            depths_w[0]=depths[0]+0.5*(depths[0]-depths[1])
            depths_w[1:-1]=0.5*(depths[1:]+depths[:-1])
            depths_w[-1]=depths[-1]+0.5*(depths[-1]-depths[-2])
            if verbo: print(depths_w)
            zr=np.zeros((len(depths_w),ny2-ny1,nx2-nx1))
            zw=np.zeros((len(depths_w)+1,ny2-ny1,nx2-nx1))

            for iz in range(0,len(depths_w)): zr[iz,:,:]=depths_w[iz]; 
            for iz in range(1,len(depths)+1): zw[iz,:,:]=depths[iz-1]; 
            zw[0,:,:]=depths[0]-2*(depths[0]-depths_w[0])
            zw[-1,:,:]=depths[-1]-2*(depths[-1]-depths_w[-1])



        ####################################

        rho0=ncfile.rho0; g=cst().g

        if varname in ['dzv']:

            v = variables(ncfile,'v',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo,ncname=ncname)
            var = (v[1:,:,:] - v[:-1,:,:])/roms.rho2v(1*(zr[1:,:,:]-zr[:-1,:,:]))

        elif varname in ['dzu']:
            if verbo: print('dzu')
            u = variables(ncfile,'u',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo,ncname=ncname)
            var = (u[1:,:,:] - u[:-1,:,:])/roms.rho2u(1*(zr[1:,:,:]-zr[:-1,:,:]))

        elif varname in ['bvf']:

            temp_w = variables(ncfile,'temp',ny1,ny2,nx1,nx2,depths_w,infiletime,topo, verbo=verbo,ncname=ncname)
            salt_w = variables(ncfile,'salt',ny1,ny2,nx1,nx2,depths_w,infiletime,topo, verbo=verbo,ncname=ncname)
            [rho,bvf] = roms.rho_eos(temp_w,salt_w,zr,g,rho0,zw)
            var=bvf[1:-1,:,:]/f


        elif varname in ['dzbuoy']:

            temp_w = variables(ncfile,'temp',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo,ncname=ncname)
            salt_w = variables(ncfile,'salt',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo,ncname=ncname)
            print('Should we use rho or rho1?')
            rho = roms.rho1_eos(temp_w,salt_w,zr,g,rho0,zw)
            
            #var=bvf[1:-1,:,:]      
            var0 = -g*(rho/rho0-1)
            var = (var0[1:,:,:] - var0[:-1,:,:])/(zr[1:,:,:]-zr[:-1,:,:])
            
        #elif varname in ['stretching']: (unfinished)

            #temp_w = variables(ncfile,'temp',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo) 
            #salt_w = variables(ncfile,'salt',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo)
            
            #print 'Should we use rho or rho1?'
            #rho = roms.rho1_eos(temp_w,salt_w,zr,g,rho0,zw)
            #rho_mean = tools.nanmean(tools.nanmean(rho,axis=0),axis=0)
            #print 'rho_mean.shape',rho_mean.shape
            
            
            
            ##var=bvf[1:-1,:,:]      
            #var0 = -g*(rho/rho0-1)
            #var = (var0[1:,:,:] - var0[:-1,:,:])/(zr[1:,:,:]-zr[:-1,:,:])
            
        elif varname in ['dzt']:
 
            var0 = variables(ncfile,'temp',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo) 
            var = (var0[1:,:,:] - var0[:-1,:,:])/(zr[1:,:,:]-zr[:-1,:,:])
           
            
            
        else:

            temp_w = variables(ncfile,'temp',ny1,ny2,nx1,nx2,depths_w,infiletime,topo, verbo=verbo,ncname=ncname)
            salt_w = variables(ncfile,'salt',ny1,ny2,nx1,nx2,depths_w,infiletime,topo, verbo=verbo,ncname=ncname)
            u_w = variables(ncfile,'u',ny1,ny2,nx1,nx2,depths_w,infiletime,topo, verbo=verbo,ncname=ncname)
            v_w = variables(ncfile,'v',ny1,ny2,nx1,nx2,depths_w,infiletime,topo, verbo=verbo,ncname=ncname)

            
            if varname=='pv':

                if (depths[0]>0) and (len(depths)>1): [var,pv1,pv2,pv3] = roms.PV_sig(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)
                elif (depths[0]>0): var = roms.PV_sig_bot(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)
                else: var = roms.PV(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)


            elif varname=='dxpv':
                var0 = roms.PV(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)
                var = roms.diffx(var0,roms.rho2psi(pm))

            elif varname=='dypv':
                var0 = roms.PV(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)
                var = roms.diffy(var0,roms.rho2psi(pn))

            elif varname=='pvold':
                [var,pv1,pv2,pv3] = roms.PVF(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)

            elif varname=='pv1':
                if (depths[0]>0) and (len(depths)>1): [pv,var,pv2,pv3] = roms.PV_sig(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)
                else: [pv,var,pv2,pv3] = roms.PVF(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)

            elif varname=='pv2':
                if (depths[0]>0) and (len(depths)>1): [pv,pv1,var,pv3] = roms.PV_sig(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)
                else: [pv,pv1,var,pv3] = roms.PVF(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)

            elif varname=='pv3':
                if (depths[0]>0) and (len(depths)>1): [pv,pv1,pv2,var] = roms.PV_sig(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)
                else: [pv,pv1,pv2,var] = roms.PVF(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)

            elif varname=='upv':
                pv = roms.PV(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)
                var = pv * roms.rho2v(variables(ncfile,'u',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo))

            elif varname=='vpv':
                pv = roms.PV(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)
                var = pv * roms.rho2u(variables(ncfile,'v',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo))


    ################################################################################### 

    elif varname in ['rho','rho1','buoy','dxbuoy','dybuoy','dxbuoy2','dybuoysmooth','tendency','tendency_3d','tendency_geo','tendency_ageo','strg','tendD','tendE','tendF']:	

        temp = variables(ncfile,'temp',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo) 
        salt = variables(ncfile,'salt',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo)

        rho0=ncfile.rho0; g=cst().g

        if (len(depths)>1) and (min(depths)>0):
            zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[int(depths[0])-1:int(depths[-1]),:,:]

        else:
            zr=np.zeros((len(depths),ny2-ny1,nx2-nx1))
            for iz in range(0,len(depths)): zr[iz,:,:]=depths[iz]
            if (len(depths)==1) and (np.ndim(zr)==3): zr=zr[0,:,:]


        if varname=='rho1':
            var = roms.rho1_eos(temp,salt,zr,g,rho0)
        else:
            var = roms.rho1_eos(temp,salt,zr,g,rho0)


        if varname=='rho':
            var = var-1000.

        elif varname=='rho1':
            var = var-1000.

        elif varname=='buoy':
            var = -g*(var/rho0-1)

        elif varname=='dxbuoy':
            var = roms.diffx(-g*(var/rho0),pm) #/(0.5*(f[:,1:]+f[:,:-1]))

        elif varname=='dybuoy':
            var = roms.diffy(-g*(var/rho0),pn) #/(0.5*(f[1:,:]+f[:-1,:]))

        elif varname=='dybuoysmooth':
            var = roms.diffysmooth(-g*(var/rho0),pn) #/f

        elif varname=='dxbuoy2':
            var = 0.5*(roms.u2rho(roms.diffx(-g*var/rho0-1,pm))**2+roms.v2rho(roms.diffy(-g*var/rho0-1,pn)**2))
            
        elif varname in ['tendency','tendD','tendE','tendF']:   
            

            if u==None:
                u = roms.u2rho(variables(ncfile,'u',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo) )
            #else:
                #print 'using provided u, no need to recompute it'   
                
            if v==None:
                v =roms.v2rho( variables(ncfile,'v',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo))
            #else:
                #print 'using provided v, no need to recompute it'
                
            bx = roms.u2rho(roms.diffx(-g*(var/rho0-1),pm))
            by = roms.v2rho(roms.diffy(-g*(var/rho0-1),pn))
            
            if u.shape[0]==v.shape[0]:
                vx = roms.u2rho(roms.diffx(v ,pm))
                vy = roms.v2rho(roms.diffy(v ,pn))         
                ux = roms.u2rho(roms.diffx(u ,pm))
                uy = roms.v2rho(roms.diffy(u ,pn))                      
            else:
                vx = roms.psi2rho(roms.diffx(v ,roms.rho2v(pm)))
                uy = roms.psi2rho(roms.diffy(u ,roms.rho2u(pn)))
                vy = np.zeros(pm.shape)*np.nan
                vy[1:-1,:] = roms.diffy(v ,roms.rho2v(pn))  
                ux = np.zeros(pm.shape)*np.nan
                ux[:,1:-1] = roms.diffx(u ,roms.rho2u(pm))
   

            if varname=='tendency':
                var = -1*(bx * ux * bx + by * uy * bx + bx * vx * by + by * vy * by)
            elif varname=='tendD':
                var = -0.5 * (ux+vy) * ( bx**2 + by**2 )
            elif varname=='tendE':
                var = -0.5 * (ux-vy) * ( bx**2 - by**2 ) 
            elif varname=='tendF':
                var = -1. * (vx+uy)  *  bx * by
            
            
            
        elif varname=='tendency_geo':
        
            zeta = variables(ncfile,'zeta',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo) 
            
            u = -g*roms.rho2u(roms.v2rho(roms.diffy(zeta,pn))/f);
            v = g*roms.rho2v(roms.u2rho(roms.diffx(zeta,pm))/f);
            
            bx = roms.u2rho(roms.diffx(-g*(var/rho0-1),pm))
            by = roms.v2rho(roms.diffy(-g*(var/rho0-1),pn))
            vx = roms.u2rho(roms.diffx(roms.v2rho(v) ,pm))
            vy = roms.v2rho(roms.diffy(roms.v2rho(v) ,pn))         
            ux = roms.u2rho(roms.diffx(roms.u2rho(u) ,pm))
            uy = roms.v2rho(roms.diffy(roms.u2rho(u) ,pn))        

            var = -1*(bx * ux * bx + by * uy * bx + bx * vx * by + by * vy * by)
            
        elif varname=='tendency_ageo':
        
            zeta = variables(ncfile,'zeta',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo) 

            ugeo = -g*roms.rho2u(roms.v2rho(roms.diffy(zeta,pn))/f);
            vgeo = g*roms.rho2v(roms.u2rho(roms.diffx(zeta,pm))/f);
            
            u = variables(ncfile,'u',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo)-ugeo
            v = variables(ncfile,'v',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo)-vgeo       
            
            bx = roms.u2rho(roms.diffx(-g*(var/rho0-1),pm))
            by = roms.v2rho(roms.diffy(-g*(var/rho0-1),pn))
            vx = roms.u2rho(roms.diffx(roms.v2rho(v) ,pm))
            vy = roms.v2rho(roms.diffy(roms.v2rho(v) ,pn))         
            ux = roms.u2rho(roms.diffx(roms.u2rho(u) ,pm))
            uy = roms.v2rho(roms.diffy(roms.u2rho(u) ,pn))        

            var = -1*(bx * ux * bx + by * uy * bx + bx * vx * by + by * vy * by)
            
        elif varname=='tendency_3d':    
        
            #print 'beware temporary change for test!!!!'
            #var = temp
            #bx = roms.u2rho(roms.diffx(var,pm))
            #by = roms.v2rho(roms.diffy(var,pn))
            #bz = variables(ncfile,'dzt',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo)
            #if len(depths)==1: bz = bz[0,:,:]
            
            bx = roms.u2rho(roms.diffx(-g*(var/rho0-1),pm))
            by = roms.v2rho(roms.diffy(-g*(var/rho0-1),pn))
            bz = variables(ncfile,'dzbuoy',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo)
            if len(depths)==1: bz = bz[0,:,:]

            
            #compared to tendency, we add the contribution from vertical terms:
            u = variables(ncfile,'u',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo) 
            v = variables(ncfile,'v',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo)
            w = variables(ncfile,'w',ny1,ny2,nx1,nx2,depths,infiletime,topo,pm,pn,verbo=verbo)

            vx = roms.u2rho(roms.diffx(roms.v2rho(v) ,pm))
            vy = roms.v2rho(roms.diffy(roms.v2rho(v) ,pn))         
            ux = roms.u2rho(roms.diffx(roms.u2rho(u) ,pm))
            uy = roms.v2rho(roms.diffy(roms.u2rho(u) ,pn)) 
            
            wx = roms.u2rho(roms.diffx(w ,pm))
            wy = roms.v2rho(roms.diffy(w ,pn))            

            var = -1*(bx * ux * bx + by * uy * bx + bx * vx * by + by * vy * by + bz * wx * bx + bz * wy * by)
            #print 'beware temporary change for test!!!!'
            #var = -1*(bz * wx * bx + bz * wy * by)


        elif varname in ['strg']:  

            print('only surface comnputation yet')
            zeta = variables(ncfile,'zeta',ny1,ny2,nx1,nx2,[0],infiletime,topo,verbo=verbo)

            u = -1*roms.rho2u(roms.v2rho(roms.diffy(zeta,pn))); v = roms.rho2v(roms.u2rho(roms.diffx(zeta,pm)))
            var = g*roms.strains(u,v,pm,pn)/f

    ################################################################################### 
    elif varname in ['rhoF','rho1F','buoyF','dxbuoyF','dybuoyF']:	

        temp = variables(ncfile,'temp',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo) 
        salt = variables(ncfile,'salt',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo)

        rho0=ncfile.rho0; g=cst().g

        if (len(depths)>1) and (min(depths)>0):
            zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[int(depths[0])-1:int(depths[-1]),:,:]
        else:
            zr=np.zeros((len(depths),ny2-ny1,nx2-nx1))
            for iz in range(0,len(depths)): zr[iz,:,:]=depths[iz]
            #if (len(depths)==1) and (np.ndim(zr)==3): zr=zr[0,:,:]
            
        if varname=='rho1F':
            var = roms.rho1_eosF(temp,salt,zr,g,rho0) # + rho0
        else:
            var = roms.rho_eosF(temp,salt,zr,g,rho0) # + rho0
        

        if varname in ['rhoF','rho1F']:
            var = var - 1000. + rho0
      

        elif varname=='buoyF':
            var = -g*(var/rho0)

        elif varname=='dxbuoyF':
            var = roms.diffx(-g*var/rho0,pm)/(0.5*(f[:,1:]+f[:,:-1]))

        elif varname=='dybuoyF':
            var = -1*roms.diffy(-g*var/rho0,pn)/(0.5*(f[1:,:]+f[:-1,:]))



    ################################################################################### 
    elif varname in ['bvfF','pvF']:	


        rho0=ncfile.rho0; g=cst().g

        ####################################
        #Define levels where var will be computed
        ####################################

        if (len(depths)==1) and (max(depths)<=0):

            dz=1
            if (depths[0]==0): depths=[-1]; depths_w=[depths[0]-dz,depths[0]+dz];
            elif (depths[0]<0): depths_w=[depths[0]-dz,depths[0]+dz];

            zr=np.zeros((len(depths_w),ny2-ny1,nx2-nx1))
            zw=np.zeros((len(depths_w)+1,ny2-ny1,nx2-nx1))

            for iz in range(0,len(depths_w)): zr[iz,:,:]=depths_w[iz]; 
            for iz in range(1,len(depths)+1): zw[iz,:,:]=depths[iz-1]; 
            zw[0,:,:]=depths[0]-2*(depths[0]-depths_w[0])
            zw[-1,:,:]=depths[-1]-2*(depths[-1]-depths_w[-1])

        elif (len(depths)==1) and (depths[0]>0):

            depths_w=[depths[0]-1,depths[0]+1];
            zw=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[depths,:,:]
            zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[depths_w,:,:]       

        elif min(depths)>0:

            depths_w=depths
            zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[int(depths[0])-1:int(depths[-1]),:,:]
            zw=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[int(depths[0])-1:int(depths[-1])+1,:,:]

        else:

            depths_w=np.zeros((len(depths)+1))
            depths_w[0]=depths[0]+0.5*(depths[0]-depths[1])
            depths_w[1:-1]=0.5*(depths[1:]+depths[:-1])
            depths_w[-1]=depths[-1]+0.5*(depths[-1]-depths[-2])
            if verbo: print(depths_w)
            zr=np.zeros((len(depths_w),ny2-ny1,nx2-nx1))
            zw=np.zeros((len(depths_w)+1,ny2-ny1,nx2-nx1))
            for iz in range(0,len(depths_w)): zr[iz,:,:]=depths_w[iz]; 
            for iz in range(1,len(depths)+1): zw[iz,:,:]=depths[iz-1]; 
            zw[0,:,:]=depths[0]-2*(depths[0]-depths_w[0])
            zw[-1,:,:]=depths[-1]-2*(depths[-1]-depths_w[-1])


        if verbo: print(varname + ' will be computed at depths')
        #if verbo: print zw[1:-1,0,0]
        if verbo: print(depths)

        ####################################
        
        if varname=='bvfF':

            temp_w = variables(ncfile,'temp',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo) 
            salt_w = variables(ncfile,'salt',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo)

            [rho,bvf] = roms.rho_eosF(temp_w,salt_w,zr,g,rho0,zw)
            var = bvf[1:-1,:]

        elif varname=='pvF':

            temp_w = variables(ncfile,'temp',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo) 
            salt_w = variables(ncfile,'salt',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo)
            u_w = variables(ncfile,'u',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo)
            v_w = variables(ncfile,'v',ny1,ny2,nx1,nx2,depths_w,infiletime,topo,verbo=verbo)

            [var,pv1,pv2,pv3] = roms.PVF(temp_w,salt_w,u_w,v_w,zr,zw,f,g,rho0,pm,pn)
            

    ###################################################################################
    #various computed variables


    elif varname in ['dxT']:

        var = variables(ncfile,'temp',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo)
        var = roms.diffx(var,pm)

    elif varname in ['dyT']:

        var = variables(ncfile,'temp',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo)
        var = roms.diffy(var,pn)

    elif varname in ['uper']:	

        u = variables(ncfile,'u',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo) 

        var = u - roms.nanmean(u)
        
        if len(u.shape)==3:
            if verbo: print('la coin de lui')
            for iz in range (0,u.shape[0]):
                var[iz,:,:] = u[iz,:,:] - np.mean(u[iz,:,:]) #.mean(axis=2)

        else:
            for ix in range (0,u.shape[1]):
                var[:,ix] = u[:,ix] - np.mean(u[:,ix]) #u.mean(axis=1)
        
 
    elif varname in ['vper']:	

        v = variables(ncfile,'v',ny1,ny2,nx1,nx2,depths,infiletime,topo,verbo=verbo) 
        var = v*np.nan

        #if len(v.shape)==3:
            #for iy in range (0,v.shape[1]):
                #var[:,iy,:] = v[:,iy,:] - v.mean(axis=1)
        #else:
            #for iy in range (0,v.shape[0]):
                #var[iy,:] = v[iy,:] - v.mean(axis=0)

        var = v- roms.nanmean(v)


    ###################################################################################
    #Bottom Stress
    ###################################################################################

    elif varname in ['botsu','botsv','botsuv','uvbots']:	

        ubot = roms.u2rho(variables(ncfile,'u',ny1,ny2,nx1,nx2,[1],infiletime,topo=topo,verbo=verbo))
        vbot = roms.v2rho(variables(ncfile,'v',ny1,ny2,nx1,nx2,[1],infiletime,topo=topo,verbo=verbo)) 
        rho0=ncfile.rho0
        
        zw=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)

        Hz = zw[1,:,:] - zw[0,:,:]
        rdrg = ncfile.rdrg; vonKar = 0.41; Zob = 0.01; 

        cff = np.sqrt(ubot**2 + vbot**2)

        if varname=='botsu':	
            var = [rdrg + cff * (vonKar/np.log(Hz/Zob))**2] * ubot #* rho0

        elif varname=='botsv':	
            var = [rdrg + cff * (vonKar/np.log(Hz/Zob))**2] * vbot #* rho0

        elif varname=='botsuv':	
            var = [rdrg + cff * (vonKar/np.log(Hz/Zob))**2] * cff #* rho0

        elif varname=='uvbots':	

            var = [rdrg + cff * (vonKar/np.log(Hz/Zob))**2] * ubot #* rho0
            var1 = [rdrg + cff * (vonKar/np.log(Hz/Zob))**2] * vbot #* rho0

    ###################################################################################


    elif varname in ['detabuoy','dxibuoy']:	

        rho0=ncfile.rho0; g=cst().g

        if depths==[1]:
            temp = variables(ncfile,'temp',ny1,ny2,nx1,nx2,[1,2],infiletime,topo,verbo=verbo) 
            salt = variables(ncfile,'salt',ny1,ny2,nx1,nx2,[1,2],infiletime,topo,verbo=verbo)

            z_r=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[0:2,:,:]

            var = -g*(roms.rho_eos(temp,salt,z_r,g,rho0)/rho0-1)

            ######################

            if varname=='detabuoy':
                var = roms.diffeta(var,pm,z_r)

            elif varname=='dxibuoy':
                var = roms.diffxi(var,pn,z_r)


        else:
            print('not coded yet')
            sys.exit()


    ###################################################################################
    #Surface PV flux
    ###################################################################################



    elif varname in ['J1old','J2old']:	

        rho0=ncfile.rho0; g=cst().g

        temp = variables(ncfile,'temp',ny1,ny2,nx1,nx2,[0],infiletime,topo,verbo=verbo) 
        salt = variables(ncfile,'salt',ny1,ny2,nx1,nx2,[0],infiletime,topo,verbo=verbo)

        Hz=np.abs(roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[-2,:,:])

        try:
            hbls = variables(ncfile,'hbls',ny1,ny2,nx1,nx2,[0],infiletime,topo=topo,verbo=verbo,pm=pm,pn=pn,f=f,cubic=1)
            print('using hbls')
        except:
            hbls=Hz

        if varname=='J2old':

            zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[-1,:,:]

            if uwinds==None: [uwinds,vwinds] = roms.get_winds(ncfile,infiletime,ncname.wind,ny1,ny2,nx1,nx2)

            var0 = roms.rho_eos(temp,salt,zr,g,rho0)

            var = (roms.v2rho(roms.diffy(-g*var0/rho0,pn)) * roms.u2rho(uwinds)\
                 - roms.u2rho(roms.diffx(-g*var0/rho0,pm)) * roms.v2rho(vwinds))/rho0/hbls


        elif varname=='J1old':
            
            absvrt = variables(ncfile,'absvrt',ny1,ny2,nx1,nx2,[0],infiletime,topo=topo,verbo=verbo,pm=pm,pn=pn,f=f,cubic=1) 

            [stflux,ssflux] = roms.get_buoy_flux(ncfile,infiletime,ncname.frc,ny1,ny2,nx1,nx2,'m',Hz,temp,salt)

            [alpha,beta] = roms.alphabeta(temp,salt,ncfile.rho0)

            var =  g * absvrt * f * (alpha*stflux - beta*ssflux)/hbls

            del absvrt,stflux,ssflux,alpha,beta,Hz

    ###################################################################################
    #Surface PV flux
    ###################################################################################



    elif varname in ['J2','J1']:	

        rho0=ncfile.rho0; g=cst().g

        ###################

        if varname=='J2':

            temp = variables(ncfile,'temp',ny1,ny2,nx1,nx2,[0],infiletime,topo,verbo=verbo) 
            salt = variables(ncfile,'salt',ny1,ny2,nx1,nx2,[0],infiletime,topo,verbo=verbo)
            zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[-2:,:,:]
            zw=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[-2:,:,:]
            var0 = roms.rho_eos(temp,salt,zr[-1,:,:],g,rho0)
            del temp,salt


            if uwinds==None: [uwinds,vwinds] = roms.get_winds(ncfile,infiletime,ncname.wind,ny1,ny2,nx1,nx2)

            Fx = roms.u2rho(uwinds)/rho0
            Fy = roms.v2rho(vwinds)/rho0
            del uwinds,vwinds

            ###################
            Hz=np.abs(roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[-2,:,:])
            N_w = len(ncfile.dimensions['s_w'])

            try:
                AKv = variables(ncfile,'AKv',ny1,ny2,nx1,nx2,[N_w-1],infiletime,topo,verbo=verbo) 
            except:
                AKv = variables(ncfile,'AKt',ny1,ny2,nx1,nx2,[N_w-1],infiletime,topo,verbo=verbo) 

            print(AKv)
            
            #AKv[AKv>0.5] = AKv[AKv>0.5]-0.5

            if np.min(AKv)==np.max(AKv): print('no AKv present')
            u = variables(ncfile,'u',ny1,ny2,nx1,nx2,[N_w-2,N_w-1],infiletime,topo,verbo=verbo) 
            v = variables(ncfile,'v',ny1,ny2,nx1,nx2,[N_w-2,N_w-1],infiletime,topo,verbo=verbo) 

            #Fx = (Fx - AKv * roms.u2rho((u[1,:,:] - u[0,:,:]))/(zr[1,:,:]-zr[0,:,:])  ) #/  variables(ncfile,'hbls',ny1,ny2,nx1,nx2,[0],infiletime,topo=topo,verbo=verbo,pm=pm,pn=pn,f=f,cubic=1)     
            #Fy = (Fy - AKv * roms.v2rho((v[1,:,:] - v[0,:,:]))/(zr[1,:,:]-zr[0,:,:])  )# /  variables(ncfile,'hbls',ny1,ny2,nx1,nx2,[0],infiletime,topo=topo,verbo=verbo,pm=pm,pn=pn,f=f,cubic=1)   
            
            Fx = (Fx - AKv * roms.u2rho((u[1,:,:] - u[0,:,:]))/(zr[1,:,:]-zr[0,:,:])  ) /  (zw[1,:,:]-zw[0,:,:])
            Fy = (Fy - AKv * roms.v2rho((v[1,:,:] - v[0,:,:]))/(zr[1,:,:]-zr[0,:,:])  ) /  (zw[1,:,:]-zw[0,:,:])
            
            del u,v,AKv,Hz,zr

            var = roms.v2rho(roms.diffy(-g*var0/rho0,pn)) * Fx\
                 - roms.u2rho(roms.diffx(-g*var0/rho0,pm)) * Fy



        elif varname=='J1':

            N_w = len(ncfile.dimensions['s_w']); zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[-2:,:,:]

            temp = variables(ncfile,'temp',ny1,ny2,nx1,nx2,[N_w-2,N_w-1],infiletime,topo,verbo=verbo) 
            salt = variables(ncfile,'salt',ny1,ny2,nx1,nx2,[N_w-2,N_w-1],infiletime,topo,verbo=verbo) 

            zw = roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)

            Hz = zw[-1,:,:] - zw[-2,:,:]

            [stflux,ssflux] = roms.get_buoy_flux(ncfile,infiletime,ncname.frc,ny1,ny2,nx1,nx2,'m',Hz,temp[1,:,:],salt[1,:,:])

            [alpha,beta] = roms.alphabeta(temp[1,:,:],salt[1,:,:],ncfile.rho0)

            AKt = variables(ncfile,'AKt',ny1,ny2,nx1,nx2,[N_w-1],infiletime,topo,verbo=verbo) 

            dbdt = alpha*(stflux- AKt*[(temp[1,:,:] - temp[0,:,:])/(zr[1,:,:]-zr[0,:,:])]) \
                   - beta*(ssflux- AKt*[(salt[1,:,:] - salt[0,:,:])/(zr[1,:,:]-zr[0,:,:])])


            if np.min(AKt)==np.max(AKt): print('no AKt present')

            #del stflux,ssflux,AKt

            absvrt = variables(ncfile,'absvrt',ny1,ny2,nx1,nx2,[0],infiletime,topo=topo,verbo=verbo,pm=pm,pn=pn,f=f,cubic=1) 

            var =  g * absvrt * f * (dbdt)/ variables(ncfile,'hbls',ny1,ny2,nx1,nx2,[0],infiletime,topo=topo,verbo=verbo,pm=pm,pn=pn,f=f,cubic=1)

            del absvrt,Hz


    ###################################################################################
    #Bottom PV
    ###################################################################################



    elif varname in ['Jbotold']:

        rho0=ncfile.rho0; g=cst().g
        
        [botsu,botsv] = variables(ncfile,'uvbots',ny1,ny2,nx1,nx2,[1],infiletime,topo,verbo=verbo) 

        temp = variables(ncfile,'temp',ny1,ny2,nx1,nx2,[1,2],infiletime,topo,verbo=verbo) 
        salt = variables(ncfile,'salt',ny1,ny2,nx1,nx2,[1,2],infiletime,topo,verbo=verbo)

        zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[0:2,:,:]
        Hz=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[:2,:,:]

        try:
            hbbls = variables(ncfile,'hbbls',ny1,ny2,nx1,nx2,[0],infiletime,topo=topo,verbo=verbo,pm=pm,pn=pn,f=f,cubic=1)
            print('using hbbls')
        except:
            zw=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[:2,:,:]
            hbbls=zw[-1,:,:] - zw[-2,:,:]

        buoy = -g*(roms.rho_eosF(temp,salt,zr,g,rho0)/rho0-1)
        del temp,salt

        Fx = (botsu[0,:,:]) /hbbls
        Fy = (botsv[0,:,:]) /hbbls

        var = roms.u2rho(roms.diffxi(buoy,pm,zr)) * Fy  - roms.v2rho(roms.diffeta(buoy,pn,zr)) * Fx




    ###################################################################################
    #Bottom PV
    ###################################################################################



    elif varname in ['Jbot']:

        rho0=ncfile.rho0; g=cst().g
        
        [botsu,botsv] = variables(ncfile,'uvbots',ny1,ny2,nx1,nx2,[1],infiletime,topo,verbo=verbo) 

        temp = variables(ncfile,'temp',ny1,ny2,nx1,nx2,[1,2],infiletime,topo,verbo=verbo) 
        salt = variables(ncfile,'salt',ny1,ny2,nx1,nx2,[1,2],infiletime,topo,verbo=verbo)

        zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[0:2,:,:]
        Hz=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[:2,:,:]

        buoy = -g*(roms.rho_eosF(temp,salt,zr,g,rho0)/rho0-1)
        del temp,salt

        u = variables(ncfile,'u',ny1,ny2,nx1,nx2,[1,2],infiletime,topo,verbo=verbo) 
        v = variables(ncfile,'v',ny1,ny2,nx1,nx2,[1,2],infiletime,topo,verbo=verbo)

        try:
            AKv = variables(ncfile,'AKv',ny1,ny2,nx1,nx2,[2],infiletime,topo,verbo=verbo) 
        except:
            AKv = variables(ncfile,'AKt',ny1,ny2,nx1,nx2,[2],infiletime,topo,verbo=verbo) 


        Fx = (AKv * roms.u2rho((u[1,:,:] - u[0,:,:]))/(zr[1,:,:]-zr[0,:,:]) + botsu[0,:,:]) /np.abs(Hz[1,:,:]-Hz[0,:,:])  
        Fy = (AKv * roms.v2rho((v[1,:,:] - v[0,:,:]))/(zr[1,:,:]-zr[0,:,:]) + botsv[0,:,:]) /np.abs(Hz[1,:,:]-Hz[0,:,:])
        

        var = roms.u2rho(roms.diffxi(buoy,pm,zr)) * Fy  - roms.v2rho(roms.diffeta(buoy,pn,zr)) * Fx
        
        #hbbls = variables(ncfile,'hbbls',ny1,ny2,nx1,nx2,[0],infiletime,topo=topo,verbo=verbo,pm=pm,pn=pn,f=f,cubic=1)
        #var[hbbls>=topo-1] = np.nan

    ###################################################################################
    #Bottom PV
    ###################################################################################



    elif varname in ['Jbottest']:

        rho0=ncfile.rho0; g=cst().g
        
        [botsu,botsv] = variables(ncfile,'uvbots',ny1,ny2,nx1,nx2,[1],infiletime,topo,verbo=verbo) 

        temp = variables(ncfile,'temp',ny1,ny2,nx1,nx2,[1,2],infiletime,topo,verbo=verbo) 
        salt = variables(ncfile,'salt',ny1,ny2,nx1,nx2,[1,2],infiletime,topo,verbo=verbo)

        zr=roms.get_depths_zoom(ncfile,topo,infiletime,'r',ny1,ny2,nx1,nx2)[0:2,:,:]
        Hz=roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[:2,:,:]

        buoy = -g*(roms.rho_eosF(temp,salt,zr,g,rho0)/rho0-1)
        del temp,salt

        u = variables(ncfile,'u',ny1,ny2,nx1,nx2,[1,2],infiletime,topo,verbo=verbo) 
        v = variables(ncfile,'v',ny1,ny2,nx1,nx2,[1,2],infiletime,topo,verbo=verbo)

        try:
            AKv = variables(ncfile,'AKv',ny1,ny2,nx1,nx2,[2],infiletime,topo,verbo=verbo) 
        except:
            AKv = variables(ncfile,'AKt',ny1,ny2,nx1,nx2,[2],infiletime,topo,verbo=verbo) 

        Fx = (0*AKv * roms.u2rho((u[1,:,:] - u[0,:,:]))/(zr[1,:,:]-zr[0,:,:]) + botsu[0,:,:]) /np.abs(Hz[1,:,:]-Hz[0,:,:])  
        Fy = (0*AKv * roms.v2rho((v[1,:,:] - v[0,:,:]))/(zr[1,:,:]-zr[0,:,:]) + botsv[0,:,:]) /np.abs(Hz[1,:,:]-Hz[0,:,:])

        
        #hbbls = variables(ncfile,'hbbls',ny1,ny2,nx1,nx2,[0],infiletime,topo=topo,verbo=verbo,pm=pm,pn=pn,f=f,cubic=1)
        #var = hbbls-np.abs(Hz[1,:,:]-Hz[0,:,:])
        #Fx = (botsu[0,:,:]) /hbbls
        #Fy = (botsv[0,:,:]) /hbbls

        var = roms.u2rho(roms.diffxi(buoy,pm,zr)) * Fy  - roms.v2rho(roms.diffeta(buoy,pn,zr)) * Fx

    ###################################################################################
    ###################################################################################




    elif varname in ['buoyflux']:	

        rho0=ncfile.rho0; g=cst().g

        temp = variables(ncfile,'temp',ny1,ny2,nx1,nx2,[0],infiletime,topo,verbo=verbo) 
        salt = variables(ncfile,'salt',ny1,ny2,nx1,nx2,[0],infiletime,topo,verbo=verbo)

        Hz=np.abs(roms.get_depths_zoom(ncfile,topo,infiletime,'w',ny1,ny2,nx1,nx2)[-2,:,:])

        if varname=='buoyflux':

            [stflux,ssflux] = roms.get_buoy_flux(ncfile,infiletime,ncname.frc,ny1,ny2,nx1,nx2,'m',Hz,temp,salt)

            [alpha,beta] = roms.alphabeta(temp,salt,ncfile.rho0)

            var =  alpha*stflux - beta*ssflux








    if opened==True: ncfile.close()

    ###################################################################################
    #Bottom Stress
    ###################################################################################


    if varname in ['uvbots','dbuoybots']:
        return [var,var1]

    else:
        return var







###################################################################################
#Load custom colormap defined for each variable
###################################################################################

def colormap(vavar):

    #Standard BLUE-WHITE-RED scale - useful for pos-neg field (w,pv..etc...
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
    cdict7 = {'blue':  ((0.0, 0.0, 0.0),
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

    if vavar in ['temp','salt','dsalt']:
        my_cmap = py.cm.jet
    elif vavar in ['buoyflux','div','dvdz','pv','pvF','pv1','pv2','pv3','J1',\
                    'J2','Jbot','J1old','J2old','omega','absvrt','vrt','vrtdiv','pvold',\
                    'dxpv','wb','u','uper','v','tendency','tendency_3d','tendency_geo','tendency_ageo','dxbuoy2','dxbuoy','str','wb',\
                    'Tadv','THdiff','TVmix','TForc','Sadv','SHmix','SVmix','SForc','Trate','Srate',\
                    'Tadv_pert','THdiff_pert','TVmix_pert','TForc_pert','Sadv_pert','SHmix_pert','SVmix_pert','SForc_pert','Trate_pert','Srate_pert',\
                    'TXadv','TYadv','TVadv',\
                    'bvf','bvfF','dybuoy','dzt']:
        my_cmap = col.LinearSegmentedColormap('my_colormap',cdict5,256)
    elif vavar in ['vrt','ageo','ow_alter','ow_alter1','ow_alter2']:
        my_cmap = col.LinearSegmentedColormap('my_colormap',cdict6,256)
    elif vavar in ['test']:
        my_cmap = col.LinearSegmentedColormap('my_colormap',cdict7,256)
    elif vavar in ['temp']:
        my_cmap = col.LinearSegmentedColormap('my_colormap',cdict9,256)
    elif vavar in ['w','ow']:
        my_cmap = col.LinearSegmentedColormap('my_colormap',cdict10,256)
    else:
        my_cmap = py.cm.spectral 

    return my_cmap


###################################################################################
#define levels
###################################################################################

def levels(vavar, var=None , nlev=None):



    
    if vavar=='topo':
        #levelsvar = [21, 50, 100, 500, 1000, 2000, 3000, 4000, 5000]
        levelsvar = [21, 200, 600,1000, 2000, 3000, 4000, 5000]

    elif vavar=='topo2':
        levelsvar = [21]

    #elif vavar=='w':
    #    levelsvar = np.arange(-1,1.01,0.05)*2e-3


    elif nlev>10:

        minvar=np.nanmin(var)
        maxvar=np.nanmax(var)

        if vavar in ['pv','pv1','pv2','pv3','pvold','dxpv0']:
            levelsvar = np.arange(-1,1.01,0.01)*2e-8

        elif vavar in ['vrt','div','vrtdiv','ageo']:
            levelsvar=np.arange(-2,2.05,0.05)*0.5

        elif vavar in ['w','omega','ow1']:
            maxabs = 0.5*max(abs(minvar),abs(maxvar))
            levelsvar=np.arange(-1*maxabs,maxabs+maxabs/nlev,maxabs/nlev)

        elif vavar in ['w0']:
            levelsvar = np.arange(-0.025,0.025,0.001)

        elif vavar in ['wb0']:
            levelsvar = np.arange(-0.00005,0.00005,0.000001)

        elif vavar in ['salt0']:
            levelsvar = np.arange(34.5,36.8,0.01)
            #levelsvar = np.arange(36.1,36.4,0.01)
 
        elif vavar in ['rhoF','rho1F','rho']:
            levelsvar = np.arange(25.,28.,0.01)
            
        elif vavar in ['temp2']:
            levelsvar = np.arange(12,24,0.1)

        elif vavar in ['ow','ow_alter','ow_alter1','ow_alter2']:
            levelsvar = np.arange(-1,1.01,0.01)*0.5e-8

        elif vavar in ['pv','div','vrt','w','dybuoy']:
            maxabs = 0.00025 #0.5*max(abs(minvar),abs(maxvar))
            levelsvar=np.arange(-1*maxabs,maxabs+maxabs/nlev,maxabs/nlev)

        elif vavar in ['J1','J2','Jbot','J1old','J2old','bvf0','bvfF']:
            maxabs = 0.1*max(abs(minvar),abs(maxvar))
            levelsvar=np.arange(-1*maxabs,maxabs+maxabs/nlev,maxabs/nlev)

        elif vavar in ['u','v','uper','omega','pv','pvF','pvold','dxpv','wb','buoyflux','dxbuoy']:
            maxabs = 0.5*max(abs(minvar),abs(maxvar))
            levelsvar=np.arange(-1*maxabs,maxabs+maxabs/nlev,maxabs/nlev)

        elif vavar in ['sym']:
            maxabs = max(abs(minvar),abs(maxvar))
            levelsvar=np.arange(-1*maxabs,maxabs+maxabs/nlev,maxabs/nlev)
            
        elif vavar in ['temp','rho']:
            #minvar=minvar+11 #for jim.figure1
            #maxvar=maxvar-1 #for jim.figure2
            #nlev=100
            #minvar=minvar+0.5 #for jim.figure2
            #maxvar=maxvar+0.5 #for jim.figure2
            #minvar=0.; maxvar=25.; #for GULZ.figure1

            #minvar=minvar+5 #for latmix.nefro

            levelsvar=np.arange(minvar,maxvar+(maxvar-minvar)/nlev,(maxvar-minvar)/nlev)

        elif vavar in ['dsalt']:
            levelsvar = np.arange(-0.5,0.5,0.01)

        elif vavar in ['absvrt']:
            maxabs = 5.
            levelsvar=np.arange(-1*maxabs,1*maxabs+maxabs/nlev,maxabs/nlev)

        elif vavar in ['dxbuoy2']:           
            levelsvar=np.arange(minvar,maxvar+(maxvar-minvar)/nlev,(maxvar-minvar)/nlev)
            
        elif vavar in ['tendency','tendency_3d','tendency_geo','tendency_ageo']:           
            levelsvar=np.arange(-1.,1.,0.01)*1e-15
             
        elif vavar in ['strg']:           
            levelsvar=np.arange(0,1.,0.01)*0.01
            
        elif vavar in ['hbls0']:           
            levelsvar=np.arange(50,100.,1.)
        elif vavar in ['dzt']:
            levelsvar=np.arange(0,0.03,0.0005)

        else:
            levelsvar=np.arange(minvar,maxvar+(maxvar-minvar)/nlev,(maxvar-minvar)/nlev)


    else:

        minvar=np.nanmin(var)
        maxvar=np.nanmax(var)

        if vavar in ['pv0','pvF','div','vrt0','dxpv','w']:
            maxabs = 0.05*max(abs(minvar),abs(maxvar))
            levelsvar=[-1*maxabs,maxabs]

        elif vavar in ['vrt']:
            levelsvar=[-2.,2.]

        elif vavar in ['wb']:
            maxabs = 0.5*max(abs(minvar),abs(maxvar))
            levelsvar=[-1*maxabs,-0.5*maxabs,0.5*maxabs,maxabs]

        elif vavar in ['str']:
            maxabs = max(abs(minvar),abs(maxvar))
            levelsvar=[-0.8*maxabs,-0.5*maxabs,0.5*maxabs,0.8*maxabs]

        elif vavar in ['dsalt']:
            levelsvar=[-10, -0.1]

        elif vavar in ['absvrt0']:
            levelsvar=[-10, -1]

        elif vavar in ['vrt']:
            levelsvar=[-2, 2]

        elif vavar in ['ageo']:
            levelsvar=[-1e4, 0, 1e4]

        elif vavar in ['J1','J2','Jbot','J1old','J2old']:
            maxabs = 0.1*max(abs(minvar),abs(maxvar))
            levelsvar=[-1*maxabs, maxabs]

        elif vavar in ['z','zw']:
            levelsvar=np.arange(0,100,1)
            #levelsvar=np.arange(minvar,maxvar,(maxvar-minvar)/10)

        elif vavar in ['temp0']:
            levelsvar=[17.5, 18.5, 19.5, 20.5]
            
        elif vavar in ['rho0','rhoF0','rho1F0']:
            levelsvar=[26.4, 26.6]
            
        elif vavar in ['botsuv']:
            levelsvar=np.arange(0.005,0.01,0.001)/10

        elif vavar in ['buoyF0']:
            levelsvar=np.arange(-1,1,0.05)/100

        elif vavar in ['rhoF0']:
            levelsvar=np.arange(-10,100,0.1)
                    
        elif vavar in ['dxbuoy2']:
           levelsvar=np.arange(0,maxvar+(maxvar)/10,(maxvar)/10)
            
        else:


            levelsvar=np.arange(minvar,maxvar+(maxvar-minvar)/10,(maxvar-minvar)/10)

    return levelsvar



###################################################################################
#define colorabar levels using levelsvar
###################################################################################

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
# define var and axis for plots (if horizontal or vertical)
###################################################################################

def getvarxy(var,x,y,depths,pm,pn,z_r=None,z_w=None,momoy=0):

    if type(var) is int:

        varout=var
        xvar=0
        yvar=0

    #horizontal plot
    elif (len(depths)==1):

        if np.ndim(var)==3: var = var[0,:,:]

        if (var.shape[0]>5) and (var.shape[1]>5):

            varout = var

            if x.shape==varout.shape: 
                xvar=x; yvar=y
            elif (x.shape[1]==varout.shape[1]) and (x.shape[0]==varout.shape[0]+1):
                xvar=roms.rho2v(x); yvar=roms.rho2v(y)
            elif (x.shape[0]==varout.shape[0]) and (x.shape[1]==varout.shape[1]+1):
                xvar=roms.rho2u(x); yvar=roms.rho2u(y)      
            elif (x.shape[0]==varout.shape[0]+1) and (x.shape[1]==varout.shape[1]+1):
                xvar=roms.rho2psi(x); yvar=roms.rho2psi(y)  
            '''
            #downgrade data if too many points for plotting
            if (var.shape[0]*var.shape[1]>500000):
                var = var[::var.shape[0]/500,::var.shape[1]/500]
                xvar=xvar[::var.shape[0]/500,::var.shape[1]/500]
                yvar=yvar[::var.shape[0]/500,::var.shape[1]/500]
            '''
        #xvar=(xvar-np.min(xvar))/pm[0,0]/1000
        #yvar=(yvar-np.min(yvar))/pn[0,0]/1000

    #linear (only 1D)
        else:

            #x slice
            if var.shape[0]>5:
                varout =  var[:,1]
                if y.shape[0]==varout.shape[0]: xvar = y[:,1]
                elif y.shape[0]==varout.shape[0]+1: xvar = roms.rho2v(y)[:,1]
                yvar = None

            #y slice
            elif var.shape[1]>5:
                varout =  var[1,:]
                if x.shape[1]==varout.shape[0]: xvar = x[1,:]
                elif x.shape[1]==varout.shape[0]+1: xvar = roms.rho2u(x)[1,:]
                yvar = None


    #vertical plot on sigma-levels
    elif (min(depths)>0) and ((var.shape[1]<=5) or (var.shape[2]<=5)):

        if var.shape[0]==len(depths): zvar=z_r
        elif var.shape[0]==len(depths)+1: zvar=z_w
        elif var.shape[0]==len(depths)-1: zvar=z_w[1:-1,:,:]
        else: print('HOUSTON WE HAVE A PROBLEM')

        #x slice
        if var.shape[1]>5:
            varout =  np.transpose(var[:,:,1])
            #zivar,yivar=np.mgrid[0:var.shape[0],0:var.shape[1]]
            yivar,zivar=np.meshgrid(y[:,1],list(range(var.shape[0])))
            if var.shape[1]==x.shape[0]-1: yivar=yivar+0.5; zivar=roms.rho2v(zivar)
            [xvar,yvar]=[np.transpose(yivar),np.transpose(zivar)]
            yvar = np.transpose(zvar[:,:,1])
                  
            if var.shape[1]==x.shape[0]-1: xvar,yvar = roms.rho2v(xvar),roms.rho2v(yvar)
            print('hum hum',var.shape,varout.shape,xvar.shape,yvar.shape)
            
        #y slice
        elif var.shape[2]>5:
            varout =  np.transpose(var[:,1,:]); 
            #zivar,xivar=np.mgrid[0:var.shape[0],0:var.shape[2]]
            xivar,zivar=np.meshgrid(x[1,:],list(range(var.shape[0])))
            if var.shape[2]==x.shape[1]-1: xivar=xivar+0.5; zivar=roms.rho2u(zivar)
            [xvar,yvar]=[np.transpose(xivar),np.transpose(zivar)]
            yvar = np.transpose(zvar[:,1,:])
            
            if var.shape[2]==x.shape[1]-1: xvar,yvar = roms.rho2u(xvar),roms.rho2u(yvar)
            print('hum',var.shape,varout.shape,xvar.shape,yvar.shape)


    #vertical plot on z-levels
    elif (len(depths)>1) and ((var.shape[1]<=5) or (var.shape[2]<=5)):
        #x slice
        if var.shape[1]>5:
            if var.shape[2]==x.shape[1]:
                varout =  np.transpose(var[:,:,1]); 
                #varout =  np.transpose(np.mean(var,2)); 
            else:
                varout =  np.transpose(roms.u2rho(var)[:,: ,1]); 
                #varout =  np.transpose(np.mean(roms.u2rho(var),2)); 
            #zivar,yivar=np.mgrid[0:var.shape[0],0:var.shape[1]]
            yivar,zivar=np.meshgrid(y[:,1],list(range(var.shape[0])))
            if var.shape[1]==x.shape[0]-1:
                yivar=roms.rho2u(yivar)
                zivar=roms.rho2u(zivar)
            [xvar,yvar]=[np.transpose(yivar),np.transpose(zivar)]

        #y slice
        elif var.shape[2]>5:
            if var.shape[1]==x.shape[0]:
                varout =  np.transpose(var[:,1,:]); 
            else:
                varout =  np.transpose(roms.v2rho(var)[:,1,:]); 
            #zivar,xivar=np.mgrid[0:var.shape[0],0:var.shape[2]]
            xivar,zivar=np.meshgrid(x[1,:],list(range(var.shape[0])))
            if var.shape[2]==x.shape[1]-1: 
                xivar=roms.rho2u(xivar)
                zivar=roms.rho2u(zivar)
            [xvar,yvar]=[np.transpose(xivar),np.transpose(zivar)]


        if (var.shape[0]==len(depths)): 
            for iz in range(0,var.shape[0]): yvar[:,iz] = depths[iz]
        elif (var.shape[0]==len(depths)-1): 
            for iz in range(0,var.shape[0]): yvar[:,iz] = 0.5*(depths[iz]+depths[iz+1])
        elif (var.shape[0]==len(depths)+1): 
            for iz in range(1,var.shape[0]-1): 
                yvar[:,iz] = 0.5*(depths[iz]+depths[iz+1])
            yvar[:,0] = depths[0]-0.5*(depths[0]-depths[1]); 
            yvar[:,-1] = depths[-1]+0.5*(depths[iz]-depths[iz+1]); 



    #averaged vertical plot on z-levels
    elif (len(depths)>1) and (momoy==1):
        if var.shape[2]==x.shape[1]:
            varout =  np.transpose(np.mean(var[:,1:-1,:],1)); 
        else:
            varout =  np.transpose(roms.nanmean(roms.u2rho(var)[:,1:-1,:],1)); 

        zivar,yivar=np.mgrid[0:var.shape[0],0:varout.shape[0]]

        if varout.shape[0]==x.shape[1]-1: yivar=yivar+0.5
        [xvar,yvar]=[np.transpose(yivar),np.transpose(zivar)]


        if (var.shape[0]==len(depths)): 
            for iz in range(0,var.shape[0]): yvar[:,iz] = depths[iz]
        elif (var.shape[0]==len(depths)-1): 
            for iz in range(0,var.shape[0]): yvar[:,iz] = 0.5*(depths[iz]+depths[iz+1])
        elif (var.shape[0]==len(depths)+1): 
            for iz in range(1,var.shape[0]-1): 
                yvar[:,iz] = 0.5*(depths[iz]+depths[iz+1])
            yvar[:,0] = depths[0]-0.5*(depths[0]-depths[1]); 
            yvar[:,-1] = depths[-1]+0.5*(depths[iz]-depths[iz+1]); 


    #3D plot
    else:
        varout=var
        xvar=None
        yvar=None





    return [varout,xvar,yvar]


###################################################################################
# define var and axis for u,v plots (if horizontal or vertical)
###################################################################################

def getvarxy_uv(u,v,x,y,depths,pm,pn,dens=0,momoy=0):



    if (len(depths)==1):

        if np.ndim(u)==3: u = u[0,:,:]; v = v[0,:,:]

        if u.shape!=v.shape:
            u = roms.u2rho(u); v = roms.v2rho(v)

        if dens==1:
            [uvar,vvar,xvar,yvar] = [u,v,x,y]
        elif dens==0:
            nvectx=max(u.shape[-1]/25,1)*2*2
            nvecty=max(u.shape[-2]/25,1)*2*2*2
            [uvar,vvar,xvar,yvar] = [u[2::nvecty,1::nvectx],v[2::nvecty,1::nvectx],x[2::nvecty,1::nvectx],y[2::nvecty,1::nvectx]]
        else:
            print('density is ',dens)
            nvectx=dens
            nvecty=dens
            [uvar,vvar,xvar,yvar] = [u[2::nvecty,1::nvectx],v[2::nvecty,1::nvectx],x[2::nvecty,1::nvectx],y[2::nvecty,1::nvectx]]

        coefuv=1

    else:

        if (u.shape[1]==v.shape[1]-1): u = roms.v2rho(u)
        elif (u.shape[2]==v.shape[2]-1): u = roms.u2rho(u)
        
        if dens==0:
            nvectx=max(len(depths)/8/2.,1)
            nvecty=max(max(u.shape[1],u.shape[2])/25*3/4.,1)
            coefuv=10 #20 #20 #0.5
            print('nvectx, nvecty',  nvectx, nvecty)
            
        else:
            nvectx=dens/3*2.
            nvecty=dens
            coefuv=20 #20 #0.5
            
        if u.shape==v.shape:
            [uvar,xvar,yvar]= getvarxy(u,x,y,depths,pm,pn,None,None,momoy)
            [vvar,xvar,yvar]= getvarxy(coefuv*v,x,y,depths,pm,pn,None,None,momoy)


        [uvar,vvar,xvar,yvar] = [uvar[2::nvecty,1::nvectx],vvar[2::nvecty,1::nvectx],xvar[2::nvecty,1::nvectx],yvar[2::nvecty,1::nvectx]]
        #remove zonal mean


    return  [uvar,vvar,xvar,yvar,coefuv]


###################################################################################
# define var and axis for 3D plots
###################################################################################

def getvarxyz(var,varname):

    var = np.swapaxes(var, 0, 2)
    zi,yi,xi=np.mgrid[0:var.shape[0],0:var.shape[1],0:var.shape[2]]
    zi=1.*zi/var.shape[0]; xi=1.*xi/var.shape[2]; yi=1.*yi/var.shape[1];

    #xi,yi,zi=np.mgrid[0:var.shape[0],0:var.shape[1],0:var.shape[2]]
    #zi=1.*zi/var.shape[2]; xi=1.*xi/var.shape[0]; yi=1.*yi/var.shape[1];


    if varname in ['pv','pvold']:
        varmin=-0.25e-8; varmax=-1*varmin
    else:
        varmin=np.nanmin(var)+(np.nanmax(var)-np.nanmin(var))/10
        varmax=np.nanmax(var)-(np.nanmax(var)-np.nanmin(var))/10

    return [var,zi,yi,xi,varmin,varmax]





###################################################################################
# get ocean date from file and convert it to '%m/%d - %H:00'
###################################################################################

def date(ncfile,infiletime,opened=False):

    if isinstance(ncfile, str): ncfile  = Dataset(ncfile, 'r', format='NETCDF3_CLASSIC'); opened=True

    oceantime = int(np.array(ncfile.variables['ocean_time'][infiletime]))%(360*24*3600)
    oceandate = tm.strftime("%B %d - %H:00", tm.gmtime(oceantime))

    if opened==True: ncfile.close()

    return oceandate

###################################################################################


def p(u,x=None,y=None,N=10):
    levels=np.arange(np.nanmin(u),np.nanmax(u)+(np.nanmax(u)-np.nanmin(u))/N,(np.nanmax(u)-np.nanmin(u))/N)
    if x==None: py.contourf(u,levels,extend='both'); py.colorbar(); py.show()
    else: py.contourf(x,y,u,levels,extend='both'); py.colorbar(); py.show()


def vec(x,y,u,v,nn,nny=None):
    if nny==None: nny=nn
    py.quiver(x[::nny,::nn],y[::nny,::nn],u[::nny,::nn],v[::nny,::nn]); py.axis('equal');py.show()


def p4(u,u2,u3,u4,x=None,y=None,N=10,samelev=False):
    levels=np.arange(np.nanmin(u),np.nanmax(u)+(np.nanmax(u)-np.nanmin(u))/N,(np.nanmax(u)-np.nanmin(u))/N)/2.

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














