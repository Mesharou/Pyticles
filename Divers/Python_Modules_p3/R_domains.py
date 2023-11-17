
# some ROMS tools written in python
import R_tools as tools

#import matplotlib.pyplot as py

import numpy as np

#copy data
from copy import copy

############################################################
'''
18/02/02: Add mask condition in find_larger_domain_possible
'''

############################################################
# R_DOMAINS
###########################################################"""


def find_larger_domain_possible(simul,domain,depth,max_domain):
    starting_point = domain[1],domain[0]
    inc=50.; ind = list(range(4))
    L = np.array([-inc/2, inc/2, -inc, inc],int); count=0
    newdomain = [starting_point[1]-1,starting_point[1]+1,starting_point[0]-1,starting_point[0]+1]
    #py.clf(); py.contourf(simul.topo.T,np.arange(0,-depth,2));py.contour(simul.topo.T,[-depth]); py.ion(); py.savefig('domains_test.png')
    while len(ind)>0 and count<100:
        for i in ind:
            oldval = newdomain[i]
            if i<2:
                newdomain[i]=int(np.nanmin([np.nanmax([newdomain[i]+L[i],max_domain[0]]),max_domain[1]]))
            else:
                newdomain[i]=int(np.nanmin([np.nanmax([newdomain[i]+L[i],max_domain[2]]),max_domain[3]]))
            if np.nanmin(simul.topo[newdomain[2]:newdomain[3],newdomain[0]:newdomain[1]])<-depth or np.isnan(simul.mask[newdomain[2]:newdomain[3],newdomain[0]:newdomain[1]].min()):
                newdomain[i]=oldval; L[i]=L[i]/2;
            if (L[i] in [0,-1]) or (newdomain[i]==max_domain[i]):
                del ind[ind.index(i)]
            #py.plot([newdomain[2], newdomain[3], newdomain[3],newdomain[2],newdomain[2]],\
            #[newdomain[0], newdomain[0], newdomain[1],newdomain[1], newdomain[0]] , 'g-', linewidth=1); py.draw()
            count+=1
    return newdomain



def my_domain(simul, domname,depth0):
    '''
    Define domain used to compute spectra
    
    '''
    # define center xy
    if 'lucky' in simul.simul:
        xy = tools.find_points(simul.x,simul.y,-33.,38.)
    elif 'chab' in simul.simul or 'sarg' in simul.simul or simul.simul in ['nwat','atlbig']:
        if 'sargasso' in domname:
            xy = tools.find_points(simul.x,simul.y,-76.,29.)
        elif 'bump' in domname:
            xy = tools.find_points(simul.x,simul.y,-78.23-0.5,31.56-0.8)
            #xy = tools.find_points(simul.x,simul.y,-75.7,28.499)
        else:
            xy = tools.find_points(simul.x,simul.y,-77,33.)
    elif 'nese' in simul.simul:
        xy = tools.find_points(simul.x,simul.y,-65,39.)
    elif 'cuc' in domname:
        #xy = tools.find_points(simul.x,simul.y,-125+360.,36.)
        if 'midcal' in simul.simul:
            #cuc
            xy = tools.find_points(simul.x,simul.y,237.5-360.,36.5)
        else:
            xy = tools.find_points(simul.x,simul.y,237.5,36.5)
    elif 'midcalL3' in domname:
        xy = tools.find_points(simul.x,simul.y,-120.64339,34.31438)
    elif 'midcalL2' in domname:
        xy = tools.find_points(simul.x,simul.y,-121.19,34.2)
    elif 'midcal' in domname:
        xy = tools.find_points(simul.x,simul.y,-120.723,34.23)
    elif 'nbimin' in domname:
        xy = tools.find_points(simul.x,simul.y,-79.34,25.8)
    elif 'bimin' in domname:
        xy = tools.find_points(simul.x,simul.y,-79.4,25.9)
    elif 'eurec4a' in domname:
        xy = tools.find_points(simul.x,simul.y,-54,10)
        
    ################
    res = np.mean(simul.pm)

    max_domain=[simul.coord[1]/12,simul.coord[1]-simul.coord[1]/12,simul.coord[3]/12,simul.coord[3]-simul.coord[3]/12] # Limit domain to exclude sponges
    if 'chab' in simul.simul or simul.simul in ['nwat','atlbig']:
        dm = np.int(130000*res)
        domain1= [xy[1][0]-dm,xy[1][0]+dm,xy[0][0]-dm,xy[0][0]+dm]
    elif 'nese' in simul.simul:
        dm = np.int(200000*res)
        domain1= [xy[1][0]-dm,xy[1][0]+dm,xy[0][0]-dm,xy[0][0]+dm]
    elif 'midcalL3' in domname:
        dm = np.int(40000/2*res)
        domain1= [xy[1][0]-dm,xy[1][0]+dm,xy[0][0]-dm,xy[0][0]+dm]
    elif 'midcalL2' in domname:
        dm = np.int(128000/2*res)
        domain1= [xy[1][0]-dm,xy[1][0]+dm,xy[0][0]-dm,xy[0][0]+dm]
    else:
        domain1 = find_larger_domain_possible(simul,[xy[1][0],xy[0][0]],depth0,max_domain);

    ################
    #make a square domain centered on xy

    diff = (domain1[3]-domain1[2])-(domain1[1]-domain1[0])
    domain_save = copy(domain1)
    if diff<0:
        dm = domain1[3]-domain1[2]
        domain1[0],domain1[1] = xy[1][0]-dm/2,xy[1][0]+dm/2+dm%2
        if domain1[0]<domain_save[0]: domain1[0],domain1[1] = domain_save[0],domain1[1]-domain1[0]+domain_save[0]
        if domain1[1]>domain_save[1]: domain1[0],domain1[1] = domain1[0]-(domain1[1]-domain_save[1]),domain_save[1]
    else:
        dm = domain1[1]-domain1[0]
        domain1[2],domain1[3] = xy[0][0]-dm/2,xy[0][0]+dm/2+dm%2
        if domain1[2]<domain_save[2]: domain1[2],domain1[3] = domain_save[2],domain1[3]-domain1[2]+domain_save[2]
        if domain1[3]>domain_save[3]: domain1[2],domain1[3] = domain1[2]-(domain1[3]-domain_save[3]),domain_save[3]

    ################

    return domain1

