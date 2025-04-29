###################################################################################
# R_SMOOTH
###################################################################################
"""





"""
from __future__ import print_function



#from builtins import zip
#from builtins import range
import numpy as np



#import R_tools as tools


#######################################################
# 1D- SMOOTHING (used in instability.py)
#######################################################



def smooth(x,window_len=11,window='hanning',cycle=0):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    #s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    
    if cycle==1:
        s=np.r_[x[-1-window_len+2:],x,x[:window_len-1]]
    else:
        s=np.r_[2*x[0]-x[window_len-1:0:-1],x,2*x[-1]-x[-1:-window_len:-1]]


    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')

    #return y[window_len:-window_len+1]
    return y[window_len//2:-window_len//2+1]


    
    


#######################################################
# 2D Smoothing
#######################################################

from scipy import signal

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-2*size:2*size+1, -2*sizey:2*sizey+1]
    g = np.exp(-(x**2/(2*float(size/2)**2)+y**2/(2*float(sizey/2)**2)))
    return g / g.sum()



#######################################################

def smooth_2d_fast(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    '''
    import cv2
    n = n + (n+1)%2
    im[np.isnan(im)] = 0
    improc = cv2.GaussianBlur(im,(n,n),0)
    return(improc)
    '''
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='same')
    return(improc)



#######################################################

def smooth_2d(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """

    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='same')
    return(improc)


#######################################################

def nansmooth2d(im, ng) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    if ng==0:
        
        return im
    
    else:

        if ng%2==1: ng=ng+1

        n=2*ng; ny=n

        g = gauss_kern(ng, sizey=ng)

        [M,L] = im.shape

        im_large =  np.zeros((M+int(n)*2,L+int(ny)*2))*np.nan
        im_large[n:-n,ny:-ny] = im

        im_mask =  np.ones((M+int(n)*2,L+int(ny)*2))
        im_mask[np.isnan(im_large)] = np.nan

        im_out = np.zeros((M,L))

        for i in range(M):
            for j in range(L):

                if not np.isnan(im[i,j]):
                    nang = g * im_mask[i:i+2*n+1,j:j+2*ny+1]
                    nang = nang / np.nansum(nang)
                    im_out[i,j] = np.nansum(nang * im_large[i:i+2*n+1,j:j+2*ny+1])

                else:

                    im_out[i,j] = np.nan


    return im_out

#######################################################
  

  
  
    
    
###################################################################################
# Small Functions used to find indices
###################################################################################


def find_val(field,value,init=0,smoothing=0,floor=False):
    """
    Find function y(x) corresponding to field(x,y(x))=value
    
    In case of multiple solutions , choose the one with the smallest dy/dx
    
    floor= True: find closer field value smaller than value
    floor= False: find closer field value
    """
    
    jmean = np.zeros(field.shape[0])
    
    jmean[0] = find_sup(field[0,:],value,init=init)
    
    for i in range(1,field.shape[0],1):
        jmean[i] = find_sup(field[i,:],value,init=jmean[i-1])
    
    if floor:
        for i in range(1,field.shape[0],1):
            if jmean[i]>=0:
                if field[i,np.int(jmean[i])]>value: jmean[i]-=1
            else:
                jmean[i]=0

    #Smoothing
    if smoothing==1: jmean = smooth(jmean,min(20,jmean.shape[0]-1),'flat') 
    
    return jmean
    
    
#####################    
    
    
def find_sup(line,val,init=0):

    line0 = np.nan; i1=np.nan
    if np.isnan(init): init=0
        
    for i in range(int(init),line.shape[0],1):
        if ((line[i]-val)*(line0-val))<0.: 
            i1=i
            break
        else:
            line0=line[i]
                        
    for i in range(int(init),int(np.nanmax([0,2*init-i1])),-1):
        if ((line[i]-val)*(line0-val))<0.: 
            i1=i
            break
        else:
            line0=line[i]  
                        
    return i1
            
        

    
####################################################################################

def find_max(field):
    
    """
    Find function y(x) corresponding to field(x,y(x))=max(field(x,:))
        
    """
    
    jmean = np.zeros(field.shape[0])
    for i in range(1,field.shape[0],1):
        jmean[i] = np.nanargmax(field[i,:])
        
    return jmean  
    
    


        
#######################################################
#Rotation
#######################################################

#the function
def rotate_2d(pts,cnt,ang=np.pi/4):
    '''pts = {} Rotates points(nx2) about center cnt(2) by angle ang(1) in radian'''
    return np.dot(pts-cnt,np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]]))+cnt



#######################################################
#Rotate winds or u,v to lat,lon coord
#######################################################

def rotuv(gname,u,v,ny1=None,ny2=None,nx1=None,nx2=None):

    if isinstance(gname, str):
        print('are you sure you want to use angle?')
        nc = Dataset(gname, 'r', format='NETCDF3_CLASSIC')
        angle = np.array(nc.variables['angle'][ny1:ny2,nx1:nx2])
        nc.close()
    else:
        angle = gname 

    'rotate vectors by geometric angle'
    urot = u*np.cos(angle) - v*np.sin(angle)
    vrot = u*np.sin(angle) + v*np.cos(angle)

    return [urot,vrot]
    
    
    
###################################################################################
# Define transformation matrix for local rotated axes (equivalent to interp.griddata but much faster)
###################################################################################


def transform_grid(jmean,angle,old_Ly,Lya, Lyb=0.):
    '''
    In this version  x-dimensions are supposed to be the same
    '''
    
    if Lyb==0.: Lyb=Lya

    xrot=np.arange(len(jmean)); yrot= np.arange(-Lya,Lyb+1,1.)

    Xrot,Yrot = np.meshgrid(xrot,yrot)
    Xrot,Yrot = Xrot.T,Yrot.T

    for ix in range(Xrot.shape[0]): Yrot[ix,:] = Yrot[ix,:] + jmean[ix]

    newXrot,newYrot = 0.*Xrot,0.*Yrot

    for ix in range(Xrot.shape[0]):
        xycenter = [xrot[ix],jmean[ix]]
        for iy in range(len(yrot)):
            [newXrot[ix,iy],newYrot[ix,iy]] = rotate_2d([Xrot[ix,iy],Yrot[ix,iy]],np.array(xycenter),-1*angle[ix])

    x=np.arange(len(jmean)); y=np.arange(old_Ly); X,Y = np.meshgrid(x,y);  X,Y =X.T,Y.T

    [elem,coef] = get_tri_coef(X,Y,newXrot,newYrot)


    return [elem,coef,newXrot,newYrot]


###################################################################################
# Define transformation matrix for local rotated axes (equivalent to interp.griddata but much faster)
###################################################################################


def transform_grid_zoom(imean,jmean,angle, xxx_todo_changeme, Lya, Lyb=0.):

    (old_Lx,old_Ly) = xxx_todo_changeme
    if Lyb==0.: Lyb=Lya
    
    xrot=imean; yrot= np.arange(-Lya,Lyb+1,1.)

    Xrot,Yrot = np.meshgrid(xrot,yrot)
    Xrot,Yrot = Xrot.T,Yrot.T

    for ix in range(Xrot.shape[0]): Yrot[ix,:] = Yrot[ix,:] + jmean[ix]

    newXrot,newYrot = 0.*Xrot,0.*Yrot

    for ix in range(Xrot.shape[0]):
        xycenter = [xrot[ix],jmean[ix]]
        for iy in range(len(yrot)):
            [newXrot[ix,iy],newYrot[ix,iy]] = rotate_2d([Xrot[ix,iy],Yrot[ix,iy]],np.array(xycenter),-1*angle[ix])

    x=np.arange(old_Lx); y=np.arange(old_Ly); X,Y = np.meshgrid(x,y);  X,Y =X.T,Y.T

    [elem,coef] = get_tri_coef(X,Y,newXrot,newYrot)


    return [elem,coef,newXrot,newYrot]




#######################################################
# Interpolation
#######################################################

from scipy.spatial import Delaunay
import time as tm


def get_tri_coef(X,Y,newX,newY,verbose=0):

    """
    Inputs:
        origin lon and lat 2d arrays (X,Y)
        child lon and lat 2d or 1d arrays (newX,newY)

    Ouputs:
        elem - pointers to 2d gridded data (at lonp,latp locations) from
            which the interpolation is computed (3 for each child point)
        coef - linear interpolation coefficients
    Use:
        To subsequently interpolate data from Fp to Fc, the following
        will work:      Fc  = sum(coef.*Fp(elem),3);  This line  should come in place of all
        griddata calls. Since it avoids repeated triangulations and tsearches (that are done
        with every call to griddata) it should be much faster.
    """
    

    Xp = np.array([X.ravel(), Y.ravel()]).T
    Xc = np.array([newX.ravel(), newY.ravel()]).T


    #Compute Delaunay triangulation
    if verbose==1: tstart = tm.time()
    tri = Delaunay(Xp)
    if verbose==1: print('Delaunay Triangulation', tm.time()-tstart)

    #Compute enclosing simplex and barycentric coordinate (similar to tsearchn in MATLAB)
    npts = Xc.shape[0]
    p = np.zeros((npts,3))

    points = tri.points[tri.vertices[tri.find_simplex(Xc)]]

    if verbose==1: tstart = tm.time()
    for i in range(npts):

        if verbose==1: print(np.float(i)/npts)

        if tri.find_simplex(Xc[i])==-1:  #Point outside triangulation
             p[i,:] = p[i,:] * np.nan

        else:

            if verbose==1: tstart = tm.time()
            A = np.append(np.ones((3,1)),points[i] ,axis=1)
            if verbose==1: print('append A', tm.time()-tstart)

            if verbose==1: tstart = tm.time()
            B = np.append(1., Xc[i])
            if verbose==1: print('append B', tm.time()-tstart)

            if verbose==1: tstart = tm.time()
            p[i,:] = np.linalg.lstsq(A.T,B.T)[0]
            if verbose==1: print('solve', tm.time()-tstart)




    if verbose==1: print('Coef. computation 1', tm.time()-tstart)

    if verbose==1: tstart = tm.time()
    
    if len(newY.shape)==2:
        elem = np.reshape(tri.vertices[tri.find_simplex(Xc)],(newX.shape[0],newY.shape[1],3))
        coef = np.reshape(p,(newX.shape[0],newY.shape[1],3))
    elif len(newY.shape)==1:
        elem = np.reshape(tri.vertices[tri.find_simplex(Xc)],(newX.shape[0],3))
        coef = np.reshape(p,(newX.shape[0],3))
        
    if verbose==1: print('Coef. computation 2', tm.time()-tstart)

    return [elem,coef]
    
    
   
   
   
   
   
   
   
   
 
#######################################################


def angle(jmean,imean=None,smoothing=1):
    '''   
    Get angle for local stream following axe
    '''
    
    angle = np.zeros(len(jmean))
    
    if imean==None: imean = np.arange(len(jmean))

    for i in range(1,len(jmean)-1):
        angle[i] = -1*np.arctan((jmean[i+1]-jmean[i-1])/(imean[i+1]-imean[i-1]))

    angle[0] =angle[1]
    angle[-1] = angle[-2]
    
    if smoothing==1: angle = smooth(angle,min(20,jmean.shape[0]-1),'flat')
    
    return angle
        

    

#######################################################


def proj(u,v,mask,nsmooth=20):
    '''
    project a vector field on the vector normal to a contour 
    '''

    if mask.shape[0]==u.shape[0]+1: mask=rho2psi(mask)
    
    '''
    #Get coordinates of the contour
    joe = np.zeros(mask.shape); joe[mask==1] = 1.;
    cs = py.contour(joe.T,[0.])
    cont = cs.collections[0].get_paths()[0].vertices.astype(np.int)
    '''
    
    #u = rho2v(u); v = rho2u(v); # put on psi grid
    
    #Get coordinates of the contour (1st point inside)
    mask_temp = rho2psi(rho2psi(mask))
    joe = np.zeros(mask_temp.shape); joe[mask_temp==1] = 1.;
    cs = py.contour(joe.T,[0.])
    cont = cs.collections[0].get_paths()[0].vertices.astype(np.int) + 1




    if len(u.shape)==3:

        un = np.zeros((cont.shape[0],u.shape[2]))
        ut = np.zeros((cont.shape[0],u.shape[2]))

        for iz in range(u.shape[2]):
            [un[:,iz],ut[:,iz],cont,angle,orientation] = proj_2d(u[:,:,iz],v[:,:,iz],cont,joe,nsmooth)

    else:

        un = np.zeros(cont.shape[0])
        ut = np.zeros(cont.shape[0])

        [un,ut,cont,angle,orientation] = proj_2d(u,v,cont,joe,nsmooth)


    return [un,ut,cont,angle,orientation]


#######################################################

def proj_2d(u,v,cont,mask,nsmooth=20):
    '''
    uco = u[cont[:,0],cont[:,1]]
    vco = v[cont[:,0],cont[:,1]]

    #Compute angle at each contour point
    angle = np.zeros(cont.shape[0])
    un = np.zeros(cont.shape[0])
    ut = np.zeros(cont.shape[0])
    
    X,Y = np.array(np.meshgrid(np.arange(u.shape[0]),np.arange(u.shape[1])), dtype=float)

    x = X.T[cont[:,0],cont[:,1]]
    y = Y.T[cont[:,0],cont[:,1]]
    orientation = np.zeros(cont.shape[0])

    for i in range(0,cont.shape[0]):

        #if i==0: angle[0] = -1*np.arctan((x[1]-x[-1])/(y[1]-y[-1]))
        #elif i==cont.shape[0]-1: angle[-1] = -1*np.arctan((x[0]-x[-2])/(y[0]-y[-2]))
        #else: angle[i] = -1*np.arctan((x[i+1]-x[i-1])/(y[i+1]-y[i-1]))
        if i==0: angle[0] = -1*np.arctan((y[1]-y[-1])/(x[1]-x[-1]))
        elif i==cont.shape[0]-1: angle[-1] = -1*np.arctan((y[0]-y[-2])/(x[0]-x[-2]))
        else: angle[i] = -1*np.arctan((y[i+1]-y[i-1])/(x[i+1]-x[i-1]))

        if mask[np.min([np.ceil(x[i]-np.sin(angle[i])),mask.shape[0]-1]),np.min([np.ceil(y[i]+np.cos(angle[i])),mask.shape[1]-1])]==1: orientation[i]=-1.
        else: orientation[i] =1.

        un[i] = orientation[i] *(vco[i]*np.cos(angle[i])) - uco[i]*np.sin(angle[i])
        ut[i] = orientation[i] *(vco[i]*np.sin(angle[i])) + uco[i]*np.cos(angle[i])
        #print uco[i], vco[i], angle[i]*180./np.pi, un[i]


    #if len(u.shape)==2: [un,vn]=proj_3d(u,v,mask)
    #elif len(u.shape)==3: v=proj_3d(u,v,mask)

    return [un,ut,cont,angle,orientation]
    '''

    uco = u[cont[:,0],cont[:,1]]
    vco = v[cont[:,0],cont[:,1]]

    #Compute angle at each contour point
    angle = np.zeros(cont.shape[0])
    un = np.zeros(cont.shape[0])
    ut = np.zeros(cont.shape[0])

    X,Y = np.array(np.meshgrid(np.arange(u.shape[0]),np.arange(u.shape[1])), dtype=float)

    if nsmooth>0:
        x = smooth(X.T[cont[:,0],cont[:,1]],nsmooth,'flat',cycle=1)
        y = smooth(Y.T[cont[:,0],cont[:,1]],nsmooth,'flat',cycle=1)
    else:
        x = X.T[cont[:,0],cont[:,1]]
        y = Y.T[cont[:,0],cont[:,1]]     

    orientation = np.zeros(cont.shape[0])

    for i in range(0,cont.shape[0]):

        #if i==0: angle[0] = -1*np.arctan((x[1]-x[-1])/(y[1]-y[-1]))
        #elif i==cont.shape[0]-1: angle[-1] = -1*np.arctan((x[0]-x[-2])/(y[0]-y[-2]))
        #else: angle[i] = -1*np.arctan((x[i+1]-x[i-1])/(y[i+1]-y[i-1]))
        if i==0: angle[0] = -1*np.arctan((y[1]-y[-1])/(x[1]-x[-1]))
        elif i==cont.shape[0]-1: angle[-1] = -1*np.arctan((y[0]-y[-2])/(x[0]-x[-2]))
        else: angle[i] = np.arctan((y[i+1]-y[i-1])/(x[i+1]-x[i-1]))
        test=3 #distance used to test orientation
        if mask[np.min([np.ceil(x[i]-test*np.sin(angle[i])),mask.shape[0]-1]),np.min([np.ceil(y[i]+test*np.cos(angle[i])),mask.shape[1]-1])]==1: orientation[i]=1.
        else: orientation[i] =-1.

        un[i] = orientation[i] *(vco[i]*np.cos(angle[i]) - uco[i]*np.sin(angle[i]))
        ut[i] = (vco[i]*np.sin(angle[i])) + uco[i]*np.cos(angle[i])

    return [un,ut,cont,angle,orientation]
    




  
#######################################################

def moy_h(var, nh) :
    """

    """

    
    if len(var.shape)==3: 
    
        [N,M,L] = var.shape
        varout = np.zeros([N,M,L])[::nh,::nh,:]
        for i in range(varout.shape[0]):
            for j in range(varout.shape[1]):
                for k in range(varout.shape[2]):               
                    varout[i,j,k] = np.mean(var[i*nh:(i+1)*nh-1,j*nh:(j+1)*nh-1,k])
                    
        #varout = var[::nh,::nh,:]
    else:
        
        [N,M] = var.shape
        varout = np.zeros([N,M])[::nh,::nh]
        for i in range(varout.shape[0]):
            for j in range(varout.shape[1]):
                varout[i,j] = np.mean(var[i*nh:(i+1)*nh-1,j*nh:(j+1)*nh-1])
        #varout = var[::nh,::nh]
        
       

    
    return varout


#######################################################


def azimuthalAverage(image, center=None, stddev=False, returnradii=False, return_nr=False,
        binsize=1., weights=None, steps=False, interpnan=False, left=None, right=None):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is
             None, which then uses the center of the image (including
             fractional pixels).
    stddev - if specified, return the azimuthal standard deviation instead of the average
    returnradii - if specified, return (radii_array,radial_profile)
    return_nr   - if specified, return number of pixels per radius *and* radius
    binsize - size of the averaging bin.  Can lead to strange results if
        non-binsize factors are used to specify the center and the binsize is
        too large
    weights - can do a weighted average instead of a simple average if this keyword parameter
        is set.  weights.shape must = image.shape.  weighted stddev is undefined, so don't
        set weights and stddev.
    steps - if specified, will return a double-length bin array and radial
        profile so you can plot a step-form radial profile (which more accurately
        represents what's going on)
    interpnan - Interpolate over NAN values, i.e. bins where there is no data?
        left,right - passed to interpnan; they set the extrapolated values

    If a bin contains NO DATA, it will have a NAN value because of the
    divide-by-sum-of-weights component.  I think this is a useful way to denote
    lack of data, but users let me know if an alternative is prefered...
   
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)
    if center is None:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    r = np.hypot(x - center[0], y - center[1])
    if weights is None:
        weights = np.ones(image.shape)
    elif stddev:
        raise ValueError("Weighted standard deviation is not defined.")
    # the 'bins' as initially defined are lower/upper bounds for each bin
    # so that values will be in [lower,upper)  
    nbins = int(np.round(r.max() / binsize)+1)
    maxbin = nbins * binsize
    bins = np.linspace(0,maxbin,nbins+1)
    # but we're probably more interested in the bin centers than their left or right sides...
    bin_centers = (bins[1:]+bins[:-1])/2.0
    # Find out which radial bin each point in the map belongs to
    whichbin = np.digitize(r.flat,bins)
    # how many per bin (i.e., histogram)?
    # there are never any in bin 0, because the lowest index returned by digitize is 1
    nr = np.bincount(whichbin)[1:]
    # recall that bins are from 1 to nbins (which is expressed in array terms by arange(nbins)+1 or xrange(1,nbins+1) )
    # radial_prof.shape = bin_centers.shape
    if stddev:
        radial_prof = np.array([image.flat[whichbin==b].std() for b in range(1,nbins+1)])
    else:
        # add renormalization (2pi/k/n)
        #norm = np.array([weights.flat[whichbin==b].sum() for b in xrange(1,np.int(nbins)+1)])  
        #norm_k = np.array([2*np.pi*(b) for b in xrange(1,np.int(nbins)+1)])       
        radial_prof = np.array([(image*weights).flat[whichbin==b].sum() / weights.flat[whichbin==b].sum() * 2*np.pi*b for b in range(1,np.int(nbins)+1)])  
        #or not
        #radial_prof = np.array([(image*weights).flat[whichbin==b].sum() for b in xrange(1,np.int(nbins)+1)])        
    #import pdb; pdb.set_trace()
    if interpnan:
        radial_prof = np.interp(bin_centers,bin_centers[radial_prof==radial_prof],radial_prof[radial_prof==radial_prof],left=left,right=right)
    if steps:
        xarr = np.array(list(zip(bins[:-1],bins[1:]))).ravel()
        yarr = np.array(list(zip(radial_prof,radial_prof))).ravel()
        return xarr,yarr
    elif returnradii:
        return bin_centers,radial_prof
    elif return_nr:
        return nr,bin_centers,radial_prof
    else:
        return radial_prof




  





