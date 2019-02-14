#####################

#ipart = [1000]
#filetime = time_int[0]

#####################

font = {'size'   : 8}
plt.rc('font', **font)

#####################

marksize=5.; s=1   
colorpoint = 'black'; colorpointedge='white'
lw = 0.1/np.sqrt(len(ipart))*10

base = -0.5

################################

simul.update(filetime); 
itime=np.abs(simul.dtime) * np.abs(time_int[0] - filetime); 

depth = -15 #np.min([tools.nanmean(pdepth[ipart,itime]),0]);
x0 = tools.nanmean(px[ipart,itime])-base;
y0 = tools.nanmean(py[ipart,itime])-base;
z0 = tools.nanmean(pz[ipart,itime])+0.5;

coord = part.subsection(px[ipart,itime],py[ipart,itime],ny=simul.coordmax[1],nx=simul.coordmax[3],offset=10); 


##########################################################

varname='rho'
varname2='rho'
varname='temp'
###########

if varname=='salt':
    varname_long = r'$S$'
elif varname=='temp':
    varname_long = r'$T$'
elif varname=='pv':
    varname_long = r'$PV$'    
elif varname=='rho':
    varname_long = r'$\rho$'      
elif varname=='w':
    varname_long = 'w'  




print ' ##########################################################'


i0=0
j0=0
k0=0


if np.int(z0)>=60:
    depths=[59,60]
else:
    depths=range(np.int(z0),np.int(z0)+2)

print depths

salt = var('temp',simul,coord=coord, depths=depths).data
u = var('u',simul,coord=coord, depths=depths).data
v = var('v',simul,coord=coord, depths=depths).data

print '__________________'
#temp = var(varname,simul).data
#print 'interp t', partF.interp_3d(px[ipart,itime],py[ipart,itime],pz[ipart,itime],temp,0.,px.shape[0],i0,j0,k0)
print 'pt',pt[ipart,itime]
#print 'pt',np.nansum(salt[np.int(x0)-coord[2]:np.int(x0)+2-coord[2],np.int(y0)-coord[0]:np.int(y0)+2-coord[0],:]*partF.linear_3d(x0-np.int(x0),y0-np.int(y0),z0-depths[0]))
print '__________________'
#utot = var('u',simul).data
#print 'interp_u', partF.interp_3d_u(px[ipart,itime],py[ipart,itime],pz[ipart,itime],utot,0.,px.shape[0],i0,j0,k0)
#print np.int(x0-0.5)-coord[2],np.int(y0)-coord[0],np.int(z0)-1
#print np.int(x0-0.5)+1,np.int(y0)+1,depths[0]
#print 'u de ca',u[np.int(x0-0.5)-coord[2],np.int(y0)-coord[0],0]
print 'pu',pu[ipart,itime]
#print 'pu',np.nansum(u[np.int(x0-0.5)-coord[2]:np.int(x0-0.5)+2-coord[2],np.int(y0)-coord[0]:np.int(y0)+2-coord[0],:]*partF.linear_3d(x0-0.5-np.int(x0-0.5),y0-np.int(y0),z0-depths[0]))
print '__________________'
#vtot = var('v',simul).data
#print 'interp_v', partF.interp_3d_v(px[ipart,itime],py[ipart,itime],pz[ipart,itime],vtot,0.,px.shape[0],i0,j0,k0)
print 'pv',pv[ipart,itime]
#print 'pv',np.nansum(v[np.int(x0)-coord[2]:np.int(x0)+2-coord[2],np.int(y0-0.5)-coord[0]:np.int(y0-0.5)+2-coord[0],:]*partF.linear_3d(x0-np.int(x0),y0-0.5-np.int(y0-0.5),z0-depths[0]))
print ' ##########################################################'

dpx,dpy = 1.,1.
partF.check_mask(mymask,px[ipart,itime],py[ipart,itime],dpx,dpy,0.,10.,i0,j0)
print 'dpx,dpy',dpx,dpy

###########################################################################################################
#Plot horizontal section at the particule position
###########################################################################################################
ax1 = plt.subplot(2,2,1)
##################

print "compute var", simul.simul, depth, coord, varname
salt = var(varname,simul,coord=coord, depths=[depth]).data * simul.mask[coord[2]:coord[3],coord[0]:coord[1]]

#buoy = var('rho1',simul,coord=coord, depths=[depth]).data
u = var('u',simul,coord=coord, depths=[depth]).data
v = var('v',simul,coord=coord, depths=[depth]).data
#w = var('w',simul,coord=coord, depths=[depth]).data


print 'u', u[np.int(x0-0.5)-coord[2]:np.int(x0-0.5)+2-coord[2],np.int(y0)-coord[0]:np.int(y0)+2-coord[0]]
print 'v', v[np.int(x0)-coord[2]:np.int(x0)+2-coord[2],np.int(y0-0.5)-coord[0]:np.int(y0-0.5)+2-coord[0]]
maskzoom = simul.mask[coord[2]:coord[3],coord[0]:coord[1]]
print 'mask', maskzoom[np.int(x0)-coord[2]:np.int(x0)+2-coord[2],np.int(y0)-coord[0]:np.int(y0)+2-coord[0]]
print 'salt', salt[np.int(x0)-coord[2]:np.int(x0)+2-coord[2],np.int(y0)-coord[0]:np.int(y0)+2-coord[0]]

##########################################################

coefx = 1. #/np.mean(simul.pm)/1000.
x =(np.arange(coord[2],coord[3]))*coefx
x_u = (x[1:] + x[:-1] ) *0.5
y =(np.arange(coord[0],coord[1]))*coefx
[x2d,y2d] = np.meshgrid(x,y)
##################

#levelsvar = oldsim.levels('w',salt,100)
#my_cmap=plot.ncview_colormap('jaisnd')

levelsvar = np.linspace(-0.2,0.2,100); #oldsim.levels('lo',u,100)
levelsvaru = levelsvar
my_cmap=plot.ncview_colormap('blu_red')
my_cmap.set_bad('gray', 1)


##################

#plt.contourf(x,y,salt.T,levelsvar,cmap=my_cmap,extend='both'); plt.colorbar(ticks=oldsim.clabels(levelsvar)); 
#plt.pcolormesh(x,y,ma.masked_invalid(salt.T),vmin=levelsvar.min(),vmax=levelsvar.max(),cmap=my_cmap); plt.colorbar(ticks=oldsim.clabels(levelsvar)); 

#u[u==0] = np.nan
#plt.pcolormesh(x_u,y,ma.masked_invalid(u.T),vmin=levelsvar.min(),vmax=levelsvar.max(),cmap=my_cmap); plt.colorbar(ticks=oldsim.clabels(levelsvar)); 

u0 = tools.u2rho(u)
u0[u0==0] = np.nan
plt.pcolormesh(x-0.5,y-0.5,ma.masked_invalid(u0.T),vmin=levelsvar.min(),vmax=levelsvar.max(),cmap=my_cmap); plt.colorbar(ticks=oldsim.clabels(levelsvar)); 

nn=1. #u.shape[0]/10; x[::nn]
plt.quiver(x2d[::nn,::nn],y2d[::nn,::nn],tools.u2rho(u).T[::nn,::nn],tools.v2rho(v).T[::nn,::nn], pivot='mid', color='k', scale=20, width=0.1,
    headwidth=5,headlength=10,headaxislength=6,
    minlength=0,minshaft=2);

#plot particules
if len(ipart)==1:
    plt.plot((px[ipart,:itime]-base)*coefx,(py[ipart,:itime]-base)*coefx,'o', markersize=2., markerfacecolor=colorpoint, markeredgecolor=colorpointedge); 
    plt.plot((px[ipart,itime+1:]-base)*coefx,(py[ipart,itime+1:]-base)*coefx,'o', markersize=1., markerfacecolor=colorpoint, markeredgecolor=colorpointedge);

plt.plot((px[ipart,itime]-base)*coefx,(py[ipart,itime]-base)*coefx,'o', markersize=marksize, markerfacecolor=colorpoint, markeredgecolor=colorpointedge); 

##################
plt.grid(True, axis='both', linestyle='-', color='k')
plt.axis('scaled')
plt.axis([np.min(x), np.max(x), np.min(y), np.max(y)])
#t1 = plt.title( 'particule no : ' + format(ipart[0]) + ' at ' + format([px[ipart,0][0], px[ipart,0][0], pz[ipart,0][0]]), fontsize=4); 
t1 = plt.title('u' + ' at depth ' + r'$z_0$' + ' and time ' + r'$t_0$'); 
#plt.xlabel(r'$x\,(km)$',fontsize='12')
#plt.ylabel(r'$y\,(km)$',fontsize='12')

del u,v



###########################################################################################################
#Plot horizontal section at the particule position
###########################################################################################################
ax2 = plt.subplot(2,2,2)
##################

#u = var('u',simul,coord=coord, depths=[depth]).data
v = var('v',simul,coord=coord, depths=[depth]).data

##########################################################

coefx = 1. #/np.mean(simul.pm)/1000.
x =(np.arange(coord[2],coord[3]))*coefx
y =(np.arange(coord[0],coord[1]))*coefx

x_u = (x[1:] + x[:-1]) *0.5 
y_v = (y[1:] + y[:-1]) *0.5 


##################
#v[v==0] = np.nan
#pv[v==0] = np.nanlt.pcolormesh(x,y_v,ma.masked_invalid(v.T),vmin=levelsvaru.min(),vmax=levelsvaru.max(),cmap=my_cmap); plt.colorbar(ticks=oldsim.clabels(levelsvaru)); 

v = tools.v2rho(v)
v[v==0] = np.nan
plt.pcolormesh(x-0.5,y-0.5,ma.masked_invalid(v.T),vmin=levelsvaru.min(),vmax=levelsvaru.max(),cmap=my_cmap); plt.colorbar(ticks=oldsim.clabels(levelsvaru)); 

#plot particules
if len(ipart)==1:
    plt.plot((px[ipart,:itime]-base)*coefx,(py[ipart,:itime]-base)*coefx,'o', markersize=2., markerfacecolor=colorpoint, markeredgecolor=colorpointedge); 
    plt.plot((px[ipart,itime+1:]-base)*coefx,(py[ipart,itime+1:]-base)*coefx,'o', markersize=1., markerfacecolor=colorpoint, markeredgecolor=colorpointedge);

plt.plot((px[ipart,itime]-base)*coefx,(py[ipart,itime]-base)*coefx,'o', markersize=marksize, markerfacecolor=colorpoint, markeredgecolor=colorpointedge); 

##################
plt.grid(True, axis='both', linestyle='-', color='k')
plt.axis('scaled')
plt.axis([np.min(x), np.max(x), np.min(y), np.max(y)])
#t1 = plt.title( 'particule no : ' + format(ipart[0]) + ' at ' + format([px[ipart,0][0], px[ipart,0][0], pz[ipart,0][0]]), fontsize=4); 
t1 = plt.title('v' + ' at depth ' + r'$z_0$' + ' and time ' + r'$t_0$'); 
#plt.xlabel(r'$x\,(km)$',fontsize='12')
#plt.ylabel(r'$y\,(km)$',fontsize='12')

del v


###########################################################################################################
def linterpy(var,y):
    coef = y - np.int(y)
    varout = coef*var[:,2,:] + (1-coef)*var[:,1,:]
    return varout

###########################################################################################################
def linterpx(var,x):
    coef = x - np.int(x)
    varout = coef*var[2,:,:] + (1-coef)*var[1,:,:]
    return varout

###########################################################################################################
###########################################################################################################
#Plot vertical section at the particule position
###########################################################################################################

plt.subplot(2,2,3)

###############################
#COORDINATES FOR VERTICAL SECTION (panel 2):

dx=10 #x-half interval
dy0=100 #we are plotting particules within y0+dy0, y0-dy0

x01, x02 = np.int(np.nanmin(px[ipart,itime])-dx), np.int(np.nanmin([np.nanmax(px[ipart,itime])+dx,coord[3]]))
xvar = np.arange(x01, x02,1)

depths = np.arange(np.nanmin(pdepth[ipart,:])-20,np.nanmin([np.nanmax(pdepth[ipart,:])+50, 1]),3)
depths0 = np.arange(np.nanmin(pdepth[ipart,itime])-20,np.nanmin([np.nanmax(pdepth[ipart,itime])+50, 1]),3)
coord_z = [np.int(y0)-1, np.int(y0)+3, x01, x02]

###############################
# section along x direction   
###############################


w = linterpy(tools.rho2u(var('w',simul,coord=coord_z, depths=depths).data),y0)
u = linterpy(var('u',simul,coord=coord_z, depths=depths).data,y0)
xvar_u = (xvar[1:]+xvar[:-1])/2.

#salt = linterpy(var(varname,simul,coord=coord_z, depths=depths).data,y0)
#plt.contourf(xvar*coefx,depths,salt.T,levelsvar,cmap=my_cmap,extend='both'); 
#plt.pcolormesh(xvar*coefx,depths,ma.masked_invalid(salt.T),vmin=levelsvar.min(),vmax=levelsvar.max(),cmap=my_cmap);

u[u==0] = np.nan
plt.pcolormesh((xvar_u-0.5)*coefx,depths,ma.masked_invalid(u.T),vmin=levelsvaru.min(),vmax=levelsvaru.max(),cmap=my_cmap);

#####################
# Plot vectors
coefw = 0.5*(xvar[-1]-xvar[0])/np.mean(simul.pm)/(np.max(depths)-np.min(depths))
nn=1. #u.shape[0]/1; 
nnz=1. #u.shape[1]/5; 

CS1u=plt.quiver(xvar_u[::nn]*coefx,depths[::nnz],u.T[::nnz,::nn],coefw*w.T[::nnz,::nn],pivot='mid',color='k',scale=2,width=0.1,
    headwidth=5,headlength=10,headaxislength=6,
    minlength=0,minshaft=2);

wname = r'$w = ' + '{0:04}'.format(round(0.5/coefw,3)) + '\,m.s^{-1}$'; uname=r'$ u^{\prime} = 0.5\,m.s^{-1},$'; uwname = uname + wname
plt.quiverkey(CS1u,0.8,1.02,0.5,uwname,fontproperties={'size': '6' }, coordinates='axes',color='k')

#####################
if len(ipart)==1:
    plt.plot((px[ipart,:itime]-base)*coefx,pdepth[ipart,:itime]-base,'o', markersize=1., markerfacecolor=colorpoint, markeredgecolor=colorpoint);
    plt.plot((px[ipart,itime+1:]-base)*coefx,pdepth[ipart,itime+1:]-base,'o', markersize=1., markerfacecolor=colorpoint, markeredgecolor=colorpoint);
    plt.plot((px[ipart,itime]-base)*coefx,pdepth[ipart,itime]-base,'o', markersize=marksize, markerfacecolor=colorpoint, markeredgecolor=colorpointedge);  
else: 
    for iq in ipart:
        if np.abs(py[iq,itime]-base-y0)<dy0:
            plt.plot((px[iq,itime]-base)*coefx,pdepth[iq,itime]-base,'o', markersize=marksize, markerfacecolor=colorpoint, markeredgecolor=colorpoint); 
#####################
    
#annotations
plt.plot([np.min(xvar)*coefx, np.max(xvar)*coefx], [depth, depth],'k--',lw=1.5); 
plt.text(np.min(xvar)*coefx, depth+10,r'$z_0$',fontsize='16',color='k'); 
  
#plt.xlabel(r'$x\,(km)$',fontsize='12')
#plt.ylabel(r'$z\,(m)$',fontsize='12',horizontalalignment='center')

t1 = plt.title(varname_long + ' at ' + r'$y^{\prime}_0$' + ' and time ' + r'$t_0$',horizontalalignment='right' ); 

plt.grid(True, axis='both', linestyle='-', color='k')
plt.axis([np.min(xvar)*coefx, np.max(xvar)*coefx, np.min(depths0), np.max(depths0)+5.])
#plt.axis([285, 305, -150,0])



#####################
plt.subplot(2,2,1)

plt.plot([np.min(x), np.max(x)], [y0*coefx, y0*coefx],'k--',lw=1.5); 
plt.plot([np.min(x), np.max(x)], [(y0-dy0)*coefx, (y0-dy0)*coefx],'w--',lw=0.5); 
plt.plot([np.min(x), np.max(x)], [(y0+dy0)*coefx, (y0+dy0)*coefx],'w--',lw=0.5); 
plt.text(np.min(x), y0*coefx-10,r'$y_0$',fontsize='16',color='k'); 

#####################









###########################################################################################################
###########################################################################################################
#Plot vertical section at the particule position
###########################################################################################################

plt.subplot(2,2,4)

###############################
#COORDINATES FOR VERTICAL SECTION (panel 2):
plt.text(np.min(x), y0*coefx-10,r'$y_0$',fontsize='16',color='k'); 

dy=10 #x-half interval
dx0=100 #we are plotting particules within y0+dy0, y0-dy0



y01, y02 = np.int(np.nanmin(py[ipart,itime])-dy), np.int(np.nanmin([np.nanmax(py[ipart,itime])+dy,coord[3]]))
yvar = np.arange(y01, y02,1)
yvar_v = (yvar[1:]+yvar[:-1])/2.


coord_z = [y01, y02, np.int(x0)-1, np.int(x0)+3]

###############################
# section along y direction   
###############################

#salt = linterpx(var(varname,simul,coord=coord_z, depths=depths).data,x0)
w = linterpx(tools.rho2v(var('w',simul,coord=coord_z, depths=depths).data),x0)
v = linterpx(var('v',simul,coord=coord_z, depths=depths).data,x0)

#plt.contourf(yvar*coefx,depths,salt.T,levelsvar,cmap=my_cmap,extend='both'); 
#plt.pcolormesh(yvar*coefx,depths,ma.masked_invalid(salt.T),vmin=levelsvar.min(),vmax=levelsvar.max(),cmap=my_cmap);
v[v==0] = np.nan
plt.pcolormesh((yvar_v-0.5)*coefx,depths,ma.masked_invalid(v.T),vmin=levelsvaru.min(),vmax=levelsvaru.max(),cmap=my_cmap);

#####################
# Plot vectors
coefw = 1.*(yvar[-1]-yvar[0])/np.mean(simul.pn)/(np.max(depths)-np.min(depths))
nn=1. #u.shape[0]/1; 


CS1u=plt.quiver(yvar_v[::nn]*coefx,depths[::nnz],v.T[::nnz,::nn],coefw*w.T[::nnz,::nn],pivot='mid',color='k',scale=5,width=0.1,
    headwidth=5,headlength=10,headaxislength=6,
    minlength=0,minshaft=2);

wname = r'$w = ' + '{0:04}'.format(round(0.5/coefw,3)) + '\,m.s^{-1}$'; uname=r'$ v^{\prime} = 0.5\,m.s^{-1},$'; uwname = uname + wname
plt.quiverkey(CS1u,0.8,1.02,0.5,uwname,fontproperties={'size': '6' }, coordinates='axes',color='k')

#####################
if len(ipart)==1:
    plt.plot((py[ipart,:itime]-base)*coefx,pdepth[ipart,:itime]-base,'o', markersize=1., markerfacecolor=colorpoint, markeredgecolor=colorpoint);
    plt.plot((py[ipart,itime+1:]-base)*coefx,pdepth[ipart,itime+1:]-base,'o', markersize=1., markerfacecolor=colorpoint, markeredgecolor=colorpoint);
    plt.plot((py[ipart,itime]-base)*coefx,pdepth[ipart,itime]-base,'o', markersize=marksize, markerfacecolor=colorpoint, markeredgecolor=colorpointedge);  
else: 
    for iq in ipart:
        if np.abs(px[iq,itime]-base-y0)<dy0:
            plt.plot((py[iq,itime]-base)*coefx,pdepth[iq,itime]-base,'o', markersize=marksize, markerfacecolor=colorpoint, markeredgecolor=colorpoint); 
        
#####################
    
#annotations
plt.plot([np.min(yvar)*coefx, np.max(yvar)*coefx], [depth, depth],'k--',lw=1.5); 
plt.text(np.min(yvar)*coefx, depth+10,r'$z_0$',fontsize='16',color='k'); 
  
#plt.xlabel(r'$x\,(km)$',fontsize='12')
#plt.ylabel(r'$z\,(m)$',fontsize='12',horizontalalignment='center')

t1 = plt.title(varname_long + ' at ' + r'$y^{\prime}_0$' + ' and time ' + r'$t_0$',horizontalalignment='right' ); 

plt.grid(True, axis='both', linestyle='-', color='k')
plt.axis([np.min(yvar)*coefx, np.max(yvar)*coefx, np.min(depths0), np.max(depths0)+5.])
#plt.axis([285, 305, -150,0])



#####################
plt.subplot(2,2,1)

plt.plot([y0*coefx, y0*coefx], [np.min(y), np.max(y)],'k--',lw=1.5); 
plt.plot( [(y0-dy0)*coefx, (y0-dy0)*coefx],[np.min(x), np.max(x)],'w--',lw=0.5); 
plt.plot( [(y0+dy0)*coefx, (y0+dy0)*coefx],[np.min(x), np.max(x)],'w--',lw=0.5); 
plt.text(y0*coefx-10, np.min(y), r'$y_0$',fontsize='16',color='k');   

#####################









##################

plt.savefig(which  + 'part_' +  '{0:05}'.format(iq) + '_time_'+ '{0:03}'.format(filetime) +'.png', size=None, figure=None, magnification='auto', dpi=250,bbox_inches='tight',bbox_extra_artists=[t1])
plt.clf()













