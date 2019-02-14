
pm=pm_s
pn=pn_s
timing=subtiming
x_periodic=x_periodic
y_periodic=y_periodic
coord=subcoord


[ny1tot,ny2tot,nx1tot,nx2tot] = simul.coord[0:4]

if timing: tstart2 = tm.time()     

nc = Dataset(simul.ncfile, 'r', format='NETCDF3_CLASSIC')
[ny1,ny2,nx1,nx2] = coord
nz = len(simul.coord[4])

################################
mask = copy(simul.mask)
mask[np.isnan(mask)]=0
################################

u = np.zeros((nx2-nx1-1,ny2-ny1,nz))*np.nan
v = np.zeros((nx2-nx1,ny2-ny1-1,nz))*np.nan

################################

nw = min(ng,nx1); ne = min(ng,nx2tot-nx1tot+2*ng-nx2)
ns = min(ng,ny1); nn = min(ng,ny2tot-ny1tot+2*ng-ny2)


# Fill inside points [if x periodic shift one index right in netcdf file]
if x_periodic: iper=1
else: iper = 0

u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1-ns:ny2-2*ng+nn,nx1+iper-nw:nx2-1+iper-2*ng+ne]))'))

u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn,:] = (u[ng-nw:nx2-nx1-1-ng+ne,ng-ns:ny2-ny1-ng+nn,:].T * (mask[nx1+1-nw:nx2-2*ng+ne,ny1-ns:ny2-2*ng+nn]*mask[nx1-nw:nx2-1-2*ng+ne,ny1-ns:ny2-2*ng+nn]).T).T

# Fill inside points [if y periodic shift one index north in netcdf file]
if y_periodic: jper=1
else: jper = 0

v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny1-ns+jper:ny2-1+jper-2*ng+nn,nx1-nw:nx2-2*ng+ne]))

v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn,:] = (v[ng-nw:nx2-nx1-ng+ne,ng-ns:ny2-ny1-1-ng+nn,:].T * (mask[nx1-nw:nx2-2*ng+ne,ny1+1-ns:ny2-2*ng+nn]*mask[nx1-nw:nx2-2*ng+ne,ny1-ns:ny2-1-2*ng+nn]).T).T

################################
# Filling Ghost points
################################

if nw<ng and x_periodic:
    u[ng-nw-1,ng-ns:ny2-ny1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1-ns:ny2-2*ng+nn,nx1tot]))
    for i in range(1,ng):
        u[ng-nw-1-i,ng-ns:ny2-ny1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1-ns:ny2-2*ng+nn,nx2tot-i]))
    for i in range(ng):
        v[ng-nw-1-i,ng-ns:ny2-ny1-1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny1-ns+jper:ny2-1+jper-2*ng+nn,nx2tot-1-i]))
    nw=ng 

if ne<ng and x_periodic:
    for i in range(ng):
        u[nx2-nx1-1-ng+ne+i,ng-ns:ny2-ny1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1-ns:ny2-2*ng+nn,nx1tot+i]))
    for i in range(ng):
        v[nx2-nx1-ng+ne+i,ng-ns:ny2-ny1-1-ng+nn,:] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny1-ns+jper:ny2-1+jper-2*ng+nn,nx1tot+i]))
    ne=ng

if ns<ng and y_periodic:
    v[ng-nw:nx2-nx1-ng+ne,ng-ns-1,:] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny1tot,nx1-nw:nx2-2*ng+ne]))
    for i in range(1,ng):
        v[ng-nw:nx2-nx1-ng+ne,ng-ns-1-i,:] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny2tot-i,nx1-nw:nx2-2*ng+ne]))
    for i in range(1,ng):
        u[ng-nw:nx2-nx1-1-ng+ne,ng-ns-1-i,:] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny2tot-1-i,nx1+iper-nw:nx2-1+iper-2*ng+ne]))

if nn<ng and y_periodic:
    for i in range(ng):
        v[ng-nw:nx2-nx1-ng+ne,ny2-ny1-1-ng+nn+i,:] = simul.Forder(np.squeeze(nc.variables['v'][simul.infiletime,:,ny1tot+i,nx1-nw:nx2-2*ng+ne]))
    for i in range(1,ng):
        u[ng-nw:nx2-nx1-1-ng+ne,ny2-ny1-ng+nn+i,:] = simul.Forder(np.squeeze(nc.variables['u'][simul.infiletime,:,ny1tot+i,nx1+iper-nw:nx2-1+iper-2*ng+ne]))


################################

if timing: print 'get u,v from file....', tm.time()-tstart2
if timing: tstart2 = tm.time()    

################################


print 'no omega in file, computing'
[z_r,z_w] = part.get_depths(simul,coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng)
w = partF.get_omega(u,v,z_r,z_w,pm,pn)
#################

plt.contourf(w[:,:,-2].T,levels=np.arange(-1,1.1,0.1)*1e-6,cmap=plt.cm.RdYlBu_r,extend='both'); plt.savefig('test.png'); plt.clf()

ip=20 
w[px[ip]-coord[2],py[ip]-coord[0],-2]


plt.pcolormesh(z_w[px[ip]-coord[2]-20:px[ip]-coord[2]+20,py[ip]-coord[0]-20:py[ip]-coord[0]+20,-2].T,cmap=plt.cm.RdYlBu_r,vmin=-10.,vmax=0.); plt.savefig('test.png'); plt.clf()

xx,yy = np.arange(px[ip]-coord[2]-5,px[ip]-coord[2]+5),np.arange(py[ip]-coord[0]-5,py[ip]-coord[0]+5)

plt.contour(xx+0.5,yy,z_w[px[ip]-coord[2]-5:px[ip]-coord[2]+5,py[ip]-coord[0]-5:py[ip]-coord[0]+5,0].T,cmap=plt.cm.RdYlBu_r,levels=np.linspace(-1000,0,10));

plt.pcolormesh(xx,yy,u[px[ip]-coord[2]-5:px[ip]-coord[2]+5,py[ip]-coord[0]-5:py[ip]-coord[0]+5,-2].T,cmap=plt.cm.RdYlBu_r,vmin=-0.5,vmax=.5); plt.colorbar(); plt.savefig('test.png'); plt.clf()

#################
w[np.isnan(w)] = 0.
if x_periodic and nx1<ng and nx2>nx2tot-nx1tot+ng: 
    w[0,:,:] = w[nx2tot,:,:]
    w[-1,:,:] = w[nx1tot+2*ng-1,:,:]
if y_periodic and ny1<ng and ny2>ny2tot-ny1tot+ng: 
    w[:,0,:] = w[:,ny2tot,:]
    w[:,-1,:] = w[:,ny1tot+2*ng-1,:]

nc.close()

if timing: print 'get w from file....', tm.time()-tstart2
if timing: tstart2 = tm.time()    

return u,v,w


