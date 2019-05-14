'''
Compute ptopo corresponding to topography floor depth at particles positions

Input: 
    simul, dx, dy, dz, nx, ny, nz
    px, py, pz
    ...

Outputs:
    ptopo
    
'''


###################################################################################
#@profile
def interp_2d(subrange):
    
    global px,py,topo,ng,nq,i0,j0,ptopo,coord

    [ny1,ny2,nx1,nx2] = np.array(coord)

    ptopo[subrange]=partF.interp_2d(px[subrange],py[subrange],simul.topo[nx1:nx2,ny1:ny2],0,nq,i0,j0)
    
###################################################################################

# no need to periodize topo, but force ng=0
# topo = part.periodize2d_fromvar(simul,simul.topo,coord=coord,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng) 

###############################################################################
# Get topo at particles positions
###############################################################################

nslice = nq//nproc + 1
subranges = []
procs = []

for i in range(nproc):
    subranges.append(list(range(i*nslice,np.nanmin([(i+1)*nslice,nq]))))
    procs.append(mp.Process(target=interp_2d, args=(subranges[i],)))

for p in procs: p.start()
for p in procs: p.join()   


###################################################################################

del topo

