'''
Compute plon,plat corresponding to particles positions

Input: 
    simul, dx, dy, nx, ny
    px, py
    ...

Outputs:
    plon, plat
    
'''


###################################################################################
#@profile
def interp_2d(subrange):
    
    global px,py,lon,lat,ng,nq,i0,j0,plon,plat,coord

    [ny1,ny2,nx1,nx2] = np.array(coord)

    # no need to periodize x,y, but force ng=0
    plon[subrange] = partF.interp_2d(px[subrange],py[subrange],simul.x[nx1:nx2,ny1:ny2],0,nq,i0,j0)
    plat[subrange] = partF.interp_2d(px[subrange],py[subrange],simul.y[nx1:nx2,ny1:ny2],0,nq,i0,j0)


###############################################################################
# Get T,S at particles positions
###############################################################################

nslice = nq/nproc+1; subranges=[]
for i in range(nproc): subranges.append(range(i*nslice,np.nanmin([(i+1)*nslice,nq])))

procs = [mp.Process(target=interp_2d, args=([subranges[i]])) for i in range(nproc)]

for p in procs: p.start()
for p in procs: p.join()   

###################################################################################


