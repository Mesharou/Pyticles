'''
Compute pdepth corresponsing to particles depths

Input: 
    simul, dx, dy, dz, nx, ny, nz
    px, py, pz
    ...

Outputs:
    pdepth
    
'''

###################################################################################
#@profile
def interp_3d_depth(subrange):
    
    global px,py,pz,z_w,nq,ng,i0,j0,k0,pdepth

    pdepth[subrange] = partF.interp_3d_w(px[subrange],py[subrange],pz[subrange],
                                         z_w,ng,nq,i0,j0,k0)

###################################################################################
#create shared arrays

z_w = shared_array(nx_s, ny_s, nz+1);

#Load depth

z_w[:] = part.get_depths_w(simul, x_periodic=x_periodic, y_periodic=y_periodic,
                           ng=ng, coord=coord)

# JC dfile
if dfile > 0:
    prev_time = np.floor(time)
    next_time = np.ceil(time)
else:
    prev_time = np.ceil(time)
    next_time = np.floor(time)

if not meanflow and alpha_time != 0:
    simul.update(next_time)
    z_w2 = part.get_depths_w(simul, x_periodic=x_periodic,
                             y_periodic=y_periodic, ng=ng, coord=coord)
    simul.update(prev_time)
    z_w[:] = linear(z_w[:], z_w2, alpha_time)
    del z_w2

###############################################################################
# Get depth at particles positions
###############################################################################

nslice = nq // nproc
remain = nq - nslice * nproc
i_shift = remain * (nslice + 1)

index = np.arange(nq)

subranges = []
procs = []

for i in range(nproc):
    if i < remain:
        subranges.append(index[i*(nslice+1) : np.nanmin([(i+1)*(nslice+1), nq])])
    else:
        j = i - remain
        subranges.append(index[i_shift + j * nslice : np.nanmin([i_shift + (j+1) * nslice, nq])])

    procs.append(mp.Process(target=interp_3d_depth, args=(subranges[i],)))

for p in procs: p.start()
for p in procs: p.join()   

###################################################################################

del z_w

