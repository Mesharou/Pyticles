'''
Compute ptemp,psalt corresponsing to T,S at particles positions

Input: 
    simul, dx, dy, dz, nx, ny, nz
    px, py, pz
    ...

Outputs:
    ptemp, psalt
    
'''


###################################################################################
#
def interp_3d_ts(subrange):
    
    global px,py,pz,temp,salt,nq,i0,j0,k0,ptemp,psalt
    ptemp[subrange],psalt[subrange]=partF.interp_3d_ts(px[subrange],py[subrange],pz[subrange],temp,salt,ng,nq,i0,j0,k0)

###################################################################################
#@profile
def interp_2d_ts(subrange):
    
    global px,py,temp,salt,nq,i0,j0,ptemp,psalt

    ptemp[subrange],psalt[subrange]=partF.interp_2d_ts(px[subrange],py[subrange],temp,salt,ng,nq,i0,j0)

###################################################################################
#create shared arrays
if adv3d:
    temp = shared_array(nx_s,ny_s,nz); salt = shared_array(nx_s,ny_s,nz)
else:
    temp = shared_array(nx_s,ny_s); salt = shared_array(nx_s,ny_s)

#Load T,S on sigma levels (much faster than accessing nc files through subprocesses)
if adv3d:
    [temp[:],salt[:]] = part.get_ts_io(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng,coord=coord)
    if not meanflow and alpha_time != 0:
        simul.update(np.ceil(time))
        [temp2,salt2] = part.get_ts_io(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng,coord=coord)
        simul.update(np.floor(time))
        temp[:] = linear(temp[:],temp2,alpha_time)
        salt[:] = linear(salt[:],salt2,alpha_time)

elif simul.simul[-4:]=='surf':
    [temp[:],salt[:]] = part.get_ts_io_surf(simul, x_periodic=x_periodic,
            y_periodic=y_periodic, ng=ng, coord=coord)
    if not meanflow and alpha_time != 0:
        simul.update(np.ceil(time))
        [temp2,salt2] = part.get_ts_io_surf(simul, x_periodic=x_periodic,
                y_periodic=y_periodic, ng=ng, coord=coord)
        simul.update(np.floor(time))
        temp[:] = linear(temp[:],temp2,alpha_time)
        salt[:] = linear(salt[:],salt2,alpha_time)

else:
    [temp[:],salt[:]] = part.get_ts_io_2d(simul, x_periodic=x_periodic, 
            y_periodic=y_periodic, ng=ng, advdepth=advdepth, coord=coord)
    if not meanflow and alpha_time != 0:
        simul.update(np.ceil(time))
        [temp2,salt2] = part.get_ts_io_2d(simul, x_periodic=x_periodic,
                y_periodic=y_periodic, ng=ng, advdepth = advdepth,coord=coord)
        simul.update(np.floor(time))
        temp[:] = linear(temp[:],temp2,alpha_time)
        salt[:] = linear(salt[:],salt2,alpha_time)

###############################################################################
# Get T,S at particles positions
###############################################################################

nslice = nq//nproc + 1
subranges=[]
procs = []

for i in range(nproc):
    subranges.append(list(range(i*nslice,np.nanmin([(i+1)*nslice,nq]))))
    if adv3d:
        procs.append(mp.Process(target=interp_3d_ts, args=(subranges[i],)))
    else:
        procs.append(mp.Process(target=interp_2d_ts, args=(subranges[i],)))

for p in procs: p.start()
for p in procs: p.join()   


###################################################################################

del temp,salt

