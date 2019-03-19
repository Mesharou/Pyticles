'''
Compute ptemp corresponsing to T at particles positions

Input: 
    simul, dx, dy, dz, nx, ny, nz
    px, py, pz
    ...

Outputs:
    ptemp
    
'''


###################################################################################
#@profile
def interp_3d_t(subrange):
    
    global px,py,pz,temp,nq,i0,j0,k0,ptemp,psalt

    '''
    px_F = np.asfortranarray(px[subrange])
    py_F = np.asfortranarray(py[subrange])
    pz_F = np.asfortranarray(pz[subrange])
    '''

    ptemp[subrange]=partF.interp_3d(px[subrange],py[subrange],pz[subrange],temp,ng,nq,i0,j0,k0)

    
###################################################################################
#create shared arrays
temp = shared_array(nx_s,ny_s,nz)


#Load T,S on sigma levels (much faster than accessing nc files through subprocesses)
temp[:] = part.get_t_io(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng,coord=coord)
if not meanflow and alpha_time != 0:
    simul.update(np.ceil(time))
    temp2 = part.get_t_io(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng,coord=coord)
    simul.update(np.floor(time))
    temp[:] = linear(temp[:],temp2,alpha_time)


###############################################################################
# Get T,S at particles positions
###############################################################################


nslice = nq/nproc+1; subranges=[]
for i in range(nproc): subranges.append(list(range(i*nslice,np.nanmin([(i+1)*nslice,nq]))))

procs = [mp.Process(target=interp_3d_t, args=([subranges[i]])) for i in range(nproc)]

for p in procs: p.start()
for p in procs: p.join()   



###################################################################################

del temp

