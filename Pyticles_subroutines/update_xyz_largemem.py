# - 2021-04 Added cartesian = False parameter to force advection on sigma levels
###################################################################################
#Read velocities
###################################################################################
# variables for CFL condition
from input_file import inline_cfl, umax, vmax, wmax, dzmin, cmax, klim, remove,\
                       below, debug, horizontal_diffusion, Khdiff, vertical_diffusion, Kzdiff


debug_time = False

if sedimentation or sedimentation_only: 
    #w_sed0= -25 not supposed to be defined here but in pyticles
    w_sed = w_sed0/(3600.*24.)
    #print(' ')
    #print(' ===========> Vitesse de sedimentation :')
    #print((' w(m/d), w(m/sec) = ',w_sed0, w_sed))
    #print(' ')

timing = True #timing for this subroutine

tstart = tm.time()  

###################################################################################
# Create variables here to be able to release memory completely

nx_s, ny_s = subcoord[3]-subcoord[2], subcoord[1]-subcoord[0]
i0 = subcoord[2]
j0 = subcoord[0]

if adv3d:
    u = shared_array(nx_s-1, ny_s, nz, 2)
    v = shared_array(nx_s, ny_s-1, nz, 2)  
    w = shared_array(nx_s, ny_s, nz+1, 2)  
    dz = shared_array(nx_s, ny_s, nz, 2)
else:
    u = shared_array(nx_s-1, ny_s, 2)
    v = shared_array(nx_s, ny_s-1, 2)  

pm_s = shared_array(nx_s,ny_s)
pn_s = shared_array(nx_s,ny_s)
mask_s = shared_array(nx_s,ny_s)

print(('Create u,v,w shared arrays...............', tm.time()-tstart))
tstart = tm.time()   

#######################

subtiming = True

# Total number of time steps:
istep = shared_array(1,prec='int',value=-1)

# Index of the  previous (itim[0]) and next(itim[1]) time-step for u,v,w,dz
itim = shared_array(2, prec='int', value=[0, 1])


###################################################################################
# LOADING VELOCITY FIELDS
###################################################################################

tim0 = simul.oceantime

if debug_time:
    print('-->')
    print('u0 time :', simul.time)

################
# Initialization
################

if np.isnan(pm_s[0,0]):
    pm_s[:] = part.periodize2d_fromvar(simul, simul.pm, coord=subcoord,
            x_periodic=x_periodic, y_periodic=y_periodic, ng=ng)
    pn_s[:] = part.periodize2d_fromvar(simul, simul.pn, coord=subcoord,
            x_periodic=x_periodic, y_periodic=y_periodic, ng=ng)
    mask_s[:] = part.periodize2d_fromvar(simul, maskrho, coord=subcoord,
            x_periodic=x_periodic, y_periodic=y_periodic, ng=ng)
    if subtiming: print(('Get pm,pn .....................', tm.time()-tstart))
    if subtiming: tstart = tm.time()
    
    ############################################################################
    # Load (u,v,w) on sigmal-levels at time-step t
    if adv3d:
        [u[:,:,:,itim[0]], v[:,:,:,itim[0]], w[:,:,:,itim[0]]] = \
                part.get_vel_io(simul, pm=pm_s, pn=pn_s, timing=subtiming,
                x_periodic=x_periodic, y_periodic=y_periodic, ng=ng,
                coord=subcoord, cartesian=False, u_name = u_name, v_name = v_name)
        if sedimentation:
            w[:, :, :, itim[0]] = w[:, :, :, itim[0]] + w_sed 
        elif sedimentation_only:
            w[:, :, :, itim[0]] = w[:, :, :, itim[0]] * 0 +  w_sed 
    
    elif simul.simul[-4:]=='surf' or advsurf:
        if source == 'analytical':
            [u[:,:,itim[0]],v[:,:,itim[0]]] = part.ana_vel_surf(simul,flow=flow,\
                                                               coord=subcoord,\
                                                        x_periodic=x_periodic,\
                                                        y_periodic=y_periodic,\
                                                                         ng=ng)
        elif source == 'fluid2d':
            [u[:,:,itim[0]],v[:,:,itim[0]]] = part.fluid2d_vel(simul,\
                                                 coord=subcoord,config=config,\
                                                        x_periodic=x_periodic,\
                                                        y_periodic=y_periodic,\
                                                                         ng=ng)
        else:
            [u[:,:,itim[0]], v[:,:,itim[0]]] = part.get_vel_io_surf(simul, pm=pm_s,
                pn=pn_s, timing=subtiming, x_periodic=x_periodic,
                y_periodic=y_periodic, ng=ng, coord=subcoord, u_name = u_name, v_name = v_name)
    ### JC 
    elif advzavg:
        [u[:,:,itim[0]], v[:,:,itim[0]]] = part.get_vel_io_2d_zavg(simul, pm=pm_s,
                pn=pn_s, timing=subtiming, x_periodic=x_periodic,
                y_periodic=y_periodic, ng=ng, advdepth = advdepth, z_thick=z_thick,
                coord=subcoord, u_name = u_name, v_name = v_name)
    ### JC
    else:
        [u[:,:,itim[0]], v[:,:,itim[0]]] = part.get_vel_io_2d(simul, pm=pm_s,
                pn=pn_s, timing=subtiming, x_periodic=x_periodic,
                y_periodic=y_periodic, ng=ng, advdepth = advdepth,
                coord=subcoord, u_name = u_name, v_name = v_name)

    if subtiming: print(('Computing velocity at t1.......', tm.time()-tstart))
    if subtiming: tstart = tm.time()

    if adv3d:
        z_w = part.get_depths_w(simul, x_periodic=x_periodic,
                y_periodic=y_periodic, ng=ng, coord=subcoord)
        dz[:,:,:,itim[0]] = z_w[:,:,1:] - z_w[:,:,:-1]
    
    if subtiming: print(('Computing dz at t1.............', tm.time()-tstart))
    if subtiming: tstart = tm.time()


#######################################################
# Update simul at t+1 - and get u,v,w at time-step t+1
#######################################################

if meanflow:

    try:
        tim0 = tim1
    except:
        pass
    
    tim1 = tim0 + delt
    # No need to update
    if not adv3d:
        [u[:,:,itim[1]],v[:,:,itim[1]]] = [u[:,:,itim[0]],v[:,:,itim[0]]]
    else:
        [u[:,:,:,itim[1]],v[:,:,:,itim[1]],w[:,:,:,itim[1]]]\
         = [u[:,:,:,itim[0]],v[:,:,:,itim[0]],w[:,:,:,itim[0]]]
        dz[:,:,:,itim[1]] = dz[:,:,:,itim[0]]

else:
    if dfile > 0: 
        simul.update(int(np.floor(time) + simul.dtime))
    else:
        simul.update(int(np.ceil(time) + simul.dtime))
    
    if debug_time:
        print('u1 time :', simul.time)
        print('-->')

    # JC
    tim1 = simul.oceantime

    if subtiming: print(('Update simulation..............', tm.time()-tstart))
    if subtiming: tstart = tm.time()

    if adv3d:
        [u[:,:,:,itim[1]], v[:,:,:,itim[1]], w[:,:,:,itim[1]]] = part.get_vel_io(simul,
                pm=pm_s, pn=pn_s, x_periodic=x_periodic, y_periodic=y_periodic, ng=ng,
                coord=subcoord, cartesian=False, u_name = u_name, v_name = v_name)
        if sedimentation:
            w[:,:,:,itim[1]] = w[:,:,:,itim[1]] + w_sed
        elif sedimentation_only:
                w[:, :, :, itim[1]] = w[:,:,:,itim[1]] * 0 + w_sed 

    elif simul.simul[-4:]=='surf' or advsurf:
        
        if source == 'analytical':
            [u[:,:,itim[1]],v[:,:,itim[1]]] = part.ana_vel_surf(simul,flow=flow,\
                                                                coord=subcoord,\
                                                        x_periodic=x_periodic,\
                                                        y_periodic=y_periodic,\
                                                                         ng=ng)
        elif source == 'fluid2d':
            [u[:,:,itim[1]],v[:,:,itim[1]]] = part.fluid2d_vel(simul,\
                                                 coord=subcoord,config=config,\
                                                        x_periodic=x_periodic,\
                                                        y_periodic=y_periodic,\
                                                                         ng=ng)
        else:
            [u[:,:,itim[1]], v[:,:,itim[1]]] = part.get_vel_io_surf(simul, pm=pm_s, pn=pn_s, 
                                        timing=subtiming, x_periodic=x_periodic,
                                            y_periodic=y_periodic, ng=ng, coord=subcoord,\
                                                u_name = u_name, v_name = v_name)
    ### JC 
    elif advzavg:
        [u[:,:,itim[1]], v[:,:,itim[1]]] = part.get_vel_io_2d_zavg(simul, pm=pm_s,
                pn=pn_s, timing=subtiming, x_periodic=x_periodic,
                y_periodic=y_periodic, ng=ng, advdepth = advdepth, z_thick=z_thick,
                coord=subcoord, u_name = u_name, v_name = v_name)
    ### JC
    else:
        [u[:,:,itim[1]], v[:,:,itim[1]]] = part.get_vel_io_2d(simul, pm=pm_s,
                                        pn=pn_s, timing=subtiming,
                                        x_periodic=x_periodic, y_periodic=y_periodic,
                                        ng=ng, advdepth = advdepth, coord=subcoord,\
                                            u_name = u_name, v_name = v_name)

    if subtiming: print(('Computing velocity at t2.......', tm.time()-tstart))
    if subtiming: tstart = tm.time()

    if adv3d:
        z_w = part.get_depths_w(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng,coord=subcoord)
        dz[:,:,:,itim[1]] = z_w[:,:,1:] - z_w[:,:,:-1]; del z_w

    if subtiming: print(('Computing dz at t2.............', tm.time()-tstart))
    if subtiming: tstart = tm.time()


################################################################################### 
# Compute maxvel

if subtiming: tstart = tm.time()
maxvel[0] = np.nanmax(np.abs(u))*1.5
maxvel[1] = np.nanmax(np.abs(v))*1.5
if adv3d:
    maxvel[2] = np.nanmax(np.abs(w))*1.5 
print(('maxvel is', maxvel))
if subtiming: print(('Computing maxvel.............', tm.time()-tstart))

###################################################################################
# set subtime step between two ROMS frames for integration
if not meanflow: delt[0] = simul.dt * np.sign(dfile)
if inline_cfl:
    if umax is None:
       umax = maxvel[0]
    if vmax is None:
       vmax = maxvel[1]
    if wmax is None:
       wmax = maxvel[2]
    
    nsub_steps = part.get_nsub_steps(simul=simul, cmax=cmax, umax=umax, vmax=vmax,
                                     wmax=wmax, dzmin=dzmin, coord=subcoord)
    subtstep = int(nsub_steps * np.abs(dfile))

print(f"--> nsub_steps {nsub_steps}")
dt = delt[0] / nsub_steps
dfct = 1. / nsub_steps

###################################################################################
# Multiprocess for the advance_3d part   
###################################################################################
#@profile

def advance_3d(subrange,out,step):
    
    global px, py, pz, u, v, w, pm_s, pn_s, mask_s, dz, dt, dfct, ng, nq, i0, \
    j0, k0, tim0, delt, subtstep, nxi, nyi, nz, istep, iab, itim, debug_time, \
    remove, klim, below, debug
    
    #from input_file import remove
    #if remove:
        #from input_file import klim
    # If using a Adams-Bashforth method we need to have access to previous vel. values
    if timestep[:2]=='AB': global dpx,dpy,dpz,iab
         
    px_F = np.asfortranarray(px[subrange])
    py_F = np.asfortranarray(py[subrange])
    pz_F = np.asfortranarray(pz[subrange])
    istep_F = istep[0]
    subtime = tim0 + alpha_time * delt[0] 

    if debug_time:
        print("tim0, alpha_time, delt[0]", tim0, alpha_time, delt[0])
        print("subtime", subtime)
 
    if timestep[:2]=='AB': 
        dpx_F = np.asfortranarray(dpx[subrange,:])
        dpy_F = np.asfortranarray(dpy[subrange,:])
        dpz_F = np.asfortranarray(dpz[subrange,:])    
        iab_F = np.asfortranarray(iab[:])

    for it in range(subtstep):
        
        fct = (subtime-tim0) / delt[0]  
        if debug_time:
            print('')
            print('-->')
            print('debug it fct dfct', it, fct, dfct)
            print('debug, dt, subtime', dt, subtime)
            print('')

        istep_F += 1; #print 'istep is', istep, istep_F
        
        ########################
                
        if timestep=='FE':
            partF.timestep_fe(px_F,py_F,pz_F,u,v,w,itim,fct,dfct,pm_s,pn_s,
                              mask_s,dz,dt,ng,nq,i0,j0,k0)
            
        elif timestep=='RK2':  
            partF.timestep_rk2(px_F,py_F,pz_F,u,v,w,itim,fct,dfct,pm_s,pn_s,
                               mask_s,dz,dt,ng,nq,i0,j0,k0)

        elif timestep=='RK4':
            if debug_time: 
                #print('debug px0', px_F[:3])
                #print('debug py0', py_F[:3])
                print('debug others itim,fct,dfct,dt,ng,nq,i0,j0,k0')
                print('debug others', itim,fct,dfct,dt,ng,nq,i0,j0,k0)
            
            partF.timestep_rk4(px_F,py_F,pz_F,u,v,w,itim,fct,dfct,pm_s,pn_s,
            	                   mask_s,dz,dt,ng,nq,i0,j0,k0)

        ########################
        elif timestep[:2]=='AB':
             
            # Test if istep is larger than scheme order (means we have previous values in memory) 
            if istep_F>=int(timestep[-1])-1:
                
                if timestep=='AB2':
                    partF.timestep_ab2(px_F,py_F,pz_F,dpx_F,dpy_F,dpz_F,iab_F,\
                                    u,v,w,itim,fct,dfct,pm_s,pn_s,mask_s,dz,dt,ng,nq,i0,j0,k0)
                    
                elif timestep=='AB3':
                    partF.timestep_ab3(px_F,py_F,pz_F,dpx_F,dpy_F,dpz_F,iab_F,\
                                    u,v,w,itim,fct,dfct,pm_s,pn_s,mask_s,dz,dt,ng,nq,i0,j0,k0)

                elif timestep=='AB4':
                    partF.timestep_ab4(px_F,py_F,pz_F,dpx_F,dpy_F,dpz_F,iab_F,\
                                    u,v,w,itim,fct,dfct,pm_s,pn_s,mask_s,dz,dt,ng,nq,i0,j0,k0)

                elif timestep=='ABM4':
                    partF.timestep_abm4(px_F,py_F,pz_F,dpx_F,dpy_F,dpz_F,iab_F,\
                                    u,v,w,itim,fct,dfct,pm_s,pn_s,mask_s,dz,dt,ng,nq,i0,j0,k0)

            # If not use a RK4 scheme for the first 2 or 4 time steps
            else:
                if subrange[0]==0: print(('istep is', istep_F, ' using RK4 for initialization'))
                (dpx_F[:,iab_F[-1]],dpy_F[:,iab_F[-1]],dpz_F[:,iab_F[-1]]) = \
                                      partF.timestep_rk4(px_F,py_F,pz_F,u,v,w,itim,fct,dfct,\
                                                         pm_s,pn_s,mask_s,dz,dt,ng,nq,i0,j0,k0)         

            iab_F = (iab_F+1)%ab_order  

        ########################  
        
        else:
            raise Exception("no time-stepping scheme specified")

        #add some diffusion:
        if horizontal_diffusion:
            [px_F, py_F] = part.horizontal_diffusion(px_F, py_F, Khdiff, dt, pm_s, pn_s)

        if vertical_diffusion:
            [pdz] = partF.interp_3d(px_F,py_F,pz_F,dz,ng,nq,i0,j0,k0)
            [pz_F] = part.vertical_diffusion(pz_F, Kzdiff, dt, pdz)

        #Remove particles exiting the domain:
        [px_F, py_F, pz_F] = part.cull(px_F, py_F, pz_F, nxi, nyi, nz,
                           x_periodic=x_periodic, y_periodic=y_periodic)
                           
        #Give a kick to particles trapped at surface/bottom       
        [pz_F] = part.kick(pz_F, nz)
        
        
        subtime += dt
        
    # Remove particles below/above (below=True/False) a certain sigma level (klim):
    if remove:
        [px_F, py_F, pz_F] = part.remove_depth(px_F, py_F, pz_F, klim,
                                               below=below, debug=debug)

    # Update the shared arrays 
    if timestep[:2]=='AB': 
        dpx[subrange,:], dpy[subrange,:], dpz[subrange,:] = dpx_F, dpy_F, dpz_F
        out.put(iab_F)
        
    step.put(istep_F)     
        
    px[subrange], py[subrange], pz[subrange] = px_F, py_F, pz_F


###################################################################################
# Multiprocess for the advance_2d part   
###################################################################################
#@profile
#FIXME not sure 2D was tested with dfile = -1/N
def advance_2d(subrange,out,step):
    
    global px, py, u, v, pm_s, pn_s, mask_s, dt, dfct, ng, nq, i0, j0, tim0, delt, subtstep, nxi, nyi, istep, iab, itim
    
    # If using a Adams-Bashforth method we need to have access to previous vel. values
    if timestep[:2]=='AB': global dpx,dpy,iab
    
    px_F = np.asfortranarray(px[subrange])
    py_F = np.asfortranarray(py[subrange])
    istep_F = istep[0]
    subtime = tim0


        
    for it in range(subtstep):
        
        fct = (subtime-tim0)/delt[0]
        istep_F += 1; #print 'istep is', istep, istep_F

        ########################
        if timestep=='RK4':
            partF.timestep2d_rk4(px_F,py_F,u,v,itim,fct,dfct,pm_s,pn_s,mask_s,\
                               dt,ng,nq,i0,j0)  
        ########################  
        else:
            raise Exception("time-stepping scheme not implemented for 2D yet")

        ########################  
        #Remove particles exiting the domain:
        
        [px_F, py_F] = part.cull2d(px_F, py_F, nxi, nyi, x_periodic=x_periodic,
                                   y_periodic=y_periodic)

        subtime += dt


    step.put(istep_F) 
    px[subrange],py[subrange]=px_F,py_F

###################################################################################
# ADVANCE_3D
###################################################################################

nslice = len(subsubrange)//nproc
remain = nq - nslice * nproc
i_shift = remain * (nslice + 1)

subranges=[]
nprocs=[]

for i in range(nproc):
    if i < remain:
        subranges.append(subsubrange[i*(nslice+1) : np.nanmin([(i+1)*(nslice+1), nq])])
    else:
        j = i - remain
        subranges.append(subsubrange[i_shift + j * nslice : np.nanmin([i_shift + (j+1) * nslice,
                         nq])])
    
    if len(subranges[-1]) > 0: nprocs.append(i)

################################################################################
# Run nproc simultaneous processes

procs = []
if len(nprocs)>0:
    out = mp.Queue(); step = mp.Queue()
    for i in range(nproc):
        if adv3d:
            procs.append(mp.Process(target=advance_3d,
                         args=(subranges[i], out, step,)))
        else:
            procs.append(mp.Process(target=advance_2d,
                         args=(subranges[i], out, step,)))
    for p in procs: p.start()
    for p in procs: p.join()   
    if timestep[:2]=='AB': iab[:] = out.get_nowait()
    try:
        istep[0] = step.get_nowait()
    except Queue.Empty:
        raise Exception("advance_3d did not complete")


    ############################################################################
    if timing: print(('Integration between 2 frames...', tm.time()-tstart))
    if timing: tstart = tm.time()

    ############################################################################


