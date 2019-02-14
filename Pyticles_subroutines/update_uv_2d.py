'''
Compute pu,pv corresponsing to u,v at particles positions for 2d surf 

Input: 
    simul, dx, dy, nx, ny
    px, py
    ...

Outputs:
    pu, pv
    
'''

###################################################################################
#Load u,v

if simul.simul[-4:]=='surf':
    [u,v] = part.get_vel_io_surf(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng,coord=coord)
    if not meanflow and alpha_time != 0:
        simul.update(np.ceil(time))
        [u2,v2] = part.get_vel_io_surf(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng,coord=coord)
        simul.update(np.floor(time))
        u = linear(u,u2,alpha_time)
        v = linear(v,v2,alpha_time)
else:
    [u,v] = part.get_vel_io_2d(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng, advdepth = advdepth,coord=coord)
    if not meanflow and alpha_time != 0:
        simul.update(np.ceil(time))
        [u2,v2] = part.get_vel_io_2d(simul,x_periodic=x_periodic,y_periodic=y_periodic,ng=ng, advdepth = advdepth,coord=coord)
        simul.update(np.floor(time))
        u = linear(u,u2,alpha_time)
        v = linear(v,v2,alpha_time)

###############################################################################
# Get u,v at particles positions

pu[:]=partF.interp_2d_u(px,py,u,ng,nq,i0,j0)
pv[:]=partF.interp_2d_v(px,py,v,ng,nq,i0,j0)

###################################################################################

del u,v

