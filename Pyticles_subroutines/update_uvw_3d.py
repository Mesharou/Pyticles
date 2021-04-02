'''
04-2021 JC - added cartesian parameter key to determine whether pw is omega or
     cartesian vertical velocity

JC added 3D uvw in no surface case
 Question why deleting only u,v,w and not u2,v2,w2

Compute pu,pv,pw corresponsing to u,v,w at particles positions for 2d surf 

- add w_sed0 to pw 

Input: 
    simul, dx, dy, dz, nx, ny, nz
    px, py, pz
    ...

Outputs:
    pu, pv, pw
    
'''

###################################################################################
#Load u,v
# JC dfile
if dfile > 0:
    prev_time = np.floor(time)
    next_time = np.ceil(time)
else:
    prev_time = np.ceil(time)
    next_time = np.floor(time)

if simul.simul[-4:]=='surf':
    [u, v] = part.get_vel_io_surf(simul, x_periodic=x_periodic,
                                 y_periodic=y_periodic, ng=ng, coord=coord)
    if not meanflow and alpha_time != 0:
        simul.update(next_time)
        [u2, v2] = part.get_vel_io_surf(simul, x_periodic=x_periodic, 
                                       y_periodic=y_periodic, ng=ng, coord=coord)
        simul.update(prev_time)
        u = linear(u, u2, alpha_time)
        v = linear(v, v2, alpha_time)
else:
    [u, v, w] = part.get_vel_io(simul, x_periodic=x_periodic, 
                                y_periodic=y_periodic, ng=ng, coord=coord, 
                                cartesian=cartesian)
    if not meanflow and alpha_time != 0:
        simul.update(next_time)
        [u2, v2, w2] = part.get_vel_io(simul, x_periodic=x_periodic,
                        y_periodic=y_periodic, ng=ng, coord=coord,
                        cartesian=cartesian)
        simul.update(prev_time)
        u = linear(u, u2, alpha_time)
        v = linear(v, v2, alpha_time)
        w = linear(w, w2, alpha_time)
###############################################################################
# Get u,v at particles positions

pu[:] = partF.interp_3d_u(px, py, pz, u, ng, nq, i0, j0, k0)
pv[:] = partF.interp_3d_v(px, py, pz, v, ng, nq, i0, j0, k0)
if cartesian:
    # w vertical velocity in cartesian coordinates
    pw[:] = partF.interp_3d(px, py, pz, w, ng, nq, i0, j0, k0) \
            + w_sed0 / 3600 / 24
else:
    # omega vertical velocity
    pw[:] = partF.interp_3d_w(px, py, pz, w, ng, nq, i0, j0, k0)\
            + w_sed0 / 3600 / 24


###################################################################################

del u,v,w

