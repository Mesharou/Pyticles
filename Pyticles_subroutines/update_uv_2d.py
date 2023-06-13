'''
- add 10-05-2016 JCollin: if advzavg: pu, and pv are vertically averaged over z_thick
                          around advdepth
- 
Compute pu,pv corresponding to u,v at particles positions for 2d surf 

Input: 
    simul, dx, dy, nx, ny
    px, py
    ...

Outputs:
    pu, pv
    
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

if simul.simul[-4:]=='surf' or advsurf:
    [u, v] = part.get_vel_io_surf(simul, x_periodic=x_periodic,\
             y_periodic=y_periodic, ng=ng, coord=coord,\
             u_name = u_name, v_name = v_name)
    if not meanflow and alpha_time != 0:
        simul.update(next_time)
        [u2, v2] = part.get_vel_io_surf(simul, x_periodic=x_periodic,\
                y_periodic=y_periodic, ng=ng, coord=coord,\
                u_name = u_name, v_name = v_name)
        simul.update(prev_time)
        u = linear(u, u2, alpha_time)
        v = linear(v, v2, alpha_time)
### JC 
elif advzavg:
    
    [u, v] = part.get_vel_io_2d_zavg(simul, x_periodic=x_periodic,
            y_periodic=y_periodic, ng=ng, advdepth=advdepth, z_thick=z_thick,
            coord=coord, u_name = u_name, v_name = v_name)
    if not meanflow and alpha_time != 0:
        simul.update(next_time)
        [u2, v2] = part.get_vel_io_2d_zavg(simul, x_periodic=x_periodic,
                  y_periodic=y_periodic, ng=ng, coord=coord, advdepth=advdepth,
                  z_thick=z_thick, u_name = u_name, v_name = v_name)
        simul.update(prev_time)
        u = linear(u, u2, alpha_time)
        v = linear(v, v2, alpha_time)
        if debug_zavg:
            print('-------------------------------------')
            print(f'max u = {np.max(u)}')
            print(f'max v = {np.max(v)}')
### JC
else:
    [u,v] = part.get_vel_io_2d(simul, x_periodic=x_periodic, y_periodic=y_periodic,
            ng=ng, advdepth=advdepth, coord=coord, u_name = u_name, v_name = v_name)
    if not meanflow and alpha_time != 0:
        simul.update(next_time)
        [u2,v2] = part.get_vel_io_2d(simul, x_periodic=x_periodic,
                  y_periodic=y_periodic, ng=ng, advdepth = advdepth,coord=coord,\
                  u_name = u_name, v_name = v_name)
        simul.update(prev_time)
        u = linear(u, u2, alpha_time)
        v = linear(v, v2, alpha_time)

###############################################################################
# Get u,v at particles positions

pu[:] = partF.interp_2d_u(px, py, u, ng, nq, i0, j0)
pv[:] = partF.interp_2d_v(px, py, v, ng, nq, i0, j0)

###################################################################################

del u, v

