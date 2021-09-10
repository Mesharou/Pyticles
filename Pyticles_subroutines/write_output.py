# Create or open file
if not os.path.isfile(newfile):
# If file doesn't exist: Create it.

    # Create file
    nc = Dataset(newfile, 'w', format = 'NETCDF4')

    # Write attributes
    nc.w_sed0 = w_sed0
    if not adv3d: nc.depth =  advdepth

    if meanflow: nc.meanflow = 1
    else: nc.meanflow = 0
    if initial_depth: nc.initial_depth = 1
    else: nc.initial_depth = 0

    nc.dfile = dfile

    # particles seeding
    if preserved_meter:
        nc.dx_box = dx_box
        nc.nx_box = nx_box
        nc.ny_box = ny_box
    else:
        nc.dx0 = dx0
        nc.iwd = iwd
        nc.jwd = jwd
        nc.nnx = nnx
        nc.nny = nny 
        try:
            nc.dx_m = dx_m
        except:pass 
    nc.nqmx = nqmx
    nc.nnlev = nnlev
    nc.depths0 = depths0

    nc.description = 'particles tracking'
    nc.simulation = parameters
    nc.sub =  subtstep
    nc.base =  0
    nc.ng =  ng
    if x_periodic: nc.x_periodic =  1
    else: nc.x_periodic =  0
    if y_periodic: nc.y_periodic =  1
    else: nc.y_periodic =  0

    # Options
    nc.restart = int(restart)
    nc.dfile = dfile
    nc.write_uv = int(write_uv)
    nc.adv3d = int(adv3d)
    nc.write_lonlat = int(write_lonlat)
    nc.write_depth = int(write_depth)
    nc.write_topo = int(write_topo)
    nc.write_t = int(write_t)
    nc.write_ts = int(write_ts)
    nc.initial_cond = int(initial_cond)
    nc.initial_depth = int(initial_depth)
    nc.meanflow = int(meanflow)
    nc.x_periodic = int(x_periodic)
    nc.y_periodic = int(y_periodic)
    nc.timestep = timestep
    nc.sedimentation = int(sedimentation)
    if continuous_injection : nc.dt_injection = dt_injection

    # particles seeding 
    nc.ic = ic
    nc.jc = jc
    nc.lev0 = lev0
    nc.lev1 = lev1

    nc.nqmx = nqmx
    if preserved_meter:
        nc.dx_box = dx_box
        nc.nx_box = nx_box
        nc.ny_box = ny_box
    else:
        nc.dx0 = dx0
        nc.iwd = iwd
        nc.jwd = jwd
        nc.nnx = nnx
        nc.nny = nny
    nc.nnlev = nnlev
    nc.depth0 = depths0

    #Dimensions
    nc.createDimension('time', None) 
    nc.createDimension('nq', nq) 
    # Time variables
    nc.createVariable('ocean_time', 'f', ('time',) )
    nc.createVariable('time', 'f', ('time',) )

    # Time variables
    nc.createVariable('px','d',('time','nq',))
    nc.createVariable('py','d',('time','nq',))
    if adv3d: nc.createVariable('pz','d',('time','nq',))
    
    if write_ts:
        nc.createVariable('pt','d',('time','nq',))
        nc.createVariable('ps','d',('time','nq',))
    elif write_t:
        nc.createVariable('pt','d',('time','nq',))

    if write_lonlat:    
        nc.createVariable('plon','d',('time','nq',))
        nc.createVariable('plat','d',('time','nq',))

    if write_depth:    
        nc.createVariable('pdepth','d',('time','nq',))

    if write_topo:    
        nc.createVariable('ptopo','d',('time','nq',))


    if write_uv:
        nc.createVariable('pu','d',('time','nq',))
        nc.createVariable('pv','d',('time','nq',))

    if write_uvw:
        nc.createVariable('pu','d',('time','nq',))
        nc.createVariable('pv','d',('time','nq',))
        nc.createVariable('pw','d',('time','nq',))

    #nc.createVariable('ptemp','d',('time','nq',))
    #nc.createVariable('psalt','d',('time','nq',))
    #nc.createVariable('prho1','d',('time','nq',))


else:
    # If file existes: Open it.
    nc = Dataset(newfile, 'a')


#########################################################
# Write Variables into file
try:
    if  dfile > 0:
        tcf = time - np.floor(time)
    else:
        tcf = time - np.ceil(time)
    
    nc.variables['ocean_time'][itime]= simul.oceantime + tcf * simul.dt 

except:
    print('no simul.oceantime')
    nc.variables['ocean_time'][itime]= time * delt

nc.variables['time'][itime] = time
nc.variables['px'][itime, :] = px
nc.variables['py'][itime, :] = py
if adv3d: nc.variables['pz'][itime, :] = pz

if write_ts:    
    nc.variables['pt'][itime, :] = ptemp
    nc.variables['ps'][itime, :] = psalt
elif write_t:    
    nc.variables['pt'][itime, :] = ptemp

if write_lonlat:    
    nc.variables['plon'][itime, :] = plon
    nc.variables['plat'][itime, :] = plat

if write_depth:    
    nc.variables['pdepth'][itime, :] = pdepth

if write_topo:    
    nc.variables['ptopo'][itime, :] = ptopo

if write_uv:
    nc.variables['pu'][itime, :] = pu
    nc.variables['pv'][itime, :] = pv

if write_uvw:
    nc.variables['pu'][itime, :] = pu
    nc.variables['pv'][itime, :] = pv
    nc.variables['pw'][itime, :] = pw

# Close netcdf file
nc.close()
    
