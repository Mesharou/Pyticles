
if not os.path.isfile(newfile):
# If file doesn't exist: Create it.

    nc = Dataset(newfile, 'w', format = 'NETCDF4')

    #Some parameters
    nc.description = 'particles tracking'
    nc.simulation = parameters
    nc.sub =  subtstep
    nc.base =  0
    nc.ng =  ng
    if x_periodic: nc.x_periodic =  1
    else: nc.x_periodic =  0
    if y_periodic: nc.y_periodic =  1
    else: nc.y_periodic =  0

    if not adv3d: nc.depth =  advdepth

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

    #nc.createVariable('ptemp','d',('time','nq',))
    #nc.createVariable('psalt','d',('time','nq',))
    #nc.createVariable('prho1','d',('time','nq',))


else:
    # If file existes: Open it.
    nc = Dataset(newfile, 'a')
    
#########################################################
    
# Write Variables into file
try:
    nc.variables['ocean_time'][itime]= simul.oceantime + (time - np.floor(time)) * simul.dt
except:
    print('no simul.oceantime')
    nc.variables['ocean_time'][itime]= time * delt

nc.variables['time'][itime]=time
nc.variables['px'][itime,:]=px
nc.variables['py'][itime,:]=py
if adv3d: nc.variables['pz'][itime,:]=pz

if write_ts:    
    nc.variables['pt'][itime,:]=ptemp
    nc.variables['ps'][itime,:]=psalt
elif write_t:    
    nc.variables['pt'][itime,:]=ptemp

if write_lonlat:    
    nc.variables['plon'][itime,:]=plon
    nc.variables['plat'][itime,:]=plat

if write_depth:    
    nc.variables['pdepth'][itime,:]=pdepth

if write_topo:    
    nc.variables['ptopo'][itime,:]=ptopo


if write_uv:
    nc.variables['pu'][itime,:]=pu
    nc.variables['pv'][itime,:]=pv



# Close netcdf file
nc.close()
    
