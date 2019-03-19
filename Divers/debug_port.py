# Takes all output from 2 simulation and compare them
# Results are plotted in folder_save  
# !!! nc_file_p2 is a reference file !!!

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

folder_root = '/home/jeremy/Bureau/Data/Pyticles/'
folder_save = folder_root + '/ADV_2D_NO_UV'

adv3d = False

nc_file_p2 = folder_root + 'Port_Test_P2/Case_1_Port_Test_P2_12_1550.nc'
nc_file_p3 = folder_root + '/ADV_2D_NO_UV/Case_1_ADV_2D_NO_UV_adv0000m_12_1550.nc' 
#nc_file_p3 = folder_root + 'Port_Test_P3/Case_1_Port_Test_P3_12_1550.nc'

nc = Dataset(nc_file_p2, 'r')
var_list = nc.variables
px0 = nc.variables['px'][:]
py0 = nc.variables['py'][:]
if adv3d:
    pz0 = nc.variables['pz'][:]

ocean_time0 = nc.variables['ocean_time'][:]
time0 = nc.variables['time'][:]
pt0 = nc.variables['pt'][:]
ps0 = nc.variables['ps'][:]
plon0 = nc.variables['plon'][:]
plat0 = nc.variables['plat'][:]
pdepth0 = nc.variables['pdepth'][:]
ptopo0 = nc.variables['ptopo'][:]
nc.close()


nc = Dataset(nc_file_p3, 'r')
px3 = nc.variables['px'][:]
py3 = nc.variables['py'][:]
if adv3d:
    pz3 = nc.variables['pz'][:]

ocean_time3 = nc.variables['ocean_time'][:]
time3 = nc.variables['time'][:]
pt3 = nc.variables['pt'][:]
ps3 = nc.variables['ps'][:]
plon3 = nc.variables['plon'][:]
plat3 = nc.variables['plat'][:]
pdepth3 = nc.variables['pdepth'][:]
ptopo3 = nc.variables['ptopo'][:]
nc.close()

# ========================================================


# =========================================================
# PLOTS
# PT

def diag_error(var3, var0, save_plot=False, save_name='', npart=10, title=''):
    '''compute max difference between var1 and var2
       plot some time series of both var1 and var2
       plotting npart first partciles

       var 3 : output from python 3
       var 0 : output from python 2
       npart : number of particles to plot
       save_plot : to save the figure
       save_name : name a the figure for saving
    '''
    err_max = (abs(var3 - var0)).max()
    if err_max != 0:
        print(f'Error simuation are nor identical !')
        print(f'Error is {err_max}')
    else:
        pass

    plt.figure
    plt.subplot(211)
    plt.plot(time3, var3[:,:npart])
    plt.title(f'{title}      Error Max = {err_max}')
    plt.ylabel('python 3.6')
    plt.subplot(212)
    plt.plot(time0, var0[:,:npart])
    plt.ylabel('python 2.7')
    plt.xlabel('time')
    if save_plot:
        plt.savefig(save_name)
    else:
        pass
    plt.show()
    return


#TIME
#save_name = folder_root + 'err_time.png'
#diag_error(time3, time0, save_plot=True, save_name=save_name, npart=10)

#OCEAN_TIME
#save_name = folder_root + 'err_ocean_time.png'
#diag_error(ocean_time3, ocean_time0, save_plot=True, save_name=save_name, npart=10)

# PT 
save_name = folder_root + 'err_pt.png' 
diag_error(pt3, pt0, save_plot=True, save_name=save_name, npart=10, title='temp')

# PS 
save_name = folder_root + 'err_ps.png'
diag_error(ps3, ps0, save_plot=True, save_name=save_name, npart=10, title='salinity')

# PX 
save_name = folder_root + 'err_px.png'
diag_error(px3, px0, save_plot=True, save_name=save_name, npart=10, title='px')

# PY 
save_name = folder_root + 'err_py.png'
diag_error(py3, py0, save_plot=True, save_name=save_name, npart=10, title='py')

# PZ
if adv3d:
    save_name = folder_root + 'err_pz.png'
    diag_error(pz3, pz0, save_plot=True, save_name=save_name, npart=10, title='pz')

# PLON 
save_name = folder_root + 'err_plon.png'
diag_error(plon3, plon0, save_plot=True, save_name=save_name, npart=10, title='plon')

# PLAT
save_name = folder_root + 'err_plat.png'
diag_error(plat3, plat0, save_plot=True, save_name=save_name, npart=10, title='plat')

# PDEPTH 
save_name = folder_root + 'err_pdepth.png'
diag_error(pdepth3, pdepth0, save_plot=True, save_name=save_name, npart=10, title='pdepth')

# PTOPO
save_name = folder_root + 'err_ptopo.png'
diag_error(ptopo3, ptopo0, save_plot=True, save_name=save_name, npart=10, title='ptopo')






