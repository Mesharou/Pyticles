# Takes all output from 2 simulation and compare them
# Results are plotted in folder_save  
# !!! nc_file_p2 is a reference file !!!

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
# ========================================================
# USER PARAMETERS 

config_3 = 'ADV_2D_UV'
config_2 = 'ADV_2D_UV'

folder_root = '/home/jeremy/Bureau/Data/Pyticles/'
folder_save = folder_root + config_3
generic = '/sig_vs_500m_' # name for all figs
save_name = folder_save + generic + 'err_px.png'


write_uv = True
adv3d = False
save_plot = True

if not adv3d:
    adv_title3 = '_adv0500m' 
    adv_title2 = '_adv0049sig'
else:
    adv_title3 = ''
    adv_title2 = ''


ncfile_p2 = folder_root + config_2 + '/Case_1_' + config_2 + adv_title3 + '_12_1550.nc'
ncfile_p3 = folder_root + config_3 + '/Case_1_' + config_3 + adv_title2 + '_12_1550.nc' 
#nc_file_p3 = folder_root + 'Port_Test_P3/Case_1_Port_Test_P3_12_1550.nc'

# =========== Fonctions to be in a module ================
def get_var(var, ncfile):
    '''
    Returns the netcdf variable 'var' from ncfilei
    var is a string : name of pyticiles  variable
    return np.array py_var
    use module netCDF4
    use module numpy as np
    '''
    nc = Dataset(ncfile, 'r')
    if var in nc.variables:
         py_var = nc.variables[var][:]
    else:
        py_var = []
        print(f'Error {var} is not found in file {ncfile}')
    nc.close()
    return py_var


def plot_diff(var3, var2, time0=[], save_plot=False, save_name='',
        npart=10, title=''):
    '''compute max difference between var1 and var2
       plot some time series of both var1 and var2
       plotting npart first partciles

       var 3 : output from python 3
       var 2 : output from python 2
       time : time_axis 
       npart : number of particles to plot
       save_plot : to save the figure
       save_name : figure name
    '''
    err_max = (abs(var3 - var2)).max()
    if err_max != 0:
        print(f'Error simuation are nor identical !')
        print(f'Error is {err_max}')
    else:
        pass
    
    if len(time0) == 0:
        time0 = np.arange(len(var3[:,0]))
    else:
        pass
    
    plt.figure
    plt.subplot(211)
    plt.plot(time0, var3[:,:npart])
    plt.title(f'{title}      Error Max = {err_max}')
    plt.ylabel('python 3.6')
    plt.subplot(212)
    plt.plot(time0, var2[:,:npart])
    plt.ylabel('python 2.7')
    plt.xlabel('time')
    if save_plot:
        plt.savefig(save_name)
    else:
        pass
    plt.show()
    return

# ====================================================================

#PX
px3 = get_var('px', ncfile_p3)
px2 = get_var('px', ncfile_p2)
time0 = get_var('time', ncfile_p3)

plot_diff(px3, px2, time0=time0, save_plot=save_plot, save_name='',
        npart=10, title='px')

# PU / PV 
if write_uv:
    pu2 = get_var('pu', ncfile_p2)
    pu3 = get_var('pu', ncfile_p3)
    pv2 = get_var('pv', ncfile_p2)
    pv3 = get_var('pv', ncfile_p3)

    plot_diff(pv3, pv2, time0=time0, save_plot=False, save_name='err_pv', npart=10, title='pv')
    plot_diff(pu3, pu2, time0=time0, save_plot=False, save_name='err_pu', npart=10, title='pu')

#TIME
time3 = get_var('time', ncfile_p3)
time2 = get_var('time', ncfile_p2)

save_name = folder_save+ 'err_time.png'
plt.plot(time3, time2)
plt.show()

# PT 
pt3 = get_var('pt', ncfile_p3)
pt2 = get_var('pt', ncfile_p2)

save_name = folder_save + 'err_pt.png'
plot_diff(pt3, pt2, time0=time0, save_plot=save_plot, save_name=save_name,
        npart=10, title='pt')

# PS 
ps3 = get_var('ps', ncfile_p3)
ps2 = get_var('ps', ncfile_p2)

save_name = folder_save + 'err_ps.png'
plot_diff(ps3, ps2, time0=time0, save_plot=save_plot, save_name=save_name,
        npart=10, title='ps')

# PY 
py3 = get_var('py', ncfile_p3)
py2 = get_var('py', ncfile_p2)

save_name = folder_save + 'err_py.png'
plot_diff(py3, py2, time0=time0, save_plot=save_plot, save_name=save_name,
        npart=10, title='py')

# PZ
if adv3d:
    pz3 = get_var('pz', ncfile_p3)
    pz2 = get_var('pz', ncfile_p2)

    save_name = folder_save + 'err_pz.png'
    plot_diff(pz3, pz2, time0=time0, save_plot=save_plot, save_name=save_name,
            npart=10, title='pz')

# PLON 
plon3 = get_var('plon', ncfile_p3)
plon2 = get_var('plon', ncfile_p2)

save_name = folder_save + 'err_plon.png'
plot_diff(plon3, plon2, time0=time0, save_plot=save_plot, save_name=save_name,
        npart=10, title='plon')

# PLAT
plat3 = get_var('plat', ncfile_p3)
plat2 = get_var('plat', ncfile_p2)

save_name = folder_save + 'err_plat.png'
plot_diff(plat3, plat2, time0=time0, save_plot=save_plot, save_name=save_name,
        npart=10, title='plat')

# PDEPTH
if adv3d:
    pdepth3 = get_var('pdepth', ncfile_p3)
    pdepth2 = get_var('pdepth', ncfile_p2)

    save_name = folder_save + 'err_pdepth.png'
    plot_diff(pdepth3, pdepth2, time0=time0, save_plot=save_plot,
            save_name=save_name, npart=10, title='pdepth')

# PTOPO
ptopo3 = get_var('ptopo', ncfile_p3)
ptopo2 = get_var('ptopo', ncfile_p2)

save_name = folder_save + 'err_ptopo.png'
plot_diff(ptopo3, ptopo2, time0=time0, save_plot=save_plot,
        save_name=save_name, npart=10, title='ptopo')







