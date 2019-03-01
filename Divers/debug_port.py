from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

folder_root = '/home/jeremy/Bureau/Data/Pyticles/'

nc_file_p2 = folder_root + 'Port_Test_P2/Case_1_Port_Test_P2_12_1550.nc'
nc_file_p3 = folder_root + 'Port_Test_P3/Case_1_Port_Test_P2_12_1550.nc'

nc = Dataset(nc_file_p2, 'r')
px0 = nc.variables['px'][:]
py0 = nc.variables['py'][:]
pz0 = nc.variables['pz'][:]
nc.close()


nc = Dataset(nc_file_p2, 'r')
px0_3 = nc.variables['px'][:]
py0_3 = nc.variables['py'][:]
pz0_3 = nc.variables['pz'][:]
nc.close()

print(px0.shape)

plt.plot(px0[:,0] - px0_3[:,0])
#plt.show()


diff_px0 = px0_3 - px0 
diff_py0 = py0_3 - py0 
diff_pz0 = pz0_3 - pz0

Diff_PX = []
Diff_PY = []
Diff_PZ = []

for i in range(100):
    Diff_PX.append(max(diff_px0[:,i]))
    Diff_PY.append(max(diff_py0[:,i]))
    Diff_PZ.append(max(diff_pz0[:,i]))

plt.plot(Diff_PX)
plt.plot(Diff_PY)
plt.plot(Diff_PZ)
plt.show()


for ip in range(100):
    plt.plot(pz0_3[:,ip])
    plt.show()



print(f'Max Error PX = {max(diff_px0))} ') 









