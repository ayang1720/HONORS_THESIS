'''
Alex Yang Honors Thesis Task 5 SU26
/users/PAS3252/ayang1720/HONORS_THESIS
/fs/ess/PAS2635/LandAir_Predictability
/fs/scratch/PAS3252/yang/HONORS_THESIS
'''

#imports
import pickle
import numpy as np
from netCDF4 import Dataset

#settings
number_of_ensembles=100

#script

#constructing the soil moisture array for the entire globe
paths=np.empty(number_of_ensembles,dtype='object')
for i in range(number_of_ensembles):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-07-14_21.00.00.nc'
smois_array=np.empty(number_of_ensembles,dtype='object')
for i in range(number_of_ensembles):
    fname=paths[i]
    nc=Dataset(fname)
    smois_array[i]=(np.squeeze(np.array(nc['smois'])))[:,0]
    nc.close()
    #print(np.shape(smois_array[i]))

#constructing the CONUS mask and land mask
fname='/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform/history.2025-10-01_06.00.00.nc'
nc=Dataset(fname)
latCell=np.array(nc['latCell']) #in radians
lonCell=np.array(nc['lonCell']) #in radians
latCell*=(180/np.pi) #in degrees
lonCell*=(180/np.pi) #in degrees
conus_mask=(latCell>=24.5)*(latCell<=49.4)*(lonCell>=(360-124.8))*(lonCell<=(360-66.9))
fname='/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform/x1.40962.init.nc'
nc=Dataset(fname)
land_mask=np.squeeze(np.array(nc['landmask']))
#print(np.sum(conus_mask))
#print(np.sum(land_mask))
#print(np.shape(conus_mask))
#print(np.shape(smois_array))
#print(np.shape(land_mask))

for i in smois_array: #smois_array has 100 elements (1 for each ensemble), each containing a 1d array
    i*=(land_mask*conus_mask) #restrcting soil moisture array to land and CONUS
    #print(np.sum(i!=0))
    i=i[i!=0].mean() #averages each ensemble over the entire conus domain
    print(i)