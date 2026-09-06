#NOTE: testing has revealed that all ocean smois values are equal to 1

from netCDF4 import Dataset
import numpy as np

ncells=40962
number_of_ensembles=500

layer=0
location=1

fname='/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform/history.2025-10-01_06.00.00.nc'
nc=Dataset(fname)
latCell=np.array(nc['latCell']) #in radians
lonCell=np.array(nc['lonCell']) #in radians
latCell*=(180/np.pi) #in degrees
lonCell*=(180/np.pi) #in degrees
great_plains_mask=(latCell>=37)*(latCell<=44)*(lonCell>=360-104)*(lonCell<=360-97)
california_mask=(latCell>=32.5)*(latCell<=42)*(lonCell>=360-124.5)*(lonCell<=360-114.1)
florida_mask=(latCell>=24.4)*(latCell<=31)*(lonCell>=360-87.6)*(lonCell<=360-80)
masks=[great_plains_mask,california_mask,florida_mask]
smois_mask=masks[location]

paths=np.empty(number_of_ensembles,dtype='object')
smois_array=np.zeros((number_of_ensembles,ncells))
for i in range(1):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-07-14_00.00.00.nc'
    fname=paths[i]
    nc=Dataset(fname)
    smois_array[i,:]=(np.squeeze(np.array(nc['smois'])))[:,layer]
    if i==0:
        print(smois_array[i,:])
        print(np.shape(smois_array[i,:]))
        print(smois_array[i,:][smois_mask])
        print(np.shape(smois_array[i,:][smois_mask]))
        mask=smois_array[i,:][smois_mask]!=1
        print(smois_array[i,:][smois_mask][mask])
        print(np.shape(smois_array[i,:][smois_mask][mask]))
    nc.close()
averaged_smois_array=np.zeros(number_of_ensembles)
for i in range(number_of_ensembles):
    averaged_smois_array[i]=np.mean(smois_array[i,:][smois_mask])

location=2
smois_mask=masks[location]

paths=np.empty(number_of_ensembles,dtype='object')
smois_array=np.zeros((number_of_ensembles,ncells))
for i in range(1):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-07-14_00.00.00.nc'
    fname=paths[i]
    nc=Dataset(fname)
    smois_array[i,:]=(np.squeeze(np.array(nc['smois'])))[:,layer]
    if i==0:
        print(smois_array[i,:])
        print(np.shape(smois_array[i,:]))
        print(smois_array[i,:][smois_mask])
        print(np.shape(smois_array[i,:][smois_mask]))
        mask=smois_array[i,:][smois_mask]!=1
        print(smois_array[i,:][smois_mask][mask])
        print(np.shape(smois_array[i,:][smois_mask][mask]))
    nc.close()
averaged_smois_array=np.zeros(number_of_ensembles)
for i in range(number_of_ensembles):
    averaged_smois_array[i]=np.mean(smois_array[i,:][smois_mask])
