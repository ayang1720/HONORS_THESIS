'''
Alex Yang's Honors Thesis
Top Layer vs. Bottom Layer Soil
Great Plains, California, Florida
'''

'''
edit log:
9/6/26: just finished writing the script that turns inputs into smois and russia arrays for given layer/location
'''

#imports
from netCDF4 import Dataset
import numpy as np
import sys

#settings
layers=[0,-1]
layer=layers[sys.argv[1]] #a number 0 or -1; 0 for top layer, -1 for bottom layer
locations=['great plains','california','florida']
location=sys.argv[2] #a number 0, 1, or 2; corresponds to the array indices above

ncells=40962
number_of_ensembles=500
number_of_days=19

#lat,lon, and masks
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

russia_mask=(latCell>=70)*(latCell<=90)*(lonCell>=90)*(lonCell<=120)
smois_mask=masks[location]

#constructing the soil moisture and russia z500 arrays on july 14th
paths=np.empty(number_of_ensembles,dtype='object')
smois_array=np.zeros((number_of_ensembles,ncells))
russia_array=np.zeros((number_of_ensembles,ncells))
for i in range(number_of_ensembles):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-07-14_00.00.00.nc'
    fname=paths[i]
    nc=Dataset(fname)
    smois_array[i,:]=(np.squeeze(np.array(nc['smois'])))[:,layer]
    russia_array[i,:]=(np.squeeze(np.array(nc['height_500hPa'])))
    nc.close()
averaged_smois_array=np.zeros(number_of_ensembles)
averaged_russia_array=np.zeros(number_of_ensembles)
for i in range(number_of_ensembles):
    water_mask=smois_array[i,:][smois_mask]!=1 #california and florida lat/lon ranges contain ocean, where smois==1
    averaged_smois_array[i]=np.mean(smois_array[i,:][smois_mask][water_mask])
    averaged_russia_array[i]=np.mean(russia_array[i,:][russia_mask])