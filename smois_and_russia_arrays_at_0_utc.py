'''
Alex Yang's Honors Thesis
Edit log:
7/30: created .py file, imports, and set up layout
8/2: created masks for all 5 regions of interest
set up the first loop, but still need to edit it to pull the correct arrays from
all the new 500 members -> this time around, only one variable from each region of
interest; the variable changes depending on which of the 5 regions you're examining
8/3: investigated the contribution of zero values of smois to the smois array, and
started creating the arrays of x1 and x2 (single date used for correlation)
8/4: created the loop that created arrays for 14 days for y1, y2, y3 and dumped to pickle
8/5: created the piclke file, saved to scratch directory
'''

#copied the pickle file produced to the scratch directory

'''
Notes to self:
-area average the variables (z500, z250) in those boxes; correlations with 500 ensemble members
-for the case of smois in great plains, continue using the glace mask from before, pull ensemble
data for these three regions and compute correlations as before to create correlation coefficients
'''

#imports
from netCDF4 import Dataset
import numpy as np
import pickle

#settings
number_of_ensembles=500
ncells=40962
variables_of_interest=['x1','x2','y1','y2','y3']
height='height_500hPa' #change this depending on which version of the code you're running
#height='height_250hPa'

#lat,lon, and masks
fname='/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform/history.2025-10-01_06.00.00.nc'
nc=Dataset(fname)
latCell=np.array(nc['latCell']) #in radians
lonCell=np.array(nc['lonCell']) #in radians
latCell*=(180/np.pi) #in degrees
lonCell*=(180/np.pi) #in degrees
smois_mask=(latCell>=37)*(latCell<=44)*(lonCell>=360-104)*(lonCell<=360-97)
russia_mask=(latCell>=70)*(latCell<=90)*(lonCell>=90)*(lonCell<=120) #mask for russia, newly edited

#script

#constructing the soil moisture and russia z500 arrays on july 14th
paths=np.empty(number_of_ensembles,dtype='object')
smois_array=np.zeros((number_of_ensembles,ncells))
russia_array=np.zeros((number_of_ensembles,ncells))
for i in range(number_of_ensembles):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-07-14_00.00.00.nc'
    fname=paths[i]
    nc=Dataset(fname)
    smois_array[i,:]=(np.squeeze(np.array(nc['smois'])))[:,0]
    russia_array[i,:]=(np.squeeze(np.array(nc[height])))
    nc.close()
    #array of 500 soil moistures, each with ncells=40962
    #array of 500 russia z500s, each with ncells=40962
#creates an array of size 500 for area averaged smois in CONUS and z500 in russia
averaged_smois_array=np.zeros(number_of_ensembles)
averaged_russia_array=np.zeros(number_of_ensembles)
for i in range(number_of_ensembles):
    averaged_smois_array[i]=np.mean(smois_array[i,:][smois_mask])
    averaged_russia_array[i]=np.mean(russia_array[i,:][russia_mask])
    #NOTE: make sure you use the right mask for all variables!

esa_dict={}
esa_dict['x1']=averaged_smois_array #(ensembles)
esa_dict['x2']=averaged_russia_array #(ensembles)
pickle.dump(esa_dict,open('/fs/scratch/PAS3252/yang/HONORS_THESIS/x1_and_x2_at_0_utc.pkl','wb'))