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
ncells=40962

#script

#accessing the atmospheric state variables
avg_dict=pickle.load(open('/fs/scratch/PAS3252/yang/HONORS_THESIS/averages_100.pkl','rb'))
avg_dict_2=pickle.load(open('/fs/scratch/PAS3252/yang/HONORS_THESIS/averages_100_2.pkl','rb'))
#variance of ensemble, averaged over space, over 729 increments of 3 hours
Z500_avg=avg_dict['Z500']
Z250_avg=avg_dict['Z250']
rainc_avg=avg_dict_2['rainc']
t2m_avg=avg_dict_2['t2m']
q2_avg=avg_dict_2['q2']
'''
print(np.shape(Z500_avg))
print(np.shape(Z250_avg))
print(np.shape(rainc_avg))
print(np.shape(t2m_avg))
print(np.shape(q2_avg))
'''

#constructing the soil moisture array for the entire globe on july 14th
paths=np.empty(number_of_ensembles,dtype='object')
for i in range(number_of_ensembles):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-07-14_21.00.00.nc'
smois_array=np.zeros((number_of_ensembles,ncells))
for i in range(number_of_ensembles):
    fname=paths[i]
    nc=Dataset(fname)
    smois_array[i,:]=(np.squeeze(np.array(nc['smois'])))[:,0]
    nc.close()
    #array of soil moistures, each with ncells=40962

#constructing the z500 array for the entire globe on aug 1
paths=np.empty(number_of_ensembles,dtype='object')
for i in range(number_of_ensembles):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-07-14_21.00.00.nc'
smois_array=np.zeros((number_of_ensembles,ncells))
for i in range(number_of_ensembles):
    fname=paths[i]
    nc=Dataset(fname)
    smois_array[i,:]=(np.squeeze(np.array(nc['smois'])))[:,0]
    nc.close()
    #array of soil moistures, each with ncells=40962

#constructing the z250 array for the entire globe on aug 1
paths=np.empty(number_of_ensembles,dtype='object')
for i in range(number_of_ensembles):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-07-14_21.00.00.nc'
smois_array=np.zeros((number_of_ensembles,ncells))
for i in range(number_of_ensembles):
    fname=paths[i]
    nc=Dataset(fname)
    smois_array[i,:]=(np.squeeze(np.array(nc['smois'])))[:,0]
    nc.close()
    #array of soil moistures, each with ncells=40962

#constructing the rainc array for the entire globe on aug 1
paths=np.empty(number_of_ensembles,dtype='object')
for i in range(number_of_ensembles):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-07-14_21.00.00.nc'
smois_array=np.zeros((number_of_ensembles,ncells))
for i in range(number_of_ensembles):
    fname=paths[i]
    nc=Dataset(fname)
    smois_array[i,:]=(np.squeeze(np.array(nc['smois'])))[:,0]
    nc.close()
    #array of soil moistures, each with ncells=40962

#constructing the t2m for the entire globe on aug 1
paths=np.empty(number_of_ensembles,dtype='object')
for i in range(number_of_ensembles):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-07-14_21.00.00.nc'
smois_array=np.zeros((number_of_ensembles,ncells))
for i in range(number_of_ensembles):
    fname=paths[i]
    nc=Dataset(fname)
    smois_array[i,:]=(np.squeeze(np.array(nc['smois'])))[:,0]
    nc.close()
    #array of soil moistures, each with ncells=40962

#constructing the q2 array for the entire globe on aug 1
paths=np.empty(number_of_ensembles,dtype='object')
for i in range(number_of_ensembles):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-07-14_21.00.00.nc'
smois_array=np.zeros((number_of_ensembles,ncells))
for i in range(number_of_ensembles):
    fname=paths[i]
    nc=Dataset(fname)
    smois_array[i,:]=(np.squeeze(np.array(nc['smois'])))[:,0]
    nc.close()
    #array of soil moistures, each with ncells=40962

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

#smois:units = "m3 m^{-3}" ;

#conus_smois is an array of length 100 containing the area averaged smois
#for the contiguous united states, 1 for each of the 100 ensemble members
conus_smois=np.zeros(number_of_ensembles)
for i in range(number_of_ensembles):
    print(np.shape(smois_array))
    smois_array[i,:]*=land_mask*conus_mask
    conus_smois[i]=smois_array[i][smois_array[i]!=0].mean()
    print(conus_smois[i])

#print(smois_array[0][smois_array[0]!=0].sum())
#print((smois_array[0]!=0).sum()) #counts the number of nonzero values
#print(smois_array[0][smois_array[0]!=0].mean())



#input: array of some meterological variable with shape (ensemble,ncells)
#output: an area averaged array representing only the land over CONUS
def conusify(array_to_be_conusified):
    for i in range(number_of_ensembles):
        array_to_be_conusified[i,:]*=land_mask*conus_mask
        conusified_array[i]=\
        array_to_be_conusified[i][array_to_be_conusified[i]!=0].mean()
    return conusified_array