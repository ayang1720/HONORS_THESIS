'''
Alex Yang Honors Thesis Task 6 SU26
/users/PAS3252/ayang1720/HONORS_THESIS
/fs/ess/PAS2635/LandAir_Predictability
/fs/scratch/PAS3252/yang/HONORS_THESIS

most of the code was copied over from task 5, but the new mask used for soil moisure
pulls only smois data from where coupling has shown to be strongest from the glace
experiment done by koster et. al (2006) -> we're going to restrict the conus domain
to just 37-44 N lat and 97-104 W lon (great plains showed most significant in conus)

'''

#imports
import pickle
import numpy as np
from netCDF4 import Dataset

#settings
number_of_ensembles=100
ncells=40962
response_variables=['z500','z250','rainnc','t2m','q2']
response_variables=['height_500hPa','height_250hPa','rainnc','t2m','q2']

#input: array of response variables names, that get pulled from the .nc file for august 1st
#output: array containing 100 ensemble runs of 40962 grid cells for each of the (in this case) 5 response variables
def construct_array(responses):
    array_of_global_response_variables=np.zeros((len(responses),number_of_ensembles,ncells))
    paths=np.empty(number_of_ensembles,dtype='object')
    paths_minus_3_hours=np.empty(number_of_ensembles,dtype='object')
    for i in range(number_of_ensembles):
        paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-08-01_21.00.00.nc'
        paths_minus_3_hours[i]=paths[i][:-16]+'8-01_18.00.00.nc'
    for i in range(number_of_ensembles):
        fname=paths[i]
        fname2=paths_minus_3_hours[i]
        nc=Dataset(fname)
        nd=Dataset(fname2)
        for j in range(len(response_variables)):
            if (response_variables[j]=='rainnc'):
                rainnc_in_the_3_hour_range=(np.squeeze(np.array(nc[response_variables[j]])))-\
                    (np.squeeze(np.array(nd[response_variables[j]])))
                array_of_global_response_variables[j,i,:]=rainnc_in_the_3_hour_range
            else:
                array_of_global_response_variables[j,i,:]=(np.squeeze(np.array(nc[response_variables[j]])))
            #array of z500 heights, each with ncells=40962
            #array of z250 heights, each with ncells=40962
            #array of cumulative precipitation in mm, each with ncells=40962
            #array of 2 meter temperature K, each with ncells=40962
            #array of 2 meter specific humidity kg/kg, each with ncells=40962
        nc.close()
        print('ensemble #'+str(i+1))
    return array_of_global_response_variables

#NOTE: update to the conusify function. conus_mask no longer applies. conus_mask is always true.
#NOTE: updated the function to not area-average, but return a
#input: array of some meterological variable with shape (ensemble,ncells)
#output: an array of atmospheric state values not counting the ocean (ensemble,ncells)
def conusify(array_to_be_conusified,land_mask,conus_mask):
    conusified_array=np.zeros((number_of_ensembles,ncells))
    for i in range(number_of_ensembles):
        array_to_be_conusified[i,:]*=land_mask*conus_mask
        conusified_array[i]=array_to_be_conusified[i]
    return conusified_array

#input:
#X: (soil moisture array of 100 values on july 14th)
#R: (atmosphericic state of 100 values on august 1st)
#output:
#covariance: (a single scalar value representing covariance between X/R)
def cov(X,R):
    X_mean=np.mean(X)
    #print(X_mean)
    R_mean=np.mean(R)
    #print(R_mean)
    summation=0
    for i in range(number_of_ensembles):
        summation+=(X[i]-X_mean)*(R[i]-R_mean)
    #print(summation)
    covariance=summation/(number_of_ensembles-1)
    return covariance

#input: array of response variables
#output: array of correlations between the response variable and the soil moisture
#output is now a 2d array, where axis 0 is every cell in MPAS (40962 of them) and axis 1 is the response variable
def turn_response_variables_into_correlation_array(global_states_array,conus_smois,var_conus_smois,land_mask,conus_mask,number_of_ensembles):
    correlations_array=np.zeros((ncells,global_states_array.shape[0]))
    conusify_array=np.zeros((number_of_ensembles,ncells,global_states_array.shape[0]))
    for i in range(len(global_states_array)): #loop over every response variable
        conusify_array[:,:,i]=conusify(global_states_array[i],land_mask,conus_mask)
        correlations_array[:,i]=cov(conus_smois,conusify_array[:,:,i])\
        /np.sqrt(var_conus_smois)/np.sqrt(np.var(conusify_array[:,:,i],axis=0,ddof=1)) #varying over ensembles
    return correlations_array

#script

#constructing the CONUS mask and land mask
fname='/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform/history.2025-10-01_06.00.00.nc'
nc=Dataset(fname)
latCell=np.array(nc['latCell']) #in radians
lonCell=np.array(nc['lonCell']) #in radians
latCell*=(180/np.pi) #in degrees
lonCell*=(180/np.pi) #in degrees
#conus_mask=(latCell>=24.5)*(latCell<=49.4)*(lonCell>=(360-124.8))*(lonCell<=(360-66.9))
conus_mask=True #change conus_mask to just apply to the entire world
fname='/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform/x1.40962.init.nc'
nc=Dataset(fname)
land_mask=np.squeeze(np.array(nc['landmask']))

#constructing the GLACE mask
glace_mask=(latCell>=37)*(latCell<=44)*(lonCell>=(360-104))*(lonCell<=(360-97))

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

#conus_smois is an array of length 100 containing the area averaged smois
#for the contiguous united states; 1 for each of the 100 ensemble members
conus_smois=(conusify(smois_array,land_mask,glace_mask)) #NOTE: this line edited to use glace mask
area_averaged_smois_array=np.zeros((number_of_ensembles))
for i in range(number_of_ensembles):
    conus_smois_for_this_ensemble=conus_smois[i,:]
    mean_for_this_ensemble=conus_smois_for_this_ensemble[conus_smois_for_this_ensemble!=0].mean()
    area_averaged_smois_array[i]=mean_for_this_ensemble
var_conus_smois=np.var(area_averaged_smois_array,ddof=1)

#constructing the atmospheric state array for the entire globe on aug 1 (100 arrays of 40962 cells each)
#global_arrays=np.zeros((len(response_variables),number_of_ensembles,ncells))
global_arrays=construct_array(response_variables) #shape: (# of responses, 100 ensembles, 40962 cells)
correlations=\
    turn_response_variables_into_correlation_array(global_arrays,conus_smois,var_conus_smois,land_mask,conus_mask,number_of_ensembles)

for i in range(len(response_variables)):
    for j in range(ncells):
        print(correlations[j,i])

esa_dict={}
esa_dict['smois']=conus_smois
esa_dict['land_mask']=land_mask
esa_dict['conus_mask']=conus_mask
esa_dict['glace_mask']=glace_mask
for i in range(len(response_variables)):
    esa_dict[response_variables[i]]=(conusify(global_arrays[i,:,:],land_mask,conus_mask))
pickle.dump(esa_dict,open('ESA.pkl','wb'))