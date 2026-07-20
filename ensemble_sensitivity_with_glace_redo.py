'''
Alex Yang Honors Thesis Task 6 SU26
/users/PAS3252/ayang1720/HONORS_THESIS
/fs/ess/PAS2635/LandAir_Predictability
/fs/scratch/PAS3252/yang/HONORS_THESIS
'''

#imports
import pickle
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta

#settings
number_of_ensembles=100
ncells=40962
response_variables=['height_500hPa','height_250hPa','rainnc','t2m','q2']
smois_start_time=datetime.strptime('20210714210000','%Y%m%d%H%M%S') #7/14/21 at 21Z
variables_start_time=datetime.strptime('20210801210000','%Y%m%d%H%M%S') #8/1/21 at 21Z

#input: array (5) of response variables names, that get pulled from the .nc file for august 1st
#output: array (5,100,40962) of 100 ensemble runs of 40962 grid cells for each response variable
def construct_array(responses):
    array_of_global_response_variables=np.zeros((len(responses),number_of_ensembles,ncells))
    paths=np.empty(number_of_ensembles,dtype='object')
    paths_minus_3_hours=np.empty(number_of_ensembles,dtype='object')
    for i in range(number_of_ensembles):
        time=variables_start_time.strftime('%Y-%m-%d_%H.%M.%S')
        paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.'+time+'.nc'
        paths_minus_3_hours[i]=paths[i][:-22]+str(variables_start_time-timedelta(hours=3))+'.nc'
        paths_minus_3_hours[i]=paths_minus_3_hours[i].replace(" ", "_")
        paths_minus_3_hours[i]=paths_minus_3_hours[i].replace(":", ".")
    for i in range(number_of_ensembles):
        fname=paths[i]
        fname2=paths_minus_3_hours[i]
        nc=Dataset(fname)
        nd=Dataset(fname2)
        for j in range(len(response_variables)):
            time=variables_start_time.strftime('%Y-%m-%d_%H.%M.%S')
            time=time.replace(" ", "_")
            if (response_variables[j]=='rainnc'): #find just a 3 hour integrated time span
                rainnc_in_the_3_hour_range=(np.squeeze(np.array(nc[response_variables[j]])))-\
                    (np.squeeze(np.array(nd[response_variables[j]])))
                array_of_global_response_variables[j,i,:]=rainnc_in_the_3_hour_range
            else: #finding the average of the response variables over a day
                sum_over_a_day=np.zeros(ncells)
                for k in range(8):
                    hour_step=timedelta(k*3)
                    itime=datetime.strptime(time,'%Y-%m-%d_%H.%M.%S')
                    itime+=hour_step
                    jtime=itime.strftime('%Y-%m-%d_%H.%M.%S')
                    fname3=fname[:57]+jtime+'.nc'
                    ne=Dataset(fname3)
                    array_of_global_response_variables[j,i,:]=(np.squeeze(np.array(ne[response_variables[j]])))
                    sum_over_a_day+=array_of_global_response_variables[j,i,:]
                sum_over_a_day/=8
                array_of_global_response_variables[j,i,:]=sum_over_a_day
        nc.close()
        nd.close()
        ne.close()
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
    X_mean=X[X!=0].mean() #should we remove values of smois that are zero? if so, turn this back into 2d array
    #the smois array that gets passed in is shape (100,40962), so shouldn't the outputted also be (40962)?
    #the formula we use for covariance, is it referring to smois averaged over ensemble or over cells?
    #there actually should be covariance between smois and each individual cell of R, so the output needs to be 2d
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

#constructing the CONUS mask and land mask NOTE: edited land_mask to alway be true
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
land_mask=True #land_mask now no longer applies at all

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

#generate smois correlated against smois
smois_smois_correlation=cov(smois_array,smois_array)/np.sqrt(var_conus_smois)/np.sqrt(var_conus_smois)

esa_dict={}
esa_dict['smois']=conus_smois
esa_dict['land_mask']=land_mask
esa_dict['conus_mask']=conus_mask
esa_dict['glace_mask']=glace_mask
esa_dict['correlations']=correlations
esa_dict['smois_smois_correlation']=smois_smois_correlation
for i in range(len(response_variables)):
    esa_dict[response_variables[i]]=(conusify(global_arrays[i,:,:],land_mask,conus_mask))
pickle.dump(esa_dict,open('ESA.pkl','wb'))