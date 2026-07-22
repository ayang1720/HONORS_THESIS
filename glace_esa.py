'''
Alex Yang Honors Thesis Task 6 SU26
/users/PAS3252/ayang1720/HONORS_THESIS
/fs/ess/PAS2635/LandAir_Predictability
/fs/scratch/PAS3252/yang/HONORS_THESIS
NEW 7/20/26
'''

#imports
import pickle
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
from scipy.stats import beta, false_discovery_control
import math

#settings
number_of_ensembles=100
ncells=40962
response_variables=['height_500hPa','height_250hPa','rainnc','t2m','q2']
smois_start_time=datetime.strptime('20210714210000','%Y%m%d%H%M%S') #7/14/21 at 21Z
variables_start_time=datetime.strptime('20210801210000','%Y%m%d%H%M%S') #8/1/21 at 21Z

#input: array (5) of response variables names, that get pulled from the .nc file for august 1st
#input: time, the exact datetime object referring to the date/time the file produces
#input: number_of_ensembles,ncells -> 100,40962 in this case
#output: array (5,100,40962) of 100 ensemble runs of 40962 grid cells for each response variable
def construct_array(responses,time,number_of_ensembles,ncells):
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
#input:
#X: 100 ensemble members' area averaged soil moisture array july 14th (glace mask)
#X: area averaged soil moisture (fixed) should be the same for each of the 40000 cells
#R: atmospheric state of 100 values on august 1st for a single cell
#R: R changes depending on which cell you're examining, size (40962,100)

#output:
#correlation: a 1d array (40962) representing correlation between X/R -> all values between -1 and 1

def cor(X1d,R2d):
    number_of_ensembles = len(X1d)
    X1d_mean=np.mean(X1d)
    R2d_mean=np.mean(R2d, axis=-1)
    X1d_pert = X1d - X1d_mean
    R2d_pert = (R2d.T - R2d_mean).T
    covariance = np.sum( X1d_pert * R2d_pert, axis=-1) / (number_of_ensembles-1)

    std_X = np.std( X1d, ddof=1 )
    std_R = np.std( R2d, ddof=1, axis=-1 )
    correlation = covariance / ( std_X * std_R)

    return correlation

#input: 1d array of size (40000) -> correlations
#output: 1d array of size (40000) -> p values
def pvalue(correlation,number_of_ensembles):
    n=number_of_ensembles
    dist = beta(n/2 - 1, n/2 - 1, loc=-1, scale=2)
    r=correlation
    p=2*dist.cdf(-abs(r))
    p=false_discovery_control(p, axis=0, method='by')
    return p

#script

#constructing the CONUS mask and land mask
fname='/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform/history.2025-10-01_06.00.00.nc'
nc=Dataset(fname)
latCell=np.array(nc['latCell']) #in radians
lonCell=np.array(nc['lonCell']) #in radians
latCell*=(180/np.pi) #in degrees
lonCell*=(180/np.pi) #in degrees

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
    #array of 100 soil moistures, each with ncells=40962

#creates an array of size 100 for each ensemble's glaced smois value
glaced_smois_array=np.empty(number_of_ensembles,dtype='object')
for i in range(number_of_ensembles):
    glaced_smois_array[i]=smois_array[i,:][glace_mask]
    glaced_smois_array[i]=np.mean(glaced_smois_array[i][glaced_smois_array[i]!=0])

#check the following 4 lines for accuracy
#constructing the atmospheric state array for the entire globe on aug 1 (100 arrays of 40962 cells each)
#global_arrays=np.zeros((len(response_variables),number_of_ensembles,ncells))
global_arrays=construct_array(response_variables,variables_start_time,number_of_ensembles,ncells) #shape: (# of responses, 100 ensembles, 40962 cells)
correlations=np.zeros((len(response_variables),ncells)) #(5,40962) correlation between X/R, should be all values between -1,1
pvalues=np.zeros((len(response_variables),ncells))
for i in range(len(response_variables)):
    global_arrays_2d=global_arrays[i,:,:]
    transposed_global_arrays=global_arrays_2d.T
    correlations[i,:]=cor(glaced_smois_array,transposed_global_arrays)
flag_keep = np.invert( np.isnan( correlations[2,:] ) ) #keep this flag, then apply it to the correlations array
correlations[2,:] = [0 if isinstance(x, float) and math.isnan(x) else x for x in correlations[2,:]] #help from ChatGPT, replaces all nan with 0
pvalues=pvalue(correlations,number_of_ensembles) #creates an array of 40962 pvalues for a specified time/date

for i in range(len(response_variables)):
    for j in range(ncells):
        print("correlation for "+response_variables[i]+", cell #"+str(j)+":")
        print(correlations[i,j])
        print("p-value:")
        print(pvalues[i,j])

# Array corr_old contains NaNs

# Flags indicating non-NaN values

# Array of non-Nan Values
#correlations[2,:] = correlations[2,:][ flag_keep  ]
#NOTE: don't apply the flag just yet---apply before plotting

#now, the rainnc correlations have shape 39281 instead of 40962

esa_dict={}
esa_dict['flag_keep']=flag_keep
esa_dict['smois']=glaced_smois_array
esa_dict['glace_mask']=glace_mask
esa_dict['correlations']=correlations
esa_dict['p-values']=pvalues
for i in range(len(response_variables)):
    esa_dict[response_variables[i]]=(global_arrays[i,:,:]) #each with size (100,40962)
pickle.dump(esa_dict,open('ESA.pkl','wb'))