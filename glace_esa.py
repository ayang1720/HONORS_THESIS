'''
Alex Yang Honors Thesis
/users/PAS3252/ayang1720/HONORS_THESIS
/fs/ess/PAS2635/LandAir_Predictability
/fs/scratch/PAS3252/yang/HONORS_THESIS
'''

#the more recent files have 3 days integration, the ones before it had 3 hours integration, both are correlations

#did this one already -delay 50 -loop 0 /fs/scratch/PAS3252/yang/HONORS_THESIS/timeseries/*00.png /users/PAS3252/ayang1720/HONORS_THESIS/time_series_july_to_august_pvalues_3_hours.gif

#correlations 3 days do this one next! ->accidentally did 3 hours instead of 3 days for this one -> edited to fix it
#convert -delay 50 -loop 0 /fs/scratch/PAS3252/yang/HONORS_THESIS/timeseries/*ions.png /users/PAS3252/ayang1720/HONORS_THESIS/time_series_july_to_august_correlations_3_hours.gif

#correlations 3 hours, replace the text name then do that one -> this next convert does 3 days -> edited to fix it

#convert -delay 50 -loop 0 /fs/scratch/PAS3252/yang/HONORS_THESIS/timeseries/*correlations_for_3_days.png /users/PAS3252/ayang1720/HONORS_THESIS/time_series_july_to_august_correlations_3_days.gif

#did this one already convert -delay 50 -loop 0 /fs/scratch/PAS3252/yang/HONORS_THESIS/timeseries/*ays.png /users/PAS3252/ayang1720/HONORS_THESIS/time_series_july_to_august_pvalues_3_days.gif

#NOTE: this line next convert -delay 50 -loop 0 /fs/scratch/PAS3252/yang/HONORS_THESIS/timeseries/*ros.png /users/PAS3252/ayang1720/HONORS_THESIS/time_series_july_to_august_correlations_3_days_new.gif

#imports
import pickle
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
from scipy.stats import beta, false_discovery_control

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
    day=time.strftime('%m/%d')
    array_of_global_response_variables=np.zeros((len(responses),number_of_ensembles,ncells))
    paths=np.empty(number_of_ensembles,dtype='object')
    paths_minus_3_hours=np.empty(number_of_ensembles,dtype='object')
    for i in range(number_of_ensembles):
        datetime=time.strftime('%Y-%m-%d_%H.%M.%S')
        paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.'+datetime+'.nc'
        paths_minus_3_hours[i]=paths[i][:-22]+str(time-timedelta(days=3))+'.nc' #replaced with 72 hours
        paths_minus_3_hours[i]=paths_minus_3_hours[i].replace(" ", "_")
        paths_minus_3_hours[i]=paths_minus_3_hours[i].replace(":", ".")
    for i in range(number_of_ensembles):
        fname=paths[i]
        fname2=paths_minus_3_hours[i]
        nc=Dataset(fname)
        nd=Dataset(fname2)
        for j in range(len(response_variables)):
            datetime=time.strftime('%Y-%m-%d_%H.%M.%S')
            datetime=datetime.replace(" ", "_") #7/27/26 edit: replaced with 72 hour
            if (response_variables[j]=='rainnc'): #find just a 3 hour integrated time span
                rainnc_in_the_3_hour_range=(np.squeeze(np.array(nc[response_variables[j]])))-\
                    (np.squeeze(np.array(nd[response_variables[j]])))
                array_of_global_response_variables[j,i,:]=rainnc_in_the_3_hour_range
            else: #finding the average of the response variables over a day
                sum_over_a_day=np.zeros(ncells)
                for k in range(8):
                    hour_step=timedelta(k*3)
                    itime=time.strptime(datetime,'%Y-%m-%d_%H.%M.%S')
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
        print('day: '+day)
        print('ensemble #'+str(i+1))
    return array_of_global_response_variables

#input:
#X: 100 ensemble members' area averaged soil moisture array july 14th (glace mask)
#R: atmospheric state of 100 values on variable datetime for a single cell (~39281 or 40962,100)
#output:
#correlation: a 1d array representing correlation between X/R -> all values between -1 and 1
#note that for rainnc, certain values have been masked out, so array might have size ~39281

def cor(X1d,R2d):
    number_of_ensembles = len(X1d)
    X1d_mean=np.mean(X1d)
    R2d_mean=np.mean(R2d, axis=-1)
    X1d_pert = X1d - X1d_mean
    R2d_pert = (R2d.T - R2d_mean).T
    covariance = np.sum( X1d_pert * R2d_pert, axis=-1) / (number_of_ensembles-1)

    #test with axis=0 just to see what would happen

    std_X = np.std( X1d, ddof=1 )
    std_R = np.std( R2d, ddof=1, axis=-1 )
    correlation = covariance / ( std_X * std_R)

    correlation=np.array(correlation,dtype=float)

    #flag_keep=np.invert(np.isnan(correlation)) #NOTE edited out 7/27 for correlations instead of pvalues
    #correlation=correlation[flag_keep] #these two lines get commented out for correlations

    return correlation

#input: 1d array of correlations
#output: 1d array of p values
#in this case, we're passing in a 2d array and getting out a 2d array (axis 0 is the 5 state variables)
def pvalue(correlation,number_of_ensembles):
    n=number_of_ensembles
    dist = beta(n/2 - 1, n/2 - 1, loc=-1, scale=2)
    r=correlation
    p=2*dist.cdf(-abs(r))
    p=false_discovery_control(p, axis=0, method='by')
    return p

#input: p-values array, lon array, lat array
#output: p-values array, lon array, lat array, with insignificant cells masked out
def significant_cells(pvalues,lons,lats):
    pvalue_mask=(pvalues<=0.05)
    pvalues=pvalues[pvalue_mask]
    lons=lons[pvalue_mask]
    lats=lats[pvalue_mask]
    return pvalues,lons,lats

#input: glaced smois array, 1d array of only response variables, datetime of the simulation, number of ensembles, ncells
#output: 2d array of pvalues of the correlation between the glaced smois array on 7/14 and the date specified
def glace_esa(glaced_smois_array,response_variables,time,number_of_ensembles,ncells):
    global_arrays=construct_array(response_variables,time,number_of_ensembles,ncells) #shape: (# of responses, 100 ensembles, 40962 cells)
    correlations=np.zeros((len(response_variables),ncells)) #(5,40962) correlation between X/R, should be all values between -1,1
    #pvalues=np.zeros((len(response_variables),ncells)) #commented out for correlations
    for i in range(len(response_variables)):
        global_arrays_2d=global_arrays[i,:,:]
        transposed_global_arrays=global_arrays_2d.T
        ncells_for_this_response=cor(glaced_smois_array,transposed_global_arrays).size
        correlations[i,:ncells_for_this_response]=cor(glaced_smois_array,transposed_global_arrays) #remaining will be zeros

    #pvalues=pvalue(correlations,number_of_ensembles) #creates an array of pvalues for a specified time/date
    #don't compute pvalues for correlations

    return correlations #replaced with the correlations, removed the line that flags nan vlaues NOTE edit 7/27
    #return pvalues #this gets replaced when doing correlations

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
glaced_smois_array_new=np.empty(number_of_ensembles,dtype='object')
for i in range(number_of_ensembles):
    glaced_smois_array[i]=smois_array[i,:][glace_mask]
    glaced_smois_array[i]=np.mean(glaced_smois_array[i][glaced_smois_array[i]!=0])
    #NOTE: created a new array without masking out nonzero values of smois
    glaced_smois_array_new[i]=np.mean(glaced_smois_array[i])

esa_dict={}
esa_dict['smois_1d_new']=glaced_smois_array_new
esa_dict['smois_1d']=glaced_smois_array #100 smois values
esa_dict['lat_1d']=latCell #40962 lat cells
esa_dict['lon_1d']=lonCell #40962 lon cells
pickle.dump(esa_dict,open('ESA.pkl','wb'))