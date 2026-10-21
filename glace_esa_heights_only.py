'''
Alex Yang Honors Thesis
/users/PAS3252/ayang1720/HONORS_THESIS
/fs/ess/PAS2635/LandAir_Predictability
/fs/scratch/PAS3252/yang/HONORS_THESIS
'''

#imports
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
from scipy.stats import beta, false_discovery_control

#settings
number_of_ensembles=500
ncells=40962
response_variables=['height_500hPa']
j=0
smois_start_time=datetime.strptime('20210714210000','%Y%m%d%H%M%S') #7/14/21 at 21Z
variables_start_time=datetime.strptime('20210801210000','%Y%m%d%H%M%S') #8/1/21 at 21Z

#input: array (5) of response variables names, that get pulled from the .nc file for august 1st
#input: time, the exact datetime object referring to the date/time the file produces
#input: number_of_ensembles,ncells -> 500,40962 in this case
#output: array (5,100,40962) of 100 ensemble runs of 40962 grid cells for each response variable
def construct_array(time,number_of_ensembles,ncells):
    day=time.strftime('%m/%d')
    array_of_global_response_variables=np.zeros((number_of_ensembles,ncells))
    paths=np.empty(number_of_ensembles,dtype='object')
    for i in range(number_of_ensembles):
        datetime=time.strftime('%Y-%m-%d_%H.%M.%S')
        paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.'+datetime+'.nc'
    for i in range(number_of_ensembles):
        fname=paths[i]
        nc=Dataset(fname)
        sum_over_a_day=np.zeros(ncells)
        for k in range(8):
            hour_step=timedelta(hours=(k+1)*3) #added hours here
            itime=time.strptime(datetime,'%Y-%m-%d_%H.%M.%S')
            itime+=hour_step
            jtime=itime.strftime('%Y-%m-%d_%H.%M.%S')
            fname3=fname[:57]+jtime+'.nc'
            ne=Dataset(fname3)
            array_of_global_response_variables[i,:]=(np.squeeze(np.array(ne['height_500hPa'])))
            sum_over_a_day+=array_of_global_response_variables[i,:]
        sum_over_a_day/=8
        array_of_global_response_variables[i,:]=sum_over_a_day
        nc.close()
        ne.close()
        print('day: '+day)
        print('ensemble #'+str(i+1))
    return array_of_global_response_variables

#covariance across ensembles, not across cells

def cor_pvalues(X1d,R2d):
    number_of_ensembles = len(X1d)
    X1d_mean=np.mean(X1d)
    R2d_mean=np.mean(R2d, axis=-1)
    X1d_pert = X1d - X1d_mean

    #NOTE: 10/19 removed the transpose because of a broadcast error beign thrown
    #R2d_pert = (R2d.T - R2d_mean).T
    R2d_pert = (R2d.T - R2d_mean)

    R2d_pert = (R2d.T - R2d_mean).T

    covariance = np.sum( X1d_pert * R2d_pert, axis=-1) / (number_of_ensembles-1)

    std_X = np.std( X1d, ddof=1 )

    #NOTE: 10/19 changed axis of R2D to fix broadcast error
    #std_R = np.std( R2d, ddof=1, axis=-1 )
    std_R = np.std( R2d, ddof=1, axis=0 )
    std_R = np.std( R2d, ddof=1, axis=-1 )

    correlation = covariance / ( std_X * std_R)

    correlation=np.array(correlation,dtype=float)

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

#NOTE 10/21 update: the new input for time is now just the global arrays 
#input: glaced smois array, 1d array of only response variables, datetime of the simulation, number of ensembles, ncells
#output: 2d array of pvalues of the correlation between the glaced smois array on 7/14 and the date specified
def glace_esa_pvalues(glaced_smois_array,time,number_of_ensembles,ncells):
    global_arrays=time #shape: (500 ensembles, 40962 cells)
    correlations=np.zeros(ncells) #(40962)

    '''#NOTE: edit made 10/19 because a broadcast erorr occurred using a 100 ensemble array with cor_pvalues function
    transposed_global_arrays=global_arrays_2d
    if number_of_ensembles==500:
        transposed_global_arrays=global_arrays_2d.T'''

    transposed_global_arrays=global_arrays.T

    correlations[:]=cor_pvalues(glaced_smois_array,transposed_global_arrays) #remaining will be zeros
    pvalues=pvalue(correlations,number_of_ensembles) #creates an array of pvalues for a specified time/date
    return pvalues