from datetime import datetime,timedelta
import numpy as np
from netCDF4 import Dataset
import pickle
import sys

#NOTE: created two slurm scripts that run this guy
key_int=int(sys.argv[1])

def construct_array(time,number_of_ensembles,ncells,timedelta_string):
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
            if timedelta_string=='timedelta_hours':
                hour_step=timedelta(hours=(k+1)*3) #added hours here
            if timedelta_string=='timedelta_days':
                hour_step=timedelta(days=(k+1)*3) #days here
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
        print('timedelta argument: '+timedelta_string[10:])
        print('day: '+day)
        print('ensemble #'+str(i+1))
    return array_of_global_response_variables

times=np.array(["20210714000000", "20210715000000", "20210716000000", "20210717000000", "20210718000000",\
                "20210719000000", "20210720000000", "20210721000000", "20210722000000", "20210723000000",\
                "20210724000000", "20210725000000", "20210726000000", "20210727000000", "20210728000000",\
                "20210729000000", "20210730000000", "20210731000000", "20210801000000"])

i=0
datetimes=np.empty(len(times),dtype=object)
for time in times:
    datetimes[i]=datetime.strptime(time,'%Y%m%d%H%M%S')
    i+=1

if key_int==0:
    for time in datetimes:
        statgood_dict={}
        statgood_dict['global_array']=construct_array(time,100,40962,'timedelta_days')
        pickle.dump(statgood_dict,open('/fs/scratch/PAS3252/yang/HONORS_THESIS/statgood_pickle_files/global_arrays/global_array_at_time_'+time.strftime('%m-%d')+'_using_timedelta_24_days.pkl','wb'))

if key_int==1:
    for time in datetimes:
        statgood_dict={}
        statgood_dict['global_array']=construct_array(time,500,40962,'timedelta_hours')
        pickle.dump(statgood_dict,open('/fs/scratch/PAS3252/yang/HONORS_THESIS/statgood_pickle_files/global_arrays/global_array_at_time_'+time.strftime('%m-%d')+'_using_timedelta_24_hours.pkl','wb'))


#dump using 24 hours as a timedelta, and 24 days
#should be 2 pickle files created

#then use these pickle files to correlate against both the 100 member smois array and the 500 member smois array
#when using the 100 member smois array, make sure to slice the global array to only use the first 100 members

#import this pickle file into glace_esa_heights_only.py, to prevent having to cal the construct_array function every time