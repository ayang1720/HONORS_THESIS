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

#next: run loop, check the dictionary, if masks were applied to the right arrays, plot, write correlation function
'''

'''
Notes to self:
-area average the variables (z500, z250) in those boxes; correlations with 500 ensemble members
-for the case of smois in great plains, continue using the glace mask from before, pull ensemble
data for these three regions and compute correlations as before to create correlation coefficients
'''

#imports
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import pickle

#settings
number_of_ensembles=500
ncells=40962
itime=datetime.strptime('20210714210000','%Y%m%d%H%M%S') #7/14/21 at 21Z
ftime=datetime.strptime('20210801210000','%Y%m%d%H%M%S') #8/1/21 at 21Z
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
russia_mask=(latCell>=65)*(latCell<=90)*(lonCell>=65)*(lonCell<=90) #mask for russia: 65-180 e 65-90 n
Wcasp_mask=(latCell>=40)*(latCell<=50)*(lonCell>=30)*(lonCell<=50) #mask for NELC 1: 30-50 e, 40-50 n (W of caspian sea)
#NOTE: the correlations for NELC 1 and NELC 2 are opposite in sign!
SEcasp_mask=(latCell>=30)*(latCell<=45)*(lonCell>=50)*(lonCell<=70) #mask for NELC 2: 50-70 e, 30-45 n (SE of caspian sea)
Waus_mask=(latCell>=-25)*(latCell<=-15)*(lonCell>=90)*(lonCell<=110) #mask for Waus: 15-25 s, 90-110 e

#functions

def cor(X1d,R2d): #NOTE: edit this completely later
    number_of_ensembles = len(X1d)
    X1d_mean=np.mean(X1d)
    R2d_mean=np.mean(R2d, axis=-1)
    X1d_pert = X1d - X1d_mean
    R2d_pert = (R2d.T - R2d_mean).T
    covariance = np.sum( X1d_pert * R2d_pert, axis=-1) / (number_of_ensembles-1)

    std_X = np.std( X1d, ddof=1 )
    std_R = np.std( R2d, ddof=1, axis=-1 )
    correlation = covariance / ( std_X * std_R)

    correlation=np.array(correlation,dtype=float)

    return correlation

#script

#constructing the soil moisture and russia z500 arrays on july 14th
paths=np.empty(number_of_ensembles,dtype='object')
smois_array=np.zeros((number_of_ensembles,ncells))
russia_array=np.zeros((number_of_ensembles,ncells))
for i in range(number_of_ensembles):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-07-14_21.00.00.nc'
    fname=paths[i]
    nc=Dataset(fname)
    smois_array[i,:]=(np.squeeze(np.array(nc['smois'])))[:,0]
    russia_array[i,:]=(np.squeeze(np.array(nc[height])))
    nc.close()
    #array of 500 soil moistures, each with ncells=40962
    #array of 500 russia z500s, each with ncells=40962
#creates an array of size 100 for area averaged smois in CONUS and z500 in russia
averaged_smois_array=np.zeros(number_of_ensembles)
averaged_russia_array=np.zeros(number_of_ensembles)
for i in range(number_of_ensembles):
    averaged_smois_array[i]=np.mean(smois_array[i,:][smois_mask])
    averaged_russia_array[i]=np.mean(russia_array[i,:][russia_mask])
    #NOTE: make sure you use the right mask for all variables!

#make one of the inputs in the function call whether to examine z500
time=itime
period=ftime-itime
number_of_days=period.days
print(number_of_days)
days=int(number_of_days)
print(days)
Wcasp_averages=np.zeros((days,number_of_ensembles))
SEcasp_averages=np.zeros((days,number_of_ensembles))
Waus_averages=np.zeros((days,number_of_ensembles))
day_number=0
while(time<=ftime):
    day=time.strftime('%m/%d')

    paths=np.empty(number_of_ensembles,dtype='object')

    for i in range(number_of_ensembles):
        time_string=time.strftime('%Y-%m-%d_%H.%M.%S')
        paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.'+time_string+'.nc'
        fname=paths[i]
        nc=Dataset(fname)

        #finding the average of the response variables over a day
        sum_over_a_day=np.zeros((3))

        for k in range(8):
            hour_step=timedelta(hours=k*3)
            jtime=time.strptime(time_string,'%Y-%m-%d_%H.%M.%S')
            jtime+=hour_step
            ktime=jtime.strftime('%Y-%m-%d_%H.%M.%S')
            fname3=fname[:57]+ktime+'.nc'
            ne=Dataset(fname3)

            s1=(np.squeeze(np.array(ne[height])))
            
            sum_over_a_day[0]+=np.mean(s1[Wcasp_mask])
            sum_over_a_day[1]+=np.mean(s1[SEcasp_mask])
            sum_over_a_day[2]+=np.mean(s1[Waus_mask])

        sum_over_a_day/=8

        Wcasp_averages[day_number,i]=sum_over_a_day[0]
        SEcasp_averages[day_number,i]=sum_over_a_day[1]
        Waus_averages[day_number,i]=sum_over_a_day[2]

        nc.close()
        ne.close()
        print('day: '+day)
        print('ensemble #'+str(i+1))
    time+=timedelta(days=1)
    day_number+=1

esa_dict={}
esa_dict['x1']=averaged_smois_array #(ensembles)
esa_dict['x2']=averaged_russia_array #(ensembles)
esa_dict['y1']=Wcasp_averages #(days,ensembles)
esa_dict['y2']=SEcasp_averages #(days,ensembles)
esa_dict['y3']=Waus_averages #(days,ensembles)
pickle.dump(esa_dict,open('xs_and_ys.pkl','wb'))


'''
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
    return array_of_global_response_variables'''