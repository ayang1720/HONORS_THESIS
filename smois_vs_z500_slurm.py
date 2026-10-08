'''
Alex Yang's Honors Thesis
NOTE: the original smois and russia arrays still use the 21-UTC time, instead of the new 0-UTC standard
'''

#NOTE: edited the masks for the 3 locations to be a little more focused in

#imports
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import pickle
import sys

#settings
number_of_ensembles=500
ncells=40962
input_time=sys.argv[1]
ftime=datetime.strptime(input_time,'%Y%m%d%H%M%S') #inputted by user
height='height_500hPa' #change this depending on which version of the code you're running
#height='height_250hPa'

#lat,lon, and masks
fname='/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform/history.2025-10-01_06.00.00.nc'
nc=Dataset(fname)
latCell=np.array(nc['latCell']) #in radians
lonCell=np.array(nc['lonCell']) #in radians
latCell*=(180/np.pi) #in degrees
lonCell*=(180/np.pi) #in degrees
Wcasp_mask=(latCell>=40)*(latCell<=50)*(lonCell>=30)*(lonCell<=50) #mask for NELC 1: 30-50 e, 40-50 n (W of caspian sea)
Wcasp_mask=(latCell>=45)*(latCell<=50)*(lonCell>=30)*(lonCell<=40)
#NOTE: the correlations for NELC 1 and NELC 2 are opposite in sign!
SEcasp_mask=(latCell>=30)*(latCell<=45)*(lonCell>=50)*(lonCell<=70) #mask for NELC 2: 50-70 e, 30-45 n (SE of caspian sea)
SEcasp_mask=(latCell>=30)*(latCell<=40)*(lonCell>=55)*(lonCell<=65)
Waus_mask=(latCell>=-25)*(latCell<=-15)*(lonCell>=90)*(lonCell<=110) #mask for Waus: 15-25 s, 90-110 e
Waus_mask=(latCell>=-25)*(latCell<=-15)*(lonCell>=95)*(lonCell<=105)

#script

Wcasp_averages=np.zeros((number_of_ensembles))
SEcasp_averages=np.zeros((number_of_ensembles))
Waus_averages=np.zeros((number_of_ensembles))

paths=np.empty(number_of_ensembles,dtype='object')

for i in range(number_of_ensembles):
    time_string=ftime.strftime('%Y-%m-%d_%H.%M.%S')
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.'+time_string+'.nc'
    fname=paths[i]
    nc=Dataset(fname)

    #finding the average of the response variables over a day
    sum_over_a_day=np.zeros((3))

    for k in range(8):
        hour_step=timedelta(hours=k*3)
        jtime=ftime.strptime(time_string,'%Y-%m-%d_%H.%M.%S')
        jtime+=hour_step
        ktime=jtime.strftime('%Y-%m-%d_%H.%M.%S')
        fname3=fname[:57]+ktime+'.nc'
        ne=Dataset(fname3)

        s1=(np.squeeze(np.array(ne[height])))
        
        sum_over_a_day[0]+=np.mean(s1[Wcasp_mask])
        sum_over_a_day[1]+=np.mean(s1[SEcasp_mask])
        sum_over_a_day[2]+=np.mean(s1[Waus_mask])

    sum_over_a_day/=8

    Wcasp_averages[i]=sum_over_a_day[0]
    SEcasp_averages[i]=sum_over_a_day[1]
    Waus_averages[i]=sum_over_a_day[2]

    print("currently on ensemble number " + str(i+1))

    nc.close()
    ne.close()

esa_dict={}
esa_dict['y1']=Wcasp_averages #(ensembles)
esa_dict['y2']=SEcasp_averages #(ensembles)
esa_dict['y3']=Waus_averages #(ensembles)
pickle.dump(esa_dict,open('/fs/scratch/PAS3252/yang/HONORS_THESIS/xs_and_ys_at_'+input_time[4:8]+'.pkl','wb'))