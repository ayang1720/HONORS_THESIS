'''
Alex Yang's Honors Thesis SU26 Task 4+
/users/PAS3252/ayang1720/HONORS_THESIS
/fs/ess/PAS2635/LandAir_Predictability
'''

#import packages
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import pickle

#settings
number_of_ensembles=100
number_of_iterations=729
itime=datetime.strptime('20210601000000','%Y%m%d%H%M%S') #6/1/21 at 0Z
ftime=datetime.strptime('20210831000000','%Y%m%d%H%M%S') #8/13/21 at 0Z
paths=np.empty(number_of_ensembles,dtype='object')
for i in range(number_of_ensembles):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/'

#NOTE: over here, we're going to initialize empty arrays for spatial averages
fname=paths[i]+'diag.'+itime.strftime('%Y-%m-%d_%H.%M.%S')+'.nc'
nc=Dataset(fname)
Z500=np.squeeze(np.array(nc['height_500hPa']))
Z250=np.squeeze(np.array(nc['height_250hPa']))
smois=np.squeeze(np.array(nc['smois']))
smois=smois[:,0]
Z500_average=np.zeros((len(Z500))) #(40962) instead of (40962,729,100)
Z250_average=np.zeros((len(Z250)))
smois_average=np.zeros((len(smois)))
t=itime

#redoing the loop to compute the spatially-averaged variance for each timestep in each iteration
iteration_number=1
while t<=ftime:
    Z500_array=np.zeros((len(Z500),number_of_ensembles))
    Z250_array=np.zeros((len(Z250),number_of_ensembles))
    smois_array=np.zeros((len(smois),number_of_ensembles))
    for i in range(number_of_ensembles):
        
        #accessing the netCDF4 file
        #height_500hPa(Time, nCells)
        #height_250hPa(Time, nCells)
        #smois(Time, nCells, nSoilLevels)
        fname=paths[i]+'diag.'+t.strftime('%Y-%m-%d_%H.%M.%S')+'.nc'
        nc=Dataset(fname)
        Z500=np.squeeze(np.array(nc['height_500hPa']))
        Z250=np.squeeze(np.array(nc['height_250hPa']))
        smois=np.squeeze(np.array(nc['smois']))
        smois=smois[:,0]
        
        Z500_array[:,i]=Z500
        Z250_array[:,i]=Z250
        smois_array[:,i]=smois

        nc.close()
    Z500_average[(iteration_number-1)]=np.mean(np.var(Z500_array,axis=1)) #var ens, averaged over space
    Z250_average[(iteration_number-1)]=np.mean(np.var(Z250_array,axis=1)) #var ens, averaged over space
    smois_average[(iteration_number-1)]=np.mean(np.var(smois_array,axis=1)) #var ens, averaged over space

    print('We are on time step '+str(iteration_number)+'/729.')
    iteration_number+=1
    t+=timedelta(hours=3)
    #end time loop

avg_dict={}
avg_dict['Z500']=Z500_average
avg_dict['Z250']=Z250_average
avg_dict['smois']=smois_average
pickle.dump(avg_dict,open('/fs/scratch/PAS3252/yang/HONORS_THESIS/averages_100.pkl','wb'))