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

#input: array of some meterological variable with shape (ensemble,ncells)
#output: an area averaged array of only the land over CONUS (ensemble)
def conusify(array_to_be_conusified):
    conusified_array=np.zeros(number_of_ensembles)
    for i in range(number_of_ensembles):
        array_to_be_conusified[i,:]*=land_mask*conus_mask
        conusified_array[i]=\
        array_to_be_conusified[i][array_to_be_conusified[i]!=0].mean()
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
    print(summation)
    covariance=summation/(number_of_ensembles-1)
    return covariance

#script

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

#constructing the atmospheric state array for the entire globe on aug 1
paths=np.empty(number_of_ensembles,dtype='object')
for i in range(number_of_ensembles):
    paths[i]='/fs/ess/PAS2635/LandAir_Predictability/member_'+str((i+1)).zfill(5)+'/diag.2021-08-01_21.00.00.nc'
z500_array=np.zeros((number_of_ensembles,ncells))
z250_array=np.zeros((number_of_ensembles,ncells))
rainc_array=np.zeros((number_of_ensembles,ncells))
t2m_array=np.zeros((number_of_ensembles,ncells))
q2_array=np.zeros((number_of_ensembles,ncells))
for i in range(number_of_ensembles):
    fname=paths[i]
    nc=Dataset(fname)
    z500_array[i,:]=(np.squeeze(np.array(nc['height_500hPa'])))
    #array of z500 heights, each with ncells=40962
    z250_array[i,:]=(np.squeeze(np.array(nc['height_250hPa'])))
    #array of z250 heights, each with ncells=40962
    rainc_array[i,:]=(np.squeeze(np.array(nc['rainc'])))
    #array of cumulative precipitation in mm, each with ncells=40962
    t2m_array[i,:]=(np.squeeze(np.array(nc['t2m'])))
    #array of 2 meter temperature K, each with ncells=40962
    q2_array[i,:]=(np.squeeze(np.array(nc['q2'])))
    #array of 2 meter specific humidity kg/kg, each with ncells=40962
    nc.close()

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

'''
print(np.sum(conus_mask)) #1158
print(np.sum(land_mask)) #11458
print(np.shape(conus_mask)) #(40962)
print(np.shape(smois_array)) #(100, 40962)
print(np.shape(land_mask)) #(40962)
'''

#conus_smois is an array of length 100 containing the area averaged smois
#for the contiguous united states; 1 for each of the 100 ensemble members
conus_smois=(conusify(smois_array))
conus_z500=(conusify(z500_array))
conus_z250=(conusify(z250_array))
conus_rainc=(conusify(rainc_array))
conus_t2m=(conusify(t2m_array))
conus_q2=(conusify(q2_array))

var_conus_smois=np.var(conus_smois,ddof=1)
var_conus_z500=np.var(conus_z500,ddof=1)
var_conus_z250=np.var(conus_z250,ddof=1)
var_conus_rainc=np.var(conus_rainc,ddof=1)
var_conus_t2m=np.var(conus_t2m,ddof=1)
var_conus_q2=np.var(conus_q2,ddof=1)

correlation_z500=cov(conus_smois,conus_z500)\
/np.sqrt(var_conus_smois)/np.sqrt(var_conus_z500)
correlation_z250=cov(conus_smois,conus_z250)\
/np.sqrt(var_conus_smois)/np.sqrt(var_conus_z250)
correlation_rainc=cov(conus_smois,conus_rainc)\
/np.sqrt(var_conus_smois)/np.sqrt(var_conus_rainc)
correlation_t2m=cov(conus_smois,conus_t2m)\
/np.sqrt(var_conus_smois)/np.sqrt(var_conus_t2m)
correlation_q2=cov(conus_smois,conus_q2)\
/np.sqrt(var_conus_smois)/np.sqrt(var_conus_q2)

print("The correlation between smois and z500: "+str(correlation_z500))
print("The correlation between smois and z250: "+str(correlation_z250))
print("The correlation between smois and rainc: "+str(correlation_rainc))
print("The correlation between smois and t2m: "+str(correlation_t2m))
print("The correlation between smois and q2: "+str(correlation_q2))