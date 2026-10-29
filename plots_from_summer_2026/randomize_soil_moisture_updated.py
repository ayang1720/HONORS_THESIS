'''
task 2 honors project, randomizing soil moisture Alex Yang 5/22/26
/fs/scratch/PAS3252/yang/HONORS_THESIS
/users/PAS3252/ayang1720/HONORS_THESIS
/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform/
ncdump -h x1.40962.init.nc
'''
#NOTE: the files produced will be saved to the scratch directory
'''
float smois(Time, nCells, nSoilLevels) ;
        smois:units = "m3 m^{-3}" ;
        smois:long_name = "soil moisture" ;
float xland(Time, nCells) ;
        xland:units = "unitless" ;
        xland:long_name = "land-ocean mask (1=land including sea-ice ; 2=ocean)" ;
remember that with MPAS coordinates, they're flattened into 1d arrays; keep nonland unchanged
'''

#import packages
from netCDF4 import Dataset
import shutil #from chatgpt
import numpy as np

#script
fname='x1.40962.init.nc'
new_files=np.empty(100,dtype=object) #help from gemini
for i in range(100): #creates 100 file names
    new_files[i]='/fs/scratch/PAS3252/yang/HONORS_THESIS/'+fname[:8]+'_copy'+str(i)+fname[8:]
    #print(new_files[i])

#turns these variables from the netCDF file into arrays
nc=Dataset(fname)
soil_moisture=np.squeeze(np.array(nc['smois']))
soil_layer_1=soil_moisture[:,0]
soil_layer_2=soil_moisture[:,1]
soil_layer_3=soil_moisture[:,2]
soil_layer_4=soil_moisture[:,3]



'''
Fix 1 use the landmask variable instead to allow for the boolean logic mask
'''


#land_mask=np.squeeze(np.array(nc['xland']))
land_mask=np.squeeze(np.array(nc['landmask']))
sea_mask=np.logical_not(land_mask) #from chatgpt, turns it into an inverted mask



'''
print(soil_moisture.shape)
print(land_mask.shape)
print(soil_layer_1.shape)
print(soil_layer_2.shape)
print(soil_layer_3.shape)
print(soil_layer_4.shape)
'''
#iterates over every new file, and every lat/lon cell in MPAS (done with array broadcasting)
for i in range(100):
    shutil.copy(fname,new_files[i])
    nd=Dataset(new_files[i],'r+') #from chatgpt, note that nd not nc
    soil_moisture_edit=nd.variables['smois']
    #not squeezed because we're returning it as it was

    #noise will be 0 if over water
    noise=np.random.normal(0,0.01,len(land_mask))
    #noise*=sea_mask #sea mask is 0 if not land
    
    #============Fix 3 changing from sea_mask to land_mask
    
    #noise*=sea_mask
    noise*=land_mask#=================================If inverted sea mask only sea gets noise
    
    
    
    #NOTE: should each soil layer get the same noise treatment?
    soil_layer_1_edit=soil_layer_1*(noise+1)
    soil_layer_2_edit=soil_layer_2*(noise+1)
    soil_layer_3_edit=soil_layer_3*(noise+1)
    soil_layer_4_edit=soil_layer_4*(noise+1)

# =================Fix 2 apply the perturbation to the new edited varaible
    soil_moisture_edit[0,:,0] = soil_layer_1_edit
    soil_moisture_edit[0,:,1] = soil_layer_2_edit
    soil_moisture_edit[0,:,2] = soil_layer_3_edit
    soil_moisture_edit[0,:,3] = soil_layer_4_edit
    
#     soil_moisture_edit[0,:,0]=soil_layer_1
#     soil_moisture_edit[0,:,1]=soil_layer_2
#     soil_moisture_edit[0,:,2]=soil_layer_3
#     soil_moisture_edit[0,:,3]=soil_layer_4
    

    print("we are on iteration #"+str(i))

    nd.close()


    '''
    Hi Alex,
    After looking over the code I think I noticed a few bugs/errors.
   1. The biggest issue is that xland stores its data as 1 for land and 2 for ocean.  Therefore when you apply the np.logical_not this makes everything false.  1 fix is using the landmask variable
    which stores its data as 1 for land or 0 for ocean which should fix that main issue.---------Can just get rid of entirely
   2. The perturbation you apply to each soil layer is not actually being rewritten as you compute it but then use the previous soil_layer
   3. Actually the inverted sea mask isnt really needed as when you apply that to noise*=sea_mask it only applies the noise to ocean
    
    '''