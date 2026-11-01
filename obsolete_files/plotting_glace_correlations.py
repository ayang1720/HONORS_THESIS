'''
note made 10/31/26:
pulls correlations from a pickle file and plots it
the pickle file it was supposed to pull from (esa.pkl) now no longer contains the correct array
changed it to pull from the ESA with 100 members pickle file, which produces a nonsensical plot
additionally, this creates a plot of the 5 atmospheric responses, which is now no longer the focus
'''

'''
MPAS Plots for GLACE smois
'''
#import packages
from netCDF4 import Dataset
import pickle
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import scipy.stats
import math

#settings
ncells=40962
response_variables=['height_500hPa','height_250hPa','rainnc','t2m','q2','']
response_variables=np.reshape(response_variables,(3,2))

#script
esa_dict=pickle.load(open('/users/PAS3252/ayang1720/HONORS_THESIS/pickle_files/ESA_with_100_members.pkl','rb'))
correlations=esa_dict['correlations']
smois_smois_correlation=esa_dict['smois_smois_correlation']
#print(np.shape(correlations)) #40962*5

#plotting
fname='/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform/history.2025-10-01_06.00.00.nc'
nc=Dataset(fname)
latCell=np.array(nc['latCell']) #in radians
lonCell=np.array(nc['lonCell']) #in radians
lonCell=np.rad2deg(lonCell)
latCell=np.rad2deg(latCell)
lon2d,lat2d=np.meshgrid(lonCell,latCell)

fig,axs=plt.subplots(nrows=3,ncols=2,figsize=(20,15),layout='constrained',\
subplot_kw={"projection":ccrs.PlateCarree()})
for i in range(3):
    for j in range(2):
        ax=axs[i,j]
        gl=ax.gridlines(draw_labels=True,color='black',linewidth=1,linestyle=':')
        gl.xlabel_style={'fontsize':14};gl.ylabel_style={'fontsize':14}
        ax.add_feature(cfeature.LAND,edgecolor='black',linewidth=1.0,facecolor='none',zorder=100)
        ax.set_title('Smois on 7/14 vs. '+response_variables[i,j]+' on 8/1',fontsize=24)
        if (i==2) and (j==1):
            ax.set_title('Smois on 7/14 vs. Smois on 7/14',fontsize=24)

sc=axs[0,0].scatter(lonCell,latCell,c=correlations[:,0],\
                    cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[0,0])
axs[0,1].scatter(lonCell,latCell,c=correlations[:,1],\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[0,1])
axs[1,0].scatter(lonCell,latCell,c=correlations[:,2],\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[1,0])
axs[1,1].scatter(lonCell,latCell,c=correlations[:,3],\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[1,1])
axs[2,0].scatter(lonCell,latCell,c=correlations[:,4],\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[2,0])
axs[2,1].scatter(lonCell,latCell,c=smois_smois_correlation,\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[2,1])

'''#print(np.shape(correlations))
ps=np.zeros((ncells,5))
for i in range(5):
        for k in range(40962):
            r=correlations[k,i]
            if (math.isinf(r)):
                r=0
            t_stat=r*np.sqrt((100-2)/(1-r**2))
            ps[k,j]=2*(1-scipy.stats.t.cdf(np.abs(t_stat),df=100-2))
new_ps=scipy.stats.false_discovery_control(ps,method='by')
number_of_statistically_significant=np.sum(new_ps<=0.05)
print(number_of_statistically_significant)'''

plt.savefig('/fs/scratch/PAS3252/yang/HONORS_THESIS/misc/glace_plot_sample.png')
plt.close()