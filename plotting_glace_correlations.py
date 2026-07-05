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

#settings
ncells=40962
response_variables=['height_500hPa','height_250hPa','rainnc','t2m','q2','']
response_variables=np.reshape(response_variables,(3,2))

#script
esa_dict=pickle.load(open('ESA.pkl','rb'))
correlations=esa_dict['correlations']
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

sc=axs[0,0].scatter(lonCell,latCell,c=correlations[:,0],\
                    cmap='jet',vmin=0,vmax=1,s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[0,0])
axs[0,1].scatter(lonCell,latCell,c=correlations[:,1],\
                 cmap='jet',vmin=0,vmax=1,s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[0,1])
axs[1,0].scatter(lonCell,latCell,c=correlations[:,2],\
                 cmap='jet',vmin=0,vmax=1,s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[1,0])
axs[1,1].scatter(lonCell,latCell,c=correlations[:,3],\
                 cmap='jet',vmin=0,vmax=1,s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[1,1])
axs[2,0].scatter(lonCell,latCell,c=correlations[:,4],\
                 cmap='jet',vmin=0,vmax=1,s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[2,0])

number_of_statistically_significant=0
for i in range(ncells):
    if (correlations[i,2]>=0.2) and (correlations[i,2]<=100):
        number_of_statistically_significant+=1
number_of_statistically_significant/=ncells
print(number_of_statistically_significant)

plt.savefig('glace_plot_sample.png')
plt.close()