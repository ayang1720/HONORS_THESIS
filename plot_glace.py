'''
MPAS Plots for GLACE smois
NEW 7/20/26
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
esa_dict=pickle.load(open('ESA.pkl','rb'))
correlations=esa_dict['correlations']
print(np.shape(correlations)) #40962*5

for i in range(40962):
    if(((correlations[2,i])>1) or (correlations[2,1]<-1) or math.isinf(correlations[2,1]) or np.isnan(correlations[2,1])):
        print("cell #"+str(i))

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

sc=axs[0,0].scatter(lonCell,latCell,c=correlations[0],\
                    cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[0,0])
sc=axs[0,1].scatter(lonCell,latCell,c=correlations[1],\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[0,1])
sc=axs[1,0].scatter(lonCell,latCell,c=correlations[2],\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[1,0])
sc=axs[1,1].scatter(lonCell,latCell,c=correlations[3],\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[1,1])
sc=axs[2,0].scatter(lonCell,latCell,c=correlations[4],\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[2,0])

'''
ps=np.zeros((ncells,5))
for i in range(5):
        for k in range(40962):
            r=correlations[i,k]
            if (math.isinf(r)):
                r=0
            t_stat=r*np.sqrt((100-2)/(1-r**2))
            ps[k,i]=2*(1-scipy.stats.t.cdf(np.abs(t_stat),df=100-2))
new_ps=scipy.stats.false_discovery_control(ps,method='by')
number_of_statistically_significant=np.zeros((5))
for i in range(5):
    number_of_statistically_significant[i]=np.sum(new_ps[:,i]<=0.05)
    print(number_of_statistically_significant)
'''

number_of_statistically_significant=np.zeros((5))
for i in range(5):
    number_of_statistically_significant[i]=np.sum(np.abs(correlations[i,:])>=0.2)
    print(number_of_statistically_significant[i]/40962)
statistically_significant_mask=np.abs(correlations)<0.2
correlations=np.where(statistically_significant_mask,correlations,0)

plt.savefig('glace_plot_7_20_26.png')

'''
sc=axs[0,0].scatter(lonCell,latCell,c=correlations[0],\
                    cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[0,0])
sc=axs[0,1].scatter(lonCell,latCell,c=correlations[1],\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[0,1])
sc=axs[1,0].scatter(lonCell,latCell,c=correlations[2],\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[1,0])
sc=axs[1,1].scatter(lonCell,latCell,c=correlations[3],\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[1,1])
sc=axs[2,0].scatter(lonCell,latCell,c=correlations[4],\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[2,0])
plt.savefig('glace_plot_7_20_26_only_significant.png')
'''
plt.close()