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
response_variables2=['height_500hPa','height_250hPa','rainnc','t2m','q2','']

#script
esa_dict=pickle.load(open('ESA.pkl','rb'))
correlations=esa_dict['p-values'] #edited to be p-values
flag_keep=esa_dict['flag_keep'] #NOTE: have not applied the flag yet to correlations! nan still replaced with 0!

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

pvalue_mask=correlations[0]<=0.5
lonCell1=lonCell[pvalue_mask]
latCell1=latCell[pvalue_mask]
correlations0=correlations[0][pvalue_mask]
sc=axs[0,0].scatter(lonCell1,latCell1,c=correlations0,\
                    cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[0,0])
pvalue_mask=correlations[1]<=0.5
lonCell2=lonCell[pvalue_mask]
latCell2=latCell[pvalue_mask]
correlations1=correlations[1][pvalue_mask]
sc=axs[0,1].scatter(lonCell2,latCell2,c=correlations1,\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[0,1])
lonCell3=lonCell[flag_keep]
latCell3=latCell[flag_keep]
correlations2=correlations[2][flag_keep]
pvalue_mask=correlations[2]<=0.5
lonCell3=lonCell[pvalue_mask]
latCell3=latCell[pvalue_mask]
correlations2=correlations[2][pvalue_mask]
sc=axs[1,0].scatter(lonCell3,latCell3,c=correlations2,\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[1,0])
pvalue_mask=correlations[3]<=0.5
lonCell4=lonCell[pvalue_mask]
latCell4=latCell[pvalue_mask]
correlations3=correlations[3][pvalue_mask]
sc=axs[1,1].scatter(lonCell4,latCell4,c=correlations3,\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[1,1])
pvalue_mask=correlations[4]<=0.5
lonCell5=lonCell[pvalue_mask]
latCell5=latCell[pvalue_mask]
correlations4=correlations[4][pvalue_mask]
sc=axs[2,0].scatter(lonCell5,latCell5,c=correlations4,\
                 cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
fig.colorbar(sc,ax=axs[2,0])

number_of_statistically_significant=np.zeros((5))
for i in range(5):
    number_of_statistically_significant[i]=np.sum(correlations[i,:]<=0.05)
    print("for variable " + response_variables2[i] + ", there are this many significant cells:")
    print(number_of_statistically_significant[i])
    print(number_of_statistically_significant[i]/40962)

plt.savefig('glace_plot_7_20_26.png')

plt.close()