'''
MPAS Plots for Alex Yang's Honors Thesis SU26 5/20/26
/users/PAS3252/ayang1720/HONORS_THESIS
/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform
'''
#import packages
from netCDF4 import Dataset
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

#script
fname='/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform/diag.2025-10-01_06.00.00.nc'
nc=Dataset(fname)
#geometric height interpolated to 500 hPa in units of meters (Time, nCells)
height_500hPa=np.squeeze(np.array(nc['height_500hPa'])) #nCells==40962
print(np.shape(height_500hPa))
fname='/fs/ess/PAS2635/Generalized_Predictability_MPAS/120km_uniform/history.2025-10-01_06.00.00.nc'
nc=Dataset(fname)
latCell=np.array(nc['latCell']) #in radians
lonCell=np.array(nc['lonCell']) #in radians
smois=np.squeeze(np.array(nc['smois'])) #soil moisture in meters cubed per meters cubed (Time, nCells, nSoilLevels)
print(np.shape(latCell))
print(np.shape(lonCell))
print(np.shape(smois)) #(40962, 4)

#plotting
lon2d,lat2d=np.meshgrid(lonCell,latCell)
fig,axs=plt.subplots(nrows=1,ncols=1,figsize=(10,7),layout='constrained',\
subplot_kw={"projection":ccrs.PlateCarree()})
gl=axs.gridlines(draw_labels=True,color='black',linewidth=1,linestyle=':')
gl.xlabel_style={'fontsize':14};gl.ylabel_style={'fontsize':14}
axs.add_feature(cfeature.LAND,edgecolor='black',linewidth=1.0,facecolor='none',zorder=100)
axs.set_title('October 1st at 6-UTC Contour Plot',fontsize=24)
#scatter plot
lonCell=np.rad2deg(lonCell)
latCell=np.rad2deg(latCell)
cs=axs.tricontourf(lonCell,latCell,smois[:,0],cmap='jet',transform=ccrs.PlateCarree())
axs.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree()) #got help from chatgpt
cbar=plt.colorbar(cs,orientation='horizontal')
cbar.set_label('soil moisture layer 1',fontsize=24)

#tricontourf plots (i have never used this module before, so i got help from chatGPT)

plt.savefig('plot_scatter_soil1_contour.png')
plt.close()