'''
NOTE: when making different plots, remember to change the name of the files created!
Creating a time series of MPAS plots NOTE: new update 8/3, trying new smois array
'''

#NOTE: 10/4, heavily changed this script, look at past versions if want to return to them

#import packages
import pickle
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import glace_esa
from datetime import datetime
import sys

#settings
path='/fs/scratch/PAS3252/yang/HONORS_THESIS/timeseries/smois_6_cases/'
directories=['top_gp','bot_gp','top_ca','bot_ca','top_fl','bot_fl']
run_number=0

xs_dict=pickle.load(open('/fs/scratch/PAS3252/yang/HONORS_THESIS/xs_dict.pkl','rb'))
print(np.shape(xs_dict['x_t_gp']))
xs=np.array([(xs_dict['x_t_gp']),(xs_dict['x_b_gp']),(xs_dict['x_t_ca']),\
        (xs_dict['x_b_ca']),(xs_dict['x_t_fl']),(xs_dict['x_b_fl'])])

number_of_ensembles=500
ncells=40962
response_variables=['height_500hPa','height_250hPa','rainnc','t2m','q2','']
response_variables=np.reshape(response_variables,(3,2))
response_variables2=['height_500hPa','height_250hPa','rainnc','t2m','q2','']
itime=datetime.strptime('20210714210000','%Y%m%d%H%M%S') #7/14/21 at 21Z
ftime=datetime.strptime('20210801210000','%Y%m%d%H%M%S') #8/1/21 at 21Z
time=datetime.strptime(sys.argv[1],'%Y%m%d%H%M%S')

#script
esa_dict=pickle.load(open('ESA.pkl','rb'))
latCell=esa_dict['lat_1d'] #40962 lat cells
lonCell=esa_dict['lon_1d'] #40962 lon cells

fig,axs=plt.subplots(nrows=3,ncols=2,figsize=(24,15),constrained_layout=True,\
subplot_kw={"projection":ccrs.PlateCarree()})
for i in range(3):
    for j in range(2):
        ax=axs[i,j]
        gl=ax.gridlines(draw_labels=True,color='black',linewidth=1,linestyle=':')
        gl.xlabel_style={'fontsize':14};gl.ylabel_style={'fontsize':14}
        ax.add_feature(cfeature.LAND,edgecolor='black',linewidth=1.0,facecolor='none',zorder=100)

axs[0,0].set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
axs[0,1].set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
axs[1,0].set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
axs[1,1].set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
axs[2,0].set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
axs[2,1].set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

for x in xs:
    case=1
    print("currently on case: " + str(case))
    pvalues=glace_esa.glace_esa_pvalues(x,response_variables2[:-1],time,number_of_ensembles,ncells)

    #edit on 7/27 to plot correlations instead of pvalues
    lons=lonCell
    lats=latCell
    number_of_statistically_significant=""

    #NOTE: chanegd limits vmin and vmax for correlations instead of pvalues
    if case==1:
        statgood,lons,lats=glace_esa.significant_cells(pvalues[0],lonCell,latCell)
        axs[0,0].set_title('Smois on 7/14 vs. '+response_variables[0,0]+' on '+time.strftime('%m/%d')+', Case #'+str(case),fontsize=18)
        sc=axs[0,0].tricontourf(lons,lats,statgood,np.linspace(0,0.1,11),\
                        cmap='plasma',transform=ccrs.PlateCarree(),extend='both')
        fig.colorbar(sc,ax=axs[0,0])
    if case==2:
        statgood,lons,lats=glace_esa.significant_cells(pvalues[0],lonCell,latCell)
        axs[0,1].set_title('Smois on 7/14 vs. '+response_variables[0,0]+' on '+time.strftime('%m/%d')+', Case #'+str(case),fontsize=18)
        sc=axs[0,1].tricontourf(lons,lats,statgood,np.linspace(0,0.1,11),\
                        cmap='plasma',transform=ccrs.PlateCarree(),extend='both')
        fig.colorbar(sc,ax=axs[0,1])
    if case==3:
        statgood,lons,lats=glace_esa.significant_cells(pvalues[0],lonCell,latCell)
        axs[1,0].set_title('Smois on 7/14 vs. '+response_variables[0,0]+' on '+time.strftime('%m/%d')+', Case #'+str(case),fontsize=18)
        sc=axs[1,0].tricontourf(lons,lats,statgood,np.linspace(0,0.1,11),\
                        cmap='plasma',transform=ccrs.PlateCarree(),extend='both')
        fig.colorbar(sc,ax=axs[1,0])
    if case==4:
        statgood,lons,lats=glace_esa.significant_cells(pvalues[0],lonCell,latCell)
        axs[1,1].set_title('Smois on 7/14 vs. '+response_variables[0,0]+' on '+time.strftime('%m/%d')+', Case #'+str(case),fontsize=18)
        sc=axs[1,1].tricontourf(lons,lats,statgood,np.linspace(0,0.1,11),\
                        cmap='plasma',transform=ccrs.PlateCarree(),extend='both')
        fig.colorbar(sc,ax=axs[1,1])
    if case==5:
        statgood,lons,lats=glace_esa.significant_cells(pvalues[0],lonCell,latCell)
        axs[2,0].set_title('Smois on 7/14 vs. '+response_variables[0,0]+' on '+time.strftime('%m/%d')+', Case #'+str(case),fontsize=18)
        sc=axs[2,0].tricontourf(lons,lats,statgood,np.linspace(0,0.1,11),\
                        cmap='plasma',transform=ccrs.PlateCarree(),extend='both')
        fig.colorbar(sc,ax=axs[2,0])
    if case==6:
        statgood,lons,lats=glace_esa.significant_cells(pvalues[0],lonCell,latCell)
        axs[2,1].set_title('Smois on 7/14 vs. '+response_variables[0,0]+' on '+time.strftime('%m/%d')+', Case #'+str(case),fontsize=18)
        sc=axs[2,1].tricontourf(lons,lats,statgood,np.linspace(0,0.1,11),\
                        cmap='plasma',transform=ccrs.PlateCarree(),extend='both')
        fig.colorbar(sc,ax=axs[2,1])

    case+=1
plt.savefig(path+time.strftime('%Y-%m-%d_%H.%M.%S')+'.png')
plt.close()