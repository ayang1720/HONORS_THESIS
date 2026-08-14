'''
NOTE: when making different plots, remember to change the name of the files created!
Creating a time series of MPAS plots NOTE: new update 8/3, trying new smois array
'''
#import packages
import pickle
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import glace_esa
from datetime import datetime, timedelta

#settings
path='/fs/scratch/PAS3252/yang/HONORS_THESIS/timeseries/'
number_of_ensembles=100
ncells=40962
response_variables=['height_500hPa','height_250hPa','rainnc','t2m','q2','']
response_variables=np.reshape(response_variables,(3,2))
response_variables2=['height_500hPa','height_250hPa','rainnc','t2m','q2','']
itime=datetime.strptime('20210714210000','%Y%m%d%H%M%S') #7/14/21 at 21Z
ftime=datetime.strptime('20210801210000','%Y%m%d%H%M%S') #8/1/21 at 21Z
time=itime

#script
esa_dict=pickle.load(open('ESA.pkl','rb'))
glaced_smois_array=esa_dict['smois_1d'] #100 smois values
#NOTE: edited to try the new smois 1d array instead
glaced_smois_array=esa_dict['smois_1d_new']
latCell=esa_dict['lat_1d'] #40962 lat cells
lonCell=esa_dict['lon_1d'] #40962 lon cells

while(time<=ftime):
    pvalues=glace_esa.glace_esa(glaced_smois_array,response_variables2[:-1],time,number_of_ensembles,ncells)
    fig,axs=plt.subplots(nrows=3,ncols=2,figsize=(24,15),constrained_layout=True,\
    subplot_kw={"projection":ccrs.PlateCarree()})
    for i in range(3):
        for j in range(2):
            ax=axs[i,j]
            gl=ax.gridlines(draw_labels=True,color='black',linewidth=1,linestyle=':')
            gl.xlabel_style={'fontsize':14};gl.ylabel_style={'fontsize':14}
            ax.add_feature(cfeature.LAND,edgecolor='black',linewidth=1.0,facecolor='none',zorder=100)
            if (i==2) and (j==1):
                ax.set_title('Blank Plot',fontsize=18)

    #edit on 7/27 to plot correlations instead of pvalues
    lons=lonCell
    lats=latCell
    number_of_statistically_significant=""

    #NOTE: chanegd limits vmin and vmax for correlations instead of pvalues

    #statgood,lons,lats=glace_esa.significant_cells(pvalues[0],lonCell,latCell)
    statgood=pvalues[0]
    #number_of_statistically_significant=round(statgood.size/ncells*100,3)
    axs[0,0].set_title('Smois on 7/14 vs. '+response_variables[0,0]+' on '+time.strftime('%m/%d')+', Significant: '+str(number_of_statistically_significant)+'%',fontsize=18)
    sc=axs[0,0].scatter(lons,lats,c=statgood,vmin=-.5,vmax=.5,\
                        cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())#-.5,.5 for correlations
    fig.colorbar(sc,ax=axs[0,0])
    #statgood,lons,lats=glace_esa.significant_cells(pvalues[1],lonCell,latCell)
    statgood=pvalues[1]
    #number_of_statistically_significant=round(statgood.size/ncells*100,3)
    axs[0,1].set_title('Smois on 7/14 vs. '+response_variables[0,1]+' on '+time.strftime('%m/%d')+', Significant: '+str(number_of_statistically_significant)+'%',fontsize=18)
    sc=axs[0,1].scatter(lons,lats,c=statgood,vmin=-.5,vmax=.5,\
                    cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree()) #0,0.05 for pvalues
    fig.colorbar(sc,ax=axs[0,1])
    #statgood,lons,lats=glace_esa.significant_cells(pvalues[2],lonCell,latCell)
    statgood=pvalues[2]
    #number_of_statistically_significant=round(statgood.size/ncells*100,3)
    axs[1,0].set_title('Smois on 7/14 vs. '+response_variables[1,0]+' on '+time.strftime('%m/%d')+', Significant: '+str(number_of_statistically_significant)+'%',fontsize=18)
    sc=axs[1,0].scatter(lons,lats,c=statgood,vmin=-.5,vmax=.5,\
                    cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
    fig.colorbar(sc,ax=axs[1,0])
    #statgood,lons,lats=glace_esa.significant_cells(pvalues[3],lonCell,latCell)
    statgood=pvalues[3]
    #number_of_statistically_significant=round(statgood.size/ncells*100,3)
    axs[1,1].set_title('Smois on 7/14 vs. '+response_variables[1,1]+' on '+time.strftime('%m/%d')+', Significant: '+str(number_of_statistically_significant)+'%',fontsize=18)
    sc=axs[1,1].scatter(lons,lats,c=statgood,vmin=-.5,vmax=.5,\
                    cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
    fig.colorbar(sc,ax=axs[1,1])
    #statgood,lons,lats=glace_esa.significant_cells(pvalues[4],lonCell,latCell)
    statgood=pvalues[4]
    #number_of_statistically_significant=round(statgood.size/ncells*100,3)
    axs[2,0].set_title('Smois on 7/14 vs. '+response_variables[2,0]+' on '+time.strftime('%m/%d')+', Significant: '+str(number_of_statistically_significant)+'%',fontsize=18)
    sc=axs[2,0].scatter(lons,lats,c=statgood,vmin=-.5,vmax=.5,\
                    cmap='jet',s=50,marker='*',transform=ccrs.PlateCarree())
    fig.colorbar(sc,ax=axs[2,0])
    fig.colorbar(sc,ax=axs[2,1])

    axs[0,0].set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    axs[0,1].set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    axs[1,0].set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    axs[1,1].set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    axs[2,0].set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    #plt.savefig(path+'glace_plot_'+time.strftime('%Y-%m-%d_%H.%M.%S')+'_pvalues_for_3_days.png') #changed name for correlations
    plt.savefig(path+'glace_plot_'+time.strftime('%Y-%m-%d_%H.%M.%S')+'_correlations_for_3_days_updated_with_cor_function_swapped.png') #changed name for correlations
    #plt.savefig(path+'glace_plot_'+time.strftime('%Y-%m-%d_%H.%M.%S')+'_correlations.png') #changed name for correlations
    #plt.savefig(path+'glace_plot_'+time.strftime('%Y-%m-%d_%H.%M.%S')+'.png') #changed name for correlations
    #plt.savefig(path+'glace_plot_'+time.strftime('%Y-%m-%d_%H.%M.%S')+'_correlations_for_3_days.png')

    plt.close()
    time+=timedelta(days=1)