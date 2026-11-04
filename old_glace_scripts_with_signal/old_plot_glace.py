'''
10/31 edit:
changed to only look at the height_500hPA response
'''

'''
Creating a time series of MPAS plots
'''
#import packages
import pickle
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import old_glace_esa
from datetime import datetime, timedelta

#settings
path='/fs/scratch/PAS3252/yang/HONORS_THESIS/timeseries/old_plots/'
number_of_ensembles=100
ncells=40962
response_variables=['height_500hPa']
itime=datetime.strptime('20210714210000','%Y%m%d%H%M%S') #7/14/21 at 21Z
ftime=datetime.strptime('20210801210000','%Y%m%d%H%M%S') #8/1/21 at 21Z
time=itime

#script
esa_dict=pickle.load(open('old_ESA.pkl','rb'))
glaced_smois_array=esa_dict['smois_1d'] #100 smois values
latCell=esa_dict['lat_1d'] #40962 lat cells
lonCell=esa_dict['lon_1d'] #40962 lon cells

while(time<=ftime):
    pvalues=old_glace_esa.glace_esa(glaced_smois_array,response_variables,time,number_of_ensembles,ncells)
    fig,axs=plt.subplots(nrows=1,ncols=1,figsize=(15,10),layout='constrained',\
    subplot_kw={"projection":ccrs.PlateCarree()})
    gl=axs.gridlines(draw_labels=True,color='black',linewidth=1,linestyle=':')
    gl.xlabel_style={'fontsize':14};gl.ylabel_style={'fontsize':14}
    axs.add_feature(cfeature.LAND,edgecolor='black',linewidth=1.0,facecolor='none',zorder=100)
    axs.set_xlim(-180, 180)
    axs.set_ylim(-90, 90)

    statgood,lons,lats=old_glace_esa.significant_cells(pvalues[0],lonCell,latCell)
    number_of_statistically_significant=round(statgood.size/ncells*100,3)
    axs.set_title('Smois on 7/14 vs. '+response_variables[0]+' on '+time.strftime('%m/%d'),fontsize=18)
    sc=axs.scatter(lons,lats,c=statgood,\
                        cmap='jet',vmin=0,vmax=0.05,s=50,marker='*',transform=ccrs.PlateCarree())
    fig.colorbar(sc,shrink=0.7,orientation="horizontal")

    plt.savefig(path+'glace_plot_'+time.strftime('%Y-%m-%d_%H.%M.%S')+'.png')

    plt.close()
    time+=timedelta(days=1)