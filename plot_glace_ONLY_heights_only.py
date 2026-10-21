#this version predates 10/19, computes correlation for all 6 smois cases (the new version only looks at the first case)

'''
NOTE: when making different plots, remember to change the name of the files created!
Creating a time series of MPAS plots NOTE: new update 8/3, trying new smois array
'''

#10/9 gif creation:
#convert -delay 50 -loop 0 /fs/scratch/PAS3252/yang/HONORS_THESIS/timeseries/smois_6_cases/*.png /fs/scratch/PAS3252/yang/HONORS_THESIS/timeseries/smois_6_cases/6_cases_time_series_updated.gif

#NOTE: 10/4, heavily changed this script, look at past versions if want to return to them

#import packages
import pickle
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import glace_esa_heights_only
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
response_variables=['height_500hPa']
response_variables2=['height_500hPa']
itime=datetime.strptime('20210714210000','%Y%m%d%H%M%S') #7/14/21 at 21Z
ftime=datetime.strptime('20210801210000','%Y%m%d%H%M%S') #8/1/21 at 21Z
time=datetime.strptime(sys.argv[1],'%Y%m%d%H%M%S')

#script
esa_dict=pickle.load(open('ESA.pkl','rb'))
latCell=esa_dict['lat_1d'] #40962 lat cells
lonCell=esa_dict['lon_1d'] #40962 lon cells

fname='/fs/scratch/PAS3252/yang/HONORS_THESIS/statgood_pickle_files/global_arrays/global_array_at_time_'+time.strftime('%m-%d')+'_using_timedelta_24_hours.pkl'
global_arrays_dict=pickle.load(open('fname','rb'))

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

case=1
lons,lats=lonCell,latCell
statgood_dict={}

for x in xs:
    print("currently on case: " + str(case))
    #NOTE: 10/8 edit response_variables now only contains heights
    pvalues=glace_esa_heights_only.glace_esa_pvalues(x,global_arrays_dict['global_array'],number_of_ensembles,ncells)

    statgood=pvalues #pvalues will be different depending on how many times this loop has been run
    #pvalues is a 1d array containing 40962 pvalues (one for each cell in MPAS)

    if case==1:
        axs[0,0].set_title('Smois on 7/14 vs. '+response_variables[0]+' on '+time.strftime('%m/%d')+', Case #'+str(case),fontsize=18)
        sc=axs[0,0].tricontourf(lons,lats,statgood,np.linspace(0,0.1,11),\
                        cmap='plasma',transform=ccrs.PlateCarree(),extend='both')
        fig.colorbar(sc,ax=axs[0,0])
        statgood_dict['case_1']=statgood
    if case==2:
        axs[0,1].set_title('Smois on 7/14 vs. '+response_variables[0]+' on '+time.strftime('%m/%d')+', Case #'+str(case),fontsize=18)
        sc=axs[0,1].tricontourf(lons,lats,statgood,np.linspace(0,0.1,11),\
                        cmap='plasma',transform=ccrs.PlateCarree(),extend='both')
        fig.colorbar(sc,ax=axs[0,1])
        statgood_dict['case_2']=statgood
    if case==3:
        axs[1,0].set_title('Smois on 7/14 vs. '+response_variables[0]+' on '+time.strftime('%m/%d')+', Case #'+str(case),fontsize=18)
        sc=axs[1,0].tricontourf(lons,lats,statgood,np.linspace(0,0.1,11),\
                        cmap='plasma',transform=ccrs.PlateCarree(),extend='both')
        fig.colorbar(sc,ax=axs[1,0])
        statgood_dict['case_3']=statgood
    if case==4:
        axs[1,1].set_title('Smois on 7/14 vs. '+response_variables[0]+' on '+time.strftime('%m/%d')+', Case #'+str(case),fontsize=18)
        sc=axs[1,1].tricontourf(lons,lats,statgood,np.linspace(0,0.1,11),\
                        cmap='plasma',transform=ccrs.PlateCarree(),extend='both')
        fig.colorbar(sc,ax=axs[1,1])
        statgood_dict['case_4']=statgood
    if case==5:
        axs[2,0].set_title('Smois on 7/14 vs. '+response_variables[0]+' on '+time.strftime('%m/%d')+', Case #'+str(case),fontsize=18)
        sc=axs[2,0].tricontourf(lons,lats,statgood,np.linspace(0,0.1,11),\
                        cmap='plasma',transform=ccrs.PlateCarree(),extend='both')
        fig.colorbar(sc,ax=axs[2,0])
        statgood_dict['case_5']=statgood
    if case==6:
        axs[2,1].set_title('Smois on 7/14 vs. '+response_variables[0]+' on '+time.strftime('%m/%d')+', Case #'+str(case),fontsize=18)
        sc=axs[2,1].tricontourf(lons,lats,statgood,np.linspace(0,0.1,11),\
                        cmap='plasma',transform=ccrs.PlateCarree(),extend='both')
        fig.colorbar(sc,ax=axs[2,1])
        statgood_dict['case_6']=statgood
    case+=1
plt.savefig(path+time.strftime('%Y-%m-%d')+'_created_with_pickle_files_avg_over_24_hours.png')
plt.close()

#pickle.dump(statgood_dict,open('/fs/scratch/PAS3252/yang/HONORS_THESIS/statgood_pickle_files/'+time.strftime('%m-%d')+'.pkl','wb'))

#NOTE: this is the pickle file not the png
#NOTE: updated 10/21 png name, note that these are created with pickle files
#pickle.dump(statgood_dict,open('/fs/scratch/PAS3252/yang/HONORS_THESIS/statgood_pickle_files/'+time.strftime('%m-%d')+'_created_with_pickle_files_avg_over_24_hours.pkl','wb'))