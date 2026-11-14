'''
10/31 edit:
changed to only look at the height_500hPA response

11/3 edit:
this version of the script imports the pickle files including the global arrays
for each day instead of constructing them for all 100 members, thus saving time

11/10 edit:
made it an argument to call the function to plot to consider whether we're plotting
correlations or pvalues AND function call includes whether it's avg 24 days or hours

11/13 edit:
now accepts 3 settings, creating 2**3=8 possbile plots depending on combination of booleans
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
import numpy as np

#settings

#NOTE: this setting decides whether the plots being created are correlations or pvalues
#it'll change vmin, vmax, the functions being run in old_glace_esa.py, and pvalues
correlations=False
correlations=True

#NOTE: change these settings to switch which of the 2 pickle files is used for plotting
#uses the pickle file for 24 days average (old, wrong)
#averaging_over_24_days=True
#usees the pickle file for 24 hours average (new, fixed)
averaging_over_24_days=False

#NOTE: change these settings to decide whether to plot using 850 mb or 500 mb geometric heights
#height_500_hPa=True
height_500_hPa=False

#NOTE: This value should be changed to accommodate >100 members!
number_of_ensembles=100

path='/fs/scratch/PAS3252/yang/HONORS_THESIS/timeseries/old_plots/'
ncells=40962
response_variables=['height_500hPa']
itime=datetime.strptime('20210714210000','%Y%m%d%H%M%S') #7/14/21 at 21Z
ftime=datetime.strptime('20210801210000','%Y%m%d%H%M%S') #8/1/21 at 21Z
time=itime

#load the global arrays
#NOTE: this version of the pickle file includes only 100 ensemble members, but with all 40962 cells
#each array has size (number_of_ensembles,ncells) OR (100,40962)

times=np.array(["20210714000000", "20210715000000", "20210716000000", "20210717000000", "20210718000000",\
                "20210719000000", "20210720000000", "20210721000000", "20210722000000", "20210723000000",\
                "20210724000000", "20210725000000", "20210726000000", "20210727000000", "20210728000000",\
                "20210729000000", "20210730000000", "20210731000000", "20210801000000"])
i=0
datetimes=np.empty(len(times),dtype=object)

for time in times:
    datetimes[i]=datetime.strptime(time,'%Y%m%d%H%M%S')
    i+=1

global_arrays=np.empty(len(times),dtype=object)

#plotting function

def plot_glace(correlations_boolean,over_24_days_boolean,height_500_hPa,number_of_ensembles):
    if correlations_boolean==True:
        plotting_over_string='_correlations'
    if correlations_boolean==False:
        plotting_over_string='_pvalues'

    if over_24_days_boolean==True:
        period='days'
    if over_24_days_boolean==False:
        period='hours'
    if height_500_hPa==True:
        pressure_level="500 hPa heights"
    if height_500_hPa==False:
        pressure_level="850 hPa heights"

    for i in range(19):
        time=datetimes[i]

        if height_500_hPa==True:
            fname='/fs/scratch/PAS3252/yang/HONORS_THESIS/statgood_pickle_files/850hPa_global_arrays/850_hPa_global_array_at_time_'+time.strftime('%m-%d')+'_using_timedelta_24_'+period+'.pkl'
        if height_500_hPa==False:
            fname='/fs/scratch/PAS3252/yang/HONORS_THESIS/statgood_pickle_files/global_arrays/global_array_at_time_'+time.strftime('%m-%d')+'_using_timedelta_24_'+period+'.pkl'

        global_arrays_dict=pickle.load(open(fname,'rb'))
        this_global_array=global_arrays_dict['global_array']
        #NOTE: need to consider the first 100 members only, slice the first number_of_ensembles
        #each global_array has size (number_of_ensembles,ncells)
        global_arrays[i]=this_global_array[:number_of_ensembles,:]
        print("inserting global array " + str( i + 1 ) + " of 19")

    #global_arrays contains 19 elements---one element for each day of the simulation
    #that way, the while loop below can iterate through all 19 days

    #script
    esa_dict=pickle.load(open('old_ESA.pkl','rb'))
    glaced_smois_array=esa_dict['smois_1d'] #100 smois values
    latCell=esa_dict['lat_1d'] #40962 lat cells
    lonCell=esa_dict['lon_1d'] #40962 lon cells

    i=0

    time=itime #because an earlier line made time=ftime after running the loop

    while(time<=ftime):

        if correlations_boolean==True:
            pvalues=old_glace_esa.glace_esa_input_global_arrays_correlations(glaced_smois_array,global_arrays[i],number_of_ensembles,ncells)
        if correlations_boolean==False:
            pvalues=old_glace_esa.glace_esa_input_global_arrays(glaced_smois_array,global_arrays[i],number_of_ensembles,ncells)

        fig,axs=plt.subplots(nrows=1,ncols=1,figsize=(15,10),layout='constrained',\
        subplot_kw={"projection":ccrs.PlateCarree()})
        gl=axs.gridlines(draw_labels=True,color='black',linewidth=1,linestyle=':')
        gl.xlabel_style={'fontsize':14};gl.ylabel_style={'fontsize':14}
        axs.add_feature(cfeature.LAND,edgecolor='black',linewidth=1.0,facecolor='none',zorder=100)
        axs.set_xlim(-180, 180)
        axs.set_ylim(-90, 90)

        if correlations_boolean==False:
            statgood,lons,lats=old_glace_esa.significant_cells(pvalues[0],lonCell,latCell)
            vminvalue=0
            vmaxvalue=0.05
        if correlations_boolean==True:
            statgood,lons,lats=pvalues,lonCell,latCell
            vminvalue=-.5
            vmaxvalue=.5

        axs.set_title('Smois on 7/14 vs. '+pressure_level+' on '+time.strftime('%m/%d')+plotting_over_string,fontsize=18)
        
        sc=axs.scatter(lons,lats,c=statgood,\
                            cmap='jet',vmin=vminvalue,vmax=vmaxvalue,s=50,marker='*',transform=ccrs.PlateCarree())
        fig.colorbar(sc,shrink=0.7,orientation="horizontal")

        if height_500_hPa==True:
            plt.savefig(path+'glace_plot_'+time.strftime('%Y-%m-%d_%H.%M.%S')+'_using_pickle_file_averaged_over_24_'+period+plotting_over_string+'.png')
        if height_500_hPa==False:
            plt.savefig(path+'/850_hPa/glace_plot_'+time.strftime('%Y-%m-%d_%H.%M.%S')+'_using_pickle_file_averaged_over_24_'+period+plotting_over_string+'.png')

        plt.close()
        time+=timedelta(days=1)
        print("plotting day " + str ( i + 1 ) + " of 19")
        i+=1


#script

correlations,averaging_over_24_days,height_500_hPa=True,True,True
print("Plotting using correlations == " + str(correlations))
print("Plotting using height 500 hPa == " + str(height_500_hPa))
print("Plotting using number of ensembles == " + str(number_of_ensembles))
print("Plotting using averaging over 24 days == " + str(averaging_over_24_days))
plot_glace(correlations,averaging_over_24_days,height_500_hPa,number_of_ensembles)
print("case 1")

correlations,averaging_over_24_days,height_500_hPa=True,True,False
print("Plotting using correlations == " + str(correlations))
print("Plotting using height 500 hPa == " + str(height_500_hPa))
print("Plotting using number of ensembles == " + str(number_of_ensembles))
print("Plotting using averaging over 24 days == " + str(averaging_over_24_days))
plot_glace(correlations,averaging_over_24_days,height_500_hPa,number_of_ensembles)
print("case 2")

correlations,averaging_over_24_days,height_500_hPa=True,False,True
print("Plotting using correlations == " + str(correlations))
print("Plotting using height 500 hPa == " + str(height_500_hPa))
print("Plotting using number of ensembles == " + str(number_of_ensembles))
print("Plotting using averaging over 24 days == " + str(averaging_over_24_days))
plot_glace(correlations,averaging_over_24_days,height_500_hPa,number_of_ensembles)
print("case 3")

correlations,averaging_over_24_days,height_500_hPa=False,True,True
print("Plotting using correlations == " + str(correlations))
print("Plotting using height 500 hPa == " + str(height_500_hPa))
print("Plotting using number of ensembles == " + str(number_of_ensembles))
print("Plotting using averaging over 24 days == " + str(averaging_over_24_days))
plot_glace(correlations,averaging_over_24_days,height_500_hPa,number_of_ensembles)
print("case 4")

correlations,averaging_over_24_days,height_500_hPa=True,False,False
print("Plotting using correlations == " + str(correlations))
print("Plotting using height 500 hPa == " + str(height_500_hPa))
print("Plotting using number of ensembles == " + str(number_of_ensembles))
print("Plotting using averaging over 24 days == " + str(averaging_over_24_days))
plot_glace(correlations,averaging_over_24_days,height_500_hPa,number_of_ensembles)
print("case 5")

correlations,averaging_over_24_days,height_500_hPa=False,True,False
print("Plotting using correlations == " + str(correlations))
print("Plotting using height 500 hPa == " + str(height_500_hPa))
print("Plotting using number of ensembles == " + str(number_of_ensembles))
print("Plotting using averaging over 24 days == " + str(averaging_over_24_days))
plot_glace(correlations,averaging_over_24_days,height_500_hPa,number_of_ensembles)
print("case 6")

correlations,averaging_over_24_days,height_500_hPa=False,False,True
print("Plotting using correlations == " + str(correlations))
print("Plotting using height 500 hPa == " + str(height_500_hPa))
print("Plotting using number of ensembles == " + str(number_of_ensembles))
print("Plotting using averaging over 24 days == " + str(averaging_over_24_days))
plot_glace(correlations,averaging_over_24_days,height_500_hPa,number_of_ensembles)
print("case 7")

correlations,averaging_over_24_days,height_500_hPa=False,False,False
print("Plotting using correlations == " + str(correlations))
print("Plotting using height 500 hPa == " + str(height_500_hPa))
print("Plotting using number of ensembles == " + str(number_of_ensembles))
print("Plotting using averaging over 24 days == " + str(averaging_over_24_days))
plot_glace(correlations,averaging_over_24_days,height_500_hPa,number_of_ensembles)
print("case 8")