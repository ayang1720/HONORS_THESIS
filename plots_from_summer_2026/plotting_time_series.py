'''
Alex Yang's Honors Thesis SU26 Task 3 Part 2
'''

#edited 6/16/26 to use /fs/scratch/PAS3252/yang/HONORS_THESIS/averages_100.pkl

#import packages
import pickle
import matplotlib.pyplot as plt
import numpy as np

#grabbing data
avg_dict=pickle.load(open('/fs/scratch/PAS3252/yang/HONORS_THESIS/averages_100.pkl','rb'))
#accidentally created this array to hold 40962 values instead of 729
Z500_average=avg_dict['Z500']
Z250_average=avg_dict['Z250']
smois_average=avg_dict['smois']
Z500_time=np.arange(len(Z500_average))
Z250_time=np.arange(len(Z250_average))
smois_time=np.arange(len(smois_average))
Z500_time+=1
Z250_time+=1
smois_time+=1

#plotting
fig,axs=plt.subplots(nrows=3,ncols=1,constrained_layout=True,figsize=(10,16))
axs[0].plot(Z500_time,Z500_average)
axs[1].plot(Z250_time,Z250_average)
axs[2].plot(smois_time,smois_average)
axs[0].set_xlabel('Time',fontsize=24)
axs[0].set_ylabel('Variance',fontsize=24)
axs[1].set_xlabel('Time',fontsize=24)
axs[1].set_ylabel('Variance',fontsize=24)
axs[2].set_xlabel('Time',fontsize=24)
axs[2].set_ylabel('Variance',fontsize=24)
axs[0].set_title('Height at 500mb',fontsize=30)
axs[1].set_title('Height at 250mb',fontsize=30)
axs[2].set_title('Soil Moisture',fontsize=30)
plt.savefig('time_series_plots_100_edited.png')
plt.close()

#editing the pickle file