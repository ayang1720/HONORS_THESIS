'''
Alex Yang's Honors Thesis SU26 Task 5+
'''

import pickle
import matplotlib.pyplot as plt
import numpy as np

esa_dict=pickle.load(open('ESA.pkl','rb'))
smois=esa_dict['smois']
z500=esa_dict['z500']
z250=esa_dict['z250']
rainc=esa_dict['rainc']
t2m=esa_dict['t2m']
q2=esa_dict['q2']
fig,axs=plt.subplots(nrows=3,ncols=2,constrained_layout=True,figsize=(10,16))
axs[0,0].scatter(smois,z500)
axs[0,1].scatter(smois,z250)
axs[1,0].scatter(smois,rainc)
axs[1,1].scatter(smois,t2m)
axs[2,0].scatter(smois,q2)
for i in range(3):
    for j in range(2):
        axs[i,j].set_xlabel('smois',fontsize=24)
axs[0,0].set_ylabel('z500',fontsize=24)
axs[0,1].set_ylabel('z250',fontsize=24)
axs[1,0].set_ylabel('rainc',fontsize=24)
axs[1,1].set_ylabel('t2m',fontsize=24)
axs[2,0].set_ylabel('q2',fontsize=24)
plt.savefig('esa_scatter_plots.png')
plt.close()