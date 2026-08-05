'''
Alex Yang's Honors Thesis
next: check the dictionary, if masks were applied to the right arrays, plot over the globe
'''

import pickle
import numpy as np

number_of_ensembles=500
number_of_days=19

esa_dict=pickle.load(open('/fs/scratch/PAS3252/yang/HONORS_THESIS/xs_and_ys.pkl','rb'))
averaged_smois_array=esa_dict['x1'] #(ensembles)
averaged_russia_array=esa_dict['x2'] #(ensembles)
Wcasp_averages=esa_dict['y1'] #(days,ensembles)
SEcasp_averages=esa_dict['y2'] #(days,ensembles)
Waus_averages=esa_dict['y3'] #(days,ensembles)
x1=esa_dict['x1'] #(ensembles)
x2=esa_dict['x2'] #(ensembles)
y1=esa_dict['y1'] #(days,ensembles)
y2=esa_dict['y2'] #(days,ensembles)
y3=esa_dict['y3'] #(days,ensembles)

'''
print(np.shape(y1))
print(np.shape(y2))
print(np.shape(y3))
'''

esa1_dict=pickle.load(open('ESA.pkl','rb'))
latCell=esa1_dict['lat_1d'] #40962 lat cells
lonCell=esa1_dict['lon_1d'] #40962 lon cells

#returns correlation coefficient for a 1d array and a 2d array
def cor2d(X1d,R2d,number_of_ensembles):
    number_of_ensembles = len(X1d)
    X1d_mean=np.mean(X1d)
    R2d_mean=np.mean(R2d, axis=-1)
    X1d_pert = X1d - X1d_mean
    R2d_pert = (R2d.T - R2d_mean).T
    covariance = np.sum( X1d_pert * R2d_pert, axis=-1) / (number_of_ensembles-1)

    std_X = np.std( X1d, ddof=1 )
    std_R = np.std( R2d, ddof=1, axis=-1 )
    correlation = covariance / ( std_X * std_X ) #divides by the variance of the 1d array

    return correlation

#returns actual correlation for a 1d array and a 2d array
def cor2d_correlation(X1d,R2d,number_of_ensembles):
    number_of_ensembles = len(X1d)
    X1d_mean=np.mean(X1d)
    R2d_mean=np.mean(R2d, axis=-1)
    X1d_pert = X1d - X1d_mean
    R2d_pert = (R2d.T - R2d_mean).T
    covariance = np.sum( X1d_pert * R2d_pert, axis=-1) / (number_of_ensembles-1)

    std_X = np.std( X1d, ddof=1 )
    std_R = np.std( R2d, ddof=1, axis=-1 )
    correlation = covariance / ( std_X * std_R ) #divides by the variance of the 1d array

    return correlation

#returns actual correlation for a 2 1d arrays
def cor1d_correlation(X1,X2,number_of_ensembles):
    number_of_ensembles = len(X1)
    X1d_mean=np.mean(X1)
    R2d_mean=np.mean(X2)
    X1d_pert = X1 - X1d_mean
    R2d_pert = (X2 - R2d_mean)
    covariance = np.sum( X1d_pert * R2d_pert ) / (number_of_ensembles-1)

    std_X = np.std( X1, ddof=1 )
    std_R = np.std( X2, ddof=1 )
    correlation = covariance / ( std_X * std_R) #divides by the variance of x2

    return correlation

#returns correlation coefficient for 2 1d arrays
def cor1d(X1,X2,number_of_ensembles):
    number_of_ensembles = len(X1)
    X1d_mean=np.mean(X1)
    R2d_mean=np.mean(X2)
    X1d_pert = X1 - X1d_mean
    R2d_pert = (X2 - R2d_mean)
    covariance = np.sum( X1d_pert * R2d_pert ) / (number_of_ensembles-1)

    std_X = np.std( X1, ddof=1 )
    std_R = np.std( X2, ddof=1 )
    correlation = covariance / ( std_R * std_R) #divides by the variance of x2

    return correlation

#regression coefficients
alpha=cor1d(x1,x2,number_of_ensembles)
beta=cor1d(x2,y1,number_of_ensembles)
gamma=cor1d(x2,y2,number_of_ensembles)
delta=cor1d(x2,y3,number_of_ensembles)

#primes
x1_prime=x1-alpha*x2
y1_prime=y1-beta*x2
y2_prime=y2-gamma*x2
y3_prime=y3-delta*x2

print('x1 - x1 prime')
print(x1-x1_prime) #0 if the contribution from russia is zero
print('correlation of x1 prime and x2')
print(cor1d_correlation(x1_prime,x2,number_of_ensembles)) #should be 0
print('correlation of y1 prime and x2')
print(cor2d_correlation(x2,y1_prime,number_of_ensembles)) #should be 0
print('correlation of y2 prime and x2')
print(cor2d_correlation(x2,y2_prime,number_of_ensembles)) #should be 0
print('correlation of y3 prime and x2')
print(cor2d_correlation(x2,y3_prime,number_of_ensembles)) #should be 0

#plotting

#import packages
import matplotlib.pyplot as plt

#plotting

for i in range(number_of_days):
    fig,axs=plt.subplots(nrows=3,ncols=2,constrained_layout=True,figsize=(10,16))
    
    axs[0,0].plot(x1,y1[i])
    axs[0,1].plot(x1,y2[i])
    axs[1,0].plot(x1,y3[i])
    axs[1,1].plot(x2,y1[i])
    axs[2,0].plot(x2,y2[i])
    axs[2,1].plot(x2,y3[i])

    axs[0,0].set_xlabel('x1',fontsize=24)
    axs[0,0].set_ylabel('y1',fontsize=24)
    axs[1,0].set_xlabel('x1',fontsize=24)
    axs[1,0].set_ylabel('y3',fontsize=24)
    axs[2,0].set_xlabel('x2',fontsize=24)
    axs[2,0].set_ylabel('y2',fontsize=24)
    axs[0,0].set_title('smois vs wcasp, day '+str(i+1),fontsize=30)
    axs[1,0].set_title('smois vs waus, day '+str(i+1),fontsize=30)
    axs[2,0].set_title('russia vs secasp, day '+str(i+1),fontsize=30)
    axs[0,1].set_xlabel('x1',fontsize=24)
    axs[0,1].set_ylabel('y2',fontsize=24)
    axs[1,1].set_xlabel('x2',fontsize=24)
    axs[1,1].set_ylabel('y1',fontsize=24)
    axs[2,1].set_xlabel('x2',fontsize=24)
    axs[2,1].set_ylabel('y3',fontsize=24)
    axs[0,1].set_title('smois vs secasp, day '+str(i+1),fontsize=30)
    axs[1,1].set_title('russia vs wcasp, day '+str(i+1),fontsize=30)
    axs[2,1].set_title('russia vs waus, day '+str(i+1),fontsize=30)
    plt.savefig('/fs/scratch/PAS3252/yang/HONORS_THESIS/day_'+str(i+1)+'_plots_of_xs_and_ys.png')
    plt.close()

#convert -delay 50 -loop 0 /fs/scratch/PAS3252/yang/HONORS_THESIS/*_plots_of_xs_and_ys.png /fs/scratch/PAS3252/yang/HONORS_THESIS/plots_of_xs_and_ys.gif
#saves it to the scratch directory