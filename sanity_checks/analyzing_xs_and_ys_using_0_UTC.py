'''
Alex Yang's Honors Thesis
'''

import pickle
import numpy as np
from datetime import datetime, timedelta
import pandas as pd

number_of_ensembles=500
number_of_days=19
esa_dict=pickle.load(open('/fs/scratch/PAS3252/yang/HONORS_THESIS/x1_and_x2_at_0_utc.pkl','rb'))

itime=datetime.strptime('0714','%m%d')
Wcasp_averages=np.zeros((number_of_days,number_of_ensembles))
SEcasp_averages=np.zeros((number_of_days,number_of_ensembles))
Waus_averages=np.zeros((number_of_days,number_of_ensembles))
for i in range(number_of_days):
    current_day=itime.strftime('%m%d') #in the form of a string
    path='/fs/scratch/PAS3252/yang/HONORS_THESIS/xs_and_ys_at_'+current_day+'.pkl'
    current_dict=pickle.load(open(path,'rb'))
    Wcasp_averages[i,:]=current_dict['y1']
    SEcasp_averages[i,:]=current_dict['y2']
    Waus_averages[i,:]=current_dict['y3']
    itime+=timedelta(days=1)

x1=esa_dict['x1'] #(ensembles)
x2=esa_dict['x2'] #(ensembles)
y1=Wcasp_averages #(days,ensembles)
y2=SEcasp_averages #(days,ensembles)
y3=Waus_averages #(days,ensembles)

esa1_dict=pickle.load(open('ESA.pkl','rb'))
latCell=esa1_dict['lat_1d'] #40962 lat cells
lonCell=esa1_dict['lon_1d'] #40962 lon cells

#returns a scalar correlation (beteween -1 and 1) for a 2 1d arrays
def cor1d_correlation(X1,X2,number_of_ensembles):
    X1d_mean=np.mean(X1)
    R2d_mean=np.mean(X2)
    X1d_pert = X1 - X1d_mean
    R2d_pert = X2 - R2d_mean
    covariance = np.sum( X1d_pert * R2d_pert ) / (number_of_ensembles-1)
    std_X = np.std( X1, ddof=1 )
    std_R = np.std( X2, ddof=1 )
    correlation = covariance / ( std_X * std_R)
    return correlation

#returns a scalar regression coefficient (not necessarily beteween -1 and 1) for a 2 1d arrays
def cor1d_regcoef(X1,X2,number_of_ensembles):
    X1d_mean=np.mean(X1)
    R2d_mean=np.mean(X2)
    X1d_pert = X1 - X1d_mean
    R2d_pert = X2 - R2d_mean
    covariance = np.sum( X1d_pert * R2d_pert ) / (number_of_ensembles-1)
    std_R = np.std( X2, ddof=1 )
    regcoef = covariance / ( std_R * std_R)
    return regcoef

#returns scalar correlations equal to the number of days in the 2d array Y (days,ensembles)
#inputs:
#X: 2d array of size (1, ensembles)
#Y: 2d array of size (days, ensembles)
def cor2d_correlation(X,Y,number_of_ensembles):
    X_mean=np.mean(X,axis=-1)
    Y_mean=np.mean(Y,axis=-1)
    X_pert = (X.T - X_mean).T
    Y_pert = (Y.T - Y_mean).T
    covariance = np.sum( X_pert * Y_pert , axis=-1) / (number_of_ensembles-1)
    std_X = np.std( X, ddof=1 )
    std_R = np.std( Y, ddof=1 , axis=-1)
    correlation = covariance / ( std_X * std_R)
    return correlation

#returns regression coefficient equal to the number of days in the 2d array Y (days,ensembles)
#inputs:
#X: 2d array of size (1, ensembles)
#Y: 2d array of size (days, ensembles)
def cor2d_regcoef(X,Y,number_of_ensembles):
    X_mean=np.mean(X,axis=-1)
    Y_mean=np.mean(Y,axis=-1)
    X_pert = (X.T - X_mean).T
    Y_pert = (Y.T - Y_mean).T
    covariance = np.sum( X_pert * Y_pert, axis=-1) / (number_of_ensembles-1)
    std_Y = np.std( Y, ddof=1, axis=-1)
    regcoef = covariance / ( std_Y * std_Y)
    return regcoef

#regression coefficients
alpha=cor1d_regcoef(x1,x2,number_of_ensembles)
beta=cor2d_regcoef(y1,x2,number_of_ensembles)
gamma=cor2d_regcoef(y2,x2,number_of_ensembles)
delta=cor2d_regcoef(y3,x2,number_of_ensembles)

#Dataframe for x1 vs unprimed varirables
data={
    'x1 vs. y1': np.empty(number_of_days,dtype='float'),
    'x1 vs. y2': np.empty(number_of_days,dtype='float'),
    'x1 vs. y3': np.empty(number_of_days,dtype='float')
}

for i in range(number_of_days):
    data['x1 vs. y1'][i]=cor2d_correlation(x1,y1[i,:],number_of_ensembles)
    data['x1 vs. y2'][i]=cor2d_correlation(x1,y2[i,:],number_of_ensembles)
    data['x1 vs. y3'][i]=cor2d_correlation(x1,y3[i,:],number_of_ensembles)

df=pd.DataFrame(data)

df.index=['Day ' + str(i) for i in range(1,len(df)+1)] #got some help from ChatGPT for this one

print(df)

#Dataframe for x1 vs prime variables

y1_prime=np.zeros((19,number_of_ensembles))
y2_prime=np.zeros((19,number_of_ensembles))
y3_prime=np.zeros((19,number_of_ensembles))

for i in range(number_of_days):
    y1_prime[i,:]=y1[i,:]-beta[i]*x2
    y2_prime[i,:]=y2[i,:]-gamma[i]*x2
    y3_prime[i,:]=y3[i,:]-delta[i]*x2

data={
    'x1 vs. y1_prime': np.empty(number_of_days,dtype='float'),
    'x1 vs. y2_prime': np.empty(number_of_days,dtype='float'),
    'x1 vs. y3_prime': np.empty(number_of_days,dtype='float')
}

for i in range(number_of_days):
    data['x1 vs. y1_prime'][i]=cor2d_correlation(x1,y1_prime[i,:],number_of_ensembles)
    data['x1 vs. y2_prime'][i]=cor2d_correlation(x1,y2_prime[i,:],number_of_ensembles)
    data['x1 vs. y3_prime'][i]=cor2d_correlation(x1,y3_prime[i,:],number_of_ensembles)

df=pd.DataFrame(data)

df.index=['Day ' + str(i) for i in range(1,len(df)+1)] #got some help from ChatGPT for this one

print(df)

#Dataframe for x2 vs unprimed varirables
data={
    'x2 vs. y1': np.empty(number_of_days,dtype='float'),
    'x2 vs. y2': np.empty(number_of_days,dtype='float'),
    'x2 vs. y3': np.empty(number_of_days,dtype='float')
}

for i in range(number_of_days):
    data['x2 vs. y1'][i]=cor2d_correlation(x2,y1[i,:],number_of_ensembles)
    data['x2 vs. y2'][i]=cor2d_correlation(x2,y2[i,:],number_of_ensembles)
    data['x2 vs. y3'][i]=cor2d_correlation(x2,y3[i,:],number_of_ensembles)

df=pd.DataFrame(data)

df.index=['Day ' + str(i) for i in range(1,len(df)+1)] #got some help from ChatGPT for this one

print(df)

#Dataframe for x2 vs primed variables
y1_prime=np.zeros((19,number_of_ensembles))
y2_prime=np.zeros((19,number_of_ensembles))
y3_prime=np.zeros((19,number_of_ensembles))
#primes
x1_prime=x1-alpha*x2
for i in range(number_of_days):
    y1_prime[i,:]=y1[i,:]-beta[i]*x2
    y2_prime[i,:]=y2[i,:]-gamma[i]*x2
    y3_prime[i,:]=y3[i,:]-delta[i]*x2

data={
    'x2 vs. y1_prime': np.empty(number_of_days,dtype='float'),
    'x2 vs. y2_prime': np.empty(number_of_days,dtype='float'),
    'x2 vs. y3_prime': np.empty(number_of_days,dtype='float')
}

for i in range(number_of_days):
    data['x2 vs. y1_prime'][i]=cor2d_correlation(x2,y1_prime[i,:],number_of_ensembles)
    data['x2 vs. y2_prime'][i]=cor2d_correlation(x2,y2_prime[i,:],number_of_ensembles)
    data['x2 vs. y3_prime'][i]=cor2d_correlation(x2,y3_prime[i,:],number_of_ensembles)

df=pd.DataFrame(data)

df.index=['Day ' + str(i) for i in range(1,len(df)+1)] #got some help from ChatGPT for this one

print(df)

print('correlation between x1 and x1 prime')
print(cor1d_correlation(x1,x1_prime,number_of_ensembles)) #0 if the contribution from russia is zero
print('correlation of x1 prime and x2')
print(np.mean(cor1d_correlation(x1_prime,x2,number_of_ensembles))) #should be 0
print('correlation of x1 prime and x1')
print(np.mean(cor1d_correlation(x1_prime,x1,number_of_ensembles))) #should be 0

y1_prime_and_x2=cor2d_correlation(x2,y1_prime,number_of_ensembles) #should be 0
y2_prime_and_x2=cor2d_correlation(x2,y2_prime,number_of_ensembles) #should be 0
y3_prime_and_x2=cor2d_correlation(x2,y3_prime,number_of_ensembles) #should be 0

print('correlation of y1 prime and x2')
print(np.mean(y1_prime_and_x2)) #should be 0
print('correlation of y2 prime and x2')
print(np.mean(y2_prime_and_x2)) #should be 0
print('correlation of y3 prime and x2')
print(np.mean(y3_prime_and_x2)) #should be 0

#finding variances of x1 and x1 prime
print('variance of x1')
print(np.var(x1))
print('variance of x1 prime')
print(np.var(x1_prime))

print('correlation between x1 and x2')
print(cor1d_correlation(x1,x2,number_of_ensembles))

y1s=np.zeros(number_of_days)
y2s=np.zeros(number_of_days)
y3s=np.zeros(number_of_days)

for i in range(number_of_days):
    correlation1=cor2d_correlation(x1_prime,y1_prime[i,:],number_of_ensembles)
    correlation2=cor2d_correlation(x1_prime,y2_prime[i,:],number_of_ensembles)
    correlation3=cor2d_correlation(x1_prime,y3_prime[i,:],number_of_ensembles)

    y1s[i]=correlation1
    y2s[i]=correlation2
    y3s[i]=correlation3
    
    print('for day '+str(i+1)+', correlations between x1 prime and y primes:')
    print(correlation1, correlation2, correlation3)

for i in range(number_of_days):
    correlation1=cor2d_correlation(x1_prime,y1[i,:],number_of_ensembles)
    correlation2=cor2d_correlation(x1_prime,y2[i,:],number_of_ensembles)
    correlation3=cor2d_correlation(x1_prime,y3[i,:],number_of_ensembles)

    y1s[i]=correlation1
    y2s[i]=correlation2
    y3s[i]=correlation3
    
    print('for day '+str(i+1)+', correlations between x1 prime and ys:')
    print(correlation1, correlation2, correlation3)


#find the regression coefficients with contribution from smois removed as sanity check
#regression coefficients
beta_star=cor2d_regcoef(y1,x1,number_of_ensembles)
gamma_star=cor2d_regcoef(y2,x1,number_of_ensembles)
delta_star=cor2d_regcoef(y3,x1,number_of_ensembles)

#Dataframe for x2 vs primed variables
y1_star=np.zeros((19,number_of_ensembles))
y2_star=np.zeros((19,number_of_ensembles))
y3_star=np.zeros((19,number_of_ensembles))
#primes
for i in range(number_of_days):
    y1_star[i,:]=y1[i,:]-beta_star[i]*x1
    y2_star[i,:]=y2[i,:]-gamma_star[i]*x1
    y3_star[i,:]=y3[i,:]-delta_star[i]*x1

data={
    'x1 vs. y1_star': np.empty(number_of_days,dtype='float'),
    'x1 vs. y2_star': np.empty(number_of_days,dtype='float'),
    'x1 vs. y3_star': np.empty(number_of_days,dtype='float')
}

for i in range(number_of_days):
    data['x1 vs. y1_star'][i]=cor2d_correlation(x1,y1_star[i,:],number_of_ensembles)
    data['x1 vs. y2_star'][i]=cor2d_correlation(x1,y2_star[i,:],number_of_ensembles)
    data['x1 vs. y3_star'][i]=cor2d_correlation(x1,y3_star[i,:],number_of_ensembles)

df=pd.DataFrame(data)

df.index=['Day ' + str(i) for i in range(1,len(df)+1)] #got some help from ChatGPT for this one

print(df)

#inputs:
#x1 smois 1d
#x2 z500 1d
#y1 z500 2d
#y2 z500 2d
#y3 z500 2d
#type can be 1, 2, 3, or 4
#case 1: x1 vs unprimed ys
#case 2: x1 vs primed ys
#case 3: x2 vs unprimed ys
#case 4: x2 vs primed ys

#type2 can be 1 or 2
#case 1: primes (remove contribution from x2)
#case 2: stars (remove contribution from x1)

#output:
#prints a dataframe
def turn_into_dataframe(x1,x2,y1,y2,y3,type,type2):
    possibiltiies=[(False,False),(False,True),(True,False),(True,True)]
    prime_or_star=[(True),(False)]
    is_it_prime=prime_or_star[type2-1]
    is_it_x2,is_it_primed=possibiltiies[type-1]
    print('is it x2')
    print(is_it_x2)
    print('is it primed')
    print(is_it_primed)
    print('is it prime')
    print(is_it_prime)

    #Dataframe for x2 vs primed variables
    y1_star=np.zeros((19,number_of_ensembles))
    y2_star=np.zeros((19,number_of_ensembles))
    y3_star=np.zeros((19,number_of_ensembles))
    #primes
    if is_it_prime:
        #regression coefficients
        beta_star=cor2d_regcoef(y1,x2,number_of_ensembles)
        gamma_star=cor2d_regcoef(y2,x2,number_of_ensembles)
        delta_star=cor2d_regcoef(y3,x2,number_of_ensembles)
        for i in range(number_of_days):
            y1_star[i,:]=y1[i,:]-beta_star[i]*x2
            y2_star[i,:]=y2[i,:]-gamma_star[i]*x2
            y3_star[i,:]=y3[i,:]-delta_star[i]*x2
        suffix='_prime'
    if is_it_prime==False:
        #regression coefficients
        beta_star=cor2d_regcoef(y1,x1,number_of_ensembles)
        gamma_star=cor2d_regcoef(y2,x1,number_of_ensembles)
        delta_star=cor2d_regcoef(y3,x1,number_of_ensembles)
        for i in range(number_of_days):
            y1_star[i,:]=y1[i,:]-beta_star[i]*x1
            y2_star[i,:]=y2[i,:]-gamma_star[i]*x1
            y3_star[i,:]=y3[i,:]-delta_star[i]*x1
        suffix='_star'

    if is_it_x2:
        which_x='x2'
    if is_it_x2==False:
        which_x='x1'

    if is_it_primed:
        suffix=suffix
    if is_it_primed==False:
        suffix=''
    
    data={
        which_x+' vs. y1'+suffix: np.empty(number_of_days,dtype='float'),
        which_x+' vs. y2'+suffix: np.empty(number_of_days,dtype='float'),
        which_x+' vs. y3'+suffix: np.empty(number_of_days,dtype='float')
    }

    if is_it_primed==False:
        y1_star=y1
        y2_star=y2
        y3_star=y3

    if is_it_x2==False:
        for i in range(number_of_days):
            data[which_x+' vs. y1'+suffix][i]=cor2d_correlation(x1,y1_star[i,:],number_of_ensembles)
            data[which_x+' vs. y2'+suffix][i]=cor2d_correlation(x1,y2_star[i,:],number_of_ensembles)
            data[which_x+' vs. y3'+suffix][i]=cor2d_correlation(x1,y3_star[i,:],number_of_ensembles)

    if is_it_x2==True:
        for i in range(number_of_days):
            data[which_x+' vs. y1'+suffix][i]=cor2d_correlation(x2,y1_star[i,:],number_of_ensembles)
            data[which_x+' vs. y2'+suffix][i]=cor2d_correlation(x2,y2_star[i,:],number_of_ensembles)
            data[which_x+' vs. y3'+suffix][i]=cor2d_correlation(x2,y3_star[i,:],number_of_ensembles)

    df=pd.DataFrame(data)

    df.index=['Day ' + str(i) for i in range(1,len(df)+1)] #got some help from ChatGPT for this one

    print(df)

#type can be 1, 2, 3, or 4
#case 1: x1 vs unprimed ys
#case 2: x1 vs primed ys
#case 3: x2 vs unprimed ys
#case 4: x2 vs primed ys

#type2 can be 1 or 2
#case 1: primes (remove contribution from x2)
#case 2: stars (remove contribution from x1)

print('x1 vs unprimed')
turn_into_dataframe(x1,x2,y1,y2,y3,1,1)
print('x1 vs primed primes')
turn_into_dataframe(x1,x2,y1,y2,y3,2,1)
print('x2 vs unprimed')
turn_into_dataframe(x1,x2,y1,y2,y3,3,1)
print('x2 vs primed primes')
turn_into_dataframe(x1,x2,y1,y2,y3,4,1)
print('x1 vs unprimed')
turn_into_dataframe(x1,x2,y1,y2,y3,1,2)
print('x1 vs primed stars')
turn_into_dataframe(x1,x2,y1,y2,y3,2,2)
print('x2 vs unprimed')
turn_into_dataframe(x1,x2,y1,y2,y3,3,2)
print('x2 vs primed stars')
turn_into_dataframe(x1,x2,y1,y2,y3,4,2)
#uncomment the print(df) statement in the function

'''
#plotting

#import packages
import matplotlib.pyplot as plt
from datetime import datetime,timedelta

x1=np.arange(1,20)
x2=x1
y1=y1s
y2=y2s
y3=y3s

fig,axs=plt.subplots(nrows=3,ncols=2,constrained_layout=True,figsize=(15,15))

axs[0,0].scatter(x1,y1)
axs[1,0].scatter(x1,y2)
axs[2,0].scatter(x1,y3)
axs[0,1].scatter(x2,y1)
axs[1,1].scatter(x2,y2)
axs[2,1].scatter(x2,y3)

axs[0,0].set_xlabel('time',fontsize=24)
axs[1,0].set_xlabel('time',fontsize=24)
axs[2,0].set_xlabel('time',fontsize=24)

axs[0,1].set_xlabel('time',fontsize=24)
axs[1,1].set_xlabel('time',fontsize=24)
axs[2,1].set_xlabel('time',fontsize=24)

axs[0,0].set_ylabel('x1\' vs y1\'',fontsize=24)
axs[1,0].set_ylabel('x1\' vs y2\'',fontsize=24)
axs[2,0].set_ylabel('x1\' vs y3\'',fontsize=24)

axs[0,1].set_ylabel('x2 vs y1\'',fontsize=24)
axs[1,1].set_ylabel('x2 vs y2\'',fontsize=24)
axs[2,1].set_ylabel('x2 vs y3\'',fontsize=24)

axs[0,0].set_title('smois vs wcasp',fontsize=30)
axs[1,0].set_title('smois vs secasp',fontsize=30)
axs[2,0].set_title('smois vs waus',fontsize=30)

axs[0,1].set_title('russia vs wcasp',fontsize=30)
axs[1,1].set_title('russia vs secasp',fontsize=30)
axs[2,1].set_title('russia vs waus',fontsize=30)

plt.savefig('/fs/scratch/PAS3252/yang/HONORS_THESIS/time_series_primes_utc_0.png')
plt.close()'''

'''
for i in range(number_of_days):
    day=itime.strftime('%m/%d')
    fig,axs=plt.subplots(nrows=3,ncols=2,constrained_layout=True,figsize=(15,15))
    
    axs[0,0].scatter(x1,y1[i])
    axs[0,0].set_ylim(5500,6000)
    axs[1,0].scatter(x1,y2[i])
    axs[1,0].set_ylim(5500,6000)
    axs[2,0].scatter(x1,y3[i])
    axs[2,0].set_ylim(5500,6000)

    axs[0,1].scatter(x2,y1[i])
    axs[0,1].set_ylim(5500,6000)
    axs[1,1].scatter(x2,y2[i])
    axs[1,1].set_ylim(5500,6000)
    axs[2,1].scatter(x2,y3[i])
    axs[2,1].set_ylim(5500,6000)

    axs[0,0].set_xlabel('x1',fontsize=24)
    axs[1,0].set_xlabel('x1',fontsize=24)
    axs[2,0].set_xlabel('x1',fontsize=24)

    axs[0,1].set_xlabel('x2',fontsize=24)
    axs[1,1].set_xlabel('x2',fontsize=24)
    axs[2,1].set_xlabel('x2',fontsize=24)

    axs[0,0].set_ylabel('y1',fontsize=24)
    axs[1,0].set_ylabel('y2',fontsize=24)
    axs[2,0].set_ylabel('y3',fontsize=24)

    axs[0,1].set_ylabel('y1',fontsize=24)
    axs[1,1].set_ylabel('y2',fontsize=24)
    axs[2,1].set_ylabel('y3',fontsize=24)

    axs[0,0].set_title('smois vs wcasp, '+day,fontsize=30)
    axs[1,0].set_title('smois vs secasp, '+day,fontsize=30)
    axs[2,0].set_title('smois vs waus, '+day,fontsize=30)

    axs[0,1].set_title('russia vs wcasp, '+day,fontsize=30)
    axs[1,1].set_title('russia vs secasp, '+day,fontsize=30)
    axs[2,1].set_title('russia vs waus, '+day,fontsize=30)

    plt.savefig('/fs/scratch/PAS3252/yang/HONORS_THESIS/day_'+str(i+1)+'_plots_of_xs_and_ys_utc_0.png')
    plt.close()
    itime+=timedelta(days=1)

#convert -delay 50 -loop 0 /fs/scratch/PAS3252/yang/HONORS_THESIS/*ys_utc_0.png /fs/scratch/PAS3252/yang/HONORS_THESIS/plots_of_xs_and_ys_at_utc_0.gif
#saves it to the scratch directory

x1=x1_prime
y1=y1_prime
y2=y2_prime
y3=y3_prime

itime=datetime.strptime('0714','%m%d')

for i in range(number_of_days):
    day=itime.strftime('%m/%d')
    fig,axs=plt.subplots(nrows=3,ncols=2,constrained_layout=True,figsize=(15,15))
    
    axs[0,0].scatter(x1,y1[i])
    axs[1,0].scatter(x1,y2[i])
    axs[2,0].scatter(x1,y3[i])
    axs[0,1].scatter(x2,y1[i])
    axs[1,1].scatter(x2,y2[i])
    axs[2,1].scatter(x2,y3[i])

    axs[0,0].set_xlabel('x1\'',fontsize=24)
    axs[1,0].set_xlabel('x1\'',fontsize=24)
    axs[2,0].set_xlabel('x1\'',fontsize=24)

    axs[0,1].set_xlabel('x2',fontsize=24)
    axs[1,1].set_xlabel('x2',fontsize=24)
    axs[2,1].set_xlabel('x2',fontsize=24)

    axs[0,0].set_ylabel('y1\'',fontsize=24)
    axs[1,0].set_ylabel('y2\'',fontsize=24)
    axs[2,0].set_ylabel('y3\'',fontsize=24)

    axs[0,1].set_ylabel('y1\'',fontsize=24)
    axs[1,1].set_ylabel('y2\'',fontsize=24)
    axs[2,1].set_ylabel('y3\'',fontsize=24)

    axs[0,0].set_title('smois vs wcasp, '+day,fontsize=30)
    axs[1,0].set_title('smois vs secasp, '+day,fontsize=30)
    axs[2,0].set_title('smois vs waus, '+day,fontsize=30)

    axs[0,1].set_title('russia vs wcasp, '+day,fontsize=30)
    axs[1,1].set_title('russia vs secasp, '+day,fontsize=30)
    axs[2,1].set_title('russia vs waus, '+day,fontsize=30)

    plt.savefig('/fs/scratch/PAS3252/yang/HONORS_THESIS/day_'+str(i+1)+'_plots_of_xs_and_ys_primes_utc_0.png')
    plt.close()
    itime+=timedelta(days=1)

#convert -delay 50 -loop 0 /fs/scratch/PAS3252/yang/HONORS_THESIS/*primes_utc_0.png /fs/scratch/PAS3252/yang/HONORS_THESIS/plots_of_xs_and_ys_primes_at_utc_0.gif
#saves it to the scratch directory'''