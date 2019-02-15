#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 12:21:46 2018

@author: phillips
"""

import os
os.chdir('/Users/phillips/Documents/Memory project/Stan analysis of families')


import pystan
import pandas as pd
from numpy import *
import matplotlib.pyplot as plt 

#plt.rcParams.update({'font.size': 10})
#plt.rc('figure', titlesize=10)

FONT_SIZE = 10

plt.rc('font', size=FONT_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=FONT_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=FONT_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=FONT_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=FONT_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=FONT_SIZE)    # legend fontsize
plt.rc('figure', titlesize=FONT_SIZE)  # fontsize of the figure title

from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'qt')
#get_ipython().run_line_magic('matplotlib', 'inline')
plt.style.use('seaborn-white')

from scipy.stats import beta

import time

import pickle

from heapq import nlargest

#%%

Dstn = pd.read_excel('/Users/phillips/Documents/Memory project/Stan analysis of families/Dstnfamilies.xlsx')
Dstn_NumCell = Dstn['NumCell'].values
Dstn_y = Dstn['y'].values
Dstn_y = Dstn_y/mean(Dstn_y)
Dstn_FamNum = Dstn['FamNum'].values
Dstn_FamTot = 24

Jam2 = pd.read_excel('/Users/phillips/Documents/Memory project/Stan analysis of families/Jam2families.xlsx')
Jam2_NumCell = Jam2['NumCell'].values
Jam2_y = Jam2['y'].values
Jam2_y = Jam2_y/mean(Jam2_y)
Jam2_FamNum = Jam2['FamNum'].values
Jam2_FamTot = 35

pGK = pd.read_excel('/Users/phillips/Documents/Memory project/Stan analysis of families/pGKfamilies.xlsx')
pGK_NumCell = pGK['NumCell'].values
pGK_y = pGK['y'].values
pGK_y = pGK_y/mean(pGK_y)
pGK_FamNum = pGK['FamNum'].values
pGK_FamTot = 25

#%%

model1 = """
data{
    int<lower=1> N1;
    int Dstn_NumCell[N1];
    real Dstn_y[N1];    
    int Dstn_FamNum[N1];
    int Dstn_FamTot;
    int<lower=1> N2;
    int Jam2_NumCell[N2];
    real Jam2_y[N2];    
    int Jam2_FamNum[N2];
    int Jam2_FamTot;
    int<lower=1> N3;
    int pGK_NumCell[N3];
    real pGK_y[N3];    
    int pGK_FamNum[N3];
    int pGK_FamTot;
}
parameters{
    real Dstn_mean_P;
    real<lower=0> Dstn_sigma_P;
    real<lower=0> Dstn_sigma_S;
    real<lower=0> Dstn_mean_F[Dstn_FamTot];
    real Jam2_mean_P;
    real<lower=0> Jam2_sigma_P;
    real<lower=0> Jam2_sigma_S;
    real<lower=0> Jam2_mean_F[Jam2_FamTot];
    real pGK_mean_P;
    real<lower=0> pGK_sigma_P;
    real<lower=0> pGK_sigma_S;
    real<lower=0> pGK_mean_F[pGK_FamTot];
}
model{
    vector[N1] mu1;
    vector[N2] mu2;
    vector[N3] mu3;
    for ( i in 1:N1 ) {
        mu1[i] = Dstn_mean_F[Dstn_FamNum[i]];
        Dstn_y[i] ~ normal(mu1[i], Dstn_sigma_S);
    }     
    for ( i in 1:Dstn_FamTot ) {
        Dstn_mean_F[i] ~ lognormal(Dstn_mean_P, Dstn_sigma_P);
    } 
    for ( i in 1:N2 ) {
        mu2[i] = Jam2_mean_F[Jam2_FamNum[i]];
        Jam2_y[i] ~ normal(mu2[i], Jam2_sigma_S);
    }     
    for ( i in 1:Jam2_FamTot ) {
        Jam2_mean_F[i] ~ lognormal(Jam2_mean_P, Jam2_sigma_P);
    } 
    for ( i in 1:N3 ) {
        mu3[i] = pGK_mean_F[pGK_FamNum[i]];
        pGK_y[i] ~ normal(mu3[i], pGK_sigma_S);
    } 
        for ( i in 1:pGK_FamTot ) {
        pGK_mean_F[i] ~ lognormal(pGK_mean_P, pGK_sigma_P);
    } 
}
generated quantities{
    vector[N1+N2+N3] log_lik;
    for ( i in 1:N1 ) {
        log_lik[i] = normal_lpdf(Dstn_y[i] | Dstn_mean_F[Dstn_FamNum[i]], Dstn_sigma_S);
    } 
    for ( i in 1:N2 ) {
        log_lik[i+N1] = normal_lpdf(Jam2_y[i] | Jam2_mean_F[Jam2_FamNum[i]], Jam2_sigma_S);
    } 
    for ( i in 1:N3 ) {
        log_lik[i+N1+N2] = normal_lpdf(pGK_y[i] | pGK_mean_F[pGK_FamNum[i]], pGK_sigma_S);
    } 
}

"""
#%%
Dstn_mean_P = mean(log(Dstn_y))
Dstn_sigma_P = std(log(Dstn_y))
Dstn_mean_F = zeros(Dstn_FamTot)
Dstn_std = zeros(Dstn_FamTot)
for i in range(0,Dstn_FamTot):
    Dstn_mean_F[i] = mean(Dstn_y[Dstn_FamNum==i+1])
    Dstn_std[i] = std(Dstn_y[Dstn_FamNum==i+1])
Dstn_sigma_S = mean(Dstn_std)

Jam2_mean_P = mean(log(Jam2_y))
Jam2_sigma_P = std(log(Jam2_y))
Jam2_mean_F = zeros(Jam2_FamTot)
Jam2_std = zeros(Jam2_FamTot)
for i in range(0,Jam2_FamTot):
    Jam2_mean_F[i] = mean(Jam2_y[Jam2_FamNum==i+1])
    Jam2_std[i] = std(Jam2_y[Jam2_FamNum==i+1])
Jam2_sigma_S = mean(Jam2_std)

pGK_mean_P = mean(log(pGK_y))
pGK_sigma_P = std(log(pGK_y))
pGK_mean_F = zeros(pGK_FamTot)
pGK_std = zeros(pGK_FamTot)
for i in range(0,pGK_FamTot):
    pGK_mean_F[i] = mean(pGK_y[pGK_FamNum==i+1])
    pGK_std[i] = std(pGK_y[pGK_FamNum==i+1])
pGK_sigma_S = mean(pGK_std)

#%%

dat = {
    'N1' : len(Dstn_NumCell),
    'Dstn_NumCell' : Dstn_NumCell,
    'Dstn_y' : Dstn_y,   
    'Dstn_FamNum' : Dstn_FamNum,
    'Dstn_FamTot' : Dstn_FamTot,
    'N2' : len(Jam2_NumCell),
    'Jam2_NumCell' : Jam2_NumCell,
    'Jam2_y' : Jam2_y,   
    'Jam2_FamNum' : Jam2_FamNum,
    'Jam2_FamTot' : Jam2_FamTot,
    'N3' : len(pGK_NumCell),
    'pGK_NumCell' : pGK_NumCell,
    'pGK_y' : pGK_y,   
    'pGK_FamNum' : pGK_FamNum,
    'pGK_FamTot' : pGK_FamTot
}

def initfun1():
    return dict(Dstn_sigma_S=Dstn_sigma_S,Dstn_mean_F=Dstn_mean_F,Jam2_mean_P=Jam2_mean_P,Jam2_sigma_P=Jam2_sigma_P,Jam2_sigma_S=Jam2_sigma_S,Jam2_mean_F=Jam2_mean_F,pGK_mean_P=pGK_mean_P,pGK_sigma_P=pGK_sigma_P,pGK_sigma_S=pGK_sigma_S,pGK_mean_F=pGK_mean_F)

sm = pystan.StanModel(model_code=model1)
currenttime = time.time()
#fit = sm.sampling(data=dat, iter=10000, chains=4)
fit = sm.sampling(data=dat, init=initfun1,iter=10000, chains=4)
elapsed = time.time() - currenttime      
print(elapsed) 

#%%

a = fit.extract(permuted=True)

#%%

pickle_out = open("GeneSpecific.pickle","wb")
pickle.dump(a, pickle_out)
pickle_out.close()

#%%

pickle_in = open("GeneSpecific.pickle","rb")
a = pickle.load(pickle_in)

#%%

# trace plots

fig = plt.figure(figsize=(20,16))

plt.subplot(341)
param = a['Dstn_mean_P']
plt.plot(arange(0,len(param)), param)
plt.xlabel('Iterations')
#plt.ylabel('Frequency')
#plt.title('\alpha Bmal1')
plt.title(r'$\alpha$ Bmal1')

plt.subplot(342)
param = a['Dstn_sigma_P']
plt.plot(arange(0,len(param)), param)
plt.xlabel('Iterations')
#plt.ylabel('Frequency')
plt.title(r'$\epsilon$ Bmal1')
plt.ylim(0, 1)

plt.subplot(343)
param = a['Dstn_sigma_S']
plt.plot(arange(0,len(param)), param)
plt.xlabel('Iterations')
#plt.ylabel('Frequency')
plt.title(r'$\phi$ Bmal1')
#plt.ylim(-pi, pi)

plt.subplot(3,4,4)
Dstn_mean_F = a['Dstn_mean_F']
Dstn_mean_F = Dstn_mean_F[:,4]
plt.plot(arange(0,len(Dstn_mean_F)), Dstn_mean_F)
plt.xlabel('Iterations')
#plt.ylabel('Frequency')
plt.title(r'$\phi$ cell 5')

#%%
#
#plt.subplot(346)
#param = a['epsilonCry1']
#plt.plot(arange(0,len(param)), param)
#plt.xlabel('Iterations')
##plt.ylabel('Frequency')
#plt.title(r'$\epsilon$ Cry1')
#plt.ylim(0, 1)
#
#plt.subplot(347)
#param = a['phiCry1']
#plt.plot(arange(0,len(param)), param)
#plt.xlabel('Iterations')
##plt.ylabel('Frequency')
#plt.title(r'$\phi$ Cry1')
#plt.ylim(-pi, pi)
#
#
#plt.subplot(349)
#param = a['alphaNr1d1']
#plt.plot(arange(0,len(param)), param)
#plt.xlabel('Iterations')
##plt.ylabel('Frequency')
##plt.title('\alpha Bmal1')
#plt.title(r'$\alpha$ Nr1d1')
#
#plt.subplot(3,4,10)
#param = a['epsilonNr1d1']
#plt.plot(arange(0,len(param)), param)
#plt.xlabel('Iterations')
##plt.ylabel('Frequency')
#plt.title(r'$\epsilon$ Nr1d1')
#plt.ylim(0, 1)
#
#plt.subplot(3,4,11)
#param = a['phiNr1d1']
#plt.plot(arange(0,len(param)), param)
#plt.xlabel('Iterations')
##plt.ylabel('Frequency')
#plt.title(r'$\phi$ Nr1d1')
#plt.ylim(-pi, pi)
#
#plt.subplot(344)
#param = a['beta_v']
#plt.plot(arange(0,len(param)), param)
#plt.xlabel('Iterations')
##plt.ylabel('Frequency')
#plt.title('beta_v')
#
#plt.subplot(348)
#param = a['kappa']
#plt.plot(arange(0,len(param)), param)
#plt.xlabel('Iterations')
##plt.ylabel('Frequency')
#plt.title(r'$\kappa$')
#
#plt.subplot(3,4,12)
#phi__1_add = a['phi_i_add_RC']
#phi__1_add = phi__1_add[:,4]
#plt.plot(arange(0,len(phi__1_add)), phi__1_add)
#plt.xlabel('Iterations')
##plt.ylabel('Frequency')
#plt.title(r'$\phi$ cell 5')
#
##plt.subplot(3,4,12)
##param = a['lik']
##plt.plot(arange(0,len(param)), param)
##plt.xlabel('Iterations')
###plt.ylabel('Frequency')
##plt.title('phi_1')
