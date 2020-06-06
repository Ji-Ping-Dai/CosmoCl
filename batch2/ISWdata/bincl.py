# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 08:59:21 2019

@author: Ji-Ping Dai
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

alpha1=1.1
cll=np.array([631,1202.3,2290.9,4786.3,9120.1,17378.0,36307.8,69183.1,131825.7,1000000])
fac=np.zeros(9)
for i in range(9):
    fac[i]=(pow(cll[0],-alpha1)-pow(cll[1],-alpha1))/(pow(cll[i],-alpha1)-pow(cll[i+1],-alpha1))
#for i in range(3,9):
#    fac[i]=(pow(cl[0],-alpha2)-pow(cl[1],-alpha2))/((pow(cl[i],-alpha2)-pow(cl[i+1],-alpha2)))

for i in range(9):
    cov=fits.open('cov'+str(i)+'.fits')
    cov=np.array((cov[0].data))[0]
    covnew=np.zeros([511,511])
    for j in range(511):
        for k in range(511):
            covnew[j][k]=cov[j+2][k+2]*fac[i]*fac[i]
    incov=np.linalg.inv(covnew)
    np.savetxt('bin'+str(i+1)+'incov.dat',incov)
    
for i in range(9):
    cl=np.loadtxt('clccf'+str(i)+'.cl')
    clnew=np.zeros([511,2])
    clnew[:,0]=cl[2:513,0]
    clnew[:,1]=cl[2:513,1]*fac[i]
    np.savetxt('bin'+str(i+1)+'cl.dat',clnew)
    
cl1=np.loadtxt('bin1cl.dat')
cl2=np.loadtxt('bin8cl.dat')
plt.plot(cl1[:,0],cl1[:,1])
plt.plot(cl2[:,0],cl2[:,1])
plt.axhline(0,ls=':',color='r')
plt.xlim(10,100)
plt.ylim(-4e-14,4e-14)
plt.xscale('log')

#cl=np.zeros([9,511])
#for i in range(9):
#    cl[i,:]=np.loadtxt('bin'+str(i+1)+'cl.dat')[:,1]
#
#ell=np.arange(2,513)
#fig = plt.figure(figsize=[10,7])
#plt.xlabel(r'$\ell$',size=25)
#plt.ylabel(r'$\rm Power~Spectrum$',size=25)
#plt.xticks(fontsize=25)
#plt.yticks(fontsize=25)
#plt.xscale('log')
##plt.plot(kh,skew,color='b',lw=2,ls=':',label=r'$\rm Initial$') 
##plt.plot(kh,skew_gra,color='g',lw=2,ls='-.',label=r'$\rm Gravity$')
##plt.plot(kh,skew_bias,color='r',lw=2,ls='--',label=r'$\rm Nonlinear~bias$')
##plt.plot(kh,skew+skew_gra+skew_bias,color='k',lw=3,label=r'$\rm Total$')
#
#plt.xlim(10,512)
#plt.ylim(-0.5e-13,0.5e-13)
#plt.plot(ell, cl[0,:],color='b',lw=2,label=r'$\rm bin_1$')
#plt.plot(ell, cl[4,:],color='r',lw=2,label=r'$\rm bin_5$')
#plt.plot(ell, cl[8,:],color='g',lw=2,label=r'$\rm bin_9$')
#plt.legend(loc='upper center',fontsize=20)
#
#
#
