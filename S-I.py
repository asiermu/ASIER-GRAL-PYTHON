# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 15:28:55 2023

@author: asier
"""

#T desberdinen arteko desberdintasuna ikusteko programa
#Aurreko RUNGE-KUTTA programaren oso oso antzekoa
import numpy as np
import matplotlib.pyplot as plt
import math
T1,T2,d,Y0=20,300,3,np.array([3900000,100000,0])
q,A,b,tau=4,np.array([[0,0,0,0],[1/2,0,0,0],[0,1/2,0,0],[0,0,1,0]]),np.array([1/6,1/3,1/3,1/6]),np.array([0,1/2,1/2,1])
n=1000
mu=4
beta=4
eta=3
gamma=2
N=4000000
def func1(t,z,mu,beta,eta,gamma,N):
    return np.array([mu*N-mu*z[0]-beta*z[0]*z[1]/N,beta*z[0]*z[1]/N-z[1]*mu-gamma*z[1]+z[2]*eta,gamma*z[1]-z[2]*mu-z[2]*eta])

def RK(n,T,Y0,d,mu,beta,eta,gamma,N,q,A,b,tau):
    dt=T/n
    tn=np.zeros(n+1)                    
    for i in range(n+1):
        tn[i]=i*dt                 
    Y=np.zeros((d,n+1))                 
    Y[:,0]=Y0[:]
   
    rho,K=np.zeros(q),np.zeros((d,q))       
    for a in range(n):
        for j in range(q):
            rho[j]=tn[a]+dt*tau[j]
        K[:,0]=func1(rho[0],Y[:,a],mu,beta,eta,gamma,N)
        for l in range(1,q):
            suma=np.zeros(d)
            for s in range(l):
                for h in range(d):
                    suma[h]=suma[h]+A[l,s]*K[h,s]
            K[:,l]=func1(rho[l],Y[:,a]+dt*suma,mu,beta,eta,gamma,N)
        suma=np.zeros(d)
        for s in range(q):
            for h in range(d):
                suma[h]=suma[h]+b[s]*K[h,s]
        Y[:,a+1]=Y[:,a]+dt*suma
    return Y,tn
#Hemen hasten da aurrekoarekiko desberdintasun bakarra
# T bakoitzarekiko S-I lortzen da eta gero grafikatu
y1=(RK(n,T1,Y0,d,mu,beta,eta,gamma,N,q,A,b,tau)[0][1])
x1=(RK(n,T1,Y0,d,mu,beta,eta,gamma,N,q,A,b,tau)[0][0])
y2=(RK(n,T2,Y0,d,mu,beta,eta,gamma,N,q,A,b,tau)[0][1])
x2=(RK(n,T2,Y0,d,mu,beta,eta,gamma,N,q,A,b,tau)[0][0])
plt.plot(x1,y1,color='black',label='S-I T1')
plt.plot(x2,y2,color='red',label='S-I T2')
plt.legend()
plt.title('SIRI-EREDUA')
plt.xlabel('Sentiberak')
plt.ylabel('Infektatuak')
plt.show()
