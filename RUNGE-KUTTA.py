#Programa honen helburua Runge-Kuttaren ebazpena egitea izango da
#Lehenik eta behin libreriak importatu

import numpy as np
import matplotlib.pyplot as plt
import math
def RK(n,T,Y0,d,mu,beta,eta,gamma,N,q,A,b,tau):
    #diskretizazioa egin
    dt=T/n 
    tn=np.zeros(n+1)                    
    for i in range(n+1):
        tn[i]=i*dt  
    #Y zeros osaturiko bektore bat sortu, bukaerako emaitza hor sartuko da          
    Y=np.zeros((d,n+1))                 
    Y[:,0]=Y0[:]#hasierako egora sartu
    rho,K=np.zeros(q),np.zeros((d,q)) 
    #rho-n t-k sartzen joaten dira hurrengo for ekin.   
    for a in range(n):
        for j in range(q):
            rho[j]=tn[a]+dt*tau[j]#gero hau funtzioan t bezala erabiliko da.
        K[:,0]=func1(rho[0],Y[:,a],mu,beta,eta,gamma,N)#lehenengoa finkatu
        for l in range(1,q):
            suma=np.zeros(d)
            for s in range(l):
                for h in range(d):
                    suma[h]=suma[h]+A[l,s]*K[h,s]
            K[:,l]=func1(rho[l],Y[:,a]+dt*suma,mu,beta,eta,gamma,N)#Runge_kuttaren Ki lortzen da
        suma=np.zeros(d)
        for s in range(q):
            for h in range(d):
                suma[h]=suma[h]+b[s]*K[h,s]
        Y[:,a+1]=Y[:,a]+dt*suma
    return Y,tn
#Orain datuak adieraziko dira
#T hasierako egoera eta dimentsioak
T,d,Y0=7,3,np.array([3900000,100000,0])
#Runge-KUttaren datuak
q1,A1,b1,tau1=4,np.array([[0,0,0,0],[1/2,0,0,0],[0,1/2,0,0],[0,0,1,0]]),np.array([1/6,1/3,1/3,1/6]),np.array([0,1/2,1/2,1])
q,A,b,tau=4,np.array([[0,0,0,0],[1/2,0,0,0],[0,1,0,0],[0,0,1,0]]),np.array([1/6,1/3,1/3,1/6]),np.array([0,1/2,1/2,1])
n=100000
mu=4
beta=4
eta=3
gamma=2
N=4000000
# balio propioen kalkuluak print(((-2*mu+beta-gamma-eta)+(math.sqrt((-beta+gamma+eta)^2+4*beta*eta)))/2)
# Balio propioen kalkuluak print(((-2*mu+beta-gamma-eta)-(math.sqrt((-beta+gamma+eta)^2+4*beta*eta)))/2)

#SIRI ereduaren funtzioa definitzen da
#(S,I,R)=Z Bektorea
def func1(t,z,mu,beta,eta,gamma,N):
    return np.array([mu*N-mu*z[0]-beta*z[0]*z[1]/N,beta*z[0]*z[1]/N-z[1]*mu-gamma*z[1]+z[2]*eta,gamma*z[1]-z[2]*mu-z[2]*eta])

# Bi runge kutten emaitzak
x=(RK(n,T,Y0,d,mu,beta,eta,gamma,N,q,A,b,tau)[0][0])
y=(RK(n,T,Y0,d,mu,beta,eta,gamma,N,q,A,b,tau)[0][1])
z=(RK(n,T,Y0,d,mu,beta,eta,gamma,N,q,A,b,tau)[0][2])
t=(RK(n,T,Y0,d,mu,beta,eta,gamma,N,q,A,b,tau)[1])
x1=(RK(n,T,Y0,d,mu,beta,eta,gamma,N,q1,A1,b1,tau1)[0][0])
y1=(RK(n,T,Y0,d,mu,beta,eta,gamma,N,q1,A1,b1,tau1)[0][1])
z1=(RK(n,T,Y0,d,mu,beta,eta,gamma,N,q1,A1,b1,tau1)[0][2])
#Bukaerako S eta I ren adierazpena
print(x[100000])
print(y[100000])
#Bukaerako t lortzeko funtzioa
def bukaerakodenbora(x,y,z,t):
    i=0
    while i<100000 and abs(y[i]-y[i+1])>0.01 and abs(x[i]-x[i+1])>0.01 and abs(z[i]-z[i+1]):
        i+=1
    print(t[i])
#Bukaerako t
bukaerakodenbora(x,y,z,t)
#SIR ereduaren ebazpena
plt.plot(t,x,color='blue',label="Sentiberak")
plt.plot(t,z,color='red',label="Infektatuak")
plt.plot(t,z,color='black',label='Berreskuratuak')
plt.legend()
plt.title('SIRI-EREDUA')
plt.xlabel('Denbora')
plt.ylabel('Populazio kopurua')
plt.show()
#Bi runge kutta arteko desberdintasunaren grafika
plt.plot(t,abs(x-x1),color='blue',label='Sentiberen desberdintza')
plt.plot(t,abs(y-y1),color='red',label='Infektatuen desberdintza')
plt.plot(t,abs(z-z1),color='black',label='Berreskuratuen desberdintza')
plt.legend()
plt.title('BI RUNGE-KUTTA METODOEN DESBERDINTASUNA')
plt.xlabel('Denbora')
plt.ylabel('Populazio kopurua')
plt.show()
