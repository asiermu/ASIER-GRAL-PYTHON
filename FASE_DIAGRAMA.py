# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 17:15:14 2022

@author: asier
"""
#Lan honetan fase-diagramak grafikatuko dira
#Libreria importatu
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from pylab import *
from fractions import Fraction
import random
#Fase_diagrama hutsa irudikatuko da mu-ren arabera
def fase_diagrama(mu,beta,eta,gamma,N=1):
    # Lan espazioa definitzen da, dominioa
    # Minimoak eta maximoak finkatzen dira
    smax = N
    smin = 0
    imax = N
    imin = 0
    # Ardatzaren tarte kopurua finkatzen dira.
    NS = 15
    NI = 15
    #Ardatz bakoitzeko diskretizazioa sortu.
    s =np.linspace(smin,smax,NS)
    i =np.linspace(imin,imax,NI)
    #Bektore koordenatuetatik, matrize koordenatuetara pasatu meshgrid-ekin.
    S, I = np.meshgrid(s, i)
    #Laneko sisteman S' eta I' ren zeinuak aztertu behar dira horretarako Bs eta Bi definitzen dira.
    Bs = mu*(N-S)-beta*S*I/N
    Bi =beta*S*I/N-I*mu-gamma*I+eta*(N-S-I)
    #B bektore bezala jarrita Bs eta Bi eraldatzen dira, geziak hobeto ikusteko.
    B=np.array([Bs,Bi])
    Bsbis=Bs/(1+LA.norm(B))
    Bibis=Bi/(1+LA.norm(B))
    #Figuraren tamaina finkatu
    figure(figsize=[10,10])
    #Gezien kokalekua, norabidea eta kolorea jartzen dira.
    QP = quiver(S,I, Bsbis, Bibis, color='r')
    quiverkey(QP, 0.85, 1.05, 1, '1 mT', labelpos='N')
    #Ardatzen propietateak ezarri
    #ds = (smax - smin)/(NS-1)
    #di = (imax - imin)/(NI-1)
    axis([smin,smax, imin, imax])#imin-di,smin-ds,imax+di smax+ds
    #Tituluan mu aldatzeko mu karaktere kate moduan idatziko da eta  mu zatiki moduan agertuko da tituluan.
    a="mu=" 
    b=str(Fraction(str(mu)).limit_denominator())
    c=" datuko Norabide diagrama"
    title(a+b+c)
    #Ardatzen izenak finkatzen dira
    xlabel('Sentiberak')
    ylabel('Infektatuak')
    #Oreka puntuak marraztu
    sigma=beta/(mu+eta+gamma)
    plt.plot(N,0, marker="o", color="green")
    plt.axhline(0,color="black")
    plt.axvline(0,color="black")
    a=mu/((mu+eta)*sigma)
    b=((mu+eta)/beta)*(sigma-(mu/(mu+eta)))
    print(a,b)
    plt.plot(a,b, marker="o", color="black")
#Laneko bi kasu desberdinak erakutsi
#fase_diagrama(1/(52*50))
#fase_diagrama(0.3)
#Orain emaitz batzuk sartu dira. Emaitzak kalkulatzeko Euler modifikatu metodoa erabiliko da.
def eulerModif(f,t0,tn,X0,n,mu,beta,eta,gamma,N=1):
    t=np.linspace(t0,tn,n+1)
    m=len(X0)
    x=np.zeros((m,n+1))
    x[:,0]=X0 #Lehenengo zutabea
    h=(tn-t0)/n
    for i in range(1,n+1):
        x[:,i]=x[:,i-1]+h*f(t[i-1]+h/2,x[:,i-1]+h/2 *f(t[i-1],x[:,i-1],mu,beta,eta,gamma,N),mu,beta,eta,gamma,N)
    return ((t,x))
#x bektorea [S,I] izango da
#Funtzioak definitzen dira sistema edukitzeko
def f1(t,x,mu,beta,eta,gamma,N=1):
    return(mu*(N-x[0])-beta*x[0]*x[1]/N)
def f2(t,x,mu,beta,eta,gamma,N=1):
    return(beta*x[0]*x[1]/N-x[1]*(mu)-gamma*x[1]+eta*(N-x[0]-x[1]))
def f(t,x,mu,beta,eta,gamma,N=1):
    return(np.array([f1(t,x,mu,beta,eta,gamma,N),f2(t,x,mu,beta,eta,gamma,N)]))
def fase_diagramaemaitza(mu,beta,eta,gamma,N=1):
    #Zoriozko hasierako balio batzuk emango dira eta euren emaitzak fase-diagraman sartuko dira, kontuak eduki S+I<=1 dela.
    random.seed("AsierMuoz")
    # Mu diagrama egiten da
    fase_diagrama(mu,beta,eta,gamma,N)
    for i in range(100):
        S0,I0=random.random()*N,random.random()*N
        if S0+I0<=N:
            x=np.array([S0,I0])
            (t,x)=eulerModif(f,0,100,x,1000,mu,beta,eta,gamma,N)
            plt.plot(x[0],x[1])

#fase_diagramaemaitza(0.3,0.5,0.3,0.5)
#show()
#fase_diagramaemaitza(1/(50*52),0.5,0.3,0.5)
#show()
#Zoriz aukeratzen dira lau koefizienteak
#fase diagrama, puntu kritikoen egonkortasunarekin bat datorren ikusi
#random.seed("Asier Muñoz 3")   
#mu=random.random()
#beta=random.random()
#eta=random.random()
#gamma=random.random()
mu,beta,eta,gamma=0.1,4,3,2
fase_diagramaemaitza(mu,beta,eta,gamma,2090)
show()
distribu=((mu+eta)*(beta))/(mu*(mu+eta+gamma))
print(distribu)
def egonkortasuna(mu,beta,eta,gamma):
    sigma=beta/(mu+eta+gamma)
    p=mu/(mu+eta)
    if sigma<=p and p<=1:
        print(mu,beta,eta,gamma)
        print("(1,0) puntu kritiko asintotikoki egonkorra")
    else:
        print("(0.0020, 0.3745) puntu kritiko asintotikoki egonkorra")
        print(mu,beta,eta,gamma)
egonkortasuna(mu,beta,eta,gamma)
fase_diagramaemaitza(9,4,3,2)
#distribu>1 baino handiago ez da desagertuko gaixotasuna
#beta handitzean aukera gehiago distru>1 gaixotasuna ez desagertzeko aukera
#gamma handitzea aukera gehiago ditru<=1 gaixotasuan desagertzeko aukera  
#mu,eta ez dago argi egoeraren arabera
