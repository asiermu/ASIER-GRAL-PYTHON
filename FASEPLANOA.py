# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 15:33:01 2022

@author: asier
"""
#importazioak egin
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from fractions import Fraction
import random
from numpy import linalg as LA



#Programa honetan SIRI eredua baina jaoitzekin aztertuko da
#Aurreko programen gauza bera egingo da
#s eta i bektoreak sortu, s bektoreak s ardatzean dauden zenbakiak hartuko ditu eta irudia egingo du i lortzeko,
s=np.linspace(0,1,100)#s-ren diskretizazioa
i=np.zeros(100)
i2=np.zeros(100)
#S'=0 eta I'=0 direneko kurbak irudikatu 99 zenbaki harturik
def kurbenmozketa(mu,beta,eta,gamma):
    for a in range(1,100):
        i[a]=(mu*(1-s[a]))/(beta*s[a])#i isolatu
        i2[a]=(eta*(s[a]-1))/(beta*s[a]-mu-eta-gamma)#i isolatu
    plt.plot(s,i,color="green")#Sartzen baduzu, hau for aren barruan eremua koloreztatuko zenuke.
    plt.plot(s,i2,color="purple")
    plt.plot(s,1-s,color="red")
    sigma=beta/(mu+eta+gamma)
    plt.xlabel('Sentiberak')
    plt.ylabel('Infektatuak')
    a="beta="
    b=str(Fraction(str(beta)).limit_denominator())
    g=" ,"
    c="gamma="
    d=str(Fraction(str(gamma)).limit_denominator())
    p="eta="
    k=str(Fraction(str(eta)).limit_denominator())
    l=" datuen FASE-PLANOA"
    e=","
    j="mu="
    j1=str(Fraction(str(mu)).limit_denominator())
    plt.title(a+b+g+c+d+e+p+k+j+j1+l)
    plt.plot(1,0, marker="o", color="green")
    a=mu/((mu+eta)*sigma)
    b=((mu+eta)/beta)*(sigma-(mu/(mu+eta)))
    plt.plot(a,b, marker="o", color="black")
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.axhline(0)
    plt.axvline(0)

#Bi kasuak aztertuko dira

#Lehenengo kasua mu=0.1,beta=2,eta=3,gamma=4
kurbenmozketa(0.1,4,3,2)
plt.text(0.04, 0.1, "1", fontsize=10)
plt.text(0.5, 0.3, "2", fontsize=10)
plt.text(0.23, 0.7, "3", fontsize=10)
plt.text(0.01, 0.9, "4", fontsize=10)
plt.arrow(0.35, 0.3, -0.1, 0,
                 head_width = 0.02,
                  width = 0.005)
plt.arrow(0.35, 0.3, 0, 0.1,
                  head_width = 0.02,
                  width = 0.005)
plt.arrow(0.08, 0.1, 0, 0.1,
                  head_width = 0.02,
                  width = 0.005)
plt.arrow(0.08, 0.1, 0.1, 0,
                  head_width = 0.02,
                  width = 0.005)
plt.arrow(0.2, 0.77, -0.1, 0,
                 head_width = 0.02,
                  width = 0.005)
plt.arrow(0.2, 0.77, 0, -0.1,
                  head_width = 0.02,
                  width = 0.005)
plt.arrow(0.02, 0.77, 0.01, 0,
                 head_width = 0.001,
                  width = 0.0005)
plt.arrow(0.02, 0.77, 0, -0.1,
                  head_width = 0.01,
                  width = 0.005)
plt.show()    
show()

kurbenmozketa(9,4,3,2)
plt.text(0.1, 0.1, "1", fontsize=10)
plt.text(0.5, 0.4, "2", fontsize=10)
plt.arrow(0.2, 0.06, 0.1, 0,
                 head_width = 0.02,
                  width = 0.005)
plt.arrow(0.2, 0.06, 0, 0.1,
                  head_width = 0.02,
                  width = 0.005)
plt.arrow(0.4, 0.3, 0.1, 0,
                 head_width = 0.02,
                  width = 0.005)
plt.arrow(0.4, 0.3, 0, -0.1,
                  head_width = 0.02,
                  width = 0.005)


#Bi kasu daude bata bigarren puntu kritikora joango denean eta besta (1,0)puntu kritikora joango direnean emaitza denak
#Bifurkazio puntua, hurrengoa izango da:
#Bifurkazio=((mu+eta)*(beta))/(mu*(mu+eta+gamma)) Beraz distribu 1 baino handiagoa den edo ez sitemako puntu kritikoen egonkortasuna aldatzen du, horregaitik da bifurkazio
#Bifurkazioa>1 gaixotasuna ez da desageruko
#Argi dago lau konstateen arteko erlazio bat dela eta beta txikiagotzen bada, Bifurkazioa deituriko konstantea txikiago izango dela, beraz aukera gehiago gaixotasuna desagertzeko.
#Jaiotzak eta heriotzak berdinak izan harren giza talde desberdinetan gertatzen denez eragina dute gaixotasunaren eboluzioan.

#Orain aurrekoekin egin den gauza bera
#Fase_diagrama hutsa irudikatuko da koefizienteen arabera
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
    a="beta="
    b=str(Fraction(str(beta)).limit_denominator())
    g=" ,"
    c="gamma="
    d=str(Fraction(str(gamma)).limit_denominator())
    p="eta="
    k=str(Fraction(str(eta)).limit_denominator())
    l=" datuen FASE-PLANOA"
    e=","
    j="mu="
    j1=str(Fraction(str(mu)).limit_denominator())
    plt.title(a+b+g+c+d+e+p+k+j+j1+l)
    #Ardatzen izenak finkatzen dira
    xlabel('Sentiberak')
    ylabel('Infektatuak')
    #Oreka puntuak marraztu
    sigma=beta/(mu+eta+gamma)
    plt.plot(N,0, marker="o", color="green")
    plt.axhline(0,color="black")
    plt.axvline(0,color="black")
    a=mu*N/((mu+eta)*sigma)
    b=N*((mu+eta)/beta)*(sigma-(mu/(mu+eta)))
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
    import random
    random.seed("1")
    # Mu diagrama egiten da
    fase_diagrama(mu,beta,eta,gamma,N)
    for i in range(100):
        S0,I0=random.random()*N,random.random()*N
        if S0+I0<=N:
            x=np.array([S0,I0])
            (t,x)=eulerModif(f,0,100,x,1000,mu,beta,eta,gamma,N)
            plt.plot(x[0],x[1])


mu,beta,eta,gamma=0.1,4,3,2
fase_diagramaemaitza(mu,beta,eta,gamma,2090)
show()

def egonkortasuna(mu,beta,eta,gamma,N=1):
    sigma=beta/(mu+eta+gamma)
    p=mu/(mu+eta)
    a=mu*N/((mu+eta)*sigma)
    b=N*((mu+eta)/beta)*(sigma-(mu/(mu+eta)))
    if sigma<=p and p<=1:
        print("(1,0) puntu kritiko asintotikoki egonkorra")
    else:
        print(a,b)
        print("puntu kritiko asintotikoki egonkorra")
egonkortasuna(mu,beta,eta,gamma)
fase_diagramaemaitza(9,4,3,2,2090)
show()
egonkortasuna(9,4,3,2,2090)
#distribu>1 baino handiago ez da desagertuko gaixotasuna
#beta handitzean aukera gehiago distru>1 gaixotasuna ez desagertzeko aukera
#gamma handitzea aukera gehiago ditru<=1 gaixotasuan desagertzeko aukera  
#mu eta eta ren aldaketak ez daude hain garbi zer eragiten duten