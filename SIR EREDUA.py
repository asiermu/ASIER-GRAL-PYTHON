# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 14:12:33 2022

@author: asier
"""
#Importatu libreria
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from pylab import *
from fractions import Fraction
import random
#Programa honek hiru helburu izango ditu
#Lehenengo helburua SIR-ereduaren fase-planoa egitea izango da.
#Bigarrena fase_diagrama lortzea izango da. Fase_diagrama hori beta eta gamma parametroen araberakoa izango da. Fase-diagrama honek fase-planoan egin dena berretsi beharko du.
# Hirugarrena SIR ereduaren ebazpen batzuk egitea euler metodo modifikatuarekin eta fase_diagramarekin batera irudikatzea eta ikustea bien arteko erlazioa ondo dagoela.
#Fase_diagrama hutsa irudikatuko da mu-ren arabera
def fase_planoa(beta,gamma):
    #S+I=1 marraztu  gorriz,S+I<=N=1 eremua identifikatzekko
    s=np.linspace(0,1,100)
    plt.plot(s,1-s,color="red")#S+I=1 zuzena gorriz marraztu
    #S=gamma/beta zuzena bertikala marraztu beltzez
    plt.axvline(x=gamma/beta,ymin=0, ymax=1,color="black")   
    #Bi fase-plano mota daude. 
    #Lehenengoan gamma/beta>=1, hau da R0<=1 denean. S'<0 eta I'<0 izango dena
    if gamma/beta>=1:
        plt.text(0.4, 0.4, "1", fontsize=10)#eremu bakarra
        plt.arrow(0.4, 0.3, -0.1, 0,
                  head_width = 0.02,
                  width = 0.005)
        plt.arrow(0.4, 0.3, 0, -0.1,
                  head_width = 0.02,
                  width = 0.005)
    #Bigarrena gamma/beta<1 hau da R0>1 denean eta kasu honetan bi eremu desberdin daude.
    else:
        #Arazoak gamma/beta oso txikia denean izan ere geziak ez dira kabitzen
        #gamma/beta zuzena eta S+I=1 zuzenak ebakitzen duten puntua lortu geziak eta eremuak hobeto irudikatzeko.
        import sympy as sym
        s=sym.Symbol('s')
        i=sym.Symbol('i')
        puntua=sym.solve([s+i-1,s-gamma/beta], dict=True)
        i=puntua[0].get(i)
        plt.text(gamma/beta-0.1,i-0.1, "1", fontsize=10)#Lehenengo eremua I'<0 eta S'<0 izango da.
        plt.arrow(gamma/beta-0.13,i-0.13, -0.1, 0,
                  head_width = 0.02,
                  width = 0.005)
        plt.arrow(gamma/beta-0.13,i-0.13, 0, -0.1,
                  head_width = 0.02,
                  width = 0.005)
        plt.text(gamma/beta+0.07,i-0.23, "2", fontsize=10)#Bigarren eremuan I'>0  eta S'<0 izango da.
        plt.arrow(gamma/beta+0.13,i-0.3, -0.1, 0,
                  head_width = 0.02,
                  width = 0.005)
        plt.arrow(gamma/beta+0.13,i-0.3, 0, +0.1,
                  head_width = 0.02,
                  width = 0.005)        
    #Grafikaren izena jarri.
    a=r"$\beta$="
    b=str(Fraction(str(beta)).limit_denominator())
    g=" eta "
    c=r"$\gamma$="
    d=str(Fraction(str(gamma)).limit_denominator())
    e="datuen FASE-PLANOA"
    title(a+b+g+c+d +e)
    plt.xlabel('Sentiberak')
    plt.ylabel('Infektatuak')
    #S eta I ardatzak finkatu
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.axhline(0)
    plt.axvline(0)
#Lehenengo helburua egin da. Lortu diren ondorio nagusiak hurrengoak dira:
#R0<=1 denean ez da epidemia egongo. Infektatuak beti behera eta (s,0) non s<=1 puntu kritiko guztiak asintotikoki egonkorrak izango dira.
#Ondoriozta daiteke gaixotasuna desagertuko dela, baina oraingoz ezin da jakin zenbat denboran eta sentiberen eta berreskuratuen tamaina.
#R0>1 s>gamma/beta direnak (s,0)puntu kritiko egongaitza izango da. Hasierako egoeran S0>s bada gamma/beta-ra heldu harte infektatu kopurua igoko da, maximora helduz. Beraz antzeman daiteke R0>1 baino bada gaixotasunaren desagerpena mantsoagoa izango dela, kurba gehiago eginik.
#R0 ren arabera gaixotasunaren jorrapena desberdina izango. Hau da mantsoago edo azkarrago desagertuko da gaixotasuna baina kasu denetan desagertuko da gaixotasuna. R0 eta hasierako egoeraren arabera, bukaerako sentibera eta berreskuratuen tamaina finkatuko dute.

#Fase_diagrama marraztuko da orain.
#N=1 suposatu da, beste N bat hartuz irudia proportzionala izango da.
def fase_diagrama(beta,gamma,N=1):
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
    Bs =-beta*S*I/N
    Bi =beta*S*I/N-gamma*I
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
    axis([smin, smax, imin, imax])#Nahiago izanez gero di eta ds kendu, bakarrik irudiaren egokierarako balio dute.
    #Titulua ezarri
    a=r"$\beta$="
    b=str(Fraction(str(beta)).limit_denominator())
    g=" eta "
    c=r"$\gamma$="
    d=str(Fraction(str(gamma)).limit_denominator())
    e="datuen FASE-DIAGRAMA"
    title(a+b+g+c+d +e)
    plt.xlabel('Sentiberak')
    plt.ylabel('Infektatuak')
    plt.axhline(0,color="black")
    plt.axvline(0,color="black")
#Hirugarren helburua;
#Orain emaitz batzuk sartu dira. Emaitzak kalkulatzeko Euler modifikatu metodoa erabiliko da.
def eulerModif(f,t0,tn,X0,n,beta,gamma,N=1):
    t=np.linspace(t0,tn,n+1)
    m=len(X0)
    x=np.zeros((m,n+1))
    x[:,0]=X0 #Lehenengo zutabea
    h=(tn-t0)/n
    for i in range(1,n+1):
        x[:,i]=x[:,i-1]+h*f(t[i-1]+h/2,x[:,i-1]+h/2 *f(t[i-1],x[:,i-1],beta,gamma,N),beta,gamma,N)
    return ((t,x))
#x bektorea [S,I] izango da
#Funtzioak definitzen dira sistema edukitzeko
def f1(t,x,beta,gamma,N=1):
    return(-beta*x[0]*x[1]/N)
def f2(t,x,beta,gamma,N):
    return(beta*x[0]*x[1]/N-gamma*x[1])
def f(t,x,beta,gamma,N=1):
    return(np.array([f1(t,x,beta,gamma,N),f2(t,x,beta,gamma,N)]))
def emaitza(beta,gamma,N=1):
    #Zoriozko hasierako balio batzuk emango dira eta euren emaitzak fase-diagraman sartuko dira, kontuak eduki S+I<=N dela.
    random.seed("ASIERMU")
    fase_diagrama(beta, gamma,N)
    #100 haiserako balio aukeratu dira, zoriz eta S+I<=N betetzen dutenak irudikatzen dira.
    for i in range(100):
        S0,I0=random.random()*N,random.random()*N
        if S0+I0<=N:
            x=np.array([S0,I0])
            (t,x)=eulerModif(f,0,100,x,1000,beta,gamma,N)
            S=x[0]
            I=x[1]
            plt.plot(S,I)
            plt.axvline(x=gamma*N/beta, ymin=0, ymax=N,color="purple")    
    R0=beta/gamma
    print(R0)
#fase_planoa(4/10,1/10)
#show()
#fase_planoa(1/10,1/10)
#show()
emaitza(4,2,4000000)
show()
#emaitza(1/10,1/10,4000000)
#show()
#R0<1 txikiago egonkorra sistema
#R0>1 baino handiagoa kurba gehiago beraz suposa dezakegu denbora gehiago tardatuko duela desagertzen baina ezin da jakin oingoz.Kurba handiagoa denez I max gehiago
#R0>1 I handiagotzen hasieran I0 eta S0 ren arabera 
#S0>1/R0 kurba egiten du 
#Egindako adibideak ondo ematen dute, fase-diagramak fase-planoan egindako berresten du eta emaitzak fase-diagraman zentsua daukate, gezien norabideak jarraitzen dute.