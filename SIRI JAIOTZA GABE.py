# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 09:30:42 2022

@author: asier
"""

#SIRI EREDUA JAIOTZA GABE
#Importazioak egin
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from pylab import *
from fractions import Fraction
import random


#Programa hau SIR ereduaren programaren ondorengoa da.
#SIRI eredua aztertuko da eta SIR programan jarraitutako pauso berdinak jarraituko dira.
#Kasu honetan helburuetako bat, aurreko SIR ereduarekiko konparaketa egitea da
#Fase_diagrama hutsa irudikatuko da parametroen arabera

#Fase-planoa egiterako orduan argi dago S'<0 dela, lehenengo ekuaziotik, beraz bakarrik bigarrena irudikatu behar da.

#S'=0 eta I'=0 direneko kurbak irudikatu 1000 zenbaki harturik horrela eremua marraztuta geldituko da.
def kurbenmozketa(beta,eta,gamma):
    s=np.linspace(0,1,1000)#s-ren diskretizazioa
    i2=np.zeros(1000)
    i1=np.zeros(1000)
    kont=0
    for a in range(1000):
        if beta*s[a]-gamma-eta<0:#Asintota dagoenenko  kasua konpondu, ze pythonek grafiko jarrai bezala irudikatzen du.
            i1[a]=(-eta+eta*s[a])/(beta*s[a]-gamma-eta)
            kont+=1
    plt.plot(s[:kont],i1[:kont],color="purple")
    for a in range(1000):
        if beta*s[a]-gamma-eta>0:
            i2[a]=(-eta+eta*s[a])/(beta*s[a]-gamma-eta)
    plt.plot(s[kont:],i2[kont:],color="purple")    
    plt.plot(s,1-s,color="red")
    a=r"$\beta$="
    b=str(Fraction(str(beta)).limit_denominator())
    g=" ,"
    c=r"$\gamma$="
    d=str(Fraction(str(gamma)).limit_denominator())
    p=r"$\eta$="
    k=str(Fraction(str(eta)).limit_denominator())
    l=" datuen FASE-PLANOA"
    e=","
    plt.title(a+b+g+c+d+e+p+k+l)
    plt.xlabel('Sentiberak')
    plt.ylabel('Infektatuak')
    #Puntu kritikoak irudikatu
    plt.plot(1,0, marker="o", color="green")
    plt.plot(0,eta/(eta+gamma), marker="o", color="black")
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.axhline(0)
    plt.axvline(0)
#Zailagoa da kasu honetan orokortzea beraz kasu konkretu baterako irudikatuko da.
kurbenmozketa(4,3,2)
plt.text(0.3, 0.4, "1", fontsize=10)#eremu bakarra
plt.arrow(0.4, 0.3, -0.1, 0,
          head_width = 0.02,
                  width = 0.005)
plt.arrow(0.4, 0.3, 0, 0.1,
                  head_width = 0.02,
                  width = 0.005)
plt.text(0.2, 0.75, "2", fontsize=10)#eremu bakarra
plt.arrow(0.2, 0.7, -0.1, 0,
          head_width = 0.02,
                  width = 0.005)
plt.arrow(0.2, 0.7, 0, -0.1,
                  head_width = 0.02,
                  width = 0.005)
show()
#Argi dago gaixotasuna ez dela desagertuko ze (1,0)da puntu kritiko bakarra I=0 baina beti doaz ezkerrerantz geziak zeren eta beti S'<0
#Beste puntu kritikoak (0,eta/(eta+gamma)), argi uzten du, horren 0<I<1 izango dela emaitza.
#eta eta gamma parametroen arabera, argi geratu da, berreskuratze tasa gero eta handiagoa denean bukaerako sentibera kopurua gero eta gutxiago.
#Bigarren pausua fase_diagrama marraztea da. SIR ereduan egin den gauza berdina egiten da.
def fase_diagrama(beta,eta,gamma,N=1):
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
    Bs = -beta*S*I/N
    Bi =beta*S*I/N-gamma*I+eta*(N-S-I)
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
    ds = (smax - smin)/(NS-1)
    di = (imax - imin)/(NI-1)
    axis([smin-ds, smax+ds, imin-di, imax+di])
    #Tituluan mu aldatzeko mu karaktere kate moduan idatziko da eta  mu zatiki moduan agertuko da tituluan.
    #Ardatzen izenak finkatzen dira
    a=r"$\beta$="
    b=str(Fraction(str(beta)).limit_denominator())
    g=","
    c=r"$\gamma$="
    d=str(Fraction(str(gamma)).limit_denominator())
    p=r"$\eta$="
    k=str(Fraction(str(eta)).limit_denominator())
    e=","
    l=" datuen FASE-DIAGRAMA"
    title(a+b+g+c+d+e+p+k+l)
    xlabel('Sentiberak')
    ylabel('Infektatuak')
    #Oreka puntuak marraztu
    plt.plot(N,0, marker="o", color="green")
    plt.plot(0,N*eta/(eta+gamma), marker="o", color="black")

#Orain emaitz batzuk sartu dira. Emaitzak kalkulatzeko Euler modifikatu metodoa erabiliko da.
def eulerModif(f,t0,tn,X0,n,beta,eta,gamma,N=1):
    t=np.linspace(t0,tn,n+1)
    m=len(X0)
    x=np.zeros((m,n+1))
    x[:,0]=X0 #Lehenengo zutabea
    h=(tn-t0)/n
    for i in range(1,n+1):
        x[:,i]=x[:,i-1]+h*f(t[i-1]+h/2,x[:,i-1]+h/2 *f(t[i-1],x[:,i-1],beta,eta,gamma,N),beta,eta,gamma,N)
    return ((t,x))
#x bektorea [S,I] izango da
#Funtzioak definitzen dira sistema edukitzeko
def f1(t,x,beta,eta,gamma,N=1):
    return(-beta*x[0]*x[1]/N)
def f2(t,x,beta,eta,gamma,N=1):
    return(beta*x[0]*x[1]/N-gamma*x[1]+eta*(N-x[0]-x[1]))
def f(t,x,beta,eta,gamma,N=1):
    return(np.array([f1(t,x,beta,eta,gamma,N),f2(t,x,beta,eta,gamma,N)]))
def emaitza(beta,eta,gamma,N=1):
    #Zoriozko hasierako balio batzuk emango dira eta euren emaitzak fase-diagraman sartuko dira, kontuak eduki S+I<=1 dela.
    random.seed("Asier Muñoz")
    #Diagrama egiten da
    fase_diagrama(beta,eta,gamma,N)
    for i in range(100):
        S0,I0=random.random()*N,random.random()*N
        if S0+I0<=N:
            x=np.array([S0,I0])
            (t,x)=eulerModif(f,0,100,x,1000,beta,eta,gamma,N)
            plt.plot(x[0],x[1])
emaitza(4,3,2,4)
#Argi geratu da sistema ez dela iñoiz desgartuko
#Beti bukatu da gaixotasuna(0,eta/(eta+gamma)) puntuan.
#Argi dago ze puntuta bukatuko duen eta berdin du bere hasierako egoera.