# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 15:16:42 2023

@author: asier
"""
#BALIO PROPIOEN KALKULUAK
import sympy as sp
from sympy.solvers import solve
mu = sp.Symbol('mu')
beta = sp.Symbol('beta')
gamma= sp.Symbol('gamma')
eta= sp.Symbol('eta')
sigma= sp.Symbol('sigma')
lamda= sp.Symbol('lamda')

#Lehenengo ekuazio karakteristikoa lortu


B=sp.Matrix([[-mu-(mu+eta)*(sigma-(mu/(mu+eta)))-lamda,-beta*mu/((mu+eta)*sigma)],[(mu+eta)*(sigma-(mu/(mu+eta)))-eta,beta*mu/((mu+eta)*sigma)-(mu+gamma+eta)-lamda]])
print(B.det())

#Balio propioak lortu
A=sp.Matrix([[-mu-(mu+eta)*(sigma-(mu/(mu+eta))),-beta*mu/((mu+eta)*sigma)],[(mu+eta)*(sigma-(mu/(mu+eta)))-eta,beta*mu/((mu+eta)*sigma)-(mu+gamma+eta)]])
M=A.eigenvals()
#Ahal baldin bada simplifikatu simplify erabiliz
la=(list(M.keys()))
#Bi balio propioak lortu dira l1 eta l2
l1=la[0]
l2=la[1]
A=(solve(l1,mu,beta,eta,gamma,sigma))
A2=(solve(l1,mu,beta,eta,gamma,sigma))