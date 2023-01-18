# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 16:23:16 2022

@author: asier
"""


from sympy import var
from sympy.solvers import solve

x, y, beta, mu, gamma, eta = var('x y beta mu gamma eta')

f1 = x*y*beta-mu*y-gamma*y+eta*(1-x-y)
f2 = mu-mu*x-beta*x*y
f3=x+y-1
sols = solve((f1, f2),x,y)
sols2=solve((f2,f3),x,y)
print(sols2)