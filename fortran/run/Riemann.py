#
# Luke McCulloch
# 'exact' 
# Solution to the Riemann Problem
# for the Shallow Water Equations
"""
This exact solution to the Riemann problem 
is found by Newton iteration

"""
import numpy as np
import matplotlib.pyplot as plt

"""
from sympy import *
# sympy way, but probably needs a good starting guess?
HL, HR, Hm, g = symbols('HL HR Hm g')
print solve((2.*(sqrt(g*HL) - sqrt(g*Hm))-(Hm-HR)*sqrt((g/2.)*((1./Hm)-(1./HR)))),Hm)
"""

def f(g,HL,Hm,HR):
    print 'inside the f function,'
    print 'Hm = {}, HR = {}'.format(Hm, HR)
    print 'fraction thing = {}'.format((1./Hm)-(1./HR))
    return 2.*(np.sqrt(g*HL) - np.sqrt(g*Hm)) - (Hm-HR) * np.sqrt((g/2.)*((1./Hm)+(1./HR)))

def df(g,HL,Hm,HR):
    dum1 = -((g*Hm)**(-.5))
    dum2 = -np.sqrt((g/2.)*((1./Hm)+(1./HR)))
    dum3 = -(Hm-HR)/2.
    dum4 = ((g/2.)*((1./Hm)+(1./HR)))**(-.5)
    dum5 = ((-g/2.)*(1./(Hm**2.)))
    
    return dum1+dum2+dum3*dum4*dum5

def df_easy(g,HL,Hm,HR,dx=.00000100001):
    """
    #df/dHm
    """
    d1 = (f(g,HL,Hm+dx,HR) - f(g,HL,Hm,HR) ) / dx
    return d1



    
HL      = 1.0
HRoHL   = 0.5
HR      = HRoHL*HL
g       = 9.81

Hmi     = .01
Hm=.4
tol=.0000000001
delta=1.
count = 0
dummy1=1.0




while(abs(dummy1)>tol):

    dummy1 = f(g,HL,Hmi,HR)
    print 'f = {}'.format(dummy1)
    dummy2 = df(g,HL,Hmi,HR)
         
    Hm      = Hmi - dummy1/dummy2
    delta   = Hm-Hmi
    Hmi=Hm
    count = count+1

print 'Newtons method complete'   
print 'Hm = {}'.format(Hm)
print 'f = {}'.format(f(g,HL,Hmi,HR))
print 'df = {}'.format(df(g,HL,Hmi,HR))
print 'count = {}'.format(count)
