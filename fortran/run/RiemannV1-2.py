#
# Luke McCulloch
# Solution to the Riemann Problem
# for the Shllow Water Equations

import numpy as np
import matplotlib.pyplot as plt
#from sympy import *
from matplotlib.font_manager import FontProperties


def f(g,HL,Hm,HR):
    #print 'inside the f function,'
    #print 'Hm = {}, HR = {}'.format(Hm, HR)
    #print 'fraction thing = {}'.format((1./Hm)-(1./HR))
    return 2.*(np.sqrt(g*HL) - np.sqrt(g*Hm))-(Hm-HR)*np.sqrt((g/2.)*((1./Hm)+(1./HR)))


def df(g,HL,Hm,HR):
    dum1 = -((g*Hm)**(-.5))
    dum2 = -np.sqrt((g/2.)*((1./Hm)+(1./HR)))
    dum3 = -(Hm-HR)/2.
    dum4 = ((g/2.)*((1./Hm)+(1./HR)))**(-.5)
    dum5 = ((-g/2.)*(1./(Hm**2.)))
    
    return dum1+dum2+dum3*dum4*dum5

def lamda1L(g,UL,HL):
    return UL - np.sqrt(g*HL)

def lamda1m(g,UM,HM):
    return UM - np.sqrt(g*HM)

def A(g,UL,HL):
    return UL+2.*np.sqrt(g*HL)

def S(HL,HM,uHL,uHM):
    return (uHL - uHM)/(HL-HM)





HL      = 1.0
HRoHL   = 0.5
HR      = HRoHL*HL
g       = 9.81

Hmi     = .01
Hm=.4
tol=.000000000000005
delta=1.
count = 0
dummy1=1.0

ini = -2.5
fini=2.5
npts=50
nsteps = 100
dx=(fini-ini)/npts
dt=.001


#HL, HR, Hm, g = symbols('HL HR Hm g')
#print solve((2.*(sqrt(g*HL) - sqrt(g*Hm))-(Hm-HR)*sqrt((g/2.)*((1./Hm)-(1./HR)))),Hm)

while(abs(dummy1)>tol):

    dummy1 = f(g,HL,Hmi,HR)
    dummy2 = df(g,HL,Hmi,HR)
         
    Hm      = Hmi - dummy1/dummy2
    delta   = Hm-Hmi
    Hmi=Hm
    count = count+1
dummy1 = f(g,HL,Hm,HR)
print 'Newtons method complete'   
#print 'count = {}, Hm = {}, f = {}'.format(count, Hm, dummy1)


#compute um:
# eq. 13.57
um = 2.*(np.sqrt(g*HL) - np.sqrt(g*Hm))
print 'count = {}, Hm = {}, f = {}, um = {}'.format(count, Hm, dummy1, um)


H=np.zeros((npts,nsteps),float)
u=np.zeros((npts,nsteps),float)
uH=np.zeros((npts,nsteps),float)

parameters=np.linspace(ini,fini,npts,endpoint=True)
t=1
nt=t+dt*nsteps
UL = 0
UR = 0
Rwave = A(g,UL,HL)
shockspeed = ((UL*HL)-(um*Hm))/(HL-Hm)
#print Rwave
#for n in range(nsteps):
    #print n
    #t=t+dt*n
n=1
t   = np.loadtxt("time.dat")
    #print 't= {}'.format(t)
Levt = lamda1L(g,UL,HL)*t
Mevt = lamda1m(g,um,Hm)*t
st   = shockspeed*t
for j,s in enumerate(parameters):
#print j,s # j is the index, s is the value
    
    if s<=Levt:
        H[j,n] = HL
        u[j,n] = 0.
    elif Levt<=s<=Mevt:
        H[j,n] = (1./(9.*g))*(((Rwave) - (s/t))**2)
        Htilda = H[j,n]
        u[j,n] = Rwave - 2.*np.sqrt(g*Htilda)
        #print 'u{},{}={}'.format(j,n,u[j,n])
        utilda = u[j,n]
    elif Mevt<=s<=st:
        H[j,n] = Hm
        u[j,n] = um
    elif st<=s:
        H[j,n] = HR
        u[j,n] = UR
    else:
        print 'Error off the grid'
        break
print 'u = {} and H = {}'.format(u[:,n],H[:,n])

#for t in range(0,100,10):
#plt.plot(parameters,H[:,t],label='H(t={})'.format(t))
solvH   = np.loadtxt("H.dat")
plt.plot(parameters,solvH[:],label='solver H(t={})'.format(t))
plt.plot(parameters,H[:,n],label='H(t={})'.format(t))
font = FontProperties(stretch='condensed', size='small')
plt.legend(prop=font, borderpad=0.015)
plt.xlabel(" x position, m")
plt.ylabel("H")
plt.title("H vs X position with {} cells".format(npts))
plt.axis(xmin=-2.5, ymin= 0.)
plt.grid()
plt.show()



vbl =  max(abs(Levt),abs(Mevt))

file = open('ev.dat', 'w+')

file.write( "Max Lamda\n" )
#fo.write( '%.4f' % (max(abs(Levt),abs(Mevt) ))
file.write("{0:.16f}".format(vbl))
file.write( "\nEnd of File\n" )
file.close()


