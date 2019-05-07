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


def lamda(g,U,H):
    return U - np.sqrt(g*H)

def lamda1L(g,UL,HL):
    return UL - np.sqrt(g*HL)

def lamda1m(g,UM,HM):
    return UM - np.sqrt(g*HM)

def lamda1R(g,UL,HL):
    return UL + np.sqrt(g*HL)

def A(g,UL,HL):
    return UL+2.*np.sqrt(g*HL)

def S(HL,HM,uHL,uHM):
    return (uHL - uHM)/(HL-HM)



fp = open("fifi", 'r')
for i, line in enumerate(fp):
    if i == 4:
        print line
        C=float(line)
        print C
    elif i > 4:
        break
fp.close()




HL      = 1.0
HRoHL   = float(np.loadtxt('out.dat'))
HR      = HRoHL*HL
g       = 9.81

# assign a starting value to Hm for newton's method
Hmi     = .01
Hm=.4
tol=.000000000000005
delta=1.
count = 0
dummy1=1.0

ini = -5.
fini=5.


#Load H
solvH   = np.loadtxt("H.dat")
npts=len(solvH)
nsteps = 100
dx=(fini-ini)/npts
dt=.001


#parameters=np.linspace(ini,fini,npts,endpoint=True)

#load x
parameters=np.loadtxt("x_pos.dat")

#HL, HR, Hm, g = symbols('HL HR Hm g')
#print solve((2.*(sqrt(g*HL) - sqrt(g*Hm))-(Hm-HR)*sqrt((g/2.)*((1./Hm)-(1./HR)))),Hm)

# Newton's method to compute Hm:
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

#parameters=np.linspace(ini,fini,npts,endpoint=True)
#Load t
t=float(np.loadtxt('time.dat'))
nt=t+dt*nsteps
UL = 0.
UR = 0.
Rwave = A(g,UL,HL)
#shockspeed = ((UL*HL)-(um*Hm))/(HL-Hm)
shockspeed = ((um*Hm)-(UR*HR))/(Hm-HR)
#shockspeed = ((um*Hm)-(UR*HR))/(Hm*HR)
#shockspeed = ((um*Hm)-(UR*HR))/(Hm*HR)

#print Rwave
#for n in range(nsteps):
    #print n
    #t=t+dt*n
n=1
#t   = float(np.loadtxt("time.dat"))


    #print 't= {}'.format(t)
Levt = lamda(g,UL,HL)*t
Mevt = lamda(g,um,Hm)*t
Rev  = lamda1R(g,UL,HR)
st   = shockspeed*t
for j,s in enumerate(parameters):
#print j,s # j is the index, s is the value
    
    if s<=Levt:
        H[j,n] = HL
        u[j,n] = UL
    elif Levt<s<=Mevt:
        H[j,n] = (1./(9.*g))*(((Rwave) - (s/t))**2)
        Htilda = H[j,n]
        u[j,n] = Rwave - 2.*np.sqrt(g*Htilda)
        #print 'u{},{}={}'.format(j,n,u[j,n])
        utilda = u[j,n]
    elif Mevt<s<=st:
        H[j,n] = Hm
        u[j,n] = um
    elif st<s:
        H[j,n] = HR
        u[j,n] = 0.#UR
    else:
        print 'Error off the grid'
        break
    uH[j,n]=u[j,n]*H[j,n]
    
#print 'u = {} and H = {}'.format(u[:,n],H[:,n],uH[:,n])

#for t in range(0,100,10):
#plt.plot(parameters,H[:,t],label='H(t={})'.format(t))

plt.plot(parameters,solvH[:],'x',label='Numeric H(t={})'.format(t))
plt.plot(parameters,H[:,n],label='Analytic H(t={})'.format(t))
font = FontProperties(stretch='condensed', size='small')
plt.legend(prop=font, borderpad=0.515)
plt.xlabel(" x position, m")
plt.ylabel("H")
plt.title("H vs X position with {} cells, Courant # = {}".format(npts, C))
plt.axis(xmin=-5., ymin= .0, xmax=5.0, ymax=1.1)
plt.grid()
#plt.show()

# Compute L2 norm for H:
H_L2 = 0.
for i in range(len(solvH)):
    H_L2 = H_L2 + ((H[i,n]-solvH[i])**2)*dx
H_L2 = np.sqrt(H_L2)
print 'L2 norm of H = {}'.format(H_L2)


#Load U
solvu   = np.loadtxt("u.dat")
plt.plot(parameters,solvu[:],'x',label='Numeric u(t={})'.format(t))
plt.plot(parameters,u[:,n],label='Analytic u(t={})'.format(t))
font = FontProperties(stretch='condensed', size='small')
plt.legend(prop=font, borderpad=0.515)
plt.xlabel(" x position, m")
plt.ylabel("u, m/s")
plt.title("u vs X position with {} cells, Courant # = {}".format(npts, C))
plt.axis(xmin=-5.1, ymin= -.1,ymax=4.1)
plt.grid()
#plt.show()

# Compute L2 norm for u:
u_L2 = 0.
for i in range(len(solvu)):
    u_L2 = u_L2 + ((solvu[i]-u[i,n])**2)*dx
u_L2 = np.sqrt(u_L2)
print 'L2 norm of u = {}'.format(u_L2)



#Load uH
solvuH   = np.loadtxt("uH.dat")
plt.plot(parameters,solvuH[:],label='Numeric uH(t={})'.format(t))
plt.plot(parameters,uH[:,n],label='Analytic uH(t={})'.format(t))
font = FontProperties(stretch='condensed', size='small')
plt.legend(prop=font, borderpad=0.515)
plt.xlabel(" x position, m")
plt.ylabel("Momentum, uH")
plt.title("uH vs X position with {} cells, Courant # = {}".format(npts, C))
plt.axis(xmin=-5.1, ymin= -.1,ymax=1.5)
plt.grid()
plt.show()

# Compute L2 norm for uH:
uH_L2 = 0.
for i in range(len(solvuH)):
    uH_L2 = uH_L2 + ((solvuH[i]-uH[i,n])**2)*dx
uH_L2 = np.sqrt(uH_L2)
print 'L2 norm of uH = {}'.format(uH_L2)


#Graph Sonic Point and Shock
#plt.plot(parameters,solvu[:],label='Numeric u(t={})'.format(t))
plt.plot(parameters,(solvu[:]/np.sqrt(1.*g*solvH[:])),label='Numeric Mach(t={})'.format(t))
font = FontProperties(stretch='condensed', size='small')
plt.legend(prop=font, borderpad=0.515)
plt.xlabel(" x position, m")
plt.ylabel("Mach Number")
plt.title("Mach Number vs X position with {} cells, Courant # = {}".format(npts, C))
plt.axis(xmin=-5.1, ymin= -.5, xmax = 5.1,ymax=2.1)
plt.grid()
plt.show()

vbl =  max(abs(Levt/t),abs(Mevt/t))

file = open('ev.dat', 'w+')

file.write( "Max Lamda\n" )
file.write("{0:.16f}".format(abs(vbl)))
file.write( "\nEnd of File\n" )
file.close()

#numWavespeed   = np.loadtxt("LeftWaveSpeed.dat")
nws  = np.loadtxt("LeftWaveSpeed.dat")\

print 'Comparison of Left wave speeds'
print 'The Analytic Left EV = {}, the Numeric Max Left EV = {}'.format(Levt/t,nws[0])
print ''
print 'Numeric - Analytic = {}'.format(nws[0]-Levt/t)
print ''
print 'The Left solution begins to slow down at {}, face {}'.format(nws[1],nws[2])
print ''
#print 'Analytic H({},~1.0)={}, Numeric H({},~1.0)={}'.format(nws[2],H[nws[2],n],nws[2],solvH[nws[2]])


"""
file = open('HL2_onehalf.dat', 'a')
file.write("{0:.16f}\n".format(H_L2))
file.close()
file = open('UL2_onehalf.dat', 'a')
file.write("{0:.16f}\n".format(u_L2))
file.close()
file = open('dx.dat', 'a')
file.write("{0:.16f}\n".format(dx))
file.close()

file = open('uHL2_onehalf.dat', 'a')
file.write("{0:.16f}\n".format(uH_L2))
file.close()
"""
vbl=0.5
inih=[]
for j,s in enumerate(parameters):
	if s<=0.:
		inih.append(1.0)
	elif s>0.:
		inih.append(vbl)
inih=np.asarray(inih)

plt.plot(parameters,inih[:])
font = FontProperties(stretch='condensed', size='small')
#plt.legend(prop=font, borderpad=0.515)
plt.xlabel(" x position, m")
plt.ylabel("H")
plt.title('$H$ at $t=0$')
plt.axis(xmin=-5.1, ymin= -.1, xmax = 5.5,ymax=5.5)
plt.grid()
#plt.show()
