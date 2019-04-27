#
# Luke McCulloch
# L2 Norms
# SWE

import numpy as np
import matplotlib.pyplot as plt
#from sympy import *
from matplotlib.font_manager import FontProperties


HL2 = np.loadtxt('HL2_onehalf.dat')
uL2 = np.loadtxt('UL2_onehalf.dat')
uHL2 = np.loadtxt('UHL2_onehalf.dat')
dx = np.loadtxt('dx.dat')

n=len(dx)


plt.plot(np.log(dx),np.log(HL2),label='Log L2(H) vs Log dx')
font = FontProperties(stretch='condensed', size='small')
plt.legend(prop=font, borderpad=0.515)
plt.xlabel(" Log(dx) ")
plt.ylabel("Log(L2 H)")
plt.title("$L2_{norm}$ of $H$")
plt.grid()
plt.show()

m=(HL2[5]-HL2[0])/(dx[5]-dx[0])
print 'slope of the log(H)-log(dx) graph = {}'.format(m)


plt.plot(np.log(dx),np.log(uL2),label='Log L2(u) vs Log dx')
font = FontProperties(stretch='condensed', size='small')
plt.legend(prop=font, borderpad=0.515)
plt.xlabel(" Log(dx) ")
plt.ylabel("Log(L2 u)")
plt.title("$L2_{norm}$ of $u$")
plt.grid()
plt.show()



plt.plot(np.log(dx),np.log(uHL2),label='Log L2(uH) vs Log dx')
font = FontProperties(stretch='condensed', size='small')
plt.legend(prop=font, borderpad=0.515)
plt.xlabel(" Log(dx) ")
plt.ylabel("Log(L2 uH)")
plt.title("$L2_{norm}$ of $uH$")
plt.grid()
plt.show()
##
##plt.plot(dx,HL2,label='L2(H) vs dx')
##font = FontProperties(stretch='condensed', size='small')
##plt.legend(prop=font, borderpad=0.515)
##plt.xlabel(" Log(dx) ")
##plt.ylabel("Log(H)")
##plt.grid()
##plt.show()
##
##
##
##plt.plot(dx,uL2,label='L2(u) vs dx')
##font = FontProperties(stretch='condensed', size='small')
##plt.legend(prop=font, borderpad=0.515)
##plt.xlabel(" Log(dx) ")
##plt.ylabel("Log(u)")
##plt.grid()
##plt.show()

m=(uL2[5]-uL2[0])/(dx[5]-dx[0])
print 'slope of the log(u)-log(dx) graph = {}'.format(m)
