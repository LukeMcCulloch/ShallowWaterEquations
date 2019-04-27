# Luke McCulloch
# Function to plot My Flow Solver Output
# October 2012

import numpy as np
import scipy as scipy
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


##class plot:
##
##    def __init__(self, label=""):
##
##        #constant
##        #self.maxexp = 1.e12
##        self.label = label
##        self.A = 0
##        self.B = 0
##
##    def plot(self):
##
##        plt.plot(self.omega, self.S, label=self.label)
##
##
##    def labelplot(self):
##
##        plt.grid()
##        font = FontProperties(stretch='condensed', size='x-small')
##        plt.legend(prop=font, borderpad=0.015)
##        plt.xlabel(r'x position $\x$  [m]')
##        plt.ylabel(r'Temperature $S_(\t)$  [C$^o\$]')
##        plt.title("Temperature vs Position")

def main():

    # Complex files
    fp = open("out.dat","r")
    lines=[x.strip() for x in fp if x.strip() and x.strip()[0] != "#" and x.strip()[0] !=0 ]
    fp.close()

    for line in lines:
        #print line
        lsplit = line.split()
        if lsplit[0]  == "ncells":
            ncells = int(lsplit[2])

    # Simple Files
    x       = np.loadtxt("x-pos.dat")
    t_n_40  = np.loadtxt("t_n_40.dat")
    t_n_80  = np.loadtxt("t_n_80.dat")
    t_n_120 = np.loadtxt("t_n_120.dat")
    
    tex40 = np.loadtxt("exact40.dat")
    tex80 = np.loadtxt("exact80.dat")
    tex120 = np.loadtxt("exact120.dat")
    



    print 'Number of cells, {}'.format(ncells)

    plt.plot(x,tex40,label='Exact 40 sec')
    plt.plot(x,t_n_40,label='40 sec')

    plt.plot(x,tex80,label='Exact 80 sec')
    plt.plot(x,t_n_80,label='80 sec')
    
    plt.plot(x,tex120,label='Exact 120 sec')
    plt.plot(x,t_n_120,label='120 sec')



    #plt.plot(x,tex40,marker='.', linestyle='-',label='Exact 40 sec')
    #plt.plot(x,t_n_40,marker='+', linestyle='--',label='40 sec')

    #plt.plot(x,tex80,marker='o', linestyle='-',label='Exact 80 sec')
    #plt.plot(x,t_n_80,marker='*', linestyle='-',label='80 sec')
    
    #plt.plot(x,tex120,marker='s', linestyle='-',label='Exact 120 sec')
    #plt.plot(x,t_n_120,marker='^', linestyle='-',label='120 sec')
    
    font = FontProperties(stretch='condensed', size='small')
    plt.legend(prop=font, borderpad=0.515)
    plt.xlabel(" x position, m")
    plt.ylabel("Temperature")
    plt.title("Temperature vs X position with {} cells".format(ncells))
    plt.axis(xmin=0., ymin= 0.)
    plt.grid()
    plt.show()






    plt.plot(x,t_n_40,label='40 sec')
    plt.plot(x,t_n_80,label='80 sec')
    plt.plot(x,t_n_120,label='120 sec')
    font = FontProperties(stretch='condensed', size='small')
    plt.legend(prop=font, borderpad=0.015)
    plt.xlabel(" x position, m")
    plt.ylabel("Temperature")
    plt.title("Temperature vs X position with {} cells".format(ncells))
    plt.axis(xmin=0., ymin= 0.)
    plt.grid()
    #plt.show()



if __name__ == '__main__':
    main()
