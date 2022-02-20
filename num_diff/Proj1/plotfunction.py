import numpy as np
from scipy import linalg as lg
import matplotlib.pyplot as plt



#x, y, xlabel, ylabel, axis scale, save bool
def plot(x, y, xlabel, ylabel, scale, save, savename = "1", plotlabel = None):
    try:
        k = x[0][1]     #Check if there are multiple arrays to be plotted
    except (TypeError, IndexError):
        if scale == "loglog":
            plt.loglog(x,y, label = plotlabel)
            plt.grid(True, which="both", axis='both')
        if scale == "semilogx":
            plt.grid(True, which='both', axis = 'both')
            plt.semilogx(x,y, label = plotlabel)
        if scale == "semilogy":
            plt.grid(True, which='both', axis = 'both')
            plt.semilogy(x,y, label = plotlabel)
        if scale == "plot":
            plt.grid(True)
            plt.plot(x,y, label = plotlabel)
    else:
        for i in range(len(x)):
            if scale[i] == "loglog":
                plt.loglog(x[i],y[i], label = plotlabel[i])
                plt.grid(True, which="both", axis='both')
            if scale[i] == "semilogx":
                plt.grid(True, which='both', axis = 'both')
                plt.semilogx(x[i],y[i], label = plotlabel[i])
            if scale == "semilogy":
                plt.grid(True, which='both', axis = 'both')
                plt.semilogy(x[i],y[i], label = plotlabel[i])
            if scale == "plot":
                plt.grid(True)
                plt.plot(x[i],y[i], label = plotlabel[i])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if plotlabel != None:
        plt.legend()
    if save == True:
        plt.savefig(savename + '.pdf')
    plt.show()




lotka(t, u) function
