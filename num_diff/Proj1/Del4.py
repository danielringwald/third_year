import numpy as np
from scipy import linalg as lg
from Del2 import adaptiveRK34
import matplotlib.pyplot as plt
import math
from scipy.integrate import solve_ivp
from plotfunction import plot

def calcvdp(my):
    def vdp(t, u):
        return np.array([u[1], my * (1 - u[0]**2) * u[1] - u[0]])
    return vdp

def del41():
    my = 100     #This is my in the van der Pol equation, task 4.1
    t0 = 0
    tf = 2*my
    y0 = [2, 0]
    tol = 1e-6

    adap = adaptiveRK34(calcvdp(my), t0, tf, y0, tol)
    x = np.zeros(len(adap[0]))
    y = np.zeros(len(adap[0]))
    for i in range(len(x)):
        x[i] = adap[1][i][0]
        y[i] = adap[1][i][1]

    #def plot(x, y, xlabel, ylabel, scale, save, savename = "1", plotlabel = None):
    plot(adap[0], y, "Time", "$y'_1$ = $y_2$", "plot", True, "41")
    #plot y2 against time
    # plt.plot(adap[0], x)
    # plt.show()

def del411():
    my = 100     #This is my in the van der Pol equation, task 4.1
    t0 = 0
    tf = 2*my
    y0 = [2, 0]
    tol = 1e-6

    adap = adaptiveRK34(calcvdp(my), t0, tf, y0, tol)
    x = np.zeros(len(adap[0]))
    y = np.zeros(len(adap[0]))
    for i in range(len(x)):
        x[i] = adap[1][i][0]
        y[i] = adap[1][i][1]

    #def plot(x, y, xlabel, ylabel, scale, save, savename = "1", plotlabel = None):
    #Phase portrait, y2 as a function of y1
    plt.title("$y_0$ = " + str(y0))
    plot(x, y, "$y_1$", "$y_2$", "plot", True, "411phase1")


def del42():
    mys = np.array([10, 15, 22, 33, 47, 68, 100, 150, 220, 330, 470])    #This is my in the van der Pol equation
    t0 = 0                                                               #Skip 1000 and 670
    tf = 0.7*mys
    y0 = [2, 0]
    tol = 1e-6

    adaps = [None]*len(mys)
    steps = np.array([None]*len(mys))
    for i in range(len(mys)):
        adaps[i] = adaptiveRK34(calcvdp(mys[i]), t0, tf[i], y0, tol)
        steps[i] = len(adaps[i][0]) - 1

    #Calculate slope, was roughly q = 2
    #dN/dmy = C*q*my^(q-1), stiffnes is prop. to my^(q-1), in this case my
    k = (math.log(steps[6])-math.log(steps[4]))/(math.log(mys[6])-math.log(mys[4])) #Slope
    plt.title("'adaptiveRK34' method")
    plot([mys, mys], [steps, mys**k * steps[-1]/(mys[-1]**k)], "$\mu$", "Amount of steps", ["loglog", "loglog"], True, "42", ["Numerical calculation", "Reference slope = "+str(round(k))])

    #OLD PLOT
    # plt.loglog(mys, steps)
    # plt.show()

def del43():
    y0 = [2, 0]
    mys = np.array([10, 15, 22, 33, 47, 68, 100, 150, 220, 330, 470, 680, 1000])
    t0 = 0                                                               #Skip 1000 and 670
    tf = 0.7*mys
    #sol = [None]*len(mys)
    steps = [None]*len(mys)
    for i in range(len(mys)):
        steps[i] = len(solve_ivp(calcvdp(mys[i]), [t0, tf[i]], y0, method='BDF').t) - 1
    plt.title("Internal solver: scipy.integrate.solve_IVP")
    plot(mys, steps, "$\mu$", "Amount of steps", "loglog", True, "43")

    #OLD PLOT
    # plt.loglog(mys, steps)
    # plt.show()

#Does it work for mu = 10000?? YES!
def del431():
    y0 = [2, 0]
    mys = 10000
    t0 = 0                                                               #Skip 1000 and 670
    tf = 2*mys
    #sol = [None]*len(mys)
    sol = solve_ivp(calcvdp(mys), [t0, tf], y0, method='BDF')
    y = sol.y[1]
    t = sol.t
    plt.title("Internal solver: scipy.integrate.solve_IVP")
    plot(t, y, "Time", "y", "plot", True, "solver10k")

    #OLD PLOT
    # plt.loglog(mys, steps)
    # plt.show()

del411()
