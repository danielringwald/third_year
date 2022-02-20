import numpy as np
from scipy import linalg as lg
import matplotlib.pyplot as plt
from Del2 import adaptiveRK34
import math
from plotfunction import plot   #plot(x, y, xlabel, ylabel, scale, save, savename = "1", plotlabel = " ")

class var:
    a = 3
    b = 9
    c = 15
    d = 15

def lotka(t, u):
    return np.array([var.a*u[0]-var.b*u[0]*u[1], var.c*u[0]*u[1]-var.d*u[1]])

def H(u):
    return var.c * u[0] + var.b * u[1] - var.d * math.log(u[0]) - var.a * math.log(u[1])

def del311():           #plot x and y against t
    u0 = [0.8, 1.3]     #Period changes with different start-values
    tol = 1e-11
    t0 = 0
    tf = 11
    adap = adaptiveRK34(lotka, t0, tf, u0, tol)
    x = np.zeros(len(adap[0]))
    y = np.zeros(len(adap[0]))
    for i in range(len(x)):
        x[i] = adap[1][i][0]
        y[i] = adap[1][i][1]
    plt.xlim(t0, tf)
    plt.title("$y_04$ = " + str(u0))
    plot([adap[0], adap[0]], [x, y], "Time", "Amount of specified animal","plot" , True, "311v2", ["Rabbits (x(t0) = "+str(u0[0])+")", "Foxes (y(t0) = "+str(u0[1])+")"])

def del312():       #plot y against x
    u0 = [0.9, 1.56]     #Period changes with different start-values
    tol = 1e-6
    t0 = 0
    tf = 10
    adap = adaptiveRK34(lotka, t0, tf, u0, tol)
    x = np.zeros(len(adap[0]))
    y = np.zeros(len(adap[0]))
    for i in range(len(x)):
        x[i] = adap[1][i][0]
        y[i] = adap[1][i][1]
    plot(x, y, "Rabbits", "Foxes", "plot", True, "312v2")

def del313():       #plot h against t
    u0 = [1.3, 0.8]     #Period changes with different start-values
    tol = 1e-8
    t0 = 0
    tf = 1000
    adap = adaptiveRK34(lotka, t0, tf, u0, tol)
    Hvec = np.zeros(len(adap[0]))
    for i in range(len(adap[0])):
        Hvec[i] = H(adap[1][i])
    plot(adap[0], abs(Hvec/H(u0) - 1), "Time", "Deviation", "semilogx", True, "HVSt")
    # plt.semilogx(adap[0], abs(Hvec/H(u0)-1))
    # plt.show()

del311()
