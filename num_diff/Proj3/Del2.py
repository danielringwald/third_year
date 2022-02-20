import numpy as np
from scipy import linalg as lg
from scipy import sparse as sc
from scipy.sparse import linalg as slg
import matplotlib.pyplot as plt
import math
import scipy as sp
from mpl_toolkits import mplot3d as p3d
import matplotlib.cm as cm


def g(grid):
    return np.exp(-100*(grid - 0.5)**2)

def LaxWen(u, amu):
    unew = amu/2*(1 + amu)*np.array((u[-1:] + u[:-1])) + (1 - amu**2)*np.array(u) - amu/2*(1 - amu)*np.array((u[1:] + u[:1]))
    return unew.tolist()

def LWsolver(N, M, a):
    deltax = 1/N
    dt = 1/M
    amu = a * dt / deltax
    xx = np.linspace(0, 1, N+1)
    U = [None]*(M+1)
    U[0] = g(xx)[:-1].tolist() #g is vectorized, need for appending in LaxWen

    for i in range(M):
        U[i+1] = LaxWen(U[i], amu)
    return U

def plot(N, M, tend, U):
    U = np.transpose(np.array(U))
    Ulast = U[0,:]
    Unew = np.vstack((U, Ulast))
    xx = np.linspace(0, 1, N+1)
    tt = np.linspace(0, tend, M+1)
    [T, X] = np.meshgrid(tt, xx)
    fig = plt.figure(figsize = (15, 15))
    ax = plt.axes(projection = '3d')
    ax.plot_surface(T, X, Unew, cmap=cm.plasma)
    ax.set_xlabel('Time')
    ax.set_ylabel('Distance')
    ax.set_zlabel('U')
    plt.show()

def task21():
    N = 100
    M = 1000
    a = 2
    tend = 5
    U = LWsolver(N, M, a)
    print( a*(1/M) / (1/N) )
    plot(N, M, tend, U)

def task22(): #Norm
    N = [100, 90]
    M = 100
    a = 1 #trams
    tend = 5
    Us = [np.array(LWsolver(N[0], M, a)), np.array(LWsolver(N[1], M, a))]
    rmsnorms = [np.zeros(M,), np.zeros(M,)]
    tt = np.linspace(0, tend, M)
    for i in range(M):
        rmsnorms[0][i] = math.sqrt(1/M) * lg.norm(Us[0][i,:])
        rmsnorms[1][i] = math.sqrt(1/M) * lg.norm(Us[1][i,:])
    plt.plot(tt, rmsnorms[0], label = '$a\Delta t/\Delta x = 1$')
    plt.title('adx/dt = 1')
    plt.ylabel('$L^2$ norm of U')
    plt.xlabel('Time')
    plt.show()

    plt.plot(tt, rmsnorms[1], label = '$a\Delta t/\Delta x = 0.9$')
    plt.title('adx/dt = 0.9')
    plt.ylabel('$L^2$ norm of U')
    plt.xlabel('Time')
    plt.show()

task21()
