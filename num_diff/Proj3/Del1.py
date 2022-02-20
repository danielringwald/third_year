import numpy as np
from scipy import linalg as lg
from scipy import sparse as sc
from scipy.sparse import linalg as slg
import matplotlib.pyplot as plt
import math
import scipy as sp
from mpl_toolkits import mplot3d as p3d
import matplotlib.cm as cm
#from plotfunction import plot

def g(grid):
    #return 500*grid**2 * (1 - grid)
    return math.sqrt(2)*np.sin(math.pi*grid)

def eulerstep(Tdx, uold, dt):
    u = uold[1:-1] + dt * np.matmul(Tdx, uold[1:-1])
    u = np.insert(u, 0, 0)
    u = np.append(u, 0)
    return u

def TRstep(Tdx, uold, dt):
    A = np.eye(len(uold)-2) - dt/2*Tdx
    B = np.matmul(np.eye(len(uold)-2) + dt/2*Tdx, uold[1:-1])
    u = lg.solve(A, B)
    u = np.insert(u, 0, 0)
    u = np.append(u, 0)
    return u

def diffsolve(N, M, tend, eul):
    deltax = 1/(N+1)
    dt = 1/M
    xx = np.linspace(0, 1, N+2) #Grid
    tt = np.linspace(0, tend, M+1)
    subdiag = 1/deltax**2 * np.array([1]*(N-1))
    mdiag = 1/deltax**2 * np.array([-2]*N)
    Tdx = np.diag(subdiag, -1) + np.diag(subdiag, 1) + np.diag(mdiag)
    [T, X] = np.meshgrid(tt, xx)
    U = [None]*(M+1)
    U[0] = g(xx).tolist()

    if eul:
        for i in range(M):
            U[i+1] = eulerstep(Tdx, np.array(U[i]), dt).tolist()
    else:
        for i in range(M):
            U[i+1] = TRstep(Tdx, np.array(U[i]), dt).tolist()
    U = np.array(U)
    U = np.transpose(U)

    fig = plt.figure(figsize = (15, 15))
    ax = plt.axes(projection = '3d')
    ax.plot_surface(T, X, U, cmap=cm.plasma)
    plt.show()

def task11():
    eul = True
    N = 10
    tend = 50

    #Not nice plot
    M = 221
    diffsolve(N, M, tend, eul) #Here we can start to see numerical instability
    #dt/dx^2 = 0.5475113122171946

    #Nice plot
    M = 500
    diffsolve(N, M, tend, eul)

def task12():
    eul = False
    N = 10
    M = 100
    tend = 25
    diffsolve(N, M, tend, eul)

task11()
