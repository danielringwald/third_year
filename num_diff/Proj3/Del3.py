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

def convdif(u, a, d, dt, dx):
    N = len(u)
    print(d * Tdx - a * Sdx)
    A = np.eye(N) - dt / 2 * (d * Tdx - a * Sdx)
    B = np.matmul(np.eye(N) + dt / 2 * (d * Tdx - a * Sdx), np.array(u))
    return lg.solve(A, B).tolist()


def CDsolver(N, M, tend, a, d):
    dx = 1 / N
    dt = tend / M
    xx = np.linspace(0, 1, N+1)
    U = [None]*(M+1)
    U[0] = g(xx[:-1]).tolist()
    createMatrices(N, dx)
    for m in range(M):
        U[m+1] = convdif(U[m], a, d, dt, dx)
    return U

def createMatrices(N, deltax):
    global Tdx
    global Sdx
    ssdiag = 1/deltax**2 * np.array([1.0]*(N-1))
    mdiag = 1/deltax**2 * np.array([-2.0]*(N))
    Tdx = np.diag(ssdiag, -1) + np.diag(mdiag) + np.diag(ssdiag, 1)
    Sdx = deltax/2 * (np.diag(ssdiag, 1) - np.diag(ssdiag, -1))
    Tdx[0, -1] = Tdx[-1, 0] = 1/deltax**2
    Sdx[0, -1] = -1 / (2*deltax)
    Sdx[-1, 0] = 1 / (2*deltax)

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
    plt.show()

def task31():
    N = 10
    dx = 1/N
    M = 100
    tend = 1
    a = -1
    d = 0.1
    peclet = abs(a/d)
    print(peclet)
    if peclet*dx >= 2:
        print(peclet*(1/N))
        N = math.ceil(peclet/2) * 2
        print(peclet*(1/N))
    U = CDsolver(N, M, tend, a, d)
    plot(N, M, tend, U)

task31()
