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
    return 6*np.exp(-100*(grid - 0.5)**2)

def LaxWen(u, dt):
    uprime = np.matmul(Sdx, u)
    ubis = np.matmul(Tdx, u)
    LW = u + dt* (-u) * uprime + dt**2 / 2 * (2 * u * uprime**2 + u**2 * ubis)
    return LW

def Burgers(u, d, dt):
    N = len(u)
    A = np.eye(N) - d * dt / 2 * Tdx
    B = LaxWen(u, dt) + d * dt / 2 * np.matmul(Tdx, u)
    unew = lg.solve(A, B)
    return unew.tolist()  #tolist for adding later u_n+1

def VBsolver(N, M, tend, d):
    #dx = 1/N
    dt = 1/M
    createMatrices(N, 1/N)
    xx = np.linspace(0, 1, N+1)
    U = [None] * (M+1)
    U[0] = g(xx[:-1]).tolist()
    for m in range(M):
        U[m+1] = Burgers(np.array(U[m]), d, dt)
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

def plot(N, M, tend, U, savename):
    U = np.transpose(np.array(U))
    Ulast = U[0,:]
    Unew = np.vstack((U, Ulast))
    xx = np.linspace(0, 1, N+1)
    tt = np.linspace(0, tend, M+1)
    [T, X] = np.meshgrid(tt, xx)
    fig = plt.figure(figsize = (15, 15))
    ax = plt.axes(projection = '3d')
    ax.plot_surface(T, X, Unew, cmap=cm.plasma)
    ax.set_title('N = 250, M = 1000, d = 0.01')
    ax.set_xlabel('Time')
    ax.set_ylabel('Distance')
    ax.set_zlabel('U')
    plt.savefig(savename + '.pdf')
    plt.show()

def task41():
    N = 250
    M = 1000
    tend = 1
    d = 0.01
    U = VBsolver(N, M, tend, d)
    plot(N, M, tend, U, "41")

task41()
