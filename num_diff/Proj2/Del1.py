import numpy as np
from scipy import linalg as lg
from scipy import sparse as sc
from scipy.sparse import linalg as slg
import matplotlib.pyplot as plt
import math
from plotfunction import plot

def calcf(grid):
    return -np.cos(grid)

def twopBVP(fvec, alpha, beta, L, N):
    deltax = L/(N+1)
    fvec[0] += -alpha/ (deltax**2)
    fvec[-1] += -beta/ (deltax**2)
    diaghl = [1]*(N-1)
    diagm = [-2]*N
    diagonals = sc.diags([diaghl, diagm, diaghl], [-1, 0, 1], format = 'csc')
    y = slg.spsolve(diagonals, fvec * (deltax**2))
    y = np.insert(y, 0, alpha)
    y = np.append(y, beta)
    return y

def errVSdeltax(alpha, beta, x0, L, N):
    rmsnorms = [None]*N
    deltax = np.array([None]*N)
    for i in range(N):
        n = 2**i
        grid = np.linspace(x0, L, n+2)
        exact = np.cos(grid)
        fvec = calcf(grid[1:-1])
        approx = twopBVP(fvec, alpha, beta, L, n)
        error = approx - exact
        deltax[i] = L / (n + 1)
        rmsnorms[i] = math.sqrt(deltax[i]) * lg.norm(error)

    #(x, y, xlabel, ylabel, scale, save, savename = "1", plotlabel = None
    slope = deltax**2 * rmsnorms[-1]/deltax[-1]**2
    plot([deltax, deltax], [slope, rmsnorms], "Step size, $\Delta$x", "RMSnorm of error", ["loglog", "loglog"], True, "Task11", ["slope = 2", "RMSnorm"])

def task11():
    x0 = 0
    L = 10
    N = 10
    alpha = math.cos(x0)
    beta = math.cos(L)
    #beta = 3
    errVSdeltax(alpha, beta, x0, L, N)

    #Flexa på Martin
    # grid = np.linspace(x0,L,N+2)
    # fvec = calcf(grid[1:-1])
    # y = twopBVP(fvec, alpha, beta, L, N)
    # plt.plot(grid, y)
    # plt.show()

def task12():
    def I(grid, L):
        return 1e-3*(3-2*np.cos(math.pi*grid/L)**12) #[m^4]
    def calcq(grid, q):
        return np.array([q]*len(grid))
    def calcubis(M, I, E):
        return M/(I*E)
    x0 = 0
    L = 10 #[m]
    M0 = 0
    Mf = 0
    u0 = 0
    uf = 0
    N = 999
    E = 1.9e11 #[Nm^-2]
    q = -50e3 #[N/m]
    grid = np.linspace(x0, L, N+2)
    qvec = calcq(grid[1:-1], q)
    Mvec = twopBVP(qvec, M0, Mf, L, N)
    Ivec = I(grid, L)
    ubisvec = calcubis(Mvec, Ivec, E)
    u = twopBVP(ubisvec[1:-1], u0, uf, L, N)
    print("{:.6f}".format(min(u)*1000) + "[mm]")
    plot(grid, u, "x [m]", "Deflection [m]", "plot", True, "task12")

task12()
