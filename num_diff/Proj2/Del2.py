import numpy as np
from scipy import linalg as lg
from scipy import sparse as sc
from scipy.sparse import linalg as slg
import matplotlib.pyplot as plt
import math
from plotfunction import plot

def errVSN(N):
    steps = np.array([None]*N)
    analyical_eigs = np.real([(-(math.pi/2 + k*math.pi)**2) for k in range(3)])   #Exact values
    errors = [[None]*N, [None]*N, [None]*N]
    #Beräkna alla egenvärden för olika N
    #Spara
    for i in range(3, N+3):
        step = 2**i
        numerical_eigs = slsolver(step)[0][0:3] #First 3 numerical eigenvalues
        errors[0][i-3] = numerical_eigs[0] - analyical_eigs[0]
        errors[1][i-3] = numerical_eigs[1] - analyical_eigs[1]
        errors[2][i-3] = numerical_eigs[2] - analyical_eigs[2]
        steps[i-3] = step

    slope = steps**(-2)
    plot([steps, steps, steps, steps], [slope, errors[0], errors[1], errors[2]], "Amount of steps", "Error", ["loglog", "loglog", "loglog", "loglog"], True, "Task21a", ["slope = -2", "First eigenvalue", "Second eigenvalue", "Third eigenvalue"])

def slsolver(N):
    subdiag = [1.0]*(N-1)
    subdiag[-1] = 2 #Change S with approx for y_N+1
    deltax = 1/N
    diags = [subdiag, [1.0]*(N-1)]
    S = 1/(deltax**2) * sc.diags(diags, [-1, 1], format = 'csc')
    [numerical_eigs, eigmodes] = lg.eig(S.toarray())
    numerical_eigs = np.add(numerical_eigs, -2 / deltax**2)
    #Numerical eigenvalues
    ind = np.argsort(abs(numerical_eigs))
    numerical_eigs = numerical_eigs[ind]
    eigmodes = eigmodes[:, ind]
    eigmodes = [(np.insert(eigmodes[:, i], 0, 0)) for i in range(N-2)]
    return [numerical_eigs, eigmodes]

def Vcalc(grid):
    #return 700 * (0.5 - abs(grid - 0.5))
    return 800 * np.sin(math.pi*grid)**2    #Double potenital well
    #return np.array([0] * grid)
    #return 800 * np.sin(2*math.pi*grid)**2

def schrodinger(V, N):
    deltax = 1/(N+1)
    subarray = np.array([1.0]*(N-1))
    mainarray = np.array([-2.0]*N)
    diags = 1 / deltax**2 * np.array([subarray, mainarray, subarray])
    diags[1] = np.subtract(diags[1], V)
    A = sc.diags(diags, [-1, 0, 1], format = 'csc')
    [eigs, eigvecs] = lg.eig(A.toarray())
    #[eigs, eigvecs] = slg.eigs(A, k = N - 2) #kanske fråga???
    #Sorting corresponding eigval to eigvec
    ind = np.argsort(abs(eigs))
    eigs = eigs[ind]
    eigvecs = eigvecs[:, ind]
    eigvecs = [(np.insert(eigvecs[:, i], 0, 0)) for i in range(N-2)]
    eigvecs =  [(np.append(eigvecs[i], 0)) for i in range(N-2)]
    return [eigs, eigvecs]

def task21a():
    N = 7
    errVSN(N)

def task21b(): #Plot of the eigenmodes
    print(slsolver(499)[0][0:3]) #[-2.46739906, -22.20644486, -61.68375408]
    [eigs, eigmodes] = slsolver(499)
    grid = np.linspace(0, 1, 500)
    eigenmode = ["First eigenmode", "Second eigenmode", "Third eigenmode"]
    [plt.plot(grid, eigmodes[i], label = eigenmode[i]) for i in range(3)]
    plt.xlabel("x")
    plt.ylabel("u")
    plt.grid(True)
    plt.legend()
    plt.savefig('2.1b.pdf')
    plt.show()

def task22():
    N = 100
    k = 7
    grid = np.linspace(0, 1, N+2)
    V = Vcalc(grid[1:-1])
    plt.plot(grid[1:-1], V, label = "Potential")
    [eigs, eigvecs] = schrodinger(V, N)[0:k]
    N = 60

    #Compute first 8 eigenvalues, double well
    #[(-253.30960830505052+0j), (-253.31175957658948+0j), (-550.1734224425977+0j), (-551.055885276949+0j), (-763.0068000736144+0j),
    #(-799.865424295113+0j), (-931.1500538768632+0j), (-1061.0488142024408+0j)]

    #Triple well eigenvalues
    #[(-167.01433429042908+0j), (-471.28187992229084+0j), (-475.64438840162734+0j), (-480.4019441917935+0j),
    #(-720.8103191653485+0j), (-857.0987037819682+0j), (-931.7621463745484+0j), (-1054.0856783418078+0j)]
    print([eigs[i] for i in range(8)])

    #Wave functions
    legend1 = ["$\Psi_1$", "$\Psi_2$", "$\Psi_3$", "$\Psi_4$", "$\Psi_5$", "$\Psi_6$", "$\Psi_7$", "$\Psi_8$"]
    [plt.plot(grid, eigvecs[i] * N - eigs[i], label = legend1[i]) for i in range(k)]
    plt.xlabel("x")
    plt.ylabel("u")
    plt.grid(True)
    plt.legend()
    plt.savefig('2.2adoublewell.pdf')
    #plt.savefig('2.2atriplewell.pdf')
    plt.show()

    #Probability
    plt.plot(grid[1:-1], V, label = "Potential")
    legend2 = ["|$\Psi_1$|$^2$", "|$\Psi_2$|$^2$", "|$\Psi_3$|$^2$", "|$\Psi_4$|$^2$", "|$\Psi_5$|$^2$", "|$\Psi_6$|$^2$", "|$\Psi_7$|$^2$", "|$\Psi_8$|$^2$"]
    [plt.plot(grid, abs(eigvecs[i])**2 * N**2 - eigs[i], label = legend2[i]) for i in range(k)]
    plt.xlabel("x")
    plt.ylabel("u")
    plt.grid(True)
    plt.legend()
    plt.savefig('2.2bdoublewell.pdf')
    #plt.savefig('2.2btriplewell.pdf')
    plt.show()

task22()
