import numpy as np
from scipy import linalg as lg
import matplotlib.pyplot as plt

def eulerstep(A, uold, h):
    return uold + h*np.matmul(A, uold)

def eulerint(A, y0, t0, tf, N):

    u = np.array([None]*(N+1))
    exact = np.array([None]*(N+1))
    u[0] = y0
    exact[0] = y0
    h = (tf-t0)/N
    for i in range(N):
        u[i+1] = eulerstep(A, u[i], h)
        exact[i+1] = np.matmul(lg.expm(A*(t0 + h*(i+1) - t0)), y0)

    err = u - exact
    return [u, err]

def errVSh(A, y0, t0, tf):
    N = 10
    norms = [None]*N
    h = [None]*N
    for n in range(N):
        error = eulerint(A, y0, t0, tf, 2**n)[1][-1]
        norms[n] = lg.norm(error)
        h[n] = (tf - t0)/2**n
    plt.loglog(h,norms)
    #plt.show()

N = 1000
t0 = 0
tf = 10
y0 = np.transpose(np.matrix([1,1]))
A = np.matrix([[-1, 10], [0, -3]])
t = np.linspace(t0, tf, N+1)

errVSh(A,y0, t0, tf)
plt.show()

# err = eulerint(A, y0, t0, tf, N)[1]
# rel = np.array([None]*(N+1))
#
# for n in range(N+1):
#     rel[n] = lg.norm(err[n])/(lg.norm(np.matmul(lg.expm(A*(t[n] - t0)), y0)))
#
# plt.semilogy(t, rel)
# plt.show()
