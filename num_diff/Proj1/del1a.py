import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

def eulerstep(A, uold, h):
    return uold + h * np.matmul(A, uold)

def eulerint(A, y0, t0, tf, N):
    u = np.array([None] * (N+1))
    exact = np.array([None] * (N+1))
    u[0] = y0
    exact[0] = y0
    h = (tf - t0) / N
    for i in range(N):
        u[i+1] = eulerstep(A, u[i], h)
        exact[i+1] = np.matmul(linalg.expm(A*(h*(i+1))), y0)

    err = u - exact

    return [u, err]

def errVSh(A, y0, t0, tf):
    Nmax = 15
    norms = [None] * Nmax
    h = [None] * Nmax
    for N in range(Nmax):
        err = eulerint(A, y0, t0, tf, 2**N)[1][-1]
        norms[N] = linalg.norm(err)
        h[N] = (tf - t0) / 2**N

    plt.loglog(h, norms)
    #plt.show()

y0 = np.transpose(np.matrix([[1, 1]]))
t0 = 0
tf = 10
A = np.matrix([[-1, 10], [0, -3]])
N = 20

t = np.linspace(t0, tf, N+1)

errVSh(A, y0, t0, tf)
plt.show()

## Del 1.5

# toterr = eulerint(A, y0, t0, tf, N)[1]
# exact = [None]*(N+1)
# errnorm = [None]*(N+1)
# j = 0
# for T in t:
#     exact[j] = linalg.norm(np.matmul(linalg.expm(A*(T - t0)), y0))
#     errnorm[j] = linalg.norm(toterr[j])
#     j += 1
#
# plt.semilogy(t, np.divide(errnorm, exact))
# plt.show()
