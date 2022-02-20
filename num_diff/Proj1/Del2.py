import numpy as np
from scipy import linalg as lg
import matplotlib.pyplot as plt
from plotfunction import plot

class RK4:
    A = np.array([0, 1/2, 1/2, 1])          #How Runge-Kutta order 4 looks like
    B = np.array([1/6, 1/3, 1/3, 1/6])
    C = np.array([0, 1/2, 1/2, 1])

class diffeq:
    Asys = [None]
    def calclinear(A):
        diffeq.Asys = A
        def linear(t, u):
            return np.matmul(A, u)   #y' = Ay, for task 2.1-4
        return linear

def RK34step(func, told, uold, h):
    yprimes = np.array([np.zeros(len(uold))] * 5)
    for n in range(4):
        yprimes[n+1] = func(told, uold + h * RK4.A[n] * yprimes[n])
    zprime = func(told, uold - h*yprimes[1] + 2*h*yprimes[2])   #Embedded RK
    err = h/6*(2*yprimes[2] + zprime - 2*yprimes[3] - yprimes[4])
    unew = uold

    for i in range(4): #Add the weigthed Y-primes to a new array which is what it will evaluate to if we take one step
        unew = unew + h*RK4.B[i]*yprimes[i+1]

    return [unew, err]

def RK34int(func, y0, t0, tf, N):
    h = (tf - t0)/N
    u = y0
    for i in range(N):
        u = RK34step(func, t0 + h*i, u, h)[0]    #Only want unew from RK4step
    exact = np.matmul(y0, lg.expm(diffeq.Asys*(tf-t0)))
    err = u - exact

    return [u, err]

def errVSh(func, y0, t0, tf, N):
    norms = [None]*N
    h = np.array([None]*N)
    for n in range(N):
        error = RK34int(func, y0, t0, tf, 2**n)[1]
        norms[n] = lg.norm(error)
        h[n] = (tf - t0)/2**n
    #plt.title("$y_0$ = " + str(y0))
    plot([h, h], [norms, h**4 * norms[0]/(h[0]**4)], "Step size", "Size of error", ["loglog", "loglog"], True, "21", ["Numeric", "Reference slope = 4"])

def newstep(tol, err, errold, hold, k):
    hnew = abs(tol/err)**(2/(3*k)) * (tol/errold)**(-1/(3*k)) * hold    #Equation (2) from theory
    return hnew

def adaptiveRK34(func, t0, tf, y0, tol): #Changed name to func to avoid duplicate names
    hnew = abs(tf-t0)*np.power(tol, 1/4) / (100*(1 + lg.norm(func(t0, y0))))
    y = [y0]
    t = [t0]
    errold = tol
    k = 4   #Order k of the error estimator

    while(hnew < (tf-t[-1])):
        step = RK34step(func, None, y[-1], hnew)
        t.append(t[-1] + hnew)
        y.append(step[0])
        err = lg.norm(step[1])
        hnew = newstep(tol, err, errold, hnew, k)
        errold = err
    hnew = tf - t[-1]
    y.append(RK34step(func, None, y[-1], hnew)[0])
    t.append(tf)
    y = np.array(y)
    # t = np.array(t)
    return [t, y]

### Test
def del21():
    t0 = 0
    tf = 1
    N = 15
    y0 = np.array([1])
    A = np.matrix([-1])
    errVSh(diffeq.calclinear(A), y0, t0, tf, N)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   del21()
