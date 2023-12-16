## Chinmay G Khasanis
import numpy as np

##------------- Euler's method -------------##
def euler_ode(X0, t0, tN, N, F):
    """
    euler_ode(x0, t0, tN, N, fun)
    x0      : x(t0)   
    [t0,tN] : Interval in which x(t) has to be found
    N       : Number of solution points
    fun     : f(x) = dx/dt
    """
    dimension  = X0.shape[0]  # size of the phase space vector
    X = np.zeros([dimension,N]) ; X[:,0]=X0  # Initialise the initial value
    t = np.linspace(t0, tN, N)
    h = (tN-t0)/N

    # 1st order Taylor expn for each step
    for i in range(N-1):
        X[:,i+1] = X[:,i] + h*F(X[:,i], t[i])
    return X,t

##------------- Three-term recursion -------------##
def three_term(x0, t0, tN, N, fun):
    """
    x0      : x(t0)   
    [t0,tN] : Interval in which x(t) has to be found
    N       : Number of solution points
    fun     : f(x) = dx/dt
    """
    dimension  = X0.shape[0]  # size of the phase space vector
    X = np.zeros([dimension,N]) ; X[:,0]=X0  # Initialise the initial value
    t = np.linspace(t0, tN, N)
    h = (tN-t0)/N
    X[:,1] =X[:,0] + h*fun(t[:,0], X[:,0])
    # 1st order Taylor expn for each step
    for i in range(1,N-1):
        X[:,i+1] = X[:,i] + h*F(X[:,i], t[i])
    return X,t


##------------- Runge-Kutta methods -------------##

##------------- 2nd order -------------##
def RK2(X0, t0, tN, N, fun):
    """
    RK2(X0, t0, tN, N, fun)
    x0      : x(t0)   
    [t0,tN] : Interval in which x(t) has to be found
    N       : Number of solution points
    fun     : f(x) = dx/dt
    """
    dimension  = X0.shape[0]  # size of the phase space vector
    X = np.zeros([dimension,N]) ; X[:,0]=X0  # Initialise the initial value
    t = np.linspace(t0, tN, N)
    h = (tN-t0)/N
    for i in range(N-1):
        K1 = h * fun(t[i], X[:,i])
        K2 = h * fun(t[i] + h, X[:,i] + K1)
        X[:,i+1] = X[:,i] + 0.5 * (K1 + K2)
    return X,t

##------------- 4th order -------------##
def RK4(X0, t0, tN, N, F):
    """
    RK4(X0, t0, tN, N, F)
    X0      : X(t0)
    [t0,tN] : Interval in which X(t) has to be found
    N       : Number of solution points
    Fun     : F(x) = dX/dt
    """
    dimension  = X0.shape[0] # size of the phase space vector
    X = np.zeros([dimension,N]) ; X[:,0]=X0  # Initialise the initial value
    t = np.linspace(t0, tN, N)
    h = (tN-t0)/N
    for i in range(N-1):
        K1 = h * F(t[i], X[:,i])
        K2 = h * F(t[i] + 0.5 * h, X[:,i] + 0.5 * K1)
        K3 = h * F(t[i] + 0.5 * h, X[:,i] + 0.5 * K2)
        K4 = h * F(t[i] + h, X[:,i] + K3)
        X[:,i+1] = X[:,i] + (K1 + 2 * K2 + 2 * K3 + K4) / 6
    return X,t

def RK4_mod(X0, t0, tN, N, F):
    """
    X0      : X(t0)
    [t0,tN] : Interval in which X(t) has to be found
    N       : Number of solution points
    Fun     : F(x) = dX/dt
    Returns the solution X(t) for those points before hitting boundary
    """
    dimension  = X0.shape[0]  # size of the phase space vector
    q = 1-int(2/dimension)
    X = np.zeros([dimension,N]) ; X[:,0]=X0  # Initialise the initial value
    t = np.linspace(t0, tN, N)
    h = (tN-t0)/N
    for i in range(N-1):
        K1 = h * F(t[i], X[:,i])
        K2 = h * F(t[i] + 0.5 * h, X[:,i] + 0.5 * K1)
        K3 = h * F(t[i] + 0.5 * h, X[:,i] + 0.5 * K2)
        K4 = h * F(t[i] + h, X[:,i] + K3)
        X[:,i+1] = X[:,i] + (K1 + 2 * K2 + 2 * K3 + K4) / 6
        if X[q,i+1] <= 0:  # Stop when it hit boundary
            break
    X = X[:,:i+1]
    t = t[:i+1]
    return X,t
