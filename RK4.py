def RK4(X0, t0, tN, N, F):
    """
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
