def RK2_ode(x0, t0, tN, N, fun):
    """
    x0      : x(t0)   
    [t0,tN] : Interval in which x(t) has to be found
    N       : Number of solution points
    fun     : f(x) = dx/dt
    """

    x = np.zeros(N) ; x[0]=x0   # Initialise the initial value
    t = np.linspace(t0, tN, N)
    h = (tN-t0)/N

    # RK-2 steps
    for i in range(N-1):
        k1 = h*fun(x[i], t[i])
        k2 = h*fun(x[i]+0.5*k1, t[i]+0.5*h)
        x[i+1] = x[i] + k2
    return x,t
