def euler_ode(x0, t0, tN, N, fun):
    """
    x0      : x(t0)   
    [t0,tN] : Interval in which x(t) has to be found
    N       : Number of solution points
    fun     : f(x) = dx/dt
    """

    x = np.zeros(N) ; x[0]=x0  # Initialise the initial value
    t = np.linspace(t0, tN, N)
    h = (tN-t0)/N

    # 1st order Taylor expn for each step
    for i in range(N-1):
        x[i+1] = x[i] + h*fun(x[i], t[i])
    return x,t
