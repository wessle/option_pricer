import numpy as np

# we're going to do figure out how to calculate puts first
def put():
    EA = int(input("Enter 0 for European, 1 for American: "))
    # CP = int(input("Enter 0 for call, 1 for put: "))
    S = float(input("Current price: "))
    K = float(input("Strike: "))
    T = float(input("Expiration: "))
    sig = float(input("Volatility: "))
    r = float(input("Risk-free interest rate: "))
    # a = float(input("Continuous dividend rate: "))
    M = int(input("Number of price steps: "))
    N = int(input("Number of time steps: "))

    # price and time step widths with 5*S as the max stock price
    mult = 2
    s = mult*S/M
    t = T/N

    # initialize the grid
    g = np.zeros((M+1, N+1))
    g[M] = K*np.ones((1, N+1))
    for i in range(1, M):
        g[i][N] = max(K-(M-i)*s, 0)

    # define the coefficients to be used in solving difference equations
    # note that we need to switch the j indices around: in my head I'm thinking
    # of lower indices corresponding to higher up on the y-axis, just like
    # in the x-y plane
    def a(j):
        return 0.5*r*(M-j)*t - 0.5*sig**2*(M-j)**2*t

    def b(j):
        return 1 + sig**2*(M-j)**2*t + r*t

    def c(j):
        return -0.5*r*(M-j)*t - 0.5*sig**2*(M-j)**2*t

    # work backwards from the righthand side, updating grid column-by-column
    for i in reversed(range(1,N+1)):
        A = np.identity(M+1)
        d = g[:,i]
        d[0] = 0
        d[M] = K

        for j in range(1, M):
            A[j][j-1] = a(j)
            A[j][j] = b(j)
            A[j][j+1] = c(j)

        f = np.linalg.solve(A, d)
        f = [max(f[k], EA*(K-(M-k)*s)) for k in range(M+1)]

        g[:,i-1] = f
        g[0,i-1] = 0
        g[M,i-1] = K

    # get arrays to print nicely
    np.set_printoptions(precision=3, suppress=2)

    print(g[M//mult][0], g[M//mult][N])
