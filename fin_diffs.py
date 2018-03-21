import numpy as np

# to calculate a European call, just use American call, they should be equal
def op():
    EA = int(input("Enter 0 for European, 1 for American: "))
    CP = int(input("Enter 0 for call, 1 for put: "))
    S = float(input("Current price: "))
    K = float(input("Strike: "))
    T = float(input("Expiration: "))
    sig = float(input("Volatility: "))
    r = float(input("Risk-free interest rate: "))
    q = float(input("Continuous dividend rate: "))
    M = int(input("Number of price steps: "))
    N = int(input("Number of time steps: "))

    # price and time step widths with mult*S as the max stock price
    # must be careful here in case the strike is higher than mult*S
    # mult can't be too large
    mult = 2

    # change M to make it easier to recover option price at the end
    M -= M % mult

    s = mult*S/M
    t = T/N

    # initialize the grid
    g = np.zeros((M+1, N+1))
    g[0] = g[0] + (1-CP)*(mult*S-K)*np.ones((1, N+1))
    g[M] = g[M] + CP*K*np.ones((1, N+1))
    for i in range(1, M):
        g[i][N] = (1-CP)*max((M-i)*s-K, 0) + CP*max(K-(M-i)*s, 0)

    g[M//mult][0] = S
    g[M//mult][N] = (1-CP)*max(S-K, 0) + CP*max(K-S, 0)

    # define the coefficients to be used in solving difference equations
    # note that we need to switch the j indices around: in my head I'm thinking
    # of lower indices corresponding to higher up on the y-axis, just like
    # in the x-y plane
    def a(j):
        return -0.5*(r-q)*(M-j)*t - 0.5*sig**2*(M-j)**2*t

    def b(j):
        return 1 + sig**2*(M-j)**2*t + r*t

    def c(j):
        return 0.5*(r-q)*(M-j)*t - 0.5*sig**2*(M-j)**2*t

    # work backwards from the righthand side, updating grid column-by-column
    for i in reversed(range(1,N+1)):
        A = np.identity(M+1)
        d = g[:,i]

        for j in range(1, M):
            A[j][j-1] = a(j)
            A[j][j] = b(j)
            A[j][j+1] = c(j)

        f = np.linalg.solve(A, d)
        f = [ ( (1-CP)*max(f[k], EA*((M-k)*s - K)) + CP*(max(f[k], EA*(K-(M-k)*s))) ) for k in range(M+1)]

        g[:,i-1] = f
        g[0,i-1] = (1-CP)*(mult*S-K)
        g[M,i-1] = CP*K

    print(g[M//mult][0])
