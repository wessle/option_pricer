import numpy as np

# price a European or American call or put option
def option():
    EA = int(input("Enter 0 for European, 1 for American: "))
    CP = int(input("Enter 0 for call, 1 for put: "))
    S = float(input("Current price: "))
    K = float(input("Strike: "))
    T = float(input("Expiration: "))
    n = int(input("Number of time steps: "))
    sigma = float(input("Volatility: "))
    r = float(input("Risk-free interest rate: "))
    a = float(input("Continuous dividend rate: "))

    # declare time step length, upward move, downward move, and upward move probability
    t = T/n

    # Cox, Ross, Rubinstein definition for u, d
    u = np.exp(sigma*np.sqrt(t))
    d = 1/u
    p = (np.exp((r-a)*t) - d)/(u - d)

    # will use these to store adjacent time steps of the tree
    l0, l1 = [], []

    # initialize the list of spots and option payoffs at maturity
    for i in range(n+1):
        spot = S*(u**(n-i))*(d**i)
        # max(spot-K, 0) since it's a call option
        l1.append([spot, (1-CP)*max(spot-K, 0) + CP*max(K-spot, 0)])

    # backward induction to arrive at current option price
    for i in range(n):
        l0 = []
        k = len(l1) - 1

        for j in range(k):
            spot = S*(u**(k-1-j))*(d**j)
            price = np.exp(-r*t)*(p*l1[j][1] + (1-p)*l1[j+1][1])
            l0.append([spot, (1-CP)*max(price, EA*(spot - K)) + CP*max(price, EA*(K - spot))])
        l1 = l0

    return l1[0][1]
