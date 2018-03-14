import numpy as np

# price a European or American call or put option
def option():
    EA = int(input("Enter 0 for European, 1 for American: "))
    CP = int(input("Enter 0 for call, 1 for put: "))
    S = float(input("Current price: "))
    K = float(input("Strike: "))
    T = float(input("Expiration: "))
    sigma = float(input("Volatility: "))
    r = float(input("Risk-free interest rate: "))
    a = float(input("Continuous dividend rate: "))
    n = int(input("Number of time steps: "))

    # declare time step length, upward move, downward move, and upward move probability
    t = T/n

    # definition os u, d, pu, pd, pm given in Hull's 'Options, Futures, and
    # Other Derivatives,' p. 409
    u = np.exp(sigma*np.sqrt(3*t))
    d = 1/u
    pu = np.sqrt(t/(12*sigma**2))*(r - a - 0.5*sigma**2) + 1/6
    pd = -np.sqrt(t/(12*sigma**2))*(r - a - 0.5*sigma**2) + 1/6
    pm = 2/3

    # Below are the definitions of u, d, pu, pd, pm given by Cox, Ross, Rubinstein
    # u = np.exp(sigma*np.sqrt(2*t))
    # d = 1/u
    # pu = ((np.exp(0.5*r*t) - np.exp(-sigma*np.sqrt(0.5*t))) \
        # / (np.exp(sigma*np.sqrt(0.5*t)) - np.exp(-sigma*np.sqrt(0.5*t))))**2
    # pd = ((np.exp(sigma*np.sqrt(0.5*t)) - np.exp(0.5*r*t)) \
        # / (np.exp(sigma*np.sqrt(0.5*t)) - np.exp(-sigma*np.sqrt(0.5*t))))**2
    # pm = 1 - pu - pd

    # initialize lists we'll use to store consecutive time steps
    l0, l1 = [], []

    # initialize the possible prices at expiration
    for i in range(n+1):
        spot = S*u**(n-i)
        l1.append([spot, (1-CP)*max(spot-K, 0) + CP*max(K-spot, 0)])

    for i in range(n):
        spot = S*d**(i+1)
        l1.append([spot, (1-CP)*max(spot-K, 0) + CP*max(K-spot, 0)])

    # backward induction to arrive at current option price
    for i in range(n):
        l0 = []

        for j in range(len(l1)-2):
            # current spot is equal to spot price immediately to the right
            # in the next time step
            spot = l1[j+1][0]
            price = np.exp(-r*t)*(pu*l1[j][1] + pm*l1[j+1][1] + pd*l1[j+2][1])
            l0.append([spot, (1-CP)*max(price, EA*(spot - K)) + CP*max(price, EA*(K - spot))])

        l1 = l0

    return l1[0][1]
