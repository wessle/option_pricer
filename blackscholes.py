from scipy.stats import norm
import numpy as np

# use the Black-Scholes formula to price a European call or put
def bs():
    CP = int(input("Enter 0 for call, 1 for put: "))
    S = float(input("Current price: "))
    K = float(input("Strike: "))
    tau = float(input("Expiration: "))
    sigma = float(input("Volatility: "))
    r = float(input("Risk-free interest rate: "))

    d1 = (1/(sigma*np.sqrt(tau)))*(np.log(S/K) + (r+(sigma**2/2))*tau)
    d2 = (1/(sigma*np.sqrt(tau)))*(np.log(S/K) + (r-(sigma**2/2))*tau)

    if CP == 0:
        price = S*norm.cdf(d1) - K*np.exp(-r*tau)*norm.cdf(d2)
    else:
        price = K*np.exp(-r*tau)*norm.cdf(-d2) - S*norm.cdf(-d1)

    return price
